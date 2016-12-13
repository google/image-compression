// Copyright 2015 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "image_compression/public/pvrtc_compressor.h"

#include <stddef.h>
#include <algorithm>
#include <vector>

#include "base/integral_types.h"
#include "base/logging.h"
#include "image_compression/internal/bit_util.h"
#include "image_compression/internal/color_util.h"
#include "image_compression/public/compressed_image.h"

namespace image_codec_compression {

// Four methods of encoding modulation bits for a block of pixels used by this
// compressor.
// See section 6. in "Texture Compression using Low-Frequency Signal
// Modulation" http://dl.acm.org/citation.cfm?id=844187 .
// Note that the PVRTC format supports a further modulation mode for
// punch-through alpha, which is not used here.
enum ModulationMode {
  // 1 bit per pixel for each pixel.
  kModulationMode1BPP,
  // The following three modes have 2 bits per pixel for every other pixel in
  // the block, in a checkerboard pattern. The mode specifies how to infer color
  // for the intervening pixels.
  kModulationModeAverage4,  // Average the 4 orthoganally connected neighbors.
  kModulationModeVertical,  // Average the 2 vertical neighbors.
  kModulationModeHorizontal  // Average the 2 horizontal neighbors.
};

// Block width and height as for 2BPP PVRTC.
static const uint32 kLog2BlockWidth = 3;
static const uint32 kLog2BlockHeight = 2;
static const uint32 kBlockWidth = (1 << kLog2BlockWidth);
static const uint32 kBlockHeight = (1 << kLog2BlockHeight);

//-----------------------------------------------------------------------------

//
// General helper functions.
//

// Little-endian write.
static unsigned char *Append32(uint32 value, uint8 *output) {
  *output++ = value >> 0;
  *output++ = value >> 8;
  *output++ = value >> 16;
  *output++ = value >> 24;
  return output;
}

// Returns true if |x| is a power of two.
static bool IsPowerOfTwo(uint32 x) {
  return (x != 0) && !(x & (x - 1));
}

// A quick measure of how different are two colors for the human eye.
// The bigger the return value, the more different.
static uint32 ColorDiff(Rgba8888 color0, Rgba8888 color1) {
  return std::abs(color0.r - color1.r) + std::abs(color0.g - color1.g)
      + std::abs(color0.b - color1.b) + std::abs(color0.a - color1.a);
}

// Calculates *|x| and *|y| for the Z-order curve value |z|.
static void FromZOrder(uint32 z, uint32 *x, uint32 *y) {
  *x = *y = 0;
  for (size_t j = 0; j < 16; j++) {
    *x |= ((z >> (j * 2 + 1)) & 1) << j;
    *y |= ((z >> (j * 2 + 0)) & 1) << j;
  }
}

// Returns the result of encoding the 8 bits of input down as |bit_depth| bits
// and then decoding back up to 8 bits.
// Encoding will simply preserve only the top |bit_depth| bits.
// Decoding will bitwise-or these bits repeatedly into the result, so that the
// output values range as smoothly as possible from 0 to 255.
static inline uint8 ApplyBitDepthReduction(uint8 input, uint32 bit_depth) {
  uint8 encoding_mask = GetMask(bit_depth) << (8 - bit_depth);
  uint8 encoded_bits = input & encoding_mask;

  DCHECK_GE(bit_depth, 3U);  // Not yet implemented for a bit_depth of 1 or 2.

  uint8 result = encoded_bits | (encoded_bits >> bit_depth);
  if (bit_depth <= 3) {
    // The encoded bits will have to be repeated again for the least significant
    // output bits.
    result |= encoded_bits >> (bit_depth * 2);
  }
  return result;
}

//-----------------------------------------------------------------------------

//
// PVRTC-specific helper functions operating on pixels and blocks of pixels.
//

// Returns the color interpolated between |color0| and |color1| as specified by
// |mod| which can range from 0 to 3:
//   0 = color0
//   1 = 5/8ths color0, 3/8ths color1
//   1 = 3/8ths color0, 5/8ths color1
//   3 = color1
static Rgba8888 ApplyModulation(Rgba8888 color0, Rgba8888 color1, uint32 mod) {
  Rgba8888 result(color0);
  DCHECK_GE(4U, mod);
  switch (mod) {
    case 0:
      // Do nothing; keep result = color0.
      break;
    case 1:
      result.r = (5 * color0.r + 3 * color1.r) / 8;
      result.g = (5 * color0.g + 3 * color1.g) / 8;
      result.b = (5 * color0.b + 3 * color1.b) / 8;
      result.a = (5 * color0.a + 3 * color1.a) / 8;
      break;
    case 2:
      result.r = (3 * color0.r + 5 * color1.r) / 8;
      result.g = (3 * color0.g + 5 * color1.g) / 8;
      result.b = (3 * color0.b + 5 * color1.b) / 8;
      result.a = (3 * color0.a + 5 * color1.a) / 8;
      break;
    case 3:
      result = color1;
      break;
  }
  return result;
}

// Returns which modulation (from 0 through 3) best represents |color| given
// the color palette |color0| and |color1|.
static uint32 BestModulation(Rgba8888 color, Rgba8888 color0, Rgba8888 color1) {
  uint32 diff = ColorDiff(color, color0);
  uint32 best_diff = diff;
  uint32 best_mod = 0;

  for (unsigned int current_mod = 1; current_mod < 4; ++current_mod) {
    Rgba8888 current_color = ApplyModulation(color0, color1, current_mod);
    diff = ColorDiff(color, current_color);
    if (diff < best_diff) {
      best_diff = diff;
      best_mod = current_mod;
    } else {
      // If it's not getting better here, it won't get better later.
      return best_mod;
    }
  }

  return best_mod;
}

// Returns a color bilinearly interpolated between the four input colors.
// |px| ranges from 0 (pure |color00| or |color01|) to
//      kBlockWidth (pure |color10| or color11|).
// |py| ranges from 0 (pure |color00| or |color10|) to
//      kBlockHeight (pure |color01| or |color11|).
static Rgba8888 Interpolate4_2BPP(Rgba8888 color00, Rgba8888 color01,
                                  Rgba8888 color10, Rgba8888 color11,
                                  uint32 px, uint32 py) {
  // Calculate the weights that should be applied to the four input colors.
  const uint32 a = (kBlockHeight - py) * (kBlockWidth - px);
  const uint32 b = (kBlockHeight - py) * px;
  const uint32 c = py * (kBlockWidth - px);
  const uint32 d = py * px;
  // Apply these weights.
  const uint32 downscale = kBlockWidth * kBlockHeight;
  return Rgba8888(
      (a * color00.r + b * color01.r + c * color10.r + d * color11.r) /
          downscale,
      (a * color00.g + b * color01.g + c * color10.g + d * color11.g) /
          downscale,
      (a * color00.b + b * color01.b + c * color10.b + d * color11.b) /
          downscale,
      (a * color00.a + b * color01.a + c * color10.a + d * color11.a) /
          downscale);
}

// Returns the color for a pixel in a bilinearly upscaled version of an input
// image. The input image is upscaled kBlockWidth in width and kBlockHeight in
// height. The bilinear interpolation wraps on all four edges of the image.
// For every block of pixels of size (kBlockWidth * kBlockHeight) in the
// upscaled image, where the top left is (0,0), the pixel at position
// (kBlockWidth / 2, kBlockHeight / 2) will use the uninterpolated
// low-frequency image colors, and the rest will be interpolated.
// |source| the raw pixel data for the input image.
// |width| width of the upscaled image.
// |height| height of the upscaled image.
// |x| and |y| the position of the pixel in the upscaled image.
// According to:
//  https://www.khronos.org/registry/gles/extensions/IMG/IMG_texture_compression_pvrtc.txt
// width and height must be power-of-two.
static Rgba8888 GetInterpolatedColor2BPP(const Rgba8888 *source,
                                         unsigned int width,
                                         unsigned int height,
                                         unsigned x, unsigned y) {
  // The left, top, right and bottom edges of the 2x2 pixel block in the source
  // image that will be used to interpolate. Note that through wrapping (for
  // example) source_left may be to the right of source_right.
  // width and height are power-of-two, so we can use '&' instead of '%'.
  const uint32 source_left =
      ((x - kBlockWidth / 2) & (width - 1)) >> kLog2BlockWidth;
  const uint32 source_top =
      ((y - kBlockHeight / 2) & (height - 1)) >> kLog2BlockHeight;
  const uint32 source_right =
      (source_left + 1) & ((width >> kLog2BlockWidth) - 1);
  const uint32 source_bottom =
      (source_top + 1) & ((height >> kLog2BlockHeight) - 1);

  // The bilinear weights to be used for interpolation.
  const uint32 x_weight = (x + kBlockWidth / 2) & (kBlockWidth - 1);
  const uint32 y_weight = (y + kBlockHeight / 2) & (kBlockHeight - 1);

  const uint32 source_width = width / kBlockWidth;
  const Rgba8888 color00 = source[source_top * source_width + source_left];
  const Rgba8888 color01 = source[source_top * source_width + source_right];
  const Rgba8888 color10 = source[source_bottom * source_width + source_left];
  const Rgba8888 color11 = source[source_bottom * source_width + source_right];

  return Interpolate4_2BPP(color00, color01, color10, color11,
                           x_weight, y_weight);
}

// An ordering for colors roughly based on brightness.
static uint32 ColorBrightnessOrder(Rgba8888 color) {
  return static_cast<uint32>(color.r) + static_cast<uint32>(color.g) +
         static_cast<uint32>(color.b) + static_cast<uint32>(color.a);
}

// Gets two colors that represent extremes of the range of colors within a block
// in a source image. A fast alternative to principal component analysis.
// This function also takes care of the wrapping of the coordinates, i.e. |x0|
// and |y0| can be outside the bounds of the image.
// |image| the source image pixel data.
// |width| source image width (must be a power of two).
// |height| source image height (must be a power of two).
// |x0| left edge of the block to be considered in pixels.
// |y0| top edge of the block to be considered in pixels.
// |out_index_0|, |out_index_1| output colors as indices into |image|.
static void GetExtremesFast(const Rgba8888 *image, uint32 width, uint32 height,
                            uint32 x0, uint32 y0,
                            uint32 *out_index_0, uint32 *out_index_1) {
  // Consider 5 different pairs; lightness, then R, G, B, A axes.
  #define PAIRS 5
  uint32 best_fitness[PAIRS][2];
  uint32 best_index[PAIRS][2];
  for (uint32 i = 0; i < PAIRS; i++) {
    // For each pair of colors, the first must have the lowest possible value
    // for the tested fitness, the second the highest possible; hence
    // initialize "best" with extreme high and low values.
    best_fitness[i][0] = static_cast<uint32>(-1);
    best_fitness[i][1] = 0;
    best_index[i][0] = 0;
    best_index[i][1] = 0;
  }

  for (uint32 y = y0; y < y0 + kBlockHeight; y++) {
    for (uint32 x = x0; x < x0 + kBlockWidth; x++) {
      uint32 x_wrapped = (x + width) & (width - 1);
      uint32 y_wrapped = (y + height) & (height - 1);
      uint32 index = y_wrapped * width + x_wrapped;
      Rgba8888 color = image[index];

      // For the first pair, use the lightness.
      uint32 lightness = (77 * color.r + 150 * color.g + 28 * color.b) / 256;
      if (lightness < best_fitness[0][0]) {
        best_fitness[0][0] = lightness;
        best_index[0][0] = index;
      }
      if (lightness > best_fitness[0][1]) {
        best_fitness[0][1] = lightness;
        best_index[0][1] = index;
      }

      // For the next 4 axes, use the R, G, B or A axis.
      for (int component = 0; component < 4; component++) {
        int output_pair = component + 1;
        const uint8 c = color[component];
        if (c < best_fitness[output_pair][0]) {
          best_fitness[output_pair][0] = c;
          best_index[output_pair][0] = index;
        }
        if (c > best_fitness[output_pair][1]) {
          best_fitness[output_pair][1] = c;
          best_index[output_pair][1] = index;
        }
      }
    }
  }

  // Choose the pair for which the color difference is biggest. This makes the
  // algorithm somewhat principal component-ish.
  uint32 best_pair_diff = 0;
  uint32 best_pair = 0;
  for (uint32 i = 0; i < PAIRS; i++) {
    uint32 diff = ColorDiff(image[best_index[i][0]], image[best_index[i][1]]);
    if (diff > best_pair_diff) {
      best_pair = i;
      best_pair_diff = diff;
    }
  }

  *out_index_0 = best_index[best_pair][0];
  *out_index_1 = best_index[best_pair][1];

  // *out_index_0 should be darker than *out_index_1 for consistency; swap if
  // not.
  if (ColorBrightnessOrder(image[*out_index_1]) <
          ColorBrightnessOrder(image[*out_index_0])) {
    uint32 temp = *out_index_0;
    *out_index_0 = *out_index_1;
    *out_index_1 = temp;
  }
}

// Returns the color that the input color will become after encoding as an "A"
// or "B" color in a PVRTC compressed image (where they are converted to
// 16-bit), and then decoding back to 32-bit.
// This helps the compressor choose correct modulation values once the "A" and
// "B" colors are chosen.
// |is_b| is true if this is the "B" color; "A" and "B" are encoded differently.
static Rgba8888 ApplyColorChannelReduction(Rgba8888 color, bool is_b) {
  if (color.a == 255) {
    color.r = ApplyBitDepthReduction(color.r, 5);
    color.g = ApplyBitDepthReduction(color.g, 5);
    color.b = ApplyBitDepthReduction(color.b, is_b ? 5 : 4);
  } else {
    color.r = ApplyBitDepthReduction(color.r, 4);
    color.g = ApplyBitDepthReduction(color.g, 4);
    color.b = ApplyBitDepthReduction(color.b, is_b ? 4 : 3);
    color.a = ApplyBitDepthReduction(color.a, 3);
  }
  return color;
}

// Encode two colors and a modulation mode into an unsigned int.
// The encoding is as follows, in the direction from MSB to LSB:
// 16 bit |colora|, 15 bit |colorb|, 1 bit |mod_mode|.
// Opaque colors are: 1 bit 1, 5 bit R, 5 bit G, 4/5 bit B.
// Translucent colors are: 1 bit 0, 3 bit A, 4 bit R, 4 bit G, 3/4 bit B.
static unsigned EncodeColors(Rgba8888 colora, Rgba8888 colorb,
                             ModulationMode mode) {
  unsigned value = 0;

  if (colora.a == 255) {
    SetBits(15, 1, 1, &value);
    SetBits(1, 4, colora.b >> 4, &value);
    SetBits(5, 5, colora.g >> 3, &value);
    SetBits(10, 5, colora.r >> 3, &value);
  } else {
    SetBits(15, 1, 0, &value);
    SetBits(1, 3, colora.b >> 5, &value);
    SetBits(4, 4, colora.g >> 4, &value);
    SetBits(8, 4, colora.r >> 4, &value);
    SetBits(12, 3, colora.a >> 5, &value);
  }

  if (colorb.a == 255) {
    SetBits(31, 1, 1, &value);
    SetBits(16, 5, colorb.b >> 3, &value);
    SetBits(21, 5, colorb.g >> 3, &value);
    SetBits(26, 5, colorb.r >> 3, &value);
  } else {
    SetBits(31, 1, 0, &value);
    SetBits(16, 4, colorb.b >> 4, &value);
    SetBits(20, 4, colorb.g >> 4, &value);
    SetBits(24, 4, colorb.r >> 4, &value);
    SetBits(28, 3, colorb.a >> 5, &value);
  }

  SetBits(0, 1, mode == kModulationMode1BPP ? 0 : 1, &value);
  return value;
}

// Works out which modulation mode to use for a given block in an image.
// |image_mod| the modulation information for the image.
// |width| and |height| image_mod pixel dimensions (must be a power of two).
// |block_x| block x coordinate, i.e. ranging from 0 to |width| / kBlockWidth.
// |block_y| block y coordinate, i.e. ranging from 0 to |height| / kBlockHeight.
static ModulationMode CalculateBlockModulationMode(
    const uint8 *image_mod, uint32 width, uint32 height,
    uint32 block_x, uint32 block_y) {
  // A count of how many pixels are best served by modulation values 2 or 3,
  // i.e. intermediate between the extremes of one color or the other.
  uint32 intermediate_value_count = 0;

  // A measure of how much variation between pixels there is horizontally.
  uint32 horizontal_count = 0;

  // A measure of how much variation between pixels there is vertically.
  uint32 vertical_count = 0;

  for (uint32 y = 0; y < kBlockHeight; y++) {
    for (uint32 x = 0; x < kBlockWidth; x++) {
      uint32 index =
          (block_y * kBlockHeight + y) * width + (block_x * kBlockWidth + x);

      if (image_mod[index] == 1 || image_mod[index] == 2)
          intermediate_value_count++;

      // Index of adjacent horizontal pixel in |image_mod|.
      uint32 index_adjacent_horizontal =
          (block_y * kBlockHeight + y) * width +
          ((block_x * kBlockWidth + x + 1) & (width - 1));

      // Index of adjacent vertical pixel in |image_mod|.
      uint32 index_adjacent_vertical =
          ((block_y * kBlockHeight + y + 1) & (height - 1)) * width +
          (block_x * kBlockWidth + x);

      horizontal_count +=
          std::abs(image_mod[index] - image_mod[index_adjacent_vertical]);
      vertical_count +=
          std::abs(image_mod[index] - image_mod[index_adjacent_horizontal]);
    }
  }

  if (intermediate_value_count <= 4)
      return kModulationMode1BPP;

  static const uint32 absolute_threshold = 10;
  static const uint32 ratio_threshold = 2;

  if (vertical_count > absolute_threshold &&
      vertical_count > horizontal_count * ratio_threshold)
      return kModulationModeVertical;
  else if (horizontal_count > absolute_threshold &&
           horizontal_count > vertical_count * ratio_threshold)
      return kModulationModeHorizontal;

  return kModulationModeAverage4;
}

// Calculates the 32 bits of modulation information to store for a given block
// in an image.
// |image_mod| the modulation information for the image.
// |width| and |height| image_mod pixel dimensions.
// |block_x| block x coordinate, i.e. ranging from 0 to |width| / kBlockWidth.
// |block_y| block y coordinate, i.e. ranging from 0 to |height| / kBlockHeight.
// |mode| which modulation mode to use.
uint32 CalculateBlockModulationData(const uint8 *image_mod,
                                    uint32 width, uint32 height,
                                    uint32 block_x, uint32 block_y,
                                    ModulationMode mode) {
  uint32 result = 0;
  uint32 bitpos = 0;
  for (unsigned y = 0; y < 4; y++) {
    for (unsigned x = 0; x < 8; x++) {
      size_t index = (block_y * 4 + y) * width + (block_x * 8 + x);

      if (mode == kModulationMode1BPP) {
        uint32 bit = image_mod[index] / 2;
        SetBits(bitpos, 1, bit, &result);
        bitpos++;
      } else {
        if ((x ^ y) & 1) continue;  // checkerboard
        uint32 bit = image_mod[index];
        // The bits at position 0 (0,0) and at position 20 (4,2) are the ones
        // that use only a single bit for the modulation value, and the other
        // bit for selecting the sub-mode.
        if (bitpos == 0) {
          // The saved bit chooses average-of-4 or "other".
          if (mode == kModulationModeAverage4)
            bit &= 2;
          else
            bit |= 1;
        } else if (bitpos == 20) {
          // The saved bit chooses vertical versus horizontal.
          if (mode == kModulationModeVertical)
            bit |= 1;
          else
            bit &= 2;
        }

        SetBits(bitpos, 2, bit, &result);
        bitpos += 2;
      }
    }
  }
  return result;
}

//-----------------------------------------------------------------------------

//
// Helper functions operating on entire images.
//

// Fills in the two low-resolution images representing the "a" and "b" colors in
// the source |image|.
static void Morph(const Rgba8888 *image, uint32 width, uint32 height,
                  Rgba8888 *outa, Rgba8888 *outb) {
  for (uint32 y = 0; y < height; y += kBlockHeight) {
    for (uint32 x = 0; x < width; x += kBlockWidth) {
      uint32 indexa = 0;
      uint32 indexb = 0;
      GetExtremesFast(&image[0], width, height, x, y, &indexa, &indexb);

      uint32 index_out =
          (y / kBlockHeight) * (width / kBlockWidth) + (x / kBlockWidth);

      outa[index_out] = ApplyColorChannelReduction(image[indexa], false);
      outb[index_out] = ApplyColorChannelReduction(image[indexb], true);
    }
  }
}

// Given a source |image| and two low-resolution "a" and "b" images, creates a
// 2-bits-per-pixel "mod" image, i.e. values between 0 and 3 for each pixel in
// |image|. Each output pixel is stored in a byte in |mod|, which is assumed
// to be pre-allocated.
static void Modulate(const Rgba8888 *image, uint32 width, uint32 height,
                     const Rgba8888 *imagea,
                     const Rgba8888 *imageb,
                     unsigned char *mod) {
  for (uint32 y = 0; y < height; y++) {
    for (uint32 x = 0; x < width; x++) {
      Rgba8888 colora =
          GetInterpolatedColor2BPP(&imagea[0], width, height, x, y);
      Rgba8888 colorb =
          GetInterpolatedColor2BPP(&imageb[0], width, height, x, y);
      *mod++ = BestModulation(*image++, colora, colorb);
    }
  }
}

// Takes the calculated "A" and "B" images, and the modulation information, and
// writes out the data in PVRTC format. Though the input modulation information
// has 2 bits per pixel of the uncompressed original image, this function will
// choose on a per-block basis which of the 4 modulation modes to use.
// |width| and |height| image dimensions.
// |imagea| the "A" image as described in pvrtc_compressor.h.
// |imageb| the "B" image as also described.
// |imagemod| One byte per pixel modulation data, containing 2-bit values.
// |pvr| output pvrtc data, assumed preallocated.
static void Encode(uint32 width, uint32 height,
                   const Rgba8888 *imagea, const Rgba8888 *imageb,
                   const uint8 *image_mod, uint8 *pvr) {
  // Loop through all output blocks.
  for (uint32 i = 0; i < width * height / (kBlockWidth * kBlockHeight); i++) {
    // The blocks are stored in Z-order; calculate the block x and y.
    uint32 block_x = 0;
    uint32 block_y = 0;
    FromZOrder(i, &block_x, &block_y);

    // Calculate which kind of encoding is worth doing for this block.
    ModulationMode mode = CalculateBlockModulationMode(image_mod, width, height,
                                                       block_x, block_y);

    // Given this mode, calculate the 32 bits that represent the block's
    // modulation information.
    uint32 mod_data = CalculateBlockModulationData(image_mod, width, height,
                                                   block_x, block_y, mode);

    // The 32 bits that represent the 2-color palette for this block and mode.
    uint32 color_data = EncodeColors(
        imagea[block_y * (width / kBlockWidth) + block_x],
        imageb[block_y * (width / kBlockWidth) + block_x],
        mode);

    // Write out this information.
    pvr = Append32(mod_data, pvr);
    pvr = Append32(color_data, pvr);
  }
}

// Compresses a given RGBA8888 image to 2BPP PVRTC RGBA.
// |image| source image data.
// |width| and |height| image dimensions.
// |pvr| output pvrtc data, assumed preallocated.
static void CompressPVRTC_RGBA_2BPP(const Rgba8888 *image,
                                    uint32 width, uint32 height,
                                    uint8 *pvr) {
  uint32 low_image_size = (width * height) / (kBlockWidth * kBlockHeight);
  std::vector<Rgba8888> imagea(low_image_size);
  std::vector<Rgba8888> imageb(low_image_size);
  std::vector<uint8> imagemod(width * height);

  Morph(image, width, height, &imagea[0], &imageb[0]);
  Modulate(image, width, height, &imagea[0], &imageb[0], &imagemod[0]);
  Encode(width, height, &imagea[0], &imageb[0], &imagemod[0], pvr);
}

//-----------------------------------------------------------------------------

//
// Public functions.
//

PvrtcCompressor::PvrtcCompressor() {
}

PvrtcCompressor::~PvrtcCompressor() {
}

bool PvrtcCompressor::SupportsFormat(CompressedImage::Format format) const {
  return format == CompressedImage::kRGBA;
}

bool PvrtcCompressor::IsValidCompressedImage(const CompressedImage &image) {
  const CompressedImage::Metadata &metadata = image.GetMetadata();
  return
      metadata.format == CompressedImage::kRGBA &&
      metadata.compressor_name == "pvrtc" &&
      metadata.uncompressed_height >= kBlockHeight &&
      metadata.uncompressed_width >= kBlockWidth &&
      metadata.compressed_width == metadata.compressed_height &&
      IsPowerOfTwo(metadata.uncompressed_height) &&
      IsPowerOfTwo(metadata.uncompressed_width) &&
      metadata.compressed_height == metadata.uncompressed_height &&
      metadata.compressed_width == metadata.uncompressed_width &&
      image.GetDataSize() == ComputeCompressedDataSize(metadata.format,
          metadata.uncompressed_height, metadata.uncompressed_width);
}

size_t PvrtcCompressor::ComputeCompressedDataSize(
    CompressedImage::Format format, uint32 height, uint32 width) {
  return width * height / 4;
}

bool PvrtcCompressor::Compress(CompressedImage::Format format,
                             uint32 height, uint32 width,
                             uint32 padding_bytes_per_row,
                             const uint8 *buffer, CompressedImage *image) {
  if (!buffer || !image || height == 0 || width == 0)
    return false;

  if (!IsPowerOfTwo(width) || !IsPowerOfTwo(height) || width != height)
    return false;

  if (padding_bytes_per_row != 0)
    return false;

  if (width % kBlockWidth != 0 || height % kBlockHeight != 0)
    return false;

  size_t data_size = ComputeCompressedDataSize(format, height, width);
  const CompressedImage::Metadata metadata(format, "pvrtc", height, width,
                                           height, width, 0);
  if (image->OwnsData()) {
    image->CreateOwnedData(metadata, data_size);
  } else {
    // Make sure the external storage has the correct size.
    if (image->GetDataSize() != data_size)
      return false;
    image->SetMetadata(metadata);
  }

  CompressPVRTC_RGBA_2BPP((const Rgba8888*)buffer, width, height,
                          image->GetMutableData());
  return true;
}

bool PvrtcCompressor::Decompress(const CompressedImage &image,
                                 std::vector<uint8> *decompressed_buffer) {
  return false;
}

bool PvrtcCompressor::Downsample(const CompressedImage &image,
                                 CompressedImage *downsampled_image) {
  return false;
}

bool PvrtcCompressor::Pad(const CompressedImage &image, uint32 padded_height,
                          uint32 padded_width, CompressedImage *padded_image) {
  return false;
}

bool PvrtcCompressor::CompressAndPad(CompressedImage::Format format,
                                     uint32 height, uint32 width,
                                     uint32 padded_height, uint32 padded_width,
                                     uint32 padding_bytes_per_row,
                                     const uint8 *buffer,
                                     CompressedImage *padded_image) {
  return false;
}

bool PvrtcCompressor::CreateSolidImage(CompressedImage::Format format,
                                       uint32 height, uint32 width,
                                       const uint8 *color,
                                       CompressedImage *image) {
  return false;
}

bool PvrtcCompressor::CopySubimage(const CompressedImage &image,
                                   uint32 start_row, uint32 start_column,
                                   uint32 height, uint32 width,
                                   CompressedImage *subimage) {
  return false;
}

}  // namespace image_codec_compression
