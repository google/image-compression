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

#include "image_compression/public/dxtc_compressor.h"

#include <algorithm>
#include <bitset>

#include "base/integral_types.h"
#include "base/logging.h"
#include "image_compression/internal/color_util.h"
#include "image_compression/internal/compressor4x4_helper.h"
#include "image_compression/internal/dxtc_const_color_table.h"
#include "image_compression/internal/pixel4x4.h"

namespace image_codec_compression {

//-----------------------------------------------------------------------------

//
// DXT-specific data structures.
//

// DXT1 data representation for a compressed 4x4 block of RGB pixels.
struct Dxt1Block {
  typedef Rgb888 ColorType;

  Dxt1Block() {
  }

  // Constructor for a solid-color block.
  explicit Dxt1Block(const Rgb565 &color) {
    uint16 c = ToUInt16(color);
    c1_lo = c0_lo = c & 0xff;
    c1_hi = c0_hi = c >> 8;
    for (int i = 0; i < 4; ++i)
      color_bits[i] = 0;
  }

  // Constructor for a general block.
  Dxt1Block(const Rgb565 &color0, const Rgb565 &color1, const uint8 bits[4]) {
    uint16  c0 = ToUInt16(color0);
    uint16  c1 = ToUInt16(color1);
    c0_lo = c0 & 0xff;
    c0_hi = c0 >> 8;
    c1_lo = c1 & 0xff;
    c1_hi = c1 >> 8;
    for (int i = 0; i < 4; ++i)
      color_bits[i] = bits[i];
  }

  // Compressed data packed into 8 bytes.
  uint8 c0_lo, c0_hi;
  uint8 c1_lo, c1_hi;
  uint8 color_bits[4];
};

// DXT5 data representation for a compressed 4x4 block of RGBA pixels.
struct Dxt5Block {
  typedef Rgba8888 ColorType;

  Dxt5Block() {
  }

  // Constructor for a block where all pixels have the same alpha value.
  Dxt5Block(const Dxt1Block &block, uint8 alpha) {
    dxt1_block = block;
    alpha1 = alpha0 = alpha;
    for (int i = 0; i < 6; ++i)
      alpha_bits[i] = 0;
  }

  // Constructor for a general block.
  Dxt5Block(const Dxt1Block &block, uint8 base_alphas[2], const uint8 bits[6]) {
    dxt1_block = block;
    alpha0 = base_alphas[0];
    alpha1 = base_alphas[1];
    for (int i = 0; i < 6; ++i)
      alpha_bits[i] = bits[i];
  }

  // Compressed data packed into 8 bytes + a Dxt1Block.
  uint8 alpha0, alpha1;
  uint8 alpha_bits[6];
  Dxt1Block dxt1_block;
};

// This class is used for processing alpha code bits for DXT5.  Alpha
// codes are 3 bits for each of the 16 pixels in a 4x4 block and are
// stored in 6 8-bit bytes. Therefore, we use a 48-bit bitset to make
// the bit addressing easier.
class Dxt5AlphaBits {
 public:
  // The default constructor sets all bits to 0.
  Dxt5AlphaBits() {
  }

  // Constructor that sets all bits from an array of 6 bytes.
  explicit Dxt5AlphaBits(const uint8 bytes[6]) {
    int cur_index = 0;
    for (int i = 0; i < 6; ++i) {
      uint8 cur_byte = bytes[i];
      for (int j = 0; j < 8; ++j) {
        bits_.set(cur_index++, cur_byte & 1);
        cur_byte >>= 1;
      }
    }
  }

  // Sets the 3-bit alpha code for pixel n (where 0 is the lowest-order pixel).
  void SetCode(int n, int code) {
    DCHECK_GE(n, 0);
    DCHECK_LT(n, 16);
    DCHECK_GE(code, 0);
    DCHECK_LT(code, 8);
    const int lsb = n * 3;
    bits_.set(lsb, code & 1);
    bits_.set(lsb + 1, code & 2);
    bits_.set(lsb + 2, code & 4);
  }

  // Returns the 3-bit alpha code for pixel n (where 0 is the
  // lowest-order pixel).
  int GetCode(int n) const {
    DCHECK_GE(n, 0);
    DCHECK_LT(n, 16);
    const int lsb = n * 3;
    return (static_cast<int>((bits_.test(lsb + 2)) << 2) |
            static_cast<int>((bits_.test(lsb + 1)) << 1) |
            static_cast<int>(bits_.test(lsb)));
  }

  // Converts the bitset back to an array of 6 bytes for use in a Dxt5Block.
  void GetBytes(uint8 bytes[6]) const {
    int cur_index = 0;
    for (int i = 0; i < 6; ++i) {
      uint8 cur_byte = 0;
      for (int j = 0; j < 8; ++j) {
        cur_byte |= (bits_.test(cur_index++) << j);
      }
      bytes[i] = cur_byte;
    }
  }

 private:
  std::bitset<48> bits_;
};

//-----------------------------------------------------------------------------

//
// General DXT block decoding functions.
//

// Decodes the 4 possible colors for a Dxt1Block.
static void DecodeColors(const Dxt1Block &block, Rgb888 colors[4],
                         bool swap_red_and_blue, bool always_4_color_case) {
  const uint16 color0_as_uint16 = block.c0_lo + block.c0_hi * 256;
  const uint16 color1_as_uint16 = block.c1_lo + block.c1_hi * 256;

  const Rgb565 color0_as_565 = ToRgb565(color0_as_uint16);
  const Rgb565 color1_as_565 = ToRgb565(color1_as_uint16);

  colors[0] = ExtendToRgb888(color0_as_565);
  colors[1] = ExtendToRgb888(color1_as_565);

  if (swap_red_and_blue) {
    SwapRedAndBlue(&colors[0]);
    SwapRedAndBlue(&colors[1]);
  }

  if (color0_as_uint16 == color1_as_uint16) {
    colors[2] = colors[3] = colors[1];
  } else if (always_4_color_case || color0_as_uint16 > color1_as_uint16) {
    colors[2] = CombineRgb888Fast<2, 1>(colors[0], colors[1]);
    colors[3] = CombineRgb888Fast<1, 2>(colors[0], colors[1]);
  } else {
    colors[2] = CombineRgb888Fast<1, 1>(colors[0], colors[1]);
    colors[3].r = colors[3].g = colors[3].b = 0x00;
  }
}

// Decodes the 8 possible alpha values for a Dxt5Block.
static void DecodeAlphaValues(const Dxt5Block &dxt5_block, uint8 alpha[8]) {
  const uint8 alpha0 = dxt5_block.alpha0;
  const uint8 alpha1 = dxt5_block.alpha1;

  alpha[0] = alpha0;
  alpha[1] = alpha1;

  if (alpha0 > alpha1) {
    alpha[2] = CombineUint8Fast<6, 1>(alpha[0], alpha[1]);
    alpha[3] = CombineUint8Fast<5, 2>(alpha[0], alpha[1]);
    alpha[4] = CombineUint8Fast<4, 3>(alpha[0], alpha[1]);
    alpha[5] = CombineUint8Fast<3, 4>(alpha[0], alpha[1]);
    alpha[6] = CombineUint8Fast<2, 5>(alpha[0], alpha[1]);
    alpha[7] = CombineUint8Fast<1, 6>(alpha[0], alpha[1]);
  } else {
    alpha[2] = CombineUint8Fast<4, 1>(alpha[0], alpha[1]);
    alpha[3] = CombineUint8Fast<3, 2>(alpha[0], alpha[1]);
    alpha[4] = CombineUint8Fast<2, 3>(alpha[0], alpha[1]);
    alpha[5] = CombineUint8Fast<1, 4>(alpha[0], alpha[1]);
    alpha[6] = 0;
    alpha[7] = 255;
  }
}

// Decodes a single Dxt1Block into a 4x4 array of RGB pixels.
void DecodeDxt1Block(const Dxt1Block &dxt1_block, bool swap_red_and_blue,
                     Rgb888 decoded_pixels[4][4]) {
  // Decode the 4 colors defined by the block.
  Rgb888 colors[4];
  DecodeColors(dxt1_block, colors, swap_red_and_blue, false);

  // Decode the color indices and use them to determine pixel colors.
  // The index of the color to use for pixel (0,0) is the lowest-order
  // 2 bits of the first entry in the bits[] array. The next 2 bits up
  // are for pixel (1,0), and so on.
  for (int y = 0; y < 4; ++y) {
    const uint8 bits = dxt1_block.color_bits[y];
    for (int x = 0; x < 4; ++x) {
      uint8 code = (bits >> (2 * x)) & 3;
      decoded_pixels[y][x] = colors[code];
    }
  }
}

// Decodes a single Dxt5Block into a 4x4 array of RGBA pixels.
static void DecodeDxt5Block(const Dxt5Block &block, bool swap_red_and_blue,
                            Rgba8888 decoded_pixels[4][4]) {
  // Decode the color part of the block a la DXT1, but always using
  // the 4-color case.
  Rgb888 colors[4];
  DecodeColors(block.dxt1_block, colors, swap_red_and_blue, true);

  // Decode the alpha values in the block.
  uint8 alpha[8];
  DecodeAlphaValues(block, alpha);

  // Decode the color indices and alpha codes. The color indexing is
  // the same as for DXT1. The alpha value to use for pixel (0,0) is
  // the lowest-order 3 bits of the first entry in the alpha_bits[]
  // array. The next 3 bits up are for pixel (1,0), and so on through
  // the other entries. The Dxt5AlphaBits class simplifies the loop.
  const Dxt5AlphaBits alpha_bits(block.alpha_bits);
  int cur_alpha_index = 0;
  for (int y = 0; y < 4; ++y) {
    const uint8 color_bits = block.dxt1_block.color_bits[y];
    for (int x = 0; x < 4; ++x) {
      const uint8 color_code = (color_bits >> (2 * x)) & 3;
      const uint8 alpha_code = alpha_bits.GetCode(cur_alpha_index);
      decoded_pixels[y][x] = Rgba8888(colors[color_code], alpha[alpha_code]);
      cur_alpha_index++;
    }
  }
}

//-----------------------------------------------------------------------------

//
// General DXT block encoding functions.
//

// Returns the size of a single compressed block used to encode the
// given format.
static size_t GetBlockSize(CompressedImage::Format format) {
  return GetNumFormatComponents(format) == 3 ?
      sizeof(Dxt1Block) : sizeof(Dxt5Block);
}

// Computes the two base colors to use for compressing the given block
// of pixels. The order of the two colors is arbitrary.
static void ComputeBaseColors(const Pixel4x4 &pixel4x4, bool swap_red_and_blue,
                              RgbInt base_colors[2]) {
  // Use the colors with the smallest and largest luminance. This is
  // a fast heuristic that produces reasonable results.
  RgbInt low_color = ToRgbOrBgrInt(pixel4x4.GetPixel(0, 0), swap_red_and_blue);
  RgbInt high_color = low_color;
  if (!pixel4x4.has_one_pixel()) {
    int low_luminance  = kint32max;
    int high_luminance = 0;
    for (int y = 0; y < 4; ++y) {
      for (int x = 0; x < 4; ++x) {
        const RgbInt &color = ToRgbOrBgrInt(pixel4x4.GetPixel(y, x),
                                            swap_red_and_blue);
        const int luminance = ComputeLuminanceFast(color);
        if (luminance < low_luminance) {
          low_luminance = luminance;
          low_color = color;
        }
        if (luminance > high_luminance) {
          high_luminance = luminance;
          high_color = color;
        }
      }
    }
  }
  base_colors[0] = low_color;
  base_colors[1] = high_color;
}

// Sets the color_bits parameter to the encoding for the given pixels
// using the two given base colors.
static void ComputeColorBits(const Pixel4x4 &pixel4x4, bool swap_red_and_blue,
                             const RgbInt base_colors[2], uint8 color_bits[4]) {
  // Special case for single-color blocks.
  if (pixel4x4.has_one_pixel()) {
    color_bits[0] = color_bits[1] = color_bits[2] = color_bits[3] = 0;
    return;
  }

  // Compute the 4 test colors.
  RgbInt test_colors[4];
  test_colors[0] = base_colors[0];
  test_colors[1] = base_colors[1];
  test_colors[2] = CombineRgbIntFast<2, 1>(test_colors[0], test_colors[1]);
  test_colors[3] = CombineRgbIntFast<1, 2>(test_colors[0], test_colors[1]);

  for (int y = 0; y < 4; ++y) {
    color_bits[y] = 0;
    for (int x = 0; x < 4; ++x) {
      const RgbInt color = ToRgbOrBgrInt(pixel4x4.GetPixel(y, x),
                                         swap_red_and_blue);
      int which_color = 0;
      int smallest_luminance_error =
          ComputeSquaredLuminanceDistanceFast(color, test_colors[0]);
      for (int c = 1; c < 4; ++c) {
        const int luminance_error =
            ComputeSquaredLuminanceDistanceFast(color, test_colors[c]);
        if (luminance_error < smallest_luminance_error) {
          smallest_luminance_error = luminance_error;
          which_color = c;
        }
      }
      color_bits[y] |= which_color << (2 * x);
    }
  }
}

// Sets the color0, color1, and color_bits parameters to the best encoding for
// the base_color.
static void ComputeConstantColorBits(const Pixel4x4 &pixel4x4,
                                     bool swap_red_and_blue,
                                     bool always_4_color_case,
                                     const RgbInt& base_color,
                                     Rgb565* color0,
                                     Rgb565* color1,
                                     uint8 color_bits[4]) {
  RgbInt use_color = ToRgbOrBgrInt(base_color, swap_red_and_blue);
  int which_color = GetBestDxtcConstColors(use_color, color0, color1,
                                           always_4_color_case);
  // Replicate the 2 bit which_color for all the color_bits.
  uint8 which_color_byte = which_color | (which_color << 2);
  which_color_byte = which_color_byte | (which_color_byte << 4);
  for (int y = 0; y < 4; ++y) {
    color_bits[y] = which_color_byte;
  }
}

// Computes the two base alpha values to use for compressing the given
// block of pixels. The two values are ordered to specify the correct
// compression scheme according to the DXT5 spec.
static void ComputeBaseAlphas(const Pixel4x4 &pixel4x4, uint8 base_alphas[2]) {
  // Special case for single-color blocks.
  if (pixel4x4.has_one_pixel()) {
    base_alphas[0] = base_alphas[1] = pixel4x4.GetAlpha(0, 0);
    return;
  }

  // Look for the smallest and largest alpha values, also keeping
  // track of how many fully transparent (0) and fully opaque (255)
  // values there are.
  int num_fully_transparent = 0;
  int num_fully_opaque = 0;
  int low_alpha  = 255;
  int high_alpha = 0;
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < 4; ++x) {
      const int alpha = pixel4x4.GetAlpha(y, x);
      if (alpha == 0) {
        num_fully_transparent++;
      } else if (alpha == 255) {
        num_fully_opaque++;
      } else {
        if (alpha < low_alpha)
          low_alpha = alpha;
        if (alpha > high_alpha)
          high_alpha = alpha;
      }
    }
  }

  // Make sure high and low alpha values are reasonable. This happens
  // only when all values are either 0 or 255.
  if (low_alpha > high_alpha) {
    low_alpha  = 0;
    high_alpha = 255;
  }

  // If there are enough fully transparent or opaque pixels, use the
  // scheme that explicitly includes 0 and 255.
  if (num_fully_transparent > 1 || num_fully_opaque > 1) {
    base_alphas[0] = low_alpha;
    base_alphas[1] = high_alpha;
  } else {
    if (num_fully_transparent > 0)
      low_alpha = 0;
    if (num_fully_opaque > 0)
      high_alpha = 255;
    base_alphas[0] = high_alpha;
    base_alphas[1] = low_alpha;
  }
}

// Sets alpha_bits to the encoding for the given pixels.
static void ComputeAlphaBits(const Pixel4x4 &pixel4x4,
                             const uint8 base_alphas[2], uint8 alpha_bits[6]) {
  // Special case for single-color blocks.
  if (pixel4x4.has_one_pixel()) {
    for (int i = 0; i < 6; ++i)
      alpha_bits[i] = 0;
    return;
  }

  int test_alphas[8];
  test_alphas[0] = base_alphas[0];
  test_alphas[1] = base_alphas[1];

  if (base_alphas[0] <= base_alphas[1]) {
    // Scheme that includes 0 and 255 explicitly:
    test_alphas[2] = CombineIntFast<4, 1>(base_alphas[0], base_alphas[1]);
    test_alphas[3] = CombineIntFast<3, 2>(base_alphas[0], base_alphas[1]);
    test_alphas[4] = CombineIntFast<2, 3>(base_alphas[0], base_alphas[1]);
    test_alphas[5] = CombineIntFast<1, 4>(base_alphas[0], base_alphas[1]);
    test_alphas[6] = 0;
    test_alphas[7] = 255;
  } else {
    // Scheme that just interpolates endpoint alphas:
    test_alphas[2] = CombineIntFast<6, 1>(base_alphas[0], base_alphas[1]);
    test_alphas[3] = CombineIntFast<5, 2>(base_alphas[0], base_alphas[1]);
    test_alphas[4] = CombineIntFast<4, 3>(base_alphas[0], base_alphas[1]);
    test_alphas[5] = CombineIntFast<3, 4>(base_alphas[0], base_alphas[1]);
    test_alphas[6] = CombineIntFast<2, 5>(base_alphas[0], base_alphas[1]);
    test_alphas[7] = CombineIntFast<1, 6>(base_alphas[0], base_alphas[1]);
  }

  // Choose the closest alpha for each pixel.
  Dxt5AlphaBits bits;
  int cur_alpha_index = 0;
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < 4; ++x) {
      const int alpha = pixel4x4.GetAlpha(y, x);
      int which_alpha = 0;
      int smallest_error = ComputeSquaredComponentDistance(alpha,
                                                           test_alphas[0]);
      for (int a = 1; a < 8; ++a) {
        const int error = ComputeSquaredComponentDistance(alpha,
                                                          test_alphas[a]);
        if (error < smallest_error) {
          smallest_error = error;
          which_alpha = a;
        }
      }
      bits.SetCode(cur_alpha_index++, which_alpha);
    }
  }
  bits.GetBytes(alpha_bits);
}

// Encodes a 4x4 block of RGB pixels into a Dxt1Block, which is returned.
static Dxt1Block EncodeDxt1Block(const Pixel4x4 &pixel4x4,
                                 bool swap_red_and_blue,
                                 bool always_4_color_case) {
  // Find the colors in the block to use as base colors for the encoding.
  RgbInt base_colors[2];
  ComputeBaseColors(pixel4x4, swap_red_and_blue, base_colors);

  // Convert pixel colors to Rgb565 form and to uint16 for comparisons.
  Rgb565 color0 = QuantizeToRgb565(base_colors[0]);
  Rgb565 color1 = QuantizeToRgb565(base_colors[1]);
  const uint16 color0_16 = ToUInt16(color0);
  const uint16 color1_16 = ToUInt16(color1);

  uint8 color_bits[4];
  if (color0_16 == color1_16) {
    // Shortcut for the case where the colors are the same, meaning the whole
    // block is the same color (or very nearly so).
    // This is common in some vector graphics images (like map tiles).
    ComputeConstantColorBits(pixel4x4, swap_red_and_blue, always_4_color_case,
                             base_colors[0], &color0, &color1, color_bits);
  } else {
    // Otherwise, make sure the colors are ordered properly.
    if (color0_16 < color1_16) {
      std::swap(base_colors[0], base_colors[1]);
      std::swap(color0, color1);
    }

    // Compute the color bit codes.
    ComputeColorBits(pixel4x4, swap_red_and_blue, base_colors, color_bits);
  }
  return Dxt1Block(color0, color1, color_bits);
}

// Encodes a 4x4 block of RGBA pixels into a Dxt5Block, which is returned.
static Dxt5Block EncodeDxt5Block(const Pixel4x4 &pixel4x4,
                                 bool swap_red_and_blue) {
  // Find the alpha values in the block to use as base values for the
  // encoding.
  uint8 base_alphas[2];
  ComputeBaseAlphas(pixel4x4, base_alphas);

  // Compute the alpha bit codes.
  uint8 alpha_bits[6];
  ComputeAlphaBits(pixel4x4, base_alphas, alpha_bits);
  return Dxt5Block(EncodeDxt1Block(pixel4x4, swap_red_and_blue, true),
                   base_alphas, alpha_bits);
}

//-----------------------------------------------------------------------------

//
// Other helper functions.
//

// Copies colors from one Dxt1Block to another, leaving the code bits alone.
static inline void CopyDxt1Colors(const Dxt1Block &from_block,
                                  Dxt1Block *to_block) {
  to_block->c0_lo = from_block.c0_lo;
  to_block->c0_hi = from_block.c0_hi;
  to_block->c1_lo = from_block.c1_lo;
  to_block->c1_hi = from_block.c1_hi;
}

// Given the color code bits for a row of a Dxt1Block, this returns
// the code bits to use to copy the column 3 bits to the other 3
// columns.
static inline uint8 CopyColumn3ColorBits(uint8 row_bits) {
  // Get the 2 bits for column 3.
  const uint8 col3_bit_code = (row_bits >> 6) & 3;
  // Replicating those 2 bits to all 4 columns is the same as
  // multiplying them by the binary number 01010101 = 0x55.
  return col3_bit_code * 0x55;
}

//-----------------------------------------------------------------------------

//
// Functors for Compressor4x4Helper. These are all specialized for
// Dxt1Block and Dxt5Block.
//

template <typename BlockType> struct DxtcEncode {
  const BlockType operator()(const Pixel4x4 &pixel4x4, bool swap_red_and_blue);
};

template <> const Dxt1Block DxtcEncode<Dxt1Block>::operator()(
    const Pixel4x4 &pixel4x4, bool swap_red_and_blue) {
  return EncodeDxt1Block(pixel4x4, swap_red_and_blue, false);
}

template <> const Dxt5Block DxtcEncode<Dxt5Block>::operator()(
    const Pixel4x4 &pixel4x4, bool swap_red_and_blue) {
  return EncodeDxt5Block(pixel4x4, swap_red_and_blue);
}

template <typename BlockType> struct DxtcDecode {
  void operator()(const BlockType &block, bool swap_red_and_blue,
                  typename BlockType::ColorType decoded_pixels[4][4]);
};

template <> void DxtcDecode<Dxt1Block>::operator()(
    const Dxt1Block &block, bool swap_red_and_blue,
    Rgb888 decoded_pixels[4][4]) {
  DecodeDxt1Block(block, swap_red_and_blue, decoded_pixels);
}

template <> void DxtcDecode<Dxt5Block>::operator()(
    const Dxt5Block &block, bool swap_red_and_blue,
    Rgba8888 decoded_pixels[4][4]) {
  DecodeDxt5Block(block, swap_red_and_blue, decoded_pixels);
}

template <typename BlockType> struct DxtcGetColumnPadBlock {
  const BlockType operator()(const BlockType &last_column_block);
};

template <> const Dxt1Block DxtcGetColumnPadBlock<Dxt1Block>::operator()(
    const Dxt1Block &last_column_block) {
  // Do this the quick way by using the same colors and just modifying
  // the bits for the last 3 columns.
  Dxt1Block pad_block;
  CopyDxt1Colors(last_column_block, &pad_block);
  for (int row = 0; row < 4; ++row)
    pad_block.color_bits[row] =
        CopyColumn3ColorBits(last_column_block.color_bits[row]);
  return pad_block;
}

template <> const Dxt5Block DxtcGetColumnPadBlock<Dxt5Block>::operator()(
    const Dxt5Block &last_column_block) {
  // This is essentially the same as the DXT1 version, with the
  // addition of dealing with the alpha bits.
  Dxt5Block pad_block;
  pad_block.dxt1_block =
      DxtcGetColumnPadBlock<Dxt1Block>()(last_column_block.dxt1_block);
  pad_block.alpha0 = last_column_block.alpha0;
  pad_block.alpha1 = last_column_block.alpha1;
  Dxt5AlphaBits bits(last_column_block.alpha_bits);
  for (int row = 0; row < 4; ++row) {
    const int row_offset = 4 * row;
    const int col3_code = bits.GetCode(row_offset + 3);
    for (int col = 0; col < 3; ++col)
      bits.SetCode(row_offset + col, col3_code);
  }
  bits.GetBytes(pad_block.alpha_bits);
  return pad_block;
}

template <typename BlockType> struct DxtcGetRowPadBlock {
  const BlockType operator()(const BlockType &last_row_block);
};

template <> const Dxt1Block DxtcGetRowPadBlock<Dxt1Block>::operator()(
    const Dxt1Block &last_row_block) {
  // Do this the quick way by using the same colors and just modifying
  // the bits for the last 3 rows.
  Dxt1Block pad_block;
  CopyDxt1Colors(last_row_block, &pad_block);
  const uint8 last_row_bits = last_row_block.color_bits[3];
  for (int row = 0; row < 4; ++row)
    pad_block.color_bits[row] = last_row_bits;
  return pad_block;
}

template <> const Dxt5Block DxtcGetRowPadBlock<Dxt5Block>::operator()(
    const Dxt5Block &last_row_block) {
  // This is essentially the same as the DXT1 version, with the
  // addition of dealing with the alpha bits.
  Dxt5Block pad_block;
  pad_block.dxt1_block =
      DxtcGetRowPadBlock<Dxt1Block>()(last_row_block.dxt1_block);
  pad_block.alpha0 = last_row_block.alpha0;
  pad_block.alpha1 = last_row_block.alpha1;
  Dxt5AlphaBits bits(last_row_block.alpha_bits);
  for (int col = 0; col < 4; ++col) {
    const int row3_code = bits.GetCode(4 * 3 + col);
    for (int row = 0; row < 3; ++row)
      bits.SetCode(4 * row + col, row3_code);
  }
  bits.GetBytes(pad_block.alpha_bits);
  return pad_block;
}

template <typename BlockType> struct DxtcGetCornerPadBlock {
  const BlockType operator()(const BlockType &last_block);
};

template <> const Dxt1Block DxtcGetCornerPadBlock<Dxt1Block>::operator()(
    const Dxt1Block &last_block) {
  // Do this the quick way by using the same colors and just modifying
  // the bits for all 16 pixels.
  Dxt1Block pad_block;
  CopyDxt1Colors(last_block, &pad_block);
  const uint8 last_bits = CopyColumn3ColorBits(last_block.color_bits[3]);
  for (int row = 0; row < 4; ++row)
    pad_block.color_bits[row] = last_bits;
  return pad_block;
}

template <> const Dxt5Block DxtcGetCornerPadBlock<Dxt5Block>::operator()(
    const Dxt5Block &last_block) {
  // This is essentially the same as the DXT1 version, with the
  // addition of dealing with the alpha bits.
  Dxt5Block pad_block;
  pad_block.dxt1_block =
      DxtcGetCornerPadBlock<Dxt1Block>()(last_block.dxt1_block);
  pad_block.alpha0 = last_block.alpha0;
  pad_block.alpha1 = last_block.alpha1;
  Dxt5AlphaBits bits(last_block.alpha_bits);
  uint8 corner_bits = bits.GetCode(15);
  for (int i = 0; i < 16; ++i)
    bits.SetCode(i, corner_bits);
  bits.GetBytes(pad_block.alpha_bits);
  return pad_block;
}

//-----------------------------------------------------------------------------

//
// Public functions.
//

DxtcCompressor::~DxtcCompressor() {
}

bool DxtcCompressor::SupportsFormat(CompressedImage::Format format) const {
  // DXTC supports all current formats.
  return true;
}

bool DxtcCompressor::IsValidCompressedImage(const CompressedImage &image) {
  const CompressedImage::Metadata &metadata = image.GetMetadata();
  return
      metadata.compressor_name == "dxtc" &&
      metadata.uncompressed_height > 0 &&
      metadata.uncompressed_width > 0 &&
      metadata.compressed_height >= metadata.uncompressed_height &&
      metadata.compressed_width >= metadata.uncompressed_width &&
      image.GetDataSize() ==
      ComputeCompressedDataSize(metadata.format, metadata.compressed_height,
                                metadata.compressed_width);
}

size_t DxtcCompressor::ComputeCompressedDataSize(CompressedImage::Format format,
                                                 uint32 height, uint32 width) {
  if (height == 0 || width == 0)
    return 0;
  const uint32 num_block_rows = Compressor4x4Helper::GetNumBlocks(height);
  const uint32 num_block_cols = Compressor4x4Helper::GetNumBlocks(width);
  return std::max(1U, num_block_rows) * std::max(1U, num_block_cols) *
      GetBlockSize(format);
}

bool DxtcCompressor::Compress(CompressedImage::Format format,
                              uint32 height, uint32 width,
                              uint32 padding_bytes_per_row,
                              const uint8 *buffer, CompressedImage *image) {
  if (!buffer || !image || height == 0 || width == 0) {
    return false;
  } else if (GetNumFormatComponents(format) == 3) {
    return Compressor4x4Helper::Compress<Dxt1Block, Rgb888>(
        DxtcEncode<Dxt1Block>(), "dxtc", format, height, width,
        padding_bytes_per_row, buffer, image);
  } else {
    return Compressor4x4Helper::Compress<Dxt5Block, Rgba8888>(
        DxtcEncode<Dxt5Block>(), "dxtc", format, height, width,
        padding_bytes_per_row, buffer, image);
  }
}

bool DxtcCompressor::Decompress(const CompressedImage &image,
                                std::vector<uint8> *decompressed_buffer) {
  if (!IsValidCompressedImage(image) || !decompressed_buffer) {
    return false;
  } else if (GetNumFormatComponents(image.GetMetadata().format) == 3) {
    return Compressor4x4Helper::Decompress<Dxt1Block, Rgb888>(
        DxtcDecode<Dxt1Block>(), image, decompressed_buffer);
  } else {
    return Compressor4x4Helper::Decompress<Dxt5Block, Rgba8888>(
        DxtcDecode<Dxt5Block>(), image, decompressed_buffer);
  }
}

bool DxtcCompressor::Downsample(const CompressedImage &image,
                                CompressedImage *downsampled_image) {
  if (!IsValidCompressedImage(image) || !downsampled_image) {
    return false;
  } else if (GetNumFormatComponents(image.GetMetadata().format) == 3) {
    return Compressor4x4Helper::Downsample<Dxt1Block, Rgb888>(
        DxtcEncode<Dxt1Block>(), DxtcDecode<Dxt1Block>(),
        image, downsampled_image);
  } else {
    return Compressor4x4Helper::Downsample<Dxt5Block, Rgba8888>(
        DxtcEncode<Dxt5Block>(), DxtcDecode<Dxt5Block>(),
        image, downsampled_image);
  }
}

bool DxtcCompressor::Pad(const CompressedImage &image, uint32 padded_height,
                         uint32 padded_width, CompressedImage *padded_image) {
  if (!IsValidCompressedImage(image) || !padded_image) {
    return false;
  } else if (GetNumFormatComponents(image.GetMetadata().format) == 3) {
    return Compressor4x4Helper::Pad<Dxt1Block>(
        DxtcGetColumnPadBlock<Dxt1Block>(),
        DxtcGetRowPadBlock<Dxt1Block>(),
        DxtcGetCornerPadBlock<Dxt1Block>(),
        image, padded_height, padded_width, padded_image);
  } else {
    return Compressor4x4Helper::Pad<Dxt5Block>(
        DxtcGetColumnPadBlock<Dxt5Block>(),
        DxtcGetRowPadBlock<Dxt5Block>(),
        DxtcGetCornerPadBlock<Dxt5Block>(),
        image, padded_height, padded_width, padded_image);
  }
}

bool DxtcCompressor::CompressAndPad(CompressedImage::Format format,
                                   uint32 height, uint32 width,
                                   uint32 padded_height, uint32 padded_width,
                                   uint32 padding_bytes_per_row,
                                   const uint8 *buffer,
                                   CompressedImage *padded_image) {
  if (!buffer || !padded_image || height == 0 || width == 0) {
    return false;
  } else if (GetNumFormatComponents(format) == 3) {
    return Compressor4x4Helper::CompressAndPad<Dxt1Block, Rgb888>(
        DxtcEncode<Dxt1Block>(), "dxtc", format, height, width,
        padded_height, padded_width, padding_bytes_per_row,
        buffer, padded_image);
  } else {
    return Compressor4x4Helper::CompressAndPad<Dxt5Block, Rgba8888>(
        DxtcEncode<Dxt5Block>(), "dxtc", format, height, width,
        padded_height, padded_width, padding_bytes_per_row,
        buffer, padded_image);
  }
}

bool DxtcCompressor::CreateSolidImage(CompressedImage::Format format,
                                      uint32 height, uint32 width,
                                      const uint8 *color,
                                      CompressedImage *image) {
  if (!image) {
    return false;
  } else if (GetNumFormatComponents(format) == 3) {
    return Compressor4x4Helper::CreateSolidImage<Dxt1Block>(
        "dxtc", format, height, width,
        Dxt1Block(QuantizeToRgb565(RgbInt(color[0], color[1], color[2]))),
        image);
  } else {
    return Compressor4x4Helper::CreateSolidImage<Dxt5Block>(
        "dxtc", format, height, width,
        Dxt5Block(
            Dxt1Block(QuantizeToRgb565(RgbInt(color[0], color[1], color[2]))),
            color[3]),
        image);
  }
}

bool DxtcCompressor::CopySubimage(const CompressedImage &image,
                                  uint32 start_row, uint32 start_column,
                                  uint32 height, uint32 width,
                                  CompressedImage *subimage) {
  if (!IsValidCompressedImage(image) || !subimage) {
    return false;
  } else if (GetNumFormatComponents(image.GetMetadata().format) == 3) {
    return Compressor4x4Helper::CopySubimage<Dxt1Block>(
        image, start_row, start_column, height, width, subimage);
  } else {
    return Compressor4x4Helper::CopySubimage<Dxt5Block>(
        image, start_row, start_column, height, width, subimage);
  }
}

}  // namespace image_codec_compression
