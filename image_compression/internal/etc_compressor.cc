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

#include "image_compression/public/etc_compressor.h"

#include <stddef.h>
#include <algorithm>
#include <vector>

#include "base/integral_types.h"
#include "base/logging.h"
#include "image_compression/internal/bit_util.h"
#include "image_compression/internal/color_util.h"
#include "image_compression/internal/compressor4x4_helper.h"
#include "image_compression/internal/pixel4x4.h"
#include "image_compression/public/compressed_image.h"

#define FRIEND_TEST(test_case_name, test_name) \
  friend class test_case_name##_##test_name##_Test

namespace image_codec_compression {

//============================================================================
//
// General information on ETC1 compression and decompression can be
// found on the web, if you're lucky.  This describes the version
// implemented here:
//    http://www.khronos.org/registry/gles/extensions/OES/
//              OES_compressed_ETC1_RGB8_texture.txt
// A 64-bit ETC1 block represents a 4x4 block of pixels, as follows:
//
// High-order 32 bits:
//       0:  "flip bit" - indicates orientation of 2x4 subblocks.
//       1:  "diff bit" - indicates color encoding method.
//     2-4:  codebook table codeword 2 for the second subblock (3 bits)
//     5-7:  codebook table codeword 1 for the first  subblock (3 bits)
//  If diffbit is 0:
//    8-11:  Base color 2 blue  (4 bits)
//   12-15:  Base color 1 blue  (4 bits)
//   16-19:  Base color 2 green (4 bits)
//   20-23:  Base color 1 green (4 bits)
//   24-27:  Base color 2 red   (4 bits)
//   28-31:  Base color 1 red   (4 bits)
//  If diffbit is 1:
//    8-10:  Color difference blue  (3 bits)
//   11-15:  Base color blue        (5 bits)
//   16-18:  Color difference green (3 bits)
//   19-23:  Base color green       (5 bits)
//   24-26:  Color difference red   (3 bits)
//   27-31:  Base color red         (5 bits)
//
// Low-order 32 bits:
//   Index for each of the 16 pixels into the codebook: 2 bits per
//   pixel. The low-order bit of the index is in bits 0-15 and the
//   high-order bit is in bits 16-31. The pixel order is column-major,
//   so that the placement of bit at pixel (x,y) within the 4x4 is at
//   (x * 4 + y).
//
//============================================================================

//
// ETC-specific data structures.
//

// This enum and struct is used to define pixel indices for a subblock
// of a 4x4 block of pixels.  The subblock is the left or right 2x4
// when the flip bit is false or the top or bottom 4x2 when the flip
// bit is true.
struct Etc1Subblock {
  // Sets the inclusive indices of the block pixels. For example, the
  // bottom subblock of a flipped block has x_min=0, x_max=3, y_min=2,
  // y_max=3.
  Etc1Subblock(int x_min_in, int x_max_in, int y_min_in, int y_max_in)
      : x_min(x_min_in),
        x_max(x_max_in),
        y_min(y_min_in),
        y_max(y_max_in) {
  }
  int x_min, x_max;
  int y_min, y_max;
};

// This class is used to help with ETC encoding and decoding. It
// encapsulates the codebook and pixel operations.
class EtcHelper {
 public:
  // Returns the ETC codebook value for the given codeword and pixel
  // index.
  static int GetCodebookModifier(int codeword, int pixel_index) {
    static const int kCodeBook[8][4] = {
      {  2,   8,  -2,   -8 },
      {  5,  17,  -5,  -17 },
      {  9,  29,  -9,  -29 },
      { 13,  42, -13,  -42 },
      { 18,  60, -18,  -60 },
      { 24,  80, -24,  -80 },
      { 33, 106, -33, -106 },
      { 47, 183, -47, -183 },
    };
    DCHECK_GE(codeword, 0);
    DCHECK_LE(codeword, 8);
    DCHECK_GE(pixel_index, 0);
    DCHECK_LT(pixel_index, 4);
    return kCodeBook[codeword][pixel_index];
  }

  // Returns the result of adding the ETC codebook value for the given
  // codeword and pixel index to the given color. This clamps the
  // resulting color to 0-255 for all components.
  static const RgbInt AddCodebookModifier(const RgbInt color,
                                          int codeword, int pixel_index) {
    const int modifier = GetCodebookModifier(codeword, pixel_index);
    return ClampRgbInt(color + RgbInt(modifier, modifier, modifier));
  }

  // Given the y and x of a pixel in a 4x4block (both x and y are 0-3,
  // y increases from top to bottom, and x increases from left to
  // right), this returns an index (from 0-16) into the ETC pixel
  // indices.
  static int GetPixelOrderIndex(int y, int x) {
    DCHECK_GE(x, 0);
    DCHECK_GE(y, 0);
    DCHECK_LT(x, 4);
    DCHECK_LT(y, 4);
    return x * 4 + y;
  }

  // Accesses the 2-bit code index for a pixel from the given 32-bit
  // encoded pixel indices (low-order word of the ETC1 encoding).  The
  // x and y of the pixel are passed in.
  static int GetPixelIndex(int y, int x, uint32 pixel_indices) {
    const int p = GetPixelOrderIndex(y, x);
    return GetBits(pixel_indices, p, 1) |
        (GetBits(pixel_indices, p + 16, 1) << 1);
  }

  // Stores the 2-bit code index for a pixel in the given 32-bit
  // encoding word. The x and y of the pixel are passed in.
  static void StorePixelIndex(int y, int x, int code, uint32 *pixel_indices) {
    DCHECK_GE(code, 0);
    DCHECK_LT(code, 4);
    const int p = GetPixelOrderIndex(y, x);
    SetBits(p, 1, code & 1, pixel_indices);
    SetBits(p + 16, 1, (code & 2) >> 1, pixel_indices);
  }

  // Splits a 64-bit block into two 32-bit words, swizzling bytes as necessary
  // to get the correct endian ordering.
  static void SplitBlock(uint64 block, uint32 *hi_word, uint32 *lo_word) {
#if defined IS_LITTLE_ENDIAN
    *hi_word = Swizzle32(block & 0xffffffff);
    *lo_word = Swizzle32(block >> 32);
#else
    *hi_word = block >> 32;
    *lo_word = block & 0xffffffff;
#endif
  }

  // Glues two 32-bit words together to form a 64-bit block, swizzling bytes as
  // necessary to get the correct byte order expected by hardware.
  static uint64 BuildBlock(uint32 hi_word, uint32 lo_word) {
#if defined IS_LITTLE_ENDIAN
    const uint32 swizzled_hi_word = Swizzle32(hi_word);
    const uint32 swizzled_lo_word = Swizzle32(lo_word);
    return (static_cast<uint64>(swizzled_lo_word) << 32) | swizzled_hi_word;
#else
    return (static_cast<uint64>(hi_word) << 32) | lo_word;
#endif
  }

  // Swizzles a 32-bit word in ETC1 order to match the host endianness.
  static uint32 Swizzle32(uint32 n) {
#if defined IS_LITTLE_ENDIAN
    const uint8 *bytes = reinterpret_cast<const uint8*>(&n);
    return
        static_cast<uint32>(bytes[0]) << 24 |
        static_cast<uint32>(bytes[1]) << 16 |
        static_cast<uint32>(bytes[2]) << 8 |
        static_cast<uint32>(bytes[3]);
#else
    return n;
#endif
  }
};

// This class is used to decode a 64-bit ETC1 block.
class Etc1BlockDecoder {
 public:
  // The constructor is passed a 64-bit block to decode.
  explicit Etc1BlockDecoder(uint64 bits);

  // Returns the color for the given pixel.
  const Rgb888 GetPixel(int y, int x) const {
    const int pixel_index = pixel_indices_[EtcHelper::GetPixelOrderIndex(y, x)];
    const bool is_first_block = flip_bit_ ? (y < 2) : (x < 2);
    const int codeword_index = is_first_block ? 0 : 1;

    // Add the codebook modifier to the correct base color.
    const RgbInt color = is_first_block ? base_color1_ : base_color2_;
    return ToRgb888(EtcHelper::AddCodebookModifier(
        color, codewords_[codeword_index], pixel_index));
  }

 private:
  bool diff_bit_;
  bool flip_bit_;
  int codewords_[2];
  int pixel_indices_[16];
  RgbInt base_color1_;
  RgbInt base_color2_;

  FRIEND_TEST(EtcCompressorTest, DecoderBits);
  FRIEND_TEST(EtcCompressorTest, ColorBits);
};

Etc1BlockDecoder::Etc1BlockDecoder(uint64 bits) {
  // Split the block into two 32-bit pieces for clarity.
  uint32 hi_word, lo_word;
  EtcHelper::SplitBlock(bits, &hi_word, &lo_word);

  flip_bit_ = (hi_word & 1) != 0;
  diff_bit_ = (hi_word & 2) != 0;

  codewords_[0] = GetBits(hi_word, 5, 3);
  codewords_[1] = GetBits(hi_word, 2, 3);

  if (diff_bit_) {
    // The first base color is extracted directly from 5-bit RGB.
    RgbInt base_color(GetBits(hi_word, 27, 5),
                      GetBits(hi_word, 19, 5),
                      GetBits(hi_word, 11, 5));
    base_color1_ = RgbInt(Extend5Bit(base_color.r),
                          Extend5Bit(base_color.g),
                          Extend5Bit(base_color.b));

    // The second base color is a color difference (signed 3-bit RGB)
    // added to the first base color (before extending to 8 bits).
    RgbInt diff(ExtendSignBit(GetBits(hi_word, 24, 3), 3),
                ExtendSignBit(GetBits(hi_word, 16, 3), 3),
                ExtendSignBit(GetBits(hi_word,  8, 3), 3));
    base_color += diff;
    base_color2_ = RgbInt(Extend5Bit(base_color.r),
                          Extend5Bit(base_color.g),
                          Extend5Bit(base_color.b));
  } else {
    // In individual coding mode, each of the two base colors is
    // extracted directly from 4-bit RGB.
    base_color1_ = RgbInt(Extend4Bit(GetBits(hi_word, 28, 4)),
                          Extend4Bit(GetBits(hi_word, 20, 4)),
                          Extend4Bit(GetBits(hi_word, 12, 4)));
    base_color2_ = RgbInt(Extend4Bit(GetBits(hi_word, 24, 4)),
                          Extend4Bit(GetBits(hi_word, 16, 4)),
                          Extend4Bit(GetBits(hi_word,  8, 4)));
  }

  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < 4; ++x) {
      const int p = EtcHelper::GetPixelOrderIndex(y, x);
      pixel_indices_[p] = EtcHelper::GetPixelIndex(y, x, lo_word);
    }
  }
}

//-----------------------------------------------------------------------------

//
// General ETC decompression helper functions.
//

// Decodes a single 64-bit ETC1 block into a 4x4 array of RGB pixels.
static void DecodeBlock(uint64 block, Rgb888 decoded_pixels[4][4]) {
  Etc1BlockDecoder decoder(block);
  for (int y = 0; y < 4; ++y) {
    for (int x = 0; x < 4; ++x) {
      decoded_pixels[y][x] = decoder.GetPixel(y, x);
    }
  }
}

//-----------------------------------------------------------------------------

//
// General ETC compression helper functions.
//

// Computes and returns the average color of a 2x4 subblock of a 4x4
// array of RGB pixels.
static const RgbInt ComputeAverageColor(const Pixel4x4 &pixel4x4,
                                        const Etc1Subblock &subblock) {
  // There should be 8 pixels in the subblock.
  DCHECK_EQ((subblock.x_max - subblock.x_min + 1) *
            (subblock.y_max - subblock.y_min + 1), 8);

  RgbInt sum(0, 0, 0);
  for (int y = subblock.y_min; y <= subblock.y_max; ++y) {
    for (int x = subblock.x_min; x <= subblock.x_max; ++x) {
      sum += pixel4x4.GetPixel(y, x);
    }
  }
  return sum / 8;
}

// Stores two base colors in normal (non-differential) mode into the
// given high-order block word.
static void StoreNormalModeColors(const RgbInt color1, const RgbInt color2,
                                  uint32 *hi_word) {
  SetBits(28, 4, color1.r, hi_word);
  SetBits(20, 4, color1.g, hi_word);
  SetBits(12, 4, color1.b, hi_word);
  SetBits(24, 4, color2.r, hi_word);
  SetBits(16, 4, color2.g, hi_word);
  SetBits(8,  4, color2.b, hi_word);
}

// Stores the base color and difference in differential mode into the
// given high-order block word.
static void StoreDiffModeColors(const RgbInt color1,
                                const RgbInt color_difference,
                                uint32 *hi_word) {
  SetBits(27, 5, color1.r, hi_word);
  SetBits(19, 5, color1.g, hi_word);
  SetBits(11, 5, color1.b, hi_word);
  SetBits(24, 3, color_difference.r, hi_word);
  SetBits(16, 3, color_difference.g, hi_word);
  SetBits(8,  3, color_difference.b, hi_word);
}

// Computes the approximate error between two color values.  This uses
// a very simple sum of squares metric (with no perceptual weighting),
// which is good enough for most images and is fast.
static uint32 ComputeColorError(const RgbInt c0, const RgbInt c1) {
  const RgbInt diff = c1 - c0;
  return diff.r * diff.r + diff.g * diff.g + diff.b * diff.b;
}

// Returns the error produced when encoding a subblock with the given
// table codeword and color. Returns the indices to use for each
// subblock pixel in pixel_indices.
static uint32 ComputeCodewordError(const Pixel4x4 &pixel4x4,
                                   const Etc1Subblock &subblock,
                                   int codeword, const RgbInt color,
                                   uint32 *pixel_indices) {
  uint32 cumulative_error = 0;
  *pixel_indices = 0;
  const RgbInt c0 = EtcHelper::AddCodebookModifier(color, codeword, 0);
  const RgbInt c1 = EtcHelper::AddCodebookModifier(color, codeword, 1);
  const RgbInt c2 = EtcHelper::AddCodebookModifier(color, codeword, 2);
  const RgbInt c3 = EtcHelper::AddCodebookModifier(color, codeword, 3);
  for (int y = subblock.y_min; y <= subblock.y_max; ++y) {
    for (int x = subblock.x_min; x <= subblock.x_max; ++x) {
      const RgbInt target_color = pixel4x4.GetPixel(y, x);
      int best_index = 0;
      uint32 smallest_error = ComputeColorError(target_color, c0);
      uint32 new_error = ComputeColorError(target_color, c1);
      if (new_error < smallest_error) {
        smallest_error = new_error;
        best_index = 1;
      }
      new_error = ComputeColorError(target_color, c2);
      if (new_error < smallest_error) {
        smallest_error = new_error;
        best_index = 2;
      }
      new_error = ComputeColorError(target_color, c3);
      if (new_error < smallest_error) {
        smallest_error = new_error;
        best_index = 3;
      }
      EtcHelper::StorePixelIndex(y, x, best_index, pixel_indices);
      cumulative_error += smallest_error;
    }
  }
  return cumulative_error;
}

// Returns the table codeword that produces the smallest error when
// encoding a subblock with the given color. Returns the indices to
// use for each subblock pixel in pixel_indices and returns the
// cumulative error for the entire subblock in cumulative_error.
static int FindBestCodeword(const Pixel4x4 &pixel4x4,
                            const Etc1Subblock &subblock, const RgbInt color,
                            uint32 *pixel_indices, uint32 *cumulative_error) {
  // Try each codeword and choose the one with the smallest error.
  int best_codeword = -1;
  *cumulative_error = kuint32max;
  for (int i = 0; i < 8; ++i) {
    uint32 indices;
    const uint32 error = ComputeCodewordError(pixel4x4, subblock, i,
                                              color, &indices);
    if (error < *cumulative_error) {
      best_codeword = i;
      *pixel_indices = indices;
      *cumulative_error = error;
    }
  }
  DCHECK_GE(best_codeword, 0);
  return best_codeword;
}

// Returns the table codeword as chosen by a heuristic for the subblock
// and given color.  Returns the indices to use for each subblock pixel in
// pixel_indices and returns the cumulative error for the entire subblock
// in cumulative_error.
static int FindCodewordHeuristic(
    const Pixel4x4 &pixel4x4, const Etc1Subblock &subblock,
    const RgbInt color, uint32 *pixel_indices, uint32 *cumulative_error) {
  // Compute absolute deviation and pick a codeword directly.
  // (Absolute deviation is faster to compute than standard deviation and
  // does not perform significantly worse.)
  RgbInt absdevsum(0, 0, 0);
  for (int y = subblock.y_min; y <= subblock.y_max; ++y) {
    for (int x = subblock.x_min; x <= subblock.x_max; ++x) {
      const RgbInt target_color = pixel4x4.GetPixel(y, x);
      const RgbInt diff = (color - target_color);
      absdevsum.r += abs(diff.r);
      absdevsum.g += abs(diff.g);
      absdevsum.b += abs(diff.b);
    }
  }
  absdevsum = absdevsum / 8;
  // Take largest component, pick the codeword of comparable magnitude.
  const int dev = std::max(absdevsum.r, std::max(absdevsum.g, absdevsum.b));
  int codeword;
  if (dev > 144) {
    codeword = 7;
  } else if (dev > 93) {
    codeword = 6;
  } else if (dev > 70) {
    codeword = 5;
  } else if (dev > 51) {
    codeword = 4;
  } else if (dev > 35) {
    codeword = 3;
  } else if (dev > 23) {
    codeword = 2;
  } else if (dev > 12) {
    codeword = 1;
  } else {
    codeword = 0;
  }
  *cumulative_error = ComputeCodewordError(pixel4x4, subblock, codeword,
                                           color, pixel_indices);
  return codeword;
}

// Returns the 64-bit ETC1 block representing the best encoding for
// the given 4x4 block with the given flip-bit setting. This returns
// the cumulative error computed for the block in cumulative_error.
static uint64 FindBestSubblockEncoding(
    const Pixel4x4 &pixel4x4, bool flip,
    EtcCompressor::CompressionStrategy strategy,
    uint32 *cumulative_error) {
  const Etc1Subblock subblock1 = flip ? Etc1Subblock(0, 3, 0, 1) :
                                 Etc1Subblock(0, 1, 0, 3);
  const Etc1Subblock subblock2 = flip ? Etc1Subblock(0, 3, 2, 3) :
                                 Etc1Subblock(2, 3, 0, 3);

  // Compute the average color of the two subblocks.
  const RgbInt avg_pixel1 = ComputeAverageColor(pixel4x4, subblock1);
  const RgbInt avg_pixel2 = ComputeAverageColor(pixel4x4, subblock2);

  // Quantize to 5 bit components.
  const RgbInt avg_pixel1_555 = QuantizeRgbFast<5>(avg_pixel1);
  const RgbInt avg_pixel2_555 = QuantizeRgbFast<5>(avg_pixel2);

  // Choose differential or regular mode. Differential mode is used
  // iff each component of the difference between the quantized
  // average colors is in the range [-4,3].
  const RgbInt difference = avg_pixel2_555 - avg_pixel1_555;
  const bool use_diff_mode = (difference.r >= -4 && difference.r <= 3 &&
                              difference.g >= -4 && difference.g <= 3 &&
                              difference.b >= -4 && difference.b <= 3);

  // Store results in the two 32-bit words of the block.
  uint32 hi_word = 0;
  uint32 lo_word = 0;

  // Set the flip bit.
  SetBits(0, 1, flip, &hi_word);

  // Colors that will be used when decoding this block.  These are
  // needed for computing the best codewords for the subblocks.
  RgbInt decoded_color1, decoded_color2;

  if (use_diff_mode) {
    SetBits(1, 1, 1, &hi_word);  // Set the diff bit.
    StoreDiffModeColors(avg_pixel1_555, difference, &hi_word);
    decoded_color1 = RgbInt(Extend5Bit(avg_pixel1_555.r),
                            Extend5Bit(avg_pixel1_555.g),
                            Extend5Bit(avg_pixel1_555.b));
    decoded_color2 = RgbInt(Extend5Bit(avg_pixel2_555.r),
                            Extend5Bit(avg_pixel2_555.g),
                            Extend5Bit(avg_pixel2_555.b));
  } else {
    // Quantize the average colors to RGB444.
    const RgbInt avg_pixel1_444 = QuantizeRgbFast<4>(avg_pixel1);
    const RgbInt avg_pixel2_444 = QuantizeRgbFast<4>(avg_pixel2);
    StoreNormalModeColors(avg_pixel1_444, avg_pixel2_444, &hi_word);
    decoded_color1 = RgbInt(Extend4Bit(avg_pixel1_444.r),
                            Extend4Bit(avg_pixel1_444.g),
                            Extend4Bit(avg_pixel1_444.b));
    decoded_color2 = RgbInt(Extend4Bit(avg_pixel2_444.r),
                            Extend4Bit(avg_pixel2_444.g),
                            Extend4Bit(avg_pixel2_444.b));
  }

  // Find the codewords that produce the smallest error.
  uint32 pixel_indices1 = 0, pixel_indices2 = 0;
  uint32 error1, error2;
  int codeword1;
  int codeword2;
  if (strategy == EtcCompressor::kHeuristic) {
    codeword1 = FindCodewordHeuristic(pixel4x4, subblock1, decoded_color1,
                                      &pixel_indices1, &error1);
    codeword2 = FindCodewordHeuristic(pixel4x4, subblock2, decoded_color2,
                                      &pixel_indices2, &error2);
  } else {
    codeword1 = FindBestCodeword(pixel4x4, subblock1, decoded_color1,
                                 &pixel_indices1, &error1);
    codeword2 = FindBestCodeword(pixel4x4, subblock2, decoded_color2,
                                 &pixel_indices2, &error2);
  }

  SetBits(5, 3, codeword1, &hi_word);
  SetBits(2, 3, codeword2, &hi_word);

  lo_word = pixel_indices1 | pixel_indices2;

  *cumulative_error = error1 + error2;
  return EtcHelper::BuildBlock(hi_word, lo_word);
}

// Encodes a 4x4 array of RGB pixels into a single 64-bit ETC1 block.
uint64 EncodeEtc1Block(const Pixel4x4 &pixel4x4,
                       EtcCompressor::CompressionStrategy strategy) {
  uint32 error_lr, error_tb;
  switch (strategy) {
    case EtcCompressor::kSplitHorizontally:
      return FindBestSubblockEncoding(pixel4x4, true, strategy, &error_tb);
    case EtcCompressor::kSplitVertically:
      return FindBestSubblockEncoding(pixel4x4, false, strategy, &error_lr);
    case EtcCompressor::kHeuristic: {
      // Find the average difference in color between the top-half and the
      // bottom-half and between the left-half and the right-half.  Choose to
      // split in the direction that has the larger difference.
      RgbInt sum1 = pixel4x4.GetPixel(0, 0) + pixel4x4.GetPixel(0, 1) +
          pixel4x4.GetPixel(1, 0) + pixel4x4.GetPixel(1, 1);
      RgbInt sum2 = pixel4x4.GetPixel(2, 0) + pixel4x4.GetPixel(2, 1) +
          pixel4x4.GetPixel(3, 0) + pixel4x4.GetPixel(3, 1);
      RgbInt sum3 = pixel4x4.GetPixel(0, 2) + pixel4x4.GetPixel(0, 3) +
          pixel4x4.GetPixel(1, 2) + pixel4x4.GetPixel(1, 3);
      RgbInt sum4 = pixel4x4.GetPixel(2, 2) + pixel4x4.GetPixel(2, 3) +
          pixel4x4.GetPixel(3, 2) + pixel4x4.GetPixel(2, 2);
      RgbInt left = (sum1 + sum2) / 8;
      RgbInt right = (sum3 + sum4) / 8;
      RgbInt top = (sum1 + sum3) / 8;
      RgbInt bottom = (sum2 + sum4) / 8;
      if (ComputeColorError(left, right) > ComputeColorError(top, bottom)) {
        return FindBestSubblockEncoding(pixel4x4, false, strategy, &error_lr);
      } else {
        return FindBestSubblockEncoding(pixel4x4, true, strategy, &error_tb);
      }
    }
    case EtcCompressor::kSmallerError:
    default: {
      // Try unflipped (left and right 4x2 subblocks) and flipped (top and
      // bottom 2x4 subblocks) modes.
      const uint64 block_lr = FindBestSubblockEncoding(pixel4x4, false,
                                                       strategy, &error_lr);
      const uint64 block_tb = FindBestSubblockEncoding(pixel4x4, true,
                                                       strategy, &error_tb);
      return error_lr <= error_tb ? block_lr : block_tb;
    }
  }
}

//-----------------------------------------------------------------------------

//
// Other helper functions.
//

// Creates an ETC1-encoded solid-color block.
static const uint64 CreateSolidBlock(const Rgb888 &color) {
  uint32 hi_word = 0;

  // The smallest entry in the codebook is 2 (row 0, column 0), so
  // subtract 2 if possible from the color before quantizing.
  static const int kSmallestCodebookEntry = 2;
  const RgbInt adjusted_color(std::max(0, color.r - kSmallestCodebookEntry),
                              std::max(0, color.g - kSmallestCodebookEntry),
                              std::max(0, color.b - kSmallestCodebookEntry));


  // Use differential mode, which has 5 bits per color component.
  SetBits(1, 1, 1, &hi_word);  // Set the diff bit.
  const RgbInt quantized_color = QuantizeRgbFast<5>(ToRgbInt(color));
  StoreDiffModeColors(quantized_color, RgbInt(0, 0, 0), &hi_word);

  // Use codeword 0.
  SetBits(5, 3, 0, &hi_word);
  SetBits(2, 3, 0, &hi_word);

  // All pixel indices are 0, so lo_word = 0.
  return EtcHelper::BuildBlock(hi_word, 0);
}

//-----------------------------------------------------------------------------

//
// Functors for Compressor4x4Helper.
//

class EtcEncode {
 public:
  explicit EtcEncode(EtcCompressor::CompressionStrategy strategy)
      : strategy_(strategy) {}
  uint64 operator()(const Pixel4x4 &pixel4x4, bool swap_red_and_blue) {
    DCHECK(!swap_red_and_blue);
    return EncodeEtc1Block(pixel4x4, strategy_);
  }
 private:
  EtcCompressor::CompressionStrategy strategy_;
};

struct EtcDecode {
  void operator()(uint64 block, bool swap_red_and_blue,
                  Rgb888 decoded_pixels[4][4]) {
    DCHECK(!swap_red_and_blue);
    DecodeBlock(block, decoded_pixels);
  }
};

class EtcGetColumnPadBlock {
 public:
  explicit EtcGetColumnPadBlock(EtcCompressor::CompressionStrategy strategy)
      : strategy_(strategy) {}
  const uint64 operator()(const uint64 &last_column_block) {
    // TODO(user): Optimize this by operating directly on encoded blocks?

    // Decode the block to get colors from it.
    const Etc1BlockDecoder decoder(last_column_block);

    // Replicate the last column in a 4x4.
    Pixel4x4 pixel4x4;
    for (int y = 0; y < 4; ++y) {
      const Rgb888 last_column_color = decoder.GetPixel(y, 3);
      for (int x = 0; x < 4; ++x) {
        pixel4x4.SetPixel(y, x, last_column_color);
      }
    }
    return EncodeEtc1Block(pixel4x4, strategy_);
  }
 private:
  EtcCompressor::CompressionStrategy strategy_;
};

struct EtcGetRowPadBlock {
 public:
  explicit EtcGetRowPadBlock(EtcCompressor::CompressionStrategy strategy)
      : strategy_(strategy) {}
  const uint64 operator()(const uint64 &last_row_block) {
    // TODO(user): Optimize this by operating directly on encoded blocks?

    // Decode the block to get colors from it.
    const Etc1BlockDecoder decoder(last_row_block);

    // Replicate the last row in a 4x4.
    Pixel4x4 pixel4x4;
    for (int x = 0; x < 4; ++x) {
      const Rgb888 last_row_color = decoder.GetPixel(3, x);
      for (int y = 0; y < 4; ++y) {
        pixel4x4.SetPixel(y, x, last_row_color);
      }
    }
    return EncodeEtc1Block(pixel4x4, strategy_);
  }
 private:
  EtcCompressor::CompressionStrategy strategy_;
};

struct EtcGetCornerPadBlock {
  const uint64 operator()(const uint64 &last_block) {
    // Create a solid block using the corner color.
    return CreateSolidBlock(Etc1BlockDecoder(last_block).GetPixel(3, 3));
  }
};

//-----------------------------------------------------------------------------

//
// Public functions.
//

EtcCompressor::EtcCompressor()
    : compression_strategy_(kSmallerError) {
}

EtcCompressor::~EtcCompressor() {
}

bool EtcCompressor::SupportsFormat(CompressedImage::Format format) const {
  // ETC does not support images with alpha or color swapping.
  // TODO(user): Implement color swapping if necessary.
  return format == CompressedImage::kRGB;
}

bool EtcCompressor::IsValidCompressedImage(const CompressedImage &image) {
  const CompressedImage::Metadata &metadata = image.GetMetadata();
  return
      metadata.format == CompressedImage::kRGB &&
      metadata.compressor_name == "etc" &&
      metadata.uncompressed_height > 0 &&
      metadata.uncompressed_width > 0 &&
      metadata.compressed_height >= metadata.uncompressed_height &&
      metadata.compressed_width >= metadata.uncompressed_width &&
      image.GetDataSize() ==
      (Compressor4x4Helper::GetNumBlocks(metadata.compressed_height) *
       Compressor4x4Helper::GetNumBlocks(metadata.compressed_width) *
       sizeof(uint64));
}

size_t EtcCompressor::ComputeCompressedDataSize(CompressedImage::Format format,
                                                uint32 height, uint32 width) {
  if (height == 0 || width == 0)
    return 0;
  if (format == CompressedImage::kRGB) {
    const uint32 num_block_rows = Compressor4x4Helper::GetNumBlocks(height);
    const uint32 num_block_cols = Compressor4x4Helper::GetNumBlocks(width);
    return std::max(1U, num_block_rows) * std::max(1U, num_block_cols) *
           sizeof(uint64);
  }
  return 0;
}

bool EtcCompressor::Compress(CompressedImage::Format format,
                             uint32 height, uint32 width,
                             uint32 padding_bytes_per_row,
                             const uint8 *buffer, CompressedImage *image) {
  if (!buffer || !image || height == 0 || width == 0 ||
      format != CompressedImage::kRGB) {
    return false;
  }
  return Compressor4x4Helper::Compress<uint64, Rgb888>(
      EtcEncode(compression_strategy_), "etc", format, height, width,
      padding_bytes_per_row, buffer, image);
}

bool EtcCompressor::Decompress(const CompressedImage &image,
                               std::vector<uint8> *decompressed_buffer) {
  if (!IsValidCompressedImage(image) || !decompressed_buffer)
    return false;
  return Compressor4x4Helper::Decompress<uint64, Rgb888>(
      EtcDecode(), image, decompressed_buffer);
}

bool EtcCompressor::Downsample(const CompressedImage &image,
                               CompressedImage *downsampled_image) {
  if (!IsValidCompressedImage(image) || !downsampled_image)
    return false;

  return Compressor4x4Helper::Downsample<uint64, Rgb888>(
      EtcEncode(compression_strategy_), EtcDecode(), image, downsampled_image);
}

bool EtcCompressor::Pad(const CompressedImage &image, uint32 padded_height,
                        uint32 padded_width, CompressedImage *padded_image) {
  if (!IsValidCompressedImage(image) || !padded_image)
    return false;
  return Compressor4x4Helper::Pad<uint64>(
      EtcGetColumnPadBlock(compression_strategy_),
      EtcGetRowPadBlock(compression_strategy_), EtcGetCornerPadBlock(),
      image, padded_height, padded_width, padded_image);
}

bool EtcCompressor::CompressAndPad(CompressedImage::Format format,
                                   uint32 height, uint32 width,
                                   uint32 padded_height, uint32 padded_width,
                                   uint32 padding_bytes_per_row,
                                   const uint8 *buffer,
                                   CompressedImage *padded_image) {
  if (!buffer || !padded_image || height == 0 || width == 0 ||
      format != CompressedImage::kRGB) {
    return false;
  }
  return Compressor4x4Helper::CompressAndPad<uint64, Rgb888>(
      EtcEncode(compression_strategy_), "etc", format, height, width,
      padded_height, padded_width, padding_bytes_per_row, buffer, padded_image);
}

bool EtcCompressor::CreateSolidImage(CompressedImage::Format format,
                                     uint32 height, uint32 width,
                                     const uint8 *color,
                                     CompressedImage *image) {
  if (!image || format != CompressedImage::kRGB)
    return false;

  return Compressor4x4Helper::CreateSolidImage<uint64>(
      "etc", format, height, width,
      CreateSolidBlock(Rgb888(color[0], color[1], color[2])),
      image);
}

bool EtcCompressor::CopySubimage(const CompressedImage &image,
                                 uint32 start_row, uint32 start_column,
                                 uint32 height, uint32 width,
                                 CompressedImage *subimage) {
  if (!IsValidCompressedImage(image) || !subimage)
    return false;

  return Compressor4x4Helper::CopySubimage<uint64>(
      image, start_row, start_column, height, width, subimage);
}

}  // namespace image_codec_compression
