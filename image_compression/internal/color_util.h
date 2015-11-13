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

#ifndef IMAGE_COMPRESSION_INTERNAL_COLOR_UTIL_H_
#define IMAGE_COMPRESSION_INTERNAL_COLOR_UTIL_H_

#include <algorithm>

#include "base/integral_types.h"
#include "base/logging.h"
#include "image_compression/internal/color_types.h"

namespace image_codec_compression {

//
// This file contains functions that help deal with colors in image compression
// and decompression functions.
//

//-----------------------------------------------------------------------------

//
// Binary arithmetic operators on integer color types.
//

inline const RgbInt operator+(const RgbInt& c0, const RgbInt& c1) {
  return RgbInt(c0.r + c1.r, c0.g + c1.g, c0.b + c1.b);
}

inline const RgbaInt operator+(const RgbaInt& c0, const RgbaInt& c1) {
  return RgbaInt(c0.r + c1.r, c0.g + c1.g, c0.b + c1.b, c0.a + c1.a);
}

inline const RgbInt operator-(const RgbInt& c0, const RgbInt& c1) {
  return RgbInt(c0.r - c1.r, c0.g - c1.g, c0.b - c1.b);
}

inline const RgbaInt operator-(const RgbaInt& c0, const RgbaInt& c1) {
  return RgbaInt(c0.r - c1.r, c0.g - c1.g, c0.b - c1.b, c0.a - c1.a);
}

inline const RgbInt operator/(const RgbInt& c, int32 div) {
  return RgbInt(c.r / div, c.g / div, c.b / div);
}

inline const RgbaInt operator/(const RgbaInt& c, int32 div) {
  return RgbaInt(c.r / div, c.g / div, c.b / div, c.a / div);
}

//-----------------------------------------------------------------------------

//
// Conversion functions.
//

// Converts an Rgb888 or Rgba8888 color to an RgbInt color. (Rgba8888
// drops the alpha.) These are overloaded to allow use in templates.
template <typename ColorType> inline const RgbInt ToRgbInt(const ColorType &c) {
  return RgbInt(c.r, c.g, c.b);
}

// Converts an Rgba8888 color to an RgbaInt color.
inline const RgbaInt ToRgbaInt(const Rgba8888 &c) {
  return RgbaInt(c.r, c.g, c.b, c.a);
}

// Converts an RgbInt color in which each component is assumed to be
// 0-255 to an Rgb888.
inline const Rgb888 ToRgb888(const RgbInt &c) {
  return Rgb888(c.r, c.g, c.b);
}

// Converts an RgbaInt color in which each component is assumed to be
// 0-255 to an Rgba8888.
inline const Rgba8888 ToRgba8888(const RgbaInt &c) {
  return Rgba8888(c.r, c.g, c.b, c.a);
}

// Converts an Rgb565 to a uint16.
inline uint16 ToUInt16(const Rgb565& c) {
  return (c.r << 11 |
          c.g << 5 |
          c.b);
}

// Converts a uint16 to an Rgb565.
inline const Rgb565 ToRgb565(const uint16 p) {
  return Rgb565(p >> 11,
                (p >> 5) & 0x3F,
                p & 0x1F);
}

// Swaps the red and blue channels in an RgbInt color. This makes it
// easier to deal with DirectX images, which store colors as BGR.
inline const RgbInt ToBgrInt(const RgbInt &c) {
  return RgbInt(c.b, c.g, c.r);
}

// Swaps the red and blue channels in an RgbaInt color. This makes it
// easier to deal with DirectX images, which store colors as BGRA.
inline const RgbaInt ToBgraInt(const RgbaInt &c) {
  return RgbaInt(c.b, c.g, c.r, c.a);
}

// Convenience function that swaps red and blue channels in an RgbInt
// color if a passed flag is true.
inline const RgbInt ToRgbOrBgrInt(const RgbInt &c, bool swap_red_and_blue) {
  return swap_red_and_blue ? ToBgrInt(c) : c;
}

// Convenience function that swaps red and blue channels in an RgbaInt
// color if a passed flag is true.
inline const RgbaInt ToRgbaOrBgraInt(const RgbaInt &c, bool swap_red_and_blue) {
  return swap_red_and_blue ? ToBgraInt(c) : c;
}

// Templated color swap function.
template <typename ColorType> inline void SwapRedAndBlue(ColorType *c) {
  std::swap(c->r, c->b);
}

//-----------------------------------------------------------------------------

//
// Basic compression/decompression functions.
//

// Quantizes a color component assumed to be in the range 0-255 to an
// N-bit unsigned value, where 0<N<8.  This just shifts the value, so
// it is fast but not super-precise.
template <int num_bits> inline const uint8 Quantize8Fast(int v) {
  DCHECK_GT(num_bits, 0);
  DCHECK_LT(num_bits, 8);
  DCHECK_GE(v, 0);
  DCHECK_LT(v, 256);
  return v >> (8 - num_bits);
}


// Quantizes a color component assumed to be in the range 0-255 to an N-bit
// unsigned value, where 0<N<8.  This results in the same value as you would get
// from using the floating point conversion: round(v / 255.0 * max_val).
// This trick comes from Jim Blinn's Dirty Pixels book, Chapter 19, "Three
// Wrongs Make a Right".
template <int num_bits> inline const uint8 Quantize8(int v) {
  DCHECK_GT(num_bits, 0);
  DCHECK_LT(num_bits, 8);
  DCHECK_GE(v, 0);
  DCHECK_LT(v, 256);
  int max_val = (1 << num_bits) - 1;
  int i = v * max_val + 128;
  return (i + (i >> 8)) >> 8;
}


// Quantizes an RGB color with 8 bits per component to N-bit components.
template <int num_bits> inline const RgbInt QuantizeRgbFast(const RgbInt &c) {
  return RgbInt(Quantize8Fast<num_bits>(c.r),
                Quantize8Fast<num_bits>(c.g),
                Quantize8Fast<num_bits>(c.b));
}

// Quantizes an RGBA color with 8 bits per component to N-bit components.
template <int num_bits>
inline const RgbaInt QuantizeRgbaFast(const RgbaInt &c) {
  return RgbaInt(Quantize8Fast<num_bits>(c.r),
                 Quantize8Fast<num_bits>(c.g),
                 Quantize8Fast<num_bits>(c.b),
                 Quantize8Fast<num_bits>(c.a));
}

// Quantizes an RgbInt color in which each component is assumed to be
// 0-255 to an Rgb565.
inline const Rgb565 QuantizeToRgb565(const RgbInt &c) {
  return Rgb565(Quantize8<5>(c.r),
                Quantize8<6>(c.g),
                Quantize8<5>(c.b));
}

// Extends a 4-bit color component to an 8-bit component by
// replicating all 4 bits. For example: '1011' -> '10111011'.
inline int Extend4Bit(int bits) {
  return (bits << 4) | bits;
}

// Extends a 5-bit color component to an 8-bit component by
// replicating the 3 highest-order bits as the new low-order bits.
// For example: '10110' -> '10110101'.
inline int Extend5Bit(int bits) {
  return (bits << 3) | ((bits >> 2) & 7);
}

// Converts an Rgb565 to an RgbInt color quickly.
// This conversion taken from:
//   http://developer.download.nvidia.com/compute/cuda/3_0/sdk/website/OpenCL/website/OpenCL/src/oclDXTCompression/doc/opencl_dxtc.pdf
// In particular it is claimed that this is how the hardware expands the 565
// colors. The approximation is off-by-1 from the standard rounded floating
// point conversion for the following values:
//  5-bit Fast  Float
//      3   24   25
//      7   57   58
//     24  198  197
//     28  231  230
//  6-bit Fast  Float
//     11   44   45
//     12   48   49
//     13   52   53
//     14   56   57
//     15   60   61
//     48  195  194
//     49  199  198
//     50  203  202
//     51  207  206
//     52  211  210
inline const RgbInt ExtendToRgbInt(const Rgb565 &c) {
  return RgbInt((c.r << 3) | (c.r >> 2),
                (c.g << 2) | (c.g >> 4),
                (c.b << 3) | (c.b >> 2));
}

// Extends an Rgb565 to an Rgb888.
inline const Rgb888 ExtendToRgb888(const Rgb565 &c) {
  return Rgb888((c.r << 3) | (c.r >> 2),
                (c.g << 2) | (c.g >> 4),
                (c.b << 3) | (c.b >> 2));
}

//-----------------------------------------------------------------------------

//
// Arithmetic operations on colors, using only integer math.  Some of
// these sacrifice a little bit of precision for speed, which is
// usually the right tradeoff during compression and decompression.
//

// Clamps an integer to 0-255 and returns the result.
inline int ClampTo8Bits(int value) {
  // This could more easily and clearly be written as:
  //   return std::max(0, std::min(255, value));
  //
  // However, this function is used in some super-time-critical
  // functions, and the sign-bit trick used here actually speeds
  // things up quite a bit.  Note that doing a similar trick to remove
  // the call to max actually resulted in slower code:
  //   const int one_if_max =
  //       (static_cast<uint32>(255 - value) & 0x80000000) >> 31;
  //   return zero_if_negative * ((value | (one_if_max * 255)) & 255);
  //
  // Also note that moving the sign-bit extraction into a separate
  // inline function resulted in slower code as well.  Go figure.
  const int zero_if_negative =
      1 - ((static_cast<uint32>(value) & 0x80000000) >> 31);
  return zero_if_negative * std::min(255, value);
}

// Clamps all components of an RgbInt to 0-255 and returns the result.
inline const RgbInt ClampRgbInt(const RgbInt &c) {
  return RgbInt(ClampTo8Bits(c.r), ClampTo8Bits(c.g), ClampTo8Bits(c.b));
}

// Clamps all components of an RgbaInt to 0-255 and returns the result.
inline const RgbaInt ClampRgbaInt(const RgbaInt &c) {
  return RgbaInt(ClampTo8Bits(c.r), ClampTo8Bits(c.g),
                 ClampTo8Bits(c.b), ClampTo8Bits(c.a));
}

// Linear combination of two integer values, using all integer math.
// The sum of the integer scales is used as the divisor for each
// product.  For example, given scales of 3 and 8, the linear
// combination is 3/11 of value0 and 8/11 of value1.
template<int scale0, int scale1>
inline int CombineIntFast(int value0, int value1) {
  DCHECK_NE(scale0 + scale1, 0);
  return (scale0 * value0 + scale1 * value1) / (scale0 + scale1);
}

// Version of CombineIntFast() that operates on uint8 values, avoiding
// overflow.
template<int scale0, int scale1>
inline uint8 CombineUint8Fast(uint8 value0, uint8 value1) {
  return static_cast<uint8>(CombineIntFast<scale0, scale1>(value0, value1));
}

// Version of CombineIntFast() that operates on Rgb888 colors.
template<int scale0, int scale1>
inline const Rgb888 CombineRgb888Fast(const Rgb888 &color0,
                                      const Rgb888 &color1) {
  return Rgb888(CombineUint8Fast<scale0, scale1>(color0.r, color1.r),
                CombineUint8Fast<scale0, scale1>(color0.g, color1.g),
                CombineUint8Fast<scale0, scale1>(color0.b, color1.b));
}

// Version of CombineIntFast() that operates on Rgba8888 colors.
template<int scale0, int scale1>
inline const Rgba8888 CombineRgba8888Fast(const Rgba8888 &color0,
                                          const Rgba8888 &color1) {
  return Rgba8888(CombineUint8Fast<scale0, scale1>(color0.r, color1.r),
                  CombineUint8Fast<scale0, scale1>(color0.g, color1.g),
                  CombineUint8Fast<scale0, scale1>(color0.b, color1.b),
                  CombineUint8Fast<scale0, scale1>(color0.a, color1.a));
}

// Version of CombineIntFast() that operates on RgbInt colors.
template<int scale0, int scale1>
inline const RgbInt CombineRgbIntFast(const RgbInt &color0,
                                      const RgbInt &color1) {
  return RgbInt(CombineIntFast<scale0, scale1>(color0.r, color1.r),
                CombineIntFast<scale0, scale1>(color0.g, color1.g),
                CombineIntFast<scale0, scale1>(color0.b, color1.b));
}

// Version of CombineIntFast() that operates on RgbaInt colors.
template<int scale0, int scale1>
inline const RgbaInt CombineRgbaIntFast(const RgbaInt &color0,
                                        const RgbaInt &color1) {
  return RgbaInt(CombineIntFast<scale0, scale1>(color0.r, color1.r),
                 CombineIntFast<scale0, scale1>(color0.g, color1.g),
                 CombineIntFast<scale0, scale1>(color0.b, color1.b),
                 CombineIntFast<scale0, scale1>(color0.a, color1.a));
}

// Averaging 4 uint8 values using only integer math, which is useful
// during image reduction.
inline uint8 Average4Uint8Fast(uint8 value0, uint8 value1,
                               uint8 value2, uint8 value3) {
  return static_cast<uint8>((static_cast<int>(value0) +
                             static_cast<int>(value1) +
                             static_cast<int>(value2) +
                             static_cast<int>(value3)) / 4);
}

// Templated averaging function for 4 colors.
template <typename ColorType>
inline const ColorType Average4ColorsFast(const ColorType &color0,
                                          const ColorType &color1,
                                          const ColorType &color2,
                                          const ColorType &color3);

// Specialized version for Rgb888 colors.
template <> inline const Rgb888 Average4ColorsFast(const Rgb888 &color0,
                                                   const Rgb888 &color1,
                                                   const Rgb888 &color2,
                                                   const Rgb888 &color3) {
  return Rgb888(Average4Uint8Fast(color0.r, color1.r, color2.r, color3.r),
                Average4Uint8Fast(color0.g, color1.g, color2.g, color3.g),
                Average4Uint8Fast(color0.b, color1.b, color2.b, color3.b));
}

// Specialized version for Rgba8888 colors.
template <> inline const Rgba8888 Average4ColorsFast(const Rgba8888 &color0,
                                                     const Rgba8888 &color1,
                                                     const Rgba8888 &color2,
                                                     const Rgba8888 &color3) {
  return Rgba8888(Average4Uint8Fast(color0.r, color1.r, color2.r, color3.r),
                  Average4Uint8Fast(color0.g, color1.g, color2.g, color3.g),
                  Average4Uint8Fast(color0.b, color1.b, color2.b, color3.b),
                  Average4Uint8Fast(color0.a, color1.a, color2.a, color3.a));
}

// Templated function to average the colors of the 2x2 block of pixels
// whose upper-left pixel is at the given row and column.
template <typename ColorType>
static ColorType ComputeAveragePixel2x2(const ColorType pixels[4][4],
                                        int row, int col) {
  return Average4ColorsFast<ColorType>(pixels[row][col],
                                       pixels[row][col + 1],
                                       pixels[row + 1][col],
                                       pixels[row + 1][col + 1]);
}

// Computes the approximate luminance of an RGB color.
inline int ComputeLuminanceFast(const RgbInt &c) {
  // The RGB-to-luminance formula is usually:
  //
  //    Luminance = 0.299 * red + 0.587 * green + 0.114 * blue
  //
  // Since we don't care about exact values, we use integer
  // approximations in roughly the same ratios, but which the compiler
  // can optimize very easily.
  static const int kRedFactor   = 4;
  static const int kGreenFactor = 8;
  static const int kBlueFactor  = 1;
  return (c.r * kRedFactor + c.g * kGreenFactor + c.b * kBlueFactor);
}

// Computes the distance between two RGB colors using a
// squared luminance-based metric.
inline int ComputeSquaredLuminanceDistanceFast(const RgbInt &c0,
                                               const RgbInt &c1) {
  const int diff = ComputeLuminanceFast(c1) - ComputeLuminanceFast(c0);
  return diff * diff;
}

// Computes the distance between two RGB colors using a
// luminance-of-difference metric.  This computes the luminance of the
// absolute difference between two colors rather than the difference of the
// luminance.  The latter is subject to confusing two iso-luminant but
// chromatically very different colors.
inline int ComputeDifferenceLuminanceFast(const RgbInt &c0,
                                          const RgbInt &c1) {
  const RgbInt cdiff(std::abs(c0.r - c1.r),
                     std::abs(c0.g - c1.g),
                     std::abs(c0.b - c1.b));
  const int diff_luminance = ComputeLuminanceFast(cdiff);
  return diff_luminance * diff_luminance;
}

// Computes the squared distance between two color component values.
inline int ComputeSquaredComponentDistance(int component0, int component1) {
  const int diff = component1 - component0;
  return diff * diff;
}

//-----------------------------------------------------------------------------

//
// Helper functions for dealing with colors in image buffers.
//

// Returns a mutable pointer to the first RGB or RGBA color value
// corresponding to a given row and column in an image buffer.
template <typename ColorType>
inline ColorType* GetMutableColorInImageBuffer(
    ColorType *image_buffer, int bytes_per_row, int row, int column) {
  uint8* row_start =
      reinterpret_cast<uint8*>(image_buffer) + (row * bytes_per_row);
  return reinterpret_cast<ColorType*>(row_start) + column;
}

// Returns a const pointer to the first RGB or RGBA color value
// corresponding to a given row and column in an image buffer.
template <typename ColorType>
inline const ColorType* GetColorInImageBuffer(
    const ColorType *image_buffer, int bytes_per_row, int row, int column) {
  const uint8* row_start =
      reinterpret_cast<const uint8*>(image_buffer) + (row * bytes_per_row);
  return reinterpret_cast<const ColorType*>(row_start) + column;
}

// Stores RGB or RGBA colors from a 4x4 array into an image buffer of
// the given width starting at the given row and column. The
// num_rows_in_block and num_cols_in_block parameters indicate how
// many rows and columns of the block should be stored.
template <typename ColorType>
inline void StoreColorsInBuffer(const ColorType block_pixels[4][4],
                                int image_width,
                                int num_rows_in_block, int num_cols_in_block,
                                int row, int col, ColorType *image_buffer) {
  for (int y = 0; y < num_rows_in_block; ++y) {
    ColorType *row_pixels = GetMutableColorInImageBuffer(
        image_buffer, image_width, row + y, col);
    for (int x = 0; x < num_cols_in_block; ++x)
      row_pixels[x] = block_pixels[y][x];
  }
}

//-----------------------------------------------------------------------------

//
// Streaming operators for test output.
//

template <typename T>
inline std::ostream &operator <<(std::ostream &out, const RgbColor<T> &c) {
  out << "[" << c.r << ", " << c.g << ", " << c.b << "]";
  return out;
}

template <typename T>
inline std::ostream &operator <<(std::ostream &out, const RgbaColor<T> &c) {
  out << "[" << c.r << ", " << c.g << ", " << c.b << ", " << c.a << "]";
  return out;
}

}  // namespace image_codec_compression

#endif  // IMAGE_COMPRESSION_INTERNAL_COLOR_UTIL_H_
