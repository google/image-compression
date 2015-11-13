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

#ifndef IMAGE_COMPRESSION_INTERNAL_PIXEL4X4_H_
#define IMAGE_COMPRESSION_INTERNAL_PIXEL4X4_H_

#include "base/integral_types.h"
#include "base/logging.h"
#include "image_compression/internal/color_util.h"

namespace image_codec_compression {

//
// The Pixel4x4 class is used to simplify access to a 4x4 block of
// pixels in an RGB or RGBA source image. The pixels in the block are
// first converted to RgbInts to provide a consistent and simple
// interface.
//
class Pixel4x4 {
 public:
  // The default constructor is used when clients will set pixel
  // values in a block directly.
  Pixel4x4()
      : has_one_pixel_(false) {
  }

  // Constructs a Pixel4x4 from Rgb888 or Rgba8888 image data. A
  // pointer to the image data is given along with the image size and
  // stride. The row and column parameters indicate the top-left
  // corner of the block, with rows increasing from top to bottom and
  // columns increasing from left to right.  The Rgb888 version leaves
  // the alpha_ member undefined.
  template <typename ColorType>
  inline Pixel4x4(const ColorType *pixels, uint32 height, uint32 width,
                  uint32 padding_bytes_per_row, int row, int column) {
    const uint32 bytes_per_row = width * sizeof(pixels[0]) +
                                 padding_bytes_per_row;

    // In the common case where the 4x4 lies completely within the
    // image, do all the construction work here.
    if (static_cast<int32>(height) - row >= 4 &&
        static_cast<int32>(width) - column >= 4) {
      for (int y = 0; y < 4; ++y) {
        const ColorType *row_pixels = GetImagePixel(pixels, bytes_per_row,
                                                    row + y, column);
        for (int x = 0; x < 4; ++x) {
          pixels_[y][x] = ToRgbInt(row_pixels[x]);
          SetAlpha(y, x, row_pixels[x]);
        }
      }
      has_one_pixel_ = false;
    } else {
      ConstructOutsideImage<ColorType>(pixels, height, width, bytes_per_row,
                                       row, column);
    }
  }

  // Stores a pixel color at the given row/column location.
  template <typename ColorType>
  inline void SetPixel(int row, int column, const ColorType &color) {
    pixels_[row][column] = ToRgbInt(color);
    SetAlpha(row, column, color);
  }

  // Stores an alpha value. The default implementation does nothing.
  // A template specialization below is used to store alpha if the
  // block contains RGBA pixels.
  template <typename ColorType>
  inline void SetAlpha(int row, int column, const ColorType &color) {
  }

  // Returns true if all pixels in the block are known to have the
  // same color value. This is a special case that is used to optimize
  // when compressing into an image with larger dimensions than the
  // original.
  inline bool has_one_pixel() const {
    return has_one_pixel_;
  }

  // Returns the color of a pixel.
  inline const RgbInt& GetPixel(int row, int column) const {
    DCHECK_GE(row, 0);
    DCHECK_GE(column, 0);
    DCHECK_LT(row, 4);
    DCHECK_LT(column, 4);
    return pixels_[row][column];
  }

  // Accesses an alpha value. Values range from 0 to 255 for RGBA
  // images, and are likely to be undefined for RGB images.
  inline int GetAlpha(int row, int column) const {
    DCHECK_GE(row, 0);
    DCHECK_GE(column, 0);
    DCHECK_LT(row, 4);
    DCHECK_LT(column, 4);
    return alpha_[row][column];
  }

 private:
  // Returns a pointer to the pixel at a given row and column in an
  // image whose first pixel and row size are provided. This is
  // templated to work with RGB and RGBA images.
  template <typename ColorType>
  static inline const ColorType* GetImagePixel(const ColorType *pixels,
                                               int bytes_per_row,
                                               int row, int column) {
    const uint8 *row_start =
        reinterpret_cast<const uint8*>(pixels) + (row * bytes_per_row);
    return reinterpret_cast<const ColorType*>(row_start) + column;
  }

  // Does the work of the constructor when the 4x4 region is not
  // totally contained within the image, replicating pixels where
  // necessary.
  template <typename ColorType>
  void ConstructOutsideImage(const ColorType *pixels,
                             uint32 image_height, uint32 image_width,
                             uint32 bytes_per_row, int row, int column);

  // Pixels as colors with integer components.
  RgbInt pixels_[4][4];

  // Alpha values of all pixels (for RGBA only). These are stored
  // separately because they are accessed separately.
  int alpha_[4][4];

  // True if the block is known to consist of a single pixel value,
  // which makes it trivial (and fast) to encode.
  bool has_one_pixel_;
};

// Specialization of SetAlpha() for Rgba8888.
template <>
inline void Pixel4x4::SetAlpha(int row, int column, const Rgba8888 &color) {
  alpha_[row][column] = color.a;
}

// Given a 4x4 block of pixels, downsamples them to a 2x2 and stores
// the results in one 2x2 corner (whose upper left pixel is at
// target_row, target_col) of the given Pixel4x4.
template <typename ColorType>
inline void StoreDownsampledPixels4x4(const ColorType pixels[4][4],
                                      int target_row, int target_col,
                                      Pixel4x4 *pixel4x4) {
  for (int row = 0; row < 2; ++row) {
    for (int col = 0; col < 2; ++col) {
      pixel4x4->SetPixel(target_row + row, target_col + col,
                         ComputeAveragePixel2x2(pixels, 2 * row, 2 * col));
    }
  }
}

}  // namespace image_codec_compression

#endif  // IMAGE_COMPRESSION_INTERNAL_PIXEL4X4_H_
