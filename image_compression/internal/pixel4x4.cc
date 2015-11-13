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

#include "image_compression/internal/pixel4x4.h"

#include <algorithm>

#include "base/integral_types.h"

namespace image_codec_compression {

template <typename ColorType>
void Pixel4x4::ConstructOutsideImage(const ColorType *pixels,
                                     uint32 image_height, uint32 image_width,
                                     uint32 bytes_per_row,
                                     int row, int column) {
  //
  // There are 4 potential rectangular pixel regions to fill in.
  //   Region 1: Within the image.
  //   Region 2: To the right of the image; filled in with replicated
  //             rightmost column of pixels.
  //   Region 3: Below the image; filled in with replicated last row
  //             of pixels.
  //   Region 4: Below and to right of the image; filled in with
  //             replicated lower-rightmost corner pixel.
  //
  int height = image_height;
  int width  = image_width;
  int rows_within_image = std::max(0, height - row);
  int columns_within_image = std::max(0, width  - column);
  int x_max = std::min(4, columns_within_image);

  for (int y = 0; y < 4; ++y) {
    // Get the row of pixels to copy from.
    const ColorType *row_pixels =
        GetImagePixel(pixels, bytes_per_row, std::min(row + y, height - 1), 0);
    // Regions 1 and 3:
    for (int x = 0; x < x_max; ++x)
      SetPixel(y, x, row_pixels[column + x]);
    // Regions 2 and 4:
    for (int x = columns_within_image; x < 4; ++x)
      SetPixel(y, x, row_pixels[std::min(column + x, width - 1)]);
  }

  // Set this for special case of single-color 4x4 block, which helps
  // optimize compression.
  has_one_pixel_ = (columns_within_image == 0 && rows_within_image == 0);
}

// Instantiate the function for supported pixel types.
template void Pixel4x4::ConstructOutsideImage<Rgb888>(
    const Rgb888 *pixels, uint32 image_height, uint32 image_width,
    uint32 bytes_per_row, int row, int column);
template void Pixel4x4::ConstructOutsideImage<Rgba8888>(
    const Rgba8888 *pixels, uint32 image_height, uint32 image_width,
    uint32 bytes_per_row, int row, int column);

}  // namespace image_codec_compression
