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

#ifndef IMAGE_COMPRESSION_INTERNAL_DXTC_CONST_COLOR_TABLE_H_
#define IMAGE_COMPRESSION_INTERNAL_DXTC_CONST_COLOR_TABLE_H_

#include "image_compression/internal/color_util.h"

namespace image_codec_compression {

// Sets color0 and color1 to the best end point colors to use to represent the
// constant block color, target_color.  Returns the 2-bit value to use to encode
// the pixels in the block.
// If always_4_color_case is false, will also use 1/2-way interpolation rules
// to reduce errors further.
int GetBestDxtcConstColors(const RgbInt& target_color,
                           Rgb565* color0, Rgb565* color1,
                           bool always_4_color_case);

}  // namespace image_codec_compression

#endif  // IMAGE_COMPRESSION_INTERNAL_DXTC_CONST_COLOR_TABLE_H_
