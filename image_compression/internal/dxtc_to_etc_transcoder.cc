// Copyright 2016 Google Inc. All Rights Reserved.
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

#include "image_compression/internal/color_util.h"
#include "image_compression/internal/color_types.h"
#include "image_compression/internal/pixel4x4.h"
#include "image_compression/public/etc_compressor.h"

namespace image_codec_compression {

// Forward declarations of objects only seen in .cc files.
struct Dxt1Block;
void DecodeDxt1Block(const Dxt1Block &dxt1_block, bool swap_red_and_blue,
                     Rgb888 decoded_pixels[4][4]);
uint64 EncodeEtc1Block(const Pixel4x4 &pixel4x4,
                       EtcCompressor::CompressionStrategy strategy);

void TranscodeDxt1ToEtc1(CompressedImage* image) {
  uint8* imgbytes = image->GetMutableData();
  Rgb888 rgb4x4[4][4];
  for (size_t i = 0, n = image->GetDataSize(); i < n; i += 8) {
    Dxt1Block* dxt1_block = reinterpret_cast<Dxt1Block*>(imgbytes);
    DecodeDxt1Block(*dxt1_block, false, rgb4x4);
    Pixel4x4 pixel4x4(&rgb4x4[0][0], 4, 4, 0, 0, 0);
    *(reinterpret_cast<uint64*>(imgbytes)) =
        EncodeEtc1Block(pixel4x4, EtcCompressor::kHeuristic);
    imgbytes += 8;
  }
}

}  // namespace image_codec_compression
