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

#include "image_compression/internal/compressor4x4_helper.h"

#include "base/integral_types.h"
#include "base/logging.h"

namespace image_codec_compression {

bool Compressor4x4Helper::SetUpCompressedImage(
    const std::string &compressor_name, size_t block_size,
    CompressedImage::Format format, uint32 height, uint32 width,
    uint32 padding_bytes_per_row, CompressedImage *image) {
  DCHECK(image);
  const uint32 num_block_rows = GetNumBlocks(height);
  const uint32 num_block_cols = GetNumBlocks(width);
  const size_t data_size = num_block_rows * num_block_cols * block_size;

  const CompressedImage::Metadata metadata(
      format, compressor_name, height, width,
      4 * num_block_rows, 4 * num_block_cols, padding_bytes_per_row);
  if (image->OwnsData()) {
    image->CreateOwnedData(metadata, data_size);
  } else {
    // Make sure the external storage has the correct size.
    if (image->GetDataSize() != data_size)
      return false;
    image->SetMetadata(metadata);
  }
  return true;
}

}  // namespace image_codec_compression
