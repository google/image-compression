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
//
// DxtcCompressor is a derived Compressor class that implements DXTC
// image compression and decompression.  It supports DXT1 compression
// for RGB images and DXT5 compression for RGBA images.
//
// Compression should work correctly on source images of any size.
// Since DXT compression operates on 4x4 blocks of pixels, the
// compressed output image always has dimensions that are multiples of
// 4.  When compressing, if either of the original image dimensions is
// not a multiple of 4, the last row or column of pixels in the
// original image is replicated to simulate the correct size in the
// compressed output.
//
// The compressed_height and compressed_width members of the
// CompressionSpec can be used to specify a different size for the
// resulting compressed image. For example, they can be used to
// compress to an image whose dimensions are (larger) powers of 2.
// Results in such cases are similar to those of the Pad() function.
//
// Decompression has similar sizing rules, allowing a subimage of a
// compressed image to be decompressed into an image of the
// appropriate size.
//
// Reduction of a compressed image requires both image dimensions be a
// multiple of 8.

#ifndef IMAGE_COMPRESSION_PUBLIC_DXTC_COMPRESSOR_H_
#define IMAGE_COMPRESSION_PUBLIC_DXTC_COMPRESSOR_H_

#include <stddef.h>
#include <vector>

#include "base/integral_types.h"
#include "image_compression/public/compressed_image.h"
#include "image_compression/public/compressor.h"

namespace image_codec_compression {

class DxtcCompressor : public Compressor {
 public:
  virtual ~DxtcCompressor();
  virtual bool SupportsFormat(CompressedImage::Format format) const;
  virtual bool IsValidCompressedImage(const CompressedImage &image);
  virtual size_t ComputeCompressedDataSize(CompressedImage::Format format,
                                           uint32 height, uint32 width);
  virtual bool Compress(CompressedImage::Format format,
                        uint32 height, uint32 width,
                        uint32 padding_bytes_per_row,
                        const uint8 *buffer, CompressedImage *image);
  virtual bool Decompress(const CompressedImage &image,
                          std::vector<uint8> *decompressed_buffer);
  virtual bool Downsample(const CompressedImage &image,
                          CompressedImage *downsampled_image);
  virtual bool Pad(const CompressedImage &image,
                   uint32 padded_height, uint32 padded_width,
                   CompressedImage *padded_image);
  virtual bool CompressAndPad(CompressedImage::Format format,
                              uint32 height, uint32 width,
                              uint32 padded_height, uint32 padded_width,
                              uint32 padding_bytes_per_row,
                              const uint8 *buffer,
                              CompressedImage *padded_image);
  virtual bool CreateSolidImage(CompressedImage::Format format,
                                uint32 height, uint32 width, const uint8 *color,
                                CompressedImage *image);
  virtual bool CopySubimage(const CompressedImage &image,
                            uint32 start_row, uint32 start_column,
                            uint32 height, uint32 width,
                            CompressedImage *subimage);
};

}  // namespace image_codec_compression

#endif  // IMAGE_COMPRESSION_PUBLIC_DXTC_COMPRESSOR_H_
