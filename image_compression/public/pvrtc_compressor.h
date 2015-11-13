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
// PvrtcCompressor implements PVRTC (PowerVR texture compression). It currently
// only supports compression of RGBA8888 to 2 bits per pixel PVRTC version 1
// RGBA, and will return false for most other method calls (such as
// decompression, padding or downsampling).
//
// PVRTC has two versions, and can compress to 2 or 4 bits per pixel, RGB or
// RGBA. The following is a description of the 2BPP RGBA format used here:
//
// The compressed image can be thought of as containing two 16 BPP images with
// 8x lower width and 4x lower height than the original uncompressed image.
// These images are bilinearly upscaled, where each pixel (x,y) in the lower
// resolution image is copied exactly into pixel (x * 8 + 4, y * 4 + 2) of the
// full resolution image, and the remaining full resolution pixels are
// interpolated. Here these are called the "A" and "B" images.
//
// The compressed image also contains a "modulation" image. For each block of
// 8x4 pixels in the original image, there is either:
// - 1 bit per pixel which tells the decompressor to either use the pixel from
//   the A or B image.
// - 2 bits for every other pixel (in a checkerboard pattern) which tells the
//   decompressor to either use 100% of the A image, 5/8ths A and 3/8ths B,
//   3/8ths A and 5/8ths B, or 100% B image. The remaining pixels in the
//   checkerboard are interpolated, either by averaging adjacent horizontal
//   pixels, or vertical pixels, or an average of all 4 orthogonally adjacent
//   pixels.
//
// The A, B, and modulation images are stored by writing 64-bit words which
// are an A image pixel, a B image pixel, and a block's worth of modulation
// information, and these are all stored in Z-order-curve order. The choice
// of 1bpp or 2bpp modulation and interplotion direction is encoded in these
// words also.
//
// The exact bit format of the 64-bit word is out of scope here; see the
// implementation.
//
// The compression's main speedup is that it chooses colors for each pixel of
// the A and B images by only considering the 8x4 block around the pixel, even
// though the A and B image pixels - through the interpolation - will affect
// 16x8 pixels in the output uncompressed image, and therefore would better be
// chosen in parallel with their neighbors. It also doesn't use principal
// component analysis to choose the A and B pixels, but an approximation.

#ifndef IMAGE_COMPRESSION_PUBLIC_PVRTC_COMPRESSOR_H_
#define IMAGE_COMPRESSION_PUBLIC_PVRTC_COMPRESSOR_H_

#include <stddef.h>
#include <vector>

#include "base/integral_types.h"
#include "image_compression/public/compressed_image.h"
#include "image_compression/public/compressor.h"

namespace image_codec_compression {

// This compressor only supports compression from RGBA8888 to 2BPP PVRTC RGBA.
// All unsupported methods will simply return false.
class PvrtcCompressor : public Compressor {
 public:
  PvrtcCompressor();
  virtual ~PvrtcCompressor();

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

#endif  // IMAGE_COMPRESSION_PUBLIC_PVRTC_COMPRESSOR_H_
