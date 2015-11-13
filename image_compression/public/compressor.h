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
// Compressor is a base interface class for classes that implement
// some form of block-based image compression and decompression, such
// as DXTC or ETC.
//
// Input images are restricted to the following:
//   - 8 bits per color component (0-255).
//   - RGB (24 bits per pixel) or RGBA (32 bits per pixel) format.
//   - Row-major order with interleaved color components.
//
// For example, the first 3 bytes of an RGB image are the RGB color
// values for the pixel in the top row and left-most column, the next
// 3 bytes are for the pixel just to the right of that one, and so on.
//
// Several of the functions return a CompressedImage through an
// out-parameter.  For each of these, the caller can pass a pointer to
// a default-constructed instance, in which the function will allocate
// the compressed data inside the instance. Alternatively, the pointer
// can be to an instance constructed using external storage, in which
// case that storage will be used for the result.  The
// ComputeCompressedDataSize() function can be used to determine how
// much storage to allocate when using this approach.

#ifndef IMAGE_COMPRESSION_PUBLIC_COMPRESSOR_H_
#define IMAGE_COMPRESSION_PUBLIC_COMPRESSOR_H_

#include <stddef.h>
#include <vector>

#include "base/integral_types.h"
#include "image_compression/public/compressed_image.h"

namespace image_codec_compression {

class Compressor {
 public:
  virtual ~Compressor() {}

  // Returns true if the compressor supports compressing images with
  // the given format.
  virtual bool SupportsFormat(CompressedImage::Format format) const = 0;

  // This can be used to validate a CompressedImage instance to ensure
  // that it can be processed by a derived instance. It returns false
  // if the compressor type does not match, if the sizes are
  // inconsistent for a compressed image of that type, or if the data
  // vector has the wrong size.
  virtual bool IsValidCompressedImage(const CompressedImage &image) = 0;

  // Returns the size in bytes of the image data that would result
  // from creating a compressed image with the given format and
  // dimensions. This allows callers to construct a CompressedImage
  // with external storage of the correct size to pass to various
  // functions.
  virtual size_t ComputeCompressedDataSize(CompressedImage::Format format,
                                           uint32 height, uint32 width) = 0;

  // Compresses an image of the given format represented by the given
  // buffer data and image sizes, storing the results in the given
  // CompressedImage instance.  The padding_bytes_per_row parameter
  // indicates how much additional padding each row of the original
  // image contains; 0 means the rows are contiguous. Returns false on
  // error.
  virtual bool Compress(CompressedImage::Format format,
                        uint32 height, uint32 width,
                        uint32 padding_bytes_per_row,
                        const uint8 *buffer, CompressedImage *image) = 0;

  // Decompresses a compressed image. The destination buffer vector
  // will be resized correctly to contain the decompressed data.
  // Returns false on error.
  virtual bool Decompress(const CompressedImage &image,
                          std::vector<uint8> *decompressed_buffer) = 0;

  // This function downsamples a compressed image to half its height
  // and half its width; the size of each resulting dimension will be
  // (oldsize + 1)/2. It may be used for example to generate mipmaps
  // of a compressed image without having to decompress it first.
  // Derived classes may impose restrictions on the types and sizes of
  // images that can be downsampled.  Returns false if downsampling is
  // not possible or if an error occurs.
  virtual bool Downsample(const CompressedImage &image,
                          CompressedImage *downsampled_image) = 0;

  // This function can be used to pad a compressed image to a new
  // height and width, replicating the last (bottom) row and last
  // (right) column of pixels as necessary.  Derived classes may
  // impose restrictions on the types and sizes of images that can be
  // padded. If a padded size is not larger than the original size, no
  // padding is done in that dimension. Returns false if padding is
  // not possible or if an error occurs.
  virtual bool Pad(const CompressedImage &image,
                   uint32 padded_height, uint32 padded_width,
                   CompressedImage *padded_image) = 0;

  // Compresses and pads in a single operation. Some compressors may
  // implement these more efficiently when done together. The same
  // rules and restrictions for Compress() and Pad() apply here.
  // Images resulting from this operation may vary slightly from those
  // from applying compression and padding separately.
  virtual bool CompressAndPad(CompressedImage::Format format,
                              uint32 height, uint32 width,
                              uint32 padded_height, uint32 padded_width,
                              uint32 padding_bytes_per_row,
                              const uint8 *buffer,
                              CompressedImage *padded_image) = 0;

  // This function can be used to create a compressed image of the
  // given format and size with a solid color represented by the given
  // color (a single-pixel buffer). Returns false if this is not
  // possible.
  virtual bool CreateSolidImage(CompressedImage::Format format,
                                uint32 height, uint32 width, const uint8 *color,
                                CompressedImage *image) = 0;

  // Given a compressed image, this copies a subimage region defined
  // by starting locations and sizes. The region must be completely
  // contained within the original image.  Derived classes may impose
  // other restrictions on the locations and sizes.  Returns false on
  // error.
  virtual bool CopySubimage(const CompressedImage &image,
                            uint32 start_row, uint32 start_column,
                            uint32 height, uint32 width,
                            CompressedImage *subimage) = 0;
};

}  // namespace image_codec_compression

#endif  // IMAGE_COMPRESSION_PUBLIC_COMPRESSOR_H_
