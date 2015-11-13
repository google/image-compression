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

#ifndef IMAGE_COMPRESSION_PUBLIC_COMPRESSED_IMAGE_H_
#define IMAGE_COMPRESSION_PUBLIC_COMPRESSED_IMAGE_H_

#include <stddef.h>
#include <cstring>
#include <string>

#include "base/integral_types.h"
#include "base/logging.h"

namespace image_codec_compression {

// This class is used to represent an image that has been compressed
// using a Compressor. It encapsulates all information about the image
// to make decompression and other operations possible. An instance
// can be created either with internally- or externally-managed
// storage.
class CompressedImage {
 public:
  // Supported image format types.
  enum Format {
    kRGB,   // Red, green, blue.
    kBGR,   // Blue, green, red (used by DirectX).
    kRGBA,  // Red, green, blue, alpha.
    kBGRA,  // Blue, green, red, alpha (used by DirectX).
  };

  // This struct encapsulates all the metadata of a compressed image.
  struct Metadata {
    // The constructor is passed all values.
    Metadata(Format format_in,
             const std::string &compressor_name_in,
             uint32 uncompressed_height_in,
             uint32 uncompressed_width_in,
             uint32 compressed_height_in,
             uint32 compressed_width_in,
             uint32 padding_bytes_per_row_in)
        : format(format_in),
          compressor_name(compressor_name_in),
          uncompressed_height(uncompressed_height_in),
          uncompressed_width(uncompressed_width_in),
          compressed_height(compressed_height_in),
          compressed_width(compressed_width_in),
          padding_bytes_per_row(padding_bytes_per_row_in) {
    }

    // Format of the compressed image.
    Format format;

    // Name of the compressor used to compress it.
    std::string compressor_name;

    // Height and width in pixels of the original uncompressed image.
    uint32 uncompressed_height;
    uint32 uncompressed_width;

    // Height and width in pixels of the compressed image. These may
    // differ from the originals due to size requirements imposed by
    // the compression technique.
    uint32 compressed_height;
    uint32 compressed_width;

    // Extra padding bytes in each row of the original uncompressed
    // image.  This is assumed to be the same for the target image
    // during decompression.
    uint32 padding_bytes_per_row;
  };

  // The default constructor creates an empty image.
  CompressedImage()
      : metadata_(kRGB, "", 0, 0, 0, 0, 0),
        data_size_(0),
        data_(NULL),
        owns_data_(true) {
  }

  // This constructor can be used when data storage is managed
  // externally.  It is assumed that the lifetime of the data is at
  // least as long as that of this instance.
  CompressedImage(size_t data_size, uint8 *external_data)
      : metadata_(kRGB, "", 0, 0, 0, 0, 0),
        data_size_(data_size),
        data_(external_data),
        owns_data_(false) {
    DCHECK(external_data);
  }

  // The destructor frees up the data storage if it is owned by this instance.
  ~CompressedImage() {
    if (OwnsData()) {
      delete[] data_;
    }
  }

  // Copies metadata and data from another instance, which must have
  // data. This instance will own its data regardless of the ownership
  // status of the other instance.
  void Duplicate(const CompressedImage &from) {
    // Duplicating from an instance to itself is permissible, but is a no-op if
    // the instance already owns its own data.
    if (&from != this || !OwnsData()) {
      // Save the current data pointer in case creating owned data changes it.
      const uint8 *source_data = from.data_;
      DCHECK(source_data);
      CreateOwnedData(from.metadata_, from.data_size_);
      std::memcpy(data_, source_data, data_size_);
    }
  }

  // Sets the instance to contain the given metadata and a data area
  // of the given size. The data will be freed when the instance is
  // destroyed.
  void CreateOwnedData(const Metadata &metadata, size_t data_size) {
    if (OwnsData())
      delete[] data_;
    metadata_ = metadata;
    data_size_ = data_size;
    data_ = new uint8[data_size_];
    owns_data_ = true;
  }

  // Sets the metadata in the instance. This should be called only in
  // the case where the instance was constructed with external data.
  // Otherwise, CreateOwnedData() should be used.
  void SetMetadata(const Metadata &metadata) {
    DCHECK(!OwnsData());
    metadata_ = metadata;
  }

  // Returns the metadata for the compressed image.
  const Metadata& GetMetadata() const {
    return metadata_;
  }

  // Returns true if the data is owned by this instance.
  bool OwnsData() const {
    return owns_data_;
  }

  // Returns the size of the data.
  size_t GetDataSize() const {
    return data_size_;
  }

  // Returns a const pointer to the data.
  const uint8* GetData() const {
    return data_;
  }

  // Returns a pointer to mutable data.
  uint8* GetMutableData() {
    return data_;
  }

 private:
  Metadata metadata_;

  // Size of the compressed data.
  size_t data_size_;

  // Pointer to compressed image data. The contents and interpretation
  // of the data is compressor-specific.
  uint8 *data_;

  // Indicates whether the data is managed by this instance.
  bool owns_data_;

  // Only the provided Copy functions should be used for copying.
  CompressedImage(const CompressedImage&);
  void operator=(const CompressedImage&);
};

// Returns the number of components in the given format.
inline int GetNumFormatComponents(CompressedImage::Format format) {
  switch (format) {
    case CompressedImage::kRGB:
    case CompressedImage::kBGR:
      return 3;
    case CompressedImage::kRGBA:
    case CompressedImage::kBGRA:
      return 4;
    default:  // To shut up the compiler.
      return 0;
  }
}

// Returns true if a format needs to have red and blue components swapped.
inline bool NeedsRedAndBlueSwapped(CompressedImage::Format format) {
  return format == CompressedImage::kBGR || format == CompressedImage::kBGRA;
}

}  // namespace image_codec_compression

#endif  // IMAGE_COMPRESSION_PUBLIC_COMPRESSED_IMAGE_H_
