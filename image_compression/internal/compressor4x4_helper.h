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
// Compressor4x4Helper is an internal helper class for derived
// Compressor classes that operate on 4x4 blocks of pixels. Most of
// the methods are templated on the type used to represent the encoded
// 4x4 block of pixels and also the type used to represent a single
// pixel.
//
//   ---------- FUNCTORS ---------
//
// The helper uses templated functors to do the compressor-specific
// work for most of the operations.  The functors and their template
// parameters are named consistently in the interface.  The signatures
// are as follows:
//
// EncodeFunctor:
//       const BlockType func(const Pixel4x4 &pixel4x4, bool swap_red_and_blue);
//   This function is used to encode a single 4x4 block represented by
//   a Pixel4x4 instance. The swap_red_and_blue parameter indicates
//   whether the format is BGR or BGRA in contrast to RGB or RGBA.
//
// DecodeFunctor:
//       void func(const BlockType &block, bool swap_red_and_blue,
//                 ColorType decoded_pixels[4][4]);
//   This function is used to decode a 4x4 block represented by an
//   encoded BlockType instance. The swap_red_and_blue parameter
//   indicates whether the format is BGR or BGRA in contrast to RGB or
//   RGBA.  The decoded pixels are to be stored in the array.
//
// ColumnPadBlockFunctor:
//       const BlockType func(const BlockType &last_column_block);
//   This function is used during calls to Pad() to create a block
//   used to pad to the right of the last block in a row of a
//   compressed image. The returned block should have all 4 pixel
//   columns equal to the rightmost column of last_column_block.
//
// RowPadBlockFunctor:
//       const BlockType func(const BlockType &last_row_block);
//   This function is used during calls to Pad() to create a block
//   used to pad below the last block in a column of a compressed
//   image.  The returned block should have all 4 pixel rows equal to
//   the bottom row of last_row_block.
//
// CornerPadBlockFunctor:
//       const BlockType func(const BlockType &last_block);
//   This function is used during calls to Pad() to create a block
//   used to pad below and to the right of the bottom-right corner of
//   a compressed image.  The returned block should have all 16 pixels
//   equal to the bottom-right pixel of last_block.
//

#ifndef IMAGE_COMPRESSION_INTERNAL_COMPRESSOR4X4_HELPER_H_
#define IMAGE_COMPRESSION_INTERNAL_COMPRESSOR4X4_HELPER_H_

#include <stddef.h>
#include <string.h>
#include <algorithm>
#include <string>
#include <vector>

#include "base/integral_types.h"
#include "base/logging.h"
#include "image_compression/internal/color_util.h"
#include "image_compression/internal/pixel4x4.h"
#include "image_compression/public/compressed_image.h"

namespace image_codec_compression {

class Compressor4x4Helper {
 public:
  // Given the number of pixels in the height or width of an image, this
  // returns the number of encoded blocks that are needed in that
  // dimension.
  static inline uint32 GetNumBlocks(uint32 num_pixels) {
    return (num_pixels + 3) / 4;
  }

  // Compresses an image.
  template <typename BlockType, typename ColorType, typename EncodeFunctor>
  static bool Compress(EncodeFunctor encode,
                       const std::string &compressor_name,
                       CompressedImage::Format format,
                       uint32 height, uint32 width,
                       uint32 padding_bytes_per_row,
                       const uint8 *buffer, CompressedImage *image);

  // Decompresses an image.
  template <typename BlockType, typename ColorType, typename DecodeFunctor>
  inline static bool Decompress(DecodeFunctor decode,
                                const CompressedImage &image,
                                std::vector<uint8> *decompressed_buffer);

  // Downsamples a compressed image. This returns false if the input
  // image is invalid or if there is an odd number of blocks in either
  // dimension.
  template <typename BlockType, typename ColorType,
            typename EncodeFunctor, typename DecodeFunctor>
  inline static bool Downsample(EncodeFunctor encode, DecodeFunctor decode,
                                const CompressedImage &image,
                                CompressedImage *downsampled_image);

  // Pads an image.
  template <typename BlockType, typename ColumnPadBlockFunctor,
            typename RowPadBlockFunctor, typename CornerPadBlockFunctor>
  inline static bool Pad(ColumnPadBlockFunctor get_column_pad_block,
                         RowPadBlockFunctor get_row_pad_block,
                         CornerPadBlockFunctor get_corner_pad_block,
                         const CompressedImage &image, uint32 padded_height,
                         uint32 padded_width, CompressedImage *padded_image);

  // Compresses and pads an image in a single operation. This is more
  // efficient than compressing and padding in separate calls.
  template <typename BlockType, typename ColorType, typename EncodeFunctor>
  inline static bool CompressAndPad(
      EncodeFunctor encode, const std::string &compressor_name,
      CompressedImage::Format format, uint32 height, uint32 width,
      uint32 padded_height, uint32 padded_width, uint32 padding_bytes_per_row,
      const uint8 *buffer, CompressedImage *padded_image);

  // Creates a solid image.
  template <typename BlockType>
  inline static bool CreateSolidImage(
      const std::string &compressor_name,
      CompressedImage::Format format, uint32 height, uint32 width,
      const BlockType &block, CompressedImage *image);

  // Copies a subimage.  Returns false if either starting index or
  // dimension is not a multiple of 4.
  template <typename BlockType>
  inline static bool CopySubimage(const CompressedImage &image,
                                  uint32 start_row, uint32 start_column,
                                  uint32 height, uint32 width,
                                  CompressedImage *subimage);

 private:
  // These are used by Downsample() to downsample a 2x2, 2x1, or 1x2 of blocks
  // into a single block.
  template <typename BlockType, typename ColorType,
            typename EncodeFunctor, typename DecodeFunctor>
  inline static const BlockType DownsampleBlocks2x2(
      EncodeFunctor encode, DecodeFunctor decode,
      const BlockType *blocks[2][2]);
  template <typename BlockType, typename ColorType,
            typename EncodeFunctor, typename DecodeFunctor>
  inline static const BlockType DownsampleBlocks2x1(
      EncodeFunctor encode, DecodeFunctor decode, const BlockType *blocks[2]);
  template <typename BlockType, typename ColorType,
            typename EncodeFunctor, typename DecodeFunctor>
  inline static const BlockType DownsampleBlocks1x2(
      EncodeFunctor encode, DecodeFunctor decode, const BlockType *blocks[2]);

  // Sets up the given CompressedImage to be used for an image with
  // the given info. This handles instances with internal or external
  // storage. Returns false on error.
  static bool SetUpCompressedImage(
      const std::string &compressor_name, size_t block_size,
      CompressedImage::Format format, uint32 height, uint32 width,
      uint32 padding_bytes_per_row, CompressedImage *image);
};

//-----------------------------------------------------------------------------

template <typename BlockType, typename ColorType, typename EncodeFunctor>
bool Compressor4x4Helper::Compress(
    EncodeFunctor encode, const std::string &compressor_name,
    CompressedImage::Format format, uint32 height, uint32 width,
    uint32 padding_bytes_per_row, const uint8 *buffer, CompressedImage *image) {
  DCHECK(buffer);
  DCHECK(image);
  if (!SetUpCompressedImage(compressor_name, sizeof(BlockType),
                            format, height, width,
                            padding_bytes_per_row, image)) {
    return false;
  }

  const ColorType *uncompressed_pixels =
      reinterpret_cast<const ColorType*>(buffer);
  BlockType *block = reinterpret_cast<BlockType*>(image->GetMutableData());
  const bool swap_red_and_blue = NeedsRedAndBlueSwapped(format);

  // TODO(user): Don't need to compress blocks after the first
  // row or col that is completely outside the original image (i.e.,
  //     (row > (uncompressed_height + 7) / 4
  //  or (col > (uncompressed_width  + 7) / 4)
  // Instead, just copy or memcpy blocks or rows of blocks.
  const uint32 num_block_rows = GetNumBlocks(height);
  const uint32 num_block_cols = GetNumBlocks(width);
  DCHECK_GT(num_block_rows, 0U);
  DCHECK_GT(num_block_cols, 0U);
  for (uint32 row = 0; row < num_block_rows; ++row) {
    for (uint32 col = 0; col < num_block_cols; ++col) {
      // Encode one 4x4 block of pixels. The Pixel4x4 constructor
      // handles the cases where the 4x4 window is no longer over
      // valid uncompressed data.
      DCHECK_LT(static_cast<const void*>(block),
                static_cast<const void*>(image->GetData() +
                                         image->GetDataSize()));
      *block++ = encode(Pixel4x4(uncompressed_pixels, height, width,
                                 padding_bytes_per_row, row * 4, col * 4),
                        swap_red_and_blue);
    }
  }
  return true;
}

template <typename BlockType, typename ColorType, typename DecodeFunctor>
bool Compressor4x4Helper::Decompress(DecodeFunctor decode,
                                     const CompressedImage &image,
                                     std::vector<uint8> *decompressed_buffer) {
  DCHECK(decompressed_buffer);

  const CompressedImage::Metadata &metadata = image.GetMetadata();
  decompressed_buffer->resize(metadata.uncompressed_height *
                              metadata.uncompressed_width * sizeof(ColorType));

  const BlockType *block = reinterpret_cast<const BlockType*>(image.GetData());
  ColorType *decompressed_pixels =
      reinterpret_cast<ColorType*>(&decompressed_buffer->at(0));

  const int height = static_cast<int>(metadata.uncompressed_height);
  const int width = static_cast<int>(metadata.uncompressed_width);

  const int num_block_rows = GetNumBlocks(metadata.uncompressed_height);
  const int num_block_cols = GetNumBlocks(metadata.uncompressed_width);
  const bool swap_red_and_blue = NeedsRedAndBlueSwapped(metadata.format);
  const int bytes_per_row = width * sizeof(decompressed_pixels[0]) +
                            metadata.padding_bytes_per_row;

  for (int block_row = 0; block_row < num_block_rows; ++block_row) {
    for (int block_col = 0; block_col < num_block_cols; ++block_col) {
      // Decode one 4x4 block of pixel colors.
      ColorType block_pixels[4][4];
      decode(*block++, swap_red_and_blue, block_pixels);

      // Store the pixels in the destination buffer, checking first
      // if the pixels are outside the destination image area.
      // TODO(user): Use separate loops if this is too slow.
      const int image_row = 4 * block_row;
      const int image_col = 4 * block_col;
      const int rows_inside_image = std::min(4, height - image_row);
      const int cols_inside_image = std::min(4, width - image_col);
      if (rows_inside_image > 0 && cols_inside_image > 0)
        StoreColorsInBuffer<ColorType>(block_pixels, bytes_per_row,
                                       rows_inside_image, cols_inside_image,
                                       image_row, image_col,
                                       decompressed_pixels);
    }
  }
  return true;
}

template <typename BlockType, typename ColorType,
          typename EncodeFunctor, typename DecodeFunctor>
bool Compressor4x4Helper::Downsample(EncodeFunctor encode, DecodeFunctor decode,
                                     const CompressedImage &image,
                                     CompressedImage *downsampled_image) {
  DCHECK(downsampled_image);
  // TODO(user): As with compression, can avoid reducing blocks outside
  // the original uncompressed image - just downsample once and copy
  // across rows and columns.

  // Downsampling requires an even number of blocks in both dimensions, except
  // for the special case of 1: if there is 1 block in a dimension, the
  // downsampled version will also have 1 block in that dimension, meaning the
  // minimum downsampled image size is 4 in each dimension.
  const CompressedImage::Metadata &metadata = image.GetMetadata();
  const int num_orig_block_rows = GetNumBlocks(metadata.uncompressed_height);
  const int num_orig_block_cols = GetNumBlocks(metadata.uncompressed_width);
  if ((num_orig_block_rows > 1 && num_orig_block_rows % 2 != 0) ||
      (num_orig_block_cols > 1 && num_orig_block_cols % 2 != 0)) {
    return false;
  }

  const uint32 orig_height = metadata.uncompressed_height;
  const uint32 orig_width = metadata.uncompressed_width;
  const uint32 downsampled_height = (orig_height + 1) / 2;
  const uint32 downsampled_width =  (orig_width + 1) / 2;
  if (!SetUpCompressedImage(metadata.compressor_name, sizeof(BlockType),
                            metadata.format, downsampled_height,
                            downsampled_width, 0, downsampled_image)) {
    return false;
  }

  const BlockType *orig_blocks =
      reinterpret_cast<const BlockType*>(image.GetData());
  BlockType *downsampled_block =
      reinterpret_cast<BlockType*>(downsampled_image->GetMutableData());
  const int num_downsampled_block_rows = num_orig_block_rows / 2;
  const int num_downsampled_block_cols = num_orig_block_cols / 2;

  // There are four possible cases, depending on whether there are multiple
  // rows and columns of blocks to downsample.
  if (num_orig_block_rows > 1 && num_orig_block_cols > 1) {
    // The common case of at least 2 blocks in each dimension.
    const BlockType *blocks_to_downsample[2][2];
    for (int row = 0; row < num_downsampled_block_rows; ++row) {
      for (int col = 0; col < num_downsampled_block_cols; ++col) {
        const int r0 = 2 * (num_orig_block_cols * row + col);
        const int r1 = r0 + num_orig_block_cols;
        blocks_to_downsample[0][0] = &orig_blocks[r0];
        blocks_to_downsample[0][1] = &orig_blocks[r0 + 1];
        blocks_to_downsample[1][0] = &orig_blocks[r1];
        blocks_to_downsample[1][1] = &orig_blocks[r1 + 1];
        *downsampled_block++ =
            DownsampleBlocks2x2<BlockType, ColorType,
                                EncodeFunctor, DecodeFunctor>(
                encode, decode, blocks_to_downsample);
      }
    }
  } else if (num_orig_block_rows > 1) {
    // At least 2 rows of blocks, only one column of blocks.
    const BlockType *blocks_to_downsample[2];
    for (int row = 0; row < num_downsampled_block_rows; ++row) {
      blocks_to_downsample[0] = &orig_blocks[2 * row];
      blocks_to_downsample[1] = &orig_blocks[2 * row + 1];
      *downsampled_block++ =
            DownsampleBlocks2x1<BlockType, ColorType,
                                EncodeFunctor, DecodeFunctor>(
                encode, decode, blocks_to_downsample);
    }
  } else if (num_orig_block_cols > 1) {
    // At least 2 columns of blocks, only one row of blocks.
    const BlockType *blocks_to_downsample[2];
    for (int col = 0; col < num_downsampled_block_cols; ++col) {
      blocks_to_downsample[0] = &orig_blocks[2 * col];
      blocks_to_downsample[1] = &orig_blocks[2 * col + 1];
      *downsampled_block++ =
          DownsampleBlocks1x2<BlockType, ColorType,
                              EncodeFunctor, DecodeFunctor>(
                encode, decode, blocks_to_downsample);
    }
  } else {
    // Only 1 row and 1 column.  This means each dimension of the original
    // block must be 4, 2, or 1.
    if (orig_height == 3 || orig_width == 3)
      return false;

    ColorType decoded_pixels[4][4];
    Pixel4x4 pixel4x4;
    decode(orig_blocks[0], false, decoded_pixels);

    // Replicate pixels if necessary to fill the full 4x4 array so that the
    // standard 4x4 block downsampler can be used.
    if (orig_width == 1) {
      for (int row = 0; row < 4; ++row) {
        decoded_pixels[row][1] = decoded_pixels[row][0];
        decoded_pixels[row][2] = decoded_pixels[row][0];
        decoded_pixels[row][3] = decoded_pixels[row][0];
      }
    } else if (orig_width == 2) {
      for (int row = 0; row < 4; ++row) {
        decoded_pixels[row][2] = decoded_pixels[row][0];
        decoded_pixels[row][3] = decoded_pixels[row][1];
      }
    }
    if (orig_height == 1) {
      for (int col = 0; col < 4; ++col) {
        decoded_pixels[1][col] = decoded_pixels[0][col];
        decoded_pixels[2][col] = decoded_pixels[0][col];
        decoded_pixels[3][col] = decoded_pixels[0][col];
      }
    } else if (orig_height == 2) {
      for (int col = 0; col < 4; ++col) {
        decoded_pixels[2][col] = decoded_pixels[0][col];
        decoded_pixels[3][col] = decoded_pixels[1][col];
      }
    }

    // Now downsample it.
    for (int row = 0; row < 2; ++row) {
      for (int col = 0; col < 2; ++col) {
        StoreDownsampledPixels4x4<ColorType>(decoded_pixels, 2 * row, 2 * col,
                                             &pixel4x4);
      }
    }
    *downsampled_block = encode(pixel4x4, false);
  }
  return true;
}

template <typename BlockType, typename ColumnPadBlockFunctor,
          typename RowPadBlockFunctor, typename CornerPadBlockFunctor>
bool Compressor4x4Helper::Pad(ColumnPadBlockFunctor get_column_pad_block,
                              RowPadBlockFunctor get_row_pad_block,
                              CornerPadBlockFunctor get_corner_pad_block,
                              const CompressedImage &image,
                              uint32 padded_height, uint32 padded_width,
                              CompressedImage *padded_image) {
  DCHECK(padded_image);
  // No padding to do if the sizes are not larger. Just copy.
  const CompressedImage::Metadata &metadata = image.GetMetadata();
  if (metadata.compressed_height >= padded_height &&
      metadata.compressed_width >= padded_width) {
    padded_image->Duplicate(image);
    return true;
  }

  if (!SetUpCompressedImage(metadata.compressor_name, sizeof(BlockType),
                            metadata.format, padded_height, padded_width,
                            0, padded_image)) {
    return false;
  }
  const int num_orig_block_rows = GetNumBlocks(metadata.compressed_height);
  const int num_orig_block_cols = GetNumBlocks(metadata.compressed_width);
  const int num_padded_block_rows = GetNumBlocks(padded_height);
  const int num_padded_block_cols = GetNumBlocks(padded_width);

  // Operate on blocks to save time.
  const BlockType *orig_block_start =
      reinterpret_cast<const BlockType*>(image.GetData());
  BlockType *padded_block_start =
      reinterpret_cast<BlockType*>(padded_image->GetMutableData());

  // Copy unpadded blocks.
  const BlockType *orig_block = orig_block_start;
  BlockType *padded_block = padded_block_start;
  const int num_orig_row_bytes = num_orig_block_cols * sizeof(BlockType);
  for (int row = 0; row < num_orig_block_rows; ++row) {
    memcpy(padded_block, orig_block, num_orig_row_bytes);

    if (num_orig_block_cols < num_padded_block_cols) {
      // Create a block to replicate the last column and copy it in
      // the space to the right of the orig image.
      const BlockType *last_col_block = padded_block + num_orig_block_cols - 1;
      BlockType pad_block(get_column_pad_block(*last_col_block));
      for (int col = num_orig_block_cols; col < num_padded_block_cols; ++col)
        padded_block[col] = pad_block;
    }

    orig_block += num_orig_block_cols;
    padded_block += num_padded_block_cols;
  }

  if (num_orig_block_rows < num_padded_block_rows) {
    // Create a row of blocks to replicate the last row (including
    // padding on the right) and copy it into the space below and to
    // the right of the original image.
    const BlockType *last_block_row =
        orig_block_start + (num_orig_block_rows - 1) * num_orig_block_cols;

    std::vector<BlockType> last_padded_row_blocks;
    last_padded_row_blocks.reserve(num_padded_block_cols);
    // Create row pad blocks for the left part of the row.
    for (int col = 0; col < num_orig_block_cols; ++col)
      last_padded_row_blocks.push_back(get_row_pad_block(last_block_row[col]));
    // Create a corner pad block for the rest of the row.
    if (num_orig_block_cols < num_padded_block_cols) {
      last_padded_row_blocks.insert(
          last_padded_row_blocks.end(),
          num_padded_block_cols - num_orig_block_cols,
          get_corner_pad_block(*(last_block_row + num_orig_block_cols - 1)));
    }
    // Copy the row into the padded buffer.
    BlockType *padded_block_row =
        padded_block_start + num_orig_block_rows * num_padded_block_cols;
    const int num_padded_row_bytes = num_padded_block_cols * sizeof(BlockType);
    for (int row = num_orig_block_rows; row < num_padded_block_rows; ++row) {
      memcpy(padded_block_row, &last_padded_row_blocks[0],
             num_padded_row_bytes);
      padded_block_row += num_padded_block_cols;
    }
  }

  return true;
}

template <typename BlockType, typename ColorType, typename EncodeFunctor>
bool Compressor4x4Helper::CompressAndPad(
    EncodeFunctor encode, const std::string &compressor_name,
    CompressedImage::Format format, uint32 height, uint32 width,
    uint32 padded_height, uint32 padded_width, uint32 padding_bytes_per_row,
    const uint8 *buffer, CompressedImage *padded_image) {
  DCHECK(buffer);
  DCHECK(padded_image);
  const uint32 final_height = std::max(height, padded_height);
  const uint32 final_width = std::max(width, padded_width);

  if (!SetUpCompressedImage(compressor_name, sizeof(BlockType),
                            format, final_height, final_width,
                            padding_bytes_per_row, padded_image)) {
    return false;
  }

  const ColorType *uncompressed_pixels =
      reinterpret_cast<const ColorType*>(buffer);
  BlockType *block =
      reinterpret_cast<BlockType*>(padded_image->GetMutableData());
  const bool swap_red_and_blue = NeedsRedAndBlueSwapped(format);

  const uint32 num_block_rows = GetNumBlocks(final_height);
  const uint32 num_block_cols = GetNumBlocks(final_width);
  DCHECK_GT(num_block_rows, 0U);
  DCHECK_GT(num_block_cols, 0U);
  for (uint32 row = 0; row < num_block_rows; ++row) {
    for (uint32 col = 0; col < num_block_cols; ++col) {
      // Encode one 4x4 block of pixels. The Pixel4x4 constructor
      // handles the cases where the 4x4 window is no longer over
      // valid uncompressed data.
      DCHECK_LT(static_cast<const void*>(block),
                static_cast<const void*>(padded_image->GetData() +
                                         padded_image->GetDataSize()));
      *block++ = encode(Pixel4x4(uncompressed_pixels, height, width,
                                 padding_bytes_per_row, row * 4, col * 4),
                        swap_red_and_blue);
    }
  }
  return true;
}

template <typename BlockType>
bool Compressor4x4Helper::CreateSolidImage(
    const std::string &compressor_name,
    CompressedImage::Format format, uint32 height, uint32 width,
    const BlockType &block, CompressedImage *image) {
  if (!SetUpCompressedImage(compressor_name, sizeof(BlockType),
                            format, height, width, 0, image)) {
    return false;
  }

  const uint32 num_block_rows = GetNumBlocks(height);
  const uint32 num_block_cols = GetNumBlocks(width);
  const uint32 num_blocks = num_block_rows * num_block_cols;

  // Copy block repeatedly.
  BlockType *block_start =
      reinterpret_cast<BlockType*>(image->GetMutableData());
  for (uint32 i = 0; i < num_blocks; ++i)
    block_start[i] = block;

  return true;
}

template <typename BlockType>
bool Compressor4x4Helper::CopySubimage(const CompressedImage &image,
                                       uint32 start_row, uint32 start_column,
                                       uint32 height, uint32 width,
                                       CompressedImage *subimage) {
  DCHECK(subimage);

  // Verify that all values are multiples of 4 and the subregion is
  // fully contained within the original image.
  const CompressedImage::Metadata &metadata = image.GetMetadata();
  if (start_row % 4 != 0 || start_column % 4 != 0 ||
      height % 4 != 0 || width % 4 != 0 ||
      start_row > metadata.compressed_height ||
      start_column > metadata.compressed_width ||
      start_row + height > metadata.compressed_height ||
      start_column + width > metadata.compressed_width) {
    return false;
  }

  if (!SetUpCompressedImage(metadata.compressor_name, sizeof(BlockType),
                            metadata.format, height, width, 0, subimage)) {
    return false;
  }

  const BlockType* orig_block_start =
      reinterpret_cast<const BlockType*>(image.GetData());
  BlockType* subimage_block_start =
      reinterpret_cast<BlockType*>(subimage->GetMutableData());

  const int num_orig_block_cols = GetNumBlocks(metadata.compressed_width);
  const int orig_start_block_col = GetNumBlocks(start_column);
  const int orig_start_block_row = GetNumBlocks(start_row);
  const int num_subimage_block_rows = GetNumBlocks(height);
  const int num_subimage_block_cols = GetNumBlocks(width);

  const BlockType* orig_block =
      orig_block_start + orig_start_block_row * num_orig_block_cols +
      orig_start_block_col;
  BlockType* subimage_block = subimage_block_start;
  for (int row = 0; row < num_subimage_block_rows; ++row) {
    memcpy(subimage_block, orig_block,
           num_subimage_block_cols * sizeof(*subimage_block));
    orig_block += num_orig_block_cols;
    subimage_block += num_subimage_block_cols;
  }

  return true;
}

template <typename BlockType, typename ColorType,
          typename EncodeFunctor, typename DecodeFunctor>
const BlockType Compressor4x4Helper::DownsampleBlocks2x2(
    EncodeFunctor encode, DecodeFunctor decode, const BlockType *blocks[2][2]) {
  ColorType decoded_pixels[4][4];
  Pixel4x4 pixel4x4;
  for (int row = 0; row < 2; ++row) {
    for (int col = 0; col < 2; ++col) {
      decode(*blocks[row][col], false, decoded_pixels);
      StoreDownsampledPixels4x4<ColorType>(decoded_pixels, 2 * row, 2 * col,
                                           &pixel4x4);
    }
  }
  return encode(pixel4x4, false);
}

template <typename BlockType, typename ColorType,
          typename EncodeFunctor, typename DecodeFunctor>
const BlockType Compressor4x4Helper::DownsampleBlocks2x1(
    EncodeFunctor encode, DecodeFunctor decode, const BlockType *blocks[2]) {
  ColorType decoded_pixels[4][4];
  Pixel4x4 pixel4x4;
  for (int row = 0; row < 2; ++row) {
    decode(*blocks[row], false, decoded_pixels);
    StoreDownsampledPixels4x4<ColorType>(decoded_pixels, 2 * row, 0, &pixel4x4);
    StoreDownsampledPixels4x4<ColorType>(decoded_pixels, 2 * row, 2, &pixel4x4);
  }
  return encode(pixel4x4, false);
}

template <typename BlockType, typename ColorType,
          typename EncodeFunctor, typename DecodeFunctor>
const BlockType Compressor4x4Helper::DownsampleBlocks1x2(
    EncodeFunctor encode, DecodeFunctor decode, const BlockType *blocks[2]) {
  ColorType decoded_pixels[4][4];
  Pixel4x4 pixel4x4;
  for (int col = 0; col < 2; ++col) {
    decode(*blocks[col], false, decoded_pixels);
    StoreDownsampledPixels4x4<ColorType>(decoded_pixels, 0, 2 * col, &pixel4x4);
    StoreDownsampledPixels4x4<ColorType>(decoded_pixels, 2, 2 * col, &pixel4x4);
  }
  return encode(pixel4x4, false);
}

}  // namespace image_codec_compression

#endif  // IMAGE_COMPRESSION_INTERNAL_COMPRESSOR4X4_HELPER_H_
