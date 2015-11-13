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

#ifndef IMAGE_COMPRESSION_INTERNAL_BIT_UTIL_H_
#define IMAGE_COMPRESSION_INTERNAL_BIT_UTIL_H_

#include "base/integral_types.h"
#include "base/logging.h"

namespace image_codec_compression {

//
// This file contains types and functions that help deal with bit
// operations in image compression and decompression functions.
//

// Returns a 32-bit mask containing num_ones 1's in the least
// significant bits.
inline uint32 GetMask(int num_ones) {
  return (1 << num_ones) - 1;
}

// Generic 32-bit word bit accessor. Example: GetBits(n, 3, 2) returns
// a 2-bit integer defined by bits 3 and 4 of n. Bits are numbered
// from 0 (LSB) to 31 (MSB).
inline int GetBits(uint32 bits, int start_bit, int num_bits) {
  // Shift right to remove bits to the right of start_bit, then
  // mask off all but num_bits bits.
  return (bits >> start_bit) & GetMask(num_bits);
}

// Generic 3-bit word bit setter. Example: SetBits(3, 2, 1, &n) stores
// the 2-bit integer 1 into bits 3 and 4 of n. Bits are numbered from
// 0 (LSB) to 31 (MSB).
inline void SetBits(int start_bit, int num_bits, int value, uint32 *bits) {
  DCHECK_LT(value, 1 << num_bits);

  // Negative numbers must be converted to unsigned so they do not
  // sign-extend into the rest of the target value. This cast takes
  // care of both positive and negative cases.
  const uint32 mask = GetMask(num_bits);
  const uint32 unsigned_value = static_cast<uint32>(value) & mask;

  // Clear any bits that are set, then set the new bits.
  *bits = (*bits & ~(mask << start_bit)) | (unsigned_value << start_bit);
}

// Converts a signed value represented in num_bits bits to a 32-bit
// integer by extending the sign bit.
inline int32 ExtendSignBit(int32 value, int num_bits) {
  DCHECK_GE(num_bits, 0);
  DCHECK_LT(num_bits, 32);
  DCHECK_LT(value, 1 << num_bits);
  // Shift left so that the sign bit is in the high-order bit, then
  // shift back, which should copy the sign bit all the way.
  const int num_to_shift = 32 - num_bits;
  return (value << num_to_shift) >> num_to_shift;
}

}  // namespace image_codec_compression

#endif  // IMAGE_COMPRESSION_INTERNAL_BIT_UTIL_H_
