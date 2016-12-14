# Copyright 2015 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Top-level gyp configuration for Image Compression.
#
# Projects that use Image Compression should depend on the target defined here.
{
  'targets': [
    {
      'target_name': 'imagecompression',
      'type': 'static_library',
      'sources': [
        'image_compression/internal/bit_util.h',
        'image_compression/internal/color_types.h',
        'image_compression/internal/color_util.h',
        'image_compression/internal/compressor4x4_helper.cc',
        'image_compression/internal/compressor4x4_helper.h',
        'image_compression/internal/dxtc_compressor.cc',
        'image_compression/internal/dxtc_const_color_table.cc',
        'image_compression/internal/dxtc_const_color_table.h',
        'image_compression/internal/dxtc_to_etc_transcoder.cc',
        'image_compression/internal/etc_compressor.cc',
        'image_compression/internal/pixel4x4.cc',
        'image_compression/internal/pixel4x4.h',
        'image_compression/internal/pvrtc_compressor.cc',
        'image_compression/public/compressed_image.h',
        'image_compression/public/compressor.h',
        'image_compression/public/dxtc_compressor.h',
        'image_compression/public/dxtc_to_etc_transcoder.h',
        'image_compression/public/etc_compressor.h',
        'image_compression/public/pvrtc_compressor.h',
      ],
      'include_dirs': ['.'],
      'defines': [
        'IS_LITTLE_ENDIAN',
      ],
      'direct_dependent_settings': {
        'defines': [
          'IS_LITTLE_ENDIAN',
        ],
      },
    },
  ],
}
