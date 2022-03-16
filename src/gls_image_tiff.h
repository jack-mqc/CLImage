// Copyright (c) 2021-2022 Glass Imaging Inc.
// Author: Fabio Riccardi <fabio@glass-imaging.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef gls_image_tiff_hpp
#define gls_image_tiff_hpp

#include <functional>
#include <string>
#include <span>

namespace gls {

typedef enum tiff_compression {
    NONE            = 1,        /* dump mode */
    LZW             = 5,        /* Lempel-Ziv & Welch */
    JPEG            = 7,        /* %JPEG DCT compression */
    PACKBITS        = 32773,    /* Macintosh RLE */
    DEFLATE         = 32946,    /* Deflate compression */
    ADOBE_DEFLATE   = 8,        /* Deflate compression, as recognized by Adobe */
} tiff_compression;

void read_tiff_file(const std::string& filename, int pixel_channels, int pixel_bit_depth,
                    std::function<bool(int width, int height)> image_allocator,
                    std::function<void(int tiff_bitspersample, int tiff_samplesperpixel, int row, int strip_height,
                                       uint8_t *tiff_buffer)> process_tiff_strip);

template <typename T>
void write_tiff_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                     tiff_compression compression, std::function<T*(int row)> row_pointer);

template <typename T>
void write_dng_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                    tiff_compression compression, std::function<T*(int row)> row_pointer);

}  // namespace gls

#endif /* gls_image_tiff_hpp */
