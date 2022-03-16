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

#include "gls_image_tiff.h"

#include <assert.h>

#include <tiffio.h>
#include <math.h>
#include <float.h>
#include <span>

#include <sys/stat.h>
#include <sys/time.h>

#include <iomanip>
#include <iostream>

#include "auto_ptr.hpp"

namespace gls {

void read_tiff_file(const std::string& filename, int pixel_channels, int pixel_bit_depth,
                    std::function<bool(int width, int height)> image_allocator,
                    std::function<void(int tiff_bitspersample, int tiff_samplesperpixel, int row, int strip_height,
                                       uint8_t *tiff_buffer)> process_tiff_strip) {
    auto_ptr<TIFF> tif(TIFFOpen(filename.c_str(), "r"),
                       [](TIFF *tif) { TIFFClose(tif); });

    if (tif) {
        uint32_t width, height;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

        uint16_t tiff_samplesperpixel;
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &tiff_samplesperpixel);

        uint16_t tiff_sampleformat;
        TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &tiff_sampleformat);

        uint16_t tiff_bitspersample=8;
        TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &tiff_bitspersample);

        if (tiff_sampleformat != SAMPLEFORMAT_UINT) {
            throw std::runtime_error("can not read sample format other than uint");
        }

        if ((tiff_bitspersample != 8 && tiff_bitspersample != 16)) {
            throw std::runtime_error("can not read sample with " + std::to_string(tiff_bitspersample) + " bits depth");
        }

        auto allocation_successful = image_allocator(width, height);
        if (allocation_successful) {
            auto_ptr<uint8_t> tiffbuf((uint8_t*)_TIFFmalloc(TIFFStripSize(tif)),
                                      [](uint8_t* tiffbuf) { _TIFFfree(tiffbuf); });
            if (tiffbuf) {
                uint32_t rowsperstrip = 0;
                TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

                for (uint32_t row = 0; row < height; row += rowsperstrip) {
                    uint32_t nrow = (row + rowsperstrip > height) ? (height - row) : rowsperstrip;
                    tstrip_t strip = TIFFComputeStrip(tif, row, 0);

                    if ((TIFFReadEncodedStrip(tif, strip, tiffbuf, -1)) < 0) {
                        throw std::runtime_error("Failed to encode TIFF strip.");
                    }
                    process_tiff_strip(tiff_bitspersample, tiff_samplesperpixel, row, nrow, tiffbuf);
                }
            } else {
                throw std::runtime_error("Error allocating memory buffer for TIFF strip.");
            }
        } else {
            throw std::runtime_error("Couldn't allocate image storage");
        }
    } else {
        throw std::runtime_error("Couldn't read tiff file.");
    }
}

template <typename T>
void write_tiff_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                     tiff_compression compression, std::function<T*(int row)> row_pointer) {
    auto_ptr<TIFF> tif(TIFFOpen(filename.c_str(), "w"),
                       [](TIFF *tif) { TIFFClose(tif); });
    if (tif) {
        uint16_t frame_width = width;
        uint16_t frame_height = height;
        uint32_t rowsperstrip = (uint32_t)-1;

        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, frame_width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, frame_height);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, compression);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, pixel_channels > 2 ? PHOTOMETRIC_RGB : PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, pixel_bit_depth);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, pixel_channels);
        rowsperstrip = TIFFDefaultStripSize(tif, rowsperstrip);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rowsperstrip);
        TIFFSetField(tif, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

        auto_ptr<T> tiffbuf((T*)_TIFFmalloc(TIFFStripSize(tif)),
                            [](T* tiffbuf) { _TIFFfree(tiffbuf); });

        if (tiffbuf) {
            for (int row = 0; (row < frame_height); row += rowsperstrip) {
                uint32_t nrow = (row + rowsperstrip) > frame_height ? nrow = frame_height - row : rowsperstrip;
                tstrip_t strip = TIFFComputeStrip(tif, row, 0);
                tsize_t bi = 0;
                for (int y = 0; y < nrow; ++y) {
                    for (int x = 0; x < frame_width; ++x) {
                        for (int c = 0; c < pixel_channels; c++) {
                            tiffbuf[bi++] = row_pointer(row + y)[pixel_channels * x + c];
                        }
                    }
                }
                if (TIFFWriteEncodedStrip(tif, strip, tiffbuf, bi * sizeof(T)) < 0) {
                    throw std::runtime_error("Failed to encode TIFF strip.");
                }
            }
        } else {
            throw std::runtime_error("Error allocating memory buffer for TIFF strip.");
        }
    } else {
        throw std::runtime_error("Couldn't write tiff file.");
    }
}

template <typename T>
void write_dng_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                     tiff_compression compression, std::function<T*(int row)> row_pointer) {
    auto_ptr<TIFF> tif(TIFFOpen(filename.c_str(), "w"),
                       [](TIFF *tif) { TIFFClose(tif); });
    if (tif) {
        const std::array<float, 9> color_matrix = {
            2.005, -0.771, -0.269,
            -0.752, 1.688, 0.064,
            -0.149, 0.283, 0.745
        };

        const std::array<float, 3> as_shot_neutral = {
            0.807133, 1.0, 0.913289
        };

//        const auto base_black_level = static_cast<float>(raw_color.black);
//        std::array<float, 4> black_level = {
//            base_black_level + static_cast<float>(raw_color.cblack[0]),
//            base_black_level + static_cast<float>(raw_color.cblack[1]),
//            base_black_level + static_cast<float>(raw_color.cblack[2]),
//            base_black_level + static_cast<float>(raw_color.cblack[3])
//        };
//
//        if (raw_color.cblack[4] == 2 && raw_color.cblack[5] == 2) {
//            for (int x = 0; x < raw_color.cblack[4]; ++x) {
//                for (int y = 0; y < raw_color.cblack[5]; ++y) {
//                    const auto index = y * 2 + x;
//                    black_level[index] = raw_color.cblack[6 + index];
//                }
//            }
//        }

        static const uint32_t white_level = 0xffff; // raw_color.maximum;


        static const short bayerPatternDimensions[] = { 2, 2 };

        TIFFSetField(tif, TIFFTAG_DNGVERSION, "\01\01\00\00");
        TIFFSetField(tif, TIFFTAG_SUBFILETYPE, 0);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_CFA);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
        TIFFSetField(tif, TIFFTAG_CFAREPEATPATTERNDIM, bayerPatternDimensions);
        TIFFSetField(tif, TIFFTAG_CFAPATTERN, 4, "\01\02\00\01"); // GBRG
        TIFFSetField(tif, TIFFTAG_MAKE, "DNG");
        TIFFSetField(tif, TIFFTAG_UNIQUECAMERAMODEL, "DNG");
        TIFFSetField(tif, TIFFTAG_COLORMATRIX1, 9, &color_matrix );
        TIFFSetField(tif, TIFFTAG_ASSHOTNEUTRAL, 3, &as_shot_neutral );
        // TIFFSetField(tif, TIFFTAG_BLACKLEVEL, 4, &black_level);
        TIFFSetField(tif, TIFFTAG_WHITELEVEL, 1, &white_level);

        for (int i = 0; i < height; i++)
            TIFFWriteScanline(tif, row_pointer(i), i, 0);

        TIFFWriteDirectory(tif);
    } else {
        throw std::runtime_error("Couldn't write tiff file.");
    }
}

template
void write_tiff_file<uint8_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                              tiff_compression compression, std::function<uint8_t*(int row)> row_pointer);

template
void write_tiff_file<uint16_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                               tiff_compression compression, std::function<uint16_t*(int row)> row_pointer);

template
void write_dng_file<uint16_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                              tiff_compression compression, std::function<uint16_t*(int row)> row_pointer);

}  // namespace gls
