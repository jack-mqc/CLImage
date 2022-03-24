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

#include "gls_dng_lossless_jpeg.hpp"

#include "gls_auto_ptr.hpp"

namespace gls {

void readTiffImageData(TIFF *tif, int height, int tiff_bitspersample, int tiff_samplesperpixel,
                       std::function<void(int tiff_bitspersample, int tiff_samplesperpixel, int row, int strip_height,
                                          uint8_t *tiff_buffer)> process_tiff_strip) {
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
}

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
        if (tiff_sampleformat != SAMPLEFORMAT_UINT) {
            throw std::runtime_error("can not read sample format other than uint");
        }

        uint16_t tiff_bitspersample=8;
        TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &tiff_bitspersample);
        if ((tiff_bitspersample != 8 && tiff_bitspersample != 16)) {
            throw std::runtime_error("can not read sample with " + std::to_string(tiff_bitspersample) + " bits depth");
        }

        auto allocation_successful = image_allocator(width, height);
        if (allocation_successful) {
            readTiffImageData(tif, height, tiff_bitspersample, tiff_samplesperpixel, process_tiff_strip);
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

#define TIFFTAG_FORWARDMATRIX1 50964
#define TIFFTAG_FORWARDMATRIX2 50965
#define TIFFTAG_TIMECODES 51043
#define TIFFTAG_FRAMERATE 51044
#define TIFFTAG_REELNAME 51081

#define TIFFTAG_PROFILENAME 50936
#define TIFFTAG_PROFILELOOKTABLEDIMS 50981
#define TIFFTAG_PROFILELOOKTABLEDATA 50982
#define TIFFTAG_PROFILELOOKTABLEENCODING 51108
#define TIFFTAG_DEFAULTUSERCROP 51125

static const TIFFFieldInfo xtiffFieldInfo[] = {
    { TIFFTAG_FORWARDMATRIX1, -1, -1, TIFF_SRATIONAL, FIELD_CUSTOM, 1, 1, "ForwardMatrix1" },
    { TIFFTAG_FORWARDMATRIX2, -1, -1, TIFF_SRATIONAL, FIELD_CUSTOM, 1, 1, "ForwardMatrix2" },
    { TIFFTAG_TIMECODES, -1, -1, TIFF_BYTE, FIELD_CUSTOM, 1, 1, "TimeCodes" },
    { TIFFTAG_FRAMERATE, -1, -1, TIFF_SRATIONAL, FIELD_CUSTOM, 1, 1, "FrameRate" },
    { TIFFTAG_REELNAME, -1, -1, TIFF_ASCII, FIELD_CUSTOM, 1, 0, "ReelName" },

    { TIFFTAG_PROFILENAME, -1, -1, TIFF_ASCII, FIELD_CUSTOM, 1, 0, "ProfileName" },
    { TIFFTAG_PROFILELOOKTABLEDIMS, 3, -1, TIFF_LONG, FIELD_CUSTOM, 1, 0, "ProfileLookTableDims" },
    { TIFFTAG_PROFILELOOKTABLEDATA, -1, -1, TIFF_FLOAT, FIELD_CUSTOM, 1, 0, "ProfileLookTableData" },
    { TIFFTAG_PROFILELOOKTABLEENCODING, -1, -1, TIFF_LONG, FIELD_CUSTOM, 1, 0, "ProfileLookTableEncoding" },
    { TIFFTAG_DEFAULTUSERCROP, 4, -1, TIFF_RATIONAL, FIELD_CUSTOM, 1, 0, "DefaultUserCrop" }
};

static TIFFExtendProc parent_extender = NULL;  // In case we want a chain of extensions

static void registerCustomTIFFTags( TIFF *tif )
{
    // Install the extended Tag field info
    TIFFMergeFieldInfo( tif, xtiffFieldInfo, sizeof( xtiffFieldInfo ) / sizeof( xtiffFieldInfo[0] ) );

    if( parent_extender )
        (*parent_extender)(tif);
}

static void augment_libtiff_with_custom_tags( void )
{
    static int first_time = 1;
    if( !first_time )
        return;
    first_time = 0;
    parent_extender = TIFFSetTagExtender( registerCustomTIFFTags );
}

void read_dng_file(const std::string& filename, int pixel_channels, int pixel_bit_depth,
                    std::function<bool(int width, int height)> image_allocator,
                    std::function<void(int tiff_bitspersample, int tiff_samplesperpixel, int row, int strip_height,
                                       uint8_t *tiff_buffer)> process_tiff_strip) {
    augment_libtiff_with_custom_tags();

    auto_ptr<TIFF> tif(TIFFOpen(filename.c_str(), "r"),
                       [](TIFF *tif) { TIFFClose(tif); });

    if (tif) {
        uint32_t width = 0, height = 0;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
        printf("width: %d, height: %d\n", width, height);

        uint16_t tiff_samplesperpixel = 0;
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &tiff_samplesperpixel);
        printf("tiff_samplesperpixel: %d\n", tiff_samplesperpixel);

        uint16_t tiff_sampleformat = SAMPLEFORMAT_UINT;
        TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &tiff_sampleformat);
        if (tiff_sampleformat != SAMPLEFORMAT_UINT) {
            throw std::runtime_error("can not read sample format other than uint: " + std::to_string(tiff_sampleformat));
        }

        uint16_t tiff_bitspersample = 0;
        TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &tiff_bitspersample);
        if ((tiff_bitspersample != 8 && tiff_bitspersample != 16)) {
            throw std::runtime_error("can not read sample with " + std::to_string(tiff_bitspersample) + " bits depth");
        }
        printf("tiff_bitspersample: %d\n", tiff_bitspersample);

        uint16_t compression = 0;
        TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &compression);
        if ((compression != COMPRESSION_NONE && compression != COMPRESSION_JPEG)) {
            throw std::runtime_error("only lossles JPEG compression is supported for DNG. (" + std::to_string(compression) + ")");
        }
        printf("compression: %d\n", compression);

        auto allocation_successful = image_allocator(width, height);
        if (allocation_successful) {
            if (compression == COMPRESSION_JPEG) {
                // DNG data is losslessly compressed

                uint32_t* stripbytecounts = 0;
                TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &stripbytecounts);
                printf("stripbytecounts: %d\n", stripbytecounts[0]);

                // We only expect one single big strip, but you never know
                uint32_t stripsize = stripbytecounts[0];

                tdata_t tiffbuf = _TIFFmalloc(stripsize);
                for (int strip = 0; strip < TIFFNumberOfStrips(tif); strip++) {
                    if (stripbytecounts[strip] > stripsize) {
                        tiffbuf = _TIFFrealloc(tiffbuf, stripbytecounts[strip]);
                        stripsize = stripbytecounts[strip];
                    }
                    if (TIFFReadRawStrip(tif, strip, tiffbuf, stripbytecounts[strip]) < 0) {
                        throw std::runtime_error("Failed to encode TIFF strip.");
                    }

                    // Used Adobe's version of libjpeg lossless codec
                    dng_stream stream((uint8_t *) tiffbuf, stripsize);
                    dng_spooler spooler;
                    uint32_t decodedSize = width * height * sizeof(uint16_t);
                    DecodeLosslessJPEG(stream, spooler,
                                       decodedSize,
                                       decodedSize,
                                       false, stripsize);

                    process_tiff_strip(tiff_bitspersample, tiff_samplesperpixel, 0, height, (uint8_t *) spooler.data());
                }
                _TIFFfree(tiffbuf);
            } else {
                // No compreession, read as a plain TIFF file

                readTiffImageData(tif, height, tiff_bitspersample, tiff_samplesperpixel, process_tiff_strip);
            }
        }
    } else {
        throw std::runtime_error("Couldn't read tiff file.");
    }
}

void write_dng_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                     tiff_compression compression, std::function<uint16_t*(int row)> row_pointer) {
    augment_libtiff_with_custom_tags();

    auto_ptr<TIFF> tiff(TIFFOpen(filename.c_str(), "w"),
                       [](TIFF *tif) { TIFFClose(tif); });

    if (tiff) {
        static const short bayerPatternDimensions[] = { 2, 2 };
        static const unsigned char bayerPattern[] = "\00\01\01\02"; // "\01\02\00\01";
        std::array<float, 9> color_matrix = {
            // 1.9435, -0.8992, -0.1936, 0.1144, 0.8380, 0.0475, 0.0136, 0.1203, 0.3553

            0.8059082031, -0.2783203125, -0.060546875,
            -0.5290527344,  1.441650391,   0.05517578125,
            -0.09350585938, 0.3002929688,  0.63671875
        };
//        std::array<float, 3> camera_multipliers = {
//            1.354894, 1.000000, 1.920234
//        };
        const std::array<float, 3> as_shot_neutral = {
//            1.f / (camera_multipliers[0] / camera_multipliers[1]),
//            1.f,
//            1.f / (camera_multipliers[2] / camera_multipliers[1])

            0.4731977819, 1, 0.7781155015
        };
        const std::array<float, 4> black_level = {0, 0, 0, 0};
        static const uint32_t white_level = 15000; // 0x0fff;

        TIFFSetField(tiff, TIFFTAG_DNGVERSION, "\01\04\00\00");
        TIFFSetField(tiff, TIFFTAG_DNGBACKWARDVERSION, "\01\04\00\00");
        TIFFSetField(tiff, TIFFTAG_SUBFILETYPE, 0);
        TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tiff, TIFFTAG_ROWSPERSTRIP, 1);
        TIFFSetField(tiff, TIFFTAG_ORIENTATION, ORIENTATION_LEFTTOP);  // Mirrored data
        TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_CFA);
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
        TIFFSetField(tiff, TIFFTAG_CFAREPEATPATTERNDIM, &bayerPatternDimensions);
        TIFFSetField(tiff, TIFFTAG_CFAPATTERN, 4, bayerPattern);

        TIFFSetField(tiff, TIFFTAG_MAKE, "Glass");
        TIFFSetField(tiff, TIFFTAG_UNIQUECAMERAMODEL, "Glass 1");

        TIFFSetField(tiff, TIFFTAG_COLORMATRIX1, 9, &color_matrix);
        TIFFSetField(tiff, TIFFTAG_CALIBRATIONILLUMINANT1, 21); // D65
        TIFFSetField(tiff, TIFFTAG_ASSHOTNEUTRAL, 3, &as_shot_neutral);

        TIFFSetField(tiff, TIFFTAG_CFALAYOUT, 1); // Rectangular (or square) layout
        TIFFSetField(tiff, TIFFTAG_CFAPLANECOLOR, 3, "\00\01\02"); // RGB https://www.awaresystems.be/imaging/tiff/tifftags/cfaplanecolor.html
        TIFFSetField(tiff, TIFFTAG_BAYERGREENSPLIT, 0);

        TIFFSetField(tiff, TIFFTAG_BLACKLEVEL, 4, &black_level);
        TIFFSetField(tiff, TIFFTAG_WHITELEVEL, 1, &white_level);

        TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);

        for (int row = 0; row < height; row++) {
            TIFFWriteScanline(tiff, row_pointer(row), row, 0);
        }

        TIFFWriteDirectory(tiff);
    } else {
        throw std::runtime_error("Couldn't open DNG file for writing.");
    }
}

template
void write_tiff_file<uint8_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                              tiff_compression compression, std::function<uint8_t*(int row)> row_pointer);

template
void write_tiff_file<uint16_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                               tiff_compression compression, std::function<uint16_t*(int row)> row_pointer);

}  // namespace gls
