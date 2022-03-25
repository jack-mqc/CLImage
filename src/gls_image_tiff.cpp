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
#include <variant>
#include <vector>
#include <map>

#include "gls_dng_lossless_jpeg.hpp"
#include "gls_auto_ptr.hpp"
#include "gls_tiff_metadata.hpp"

namespace gls {

static void readTiffImageData(TIFF *tif, int height, int tiff_bitspersample, int tiff_samplesperpixel,
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

void read_tiff_file(const std::string& filename, int pixel_channels, int pixel_bit_depth, tiff_metadata* metadata,
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
static void writeTiffImageData(TIFF *tif, int width, int height, int pixel_channels, int pixel_bit_depth,
                        std::function<T*(int row)> row_pointer) {
    uint32_t rowsperstrip = TIFFDefaultStripSize(tif, -1);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rowsperstrip);

    auto_ptr<T> tiffbuf((T*)_TIFFmalloc(TIFFStripSize(tif)),
                        [](T* tiffbuf) { _TIFFfree(tiffbuf); });

    if (tiffbuf) {
        for (int row = 0; (row < height); row += rowsperstrip) {
            uint32_t nrow = (row + rowsperstrip) > height ? nrow = height - row : rowsperstrip;
            tstrip_t strip = TIFFComputeStrip(tif, row, 0);
            tsize_t bi = 0;
            for (int y = 0; y < nrow; ++y) {
                for (int x = 0; x < width; ++x) {
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
}

template <typename T>
void write_tiff_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                     tiff_compression compression, tiff_metadata* metadata, std::function<T*(int row)> row_pointer) {
    auto_ptr<TIFF> tif(TIFFOpen(filename.c_str(), "w"),
                       [](TIFF *tif) { TIFFClose(tif); });
    if (tif) {
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, compression);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, pixel_channels > 2 ? PHOTOMETRIC_RGB : PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, pixel_bit_depth);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, pixel_channels);

        TIFFSetField(tif, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

        writeTiffImageData(tif, width, height, pixel_channels, pixel_bit_depth, row_pointer);
    } else {
        throw std::runtime_error("Couldn't write tiff file.");
    }
}

void read_dng_file(const std::string& filename, int pixel_channels, int pixel_bit_depth, gls::tiff_metadata* metadata,
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
        printf("compression: %d\n", compression);

        if (metadata) {
            getAllTIFFTags(tif, metadata);
        }

        try {
            const auto blackLevel = (*metadata)["BlackLevel"];
            std::cout << "BlackLevel: " << std::get<std::vector<float>>(blackLevel)[0] << std::endl;
        }
        catch (const std::bad_variant_access& ex) {
            std::cout << ex.what() << '\n';
        }

        auto allocation_successful = image_allocator(width, height);
        if (allocation_successful) {
            if (compression == COMPRESSION_JPEG) {
                // DNG data is losslessly compressed

                uint32_t* stripbytecounts = 0;
                TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &stripbytecounts);
                printf("stripbytecounts: %d\n", stripbytecounts[0]);

                // We only expect one single big strip, but you never know
                uint32_t stripsize = stripbytecounts[0];

                if (TIFFNumberOfStrips(tif) != 1) {
                    throw std::runtime_error("Only one TIFF strip expected for compressed DNG files.");
                }

                tdata_t tiffbuf = _TIFFmalloc(stripsize);
                for (int strip = 0; strip < TIFFNumberOfStrips(tif); strip++) {
                    if (stripbytecounts[strip] > stripsize) {
                        tiffbuf = _TIFFrealloc(tiffbuf, stripbytecounts[strip]);
                        stripsize = stripbytecounts[strip];
                    }
                    if (TIFFReadRawStrip(tif, strip, tiffbuf, stripbytecounts[strip]) < 0) {
                        throw std::runtime_error("Failed to read compressed TIFF strip.");
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

        if (metadata) {
            fetchExifMetaData(tif, metadata);
        }
    } else {
        throw std::runtime_error("Couldn't read tiff file.");
    }
}

void write_dng_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                    tiff_compression compression, tiff_metadata* metadata, std::function<uint16_t*(int row)> row_pointer) {
    if (compression != COMPRESSION_NONE &&
        compression != COMPRESSION_JPEG &&
        compression != COMPRESSION_ADOBE_DEFLATE) {
        throw std::runtime_error("Only lossles JPEG and ADOBE_DEFLATE compression schemes are supported for DNG files. (" + std::to_string(compression) + ")");
    }

    augment_libtiff_with_custom_tags();

    auto_ptr<TIFF> tif(TIFFOpen(filename.c_str(), "w"),
                       [](TIFF *tif) { TIFFClose(tif); });

    if (tif) {
        if (metadata) {
            try {
                const auto blackLevel = (*metadata)["BlackLevel"];
                std::cout << "BlackLevel: " << std::get<std::vector<float>>(blackLevel)[0] << std::endl;
            }
            catch (const std::bad_variant_access& ex) {
                std::cout << ex.what() << '\n';
            }
        }

        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);

        TIFFSetField(tif, TIFFTAG_DNGVERSION, "\01\04\00\00");
        TIFFSetField(tif, TIFFTAG_DNGBACKWARDVERSION, "\01\03\00\00");
        TIFFSetField(tif, TIFFTAG_SUBFILETYPE, 0);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, compression);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_CFA);
        TIFFSetField(tif, TIFFTAG_CFALAYOUT, 1); // Rectangular (or square) layout
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);

        TIFFSetField(tif, TIFFTAG_MAKE, "Glass");
        TIFFSetField(tif, TIFFTAG_UNIQUECAMERAMODEL, "Glass 1");

        if (metadata) {
            writeVectorMetadata<uint16_t>(tif, metadata, "CFARepeatPatternDim");
            writeVectorMetadata<uint8_t>(tif, metadata, "CFAPattern");

            writeVectorMetadata<float>(tif, metadata, "ColorMatrix1");
            writeVectorMetadata<float>(tif, metadata, "ColorMatrix2");
            writeVectorMetadata<float>(tif, metadata, "AsShotNeutral");

            writeScalarMetadata<uint16_t>(tif, metadata, "CalibrationIlluminant1");
            writeScalarMetadata<uint16_t>(tif, metadata, "CalibrationIlluminant2");

            writeVectorMetadata<uint16_t>(tif, metadata, "BlackLevelRepeatDim");
            writeVectorMetadata<float>(tif, metadata, "BlackLevel");
            writeVectorMetadata<uint32_t>(tif, metadata, "WhiteLevel");

            writeScalarMetadata<uint32_t>(tif, metadata, "BayerGreenSplit");
            writeScalarMetadata<float>(tif, metadata, "BaselineExposure");
        }

        if (compression == COMPRESSION_JPEG) {
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, height);

            std::vector<uint16_t> outputBuffer(height * width);
            dng_stream out_stream((uint8_t*) outputBuffer.data(), outputBuffer.size() * sizeof(uint16_t));

            EncodeLosslessJPEG(row_pointer(0), height, width,
                               /*srcChannels=*/ 1, /*srcBitDepth=*/ 16, // TODO: reflect the actual bit depth
                               /*srcRowStep=*/ width, /*srcColStep=*/ 1, out_stream);

            if (TIFFWriteRawStrip(tif, 0, outputBuffer.data(), out_stream.Position()) < 0) {
                throw std::runtime_error("Failed to write TIFF data.");
            }
            printf("Wrote %ld compressed bytes\n", out_stream.Position());
        } else {
            writeTiffImageData(tif, width, height, pixel_channels, pixel_bit_depth, row_pointer);
        }

        TIFFWriteDirectory(tif);
    } else {
        throw std::runtime_error("Couldn't open DNG file for writing.");
    }
}

template
void write_tiff_file<uint8_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                              tiff_compression compression, tiff_metadata* metadata, std::function<uint8_t*(int row)> row_pointer);

template
void write_tiff_file<uint16_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                               tiff_compression compression, tiff_metadata* metadata, std::function<uint16_t*(int row)> row_pointer);

}  // namespace gls
