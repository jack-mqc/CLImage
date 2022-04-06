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

inline static uint16_t swapBytes(uint16_t in) {
    return ((in & 0xff) << 8) | (in >> 8);
}

int findGcd(int a, int b) {
  while ((a % b) > 0) {
    int R = a % b;
    a = b;
    b = R;
  }
  return b;
}

inline static void unpack12BitsInto16Bits(uint16_t *out, const uint16_t *in, size_t in_size) {
    for (int i = 0; i < in_size; i += 3) {
        uint16_t in0 = swapBytes(in[i]);
        uint16_t in1 = swapBytes(in[i+1]);
        uint16_t in2 = swapBytes(in[i+2]);

        *out++ = in0 >> 4;
        *out++ = ((in0 << 8) & 0xfff) | (in1 >> 8);
        *out++ = ((in1 << 4) & 0xfff) | (in2 >> 12);
        *out++ = in2 & 0xfff;
    }
}

inline static void unpack14BitsInto16Bits(uint16_t *out, const uint16_t *in, size_t in_size) {
    for (int i = 0; i < in_size; i += 7) {
        uint16_t in0 = swapBytes(in[i]);
        uint16_t in1 = swapBytes(in[i+1]);
        uint16_t in2 = swapBytes(in[i+2]);
        uint16_t in3 = swapBytes(in[i+3]);
        uint16_t in4 = swapBytes(in[i+4]);
        uint16_t in5 = swapBytes(in[i+5]);
        uint16_t in6 = swapBytes(in[i+6]);

        *out++ = in0 >> 2;
        *out++ = ((in0 << 12) & 0x3fff) | (in1 >> 4);
        *out++ = ((in1 << 10) & 0x3fff) | (in2 >> 6);
        *out++ = ((in2 <<  8) & 0x3fff) | (in3 >> 8);
        *out++ = ((in3 <<  6) & 0x3fff) | (in4 >> 10);
        *out++ = ((in4 <<  4) & 0x3fff) | (in5 >> 12);
        *out++ = ((in5 <<  2) & 0x3fff) | (in6 >> 14);
        *out++ = (in6 & 0x3fff);
    }
}

// TODO: finish building a general unpacking function
/*
inline static void unpackTo16Bits(uint16_t *out, const uint16_t *in, int bitspersample, size_t in_size) {
    int gcd = findGcd(16, bitspersample);

    std::cout << "gcd: " << gcd << std::endl;

    int inputWords = bitspersample / gcd;
    int outputWords = 16 / gcd;

    std::cout << "inputWords: " << inputWords << ", outputWords: " << outputWords << std::endl;

    uint16_t inputBuffer[inputWords];
    uint16_t outputBuffer[outputWords];

    for (int i = 0; i < in_size; i += inputWords) {
        for (int j = 0; j < inputWords; j++) {
            inputBuffer[j] = swapBytes(in[i + j]);
        }

        int shift = 16 - bitspersample;
        for (int k = 0; k < outputWords; k++) {

        }
    }
}
*/

static void readTiffImageData(TIFF *tif, int width, int height, int tiff_bitspersample, int tiff_samplesperpixel,
                       std::function<void(int tiff_bitspersample, int tiff_samplesperpixel, int row, int strip_height,
                                          uint8_t *tiff_buffer)> process_tiff_strip) {
    size_t stripSize = TIFFStripSize(tif);
    auto_ptr<uint8_t> tiffbuf((uint8_t*)_TIFFmalloc(stripSize),
                              [](uint8_t* tiffbuf) { _TIFFfree(tiffbuf); });

    auto_ptr<uint8_t> decodedBuffer = tiff_bitspersample != 16
                                    ? auto_ptr<uint8_t>((uint8_t*)_TIFFmalloc(16 * stripSize / tiff_bitspersample),
                                                         [](uint8_t* buffer) { _TIFFfree(buffer); })
                                    : auto_ptr<uint8_t>(nullptr, [](uint8_t* buffer){ });

    printf("stripSize: %ld, width: %d\n", stripSize, width);

    if (tiffbuf) {
        uint32_t rowsperstrip = 0;
        TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

        for (uint32_t row = 0; row < height; row += rowsperstrip) {
            uint32_t nrow = (row + rowsperstrip > height) ? (height - row) : rowsperstrip;
            tstrip_t strip = TIFFComputeStrip(tif, row, 0);

            if ((TIFFReadEncodedStrip(tif, strip, tiffbuf, -1)) < 0) {
                throw std::runtime_error("Failed to encode TIFF strip.");
            }

            if (tiff_bitspersample == 12) {
                unpack12BitsInto16Bits((uint16_t*) decodedBuffer.get(), (uint16_t*) tiffbuf.get(), stripSize / sizeof(uint16_t));

                process_tiff_strip(/* tiff_bitspersample=*/ 16, tiff_samplesperpixel, row, nrow, decodedBuffer);
            } else if (tiff_bitspersample == 14) {
                unpack14BitsInto16Bits((uint16_t*) decodedBuffer.get(), (uint16_t*) tiffbuf.get(), stripSize / sizeof(uint16_t));

                process_tiff_strip(/* tiff_bitspersample=*/ 16, tiff_samplesperpixel, row, nrow, decodedBuffer);
            } else if (tiff_bitspersample == 16) {
                process_tiff_strip(/* tiff_bitspersample=*/ 16, tiff_samplesperpixel, row, nrow, tiffbuf);
            } else {
                throw std::runtime_error("tiff_bitspersample " + std::to_string(tiff_bitspersample) + " not supported.");
            }
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

        uint16_t tiff_bitspersample;
        TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &tiff_bitspersample);
        if ((tiff_bitspersample != 8 && tiff_bitspersample != 16)) {
            throw std::runtime_error("can not read sample with " + std::to_string(tiff_bitspersample) + " bits depth");
        }

        auto allocation_successful = image_allocator(width, height);
        if (allocation_successful) {
            readTiffImageData(tif, width, height, tiff_bitspersample, tiff_samplesperpixel, process_tiff_strip);
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
        if (metadata) {
            getAllTIFFTags(tif, metadata);
        }

        uint32_t subfileType = 0;
        TIFFGetField(tif, TIFFTAG_SUBFILETYPE, &subfileType);

        if (subfileType & 1) {
            // This is a preview, look for the real image

            uint16_t subifdCount;
            uint64_t* subIFD;
            TIFFGetField(tif, TIFFTAG_SUBIFD, &subifdCount, &subIFD);

            uint16_t numberOfDirectories = TIFFNumberOfDirectories(tif);
#if defined(__ANDROID__)
            printf("numberOfDirectories: %d, subfileType: %d, subifdCount: %d, subIFD: %ld\n",
                   numberOfDirectories, subfileType, subifdCount, subIFD[0]);
#else
            printf("numberOfDirectories: %d, subfileType: %d, subifdCount: %d, subIFD: %lld\n",
                   numberOfDirectories, subfileType, subifdCount, subIFD[0]);
#endif
            for (int i = 0; i < subifdCount; i++) {
                TIFFSetSubDirectory(tif, subIFD[i]);

                uint32_t subfileType = 0;
                TIFFGetField(tif, TIFFTAG_SUBFILETYPE, &subfileType);

                printf("subfile %i, subfileType: %d\n", i, subfileType);

                if ((subfileType & 1) == 0) {
                    if (metadata) {
                        getAllTIFFTags(tif, metadata);
                    }
                    break;
                }
            }
        }

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
        printf("tiff_bitspersample: %d\n", tiff_bitspersample);

        uint16_t compression = 0;
        TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &compression);
        printf("compression: %d\n", compression);

        uint32_t width = 0, height = 0;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
        printf("width: %d, height: %d\n", width, height);

        auto allocation_successful = image_allocator(width, height);
        if (allocation_successful) {
            printf("TIFFIsTiled: %d\n", TIFFIsTiled(tif));

            if (TIFFIsTiled(tif)) {
                uint32_t maxTileWidth, maxTileHeight;

                TIFFGetField(tif, TIFFTAG_TILEWIDTH, &maxTileWidth);
                TIFFGetField(tif, TIFFTAG_TILELENGTH, &maxTileHeight);

                printf("tileWidth: %d, tileHeight: %d\n", maxTileWidth, maxTileHeight);

                tmsize_t tileSize = TIFFTileSize(tif);
                uint32_t tileCount = TIFFNumberOfTiles(tif);
                uint32_t tileCountX = (width + maxTileWidth - 1) / maxTileWidth;
                auto_ptr<uint16_t> tiffbuf((uint16_t*)_TIFFmalloc(tileSize),
                                           [](uint16_t* buf) { _TIFFfree(buf); });

                auto_ptr<uint16_t> imagebuf((uint16_t*)_TIFFmalloc(width * height * sizeof(uint16_t)),
                                            [](uint16_t* buf) { _TIFFfree(buf); });

                if (compression == COMPRESSION_JPEG) {
                    for (uint32_t tile = 0; tile < tileCount; tile++) {
                        uint32_t tileX = maxTileWidth * (tile % tileCountX);
                        uint32_t tileY = maxTileHeight * (tile / tileCountX);

                        uint32_t tileWidth = std::min(tileX + maxTileWidth, width) - tileX;
                        uint32_t tileHeight = std::min(tileY + maxTileHeight, height) - tileY;

                        tmsize_t tileBytes = TIFFReadRawTile(tif, tile, (uint8_t *) tiffbuf.get(), (tsize_t) -1);
                        if (tileBytes < 0) {
                            throw std::runtime_error("Failed to read TIFF tile " + std::to_string(tile));
                        } else {
                            // Used Adobe's version of libjpeg lossless codec
                            dng_stream stream((uint8_t *) tiffbuf.get(), tileSize);
                            dng_spooler spooler;
                            uint32_t decodedSize = maxTileWidth * maxTileHeight * sizeof(uint16_t);
                            DecodeLosslessJPEG(stream, spooler,
                                               decodedSize,
                                               decodedSize,
                                               false, tileSize);

                            uint16_t *tilePixels = (uint16_t *) spooler.data();
                            for (int y = 0; y < tileHeight; y++) {
                                for (int x = 0; x < tileWidth; x++) {
                                    imagebuf[(tileY + y) * width + tileX + x] = tilePixels[y * tileWidth + x];
                                }
                            }
                        }
                    }

                    // The output of the JPEG decoder is always 16 bits
                    process_tiff_strip(/*tiff_bitspersample=*/ 16, tiff_samplesperpixel, 0, height, (uint8_t *) imagebuf.get());
                } else {
                    throw std::runtime_error("Not implemented yet...");
                }
            } else {
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

                        // The output of the JPEG decoder is always 16 bits
                        process_tiff_strip(/*tiff_bitspersample=*/ 16, tiff_samplesperpixel, 0, height, (uint8_t *) spooler.data());
                    }
                    _TIFFfree(tiffbuf);
                } else {
                    // No compreession, read as a plain TIFF file

                    readTiffImageData(tif, width, height, tiff_bitspersample, tiff_samplesperpixel, process_tiff_strip);
                }
            }
        }
        if (metadata) {
            getExifMetaData(tif, metadata);
        }
    } else {
        throw std::runtime_error("Couldn't read dng file.");
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
            setMetadata(tif, metadata, TIFFTAG_DATETIME);

            setMetadata(tif, metadata, TIFFTAG_CFAREPEATPATTERNDIM);
            setMetadata(tif, metadata, TIFFTAG_CFAPATTERN);

            setMetadata(tif, metadata, TIFFTAG_COLORMATRIX1);
            setMetadata(tif, metadata, TIFFTAG_COLORMATRIX2);
            setMetadata(tif, metadata, TIFFTAG_ASSHOTNEUTRAL);

            setMetadata(tif, metadata, TIFFTAG_CALIBRATIONILLUMINANT1);
            setMetadata(tif, metadata, TIFFTAG_CALIBRATIONILLUMINANT2);

            setMetadata(tif, metadata, TIFFTAG_BLACKLEVELREPEATDIM);
            setMetadata(tif, metadata, TIFFTAG_BLACKLEVEL);
            setMetadata(tif, metadata, TIFFTAG_WHITELEVEL);

            setMetadata(tif, metadata, TIFFTAG_BAYERGREENSPLIT);
            setMetadata(tif, metadata, TIFFTAG_BASELINEEXPOSURE);
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
