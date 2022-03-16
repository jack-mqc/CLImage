//
//  gls_image_tiff.cpp
//  CLImageTest
//
//  Created by Fabio Riccardi on 3/15/22.
//

#include "gls_image_tiff.h"

#include <assert.h>

#include <tiffio.h>
#include <math.h>
#include <float.h>
#include <span>

namespace gls {

// exception friendly pointer with destructor
template <typename T>
struct auto_ptr : std::unique_ptr<T, std::function<void(T*)>> {
   public:
    auto_ptr(T* val, std::function<void(T*)> destroyer)
        : std::unique_ptr<T, std::function<void(T*)>>(val, destroyer) {}

    // Allow the auto_ptr to be implicitly converted to the underlying pointer type
    operator T*() const { return this->get(); }
    operator T*() { return this->get(); }
};

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
                        throw std::runtime_error("invalid strip");
                    }
                    process_tiff_strip(tiff_bitspersample, tiff_samplesperpixel, row, nrow, tiffbuf);
                }
            } else {
                throw std::runtime_error("error allocating memory buffer for TIFF strip");
            }
        } else {
            throw std::runtime_error("Couldn't allocate image storage");
        }
    }
}

template <typename T>
bool write_tiff_file(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                     /* int compression, */ std::function<T*(int row)> row_pointer) {
    auto_ptr<TIFF> tif(TIFFOpen(filename.c_str(), "w"),
                       [](TIFF *tif) { TIFFClose(tif); });
    if (tif) {
        uint16_t frame_width = width;
        uint16_t frame_height = height;
        uint32_t rowsperstrip = (uint32_t)-1;

        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, frame_width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, frame_height);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE /* compression */);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, pixel_channels > 2 ? PHOTOMETRIC_RGB : PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, pixel_bit_depth);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, pixel_channels);
        rowsperstrip = TIFFDefaultStripSize(tif, rowsperstrip);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rowsperstrip);
        TIFFSetField(tif, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

        // TIFFSetField(tif,TIFFTAG_IMAGEDESCRIPTION, imgd);

        // write frame data
        // data is broken up into strips where each strip contains rowsperstrip complete rows of data
        // each stript then has a size of rowsperstrip*frame_width pixels. the last strip is possibly
        // smaller, so it is NOT padded with dummy data.
        T* const buf = (T*)_TIFFmalloc(TIFFStripSize(tif));  // data buffer for a strip of the image
        for (unsigned int row = 0; (row < frame_height); row += rowsperstrip) {
            // compute rows in this strip:
            uint32_t nrow = rowsperstrip;
            if ((row + rowsperstrip) > frame_height) {
                nrow = frame_height - row;  // this is the last strip ... and it is a bit smaller! ...
                                            // it only contains the last rows of the image
            }
            tstrip_t strip = TIFFComputeStrip(tif, row, 0);
            tsize_t bi = 0;
            // go through the fraem row-wise
            for (int y = 0; y < nrow; ++y) {
                for (int x = 0; x < frame_width; ++x) {  // go through all pixels in the current row
                    for (int c = 0; c < pixel_channels; c++) {
                        buf[bi++] = row_pointer(row + y)[pixel_channels * x + c];
                    }
                }
            }
            if (TIFFWriteEncodedStrip(tif, strip, buf, bi * sizeof(T)) < 0) {
                return false;
            }
        }
        return true;
    }
    return false;
}

template
bool write_tiff_file<uint8_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                              /* int compression, */ std::function<uint8_t*(int row)> row_pointer);

template
bool write_tiff_file<uint16_t>(const std::string& filename, int width, int height, int pixel_channels, int pixel_bit_depth,
                               /* int compression, */ std::function<uint16_t*(int row)> row_pointer);

}  // namespace gls
