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

#include "gls_tiff_metadata.hpp"

#include <iostream>

namespace gls {

struct TiffFieldInfo {
    const ttag_t        field_tag;          /* field's tag */
    const int           field_readcount;    /* read count/TIFF_VARIABLE/TIFF_SPP */
    const int           field_writecount;   /* write count/TIFF_VARIABLE */
    const TIFFDataType  field_type;         /* type of associated data */
    const char*         field_name;         /* ASCII name */
};

template<class T>
bool getMetaData(TIFF* tif, const TiffFieldInfo& fi, tiff_metadata* metadata, const std::string& key) {
    if (fi.field_readcount == TIFF_VARIABLE2 || fi.field_readcount == TIFF_VARIABLE || fi.field_readcount > 1) {

        size_t count = 0;
        T* data;

        if (fi.field_readcount == TIFF_VARIABLE) {
            uint16_t gotcount = 0;
            TIFFGetField(tif, fi.field_tag, &gotcount, &data);
            count = gotcount;
        } else if (fi.field_readcount == TIFF_VARIABLE2) {
            uint32_t gotcount = 0;
            TIFFGetField(tif, fi.field_tag, &gotcount, &data);
            count = gotcount;
        } else {
            TIFFGetField(tif, fi.field_tag, &data);
            count = fi.field_readcount;
        }

        std::vector<T> values;
        values.resize(count);
        for (unsigned i = 0; i < count; i++) {
            values[i] = data[i];
        }

        std::cout << "New metadata vector (" << values.size() << ") " << key << ": ";
        for (int i = 0; i < values.size() && i < 10; i++) {
            const auto& v = values[i];
            if (sizeof(v) == 1) {
                std::cout << (int) v;
            } else {
                std::cout << v;
            }
            if (i < 9 && i < values.size() - 1) {
                std::cout << ", ";
            } else if (i == 9 && i < values.size() - 1) {
                std::cout << "...";
            }
        }
        std::cout << std::endl;

        metadata->insert({ key, values });
        return true;
    } else if (fi.field_readcount == 1) {
        T data;
        TIFFGetField(tif, fi.field_tag, &data);

        std::cout << "New metadata scalar " << key << ": " << data << std::endl;

        metadata->insert({ key, data });
        return true;
    }
    return false;
}

bool getMetaDataString(TIFF* tif, const TiffFieldInfo& fi, tiff_metadata* metadata, const std::string& key) {
    if (fi.field_readcount > 1) {
        char* data;
        TIFFGetField(tif, fi.field_tag, &data);

        std::cout << "New metadata string " << key << ": " << data << std::endl;

        metadata->insert({ key, data });
        return true;
    }
    return false;
}

void getAllTIFFTags(TIFF* tif, tiff_metadata* metadata) {
    if (tif) {
        int cnt = TIFFGetTagListCount(tif);
        for (int i = 0; i < cnt; i++) {
            ttag_t tag = TIFFGetTagListEntry(tif, i);
            const TIFFField* tfi = TIFFFieldWithTag(tif, tag);

            TiffFieldInfo fi = {
                .field_tag = TIFFFieldTag(tfi),
                .field_readcount = TIFFFieldReadCount(tfi),
                .field_writecount = TIFFFieldWriteCount(tfi),
                .field_type = TIFFFieldDataType(tfi),
                .field_name = TIFFFieldName(tfi),
            };

            const std::string exifName = (fi.field_name != nullptr) ? std::string(fi.field_name) : std::to_string(fi.field_tag);

            bool usedMetaData = false;

            switch (fi.field_type) {
                case TIFF_BYTE: {
                    usedMetaData = getMetaData<uint8_t>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_UNDEFINED: {
                    usedMetaData = getMetaData<uint8_t>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_ASCII: {
                    usedMetaData = getMetaDataString(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_SHORT: {
                    usedMetaData = getMetaData<uint16_t>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_LONG: {
                    usedMetaData = getMetaData<uint32_t>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_SBYTE: {
                    usedMetaData = getMetaData<int8_t>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_SSHORT: {
                    usedMetaData = getMetaData<int16_t>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_SLONG: {
                    usedMetaData = getMetaData<int32_t>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_SRATIONAL:
                case TIFF_RATIONAL:
                case TIFF_FLOAT: {
                    usedMetaData = getMetaData<float>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_DOUBLE: {
                    usedMetaData = getMetaData<double>(tif, fi, metadata, exifName);
                    break;
                }
                case TIFF_IFD:
                case TIFF_IFD8:
                    std::cout << "Skipping offser field: " << fi.field_name << std::endl;
                    break;
                default:
                    throw std::runtime_error("Unknown TIFF field type: " + std::to_string(fi.field_type));
            }
        }
    }
}

void fetchExifMetaData(TIFF* tif, tiff_metadata* metadata) {
    if (tif) {
        uint32_t exif_offset;
        if (TIFFGetField(tif, TIFFTAG_EXIFIFD, &exif_offset)) {
            TIFFReadEXIFDirectory(tif, exif_offset);

            getAllTIFFTags(tif, metadata);
        }
    }
}

}  // namespace gls
