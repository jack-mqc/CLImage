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
    const TIFFDataType  field_type;         /* type of associated data */
    const char*         field_name;         /* ASCII name */
};

template<class T>
bool getMetaData(TIFF* tif, const TIFFField* tf, tiff_metadata* metadata, const std::string& key) {
    const auto field_tag = TIFFFieldTag(tf);
    const auto field_readcount = TIFFFieldReadCount(tf);

    if (field_readcount == TIFF_VARIABLE2 || field_readcount == TIFF_VARIABLE || field_readcount > 1) {
        size_t count = 0;
        T* data;

        if (field_readcount == TIFF_VARIABLE) {
            uint16_t gotcount = 0;
            TIFFGetField(tif, field_tag, &gotcount, &data);
            count = gotcount;
        } else if (field_readcount == TIFF_VARIABLE2) {
            uint32_t gotcount = 0;
            TIFFGetField(tif, field_tag, &gotcount, &data);
            count = gotcount;
        } else {
            TIFFGetField(tif, field_tag, &data);
            count = field_readcount;
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
    } else if (field_readcount == 1) {
        T data;
        TIFFGetField(tif, field_tag, &data);

        std::cout << "New metadata scalar " << key << ": " << data << std::endl;

        metadata->insert({ key, data });
        return true;
    }
    return false;
}

bool getMetaDataString(TIFF* tif, const TIFFField* tf, tiff_metadata* metadata, const std::string& key) {
    const auto field_readcount = TIFFFieldReadCount(tf);
    if (field_readcount > 1) {
        char* data;
        TIFFGetField(tif, TIFFFieldTag(tf), &data);

        std::cout << "New metadata string " << key << ": " << data << std::endl;

        metadata->insert({ key, data });
        return true;
    }
    return false;
}

void getMetadata(TIFF* tif, ttag_t field_tag, tiff_metadata* metadata) {
    const TIFFField* tf = TIFFFieldWithTag(tif, field_tag);

    const auto field_name = TIFFFieldName(tf);

    const std::string exifName = (field_name != nullptr) ? std::string(field_name) : std::to_string(field_tag);

    const auto field_type = TIFFFieldDataType(tf);
    switch (field_type) {
        case TIFF_BYTE: {
            getMetaData<uint8_t>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_UNDEFINED: {
            getMetaData<uint8_t>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_ASCII: {
            getMetaDataString(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_SHORT: {
            getMetaData<uint16_t>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_LONG: {
            getMetaData<uint32_t>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_SBYTE: {
            getMetaData<int8_t>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_SSHORT: {
            getMetaData<int16_t>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_SLONG: {
            getMetaData<int32_t>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_SRATIONAL:
        case TIFF_RATIONAL:
        case TIFF_FLOAT: {
            getMetaData<float>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_DOUBLE: {
            getMetaData<double>(tif, tf, metadata, exifName);
            break;
        }
        case TIFF_IFD:
        case TIFF_IFD8:
            std::cout << "Skipping offset field: " << field_name << std::endl;
            break;
        default:
            throw std::runtime_error("Unknown TIFF field type: " + std::to_string(field_type));
    }
}

void getAllTags(TIFF* tif, tiff_metadata* metadata) {
    if (tif) {
        int tag_count = TIFFGetTagListCount(tif);
        for (int i = 0; i < tag_count; i++) {
            ttag_t field_tag = TIFFGetTagListEntry(tif, i);

            getMetadata(tif, field_tag, metadata);
        }
    }
}

void getExifMetaData(TIFF* tif, tiff_metadata* metadata) {
    if (tif) {
        uint32_t exif_offset;
        if (TIFFGetField(tif, TIFFTAG_EXIFIFD, &exif_offset)) {
            TIFFReadEXIFDirectory(tif, exif_offset);

            getAllTags(tif, metadata);
        }
    }
}

template <typename T>
void setMetadata(TIFF* tif, const TIFFField* tf, const tiff_metadata_item& item) {
    const auto writeCount = TIFFFieldWriteCount(tf);
    if (writeCount == 1) {
        const auto value = std::get<T>(item);
        TIFFSetField(tif, TIFFFieldTag(tf), value);
    } else {
        const auto value = std::get<std::vector<T>>(item);
        if (writeCount < 0) {
            TIFFSetField(tif, TIFFFieldTag(tf), (uint16_t) value.size(), value.data());
        } else {
            if (writeCount != value.size()) {
                throw std::runtime_error("Vector size mismatch, should be: " + std::to_string(writeCount) + ", got: " + std::to_string(value.size()));
            }
            TIFFSetField(tif, TIFFFieldTag(tf), value.data());
        }
    }
}

void setMetadataString(TIFF* tif, const TIFFField* tf, const tiff_metadata_item& item) {
    const auto string = std::get<std::string>(item);
    const auto writeCount = TIFFFieldWriteCount(tf);
    if (writeCount < 0) {
        TIFFSetField(tif, TIFFFieldTag(tf), (uint16_t) string.size(), string.c_str());
    } else {
        const auto stringSize = string.size() + 1; // The string is null terminated
        if (writeCount != stringSize) {
            throw std::runtime_error("String size mismatch, should be: " + std::to_string(writeCount) + ", got: " + std::to_string(stringSize));
        }
        TIFFSetField(tif, TIFFFieldTag(tf), string.c_str());
    }
}

void setMetadata(TIFF* tif, tiff_metadata* metadata, const std::string& key) {
    const auto entry = metadata->find(key);
    if (entry != metadata->end()) {
        const TIFFField* tf = TIFFFieldWithName(tif, key.c_str());
        if (tf) {
            const auto field_type = TIFFFieldDataType(tf);

            switch (field_type) {
                case TIFF_BYTE: {
                    setMetadata<uint8_t>(tif, tf, entry->second);
                    break;
                }
                case TIFF_UNDEFINED: {
                    setMetadata<uint8_t>(tif, tf, entry->second);
                    break;
                }
                case TIFF_ASCII: {
                    setMetadataString(tif, tf, entry->second);
                    break;
                }
                case TIFF_SHORT: {
                    setMetadata<uint16_t>(tif, tf, entry->second);
                    break;
                }
                case TIFF_LONG: {
                    setMetadata<uint32_t>(tif, tf, entry->second);
                    break;
                }
                case TIFF_SBYTE: {
                    setMetadata<int8_t>(tif, tf, entry->second);
                    break;
                }
                case TIFF_SSHORT: {
                    setMetadata<int16_t>(tif, tf, entry->second);
                    break;
                }
                case TIFF_SLONG: {
                    setMetadata<int32_t>(tif, tf, entry->second);
                    break;
                }
                case TIFF_SRATIONAL:
                case TIFF_RATIONAL:
                case TIFF_FLOAT: {
                    setMetadata<float>(tif, tf, entry->second);
                    break;
                }
                case TIFF_DOUBLE: {
                    setMetadata<double>(tif, tf, entry->second);
                    break;
                }
                default:
                    throw std::runtime_error("Unknown TIFF field type: " + std::to_string(field_type));
            }
        }
    }
}

static const TIFFFieldInfo xtiffFieldInfo[] = {
    { TIFFTAG_DNG_IMAGEWIDTH, -1, -1, TIFF_LONG, FIELD_CUSTOM, 1, 0, "DNG ImageWidth" },
    { TIFFTAG_DNG_IMAGEHEIGHT, -1, -1, TIFF_LONG, FIELD_CUSTOM, 1, 0, "DNG ImageHeight" },
    { TIFFTAG_DNG_BITSPERSAMPLE, -1, -1, TIFF_LONG, FIELD_CUSTOM, 1, 0, "DNG BitsPerSample" },

    { TIFFTAG_FORWARDMATRIX1, -1, -1, TIFF_SRATIONAL, FIELD_CUSTOM, 1, 1, "ForwardMatrix1" },
    { TIFFTAG_FORWARDMATRIX2, -1, -1, TIFF_SRATIONAL, FIELD_CUSTOM, 1, 1, "ForwardMatrix2" },
    { TIFFTAG_TIMECODES, -1, -1, TIFF_BYTE, FIELD_CUSTOM, 1, 1, "TimeCodes" },
    { TIFFTAG_FRAMERATE, -1, -1, TIFF_SRATIONAL, FIELD_CUSTOM, 1, 1, "FrameRate" },
    { TIFFTAG_REELNAME, -1, -1, TIFF_ASCII, FIELD_CUSTOM, 1, 0, "ReelName" },

    { TIFFTAG_PROFILENAME, -1, -1, TIFF_ASCII, FIELD_CUSTOM, 1, 0, "ProfileName" },
    { TIFFTAG_PROFILELOOKTABLEDIMS, -1, -1, TIFF_LONG, FIELD_CUSTOM, 1, 0, "ProfileLookTableDims" },
    { TIFFTAG_PROFILELOOKTABLEDATA, -1, -1, TIFF_FLOAT, FIELD_CUSTOM, 1, 0, "ProfileLookTableData" },
    { TIFFTAG_PROFILELOOKTABLEENCODING, -1, -1, TIFF_LONG, FIELD_CUSTOM, 1, 0, "ProfileLookTableEncoding" },
    { TIFFTAG_DEFAULTUSERCROP, -1, -1, TIFF_RATIONAL, FIELD_CUSTOM, 1, 0, "DefaultUserCrop" },

    { TIFFTAG_RATING, -1, -1, TIFF_SHORT, FIELD_CUSTOM, 1, 0, "Rating" },
    { TIFFTAG_RATINGPERCENT, -1, -1, TIFF_SHORT, FIELD_CUSTOM, 1, 0, "RatingPercent" },
    { TIFFTAG_TIFFEPSTANDARDID, -1, -1, TIFF_SHORT, FIELD_CUSTOM, 1, 0, "TIFF-EP Standard ID" },
};

static TIFFExtendProc parent_extender = NULL;  // In case we want a chain of extensions

static void registerCustomTIFFTags(TIFF *tif) {
    // Install the extended Tag field info
    TIFFMergeFieldInfo( tif, xtiffFieldInfo, sizeof( xtiffFieldInfo ) / sizeof( xtiffFieldInfo[0] ) );

    if (parent_extender)
        parent_extender(tif);
}

void augment_libtiff_with_custom_tags() {
    static bool first_time = true;
    if (first_time) {
        parent_extender = TIFFSetTagExtender( registerCustomTIFFTags );
        first_time = false;
    }
}

}  // namespace gls
