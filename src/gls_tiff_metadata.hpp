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

#ifndef TiffMetadata_hpp
#define TiffMetadata_hpp

#include <variant>
#include <vector>
#include <map>
#include <string>

#include <tiffio.h>

namespace gls {

typedef std::variant<uint8_t, uint16_t, uint32_t, int8_t, int16_t, int32_t, float, double,
                     std::vector<uint8_t>, std::vector<uint16_t>, std::vector<uint32_t>,
                     std::vector<int8_t>, std::vector<int16_t>, std::vector<int32_t>,
                     std::vector<float>, std::vector<double>, std::string> tiff_metadata_item;

class tiff_metadata: public std::map<const std::string, const tiff_metadata_item> { };

void fetchExifMetaData(TIFF* tif, tiff_metadata* metadata);

void getAllTIFFTags(TIFF* tif, tiff_metadata* metadata);

template <typename T>
void writeVectorMetadata(TIFF* tif, tiff_metadata* metadata, const std::string& key) {
    const auto entry = metadata->find(key);
    if (entry != metadata->end()) {
        const auto value = std::get<std::vector<T>>(entry->second);
        const TIFFField* tf = TIFFFieldWithName(tif, key.c_str());
        if (tf) {
            TIFFSetField(tif, TIFFFieldTag(tf), value.size(), value.data());
        }
    }
}

template <typename T>
void writeScalarMetadata(TIFF* tif, tiff_metadata* metadata, const std::string& key) {
    const auto entry = metadata->find(key);
    if (entry != metadata->end()) {
        const auto value = std::get<T>(entry->second);
        const TIFFField* tf = TIFFFieldWithName(tif, key.c_str());
        if (tf) {
            TIFFSetField(tif, TIFFFieldTag(tf), value);
        }
    }
}

} // namespace gls

#endif /* TiffMetadata_hpp */
