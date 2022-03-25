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

#include <filesystem>
#include <string>

#include "gls_logging.h"
#include "gls_image.hpp"

#include "demosaic.hpp"
#include "gls_tiff_metadata.hpp"

static const char* TAG = "RawPipeline Test";

inline uint16_t clamp(int x) { return x < 0 ? 0 : x > 0xffff ? 0xffff : x; }

//static const short bayerPatternDimensions[] = { 2, 2 };
//static const unsigned char bayerPattern[] = "\00\01\01\02"; // "\01\02\00\01";
//
//std::array<float, 9> color_matrix = {
//    1.9435, -0.8992, -0.1936, 0.1144, 0.8380, 0.0475, 0.0136, 0.1203, 0.3553
//};
//std::array<float, 3> camera_multipliers = {
//    1.354894, 1.000000, 1.920234
//};
//const std::array<float, 3> as_shot_neutral = {
//    1.f / (camera_multipliers[0] / camera_multipliers[1]),
//    1.f,
//    1.f / (camera_multipliers[2] / camera_multipliers[1])
//};

int main(int argc, const char* argv[]) {
    printf("RawPipeline Test!\n");

    if (argc > 1) {
        auto input_path = std::filesystem::path(argv[1]);

        LOG_INFO(TAG) << "Processing: " << input_path.filename() << std::endl;

        // gls::tiff_metadata metadata;
        // const auto inputImage = gls::image<gls::luma_pixel_16>::read_dng_file(input_path.string(), &metadata);

        gls::tiff_metadata metadata;
        metadata.insert({ "ColorMatrix1", std::vector<float>{ 1.9435, -0.8992, -0.1936, 0.1144, 0.8380, 0.0475, 0.0136, 0.1203, 0.3553 } });
        metadata.insert({ "AsShotNeutral", std::vector<float>{ 0.7380, 1, 0.5207 } });
        metadata.insert({ "CFARepeatPatternDim", std::vector<uint16_t>{ 2, 2 } });
        metadata.insert({ "CFAPattern", std::vector<uint8_t>{ 1, 2, 0, 1 } });
        metadata.insert({ "BlackLevel", std::vector<float>{ 0 } });
        metadata.insert({ "WhiteLevel", std::vector<uint32_t>{ 0x0fff } });

        const auto inputImage = gls::image<gls::luma_pixel_16>::read_png_file(input_path.string());

        LOG_INFO(TAG) << "read inputImage of size: " << inputImage->width << " x " << inputImage->height << std::endl;

        const auto rgb_image = demosaicImage(*inputImage, metadata);

        rgb_image->write_tiff_file(input_path.replace_extension("_rgb.tiff").c_str());

        LOG_INFO(TAG) << "done with inputImage size: " << inputImage->width << " x " << inputImage->height << std::endl;

        auto output_file = input_path.replace_extension("_my.dng").c_str();

        inputImage->write_dng_file(output_file, /*compression=*/ gls::JPEG, &metadata);
    }
}
