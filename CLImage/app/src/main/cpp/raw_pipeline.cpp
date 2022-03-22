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

static const char* TAG = "RawPipeline Test";

int main(int argc, const char* argv[]) {
    printf("RawPipeline Test!\n");

    if (argc > 1) {
        auto input_path = std::filesystem::path(argv[1]);
        auto input_dir = input_path.parent_path();

        LOG_INFO(TAG) << "Processing: " << input_path.filename() << std::endl;

        // Read the input file into an image object
        auto inputImage = gls::image<gls::luma_pixel_16>::read_png_file(input_path.string());

        const auto rgb_image = demosaicImage(*inputImage);

        rgb_image->write_tiff_file(input_path.replace_extension("_rgb.tiff").c_str());

        LOG_INFO(TAG) << "done with inputImage size: " << inputImage->width << " x " << inputImage->height << std::endl;

        auto output_file = input_path.replace_extension(".dng").c_str();

        inputImage->write_dng_file(output_file);
    }
}
