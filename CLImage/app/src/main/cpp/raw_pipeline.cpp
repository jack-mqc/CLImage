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
#include <cmath>

#include "gls_logging.h"
#include "gls_image.hpp"

#include "demosaic.hpp"
#include "gls_tiff_metadata.hpp"

#include "gls_linalg.hpp"

static const char* TAG = "RawPipeline Test";

inline uint16_t clamp(int x) { return x < 0 ? 0 : x > 0xffff ? 0xffff : x; }

inline uint8_t clamp8(int x) { return x < 0 ? 0 : x > 0xff ? 0xff : x; }

// IMX492 M43-ish Sony Sensor
void IMX492Metadata(gls::tiff_metadata *metadata) {
    metadata->insert({ TIFFTAG_COLORMATRIX1, std::vector<float>{ 1.0781, -0.4173, -0.0976, -0.0633, 0.9661, 0.0972, 0.0073, 0.1349, 0.3481 } });
    metadata->insert({ TIFFTAG_ASSHOTNEUTRAL, std::vector<float>{ 1 / 2.001460, 1, 1 / 1.864002 } });
    metadata->insert({ TIFFTAG_CFAREPEATPATTERNDIM, std::vector<uint16_t>{ 2, 2 } });
    metadata->insert({ TIFFTAG_CFAPATTERN, std::vector<uint8_t>{ 1, 2, 0, 1 } });
    metadata->insert({ TIFFTAG_BLACKLEVEL, std::vector<float>{ 0 } });
    metadata->insert({ TIFFTAG_WHITELEVEL, std::vector<uint32_t>{ 0xfff } });
//    metadata->insert({ EXIFTAG_ISOSPEED, std::vector<uint16_t>{ 100 } });
//    metadata->insert({ EXIFTAG_SHUTTERSPEEDVALUE, 1.0/100.0});
}

// IMX571 APS-C Sony Sensor
void IMX571Metadata(gls::tiff_metadata *metadata) {
    metadata->insert({ TIFFTAG_COLORMATRIX1, std::vector<float>{ 2.5251, -1.3908, -0.3936, -0.5996, 1.7697, -0.1700, 0.2232, -0.2430, 1.2527 } });
    metadata->insert({ TIFFTAG_ASSHOTNEUTRAL, std::vector<float>{ 0.572128, 1.000000, 1.313796 } });
    metadata->insert({ TIFFTAG_CFAREPEATPATTERNDIM, std::vector<uint16_t>{ 2, 2 } });
    metadata->insert({ TIFFTAG_CFAPATTERN, std::vector<uint8_t>{ 2, 1, 1, 0 } });
    metadata->insert({ TIFFTAG_BLACKLEVEL, std::vector<float>{ 0 } });
    metadata->insert({ TIFFTAG_WHITELEVEL, std::vector<uint32_t>{ 0xffff } });
}

static float sigmoid(float x, float s) {
    return 0.5 * (tanh(s * x - 0.3 * s) + 1);
}

static float toneCurve(float x) {
    float s = 3.5;
    return (sigmoid(pow(x, 1/2.2), s) - sigmoid(0, s)) / (sigmoid(1, s) - sigmoid(0, s));
}

gls::image<gls::rgb_pixel>::unique_ptr applyToneCurve(const gls::image<gls::rgb_pixel_16>& image) {
    auto output_image = std::make_unique<gls::image<gls::rgb_pixel>>(image.width, image.height);

    int max_val = 0;
    int max_out = 0;
    for (int y = 0; y < image.height; y++) {
        for (int x = 0; x < image.width; x++) {
            const gls::rgb_pixel_16 &p = image[y][x];

            auto maxp = std::max(p.red, std::max(p.red, p.blue));
            if (maxp > max_val) {
                max_val = maxp;
            }

            auto op = (*output_image)[y][x] = {
                clamp8(0xff * toneCurve(0.8 * p[0] / (float) 0xffff)),
                clamp8(0xff * toneCurve(0.8 * p[1] / (float) 0xffff)),
                clamp8(0xff * toneCurve(0.8 * p[2] / (float) 0xffff))
            };

            auto maxop = std::max(op.red, std::max(op.red, op.blue));
            if (maxop > max_out) {
                max_out = maxop;
            }
        }
    }

    printf("max_val: %d, max_out: %d, tc(1): %f\n", max_val, max_out, toneCurve(1.0));

    return output_image;
}

int main(int argc, const char* argv[]) {
    printf("RawPipeline Test!\n");

    if (argc > 1) {
        auto input_path = std::filesystem::path(argv[1]);

        LOG_INFO(TAG) << "Processing: " << input_path.filename() << std::endl;

        gls::tiff_metadata metadata;
        const auto inputImage = gls::image<gls::luma_pixel_16>::read_dng_file(input_path.string(), &metadata);

        // inputImage->write_png_file((input_path.parent_path() / input_path.stem()).string() + "_raw.png");

//        gls::tiff_metadata metadata;
//        IMX492Metadata(&metadata);
//        const auto inputImage = gls::image<gls::luma_pixel_16>::read_png_file(input_path.string());

//        for (auto& p : inputImage->pixels()) {
//            p = gls::luma_pixel_16 {
//                (uint16_t) (p[0] >> 4)
//            };
//        }

        LOG_INFO(TAG) << "read inputImage of size: " << inputImage->width << " x " << inputImage->height << std::endl;

        const bool useGPU = true;
        if (useGPU) {
            const auto rgb_image = demosaicImageGPU(*inputImage, &metadata, /*auto_white_balance=*/ true);
            rgb_image->write_png_file((input_path.parent_path() / input_path.stem()).string() + "_rgb.png", /*skip_alpha=*/ true);
        } else {
            const auto rgb_image = demosaicImage(*inputImage, &metadata, /*auto_white_balance=*/ true);
            auto output_image = applyToneCurve(*rgb_image);
            LOG_INFO(TAG) << "...done with demosaicing (CPU)." << std::endl;
            output_image->write_png_file((input_path.parent_path() / input_path.stem()).string() + "_rgb.png");
        }

        auto output_file = (input_path.parent_path() / input_path.stem()).string() + "_my.dng";
        inputImage->write_dng_file(output_file, /*compression=*/ gls::JPEG, &metadata);
    }
}
