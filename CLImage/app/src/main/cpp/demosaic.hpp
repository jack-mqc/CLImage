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

#ifndef demosaic_hpp
#define demosaic_hpp

#include "gls_image.hpp"
#include "gls_tiff_metadata.hpp"

enum BayerPattern {
    grbg = 0,
    gbrg = 1,
    rggb = 2,
    bggr = 3
};

void interpolateGreen(const gls::image<gls::luma_pixel_16>& rawImage,
                      gls::image<gls::rgb_pixel_16>* rgbImage, BayerPattern bayerPattern);

void interpolateRedBlue(gls::image<gls::rgb_pixel_16>* image, BayerPattern bayerPattern);

gls::image<gls::rgb_pixel_16>::unique_ptr demosaicImage(const gls::image<gls::luma_pixel_16>& rawImage,
                                                        const gls::tiff_metadata& metadata);

gls::image<gls::rgba_pixel>::unique_ptr demosaicImageGPU(const gls::image<gls::luma_pixel_16>& rawImage,
                                                         const gls::tiff_metadata& metadata);
#endif /* demosaic_hpp */
