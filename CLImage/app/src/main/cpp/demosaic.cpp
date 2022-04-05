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

#include <limits.h>
#include <math.h>

#include "demosaic.hpp"
#include "gls_color_science.hpp"
#include "gls_linalg.hpp"

enum { red = 0, green = 1, blue = 2, green2 = 3 };

static const std::array<std::array<gls::point, 4>, 4> bayerOffsets = {{
    { gls::point{1, 0}, gls::point{0, 0}, gls::point{0, 1}, gls::point{1, 1} }, // grbg
    { gls::point{0, 1}, gls::point{0, 0}, gls::point{1, 0}, gls::point{1, 1} }, // gbrg
    { gls::point{0, 0}, gls::point{1, 0}, gls::point{1, 1}, gls::point{0, 1} }, // rggb
    { gls::point{1, 1}, gls::point{1, 0}, gls::point{0, 0}, gls::point{0, 1} }  // bggr
}};

inline uint16_t clamp(int x) { return x < 0 ? 0 : x > 0xffff ? 0xffff : x; }

void interpolateGreen(const gls::image<gls::luma_pixel_16>& rawImage,
                      gls::image<gls::rgb_pixel_16>* rgbImage, BayerPattern bayerPattern) {
    const int width = rawImage.width;
    const int height = rawImage.height;

    auto offsets = bayerOffsets[bayerPattern];
    const gls::point r = offsets[red];
    const gls::point g = offsets[green];

    // copy RAW data to RGB layer and remove hot pixels
    for (int y = 0; y < height; y++) {
        int color = (y & 1) == (r.y & 1) ? red : blue;
        int x0 = (y & 1) == (g.y & 1) ? g.x + 1 : g.x;
        for (int x = 0; x < width; x++) {
            bool colorPixel = (x & 1) == (x0 & 1);
            int channel = colorPixel ? color : green;

            int value = rawImage[y][x];
            if (x >= 2 && x < width - 2 && y >= 2 && y < height - 2) {
                int v[12];
                int n;
                if (!colorPixel) {
                    n = 8;
                    v[0] = rawImage[y - 1][x - 1];
                    v[1] = rawImage[y - 1][x + 1];
                    v[2] = rawImage[y + 1][x - 1];
                    v[3] = rawImage[y + 1][x + 1];

                    v[4] = 2 * rawImage[y - 1][x];
                    v[5] = 2 * rawImage[y + 1][x];
                    v[6] = 2 * rawImage[y][x + 1];
                    v[7] = 2 * rawImage[y][x + 1];
                } else {
                    n = 12;
                    v[0] = rawImage[y - 2][x];
                    v[1] = rawImage[y + 2][x];
                    v[2] = rawImage[y][x - 2];
                    v[3] = rawImage[y][x + 2];

                    v[4] = 2 * rawImage[y - 1][x - 1];
                    v[5] = 2 * rawImage[y - 1][x + 1];
                    v[6] = 2 * rawImage[y + 1][x - 1];
                    v[7] = 2 * rawImage[y + 1][x + 1];

                    v[8] = 2 * rawImage[y - 1][x];
                    v[9] = 2 * rawImage[y + 1][x];
                    v[10] = 2 * rawImage[y][x - 1];
                    v[11] = 2 * rawImage[y][x + 1];
                };
                bool replace = true;
                for (int i = 0; i < n; i++)
                    if (value < 2 * v[i]) {
                        replace = false;
                        break;
                    }
                if (replace) value = (v[0] + v[1] + v[2] + v[3]) / 4;
            }
            (*rgbImage)[y][x][channel] = value;
        }
    }

    // green channel interpolation

    for (int y = 2; y < height - 2; y++) {
        int color = (y & 1) == (r.y & 1) ? red : blue;
        int x0 = (y & 1) == (g.y & 1) ? g.x + 1 : g.x;

        int g_left  = (*rgbImage)[y][x0 - 1][green];
        int c_xy    = (*rgbImage)[y][x0][color];
        int c_left  = (*rgbImage)[y][x0 - 2][color];

        for (int x = x0 + 2; x < width - 2; x += 2) {
            int g_right = (*rgbImage)[y][x + 1][green];
            int g_up    = (*rgbImage)[y - 1][x][green];
            int g_down  = (*rgbImage)[y + 1][x][green];
            int g_dh    = abs(g_left - g_right);
            int g_dv    = abs(g_up - g_down);

            int c_right = (*rgbImage)[y][x + 2][color];
            int c_up    = (*rgbImage)[y - 2][x][color];
            int c_down  = (*rgbImage)[y + 2][x][color];
            int c_dh    = abs(c_left + c_right - 2 * c_xy);
            int c_dv    = abs(c_up + c_down - 2 * c_xy);

            // Minimum derivative value for edge directed interpolation (avoid aliasing)
            int dThreshold = 1200;

            // we're doing edge directed bilinear interpolation on the green channel,
            // which is a low pass operation (averaging), so we add some signal from the
            // high frequencies of the observed color channel

            int sample;
            if (g_dv + c_dv > dThreshold && g_dv + c_dv > g_dh + c_dh) {
                sample = (g_left + g_right) / 2;
                if (sample < 4 * c_xy && c_xy < 4 * sample) {
                    sample += (c_xy - (c_left + c_right) / 2) / 4;
                }
            } else if (g_dh + c_dh > dThreshold && g_dh + c_dh > g_dv + c_dv) {
                sample = (g_up + g_down) / 2;
                if (sample < 4 * c_xy && c_xy < 4 * sample) {
                    sample += (c_xy - (c_up + c_down) / 2) / 4;
                }
            } else {
                sample = (g_up + g_left + g_down + g_right) / 4;
                if (sample < 4 * c_xy && c_xy < 4 * sample) {
                    sample += (c_xy - (c_left + c_right + c_up + c_down) / 4) / 8;
                }
            }

            (*rgbImage)[y][x][green] = clamp(sample);
            g_left = g_right;
            c_left = c_xy;
            c_xy   = c_right;
        }
    }

    // get the constant component out of the reconstructed green pixels and add to it
    // the "high frequency" part of the corresponding observed color channel

    for (int y = 2; y < height - 2; y++) {
        int channel = (y & 1) == (r.y & 1) ? red : blue;
        int x0 = 2 + ((y & 1) == (g.y & 1) ? g.x + 1 : g.x);

        int g_xy            = (*rgbImage)[y][x0][green];
        int g_left          = (*rgbImage)[y][x0 - 2][green];
        int g_top_left      = (*rgbImage)[y - 2][x0 - 2][green];
        int g_bottom_left   = (*rgbImage)[y + 2][x0 - 2][green];

        int c_xy            = (*rgbImage)[y][x0][channel];
        int c_left          = (*rgbImage)[y][x0 - 2][channel];
        int c_top_left      = (*rgbImage)[y - 2][x0 - 2][channel];
        int c_bottom_left   = (*rgbImage)[y + 2][x0 - 2][channel];

        for (int x = 2; x < width - 2; x += 2) {
            int g_right         = (*rgbImage)[y][x + 2][green];
            int g_top_right     = (*rgbImage)[y - 2][x + 2][green];
            int g_bottom_right  = (*rgbImage)[y + 2][x + 2][green];
            int g_up            = (*rgbImage)[y - 2][x][green];
            int g_down          = (*rgbImage)[y + 2][x][green];

            int c_right         = (*rgbImage)[y][x + 2][channel];
            int c_top_right     = (*rgbImage)[y - 2][x + 2][channel];
            int c_bottom_right  = (*rgbImage)[y + 2][x + 2][channel];
            int c_up            = (*rgbImage)[y - 2][x][channel];
            int c_down          = (*rgbImage)[y + 2][x][channel];

            // Only work on the pixels that have a strong enough correlation between channels

            if (g_xy < 4 * c_xy && c_xy < 4 * g_xy) {
                int dh = g_xy - (g_left + g_right) / 2;
                int dv = g_xy - (g_up + g_down) / 2;
                int ne = g_xy - (g_top_left + g_bottom_right) / 2;
                int nw = g_xy - (g_top_right + g_bottom_left) / 2;

                int cdh = c_xy - (c_left + c_right) / 2;
                int cdv = c_xy - (c_up + c_down) / 2;
                int cne = c_xy - (c_top_left + c_bottom_right) / 2;
                int cnw = c_xy - (c_top_right + c_bottom_left) / 2;

                enum GradientDirection {
                    horizontal = 0,
                    vertical = 1,
                    northEast = 2,
                    northWest = 3,
                    none = 4
                };

                int gradients[4] = {
                    abs(dh) + abs(cdh), // horizontal
                    abs(dv) + abs(cdv), // vertical
                    abs(ne) + abs(cne), // northEast
                    abs(nw) + abs(cnw)  // northWest
                };

                GradientDirection minimumDirection = none;
                int minimumGradient = INT_MAX;
                for (GradientDirection g : { horizontal, vertical, northEast, northWest }) {
                    if (gradients[g] < minimumGradient) {
                        minimumDirection = g;
                        minimumGradient = gradients[g];
                    }
                }

                // Only work on parts of the image that have enough "detail"

                if (minimumDirection != none && minimumGradient > g_xy / 4) {
                    int sample;
                    switch (minimumDirection) {
                        case horizontal:
                            sample = (g_xy + (g_left + g_right) / 2 + cdh) / 2;
                            break;
                        case vertical:
                            sample = (g_xy + (g_up + g_down) / 2 + cdv) / 2;
                            break;
                        case northEast:
                            sample = (g_xy + (g_top_left + g_bottom_right) / 2 + cne) / 2;
                            break;
                        case northWest:
                            sample = (g_xy + (g_top_right + g_bottom_left) / 2 + cnw) / 2;
                            break;
                        case none:
                            // never happens, just make the compiler happy
                            sample = (*rgbImage)[y][x][green];
                            break;
                    }

                    (*rgbImage)[y][x][green] = clamp(sample);
                }
            }

            g_left          = g_xy;
            g_xy            = g_right;
            g_top_left      = g_up;
            g_up            = g_top_right;
            g_bottom_left   = g_down;
            g_down          = g_bottom_right;
            c_left          = c_xy;
            c_xy            = c_right;
            c_top_left      = c_up;
            c_up            = c_top_right;
            c_bottom_left   = c_down;
            c_down          = c_bottom_right;
        }
    }
}

void interpolateRedBlue(gls::image<gls::rgb_pixel_16>* image, BayerPattern bayerPattern) {
    const int width = image->width;
    const int height = image->height;

    auto offsets = bayerOffsets[bayerPattern];

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            for (int color : {red, blue}) {
                const gls::point c = offsets[color];

                if (((x + c.x) & 1) != (c.x & 1) || ((y + c.y) & 1) != (c.y & 1)) {
                    int sample;
                    int cg = (*image)[y + c.y][x + c.x][green];

                    if (((x + c.x) & 1) != (c.x & 1) && ((y + c.y) & 1) != (c.y & 1)) {
                        // Pixel at color location
                        int g_ne = (*image)[y + c.y + 1][x + c.x - 1][green];
                        int g_nw = (*image)[y + c.y + 1][x + c.x + 1][green];
                        int g_sw = (*image)[y + c.y - 1][x + c.x + 1][green];
                        int g_se = (*image)[y + c.y - 1][x + c.x - 1][green];

                        int c_ne = g_ne - (*image)[y + c.y + 1][x + c.x - 1][color];
                        int c_nw = g_nw - (*image)[y + c.y + 1][x + c.x + 1][color];
                        int c_sw = g_sw - (*image)[y + c.y - 1][x + c.x + 1][color];
                        int c_se = g_se - (*image)[y + c.y - 1][x + c.x - 1][color];

                        int d_ne_sw = abs(c_ne - c_sw);
                        int d_nw_se = abs(c_nw - c_se);

                        // Minimum gradient for edge directed interpolation
                        int dThreshold = 800;
                        if (d_ne_sw > dThreshold && d_ne_sw > d_nw_se) {
                            sample = cg - (c_nw + c_se) / 2;
                        } else if (d_nw_se > dThreshold && d_nw_se > d_ne_sw) {
                            sample = cg - (c_ne + c_sw) / 2;
                        } else {
                            sample = cg - (c_ne + c_sw + c_nw + c_se) / 4;
                        }
                    } else if (((x + c.x) & 1) == (c.x & 1) && ((y + c.y) & 1) != (c.y & 1)) {
                        // Pixel at green location - vertical
                        int g_up    = (*image)[y + c.y - 1][x + c.x][green];
                        int g_down  = (*image)[y + c.y + 1][x + c.x][green];

                        int c_up    = g_up - (*image)[y + c.y - 1][x + c.x][color];
                        int c_down  = g_down - (*image)[y + c.y + 1][x + c.x][color];

                        sample = cg - (c_up + c_down) / 2;
                    } else {
                        // Pixel at green location - horizontal
                        int g_left  = (*image)[y + c.y][x + c.x - 1][green];
                        int g_right = (*image)[y + c.y][x + c.x + 1][green];

                        int c_left  = g_left - (*image)[y + c.y][x + c.x - 1][color];
                        int c_right = g_right - (*image)[y + c.y][x + c.x + 1][color];

                        sample = cg - (c_left + c_right) / 2;
                    }

                    (*image)[y + c.y][x + c.x][color] = clamp(sample);
                }
            }
        }
    }
}

// sRGB -> XYZ D65 Transform: xyz_rgb * rgb_color -> xyz_color
const gls::Matrix<3, 3> xyz_rgb = {
    { 0.4124564, 0.3575761, 0.1804375 },
    { 0.2126729, 0.7151522, 0.0721750 },
    { 0.0193339, 0.1191920, 0.9503041 }
};

// XYZ D65 -> sRGB Transform: rgb_xyz * xyx_color -> rgb_color
const gls::Matrix<3, 3> rgb_xyz = {
    {  3.2404542, -1.5371385, -0.4985314 },
    { -0.9692660,  1.8760108,  0.0415560 },
    {  0.0556434, -0.2040259,  1.0572252 }
};

gls::Matrix<3, 3> cam_xyz_coeff(gls::Vector<3>& pre_mul, const gls::Matrix<3, 3>& cam_xyz) {
    // Compute sRGB -> XYZ -> Camera
    auto rgb_cam = cam_xyz * xyz_rgb;

    // Normalize rgb_cam so that rgb_cam * (1,1,1) == (1,1,1).
    // This maximizes the uint16 dynamic range and makes sure
    // that highlight clipping is white in both camera and target
    // color spaces, so that clipping doesn't turn pink

    auto cam_white = rgb_cam * gls::Vector<3>({ 1, 1, 1 });

    gls::Matrix<3, 3> mPreMul = {
        { 1 / cam_white[0], 0, 0 },
        { 0, 1 / cam_white[1], 0 },
        { 0, 0, 1 / cam_white[2] }
    };

    for (int i = 0; i < 3; i++) {
        if (cam_white[i] > 0.00001) {
            pre_mul[i] = 1 / cam_white[i];
        } else {
            throw std::range_error("");
        }
    }

    // Return Camera -> sRGB
    return inverse(mPreMul * rgb_cam);
}

gls::image<gls::rgb_pixel_16>::unique_ptr demosaicImage(const gls::image<gls::luma_pixel_16>& rawImage,
                                                        const gls::tiff_metadata& metadata) {
    const auto color_matrix1 = getVector<float>(metadata, TIFFTAG_COLORMATRIX1);
    const auto color_matrix2 = getVector<float>(metadata, TIFFTAG_COLORMATRIX2);

    // If present ColorMatrix2 is usually D65 and ColorMatrix1 is Standard Light A
    const auto& color_matrix = color_matrix2.empty() ? color_matrix1 : color_matrix2;

    const auto as_shot_neutral = getVector<float>(metadata, TIFFTAG_ASSHOTNEUTRAL);
    const auto black_level_vec = getVector<float>(metadata, TIFFTAG_BLACKLEVEL);
    const auto white_level_vec = getVector<uint32_t>(metadata, TIFFTAG_WHITELEVEL);
    const auto cfa_pattern = getVector<uint8_t>(metadata, TIFFTAG_CFAPATTERN);

    const float black_level = black_level_vec.empty() ? 0 : black_level_vec[0];
    const uint32_t white_level = white_level_vec.empty() ? 0xffff : white_level_vec[0];

    const auto bayerPattern = std::memcmp(cfa_pattern.data(), "\00\01\01\02", 4) == 0 ? BayerPattern::rggb
                            : std::memcmp(cfa_pattern.data(), "\02\01\01\00", 4) == 0 ? BayerPattern::bggr
                            : std::memcmp(cfa_pattern.data(), "\01\00\02\01", 4) == 0 ? BayerPattern::grbg
                            : BayerPattern::gbrg;

    std::cout << "as_shot_neutral: " << gls::Vector<3>(as_shot_neutral) << std::endl;

    gls::Vector<3> cam_mul = 1.0 / gls::Vector<3>(as_shot_neutral);

    // TODO: this should be CameraCalibration * ColorMatrix * AsShotWhite
    gls::Matrix<3, 3> cam_xyz = color_matrix;

    std::cout << "cam_xyz:\n" << cam_xyz << std::endl;

    gls::Vector<3> pre_mul;
    const auto rgb_cam = cam_xyz_coeff(pre_mul, cam_xyz);

    std::cout << "*** pre_mul: " << pre_mul << std::endl;
    std::cout << "*** cam_mul: " << cam_mul << std::endl;

    // Save the whitening transformation
    const auto inv_cam_white = pre_mul;

    gls::Matrix<3, 3> mCamMul = {
        { pre_mul[0] / cam_mul[0], 0, 0 },
        { 0, pre_mul[1] / cam_mul[1], 0 },
        { 0, 0, pre_mul[2] / cam_mul[2] }
    };

    // If cam_mul is available use that instead of pre_mul
    for (int i = 0; i < 3; i++) {
        pre_mul[i] = cam_mul[i];
    }

    {
        const gls::Vector<3> d65_white = { 0.95047, 1.0, 1.08883 };
        const gls::Vector<3> d50_white = { 0.9642, 1.0000, 0.8249 };

        std::cout << "XYZ D65 White -> sRGB: " << rgb_xyz * d65_white << std::endl;
        std::cout << "sRGB White -> XYZ D65: " << xyz_rgb * gls::Vector({1, 1, 1}) << std::endl;
        std::cout << "inverse(cam_xyz) * (1 / pre_mul): " << inverse(cam_xyz) * (1 / pre_mul) << std::endl;

        auto cam_white = cam_xyz * xyz_rgb * gls::Vector<3>({ 1, 1, 1 });
        std::cout << "inverse(cam_xyz) * cam_white: " << inverse(cam_xyz) * cam_white << ", CCT: " << XYZtoCorColorTemp(inverse(cam_xyz) * cam_white) << std::endl;
        std::cout << "xyz_rgb * gls::Vector<3>({ 1, 1, 1 }): " << xyz_rgb * gls::Vector<3>({ 1, 1, 1 }) << ", CCT: " << XYZtoCorColorTemp(xyz_rgb * gls::Vector<3>({ 1, 1, 1 })) << std::endl;

        std::cout << "***>> pizza: " << XYZtoCorColorTemp(xyz_rgb * rgb_cam * ((1 / pre_mul) * gls::Vector<3>({ 1, 1, 1 }))) << std::endl;

        const auto wb_out = xyz_rgb * rgb_cam * mCamMul * gls::Vector({1, 1, 1});
        std::cout << "wb_out: " << wb_out << ", CCT: " << XYZtoCorColorTemp(wb_out) << std::endl;

        const auto no_wb_out = xyz_rgb * rgb_cam * (1 / inv_cam_white);
        std::cout << "no_wb_out: " << wb_out << ", CCT: " << XYZtoCorColorTemp(no_wb_out) << std::endl;
    }

    // Scale Input Image
    auto minmax = std::minmax_element(std::begin(pre_mul), std::end(pre_mul));
    gls::Vector<4> scale_mul;
    for (int c = 0; c < 4; c++) {
        int pre_mul_idx = c == 3 ? 1 : c;
        printf("pre_mul[c]: %f, *minmax.second: %f, white_level: %d\n", pre_mul[pre_mul_idx], *minmax.second, white_level);
        scale_mul[c] = (pre_mul[pre_mul_idx] / *minmax.first) * 65535.0 / (white_level - black_level);
    }
    printf("scale_mul: %f, %f, %f, %f\n", scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);

    const auto offsets = bayerOffsets[bayerPattern];
    gls::image<gls::luma_pixel_16> scaledRawImage = gls::image<gls::luma_pixel_16>(rawImage.width, rawImage.height);
    for (int y = 0; y < rawImage.height / 2; y++) {
        for (int x = 0; x < rawImage.width / 2; x++) {
            for (int c = 0; c < 4; c++) {
                const auto& o = offsets[c];
                scaledRawImage[2 * y + o.y][2 * x + o.x] = clamp(scale_mul[c] * (rawImage[2 * y + o.y][2 * x + o.x] - black_level));
            }
        }
    }

    auto rgbImage = std::make_unique<gls::image<gls::rgb_pixel_16>>(rawImage.width, rawImage.height);

    printf("interpolating green channel...\n");
    interpolateGreen(scaledRawImage, rgbImage.get(), bayerPattern);

    printf("interpolating red and blue channels...\n");
    interpolateRedBlue(rgbImage.get(), bayerPattern);

    printf("...done with demosaicing.\n");

    // Transform to RGB space
    for (int y = 0; y < rgbImage->height; y++) {
        for (int x = 0; x < rgbImage->width; x++) {
            auto& p = (*rgbImage)[y][x];
            const auto op = rgb_cam * /* mCamMul * */ gls::Vector<3>({ (float) p[0], (float) p[1], (float) p[2] });
            p = { clamp(op[0]), clamp(op[1]), clamp(op[2]) };
        }
    }

    return rgbImage;
}
