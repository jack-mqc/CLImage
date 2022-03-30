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

        int hl = (*rgbImage)[y][x0 - 1][green];
        int cxy = (*rgbImage)[y][x0][color];
        int chl = (*rgbImage)[y][x0 - 2][color];

        for (int x = 2; x < width - 2; x++) {
            if ((x & 1) == (x0 & 1)) {
                int hr = (*rgbImage)[y][x + 1][green];
                int vu = (*rgbImage)[y - 1][x][green];
                int vd = (*rgbImage)[y + 1][x][green];
                int dh = abs(hl - hr);
                int dv = abs(vu - vd);

                int chr = (*rgbImage)[y][x + 2][color];
                int cvu = (*rgbImage)[y - 2][x][color];
                int cvd = (*rgbImage)[y + 2][x][color];
                int cdh = abs(chl + chr - 2 * cxy);
                int cdv = abs(cvu + cvd - 2 * cxy);

                // we're doing edge directed bilinear interpolation on the green channel,
                // which is a low pass operation (averaging), so we add some signal from the
                // high frequencies of the observed color channel

                int sample;
                if (dv + cdv - (dh + cdh) > 0) {
                    sample = (hl + hr) / 2;
                    if (sample < 4 * cxy && cxy < 4 * sample) sample += (cxy - (chl + chr) / 2) / 4;
                } else if (dh + cdh - (dv + cdv) > 0) {
                    sample = (vu + vd) / 2;
                    if (sample < 4 * cxy && cxy < 4 * sample) sample += (cxy - (cvu + cvd) / 2) / 4;
                } else {
                    sample = (vu + hl + vd + hr) / 4;
                    if (sample < 4 * cxy && cxy < 4 * sample)
                        sample += (cxy - (chl + chr + cvu + cvd) / 4) / 8;
                }

                (*rgbImage)[y][x][green] = clamp(sample);
                hl = hr;
                chl = cxy;
                cxy = chr;
            }
        }
    }

    // get the constant component out of the reconstructed green pixels and add to it
    // the "high frequency" part of the corresponding observed color channel

    for (int y = 2; y < height - 2; y++) {
        int channel = (y & 1) == (r.y & 1) ? red : blue;
        int x0 = 2 + ((y & 1) == (g.y & 1) ? g.x + 1 : g.x);

        int xy = (*rgbImage)[y][x0][green];
        int hl = (*rgbImage)[y][x0 - 2][green];
        int ul = (*rgbImage)[y - 2][x0 - 2][green];
        int bl = (*rgbImage)[y + 2][x0 - 2][green];

        int cxy = (*rgbImage)[y][x0][channel];
        int chl = (*rgbImage)[y][x0 - 2][channel];
        int cul = (*rgbImage)[y - 2][x0 - 2][channel];
        int cbl = (*rgbImage)[y + 2][x0 - 2][channel];

        for (int x = 2; x < width - 2; x += 2) {
            int hr = (*rgbImage)[y][x + 2][green];
            int ur = (*rgbImage)[y - 2][x + 2][green];
            int br = (*rgbImage)[y + 2][x + 2][green];
            int vu = (*rgbImage)[y - 2][x][green];
            int vd = (*rgbImage)[y + 2][x][green];

            int chr = (*rgbImage)[y][x + 2][channel];
            int cur = (*rgbImage)[y - 2][x + 2][channel];
            int cbr = (*rgbImage)[y + 2][x + 2][channel];
            int cvu = (*rgbImage)[y - 2][x][channel];
            int cvd = (*rgbImage)[y + 2][x][channel];

            // Only work on the pixels that have a strong enough correlation between channels

            if (xy < 4 * cxy && cxy < 4 * xy) {
                int dh = xy - (hl + hr) / 2;
                int dv = xy - (vu + vd) / 2;
                int ne = xy - (ul + br) / 2;
                int nw = xy - (ur + bl) / 2;

                int cdh = cxy - (chl + chr) / 2;
                int cdv = cxy - (cvu + cvd) / 2;
                int cne = cxy - (cul + cbr) / 2;
                int cnw = cxy - (cur + cbl) / 2;

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

                if (minimumDirection != none && minimumGradient > xy / 4) {
                    int sample;
                    switch (minimumDirection) {
                        case horizontal:
                            sample = (xy + (hl + hr) / 2 + cdh) / 2;
                            break;
                        case vertical:
                            sample = (xy + (vu + vd) / 2 + cdv) / 2;
                            break;
                        case northEast:
                            sample = (xy + (ul + br) / 2 + cne) / 2;
                            break;
                        case northWest:
                            sample = (xy + (ur + bl) / 2 + cnw) / 2;
                            break;
                        case none:
                            // never happens, just make the compiler happy
                            sample = (*rgbImage)[y][x][green];
                            break;
                    }

                    (*rgbImage)[y][x][green] = clamp(sample);
                }
            }

            hl = xy;
            xy = hr;
            ul = vu;
            vu = ur;
            bl = vd;
            vd = br;
            chl = cxy;
            cxy = chr;
            cul = cvu;
            cvu = cur;
            cbl = cvd;
            cvd = cbr;
        }
    }
}

void interpolateRedBlue(gls::image<gls::rgb_pixel_16>* image, BayerPattern bayerPattern) {
    const int width = image->width;
    const int height = image->height;

    auto offsets = bayerOffsets[bayerPattern];

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            for (int i = 0; i < 2; i++) {
                const int channel = i == 0 ? red : blue;
                const gls::point c = offsets[channel];

                if (((x + c.x) & 1) != (c.x & 1) || ((y + c.y) & 1) != (c.y & 1)) {
                    int sample;
                    int cg = (*image)[y + c.y][x + c.x][green];

                    if (((x + c.x) & 1) != (c.x & 1) && ((y + c.y) & 1) != (c.y & 1)) {
                        // Pixel at color location
                        int gne = (*image)[y + c.y + 1][x + c.x - 1][green];
                        int gnw = (*image)[y + c.y + 1][x + c.x + 1][green];
                        int gsw = (*image)[y + c.y - 1][x + c.x + 1][green];
                        int gse = (*image)[y + c.y - 1][x + c.x - 1][green];

                        int cne = gne - (*image)[y + c.y + 1][x + c.x - 1][channel];
                        int cnw = gnw - (*image)[y + c.y + 1][x + c.x + 1][channel];
                        int csw = gsw - (*image)[y + c.y - 1][x + c.x + 1][channel];
                        int cse = gse - (*image)[y + c.y - 1][x + c.x - 1][channel];

                        sample = cg - (cne + csw + cnw + cse) / 4;
                    } else if (((x + c.x) & 1) == (c.x & 1) && ((y + c.y) & 1) != (c.y & 1)) {
                        // Pixel at green location - vertical
                        int gu = (*image)[y + c.y - 1][x + c.x][green];
                        int gd = (*image)[y + c.y + 1][x + c.x][green];

                        int cu = gu - (*image)[y + c.y - 1][x + c.x][channel];
                        int cd = gd - (*image)[y + c.y + 1][x + c.x][channel];

                        sample = cg - (cu + cd) / 2;
                    } else {
                        // Pixel at green location - horizontal
                        int gl = (*image)[y + c.y][x + c.x - 1][green];
                        int gr = (*image)[y + c.y][x + c.x + 1][green];

                        int cl = gl - (*image)[y + c.y][x + c.x - 1][channel];
                        int cr = gr - (*image)[y + c.y][x + c.x + 1][channel];

                        sample = cg - (cl + cr) / 2;
                    }

                    (*image)[y + c.y][x + c.x][channel] = clamp(sample);
                }
            }
        }
    }
}

// XYZ -> RGB Transform
const gls::Matrix<3, 3> xyz_rgb({
    { 0.4124564, 0.3575761, 0.1804375 },
    { 0.2126729, 0.7151522, 0.0721750 },
    { 0.0193339, 0.1191920, 0.9503041 }
});

gls::Matrix<3, 3> cam_xyz_coeff(gls::Vector<4>& pre_mul, const gls::Matrix<3, 3>& cam_xyz) {
    auto cam_rgb = cam_xyz * xyz_rgb;

    // Normalize cam_rgb so that cam_rgb * (1,1,1) == (1,1,1)
    auto norm = cam_rgb * gls::Vector<3>({ 1, 1, 1 });
    for (int i = 0; i < 3; i++) {
        if (norm[i] > 0.00001) {
            cam_rgb[i] = cam_rgb[i] / norm[i];
            pre_mul[i] = 1 / norm[i];
        } else {
            cam_rgb[i] = gls::Vector<3>({ 0, 0, 0 });
            pre_mul[i] = 1;
        }
    }

    return inverse(cam_rgb);
}

template <typename T>
std::vector<T> getVector(const gls::tiff_metadata& metadata, ttag_t key) {
    const auto& entry = metadata.find(key);
    if (entry != metadata.end()) {
        return std::get<std::vector<T>>(metadata.find(key)->second);
    }
    return std::vector<T>();
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

    gls::Vector<3> cam_mul;
    for (int i = 0; i < 3; i++) {
        cam_mul[i] = 1.0 / as_shot_neutral[i];
    }

    gls::Matrix<3, 3> cam_xyz;
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            // TODO: this should be CameraCalibration * ColorMatrix * AsShotWhite
            cam_xyz[j][i] = color_matrix[3 * j + i];
        }
    }

    std::cout << "cam_xyz:\n" << cam_xyz << std::endl;

    gls::Vector<4> pre_mul;
    const auto rgb_cam = cam_xyz_coeff(pre_mul, cam_xyz);

    std::cout << "*** pre_mul: " << pre_mul << std::endl;
    std::cout << "*** cam_mul: " << cam_mul << std::endl;

    // If cam_mul is available use that instead of pre_mul
    for (int i = 0; i < 3; i++) {
        pre_mul[i] = cam_mul[i];
    }
    pre_mul[3] = pre_mul[1];

    std::cout << "rgb_cam:\n" << rgb_cam << std::endl;
    std::cout << "--> pre_mul: " << pre_mul << std::endl;

//    {
//        const gls::Matrix<1, 4> d65_white = {0.95047f, 1.0f, 1.08883f, 1};
//        const gls::Matrix<4, 1> d50_white = {0.9642, 1.0000, 0.8249, 1};
//
//        float cct = XYZtoCorColorTemp(cam_mul_xyz);
//        printf("*** Correlated Color Temperature: %f\n", cct);
//    }

    // Scale Input Image
    auto minmax = std::minmax_element(std::begin(pre_mul), std::end(pre_mul));
    gls::Vector<4> scale_mul;
    for (int c = 0; c < 4; c++) {
        printf("pre_mul[c]: %f, *minmax.second: %f, white_level: %d\n", pre_mul[c], *minmax.second, white_level);
        scale_mul[c] = (pre_mul[c] / *minmax.first) * 65535.0 / (white_level - black_level);
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
            const auto op = rgb_cam * gls::Vector<3>({ (float) p[0], (float) p[1], (float) p[2] });
            p = { clamp(op[0]), clamp(op[1]), clamp(op[2]) };
        }
    }

    return rgbImage;
}
