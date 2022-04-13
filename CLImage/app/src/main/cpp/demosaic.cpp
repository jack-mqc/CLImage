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

#include <climits>
#include <cmath>
#include <sys/types.h>

#include "demosaic.hpp"
#include "gls_color_science.hpp"
#include "gls_linalg.hpp"

#include "gls_cl.hpp"
#include "gls_logging.h"
#include "gls_cl_image.hpp"

static const char* TAG = "CLImage Pipeline";

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

gls::rectangle alignToQuad(const gls::rectangle& rect) {
    gls::rectangle alignedRect = rect;
    if (alignedRect.y & 1) {
        alignedRect.y += 1;
        alignedRect.height -= 1;
    }
    if (alignedRect.height & 1) {
        alignedRect.height -= 1;
    }
    if (alignedRect.x & 1) {
        alignedRect.x += 1;
        alignedRect.width -= 1;
    }
    if (alignedRect.width & 1) {
        alignedRect.width -= 1;
    }
    return alignedRect;
}

void colorcheck(const gls::image<gls::luma_pixel_16>& rawImage, BayerPattern bayerPattern, uint32_t black, std::array<gls::rectangle, 24> gmb_samples) {
// ColorChecker Chart under 6500-kelvin illumination
  static gls::Matrix<gmb_samples.size(), 3> gmb_xyY = {
    { 0.400, 0.350, 10.1 },        // Dark Skin
    { 0.377, 0.345, 35.8 },        // Light Skin
    { 0.247, 0.251, 19.3 },        // Blue Sky
    { 0.337, 0.422, 13.3 },        // Foliage
    { 0.265, 0.240, 24.3 },        // Blue Flower
    { 0.261, 0.343, 43.1 },        // Bluish Green
    { 0.506, 0.407, 30.1 },        // Orange
    { 0.211, 0.175, 12.0 },        // Purplish Blue
    { 0.453, 0.306, 19.8 },        // Moderate Red
    { 0.285, 0.202, 6.6  },        // Purple
    { 0.380, 0.489, 44.3 },        // Yellow Green
    { 0.473, 0.438, 43.1 },        // Orange Yellow
    { 0.187, 0.129, 6.1  },        // Blue
    { 0.305, 0.478, 23.4 },        // Green
    { 0.539, 0.313, 12.0 },        // Red
    { 0.448, 0.470, 59.1 },        // Yellow
    { 0.364, 0.233, 19.8 },        // Magenta
    { 0.196, 0.252, 19.8 },        // Cyan
    { 0.310, 0.316, 90.0 },        // White
    { 0.310, 0.316, 59.1 },        // Neutral 8
    { 0.310, 0.316, 36.2 },        // Neutral 6.5
    { 0.310, 0.316, 19.8 },        // Neutral 5
    { 0.310, 0.316, 9.0 },         // Neutral 3.5
    { 0.310, 0.316, 3.1 } };       // Black

    const auto offsets = bayerOffsets[bayerPattern];

    auto* writeRawImage = (gls::image<gls::luma_pixel_16>*) &rawImage;

    gls::Matrix<gmb_samples.size(), 4> gmb_cam;
    gls::Matrix<gmb_samples.size(), 3> gmb_xyz;

    for (int sq = 0; sq < gmb_samples.size(); sq++) {
        std::array<int, 3> count { /* zero */ };
        auto patch = alignToQuad(gmb_samples[sq]);
        for (int y = patch.y; y < patch.y + patch.height; y += 2) {
            for (int x = patch.x; x < patch.x + patch.width; x += 2) {
                for (int c = 0; c < 4; c++) {
                    const auto& o = offsets[c];
                    int val = rawImage[y + o.y][x + o.x];
                    gmb_cam[sq][c == 3 ? 1 : c] += (float) val;
                    count[c == 3 ? 1 : c]++;
                    // Mark image to identify sampled areas
                    (*writeRawImage)[y + o.y][x + o.x] = black + (val - black) / 2;
                }
            }
        }

        for (int c = 0; c < 3; c++) {
            gmb_cam[sq][c] = gmb_cam[sq][c] / (float) count[c] - (float) black;
        }
        gmb_xyz[sq][0] = gmb_xyY[sq][2] * gmb_xyY[sq][0] / gmb_xyY[sq][1];
        gmb_xyz[sq][1] = gmb_xyY[sq][2];
        gmb_xyz[sq][2] = gmb_xyY[sq][2] * (1 - gmb_xyY[sq][0] - gmb_xyY[sq][1]) / gmb_xyY[sq][1];
    }

    gls::Matrix<gmb_samples.size(), 3> inverse = pseudoinverse(gmb_xyz);

    gls::Matrix<3, 3> cam_xyz;
    for (int pass=0; pass < 2; pass++) {
        for (int i = 0; i < 3 /* colors */; i++) {
            for (int j = 0; j < 3; j++) {
                cam_xyz[i][j] = 0;
                for (int k = 0; k < gmb_samples.size(); k++)
                    cam_xyz[i][j] += gmb_cam[k][i] * inverse[k][j];
            }
        }

        gls::Vector<3> pre_mul;
        cam_xyz_coeff(pre_mul, cam_xyz);

        gls::Vector<4> balance;
        for (int c = 0; c < 4; c++) {
            balance[c] = pre_mul[c == 3 ? 1 : c] * gmb_cam[20][c];
        }
        for (int sq = 0; sq < gmb_samples.size(); sq++) {
            for (int c = 0; c < 4; c++) {
                gmb_cam[sq][c] *= balance[c];
            }
        }
    }

    float norm = 1 / (cam_xyz[1][0] + cam_xyz[1][1] + cam_xyz[1][2]);
    printf("Color Matrix: ");
    for (int c = 0; c < 3; c++) {
        for (int j = 0; j < 3; j++)
            printf("%.4f, ", cam_xyz[c][j] * norm);
    }
    printf("\n");
}

void white_balance(const gls::image<gls::luma_pixel_16>& rawImage, gls::Vector<3>* wb_mul, uint32_t white, uint32_t black, BayerPattern bayerPattern) {
    const auto offsets = bayerOffsets[bayerPattern];

    std::array<float, 8> fsum { /* zero */ };
    for (int y = 0; y < rawImage.height/2; y += 8) {
        for (int x = 0; x < rawImage.width/2; x += 8) {
            std::array<uint32_t, 8> sum { /* zero */ };
            for (int j = y; j < 8 && j < rawImage.height/2; j++) {
                for (int i = x; i < 8 && i < rawImage.width/2; i++) {
                    for (int c = 0; c < 4; c++) {
                        const auto& o = offsets[c];
                        uint32_t val = rawImage[2 * j + o.y][2 * i + o.x];
                        if (val > white - 25) {
                            goto skip_block;
                        }
                        if ((val -= black) < 0) {
                            val = 0;
                        }
                        sum[c] += val;
                        sum[c+4]++;
                    }
                }
            }
            for (int i = 0; i < 8; i++) {
                fsum[i] += (float) sum[i];
            }
            skip_block:
            ;
        }
    }
    // Aggregate green2 data to green
    fsum[1] += fsum[3];
    fsum[5] += fsum[7];

    for (int c = 0; c < 3; c++) {
        if (fsum[c] != 0) {
            (*wb_mul)[c] = fsum[c+4] / fsum[c];
        }
    }
    // Normalize with green = 1
    *wb_mul = *wb_mul / (*wb_mul)[1];
}

// Coordinates of the GretagMacbeth ColorChecker squares
// x, y, width, height from 2022-04-07-17-21-33-023.png
//static std::array<gls::rectangle, 24> gmb_samples = {{
//    { 4729, 3390, 80, 80 },
//    { 4597, 3385, 66, 81 },
//    { 4459, 3380, 71, 74 },
//    { 4317, 3373, 70, 71 },
//    { 4172, 3366, 86, 78 },
//    { 4032, 3357, 86, 76 },
//    { 4742, 3256, 70, 71 },
//    { 4599, 3251, 81, 74 },
//    { 4462, 3244, 75, 73 },
//    { 4321, 3230, 75, 83 },
//    { 4184, 3223, 71, 78 },
//    { 4040, 3218, 80, 74 },
//    { 4751, 3118, 70, 73 },
//    { 4609, 3113, 72, 69 },
//    { 4469, 3106, 72, 72 },
//    { 4330, 3098, 75, 73 },
//    { 4189, 3097, 79, 62 },
//    { 4054, 3091, 69, 63 },
//    { 4757, 2979, 71, 74 },
//    { 4616, 2975, 65, 73 },
//    { 4475, 2965, 75, 74 },
//    { 4341, 2959, 67, 73 },
//    { 4200, 2952, 71, 72 },
//    { 4061, 2942, 72, 74 },
//}};

// Coordinates of the GretagMacbeth ColorChecker squares
// x, y, width, height from 2022-04-12-10-43-56-566.dng
static std::array<gls::rectangle, 24> gmb_samples = {{
    { 4886, 2882, 285, 273 },
    { 4505, 2899, 272, 235 },
    { 4122, 2892, 262, 240 },
    { 3742, 2900, 256, 225 },
    { 3352, 2897, 258, 227 },
    { 2946, 2904, 282, 231 },
    { 4900, 2526, 274, 244 },
    { 4513, 2526, 262, 237 },
    { 4133, 2529, 235, 227 },
    { 3733, 2523, 254, 237 },
    { 3347, 2530, 245, 234 },
    { 2932, 2531, 283, 233 },
    { 4899, 2151, 283, 252 },
    { 4519, 2155, 261, 245 },
    { 4119, 2157, 269, 245 },
    { 3737, 2160, 246, 226 },
    { 3335, 2168, 261, 239 },
    { 2957, 2183, 233, 214 },
    { 4923, 1784, 265, 243 },
    { 4531, 1801, 250, 219 },
    { 4137, 1792, 234, 226 },
    { 3729, 1790, 254, 230 },
    { 3337, 1793, 250, 232 },
    { 2917, 1800, 265, 228 },
}};


void unpackRawMetadata(const gls::image<gls::luma_pixel_16>& rawImage,
                       gls::tiff_metadata* metadata,
                       BayerPattern *bayerPattern,
                       float *black_level,
                       gls::Vector<4> *scale_mul,
                       gls::Matrix<3, 3> *rgb_cam,
                       bool auto_white_balance) {
    const auto color_matrix1 = getVector<float>(*metadata, TIFFTAG_COLORMATRIX1);
    const auto color_matrix2 = getVector<float>(*metadata, TIFFTAG_COLORMATRIX2);

    // If present ColorMatrix2 is usually D65 and ColorMatrix1 is Standard Light A
    const auto& color_matrix = color_matrix2.empty() ? color_matrix1 : color_matrix2;

    auto as_shot_neutral = getVector<float>(*metadata, TIFFTAG_ASSHOTNEUTRAL);
    const auto black_level_vec = getVector<float>(*metadata, TIFFTAG_BLACKLEVEL);
    const auto white_level_vec = getVector<uint32_t>(*metadata, TIFFTAG_WHITELEVEL);
    const auto cfa_pattern = getVector<uint8_t>(*metadata, TIFFTAG_CFAPATTERN);

    *black_level = black_level_vec.empty() ? 0 : black_level_vec[0];
    const uint32_t white_level = white_level_vec.empty() ? 0xffff : white_level_vec[0];

    *bayerPattern = std::memcmp(cfa_pattern.data(), "\00\01\01\02", 4) == 0 ? BayerPattern::rggb
                  : std::memcmp(cfa_pattern.data(), "\02\01\01\00", 4) == 0 ? BayerPattern::bggr
                  : std::memcmp(cfa_pattern.data(), "\01\00\02\01", 4) == 0 ? BayerPattern::grbg
                  : BayerPattern::gbrg;

    std::cout << "as_shot_neutral: " << gls::Vector<3>(as_shot_neutral) << std::endl;

    // Uncomment this to characterize sensor
    // colorcheck(rawImage, *bayerPattern, *black_level, gmb_samples);

    gls::Vector<3> cam_mul = 1.0 / gls::Vector<3>(as_shot_neutral);

    // TODO: this should be CameraCalibration * ColorMatrix * AsShotWhite
    gls::Matrix<3, 3> cam_xyz = color_matrix;

    std::cout << "cam_xyz:\n" << cam_xyz << std::endl;

    gls::Vector<3> pre_mul;
    *rgb_cam = cam_xyz_coeff(pre_mul, cam_xyz);

    std::cout << "*** pre_mul: " << pre_mul << std::endl;
    std::cout << "*** cam_mul: " << cam_mul << std::endl;

    // Save the whitening transformation
    const auto inv_cam_white = pre_mul;

    if (auto_white_balance) {
        white_balance(rawImage, &cam_mul, white_level, *black_level, *bayerPattern);

        printf("Auto White Balance: %f, %f, %f\n", cam_mul[0], cam_mul[1], cam_mul[2]);

        // Convert cam_mul from camera to XYZ
        const auto cam_mul_xyz = cam_xyz * cam_mul;
        std::cout << "cam_mul_xyz: " << cam_mul_xyz << ", CCT: " << XYZtoCorColorTemp(cam_mul_xyz) << std::endl;

        for (int c = 0; c < 3; c++) {
            as_shot_neutral[c] = 1 / cam_mul[c];
        }
        (*metadata)[TIFFTAG_ASSHOTNEUTRAL] = as_shot_neutral;
    }

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

        auto cam_white = inverse(cam_xyz) * xyz_rgb * gls::Vector<3>({ 1, 1, 1 });
        std::cout << "inverse(cam_xyz) * cam_white: " << cam_xyz * cam_white << ", CCT: " << XYZtoCorColorTemp(cam_xyz * cam_white) << std::endl;
        std::cout << "xyz_rgb * gls::Vector<3>({ 1, 1, 1 }): " << xyz_rgb * gls::Vector<3>({ 1, 1, 1 }) << ", CCT: " << XYZtoCorColorTemp(xyz_rgb * gls::Vector<3>({ 1, 1, 1 })) << std::endl;

        const auto wb_out = xyz_rgb * *rgb_cam * mCamMul * gls::Vector({1, 1, 1});
        std::cout << "wb_out: " << wb_out << ", CCT: " << XYZtoCorColorTemp(wb_out) << std::endl;

        const auto no_wb_out = xyz_rgb * *rgb_cam * (1 / inv_cam_white);
        std::cout << "no_wb_out: " << wb_out << ", CCT: " << XYZtoCorColorTemp(no_wb_out) << std::endl;
    }

    // Scale Input Image
    auto minmax = std::minmax_element(std::begin(pre_mul), std::end(pre_mul));
    for (int c = 0; c < 4; c++) {
        int pre_mul_idx = c == 3 ? 1 : c;
        printf("pre_mul[c]: %f, *minmax.second: %f, white_level: %d\n", pre_mul[pre_mul_idx], *minmax.second, white_level);
        (*scale_mul)[c] = (pre_mul[pre_mul_idx] / *minmax.first) * 65535.0 / (white_level - *black_level);
    }
    printf("scale_mul: %f, %f, %f, %f\n", (*scale_mul)[0], (*scale_mul)[1], (*scale_mul)[2], (*scale_mul)[3]);
}

gls::image<gls::rgb_pixel_16>::unique_ptr demosaicImage(const gls::image<gls::luma_pixel_16>& rawImage,
                                                        gls::tiff_metadata* metadata, bool auto_white_balance) {
    BayerPattern bayerPattern;
    float black_level;
    gls::Vector<4> scale_mul;
    gls::Matrix<3, 3> rgb_cam;

    unpackRawMetadata(rawImage, metadata, &bayerPattern, &black_level, &scale_mul, &rgb_cam, auto_white_balance);

    LOG_INFO(TAG) << "Begin demosaicing image (CPU)..." << std::endl;

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

int scaleRawData(gls::OpenCLContext* glsContext,
                 const gls::cl_image_2d<gls::luma_pixel_16>& rawImage,
                 gls::cl_image_2d<gls::luma_pixel_float>* scaledRawImage,
                 BayerPattern bayerPattern, gls::Vector<4> scaleMul, float blackLevel) {
    try {
        // Load the shader source
        const auto demosaicProgram = glsContext->loadProgram("demosaic");

        // Bind the kernel parameters
        auto kernel = cl::KernelFunctor<cl::Image2D,  // rawImage
                                        cl::Image2D,  // scaledRawImage
                                        int,          // bayerPattern
                                        cl::Buffer,   // scaleMul
                                        float         // blackLevel
                                        >(demosaicProgram, "scaleRawData");

        cl::Buffer scaleMulBuffer(scaleMul.begin(), scaleMul.end(), true);

        // Work on one Quad (2x2) at a time
        kernel(gls::OpenCLContext::buildEnqueueArgs(scaledRawImage->width/2, scaledRawImage->height/2),
               rawImage.getImage2D(), scaledRawImage->getImage2D(), bayerPattern, scaleMulBuffer, blackLevel);
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}

int interpolateGreen(gls::OpenCLContext* glsContext,
                     const gls::cl_image_2d<gls::luma_pixel_float>& rawImage,
                     gls::cl_image_2d<gls::luma_pixel_float>* greenImage,
                     BayerPattern bayerPattern) {
    try {
        // Load the shader source
        const auto demosaicProgram = glsContext->loadProgram("demosaic");

        // Bind the kernel parameters
        auto kernel = cl::KernelFunctor<cl::Image2D,  // rawImage
                                        cl::Image2D,  // greenImage
                                        int           // bayerPattern
                                        >(demosaicProgram, "interpolateGreen");

        // Schedule the kernel on the GPU
        kernel(gls::OpenCLContext::buildEnqueueArgs(greenImage->width, greenImage->height),
               rawImage.getImage2D(), greenImage->getImage2D(), bayerPattern);
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}

int interpolateRedBlue(gls::OpenCLContext* glsContext,
                       const gls::cl_image_2d<gls::luma_pixel_float>& rawImage,
                       const gls::cl_image_2d<gls::luma_pixel_float>& greenImage,
                       gls::cl_image_2d<gls::rgba_pixel_float>* rgbImage,
                       BayerPattern bayerPattern) {
    try {
        // Load the shader source
        const auto demosaicProgram = glsContext->loadProgram("demosaic");

        // Bind the kernel parameters
        auto kernel = cl::KernelFunctor<cl::Image2D,  // rawImage
                                        cl::Image2D,  // greenImage
                                        cl::Image2D,  // rgbImage
                                        int           // bayerPattern
                                        >(demosaicProgram, "interpolateRedBlue");

        // Schedule the kernel on the GPU
        kernel(gls::OpenCLContext::buildEnqueueArgs(rgbImage->width, rgbImage->height),
               rawImage.getImage2D(), greenImage.getImage2D(), rgbImage->getImage2D(), bayerPattern);
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}

int applyKernel(gls::OpenCLContext* glsContext, const std::string& filterName,
                const gls::cl_image_2d<gls::rgba_pixel_float>& inputImage,
                gls::cl_image_2d<gls::rgba_pixel_float>* outputImage) {
    try {
        // Load the shader source
        const auto program = glsContext->loadProgram("demosaic");

        // Bind the kernel parameters
        auto kernel = cl::KernelFunctor<cl::Image2D,  // inputImage
                                        cl::Image2D   // outputImage
                                        >(program, filterName);

        // Schedule the kernel on the GPU
        kernel(gls::OpenCLContext::buildEnqueueArgs(outputImage->width, outputImage->height),
               inputImage.getImage2D(), outputImage->getImage2D());
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}

int convertTosRGB(gls::OpenCLContext* glsContext,
                  const gls::cl_image_2d<gls::rgba_pixel_float>& linearImage,
                  gls::cl_image_2d<gls::rgba_pixel>* rgbImage,
                  const gls::Matrix<3, 3>& transform) {
    try {
        // Load the shader source
        const auto demosaicProgram = glsContext->loadProgram("demosaic");

        // float4 transform[3];
        gls::Matrix<3, 4> paddedTransform;
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                paddedTransform[r][c] = transform[r][c];
            }
        }

        cl::Buffer transformBuffer(paddedTransform.span().begin(), paddedTransform.span().end(), true);

        // Bind the kernel parameters
        auto kernel = cl::KernelFunctor<cl::Image2D,  // linearImage
                                        cl::Image2D,  // rgbImage
                                        cl::Buffer    // transform
                                        >(demosaicProgram, "convertTosRGB");

        // Schedule the kernel on the GPU
        kernel(gls::OpenCLContext::buildEnqueueArgs(rgbImage->width, rgbImage->height),
               linearImage.getImage2D(), rgbImage->getImage2D(), transformBuffer);
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}

gls::image<gls::rgba_pixel>::unique_ptr demosaicImageGPU(const gls::image<gls::luma_pixel_16>& rawImage,
                                                         gls::tiff_metadata* metadata, bool auto_white_balance) {
    BayerPattern bayerPattern;
    float black_level;
    gls::Vector<4> scale_mul;
    gls::Matrix<3, 3> rgb_cam;

    unpackRawMetadata(rawImage, metadata, &bayerPattern, &black_level, &scale_mul, &rgb_cam, auto_white_balance);

    gls::OpenCLContext glsContext("");
    auto clContext = glsContext.clContext();

    LOG_INFO(TAG) << "Begin demosaicing image (GPU)..." << std::endl;

    gls::cl_image_2d<gls::luma_pixel_16> clRawImage(clContext, rawImage);
    gls::cl_image_2d<gls::luma_pixel_float> clScaledRawImage(clContext, rawImage.width, rawImage.height);

    scaleRawData(&glsContext, clRawImage, &clScaledRawImage, bayerPattern, scale_mul, black_level / 0xffff);

    gls::cl_image_2d<gls::luma_pixel_float> clGreenImage(clContext, rawImage.width, rawImage.height);

    interpolateGreen(&glsContext, clScaledRawImage, &clGreenImage, bayerPattern);

    gls::cl_image_2d<gls::rgba_pixel_float> clLinearRGBImage(clContext, rawImage.width, rawImage.height);

    interpolateRedBlue(&glsContext, clScaledRawImage, clGreenImage, &clLinearRGBImage, bayerPattern);

    applyKernel(&glsContext, "rgbToYCbCrImage", clLinearRGBImage, &clLinearRGBImage);

    gls::cl_image_2d<gls::rgba_pixel_float> clDenoisedRGBImage(clContext, rawImage.width, rawImage.height);
    applyKernel(&glsContext, "denoiseImage", clLinearRGBImage, &clDenoisedRGBImage);

    // applyKernel(&glsContext, "sharpenLumaImage", clDenoisedRGBImage, &clLinearRGBImage);

    applyKernel(&glsContext, "yCbCrtoRGBImage", clDenoisedRGBImage, &clDenoisedRGBImage);

//    gls::cl_image_2d<gls::rgba_pixel_16> clSharpenedRGBImage(clContext, rawImage.width, rawImage.height);
//    sharpenImage(&glsContext, clDenoisedRGBImage, &clSharpenedRGBImage);

    gls::cl_image_2d<gls::rgba_pixel> clsRGBImage(clContext, rawImage.width, rawImage.height);
    convertTosRGB(&glsContext, clDenoisedRGBImage, &clsRGBImage, rgb_cam);

    auto rgbaImage = clsRGBImage.toImage();

    LOG_INFO(TAG) << "...done with demosaicing (GPU)." << std::endl;

    return rgbaImage;
}
