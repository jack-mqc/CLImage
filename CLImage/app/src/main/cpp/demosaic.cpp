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

static const std::array<std::array<gls::point, 4>, 4> bayerOffsets = {
    std::array<gls::point, 4> { gls::point{1, 0}, gls::point{0, 0}, gls::point{0, 1}, gls::point{1, 1} }, // grbg
    std::array<gls::point, 4> { gls::point{0, 1}, gls::point{0, 0}, gls::point{1, 0}, gls::point{1, 1} }, // gbrg
    std::array<gls::point, 4> { gls::point{0, 0}, gls::point{1, 0}, gls::point{1, 1}, gls::point{0, 1} }, // rggb
    std::array<gls::point, 4> { gls::point{1, 1}, gls::point{1, 0}, gls::point{0, 0}, gls::point{0, 1} }  // bggr
};

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
const float xyz_rgb[3][3] = {
    {0.4124564, 0.3575761, 0.1804375},
    {0.2126729, 0.7151522, 0.0721750},
    {0.0193339, 0.1191920, 0.9503041}
};

void pseudoinverse(const float (*in)[3], float (*out)[3], int size) {
    double work[3][6], num;
    int i, j, k;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 6; j++) {
            work[i][j] = j == i + 3;
        }
        for (j = 0; j < 3; j++) {
            for (k = 0; k < size && k < 4; k++) {
                work[i][j] += in[k][i] * in[k][j];
            }
        }
    }
    for (i = 0; i < 3; i++) {
        num = work[i][i];
        for (j = 0; j < 6; j++) {
            if (fabs(num) > 0.00001f) {
                work[i][j] /= num;
            }
        }
        for (k = 0; k < 3; k++) {
            if (k == i)
                continue;
            num = work[k][i];
            for (j = 0; j < 6; j++) {
                work[k][j] -= work[i][j] * num;
            }
        }
    }
    for (i = 0; i < size && i < 4; i++) {
        for (j = 0; j < 3; j++) {
            for (out[i][j] = k = 0; k < 3; k++) {
                out[i][j] += work[j][k + 3] * in[i][k];
            }
        }
    }
}

void cam_xyz_coeff(float rgb_cam[3][3], float pre_mul[3], const float cam_xyz[3][3]) {
    float cam_rgb[3][3];
    for (int i = 0; i < 3; i++) /* Multiply out XYZ colorspace */
        for (int j = 0; j < 3; j++) {
            cam_rgb[i][j] = 0;
            for (int k = 0; k < 3; k++)
                cam_rgb[i][j] += cam_xyz[i][k] * xyz_rgb[k][j];
        }

    for (int i = 0; i < 3; i++) { /* Normalize cam_rgb so that */
        double num = 0;
        for (int j = 0; j < 3; j++) /* cam_rgb * (1,1,1) is (1,1,1,1) */
            num += cam_rgb[i][j];
        if (num > 0.00001) {
            for (int j = 0; j < 3; j++) {
                cam_rgb[i][j] /= num;
            }
            pre_mul[i] = 1 / num;
        } else {
            for (int j = 0; j < 3; j++) {
                cam_rgb[i][j] = 0.0;
            }
            pre_mul[i] = 1.0;
        }
    }

    float inverse[3][3];
    pseudoinverse(cam_rgb, inverse, 3);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rgb_cam[i][j] = inverse[j][i];
        }
    }
}

template <typename T>
std::vector<T> getVector(const gls::tiff_metadata& metadata, const std::string& key) {
    const auto& entry = metadata.find(key);
    if (entry != metadata.end()) {
        return std::get<std::vector<T>>(metadata.find(key)->second);
    }
    return std::vector<T>();
}

gls::image<gls::rgb_pixel_16>::unique_ptr demosaicImage(const gls::image<gls::luma_pixel_16>& rawImage,
                                                        const gls::tiff_metadata& metadata) {
    const auto color_matrix1 = getVector<float>(metadata, "ColorMatrix1");
    const auto color_matrix2 = getVector<float>(metadata, "ColorMatrix2");

    // If present ColorMatrix2 is usually D65 and ColorMatrix1 is Standard Light A
    const auto& color_matrix = color_matrix2.empty() ? color_matrix1 : color_matrix2;

    const auto as_shot_neutral = getVector<float>(metadata, "AsShotNeutral");
    const auto black_level_vec = getVector<float>(metadata, "BlackLevel");
    const auto white_level_vec = getVector<uint32_t>(metadata, "WhiteLevel");
    const auto cfa_pattern = getVector<uint8_t>(metadata, "CFAPattern");

    const float black_level = black_level_vec.empty() ? 0 : black_level_vec[0];
    const uint32_t white_level = white_level_vec.empty() ? 0xffff : white_level_vec[0];

    const auto bayerPattern = std::memcmp(cfa_pattern.data(), "\00\01\01\02", 4) == 0 ? BayerPattern::rggb
                            : std::memcmp(cfa_pattern.data(), "\02\01\01\00", 4) == 0 ? BayerPattern::bggr
                            : std::memcmp(cfa_pattern.data(), "\01\00\02\01", 4) == 0 ? BayerPattern::grbg
                            : BayerPattern::gbrg;

    float cam_mul[3];
    for (int i = 0; i < 3; i++) {
        cam_mul[i] = 1.0 / as_shot_neutral[i];
    }

    float cam_xyz[3][3];
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            // TODO: this should be CameraCalibration * ColorMatrix * AsShotWhite
            cam_xyz[j][i] = color_matrix[3 * j + i];
        }
    }

    float rgb_cam[3][3];
    float pre_mul[4];
    cam_xyz_coeff(rgb_cam, pre_mul, cam_xyz);

    printf("*** pre_mul: %f, %f, %f\n", pre_mul[0], pre_mul[1], pre_mul[2]);
    printf("*** cam_mul: %f, %f, %f\n", cam_mul[0], cam_mul[1], cam_mul[2]);

    // If cam_mul is available use that instead of pre_mul
    for (int i = 0; i < 3; i++) {
        pre_mul[i] = cam_mul[i];
    }
    pre_mul[3] = pre_mul[1];

    printf("rgb_cam: \n");
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            printf("%f, ", rgb_cam[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    printf("pre_mul: \n");
    for (int i = 0; i < 4; i++) {
        printf("%f, ", pre_mul[i]);
    }
    printf("\n");

//    {
//        const gls::Matrix<1, 4> d65_white = {0.95047f, 1.0f, 1.08883f, 1};
//        const gls::Matrix<4, 1> d50_white = {0.9642, 1.0000, 0.8249, 1};
//
//        gls::Matrix<4, 4> M = {
//            5,  -2,  2,  7,
//            1,   0,  0,  3,
//           -3,   1,  5,  0,
//            3,  -1, -9,  4,
//        };
//
//        gls::Matrix<4, 4> I = {
//            1, 0, 0, 0,
//            0, 1, 0, 0,
//            0, 0, 1, 0,
//            0, 0, 0, 1
//        };
//
//        const auto R = d65_white * (M + I) * d50_white;
//
//        print("M", R);
//
//        print("cofactor", gls::cofactor(M, 0, 0));
//
//        printf("determinant(M): %f\n", gls::determinant(M));
//
//        print("adjoint(M)", gls::adjoint(M));
//
//        print("inverse(M)", gls::inverse(M));
//
//        print("M / I", M / I);
//
//        print("I * 2", I * 2);
//
//        float cct = XYZtoCorColorTemp(cam_mul_xyz);
//        printf("*** Correlated Color Temperature: %f\n", cct);
//    }

    // Scale Input Image
    auto minmax = std::minmax_element(std::begin(pre_mul), std::end(pre_mul));
    std::array<float, 4> scale_mul;
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
            p = {
                clamp(p[0] * rgb_cam[0][0] + p[1] * rgb_cam[0][1] + p[2] * rgb_cam[0][2]),
                clamp(p[0] * rgb_cam[1][0] + p[1] * rgb_cam[1][1] + p[2] * rgb_cam[1][2]),
                clamp(p[0] * rgb_cam[2][0] + p[1] * rgb_cam[2][1] + p[2] * rgb_cam[2][2])
            };
        }
    }

    return rgbImage;
}

#if 0
int main(int argc, const char* argv[]) {
    printf("Hello CLImage!\n");

    if (argc > 1) {
        auto input_path = std::filesystem::path(argv[1]);
        auto input_dir = input_path.parent_path();

        std::cout << "Processing: " << input_path.filename() << std::endl;

        // Read the input file into an image object
        auto inputImage = gls::image<gls::luma_pixel_16>::read_png_file(input_path.string());

        const auto rgb_image = demosaicImage(*inputImage);

        rgb_image->write_png_file(input_path.replace_extension("_rgb.png").c_str());

        std::cout << "done with inputImage size: " << inputImage->width << " x " << inputImage->height << std::endl;

        auto output_file = input_path.replace_extension(".dng").c_str();

        inputImage->write_dng_file(output_file);
    }
}
#endif
