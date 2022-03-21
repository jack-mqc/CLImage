//
//  demosaic.cpp
//  CLImageTest
//
//  Created by Fabio Riccardi on 3/18/22.
//

#include <limits.h>
#include <math.h>

#include "demosaic.hpp"

enum { red = 0, green = 1, blue = 2 };

inline uint16_t clamp(int x) { return x < 0 ? 0 : x > 0xffff ? 0xffff : x; }

void interpolateGreen(const gls::image<gls::luma_pixel_16>& rawImage,
                      gls::image<gls::rgb_pixel_16>* rgbImage, int gx, int gy, int ry) {
    const int width = rawImage.width;
    const int height = rawImage.height;

    // copy RAW data to RGB layer and remove hot pixels
    for (int y = 0; y < height; y++) {
        int color = (y & 1) == (ry & 1) ? red : blue;
        int x0 = (y & 1) == (gy & 1) ? gx + 1 : gx;
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
        int color = (y & 1) == (ry & 1) ? red : blue;
        int x0 = (y & 1) == (gy & 1) ? gx + 1 : gx;

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
        int color = (y & 1) == (ry & 1) ? red : blue;
        int x0 = 2 + ((y & 1) == (gy & 1) ? gx + 1 : gx);

        int xy = (*rgbImage)[y][x0][green];
        int hl = (*rgbImage)[y][x0 - 2][green];
        int ul = (*rgbImage)[y - 2][x0 - 2][green];
        int bl = (*rgbImage)[y + 2][x0 - 2][green];

        int cxy = (*rgbImage)[y][x0][color];
        int chl = (*rgbImage)[y][x0 - 2][color];
        int cul = (*rgbImage)[y - 2][x0 - 2][color];
        int cbl = (*rgbImage)[y + 2][x0 - 2][color];

        for (int x = 2; x < width - 2; x += 2) {
            int hr = (*rgbImage)[y][x + 2][green];
            int ur = (*rgbImage)[y - 2][x + 2][green];
            int br = (*rgbImage)[y + 2][x + 2][green];
            int vu = (*rgbImage)[y - 2][x][green];
            int vd = (*rgbImage)[y + 2][x][green];

            int chr = (*rgbImage)[y][x + 2][color];
            int cur = (*rgbImage)[y - 2][x + 2][color];
            int cbr = (*rgbImage)[y + 2][x + 2][color];
            int cvu = (*rgbImage)[y - 2][x][color];
            int cvd = (*rgbImage)[y + 2][x][color];

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

                int gradients[4] = {
                    abs(dh) + abs(cdh), // horizontal
                    abs(dv) + abs(cdv), // vertical
                    abs(ne) + abs(cne), // north-east
                    abs(nw) + abs(cnw)  // north-west
                };

                int mind = 4, maxd = 4;
                int ming = INT_MAX;
                int maxg = 0;
                for (int i = 0; i < 4; i++) {
                    if (gradients[i] < ming) {
                        ming = gradients[i];
                        mind = i;
                    }
                    if (gradients[i] > maxg) {
                        maxg = gradients[i];
                        maxd = i;
                    }
                }

                // Only work on parts of the image that have enough "detail"

                if (mind != 4 && ming > xy / 4) {
                    int sample;
                    switch (mind) {
                        case 0:  // horizontal
                            sample = (xy + (hl + hr) / 2 + cdh) / 2;
                            break;
                        case 1:  // vertical
                            sample = (xy + (vu + vd) / 2 + cdv) / 2;
                            break;
                        case 2:  // north-east
                            sample = (xy + (ul + br) / 2 + cne) / 2;
                            break;
                        case 3:  // north-west
                            sample = (xy + (ur + bl) / 2 + cnw) / 2;
                            break;
                        case 4:  // flat
                            // nothing to do
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

void interpolateRedBlue(gls::image<gls::rgb_pixel_16>* image, int rx0, int ry0, int bx0, int by0) {
    const int width = image->width;
    const int height = image->height;

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            for (int i = 0; i < 2; i++) {
                int cx0, cy0, color;
                if (i == 0) {
                    cx0 = rx0;
                    cy0 = ry0;
                    color = red;
                } else {
                    cx0 = bx0;
                    cy0 = by0;
                    color = blue;
                }

                if (((x + cx0) & 1) != (cx0 & 1) || ((y + cy0) & 1) != (cy0 & 1)) {
                    int sample;
                    int cg = (*image)[y + cy0][x + cx0][green];

                    if (((x + cx0) & 1) != (cx0 & 1) && ((y + cy0) & 1) != (cy0 & 1)) {
                        // Pixel at color location
                        int gne = (*image)[y + cy0 + 1][x + cx0 - 1][green];
                        int gnw = (*image)[y + cy0 + 1][x + cx0 + 1][green];
                        int gsw = (*image)[y + cy0 - 1][x + cx0 + 1][green];
                        int gse = (*image)[y + cy0 - 1][x + cx0 - 1][green];

                        int cne = gne - (*image)[y + cy0 + 1][x + cx0 - 1][color];
                        int cnw = gnw - (*image)[y + cy0 + 1][x + cx0 + 1][color];
                        int csw = gsw - (*image)[y + cy0 - 1][x + cx0 + 1][color];
                        int cse = gse - (*image)[y + cy0 - 1][x + cx0 - 1][color];

                        sample = cg - (cne + csw + cnw + cse) / 4;
                    } else if (((x + cx0) & 1) == (cx0 & 1) && ((y + cy0) & 1) != (cy0 & 1)) {
                        // Pixel at green location - vertical
                        int gu = (*image)[y + cy0 - 1][x + cx0][green];
                        int gd = (*image)[y + cy0 + 1][x + cx0][green];

                        int cu = gu - (*image)[y + cy0 - 1][x + cx0][color];
                        int cd = gd - (*image)[y + cy0 + 1][x + cx0][color];

                        sample = cg - (cu + cd) / 2;
                    } else {
                        // Pixel at green location - horizontal
                        int gl = (*image)[y + cy0][x + cx0 - 1][green];
                        int gr = (*image)[y + cy0][x + cx0 + 1][green];

                        int cl = gl - (*image)[y + cy0][x + cx0 - 1][color];
                        int cr = gr - (*image)[y + cy0][x + cx0 + 1][color];

                        sample = cg - (cl + cr) / 2;
                    }

                    (*image)[y + cy0][x + cx0][color] = clamp(sample);
                }
            }
        }
    }
}
