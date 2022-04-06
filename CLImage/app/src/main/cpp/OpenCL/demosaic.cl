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

// const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

enum BayerPattern {
    grbg = 0,
    gbrg = 1,
    rggb = 2,
    bggr = 3
};

enum { raw_red = 0, raw_green = 1, raw_blue = 2, raw_green2 = 3 };

constant const int2 bayerOffsets[4][4] = {
    { {1, 0}, {0, 0}, {0, 1}, {1, 1} }, // grbg
    { {0, 1}, {0, 0}, {1, 0}, {1, 1} }, // gbrg
    { {0, 0}, {1, 0}, {1, 1}, {0, 1} }, // rggb
    { {1, 1}, {1, 0}, {0, 0}, {0, 1} }  // bggr
};

// Work on one Quad (2x2) at a time
kernel void scaleRawData(read_only image2d_t rawImage, write_only image2d_t scaledRawImage,
                         int bayerPattern, constant float scaleMul[4], float blackLevel) {
    const int2 imageCoordinates = (int2) (2 * get_global_id(0), 2 * get_global_id(1));
    for (int c = 0; c < 4; c++) {
        int2 o = bayerOffsets[bayerPattern][c];
        write_imagef(scaledRawImage, imageCoordinates + (int2) (o.x, o.y),
                     clamp(scaleMul[c] * (read_imagef(rawImage, imageCoordinates + (int2) (o.x, o.y)).x - blackLevel), 0.0f, 1.0f));
    }
}

kernel void interpolateGreen(read_only image2d_t rawImage, write_only image2d_t greenImage, int bayerPattern) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    const int x = imageCoordinates.x;
    const int y = imageCoordinates.y;

    const int2 g = bayerOffsets[bayerPattern][raw_green];
    int x0 = (y & 1) == (g.y & 1) ? g.x + 1 : g.x;

    if ((x0 & 1) == (x & 1)) {
        float g_left  = read_imagef(rawImage, (int2)(x - 1, y)).x;
        float g_right = read_imagef(rawImage, (int2)(x + 1, y)).x;
        float g_up    = read_imagef(rawImage, (int2)(x, y - 1)).x;
        float g_down  = read_imagef(rawImage, (int2)(x, y + 1)).x;
        float g_dh    = fabs(g_left - g_right);
        float g_dv    = fabs(g_up - g_down);

        float c_xy    = read_imagef(rawImage, (int2)(x, y)).x;

        float c_left  = read_imagef(rawImage, (int2)(x - 2, y)).x;
        float c_right = read_imagef(rawImage, (int2)(x + 2, y)).x;
        float c_up    = read_imagef(rawImage, (int2)(x, y - 2)).x;
        float c_down  = read_imagef(rawImage, (int2)(x, y + 2)).x;
        float c_dh    = fabs(c_left + c_right - 2 * c_xy);
        float c_dv    = fabs(c_up + c_down - 2 * c_xy);

        // Minimum derivative value for edge directed interpolation (avoid aliasing)
        float dThreshold = 0.0183; // TODO: use noise model

        // we're doing edge directed bilinear interpolation on the green channel,
        // which is a low pass operation (averaging), so we add some signal from the
        // high frequencies of the observed color channel

        float sample;
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

        write_imagef(greenImage, imageCoordinates, clamp(sample, 0.0f, 1.0f));
    } else {
        write_imagef(greenImage, imageCoordinates, read_imagef(rawImage, (int2)(x, y)).x);
    }
}

kernel void interpolateRedBlue(read_only image2d_t rawImage, read_only image2d_t greenImage,
                               write_only image2d_t rgbImage, int bayerPattern) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    const int x = imageCoordinates.x;
    const int y = imageCoordinates.y;

    const int2 r = bayerOffsets[bayerPattern][raw_red];
    const int2 g = bayerOffsets[bayerPattern][raw_green];
    const int2 b = bayerOffsets[bayerPattern][raw_blue];

    int color = (r.x & 1) == (x & 1) && (r.y & 1) == (y & 1) ? raw_red :
                (g.x & 1) == (x & 1) && (g.y & 1) == (y & 1) ? raw_green :
                (b.x & 1) == (x & 1) && (b.y & 1) == (y & 1) ? raw_blue : raw_green2;

    float green = read_imagef(greenImage, imageCoordinates).x;
    float red;
    float blue;
    switch (color) {
        case raw_red:
        case raw_blue:
        {
            float c1 = read_imagef(rawImage, imageCoordinates).x;

            float g_top_left      = read_imagef(greenImage, (int2)(x - 1, y - 1)).x;
            float g_top_right     = read_imagef(greenImage, (int2)(x + 1, y - 1)).x;
            float g_bottom_left   = read_imagef(greenImage, (int2)(x - 1, y + 1)).x;
            float g_bottom_right  = read_imagef(greenImage, (int2)(x + 1, y + 1)).x;

            float c2_top_left     = g_top_left     - read_imagef(rawImage, (int2)(x - 1, y - 1)).x;
            float c2_top_right    = g_top_right    - read_imagef(rawImage, (int2)(x + 1, y - 1)).x;
            float c2_bottom_left  = g_bottom_left  - read_imagef(rawImage, (int2)(x - 1, y + 1)).x;
            float c2_bottom_right = g_bottom_right - read_imagef(rawImage, (int2)(x + 1, y + 1)).x;

            float d_ne_sw = fabs(c2_top_left - c2_bottom_right);
            float d_nw_se = fabs(c2_top_right - c2_bottom_left);

            // Minimum gradient for edge directed interpolation
            float dThreshold = 0.0122f; // TODO: use noise model
            float c2;
            if (d_ne_sw > dThreshold && d_ne_sw > d_nw_se) {
                c2 = green - (c2_top_right + c2_bottom_left) / 2;
            } else if (d_nw_se > dThreshold && d_nw_se > d_ne_sw) {
                c2 = green - (c2_top_left + c2_bottom_right) / 2;
            } else {
                c2 = green - (c2_top_left + c2_top_right + c2_bottom_left + c2_bottom_right) / 4;
            }

            if (color == raw_red) {
                red = c1;
                blue = c2;
            } else {
                blue = c1;
                red = c2;
            }
        }
        break;

        case raw_green:
        case raw_green2:
        {
            float g_left    = read_imagef(greenImage, (int2)(x - 1, y)).x;
            float g_right   = read_imagef(greenImage, (int2)(x + 1, y)).x;
            float g_up      = read_imagef(greenImage, (int2)(x, y - 1)).x;
            float g_down    = read_imagef(greenImage, (int2)(x, y + 1)).x;

            float c1_left   = g_left  - read_imagef(rawImage, (int2)(x - 1, y)).x;
            float c1_right  = g_right - read_imagef(rawImage, (int2)(x + 1, y)).x;
            float c2_up     = g_up    - read_imagef(rawImage, (int2)(x, y - 1)).x;
            float c2_down   = g_down  - read_imagef(rawImage, (int2)(x, y + 1)).x;

            float c1 = green - (c1_left + c1_right) / 2;
            float c2 = green - (c2_up + c2_down) / 2;

            if (color == raw_green2) {
                red = c1;
                blue = c2;
            } else {
                blue = c1;
                red = c2;
            }
        }
        break;
    }

    write_imagef(rgbImage, imageCoordinates, (float4)(clamp((float3)(red, green, blue), 0.0f, 1.0f), 0));
}

static float3 sigmoid(float3 x, float s) {
    return 0.5 * (tanh(s * x - 0.3 * s) + 1);
}

static float3 toneCurve(float3 x) {
    float s = 3.5;
    return (sigmoid(native_powr(x, 1/2.2), s) - sigmoid(0, s)) / (sigmoid(1, s) - sigmoid(0, s));
}

kernel void convertTosRGB(read_only image2d_t linearImage, write_only image2d_t rgbImage,
                          constant float3 transform[3]) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    float3 linear = read_imagef(linearImage, imageCoordinates).xyz;
    float3 rgb = (float3) (dot(transform[0], linear), dot(transform[1], linear), dot(transform[2], linear));

    write_imagef(rgbImage, imageCoordinates, (float4) (toneCurve(0.8 * clamp(rgb, 0.0f, 1.0f)), 0.0f));
}
