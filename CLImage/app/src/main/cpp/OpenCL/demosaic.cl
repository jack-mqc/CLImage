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

#if defined(__QCOMM_QGPU_A3X__) || \
    defined(__QCOMM_QGPU_A4X__) || \
    defined(__QCOMM_QGPU_A5X__) || \
    defined(__QCOMM_QGPU_A6X__) || \
    defined(__QCOMM_QGPU_A7V__) || \
    defined(__QCOMM_QGPU_A7P__)

#define _fabs(a) \
   ({ __typeof__ (a) _a = (a); \
     _a >= 0 ? _a : - _a; })

#define _fmax(a, b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define _fmin(a, b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define _clamp(x, minval, maxval) \
   ({ __typeof__ (x) _x = (x); \
      __typeof__ (minval) _minval = (minval); \
      __typeof__ (maxval) _maxval = (maxval); \
     _x < _minval ? _minval : _x > _maxval ? _maxval : _x; })

#define _mix(x, y, a) \
   ({ __typeof__ (x) _x = (x); \
      __typeof__ (y) _y = (y); \
      __typeof__ (a) _a = (a); \
      _x + (_y - _x) * _a; })

// Qualcomm's default implementation of smoothstep can be really slow...

#define _smoothstep(a, b, x) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      __typeof__ (x) _x = (x), t; \
      t = _clamp(_x * (1 / (_b - _a)) - (_a / (_b - _a)), 0, 1); t * t * (3 - 2 * t); })

#define fabs _fabs
#define fmax _fmax
#define fmin _fmin
#define clamp _clamp
#define mix _mix
#define smoothstep _smoothstep

#endif

// Work on one Quad (2x2) at a time
kernel void scaleRawData(read_only image2d_t rawImage, write_only image2d_t scaledRawImage,
                         int bayerPattern, constant float scaleMul[4], float blackLevel) {
    const int2 imageCoordinates = (int2) (2 * get_global_id(0), 2 * get_global_id(1));
    for (int c = 0; c < 4; c++) {
        int2 o = bayerOffsets[bayerPattern][c];
        write_imagef(scaledRawImage, imageCoordinates + (int2) (o.x, o.y),
                     clamp(scaleMul[c] * (read_imagef(rawImage, imageCoordinates + (int2) (o.x, o.y)).x - blackLevel), 0.0, 1.0));
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

        float c_xy    = read_imagef(rawImage, (int2)(x, y)).x;

        float c_left  = read_imagef(rawImage, (int2)(x - 2, y)).x;
        float c_right = read_imagef(rawImage, (int2)(x + 2, y)).x;
        float c_up    = read_imagef(rawImage, (int2)(x, y - 2)).x;
        float c_down  = read_imagef(rawImage, (int2)(x, y + 2)).x;

        float c2_top_left       = read_imagef(rawImage, (int2)(x - 1, y - 1)).x;
        float c2_top_right      = read_imagef(rawImage, (int2)(x + 1, y - 1)).x;
        float c2_bottom_left    = read_imagef(rawImage, (int2)(x - 1, y + 1)).x;
        float c2_bottom_right   = read_imagef(rawImage, (int2)(x + 1, y + 1)).x;
        float c2_ave = (c2_top_left + c2_top_right + c2_bottom_left + c2_bottom_right) / 4;

        float g_ave = (g_left + g_right + g_up + g_down) / 4;

        float2 dv = (float2) (fabs(g_left - g_right), fabs(g_up - g_down));
        dv += (float2) (fabs(c_left + c_right - 2 * c_xy) / 2, fabs(c_up + c_down - 2 * c_xy) / 2);

        // If the gradient estimation at this location is too weak, try it on a 3x3 patch
        if (length(dv) < 0.1) {
            dv = 0;
            for (int j = -1; j <= 1; j++) {
                for (int i = -1; i <= 1; i++) {
                    float v_left  = read_imagef(rawImage, (int2)(x+i - 1, j+y)).x;
                    float v_right = read_imagef(rawImage, (int2)(x+i + 1, j+y)).x;
                    float v_up    = read_imagef(rawImage, (int2)(x+i, j+y - 1)).x;
                    float v_down  = read_imagef(rawImage, (int2)(x+i, j+y + 1)).x;

                    dv += (float2) (fabs(v_left - v_right), fabs(v_up - v_down));
                }
            }
        }

        // we're doing edge directed bilinear interpolation on the green channel,
        // which is a low pass operation (averaging), so we add some signal from the
        // high frequencies of the observed color channel

        // Estimate the whiteness of the pixel value and use that to weight the amount of HF correction
        float cMax = fmax(c_xy, fmax(g_ave, c2_ave));
        float cMin = fmin(c_xy, fmin(g_ave, c2_ave));
        float whiteness = smoothstep(0.25, 0.35, cMin/cMax);

        float sample_h = (g_left + g_right) / 2 + whiteness * (c_xy - (c_left + c_right) / 2) / 4;
        float sample_v = (g_up + g_down) / 2 + whiteness * (c_xy - (c_up + c_down) / 2) / 4;
        float sample_flat = g_ave + whiteness * (c_xy - (c_left + c_right + c_up + c_down) / 4) / 4;

        // TODO: replace eps with some multiple of sigma from the noise model
        float eps = 0.01;
        float flatness = 1 - smoothstep(eps / 2.0, eps, fabs(dv.x - dv.y));
        float sample = mix(dv.x > dv.y ? sample_v : sample_h, sample_flat, flatness);

        write_imagef(greenImage, imageCoordinates, clamp(sample, 0.0, 1.0));
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

            float2 dc = (float2) (fabs(c2_top_left - c2_bottom_right), fabs(c2_top_right - c2_bottom_left));

            float alpha = length(dc) > 0.01 ? atan2(dc.y, dc.x) / M_PI_2_F : 0.5;
            float c2 = green - mix((c2_top_right + c2_bottom_left) / 2,
                                   (c2_top_left + c2_bottom_right) / 2, alpha);

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

    write_imagef(rgbImage, imageCoordinates, (float4)(clamp((float3)(red, green, blue), 0.0, 1.0), 0));
}

/// ---- Image Denoising ----

float3 denoiseRGB(float3 inputValue, image2d_t inputImage, int2 imageCoordinates) {
    const int filterSize = 5;
    const float sigmaR = 0.0025;

    float3 filtered_pixel = 0;
    float3 kernel_norm = 0;
    for (int y = -filterSize / 2; y <= filterSize / 2; y++) {
        for (int x = -filterSize / 2; x <= filterSize / 2; x++) {
            float3 inputSample = read_imagef(inputImage, imageCoordinates + (int2)(x, y)).xyz;

            float3 inputDiff = (inputSample - inputValue);
            float sampleWeight = exp(-0.3 * dot(inputDiff, inputDiff) / (2 * sigmaR * sigmaR));

            filtered_pixel += sampleWeight * inputSample;
            kernel_norm += sampleWeight;
        }
    }
    return filtered_pixel / kernel_norm;
}

float3 rgbToYCbCr(float3 inputValue) {
    const float3 matrix[3] = {
        {  0.2126,  0.7152,  0.0722 },
        { -0.1146, -0.3854,  0.5    },
        {  0.5,    -0.4542, -0.0458 }
    };

    return (float3) (dot(matrix[0], inputValue), dot(matrix[1], inputValue), dot(matrix[2], inputValue));
}

float3 yCbCrtoRGB(float3 inputValue) {
    const float3 matrix[3] = {
        { 1.0,  0.0,      1.5748 },
        { 1.0, -0.1873,  -0.4681 },
        { 1.0,  1.8556,   0.0    }
    };

    return (float3) (dot(matrix[0], inputValue), dot(matrix[1], inputValue), dot(matrix[2], inputValue));
}

kernel void rgbToYCbCrImage(read_only image2d_t inputImage, write_only image2d_t outputImage) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    float3 outputPixel = rgbToYCbCr(read_imagef(inputImage, imageCoordinates).xyz);

    write_imagef(outputImage, imageCoordinates, (float4) (outputPixel, 0.0));
}

kernel void yCbCrtoRGBImage(read_only image2d_t inputImage, write_only image2d_t outputImage) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    float3 outputPixel = yCbCrtoRGB(read_imagef(inputImage, imageCoordinates).xyz);

    write_imagef(outputImage, imageCoordinates, (float4) (outputPixel, 0.0));
}

float3 denoiseLumaChroma(float3 inputValue, image2d_t inputImage, int2 imageCoordinates) {
    float3 filtered_pixel = 0;
    float3 kernel_norm = 0;
    const int filterSize = 7;
    const float lumaSigmaR = 0.0005;
    const float chromaSigmaR = 0.005;

    const float3 inputYCC = read_imagef(inputImage, imageCoordinates).xyz;

    for (int y = -filterSize / 2; y <= filterSize / 2; y++) {
        for (int x = -filterSize / 2; x <= filterSize / 2; x++) {
            float3 inputSampleYCC = read_imagef(inputImage, imageCoordinates + (int2)(x, y)).xyz;

            float3 inputDiff = inputSampleYCC - inputYCC;
            float lumaWeight = exp(-0.3 * dot(inputDiff.x, inputDiff.x) / (2 * lumaSigmaR * lumaSigmaR));
            float chromaWeight = exp(-0.3 * dot(inputDiff.yz, inputDiff.yz) / (2 * chromaSigmaR * chromaSigmaR));

            float3 sampleWeight = (float3) (lumaWeight, chromaWeight, chromaWeight);

            filtered_pixel += sampleWeight * inputSampleYCC;
            kernel_norm += sampleWeight;
        }
    }
    return filtered_pixel / kernel_norm;
}

kernel void denoiseImage(read_only image2d_t inputImage, write_only image2d_t denoisedImage) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    float3 inputPixel = read_imagef(inputImage, imageCoordinates).xyz;
    float3 denoisedPixel = denoiseLumaChroma(inputPixel, inputImage, imageCoordinates);

    write_imagef(denoisedImage, imageCoordinates, (float4) (denoisedPixel, 0.0));
}

/// ---- Tone Curve ----

float3 sigmoid(float3 x, float s) {
    return 0.5 * (tanh(s * x - 0.3 * s) + 1);
}

// This tone curve is designed to mostly match the default curve from DNG files
// TODO: it would be nice to have separate control on highlights and shhadows contrast

float3 toneCurve(float3 x) {
    float s = 3.5;
    return (sigmoid(native_powr(0.95 * x, 0.5), s) - sigmoid(0, s)) / (sigmoid(1, s) - sigmoid(0, s));
}

float3 saturationBoost(float3 value, float saturation) {
    // Saturation boost with highlight protection
    const float luma = 0.2126 * value.x + 0.7152 * value.y + 0.0722 * value.z; // BT.709-2 (sRGB) luma primaries
    const float3 clipping = smoothstep(0.75, 2.0, value);
    return mix(luma, value, mix(saturation, 1.0, clipping));
}

float3 contrastBoost(float3 value, float contrast) {
    const float gray = 0.2;
    const float3 clipping = smoothstep(0.9, 2.0, value);
    return mix(gray, value, mix(contrast, 1.0, clipping));
}

float3 gaussianBlur(image2d_t inputImage, int2 imageCoordinates) {
    // Compute a blurred version of the image
    float3 blurred_pixel = 0;
    float3 kernel_norm = 0;
    const int filterSize = 7;
    const float sigmaS = (float) filterSize / 3.0;
    for (int y = -filterSize / 2; y <= filterSize / 2; y++) {
        for (int x = -filterSize / 2; x <= filterSize / 2; x++) {
            float kernelWeight = native_exp(-((float)(x * x + y * y) / (2 * sigmaS * sigmaS)));
            blurred_pixel += kernelWeight * read_imagef(inputImage, imageCoordinates + (int2)(x, y)).xyz;
            kernel_norm += kernelWeight;
        }
    }
    return blurred_pixel / kernel_norm;
}

float gaussianBlurLuma(image2d_t inputImage, int2 imageCoordinates) {
    // Compute a blurred version of the image
    float blurred_pixel = 0;
    float kernel_norm = 0;
    const int filterSize = 7;
    const float sigmaS = (float) filterSize / 3.0;
    for (int y = -filterSize / 2; y <= filterSize / 2; y++) {
        for (int x = -filterSize / 2; x <= filterSize / 2; x++) {
            float kernelWeight = native_exp(-((float)(x * x + y * y) / (2 * sigmaS * sigmaS)));
            blurred_pixel += kernelWeight * read_imagef(inputImage, imageCoordinates + (int2)(x, y)).x;
            kernel_norm += kernelWeight;
        }
    }
    return blurred_pixel / kernel_norm;
}

float3 sharpen(float3 pixel_value, image2d_t inputImage, int2 imageCoordinates) {
    float3 dx = read_imagef(inputImage, imageCoordinates + (int2)(1, 0)).xyz - pixel_value;
    float3 dy = read_imagef(inputImage, imageCoordinates + (int2)(0, 1)).xyz - pixel_value;

    const float amount = 1.5;

    // Smart sharpening
    float3 sharpening = amount * smoothstep(0.0, 0.03, length(dx) + length(dy))       // Gradient magnitude thresholding
                               * (1.0 - smoothstep(0.95, 1.0, pixel_value))          // Highlight ringing protection
                               * (0.6 + 0.4 * smoothstep(0.0, 0.1, pixel_value));   // Shadows ringing protection

    float3 blurred_pixel = gaussianBlur(inputImage, imageCoordinates);

    return mix(blurred_pixel, pixel_value, fmax(sharpening, 1.0));
}

float3 sharpenLuma(image2d_t inputImage, int2 imageCoordinates) {
    float3 inputPixel = read_imagef(inputImage, imageCoordinates).xyz;

    float dx = read_imagef(inputImage, imageCoordinates + (int2)(1, 0)).x - inputPixel.x;
    float dy = read_imagef(inputImage, imageCoordinates + (int2)(0, 1)).x - inputPixel.x;

    const float amount = 1.5;

    // Smart sharpening
    float sharpening = amount * smoothstep(0.0, 0.03, length(dx) + length(dy))       // Gradient magnitude thresholding
                              * (1.0 - smoothstep(0.95, 1.0, inputPixel.x))         // Highlight ringing protection
                              * (0.6 + 0.4 * smoothstep(0.0, 0.1, inputPixel.x));  // Shadows ringing protection

    float blurred_pixel = gaussianBlurLuma(inputImage, imageCoordinates);

    return (float3) (mix(blurred_pixel, inputPixel.x, fmax(sharpening, 1.0)), inputPixel.yz);
}

kernel void sharpenImage(read_only image2d_t inputImage, write_only image2d_t sharpenedImage) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    float3 inputPixel = read_imagef(inputImage, imageCoordinates).xyz;
    float3 sharpenedPixel = sharpen(inputPixel, inputImage, imageCoordinates);

    write_imagef(sharpenedImage, imageCoordinates, (float4) (sharpenedPixel, 0.0));
}

kernel void sharpenLumaImage(read_only image2d_t inputImage, write_only image2d_t sharpenedImage) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    float3 sharpenedPixel = sharpenLuma(inputImage, imageCoordinates);

    write_imagef(sharpenedImage, imageCoordinates, (float4) (sharpenedPixel, 0.0));
}

kernel void convertTosRGB(read_only image2d_t linearImage, write_only image2d_t rgbImage,
                          constant float3 transform[3]) {
    const int2 imageCoordinates = (int2) (get_global_id(0), get_global_id(1));

    float3 linear = read_imagef(linearImage, imageCoordinates).xyz;

    linear = sharpen(linear, linearImage, imageCoordinates);

    linear = contrastBoost(linear, 1.2);

    float3 rgb = (float3) (dot(transform[0], linear), dot(transform[1], linear), dot(transform[2], linear));

    write_imagef(rgbImage, imageCoordinates, (float4) (toneCurve(clamp(rgb, 0.0, 1.0)), 0.0));
}
