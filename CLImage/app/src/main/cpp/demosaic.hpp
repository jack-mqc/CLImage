//
//  demosaic.hpp
//  CLImageTest
//
//  Created by Fabio Riccardi on 3/18/22.
//

#ifndef demosaic_hpp
#define demosaic_hpp

#include "gls_image.hpp"

enum BayerPattern {
    grbg = 0,
    gbrg = 1,
    rggb = 2,
    bggr = 3
};

void interpolateGreen(const gls::image<gls::luma_pixel_16>& rawImage, gls::image<gls::rgb_pixel_16>* rgbImage, BayerPattern bayerPattern);

void interpolateRedBlue(gls::image<gls::rgb_pixel_16>* image, BayerPattern bayerPattern);

gls::image<gls::rgb_pixel_16>::unique_ptr demosaicImage(const gls::image<gls::luma_pixel_16>& rawImage);

#endif /* demosaic_hpp */
