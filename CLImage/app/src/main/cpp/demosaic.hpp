//
//  demosaic.hpp
//  CLImageTest
//
//  Created by Fabio Riccardi on 3/18/22.
//

#ifndef demosaic_hpp
#define demosaic_hpp

#include "gls_image.hpp"

void interpolateGreen(const gls::image<gls::luma_pixel_16>& rawImage, gls::image<gls::rgb_pixel_16>* rgbImage, int gx, int gy, int ry );

void interpolateRedBlue(gls::image<gls::rgb_pixel_16>* image, int rx0, int ry0, int bx0, int by0 );

#endif /* demosaic_hpp */
