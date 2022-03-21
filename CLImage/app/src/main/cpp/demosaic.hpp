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

enum { red = 0, green = 1, blue = 2 };

inline static std::array<gls::point, 3> bayerOffsets(BayerPattern bayerPattern) {
    switch (bayerPattern) {
        case grbg:
            return { gls::point{1, 0}, gls::point{0, 0}, gls::point{0, 1} };
        case gbrg:
            return { gls::point{0, 1}, gls::point{0, 0}, gls::point{1, 0} };
        case rggb:
            return { gls::point{0, 0}, gls::point{1, 0}, gls::point{1, 1} };
        case bggr:
            return { gls::point{1, 1}, gls::point{1, 0}, gls::point{0, 0} };
    }
}

void interpolateGreen(const gls::image<gls::luma_pixel_16>& rawImage, gls::image<gls::rgb_pixel_16>* rgbImage, BayerPattern bayerPattern);

void interpolateRedBlue(gls::image<gls::rgb_pixel_16>* image, BayerPattern bayerPattern);

#endif /* demosaic_hpp */
