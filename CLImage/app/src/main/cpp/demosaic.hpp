//
//  demosaic.hpp
//  CLImageTest
//
//  Created by Fabio Riccardi on 3/18/22.
//

#ifndef demosaic_hpp
#define demosaic_hpp

void interpolateGreen(unsigned short *srcData, unsigned short *destData,
                      int width, int height, int srcLineStride, int destLineStride,
                      int srcOffset, int rOffset, int gOffset, int bOffset,
                      int gx, int gy, int ry);

void interpolateRedBlue(unsigned short *data, int width, int height, int lineStride,
                        int rOffset, int gOffset, int bOffset,
                        int rx0, int ry0, int bx0, int by0);

#endif /* demosaic_hpp */
