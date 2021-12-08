//
// Created by Fabio Riccardi on 12/8/21.
//

#ifndef CLIMAGE_CL_PIPELINE_H
#define CLIMAGE_CL_PIPELINE_H

#include "climage/gls_cl_image.hpp"

// Simple image processing with opencl.hpp, using cl_image to pass data to and from the GPU
int blur(const gls::cl_image_2d<gls::rgba_pixel>& input, gls::cl_image_2d<gls::rgba_pixel>* output);

#endif //CLIMAGE_CL_PIPELINE_H
