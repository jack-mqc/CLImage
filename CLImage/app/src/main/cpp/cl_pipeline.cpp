//
// Created by Fabio Riccardi on 12/8/21.
//

#include "cl_pipeline.h"

#include "climage/gls_logging.h"

static const char* TAG = "CLImage Pipeline";

// Simple image processing with opencl.hpp, using cl_image to pass data to and from the GPU
int blur(const gls::cl_image_2d<gls::rgba_pixel>& input, gls::cl_image_2d<gls::rgba_pixel>* output) {
    try {
        // Load the shader source
        const auto blurProgram = gls::loadOpenCLProgram("blur");

        // Bind the kernel's parameters
        auto blurKernel = cl::KernelFunctor<cl::Image2D,  // input
                cl::Image2D   // output
        >(*blurProgram, "blur");

        // Schedule the kernel on the GPU
        blurKernel(gls::buildEnqueueArgs(output->width, output->height), input.getImage2D(), output->getImage2D());
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}
