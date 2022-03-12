/*******************************************************************************
 * Copyright (c) 2021-2022 Glass Imaging Inc.
 * Author: Fabio Riccardi <fabio@glass-imaging.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

#include <string>

#include "gls_cl.hpp"
#include "gls_logging.h"
#include "gls_cl_image.hpp"

#include "cl_pipeline.h"

static const char* TAG = "CLImage Test";

int main(int argc, const char* argv[]) {
    printf("Hello CLImage!\n");

    if (argc > 1) {
        // Initialize the OpenCL environment and get the context
        gls::OpenCLContext glsContext("");
        auto clContext = glsContext.clContext();

        // Read the input file into an image object
        auto inputImage = gls::image<gls::rgba_pixel>::read_png_file(argv[1]);

        LOG_INFO(TAG) << "inputImage size: " << inputImage->width << " x " << inputImage->height << std::endl;

        // Load image data in OpenCL image texture
        gls::cl_image_2d<gls::rgba_pixel> clInputImage(clContext, *inputImage);

        // Output Image from OpenCL processing
        gls::cl_image_2d<gls::rgba_pixel> clOutputImage(clContext, clInputImage.width, clInputImage.height);

        // Execute OpenCL Blur algorithm
        if (blur(&glsContext, clInputImage, &clOutputImage) == 0) {
            LOG_INFO(TAG) << "All done with Blur" << std::endl;
        } else {
            LOG_ERROR(TAG) << "Something wrong with the Blur." << std::endl;
        }

        // Use OpenCL's memory mapping for zero-copy image Output
        auto outputImage = clOutputImage.mapImage();
        outputImage.write_png_file("output.png");
        clOutputImage.unmapImage(outputImage);
    }
}
