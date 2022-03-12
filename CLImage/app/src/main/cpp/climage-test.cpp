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

#include <jni.h>
#include <string>

#include "gls_android_support.h"

#include "gls_cl.hpp"
#include "gls_logging.h"
#include "gls_cl_image.hpp"

#include "cl_pipeline.h"

static const char* TAG = "CLImage Test";

extern "C" JNIEXPORT jint JNICALL
Java_com_glassimaging_climage_TestCLImage_testCLImage(
        JNIEnv* env,
        jobject /* this */,
        jobject assetManager,
        jobject inputBitmap,
        jobject outputBitmap) {
    // Bind input and output Bitmaps
    gls::AndroidBitmap input(env, inputBitmap);
    gls::AndroidBitmap output(env, outputBitmap);

    // Initialize the OpenCL environment and get the context
    gls::OpenCLContext glsContext("");
    auto clContext = glsContext.clContext();

    // Load OpenCL Shaders stored in App's assets
    gls::loadOpenCLShaders(env, assetManager, glsContext.getShadersMap());
    gls::loadOpenCLBytecode(env, assetManager, glsContext.getBytecodeMap());

    try {
        // Input Image for OpenCL processing
        gls::cl_image_2d<gls::rgba_pixel>::unique_ptr clInputImage;
        {
            // Input image wrapped around the Android Bitmap
            gls::image<gls::rgba_pixel> inputImage(input.info().width, input.info().height,
                                                    input.lockPixels<gls::rgba_pixel>());
            // Move data from the Input Bitmap to the OpenCL Image
            clInputImage = std::make_unique<gls::cl_image_2d<gls::rgba_pixel>>(clContext, inputImage);

            // Release Bitmap Object
            input.unLockPixels();
        }

        // Output Image from OpenCL processing
        gls::cl_image_2d<gls::rgba_pixel> clOutputImage(clContext, clInputImage->width, clInputImage->height);

        // Execute OpenCL Blur algorithm
        if (blur(&glsContext, *clInputImage, &clOutputImage) == 0) {
            LOG_INFO(TAG) << "All done with Blur" << std::endl;
        } else {
            LOG_ERROR(TAG) << "Something wrong with the Blur." << std::endl;
        }

        // Output image wrapped around the Android Bitmap
        gls::image<gls::rgba_pixel> outputImage(clOutputImage.width, clOutputImage.height, output.lockPixels<gls::rgba_pixel>());

        // Wait for GPU to complete execution and move data to the output Bitmap
        clOutputImage.copyPixelsTo(&outputImage);

        // Release Bitmap object
        output.unLockPixels();

        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what())
                       << " - " << gls::clStatusToString(err.err()) << std::endl;
        return -1;
    }
}
