/*******************************************************************************
 * Copyright (c) 2021 Glass Imaging Inc.
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

#include "android_support.h"

#include "gls_cl.hpp"
#include "gls_logging.h"
#include "gls_cl_image.hpp"

static const char* TAG = "CLImage Test";

int blur(const gls::cl_image_2d<gls::rgba_pixel>& input, gls::cl_image_2d<gls::rgba_pixel>* output) {
    try {
        const auto blurProgram = gls::loadOpenCLProgram("blur");

        auto blurKernel = cl::KernelFunctor<cl::Image2D,  // input
                                            cl::Image2D   // output
                                            >(*blurProgram, "blur");

        blurKernel(gls::buildEnqueueArgs(output->width, output->height), input.getImage2D(), output->getImage2D());
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}

int initializeOpenCL() {
    try {
        cl::Context context = gls::getContext();
        cl::Context::setDefault(context);
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}

extern "C" JNIEXPORT jint JNICALL
Java_com_glassimaging_climage_TestCLImage_testCLImage(
        JNIEnv* env,
        jobject /* this */,
        jobject assetManager,
        jobject inputBitmap,
        jobject outputBitmap) {

    // Load and bind OpenCL library
    initializeOpenCL();

    // Load OpenCL Shader resources
    loadOpenCLShaders(env, assetManager, gls::getShadersMap());

    // Bind input and output Bitmaps
    AndroidBitmap input(env, inputBitmap);
    AndroidBitmap output(env, outputBitmap);

    // Load and bind OpenCL library
    try {
        const auto context = gls::getContext();

        // Input Image
        gls::cl_image_2d<gls::rgba_pixel>::unique_ptr inputImage;
        {
            gls::image<gls::rgba_pixel> inputCPU(input.info().width, input.info().height,
                                                 input.lockPixels<gls::rgba_pixel>());
            inputImage = std::make_unique<gls::cl_image_2d<gls::rgba_pixel>>(context, inputCPU);
            input.unLockPixels();
        }

        gls::cl_image_2d<gls::rgba_pixel> outputImage(context, inputImage->width, inputImage->height);

        if (blur(*inputImage, &outputImage) == 0) {
            LOG_INFO(TAG) << "Blur went fine" << std::endl;
        }

        gls::image<gls::rgba_pixel> outputImageCPU(outputImage.width, outputImage.height, output.lockPixels<gls::rgba_pixel>());

        for (auto& p : outputImageCPU.pixels()) {
            p = gls::rgba_pixel(0, 255, 0, 255);
        }

        outputImage.copyPixelsTo(&outputImageCPU);
        output.unLockPixels();

        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what())
                       << " - " << gls::clStatusToString(err.err()) << std::endl;
        return -1;
    }
}
