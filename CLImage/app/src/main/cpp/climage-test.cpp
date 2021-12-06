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

#define CL_TARGET_OPENCL_VERSION 200
#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS true

// Include cl_icd_wrapper.h before <CL/*>
#include "cl_icd_wrapper.h"

#include "CL/opencl.hpp"

#include "logging.h"
#include "cl_error.h"

static const char* TAG = "CLImage Test";

cl::Context getContext() {
    static bool initialized = false;

    if (!initialized) {
        CL_WRAPPER_NS::bindOpenCLLibrary();

        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        cl::Platform platform;
        for (auto& p : platforms) {
            std::string version = p.getInfo<CL_PLATFORM_VERSION>();
            if (version.find("OpenCL 2.") != std::string::npos) {
                platform = p;
            }
        }
        if (platform() == 0) {
            throw cl::Error(-1, "No OpenCL 2.0 platform found.");
        }

        cl::Platform defaultPlatform = cl::Platform::setDefault(platform);
        if (defaultPlatform != platform) {
            throw cl::Error(-1, "Error setting default platform.");
        }

        cl_context_properties properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(platform)(), 0};
        cl::Context context(CL_DEVICE_TYPE_ALL, properties);

        cl::Device d = cl::Device::getDefault();
        LOG_INFO(TAG) << "- Device: " << d.getInfo<CL_DEVICE_NAME>() << std::endl;
        LOG_INFO(TAG) << "- Device Version: " << d.getInfo<CL_DEVICE_VERSION>() << std::endl;
        LOG_INFO(TAG) << "- Driver Version: " << d.getInfo<CL_DRIVER_VERSION>() << std::endl;
        LOG_INFO(TAG) << "- OpenCL C Version: " << d.getInfo<CL_DEVICE_OPENCL_C_VERSION>() << std::endl;
        LOG_INFO(TAG) << "- Compute Units: " << d.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
        LOG_INFO(TAG) << "- CL_DEVICE_EXTENSIONS: " << d.getInfo<CL_DEVICE_EXTENSIONS>() << std::endl;

        LOG_INFO(TAG) << "- CL_DEVICE_MAX_WORK_GROUP_SIZE: " << d.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << std::endl;

        cl::Context::setDefault(context);

        initialized = true;

        return context;
    } else
        return cl::Context::getDefault();
}

extern "C" JNIEXPORT jstring JNICALL
Java_com_glassimaging_climage_TestCLImage_testCLImage(
        JNIEnv* env,
        jobject /* this */) {

    // Load and bind OpenCL library
    try {
        cl::Context context = getContext();
        cl::Context::setDefault(context);
        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << getOpenCLError(err.err())
                       << std::endl;
    }

    std::string hello = "Hello from C++";
    return env->NewStringUTF(hello.c_str());
}
