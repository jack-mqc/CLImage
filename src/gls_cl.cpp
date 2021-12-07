//
// Created by Fabio Riccardi on 12/7/21.
//

#include "gls_cl.hpp"
#include "gls_logging.h"

namespace gls {

static const char* TAG = "CLImage";

#ifdef __APPLE__

cl::Context getContext() {
    cl::Context context = cl::Context::getDefault();

    static bool initialized = false;
    if (!initialized) {
        std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

        // Macs have multiple GPUs, select the one with most compute units
        int max_compute_units = 0;
        cl::Device best_device;
        for (const auto& d : devices) {
            int device_compute_units = d.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
            if (device_compute_units > max_compute_units) {
                max_compute_units = device_compute_units;
                best_device = d;
            }
        }
        cl::Device::setDefault(best_device);

        cl::Device d = cl::Device::getDefault();
        LOG_INFO(TAG) << "OpenCL Default Device: " << d.getInfo<CL_DEVICE_NAME>() << std::endl;
        LOG_INFO(TAG) << "- Device Version: " << d.getInfo<CL_DEVICE_VERSION>() << std::endl;
        LOG_INFO(TAG) << "- Driver Version: " << d.getInfo<CL_DRIVER_VERSION>() << std::endl;
        LOG_INFO(TAG) << "- OpenCL C Version: " << d.getInfo<CL_DEVICE_OPENCL_C_VERSION>() << std::endl;
        LOG_INFO(TAG) << "- Compute Units: " << d.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
        LOG_INFO(TAG) << "- CL_DEVICE_EXTENSIONS: " << d.getInfo<CL_DEVICE_EXTENSIONS>() << std::endl;

        LOG_INFO(TAG) << "- CL_DEVICE_MAX_WORK_GROUP_SIZE: " << d.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << std::endl;

        initialized = true;
    }
    return context;
}

#elif __ANDROID__

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
        if (platform() == nullptr) {
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
#endif

}  // namespace gls
