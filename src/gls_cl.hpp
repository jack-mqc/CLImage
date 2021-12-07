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

#ifndef GLS_CL_HPP
#define GLS_CL_HPP

#define CL_HPP_ENABLE_EXCEPTIONS

#ifdef __APPLE__
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_USE_CL_IMAGE2D_FROM_BUFFER_KHR true

#include <OpenCL/cl_ext.h>

#include "CL/opencl.hpp"
#elif __ANDROID__
#include <map>

#define CL_TARGET_OPENCL_VERSION 200
#define CL_HPP_TARGET_OPENCL_VERSION 200

// Include cl_icd_wrapper.h before <CL/*>
#include "gls_icd_wrapper.h"

#include <CL/cl_ext.h>

#include <CL/opencl.hpp>
#endif

namespace gls {
    std::string clStatusToString(cl_int status);
}  // namespace gls
#endif /* GLS_CL_HPP */
