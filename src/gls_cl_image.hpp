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

#ifndef CL_IMAGE_H
#define CL_IMAGE_H

#include "gls_cl.hpp"
#include "gls_image.hpp"

namespace gls {

template <typename T>
class cl_image : public basic_image<T> {
   public:
    typedef std::unique_ptr<cl_image<T>> unique_ptr;

    cl_image(int _width, int _height) : basic_image<T>(_width, _height) {}

    static inline cl::ImageFormat ImageFormat();
};

template <typename T>
class cl_image_2d : public cl_image<T> {
   protected:
    struct payload {
        const cl::Image2D image;
    };

    const std::unique_ptr<payload> _payload;

    cl_image_2d(cl::Context context, int _width, int _height, std::unique_ptr<payload> payload)
        : cl_image<T>(_width, _height), _payload(std::move(payload)) {}

   public:
    typedef std::unique_ptr<cl_image_2d<T>> unique_ptr;

    cl_image_2d(cl::Context context, int _width, int _height)
        : cl_image<T>(_width, _height), _payload(buildPayload(context, _width, _height)) {}

    cl_image_2d(cl::Context context, const gls::image<T>& other)
        : cl_image<T>(other.width, other.height),
          _payload(buildPayload(context, other.width, other.height, other.pixels().data())) {}

    static inline std::unique_ptr<payload> buildPayload(cl::Context context, int width, int height, T* data = nullptr) {
        cl_mem_flags mem_flags = CL_MEM_READ_WRITE | (data ? CL_MEM_COPY_HOST_PTR : 0);
        return std::make_unique<payload>(payload{cl::Image2D(context, mem_flags, cl_image<T>::ImageFormat(), width,
                                                             height, 0, data)});
    }

    static inline unique_ptr fromImage(cl::Context context, const gls::image<T>& other) {
        return std::make_unique<cl_image_2d<T>>(context, other);
    }

    inline typename gls::image<T>::unique_ptr toImage() const {
        auto image = std::make_unique<gls::image<T>>(gls::image<T>::width, gls::image<T>::height);
        copyPixelsTo(image.get());
        return image;
    }

    void copyPixelsFrom(const image<T>& other) const {
        assert(other.width == image<T>::width && other.height == image<T>::height);
        cl::enqueueWriteImage(_payload->image, true, {0, 0, 0}, {(size_t)image<T>::width, (size_t)image<T>::height, 1},
                              image<T>::pixel_size * image<T>::width, 0, other.pixels().data());
    }

    void copyPixelsTo(image<T>* other) const {
        assert(other->width == image<T>::width && other->height == image<T>::height);
        cl::enqueueReadImage(_payload->image, CL_TRUE, {0, 0, 0},
                             {(size_t)image<T>::width, (size_t)image<T>::height, 1},
                             image<T>::pixel_size * image<T>::width, 0, other->pixels().data());
    }

    virtual image<T> mapImage() {
        size_t row_pitch;
        size_t slice_pitch;
        cl::CommandQueue queue = cl::CommandQueue::getDefault();
        T* image_data =
            (T*)queue.enqueueMapImage(_payload->image, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, {0, 0, 0},
                                      {(size_t)image<T>::width, (size_t)image<T>::height, 1}, &row_pitch, &slice_pitch);
        assert(image_data != nullptr);

        size_t stride = row_pitch / image<T>::pixel_size;
        size_t data_size = stride * image<T>::height;
        return gls::image(image<T>::width, image<T>::height, (int) stride, std::span<T>(image_data, data_size));
    }

    void unmapImage(const image<T>& mappedImage) { cl::enqueueUnmapMemObject(_payload->image, (void*)mappedImage[0]); }

    cl::Image2D getImage2D() const { return _payload->image; }
};

template <typename T>
class cl_image_buffer_2d : public cl_image_2d<T> {
   private:
    struct payload : public cl_image_2d<T>::payload {
        const cl::Buffer buffer;
    };

   public:
    cl_image_buffer_2d(cl::Context context, int _width, int _height)
        : cl_image_2d<T>(context, _width, _height, buildPayload(context, _width, _height)) {}

    static inline std::unique_ptr<payload> buildPayload(cl::Context context, int width, int height) {
        const int pitch_alignment = cl::Device::getDefault().getInfo<CL_DEVICE_IMAGE_PITCH_ALIGNMENT>();
        const int pitch = pitch_alignment * ((width * sizeof(T) + pitch_alignment - 1) / pitch_alignment) / sizeof(T);

        // TODO: this is too strong of an assumption
        assert(width == pitch);

        auto buffer = cl::Buffer(context, CL_MEM_READ_WRITE, pitch * height * sizeof(T));
        auto image = cl::Image2D(context, cl_image<T>::ImageFormat(), buffer, width, height, pitch * sizeof(T));
        return std::make_unique<payload>(payload{{image}, buffer});
    }

    image<T> mapImage() override {
        size_t pixel_count = image<T>::width * image<T>::height;
        T* image_data =
            cl::enqueueMapBuffer(getBuffer(), true, CL_MAP_READ | CL_MAP_WRITE, 0, image<T>::pixel_size * pixel_count);

        return gls::image(image<T>::width, image<T>::height, std::span<T>(image_data, pixel_count));
    }

    cl::Buffer getBuffer() const { return static_cast<const payload*>(this->_payload.get())->buffer; }
};

template <typename T>
class cl_image_3d : public cl_image<T> {
    const cl::Image3D _image;

   public:
    typedef std::unique_ptr<cl_image_3d<T>> unique_ptr;

    cl_image_3d(cl::Context context, int _width, int _height)
        : cl_image<T>(_width, _height), _image(buildImage(context, _width, _height)) {}

    cl_image_3d(cl::Context context, const gls::image<T>& other) : cl_image_3d(context, other.width, other.height) {
        copyPixelsFrom(other);
    }

    inline typename gls::image<T>::unique_ptr toImage() const {
        auto image = std::make_unique<gls::image<T>>(gls::image<T>::width, gls::image<T>::height);
        copyPixelsTo(image.get());
        return image;
    }

    void copyPixelsFrom(const image<T>& other) const {
        assert(other.width == image<T>::width && other.height == image<T>::height);
        size_t depth = image<T>::height / image<T>::width;
        cl::enqueueWriteImage(_image, true, {0, 0, 0}, {(size_t)image<T>::width, (size_t)image<T>::width, depth},
                              image<T>::pixel_size * image<T>::width,
                              image<T>::width * image<T>::width * image<T>::pixel_size, other.pixels().data());
    }

    void copyPixelsTo(image<T>* other) const {
        assert(other->width == image<T>::width && other->height == image<T>::height);
        size_t depth = image<T>::height / image<T>::width;
        cl::enqueueReadImage(_image, true, {0, 0, 0}, {(size_t)image<T>::width, (size_t)image<T>::width, depth},
                             image<T>::pixel_size * image<T>::width,
                             image<T>::width * image<T>::width * image<T>::pixel_size, other->pixels().data());
    }

    static inline cl::Image3D buildImage(cl::Context context, int width, int height) {
        size_t depth = height / width;
        return cl::Image3D(context, CL_MEM_READ_WRITE, cl_image<T>::ImageFormat(), width, width, depth);
    }

    cl::Image3D getImage3D() const { return _image; }
};

template <typename T>
class cl_image_2d_array : public cl_image<T> {
    const cl::Image2DArray _image;

   public:
    typedef std::unique_ptr<cl_image_2d_array<T>> unique_ptr;

    cl_image_2d_array(cl::Context context, int _width, int _height)
        : cl_image<T>(_width, _height), _image(buildImage(context, _width, _height)) {}

    cl_image_2d_array(cl::Context context, const gls::image<T>& other) : cl_image_2d_array(other.width, other.height) {
        copyPixelsFrom(other);
    }

    inline typename gls::image<T>::unique_ptr toImage() const {
        auto image = std::make_unique<gls::image<T>>(gls::image<T>::width, gls::image<T>::height);
        copyPixelsTo(image.get());
        return image;
    }

    void copyPixelsFrom(const image<T>& other) const {
        assert(other.width == image<T>::width && other.height == image<T>::height);
        size_t depth = image<T>::height / image<T>::width;
        cl::enqueueWriteImage(_image, true, {0, 0, 0}, {(size_t)image<T>::width, (size_t)image<T>::width, depth},
                              image<T>::pixel_size * image<T>::width,
                              image<T>::width * image<T>::width * image<T>::pixel_size, other.pixels().data());
    }

    void copyPixelsTo(image<T>* other) const {
        assert(other->width == image<T>::width && other->height == image<T>::height);
        size_t depth = image<T>::height / image<T>::width;
        cl::enqueueReadImage(_image, true, {0, 0, 0}, {(size_t)image<T>::width, (size_t)image<T>::width, depth},
                             image<T>::pixel_size * image<T>::width,
                             image<T>::width * image<T>::width * image<T>::pixel_size, other->pixels().data());
    }

    static inline cl::Image2DArray buildImage(cl::Context context, int width, int height) {
        size_t depth = height / width;
        return cl::Image2DArray(context, CL_MEM_READ_WRITE, cl_image<T>::ImageFormat(), depth, width, width, width,
                                width * width);
    }

    cl::Image2DArray getImage2DArray() const { return _image; }
};

template <>
inline cl::ImageFormat cl_image<gls::luma_pixel>::ImageFormat() {
    return cl::ImageFormat(CL_R, CL_UNORM_INT8);
}

template <>
inline cl::ImageFormat cl_image<gls::luma_alpha_pixel>::ImageFormat() {
    return cl::ImageFormat(CL_RG, CL_UNORM_INT8);
}

template <>
inline cl::ImageFormat cl_image<gls::rgba_pixel>::ImageFormat() {
    return cl::ImageFormat(CL_RGBA, CL_UNORM_INT8);
}

template <>
inline cl::ImageFormat cl_image<gls::luma_pixel_16>::ImageFormat() {
    return cl::ImageFormat(CL_R, CL_UNORM_INT16);
}

template <>
inline cl::ImageFormat cl_image<gls::luma_alpha_pixel_16>::ImageFormat() {
    return cl::ImageFormat(CL_RG, CL_UNORM_INT16);
}

template <>
inline cl::ImageFormat cl_image<gls::rgba_pixel_16>::ImageFormat() {
    return cl::ImageFormat(CL_RGBA, CL_UNORM_INT16);
}

template <>
inline cl::ImageFormat cl_image<gls::luma_pixel_fp32>::ImageFormat() {
    return cl::ImageFormat(CL_R, CL_FLOAT);
}

template <>
inline cl::ImageFormat cl_image<gls::luma_alpha_pixel_fp32>::ImageFormat() {
    return cl::ImageFormat(CL_RG, CL_FLOAT);
}

template <>
inline cl::ImageFormat cl_image<gls::rgba_pixel_fp32>::ImageFormat() {
    return cl::ImageFormat(CL_RGBA, CL_FLOAT);
}

#ifndef __APPLE__
template <>
inline cl::ImageFormat cl_image<gls::luma_pixel_fp16>::ImageFormat() {
    return cl::ImageFormat(CL_R, CL_HALF_FLOAT);
}

template <>
inline cl::ImageFormat cl_image<gls::luma_alpha_pixel_fp16>::ImageFormat() {
    return cl::ImageFormat(CL_RG, CL_HALF_FLOAT);
}

template <>
inline cl::ImageFormat cl_image<gls::rgba_pixel_fp16>::ImageFormat() {
    return cl::ImageFormat(CL_RGBA, CL_HALF_FLOAT);
}
#endif

}  // namespace gls

#endif /* CL_IMAGE_H */
