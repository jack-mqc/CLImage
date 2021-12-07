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

#ifndef gls_image_h
#define gls_image_h

#include <sys/types.h>

#include <functional>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gls_image_jpeg.h"
#include "gls_image_png.h"

namespace gls {

struct point {
    int x;
    int y;

    point(int _x, int _y) : x(_x), y(_y) {}
};

template <typename T>
struct basic_luma_pixel {
    constexpr static int channels = 1;
    constexpr static int bit_depth = 8 * sizeof(T);

    union {
        struct {
            T luma;
        };
        T v[channels];
    };

    basic_luma_pixel() {}
    basic_luma_pixel(T _luma) : luma(_luma) {}

    T& operator[](int c) { return v[c]; }
    const T& operator[](int c) const { return v[c]; }
};

template <typename T>
struct basic_luma_alpha_pixel {
    constexpr static int channels = 2;
    constexpr static int bit_depth = 8 * sizeof(T);

    union {
        struct {
            T luma;
            T alpha;
        };
        T v[channels];
    };

    basic_luma_alpha_pixel() {}
    basic_luma_alpha_pixel(T _luma, T _alpha) : luma(_luma), alpha(_alpha) {}

    T& operator[](int c) { return v[c]; }
    const T& operator[](int c) const { return v[c]; }
};

template <typename T>
struct basic_rgb_pixel {
    constexpr static int channels = 3;
    constexpr static int bit_depth = 8 * sizeof(T);

    union {
        struct {
            T red;
            T green;
            T blue;
        };
        T v[channels];
    };

    basic_rgb_pixel() {}
    basic_rgb_pixel(T _red, T _green, T _blue) : red(_red), green(_green), blue(_blue) {}

    T& operator[](int c) { return v[c]; }
    const T& operator[](int c) const { return v[c]; }
};

template <typename T>
struct basic_rgba_pixel {
    constexpr static int channels = 4;
    constexpr static int bit_depth = 8 * sizeof(T);

    union {
        struct {
            T red;
            T green;
            T blue;
            T alpha;
        };
        T v[channels];
    };

    basic_rgba_pixel() {}
    basic_rgba_pixel(T _red, T _green, T _blue, T _alpha) : red(_red), green(_green), blue(_blue), alpha(_alpha) {}

    T& operator[](int c) { return v[c]; }
    const T& operator[](int c) const { return v[c]; }
};

typedef basic_luma_pixel<uint8_t> luma_pixel;
typedef basic_luma_alpha_pixel<uint8_t> luma_alpha_pixel;
typedef basic_rgb_pixel<uint8_t> rgb_pixel;
typedef basic_rgba_pixel<uint8_t> rgba_pixel;

typedef basic_luma_pixel<uint16_t> luma_pixel_16;
typedef basic_luma_alpha_pixel<uint16_t> luma_alpha_pixel_16;
typedef basic_rgb_pixel<uint16_t> rgb_pixel_16;
typedef basic_rgba_pixel<uint16_t> rgba_pixel_16;

typedef basic_luma_pixel<float> luma_pixel_fp32;
typedef basic_luma_alpha_pixel<float> luma_alpha_pixel_fp32;
typedef basic_rgb_pixel<float> rgb_pixel_fp32;
typedef basic_rgba_pixel<float> rgba_pixel_fp32;

#ifndef __APPLE__
typedef __fp16 float16_t;
typedef basic_luma_pixel<float16_t> luma_pixel_fp16;
typedef basic_luma_alpha_pixel<float16_t> luma_alpha_pixel_fp16;
typedef basic_rgb_pixel<float16_t> rgb_pixel_fp16;
typedef basic_rgba_pixel<float16_t> rgba_pixel_fp16;
#endif

#ifdef __APPLE__
typedef float float_type;
#else
typedef float16_t float_type;
#endif
typedef basic_luma_pixel<float_type> luma_pixel_float;
typedef basic_luma_alpha_pixel<float_type> luma_alpha_pixel_float;
typedef basic_rgb_pixel<float_type> rgb_pixel_float;
typedef basic_rgba_pixel<float_type> rgba_pixel_float;

template <typename T>
class basic_image {
   public:
    const int width;
    const int height;

    static const constexpr int pixel_bit_depth = T::bit_depth;
    static const constexpr int pixel_channels = T::channels;
    static const constexpr int pixel_size = sizeof(T);

    typedef std::unique_ptr<basic_image<T>> unique_ptr;

    basic_image(int _width, int _height) : width(_width), height(_height) {}
};

template <typename T>
class image : public basic_image<T> {
   public:
    typedef std::unique_ptr<image<T>> unique_ptr;

   protected:
    const std::unique_ptr<std::vector<T>> _data_store = nullptr;
    const std::span<T> _data;

    // Data is owned by caller, the image is only a wrapper around it
    image(int _width, int _height, std::unique_ptr<std::vector<T>> data_store)
        : basic_image<T>(_width, _height),
          _data_store(std::move(data_store)),
          _data(_data_store->data(), _data_store->size()) {
        assert(_width * _height <= _data.size());
    }

   public:
    // Data is owned by the image and retained by _data_store
    image(int _width, int _height)
        : basic_image<T>(_width, _height),
          _data_store(std::make_unique<std::vector<T>>(_width * _height)),
          _data(_data_store->data(), _data_store->size()) {}

    // Data is owned by caller, the image is only a wrapper around it
    image(int _width, int _height, std::span<T> data) : basic_image<T>(_width, _height), _data(data) {
        assert(_width * _height <= data.size());
    }

    // row access
    T* operator[](int row) { return &_data[basic_image<T>::width * row]; }

    const T* operator[](int row) const { return &_data[basic_image<T>::width * row]; }

    const std::span<T> pixels() const { return _data; }

    const size_t size_in_bytes() { return _data.size() * basic_image<T>::pixel_size; }

    // image factory from PNG file
    static unique_ptr read_png_file(const std::string& filename) {
        unique_ptr image = nullptr;

        auto image_allocator = [&image](int width, int height, std::vector<uint8_t*>* row_pointers) -> bool {
            if ((image = std::make_unique<gls::image<T>>(width, height)) == nullptr) {
                return false;
            }
            for (int i = 0; i < height; ++i) {
                (*row_pointers)[i] = (uint8_t*)(*image)[i];
            }
            return true;
        };

        gls::read_png_file(filename, T::channels, T::bit_depth, image_allocator);

        return image;
    }

    // Write image to PNG file
    int write_png_file(const std::string& filename) const {
        auto row_pointer = [this](int row) -> uint8_t* { return (uint8_t*)(*this)[row]; };
        return gls::write_png_file(filename, basic_image<T>::width, basic_image<T>::height, T::channels, T::bit_depth, row_pointer);
    }

    // image factory from JPEG file
    static unique_ptr read_jpeg_file(const std::string& filename) {
        unique_ptr image = nullptr;

        auto image_allocator = [&image](int width, int height) -> std::span<uint8_t> {
            if ((image = std::make_unique<gls::image<T>>(width, height)) == nullptr) {
                return std::span<uint8_t>();
            }
            return std::span<uint8_t>((uint8_t*)(*image)[0], sizeof(T) * width * height);
        };

        gls::read_jpeg_file(filename, T::channels, T::bit_depth, image_allocator);

        return image;
    }

    // Write image to JPEG file
    int write_jpeg_file(const std::string& filename, int quality) const {
        auto image_data = [this]() -> std::span<uint8_t> {
            return std::span<uint8_t>((uint8_t*)this->_data.data(), sizeof(T) * this->_data.size());
        };
        return gls::write_jpeg_file(filename, basic_image<T>::width, basic_image<T>::height, T::channels, T::bit_depth, image_data, quality);
    }
};

template <typename T>
class image_3d : public image<T> {
   public:
    const int depth;

    typedef std::unique_ptr<image_3d<T>> unique_ptr;

    // Data is owned by the image and retained by _data_store
    image_3d(int _width, int _height, int _depth)
        : image<T>(_width, _height, std::make_unique<std::vector<T>>(_width * _height * _depth)), depth(_depth) {
        assert(image<T>::width * image<T>::height * depth == image<T>::data.size());
    }

    // Data is owned by caller, the image is only a wrapper around it
    image_3d(int _width, int _height, int _depth, std::span<T> data)
        : image<T>::width(_width), image<T>::height(_height), image<T>::_data(data), depth(_depth) {
        assert(image<T>::width * image<T>::height * depth == image<T>::data.size());
    }

    // plane access
    image<T> operator[](int plane) {
        return image<T>(
            image<T>::width, image<T>::height,
            std::span(image<T>::_data[image<T>::width * image<T>::height * plane], image<T>::width * image<T>::height));
    }

    const image<T> operator[](int plane) const {
        return image<T>(
            image<T>::width, image<T>::height,
            std::span(image<T>::_data[image<T>::width * image<T>::height * plane], image<T>::width * image<T>::height));
    }

    // image factory from PNG file
    static unique_ptr read_png_file(const std::string& filename, int planes) {
        unique_ptr result = nullptr;
        std::unique_ptr<image<T>> image_2d = nullptr;

        auto image_allocator = [&result, &image_2d, planes](int width, int height,
                                                            std::vector<uint8_t*>* row_pointers) -> bool {
            if (height % planes != 0) {
                // We should fail harder here...
                return false;
            }
            if ((result = std::make_unique<gls::image_3d<T>>(width, height / planes, planes)) == nullptr) {
                return false;
            }
            if ((image_2d = std::make_unique<gls::image<T>>(width, height, result->pixels())) == nullptr) {
                return false;
            }
            for (int i = 0; i < height; ++i) {
                (*row_pointers)[i] = (uint8_t*)(*image_2d)[i];
            }
            return true;
        };

        gls::read_png_file(filename, T::channels, T::bit_depth, image_allocator);

        return result;
    }

    // Write image to PNG file
    int write_png_file(const std::string& filename) const {
        image<T> image_2d(image<T>::width, image<T>::height * depth, image<T>::pixels());
        return image_2d.write_png_file(filename);
    }
};

}  // namespace gls

#endif /* gls_image_h */