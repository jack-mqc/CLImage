//
// Created by Fabio Riccardi on 8/23/21.
//

#ifndef LEICAPHOTOAPPX_ANDROID_SUPPORT_H
#define LEICAPHOTOAPPX_ANDROID_SUPPORT_H

#include <span>
#include <map>
#include <android/bitmap.h>

#include <exception>

class exception : public std::exception {
    const std::string _message;

public:
    explicit exception(std::string message) : _message(std::move(message)) {}

    [[nodiscard]] const char* what() const noexcept override { return _message.c_str(); }
};

std::string toString(JNIEnv *env, jstring jStr);

void loadResourceData(JNIEnv *env, jobject assetManager,
                      std::vector<std::byte> *modelData,
                      const std::string &resourceName);

void loadOpenCLShaders(JNIEnv *env, jobject assetManager,
                       std::map<std::string, std::string> *shaders);

void loadOpenCLBytecode(JNIEnv *env, jobject assetManager,
                        std::map<std::string, std::vector<unsigned char>> *bytecodes);

template<class T>
class JavaArray : public std::vector<T> {
public:
    JavaArray(JNIEnv *env, jfloatArray array) {
        jsize arrayLength = env->GetArrayLength(array);
        auto *arrayData = (T *) env->GetPrimitiveArrayCritical(array, nullptr);
        this->assign(arrayData, arrayData + arrayLength);
        env->ReleasePrimitiveArrayCritical(array, arrayData, 0);
    }
};

template<class T>
class JavaArrayCritical : public std::span<T> {
    JNIEnv *_env;
    jfloatArray _array;
public:
    JavaArrayCritical(JNIEnv *env, jfloatArray array) : _env(env), _array(array), std::span<T>(
            (T *) env->GetPrimitiveArrayCritical(array, nullptr),
            env->GetArrayLength(array)
    ) {}

    ~JavaArrayCritical() {
        _env->ReleasePrimitiveArrayCritical(_array, this->data(), 0);
    }
};

struct AndroidBitmap {
    JNIEnv *_env;
    jobject _bitmap;
    AndroidBitmapInfo _info;

    AndroidBitmap(JNIEnv *env, jobject bitmap) : _env(env), _bitmap(bitmap) {
        int status = AndroidBitmap_getInfo(env, bitmap, &this->_info);
        if (status != ANDROID_BITMAP_RESULT_SUCCESS) {
            throw exception("Failed accessing Android Bitmap object: " + std::to_string(status));
        }
    }

    [[nodiscard]] const AndroidBitmapInfo &info() const {
        return _info;
    }

    template <typename T>
    [[nodiscard]] std::span<T> lockPixels() const {
        uint8_t *bitmapData = nullptr;
        int status = AndroidBitmap_lockPixels(_env, _bitmap, (void **) &bitmapData);
        if (status != ANDROID_BITMAP_RESULT_SUCCESS) {
            throw exception("AndroidBitmapInfo::lockPixels failure: " + std::to_string(status));
        }
        int pixelSize = 0;
        switch (_info.format) {
            case ANDROID_BITMAP_FORMAT_RGBA_8888:
                pixelSize = 4;
                break;
            case ANDROID_BITMAP_FORMAT_RGB_565:
            case ANDROID_BITMAP_FORMAT_RGBA_4444:
                pixelSize = 2;
                break;
            case ANDROID_BITMAP_FORMAT_A_8:
                pixelSize = 1;
                break;
            case ANDROID_BITMAP_FORMAT_RGBA_F16:
                pixelSize = 8;
                break;
            default:
                throw exception("Unexpected Bitmap format: " + std::to_string(_info.format));
        }
        assert(pixelSize * _info.width == _info.stride);
        return std::span((T *) bitmapData, pixelSize * _info.width * _info.height / sizeof(T));
    }

    void unLockPixels() const {
        int status = AndroidBitmap_unlockPixels(_env, _bitmap);
        if (status != ANDROID_BITMAP_RESULT_SUCCESS) {
            throw exception("AndroidBitmapInfo::unLockPixels failure: " + std::to_string(status));
        }
    }
};

#endif //LEICAPHOTOAPPX_ANDROID_SUPPORT_H
