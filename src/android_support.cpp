//
// Created by Fabio Riccardi on 9/9/21.
//

#include <string>
#include <vector>
#include <map>

#include <jni.h>

#include <android/asset_manager.h>
#include <android/asset_manager_jni.h>

#include "android_support.h"
#include "gls_logging.h"

static const char *TAG = "AndroidSupport";

std::string toString(JNIEnv *env, jstring jStr) {
    const char *cStr = env->GetStringUTFChars(jStr, nullptr);
    std::string str = std::string(cStr);
    env->ReleaseStringUTFChars(jStr, cStr);
    return str;
}

void loadOpenCLShaders(JNIEnv *env, jobject assetManager,
                       std::map<std::string, std::string> *shaders) {
    AAssetManager *mgr = AAssetManager_fromJava(env, assetManager);
    AAssetDir *assetDir = AAssetManager_openDir(mgr, "");
    const char *asset_filename;
    while ((asset_filename = AAssetDir_getNextFileName(assetDir)) != nullptr) {
        std::string filename(asset_filename);
        if (filename.ends_with(".cl")) {
            LOG_INFO(TAG) << "Loading OpenCL shader: " << filename << std::endl;
            std::string filePath = /* "OpenCL/" + */ filename;
            AAsset *asset = AAssetManager_open(mgr, filePath.c_str(), AASSET_MODE_BUFFER);
            off_t assetLength = AAsset_getLength(asset);
            const void *assetBuffer = AAsset_getBuffer(asset);
            std::string source((char *) assetBuffer, (char *) assetBuffer + assetLength);
            (*shaders)[filename] = source;
            AAsset_close(asset);
        }
    }
    AAssetDir_close(assetDir);
}

void loadOpenCLBytecode(JNIEnv *env, jobject assetManager,
                       std::map<std::string, std::vector<unsigned char>> *bytecodes) {
    AAssetManager *mgr = AAssetManager_fromJava(env, assetManager);
    AAssetDir *assetDir = AAssetManager_openDir(mgr, "");
    const char *asset_filename;
    while ((asset_filename = AAssetDir_getNextFileName(assetDir)) != nullptr) {
        std::string filename(asset_filename);
        if (filename.ends_with(".o")) {
            LOG_INFO(TAG) << "Loading OpenCL binary shader: " << filename << std::endl;
            std::string filePath = /* "OpenCL/" + */ filename;
            AAsset *asset = AAssetManager_open(mgr, filePath.c_str(), AASSET_MODE_BUFFER);
            off_t assetLength = AAsset_getLength(asset);
            const void *assetBuffer = AAsset_getBuffer(asset);
            std::vector<unsigned char> bytecode((char *) assetBuffer, (char *) assetBuffer + assetLength);
            (*bytecodes)[filename] = bytecode;
            AAsset_close(asset);
        }
    }
    AAssetDir_close(assetDir);
}

void loadResourceData(JNIEnv *env, jobject assetManager,
                      std::vector <std::byte> *modelData,
                      const std::string &resourceName) {
    AAssetManager *mgr = AAssetManager_fromJava(env, assetManager);
    AAssetDir *assetDir = AAssetManager_openDir(mgr, "");
    const char *filename;
    while ((filename = AAssetDir_getNextFileName(assetDir)) != nullptr) {
        if (std::string(filename) == resourceName) {
            AAsset *asset = AAssetManager_open(mgr, filename, AASSET_MODE_BUFFER);
            off_t assetLength = AAsset_getLength(asset);
            const void *assetBuffer = AAsset_getBuffer(asset);
            *modelData = std::vector<std::byte>((std::byte *) assetBuffer,
                                                (std::byte *) assetBuffer + assetLength);
            AAsset_close(asset);
            break;
        }
    }
    AAssetDir_close(assetDir);
}
