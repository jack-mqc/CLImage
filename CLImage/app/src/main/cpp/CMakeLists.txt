# For more information about using CMake with Android Studio, read the
# documentation: https://d.android.com/studio/projects/add-native-code.html

# Sets the minimum version of CMake required to build the native library.

cmake_minimum_required(VERSION 3.10.2)

# Declares and names the project.

project("climage")

add_library( libz STATIC IMPORTED )
set_target_properties( libz PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/jni/${ANDROID_ABI}/libz.a )
add_library( libpng STATIC IMPORTED )
set_target_properties( libpng PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/jni/${ANDROID_ABI}/libpng.a )
add_library( libjpg STATIC IMPORTED )
set_target_properties( libjpg PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/jni/${ANDROID_ABI}/libjpeg-turbo.a )

add_library( libtiff STATIC IMPORTED )
set_target_properties( libtiff PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/jni/${ANDROID_ABI}/libtiff.a )
add_library( libtiffxx STATIC IMPORTED )
set_target_properties( libtiffxx PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/jni/${ANDROID_ABI}/libtiffxx.a )

# Creates and names a library, sets it as either STATIC
# or SHARED, and provides the relative paths to its source code.
# You can define multiple libraries, and CMake builds them for you.
# Gradle automatically packages shared libraries with your APK.

add_library( # Sets the name of the library.
        climage

        # Sets the library as a shared library.
        SHARED

        # Provides a relative path to your source file(s).
        climage/gls_icd_wrapper.cpp
        climage/gls_cl_error.cpp
        climage/gls_logging.cpp
        climage/gls_cl.cpp
        climage/gls_android_support.cpp
        cl_pipeline.cpp
        climage-test.cpp)

target_compile_options( climage PRIVATE -Wall -Werror -DUSE_ASSET_MANAGER )

# Add header files for opencl
target_include_directories( climage PRIVATE ${CMAKE_SOURCE_DIR}/headers )
target_include_directories( climage PRIVATE ${CMAKE_SOURCE_DIR}/climage )

# Searches for a specified prebuilt library and stores the path as a
# variable. Because CMake includes system libraries in the search path by
# default, you only need to specify the name of the public NDK library
# you want to add. CMake verifies that the library exists before
# completing its build.

find_library( # Sets the name of the path variable.
        log-lib

        # Specifies the name of the NDK library that
        # you want CMake to locate.
        log)

# Specifies libraries CMake should link to your target library. You
# can link multiple libraries, such as libraries you define in this
# build script, prebuilt third-party libraries, or system libraries.

target_link_libraries( # Specifies the target library.
        climage

        # Links the target library to the log library
        # included in the NDK.
        ${log-lib}
        -ljnigraphics
        -landroid)

add_executable(
        CLImageTest
        climage/gls_image_png.cpp
        climage/gls_image_jpeg.cpp
        climage/gls_image_tiff.cpp
        climage/gls_icd_wrapper.cpp
        climage/gls_cl_error.cpp
        climage/gls_logging.cpp
        climage/gls_cl.cpp
        cl_pipeline.cpp
        main.cpp
)

target_link_libraries( # Specifies the target library.
        CLImageTest
        libjpg
        libpng
        libz
        libtiff
        libtiffxx
        ${log-lib})

target_compile_options( CLImageTest PRIVATE -Wall -Werror -DUSE_IOSTREAM_LOG )

target_include_directories( CLImageTest PRIVATE ${CMAKE_SOURCE_DIR}/headers )
target_include_directories( CLImageTest PRIVATE ${CMAKE_SOURCE_DIR}/climage )

add_executable(
        RawPipeline
        climage/gls_image_png.cpp
        climage/gls_image_jpeg.cpp
        climage/gls_image_tiff.cpp
        climage/gls_tiff_metadata.cpp
        climage/gls_dng_lossless_jpeg.cpp
        climage/gls_icd_wrapper.cpp
        climage/gls_cl_error.cpp
        climage/gls_logging.cpp
        climage/gls_cl.cpp
        climage/gls_color_science.cpp
        cl_pipeline.cpp
        demosaic.cpp
        raw_pipeline.cpp
        filesystem/operations.cpp
        filesystem/directory_iterator.cpp
)

target_link_libraries( # Specifies the target library.
        RawPipeline
        libjpg
        libpng
        libz
        libtiff
        libtiffxx
        ${log-lib})

target_compile_options( RawPipeline PRIVATE -Wall -Werror -DUSE_IOSTREAM_LOG )

target_include_directories( RawPipeline PRIVATE ${CMAKE_SOURCE_DIR}/headers )
target_include_directories( RawPipeline PRIVATE ${CMAKE_SOURCE_DIR}/climage )

add_executable(
        ShaderCompiler
        climage/gls_icd_wrapper.cpp
        climage/gls_cl_error.cpp
        climage/gls_logging.cpp
        climage/gls_cl.cpp
        ShaderCompiler.cpp
)

target_compile_options( ShaderCompiler PRIVATE -Wall -Werror -DUSE_IOSTREAM_LOG )

target_include_directories( ShaderCompiler PRIVATE ${CMAKE_SOURCE_DIR}/headers )
target_include_directories( ShaderCompiler PRIVATE ${CMAKE_SOURCE_DIR}/climage )
