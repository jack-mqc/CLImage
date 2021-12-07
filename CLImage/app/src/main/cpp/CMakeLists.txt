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
        climage/cl_error.cpp
        climage/gls_logging.cpp
        climage-test.cpp)

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
        ${log-lib})

add_executable(
        CLImageTest
        climage/gls_image.hpp
        climage/gls_image_jpeg.cpp
        main.cpp
)

target_link_libraries( # Specifies the target library.
        CLImageTest
        libjpg
        libpng
        libz
        ${log-lib}
        )

target_compile_options( CLImageTest PRIVATE -Wall -Werror -DUSE_IOSTREAM_LOG -DCOMMAND_LINE_TOOL )

target_include_directories( CLImageTest PRIVATE ${CMAKE_SOURCE_DIR}/headers )
target_include_directories( CLImageTest PRIVATE ${CMAKE_SOURCE_DIR}/climage )
