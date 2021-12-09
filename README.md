# CLImage
Simple OpenCL Image Processing with C++ on Android and macOS

CLImage provides a modern no-fuss C++ approach to using OpenCL on Android and Unix based systems.

On Android it also provides a bridge to the OpenCL libraries installed by the device manufacturer,
allowing the use of powerful header based wrappers such as opencl.hpp

Two example projects are provided to demonstrate climage's use on both Android and macOS.

To accompany the OpenCL/opencl.hpp programming environment, a set of image management classes
allows to easily create, destroy, modify image objects on the CPU and on the GPU.

As an example:

    auto inputImage = gls::image<gls::rgba_pixel>::read_png_file(filePath);

reads a PNG image from the file system and returns a uinique_ptr to an in memory RGBA representation
of the image.

A variety of image formats are available for the most common image layouts (Luma, LumaAlpha, RGB, RGBA)
and data formats (uint8, uint16, float32 and float16), image formats and pixel data types are defined in
gls_image.hpp

To create an image that can be directly accessed by OpenCL we can use:

    gls::cl_image_2d<gls::rgba_pixel> clInputImage(context, *inputImage);

This will create an image-like object that acts as a wrapper to the OpenCL texture image, initialized with
the geometry, data type and data content of inputImage.

All image objects are typed, making it possible to use the C++ type system to robustly interface between C++
and OpenCL/opencl.hpp

To get the underlying OpenCL Image2D object, the getImage2D() accessor is provided. The following code fragmnent
illustrates how to use the cl_image types with opencl.hpp

    int blur(const gls::cl_image_2d<gls::rgba_pixel>& input, gls::cl_image_2d<gls::rgba_pixel>* output) {
        try {
            // Load the shader source
            const auto blurProgram = gls::loadOpenCLProgram("blur");

            // Bind the kernel's parameters
            auto blurKernel = cl::KernelFunctor<cl::Image2D,  // input
                                                cl::Image2D   // output
                                                >(*blurProgram, "blur");

            // Schedule the kernel on the GPU
            blurKernel(gls::buildEnqueueArgs(output->width, output->height), input.getImage2D(), output->getImage2D());
            return 0;
        } catch (cl::Error& err) {
            LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                           << std::endl;
            return -1;
        }
    }

To retrieve data from a cl_image, both memory mapping and memory copy are available, for instance:

    auto outputImage = clOutputImage.mapImage();
    outputImage.write_png_file("output.png");
    clOutputImage.unmapImage(outputImage);

Creates a zero-copy gls::image instance wrapped around the cl_image OpenCL representation, the image is written to the
file system as a PNG file and then unmapped from memory.
