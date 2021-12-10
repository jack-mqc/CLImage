//
// Created by Fabio Riccardi on 9/10/21.
//

#include <fstream>
#include <sstream>

#include "cl_pipeline.h"
#include "gls_cl.hpp"
#include "gls_logging.h"

static const char* TAG = "ShaderCompiler";

extern "C" int __android_log_print(int prio, const char *tag, const char *fmt, ...) {
    return 0;
}

int main(int argc, const char *argv[])
{
    std::cout << "OpenCL Shader Compiler." << std::endl;

    if (argc > 1 && strcmp(argv[1], "-help") == 0) {
        std::cout << "Usage: cat shader.cl | " << argv[0] << " outfile" << std::endl;
        return 0;
    }

    // Read shader source from stdin
    std::stringbuf buffer;
    std::ostream os (&buffer);
    for (std::string line; std::getline(std::cin, line);) {
        os << line << std::endl;
    }

    try {
        // initializeOpenCL();
        cl::Context context = gls::getContext();
        cl::Program program(buffer.str());
        if (gls::buildProgram(program) != 0) {
            return -1;
        }

        std::vector<std::vector<unsigned char>> binaries;
        cl_int result = program.getInfo(CL_PROGRAM_BINARIES, &binaries);
        if (result != CL_SUCCESS) {
            LOG_INFO(TAG) << "CL_PROGRAM_BINARIES returned: " << gls::clStatusToString(result) << std::endl;
            return -1;
        }

        gls::SaveBinaryFile(argc > 1 ? argv[1] : "binaryShader.o", binaries[0]);

        return 0;
    } catch (cl::Error& err) {
        LOG_ERROR(TAG) << "Caught Exception: " << std::string(err.what()) << " - " << gls::clStatusToString(err.err())
                       << std::endl;
        return -1;
    }
}
