#!/bin/sh

SDK_DIR=$1
SRC_DIR=$2
DST_DIR=$3

ADB=${SDK_DIR}/platform-tools/adb
ShaderCompiler='app/build/intermediates/cmake/arm8Release/obj/arm64-v8a/ShaderCompiler'

$ADB push $ShaderCompiler /data/local/tmp > /dev/null

for path in "${SRC_DIR}"/*.cl
do
    filename="${path##*/}"
    shader_name="${filename%.*}"
    object_file="/data/local/tmp/${shader_name}.o"

    echo Compiling OpenCL shader file: "${filename}"

    ("${ADB}" shell /data/local/tmp/ShaderCompiler "${object_file}") < "${path}"
    ${ADB} pull "${object_file}" ${DST_DIR}/"${shader_name}.o" > /dev/null
done
