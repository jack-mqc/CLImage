#!/bin/sh

SDK_DIR=$1
SRC_DIR=$2
DST_DIR=$3

ADB=${SDK_DIR}/platform-tools/adb
ShaderCompiler='app/build/intermediates/cmake/arm8Release/obj/arm64-v8a/ShaderCompiler'

$ADB push $ShaderCompiler /data/local/tmp

for path in "${SRC_DIR}"/*.cl
do
	("${ADB}" shell /data/local/tmp/ShaderCompiler /data/local/tmp/compiled.o) < "${path}"
	filename="${path##*/}"
	shader_name="${filename%.*}"
	${ADB} pull /data/local/tmp/compiled.o ${DST_DIR}/"${shader_name}.o"
done
