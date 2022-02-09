cmake_minimum_required(VERSION 3.9.1)
project(libavc C CXX)
enable_language(ASM)

set(AVC_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
set(AVC_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}")

if("${AVC_ROOT}" STREQUAL "${AVC_CONFIG_DIR}")
  message(
    FATAL_ERROR
      "Building from within the libavc source tree is not supported.\n"
      "Hint: Run these commands\n"
      "$ rm -rf CMakeCache.txt CMakeFiles\n"
      "$ mkdir -p ./build\n"
      "$ cd ./build\n"
      "And re-run CMake from the build directory.")
endif()

set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

set(CMAKE_STATIC_LIBRARY_PREFIX "")

if(${CMAKE_CROSSCOMPILING})
  set(THREADS_PTHREAD_ARG
      "2"
      CACHE STRING
            "Refer https://gitlab.kitware.com/cmake/cmake/-/issues/16920" FORCE)
endif()

include("${AVC_ROOT}/cmake/utils.cmake")

libavc_add_compile_options()
libavc_add_definitions()
libavc_set_link_libraries()

include("${AVC_ROOT}/common/common.cmake")
include("${AVC_ROOT}/decoder/libavcdec.cmake")
include("${AVC_ROOT}/encoder/libavcenc.cmake")

include("${AVC_ROOT}/test/decoder/avcdec.cmake")
include("${AVC_ROOT}/test/encoder/avcenc.cmake")

include("${AVC_ROOT}/fuzzer/avc_dec_fuzzer.cmake")
include("${AVC_ROOT}/fuzzer/avc_enc_fuzzer.cmake")
