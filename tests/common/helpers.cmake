#******************************************************************************
#
# Copyright (C) 2026 Ittiam Systems Pvt Ltd, Bangalore
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at:
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#******************************************************************************

# Determine MD5 implementation based on target
if(("${SYSTEM_PROCESSOR}" STREQUAL "aarch64") OR
   ("${SYSTEM_PROCESSOR}" STREQUAL "aarch32") OR
   ("${SYSTEM_PROCESSOR}" STREQUAL "arm64"))
  set(USE_LIGHTWEIGHT_MD5 TRUE)
else()
  find_package(OpenSSL REQUIRED)
  set(USE_OPENSSL_MD5 TRUE)
endif()

if(USE_LIGHTWEIGHT_MD5)
  # Fetch lightweight MD5 implementation (mackron/md5)
  include(FetchContent)
  FetchContent_Declare(
    mackron_md5
    GIT_REPOSITORY https://github.com/mackron/md5.git
    GIT_TAG        master
  )
  FetchContent_MakeAvailable(mackron_md5)
endif()

# Define the common test helper library
add_library(libavc_test_common STATIC
  ${AVC_ROOT}/tests/common/BitsBuf.cpp
  ${AVC_ROOT}/tests/common/RawBuf.cpp
  ${AVC_ROOT}/tests/common/Md5Wrapper.cpp
  ${AVC_ROOT}/tests/common/BitsFile.cpp
  ${AVC_ROOT}/tests/common/RawFile.cpp
)

# Export include directories so other libraries/tests can find headers
target_include_directories(libavc_test_common PUBLIC
  ${AVC_ROOT}/tests/common
)

if(USE_LIGHTWEIGHT_MD5)
  target_compile_definitions(libavc_test_common PUBLIC USE_LIGHTWEIGHT_MD5)
  target_include_directories(libavc_test_common PUBLIC ${mackron_md5_SOURCE_DIR})
endif()

if(USE_OPENSSL_MD5)
  target_link_libraries(libavc_test_common PRIVATE OpenSSL::Crypto)
endif()

target_link_libraries(libavc_test_common PUBLIC gtest)
