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

# Define the decoder test helper library
add_library(libavc_dec_helper STATIC
  ${AVC_ROOT}/tests/decoder/DecHelper.cpp
)

# Export include directories
target_include_directories(libavc_dec_helper PUBLIC
  ${AVC_ROOT}/tests/decoder
  ${AVC_ROOT}/tests/common
  ${AVC_ROOT}/common
  ${AVC_ROOT}/decoder
)

# Link against the native decoder, common test helper, and GTest/GMock
target_link_libraries(libavc_dec_helper PUBLIC
  libavcdec
  libavc_test_common
  GTest::gtest
  GTest::gmock
)

# Define the decoder regression test executable
libavc_add_executable(
  avc_dec_tests
  libavc_dec_helper
  SOURCES ${AVC_ROOT}/tests/decoder/DecTests.cpp
  LIBS GTest::gtest_main
)
target_compile_definitions(avc_dec_tests PRIVATE "PROJECT_DIR=\"${AVC_ROOT}\"")
gtest_discover_tests(
  avc_dec_tests
  ${AVC_GTEST_DISCOVERY_MODE}
  PROPERTIES ENVIRONMENT "AVC_TEST_DIR=${avc_test_data_SOURCE_DIR}"
)
