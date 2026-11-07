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

# Define the encoder test helper library
add_library(libavc_enc_helper STATIC
  ${AVC_ROOT}/tests/encoder/EncHelper.cpp
)

# Export include directories
target_include_directories(libavc_enc_helper PUBLIC
  ${AVC_ROOT}/tests/encoder
  ${AVC_ROOT}/tests/common
  ${AVC_ROOT}/common
  ${AVC_ROOT}/encoder
)

# Link against native encoder and common test helper library
target_link_libraries(libavc_enc_helper PUBLIC
  libavcenc
  libavc_test_common
)

# Define the encoder regression test executable
libavc_add_executable(
  avc_enc_tests
  libavc_enc_helper
  SOURCES ${AVC_ROOT}/tests/encoder/EncTests.cpp
  LIBS libavc_dec_helper GTest::gtest_main
)
target_compile_definitions(avc_enc_tests PRIVATE "PROJECT_DIR=\"${AVC_ROOT}\"")
gtest_discover_tests(
  avc_enc_tests
  ${AVC_GTEST_DISCOVERY_MODE}
  PROPERTIES ENVIRONMENT "AVC_TEST_DIR=${avc_test_data_SOURCE_DIR}"
)
