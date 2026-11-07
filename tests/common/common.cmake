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

if(ENABLE_TESTS)
  enable_testing()
  include(GoogleTest)
  if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR AND CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
    set(AVC_GTEST_DISCOVERY_MODE DISCOVERY_MODE PRE_TEST)
  else()
    set(AVC_GTEST_DISCOVERY_MODE)
  endif()
  include("${AVC_ROOT}/tests/common/helpers.cmake")

  include(FetchContent)
  FetchContent_Declare(
    avc_test_data
    URL "https://storage.googleapis.com/android_media/external/libavc/tests/AvcTestRes-1.1.zip"
  )
  FetchContent_MakeAvailable(avc_test_data)

  include("${AVC_ROOT}/tests/decoder/DecHelper.cmake")
  include("${AVC_ROOT}/tests/encoder/EncHelper.cmake")
endif()
