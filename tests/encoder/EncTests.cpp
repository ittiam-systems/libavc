/******************************************************************************
 *
 * Copyright (C) 2026 Ittiam Systems Pvt Ltd, Bangalore
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 ******************************************************************************/

#include <gtest/gtest.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "DecHelper.h"
#include "EncHelper.h"
#include "Md5Wrapper.h"
#include "TestEnums.h"

namespace libavc {
namespace test {

struct EncodeStreamConfig {
  std::string filename;
  size_t width;
  size_t height;
  Format format;
};

void PrintTo(const EncodeStreamConfig& config, ::std::ostream* os) {
  *os << config.filename;
}

// Helper to resolve absolute file path for test assets
static std::string getFullPath(const std::string& filename) {
  if (const char* envPath = std::getenv("AVC_TEST_DIR")) {
    return std::string(envPath) + "/" + filename;
  }
  return filename;
}

class EncTestFixture
    : public ::testing::TestWithParam<std::tuple<EncodeStreamConfig, size_t>> {
 protected:
  EncTestFixture() = default;
  ~EncTestFixture() override = default;
};

TEST_P(EncTestFixture, EncodeVerify) {
  auto [config, cores] = GetParam();

  std::string inputPath = getFullPath(config.filename);
  std::string outBitstreamPath = ::testing::TempDir() + "temp_output.h264";
  std::string reconMd5Path = ::testing::TempDir() + "temp_recon_md5.txt";

  std::ifstream f(inputPath);
  if (!f.good()) {
    GTEST_FAIL() << "Input file not found: " << inputPath;
  }
  f.close();

  // Clean up any stale outputs
  std::remove(outBitstreamPath.c_str());

  auto encBuilder = EncHelper::Builder()
                        .setInputFilePath(inputPath)
                        .setOutFilePath(outBitstreamPath)
                        .setOutputReconMd5Path(reconMd5Path)
                        .setFormat(config.format)
                        .setResolution(config.width, config.height)
                        .setFrameRate(30000, 1000)
                        .setCores(cores);

  std::unique_ptr<EncHelper> helper = encBuilder.build();
  ASSERT_NE(helper, nullptr)
      << "Failed to build EncHelper for input: " << inputPath;
  EXPECT_TRUE(helper->encodeFile())
      << "Encoding failed for input: " << inputPath;

  auto decBuilder = DecHelper::Builder()
                        .setInputFilePath(outBitstreamPath)
                        .setRefMd5Path(reconMd5Path)
                        .setFormat(libavc::test::Format::yuv420p)
                        .setCores(cores);
  std::unique_ptr<DecHelper> decHelper = decBuilder.build();
  ASSERT_NE(decHelper, nullptr)
      << "Failed to build DecHelper for input: " << outBitstreamPath;
  EXPECT_TRUE(decHelper->decodeFile())
      << "Decoding failed for input: " << outBitstreamPath;
}

static const std::vector<EncodeStreamConfig> kEncodeStreams = {
    {"bbb_176x144_yuv420.yuv", 176, 144, Format::yuv420p}};

INSTANTIATE_TEST_SUITE_P(EncoderRegression, EncTestFixture,
                         ::testing::Combine(::testing::ValuesIn(kEncodeStreams),
                                            ::testing::Values(1, 2, 3, 4)));

}  // namespace test
}  // namespace libavc
