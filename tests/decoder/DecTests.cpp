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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "DecHelper.h"
#include "TestEnums.h"

namespace libavc {
namespace test {

struct DecodeStreamConfig {
  std::string inputFile;
  Format format;
  std::string refMd5File;
};

void PrintTo(const DecodeStreamConfig& config, ::std::ostream* os) {
  *os << config.inputFile;
}

// Helper to resolve absolute file path for test assets
static std::string getFullPath(const std::string& inputFile) {
  if (const char* envPath = std::getenv("AVC_TEST_DIR")) {
    return std::string(envPath) + "/" + inputFile;
  }
  return inputFile;
}

class DecTestFixture
    : public ::testing::TestWithParam<std::tuple<DecodeStreamConfig, size_t>> {
 protected:
  DecTestFixture() = default;
  ~DecTestFixture() override = default;
};

TEST_P(DecTestFixture, DecodeVerify) {
  auto [config, cores] = GetParam();

  std::string inputPath = getFullPath(config.inputFile);
  std::string refMd5Path = getFullPath(config.refMd5File);

  auto builder = DecHelper::Builder()
                     .setInputFilePath(inputPath)
                     .setRefMd5Path(refMd5Path)
                     .setFormat(config.format)
                     .setCores(cores);

  std::unique_ptr<DecHelper> helper = builder.build();
  ASSERT_NE(helper, nullptr)
      << "Failed to build DecHelper for input: " << inputPath;
  EXPECT_TRUE(helper->decodeFile())
      << "Decoding failed for input: " << inputPath;
}

static const std::vector<DecodeStreamConfig> kDecodeStreams = {
    {"bbb_176x144_baseline.h264", Format::yuv420p, "bbb_176x144_baseline_md5.txt"},
    {"bbb_176x144_main.h264", Format::yuv420p, "bbb_176x144_main_md5.txt"},
    {"bbb_176x144_high.h264", Format::yuv420p, "bbb_176x144_high_md5.txt"},
};

class DecStreamTestFixture
    : public ::testing::TestWithParam<DecodeStreamConfig> {
 protected:
  DecStreamTestFixture() = default;
  ~DecStreamTestFixture() override = default;
};

TEST_P(DecStreamTestFixture, DecodeDisableRefRetention) {
  struct TempFileGuard {
    std::string path;
    ~TempFileGuard() { std::remove(path.c_str()); }
  };

  constexpr size_t kNumLoops = 2;
  constexpr size_t kNumCores = 3;
  const auto& config = GetParam();
  std::string inputPath = getFullPath(config.inputFile);
  std::string sanitizedName = config.inputFile;
  for (char& c : sanitizedName) {
    if (c == '/' || c == '\\') {
      c = '_';
    }
  }
  std::string intraBitstreamPath =
      ::testing::TempDir() + "disable_ref_retention_intra_" + sanitizedName;
  std::string intraMd5Path =
      ::testing::TempDir() + "disable_ref_retention_intra_" + sanitizedName +
      ".md5";
  TempFileGuard intraBitstreamGuard{intraBitstreamPath};
  TempFileGuard intraMd5Guard{intraMd5Path};

  {
    auto extractHelper = DecHelper::Builder()
                             .setInputFilePath(inputPath)
                             .setOutputSyncFramesPath(intraBitstreamPath)
                             .setFormat(config.format)
                             .build();
    ASSERT_NE(extractHelper, nullptr)
        << "Failed to build DecHelper for intra extraction: " << inputPath;
    ASSERT_TRUE(extractHelper->decodeFile())
        << "Failed to extract intra frames for input: " << inputPath;
  }
  {
    auto md5Helper = DecHelper::Builder()
                         .setInputFilePath(intraBitstreamPath)
                         .setOutputMd5Path(intraMd5Path)
                         .setFormat(config.format)
                         .setNumLoops(kNumLoops)
                         .build();
    ASSERT_NE(md5Helper, nullptr)
        << "Failed to build DecHelper for intra MD5 generation: " << inputPath;
    ASSERT_TRUE(md5Helper->decodeFile())
        << "Failed to generate intra MD5 file for input: " << inputPath;
  }

  auto builder = DecHelper::Builder()
                     .setInputFilePath(inputPath)
                     .setRefMd5Path(intraMd5Path)
                     .setFormat(config.format)
                     .setCores(kNumCores)
                     .setDisableRefRetention(true)
                     .setNumLoops(kNumLoops);

  std::unique_ptr<DecHelper> helper = builder.build();
  ASSERT_NE(helper, nullptr)
      << "Failed to build DecHelper for input: " << inputPath;
  EXPECT_TRUE(helper->decodeFile())
      << "Decoding failed with disable_ref_retention enabled for input: "
      << inputPath;
  EXPECT_GT(helper->numRefMd5s(), 0u)
      << "Expected intra reference MD5s for input: " << inputPath;
  EXPECT_EQ(helper->numFramesDecoded(), helper->numRefMd5s())
      << "Decoded intra frame count mismatch for input: " << inputPath;
}

INSTANTIATE_TEST_SUITE_P(DecoderRegression, DecTestFixture,
                         ::testing::Combine(::testing::ValuesIn(kDecodeStreams),
                                            ::testing::Values(1, 2, 3)));

INSTANTIATE_TEST_SUITE_P(DecoderRegression, DecStreamTestFixture,
                         ::testing::ValuesIn(kDecodeStreams));

}  // namespace test
}  // namespace libavc
