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

#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "BitsBuf.h"
#include "BitsFile.h"
#include "RawBuf.h"
#include "RawFile.h"
#include "TestEnums.h"

class DecHelper {
 public:
  class Builder {
   public:
    Builder() = default;
    ~Builder() = default;

    Builder& setFormat(libavc::test::Format format) {
      mFormat = format;
      return *this;
    }
    Builder& setCores(size_t cores) {
      mCores = cores;
      return *this;
    }
    Builder& setArch(libavc::test::Arch arch) {
      mArch = arch;
      return *this;
    }
    Builder& setInputFilePath(const std::string& path) {
      mInputFilePath = path;
      return *this;
    }
    Builder& setOutputFilePath(const std::string& path) {
      mOutputFilePath = path;
      return *this;
    }
    Builder& setRefMd5Path(const std::string& path) {
      mRefMd5Path = path;
      return *this;
    }
    Builder& setOutputMd5Path(const std::string& path) {
      mOutputMd5Path = path;
      return *this;
    }
    Builder& setOutputSyncFramesPath(const std::string& path) {
      mOutputSyncFramesPath = path;
      return *this;
    }
    Builder& setDisableRefRetention(bool disableRefRetention) {
      mDisableRefRetention = disableRefRetention;
      return *this;
    }
    Builder& setNumLoops(size_t numLoops) {
      mNumLoops = numLoops;
      return *this;
    }

    // Validates parameters and returns an initialized DecHelper instance
    std::unique_ptr<DecHelper> build();

   private:
    libavc::test::Format mFormat = libavc::test::Format::yuv420p;
    size_t mCores = 1;
    libavc::test::Arch mArch = libavc::test::Arch::defaultArch;
    std::optional<std::string> mInputFilePath;
    std::optional<std::string> mOutputFilePath;
    std::optional<std::string> mRefMd5Path;
    std::optional<std::string> mOutputMd5Path;
    std::optional<std::string> mOutputSyncFramesPath;
    bool mDisableRefRetention = false;
    size_t mNumLoops = 1;
    friend class DecHelper;
  };

  ~DecHelper() { deinitDecoder(); }

  // Executes the entire decoding process using the configured input/output
  // files. Uses inline GTest assertions to validate frame MD5 checksums.
  bool decodeFile();

  // Low-level chunk-based frame decoding
  bool decodeHeader(size_t& bytesConsumed);
  bool decodeFrame(size_t& bytesConsumed, bool& frameReady, size_t& frameIndex,
                   bool noInput = false);
  bool flushDecoder(size_t& frameIndex);
  void deinitDecoder();

  // Computes plane-by-plane MD5 checksum of a raw frame (ignoring stride
  // padding row-by-row)
  std::string computeFrameMd5(const RawBuf& buf);

  size_t width() const { return mWidth; }
  size_t height() const { return mHeight; }
  bool isHeaderDecoded() const { return mHeaderDecoded; }
  size_t numFramesDecoded() const { return mFramesDecoded; }
  size_t numRefMd5s() const { return mRefMd5s.size(); }

 private:
  DecHelper() = default;  // Private constructor, instantiated only via Builder
  bool initDecoder();
  void processDecodedFrame(RawBuf& buf, size_t& frameIndex);

  void* mCodec = nullptr;  // Opaque pointer to internal iv_obj_t
  bool mHeaderDecoded = false;
  bool mInFlushMode = false;
  size_t mWidth = 0;
  size_t mHeight = 0;
  size_t mBitDepth = 8;
  libavc::test::Format mFormat = libavc::test::Format::yuv420p;
  size_t mCores = 1;
  libavc::test::Arch mArch = libavc::test::Arch::defaultArch;

  std::optional<std::string> mInputFilePath;
  std::optional<std::string> mOutputFilePath;
  std::optional<std::string> mRefMd5Path;
  std::vector<std::string> mRefMd5s;
  std::ofstream mOutputMd5File;
  std::ofstream mOutputSyncFramesFile;
  bool mDisableRefRetention = false;
  size_t mNumLoops = 1;
  size_t mFramesDecoded = 0;

  BitsFile mBitsFile;
  RawFile mOutFile;

  BitsBuf mInputBuf;
  RawBuf mOutputBuf;
};
