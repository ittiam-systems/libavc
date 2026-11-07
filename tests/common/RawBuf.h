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
#include <vector>

#include "TestEnums.h"

class RawBuf {
 public:
  RawBuf() = default;
  ~RawBuf() = default;

  // Allocates buffer with aligned strides for each plane based on video format.
  bool allocBuffer(size_t width, size_t height, size_t bitDepth,
                   libavc::test::Format format, size_t alignment = 1);
  void freeBuffer();

  uint8_t* planeData(size_t planeIdx);
  const uint8_t* planeData(size_t planeIdx) const;

  size_t stride(size_t planeIdx) const { return mStrides[planeIdx]; }
  size_t planeWidth(size_t planeIdx) const;
  size_t planeHeight(size_t planeIdx) const;
  size_t bitDepth() const { return mBitDepth; }
  size_t numPlanes() const;

  size_t width() const { return mWidth; }
  size_t height() const { return mHeight; }
  libavc::test::Format format() const { return mFormat; }

 private:
  std::vector<uint8_t> mStorage;
  size_t mWidth = 0;
  size_t mHeight = 0;
  size_t mBitDepth = 0;
  libavc::test::Format mFormat = libavc::test::Format::yuv420p;
  size_t mStrides[3] = {0, 0, 0};
  size_t mOffsets[3] = {0, 0, 0};
};
