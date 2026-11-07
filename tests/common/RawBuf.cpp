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

#include "RawBuf.h"

#include <stdexcept>

size_t RawBuf::numPlanes() const {
  switch (mFormat) {
    case libavc::test::Format::yuv400p:
      return 1;
    case libavc::test::Format::yuv420sp_uv:
    case libavc::test::Format::yuv420sp_vu:
      return 2;
    case libavc::test::Format::yuv420p:
    case libavc::test::Format::yuv422p:
    case libavc::test::Format::yuv444p:
      return 3;
  }
  return 0;
}

size_t RawBuf::planeWidth(size_t planeIdx) const {
  if (planeIdx == 0) {
    return mWidth;
  }
  if (planeIdx >= numPlanes()) {
    return 0;
  }
  switch (mFormat) {
    case libavc::test::Format::yuv400p:
      return 0;
    case libavc::test::Format::yuv420p:
    case libavc::test::Format::yuv422p:
      return mWidth / 2;
    case libavc::test::Format::yuv420sp_uv:
    case libavc::test::Format::yuv420sp_vu:
      return mWidth;  // Interleaved UV has same byte-width as Y luma
    case libavc::test::Format::yuv444p:
      return mWidth;
  }
  return 0;
}

size_t RawBuf::planeHeight(size_t planeIdx) const {
  if (planeIdx == 0) {
    return mHeight;
  }
  if (planeIdx >= numPlanes()) {
    return 0;
  }
  switch (mFormat) {
    case libavc::test::Format::yuv400p:
      return 0;
    case libavc::test::Format::yuv420p:
    case libavc::test::Format::yuv420sp_uv:
    case libavc::test::Format::yuv420sp_vu:
      return mHeight / 2;
    case libavc::test::Format::yuv422p:
    case libavc::test::Format::yuv444p:
      return mHeight;
  }
  return 0;
}

bool RawBuf::allocBuffer(size_t width, size_t height, size_t bitDepth,
                         libavc::test::Format format, size_t alignment) {
  mWidth = width;
  mHeight = height;
  mBitDepth = bitDepth ? bitDepth : 8;
  mFormat = format;

  size_t planesCount = numPlanes();
  size_t totalSize = 0;

  for (size_t i = 0; i < 3; ++i) {
    mStrides[i] = 0;
    mOffsets[i] = 0;
  }

  size_t pixelSize = (mBitDepth + 7) / 8;
  for (size_t i = 0; i < planesCount; ++i) {
    size_t w = planeWidth(i);
    size_t h = planeHeight(i);

    // Align stride to alignment boundary
    size_t alignedStride = (w + alignment - 1) & ~(alignment - 1);
    mStrides[i] = alignedStride;
    mOffsets[i] = totalSize;
    totalSize += alignedStride * h * pixelSize;
  }

  mStorage.resize(totalSize);
  return true;
}

void RawBuf::freeBuffer() {
  mStorage.clear();
  mStorage.shrink_to_fit();
  mWidth = 0;
  mHeight = 0;
  for (size_t i = 0; i < 3; ++i) {
    mStrides[i] = 0;
    mOffsets[i] = 0;
  }
}

uint8_t* RawBuf::planeData(size_t planeIdx) {
  if (planeIdx >= numPlanes() || mStorage.empty()) {
    return nullptr;
  }
  return mStorage.data() + mOffsets[planeIdx];
}

const uint8_t* RawBuf::planeData(size_t planeIdx) const {
  if (planeIdx >= numPlanes() || mStorage.empty()) {
    return nullptr;
  }
  return mStorage.data() + mOffsets[planeIdx];
}
