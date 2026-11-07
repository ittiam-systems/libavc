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

#include "RawFile.h"

RawFile::RawFile() : mFile(nullptr), mReadMode(true) {}

RawFile::~RawFile() { close(); }

bool RawFile::open(const std::string& path, bool readMode) {
  close();
  mReadMode = readMode;
  mFile = fopen(path.c_str(), mReadMode ? "rb" : "wb");
  return mFile != nullptr;
}

void RawFile::close() {
  if (mFile) {
    fclose(mFile);
    mFile = nullptr;
  }
}

bool RawFile::read(RawBuf& buf) {
  if (!mFile || !mReadMode) {
    return false;
  }

  size_t planesCount = buf.numPlanes();
  for (size_t i = 0; i < planesCount; ++i) {
    uint8_t* planePtr = buf.planeData(i);
    size_t planeStride = buf.stride(i);
    size_t planeW = buf.planeWidth(i);
    size_t planeH = buf.planeHeight(i);

    if (!planePtr) {
      return false;
    }

    for (size_t row = 0; row < planeH; ++row) {
      size_t readBytes =
          fread(planePtr + (row * planeStride), 1, planeW, mFile);
      if (readBytes != planeW) {
        return false;  // EOF or read error
      }
    }
  }
  return true;
}

bool RawFile::write(const RawBuf& buf) {
  if (!mFile || mReadMode) {
    return false;
  }

  size_t planesCount = buf.numPlanes();
  for (size_t i = 0; i < planesCount; ++i) {
    const uint8_t* planePtr = buf.planeData(i);
    size_t planeStride = buf.stride(i);
    size_t planeW = buf.planeWidth(i);
    size_t planeH = buf.planeHeight(i);

    if (!planePtr) {
      return false;
    }

    for (size_t row = 0; row < planeH; ++row) {
      size_t writtenBytes =
          fwrite(planePtr + (row * planeStride), 1, planeW, mFile);
      if (writtenBytes != planeW) {
        return false;  // Write error
      }
    }
  }
  fflush(mFile);
  return true;
}

bool RawFile::isEof() const {
  if (!mFile) {
    return true;
  }
  return feof(mFile) != 0;
}
