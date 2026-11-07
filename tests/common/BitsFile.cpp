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

#include "BitsFile.h"

#include <algorithm>

BitsFile::BitsFile() : mFile(nullptr), mReadMode(true) {}

BitsFile::~BitsFile() { close(); }

bool BitsFile::open(const std::string& path, bool readMode) {
  close();
  mReadMode = readMode;
  mFile = fopen(path.c_str(), mReadMode ? "rb" : "wb");
  return mFile != nullptr;
}

size_t BitsFile::seek(size_t offset) {
  if (!mFile) {
    return 0;
  }
  return fseek(mFile, offset, SEEK_SET);
}

void BitsFile::close() {
  if (mFile) {
    fclose(mFile);
    mFile = nullptr;
  }
}

size_t BitsFile::read(BitsBuf& buf, size_t bytesToRead) {
  if (!mFile || !mReadMode) {
    return 0;
  }
  size_t toRead = std::min(bytesToRead, buf.capacity());
  size_t readBytes = fread(buf.data(), 1, toRead, mFile);
  buf.setSize(readBytes);
  buf.setOffset(0);
  return readBytes;
}

size_t BitsFile::write(const BitsBuf& buf, size_t bytesToWrite) {
  if (!mFile || mReadMode) {
    return 0;
  }
  size_t toWrite = std::min(bytesToWrite, buf.size());
  size_t written = fwrite(buf.data() + buf.offset(), 1, toWrite, mFile);
  fflush(mFile);
  return written;
}

bool BitsFile::isEof() const {
  if (!mFile) {
    return true;
  }
  return feof(mFile) != 0;
}
