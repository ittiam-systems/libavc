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

class BitsBuf {
 public:
  BitsBuf() = default;
  ~BitsBuf() = default;

  bool allocBuffer(size_t capacity);
  void freeBuffer();

  uint8_t* data() { return mStorage.data(); }
  const uint8_t* data() const { return mStorage.data(); }
  size_t size() const { return mSize; }
  size_t capacity() const { return mStorage.size(); }
  size_t offset() const { return mOffset; }

  void setSize(size_t size) { mSize = size; }
  void setOffset(size_t offset) { mOffset = offset; }

 private:
  std::vector<uint8_t> mStorage;
  size_t mSize = 0;
  size_t mOffset = 0;
};
