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
#include <cstdio>
#include <string>

#include "BitsBuf.h"

class BitsFile {
 public:
  BitsFile();
  ~BitsFile();

  bool open(const std::string& path, bool readMode);
  void close();

  size_t read(BitsBuf& buf, size_t bytesToRead);
  size_t write(const BitsBuf& buf, size_t bytesToWrite);
  size_t seek(size_t offset);

  bool isEof() const;
  bool isValid() const { return mFile != nullptr; }

 private:
  FILE* mFile = nullptr;
  bool mReadMode = true;
};
