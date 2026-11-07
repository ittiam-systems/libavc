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

#include "RawBuf.h"

class RawFile {
 public:
  RawFile();
  ~RawFile();

  bool open(const std::string& path, bool readMode);
  void close();

  // Reads a single raw frame from disk, mapping contiguous rows into
  // stride-aligned RawBuf planes.
  bool read(RawBuf& buf);

  // Writes a single raw frame to disk, stripping stride padding to write
  // contiguous rows.
  bool write(const RawBuf& buf);

  bool isEof() const;
  bool isValid() const { return mFile != nullptr; }

 private:
  FILE* mFile = nullptr;
  bool mReadMode = true;
};
