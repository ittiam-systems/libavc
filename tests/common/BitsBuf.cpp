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

#include "BitsBuf.h"

bool BitsBuf::allocBuffer(size_t capacity) {
  mStorage.resize(capacity);
  mSize = 0;
  mOffset = 0;
  return true;
}

void BitsBuf::freeBuffer() {
  mStorage.clear();
  mStorage.shrink_to_fit();
  mSize = 0;
  mOffset = 0;
}
