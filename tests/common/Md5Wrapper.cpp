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

#include "Md5Wrapper.h"

#if defined(USE_LIGHTWEIGHT_MD5)
#define MD5_IMPLEMENTATION
#include "md5.h"
typedef md5_context MD5_CTX;
#define MD5_Init md5_init
#define MD5_Update md5_update
#define MD5_Final(digest, ctx) md5_finalize(ctx, digest)
#define MD5_DIGEST_LENGTH MD5_SIZE
#else
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#include <openssl/md5.h>
#endif

#include <iomanip>
#include <sstream>

Md5Wrapper::Md5Wrapper() { mContext = new MD5_CTX(); }

Md5Wrapper::~Md5Wrapper() {
  if (mContext) {
    delete static_cast<MD5_CTX*>(mContext);
  }
}

void Md5Wrapper::init() { MD5_Init(static_cast<MD5_CTX*>(mContext)); }

void Md5Wrapper::update(const uint8_t* data, size_t size) {
  MD5_Update(static_cast<MD5_CTX*>(mContext), data, size);
}

std::string Md5Wrapper::finalize() {
  unsigned char digest[MD5_DIGEST_LENGTH];
  MD5_Final(digest, static_cast<MD5_CTX*>(mContext));

  std::stringstream ss;
  for (int i = 0; i < MD5_DIGEST_LENGTH; ++i) {
    ss << std::hex << std::setw(2) << std::setfill('0')
       << static_cast<int>(digest[i]);
  }
  return ss.str();
}

#if !defined(USE_LIGHTWEIGHT_MD5) && (defined(__GNUC__) || defined(__clang__))
#pragma GCC diagnostic pop
#endif
