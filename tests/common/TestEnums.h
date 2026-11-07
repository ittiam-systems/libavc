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

namespace libavc {
namespace test {

enum class Format {
  yuv400p,      // Planar YUV Luma only
  yuv420p,      // Planar YUV format
  yuv420sp_uv,  // Semi-planar YUV format
  yuv420sp_vu,  // Semi-planar YUV alternative format
  yuv422p,      // Planar YUV alternative format
  yuv444p,      // Planar YUV format
};

enum class Arch {
  defaultArch,
  generic,
  armv7,
  armv8,
  x86_sse42,
  x86_avx2,
};

}  // namespace test
}  // namespace libavc
