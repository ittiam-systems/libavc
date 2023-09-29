/******************************************************************************
 *
 * Copyright (C) 2015 The Android Open Source Project
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
 *****************************************************************************
 * Originally developed and contributed by Ittiam Systems Pvt. Ltd, Bangalore
*/

#ifndef _RC_TRACE_SUPPORT_H_
#define _RC_TRACE_SUPPORT_H_

#if DEBUG_RC
#define TRACE_PRINTF(...)                                                   \
{                                                                           \
    printf("\n[RC DBG] %s/%d:: ", __FUNCTION__, __LINE__);                  \
    printf(__VA_ARGS__);                                                    \
}
#else
#define TRACE_PRINTF(...) {}
#endif

#endif // _RC_TRACE_SUPPORT_H_
