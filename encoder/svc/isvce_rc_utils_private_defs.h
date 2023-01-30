/******************************************************************************
 *
 * Copyright (C) 2022 The Android Open Source Project
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

#ifndef _ISVCE_RC_UTILS_PRIVATE_DEFS_H_
#define _ISVCE_RC_UTILS_PRIVATE_DEFS_H_

#include "ih264_typedefs.h"
#include "isvc_structs.h"
#include "isvce_rc_utils.h"

/* Macros */
#define WT_LUMA_GPP 4

#define WT_TOTAL_GPP 6

/* Typedefs */
typedef DOUBLE FT_GET_GPP(yuv_buf_props_t *ps_input_buf);

/* Structs */
typedef struct
{
    /**
     * function pointer to the leaf level function for get gpp
     */
    FT_GET_GPP *pf_get_gpp;

} svc_rc_utils_state_t;

/* SSE42 Declarations */
extern FT_GET_GPP isvce_get_gpp_sse42;

/* NEON Declarations */
extern FT_GET_GPP isvce_get_gpp_neon;

#endif
