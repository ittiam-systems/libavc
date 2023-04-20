/******************************************************************************
 *
 * Copyright (C) 2021 The Android Open Source Project
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

#ifndef _IMVCD_DEFS_H_
#define _IMVCD_DEFS_H_

#include <stdint.h>

#include "ih264_typedefs.h"
#include "imvc_defs.h"
#include "ih264d_defs.h"

#define MVC_MAX_REF_PICS MAX(16 * LOG2_MAX_NUM_VIEWS, 2 * H264_MAX_REF_PICS)

/* Set to identify between actual ref pic with valid u1_pic_buf_id
   and replicated IVP ref pic! Ensure that MAX_VAL_PIC_BUF_ID-MAX_NUM_VIEWS
   is still greater than any possible value of u1_pic_buf_id */
#define IVP_PIC_BUF_ID UINT8_MAX

/* In FGC SEI
   - Worst-case bits for all elements before 'num_intensity_intervals_minus1' = 47
   - Worst-case bits for all elements before 'film_grain_characteristics_repetition_period', not
     including elements from previous line = 3 * (8 + 3 + 256 * (8 + 8 + 8 * 16)) = 110625
   - Worst-case bits for 'film_grain_characteristics_repetition_period' = 30
   Total of (47 + 110625 + 30) = 110702 byte */
#define MAX_FGC_SEI_SIZE 110702

#define MIN_BITSTREAMS_BUF_SIZE (MAX_FGC_SEI_SIZE + 256000)

#endif
