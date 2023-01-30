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

#ifndef _ISVCE_DOWNSCALER_PRIVATE_DEFS_H_
#define _ISVCE_DOWNSCALER_PRIVATE_DEFS_H_
#include "ih264_typedefs.h"
#include "isvc_macros.h"
#include "ih264_debug.h"
#include "isvc_structs.h"
#include "isvce_downscaler.h"

/* Macros */
#define DOWNSCALER_Q 16

#define FILTER_COEFF_Q 7

#define NUM_SCALER_FILTER_TAPS 8

#define NUM_SCALER_FILTER_PHASES 8

/* Typedefs */
typedef WORD8 (*FILTER_COEFF_ARRAY)[NUM_SCALER_FILTER_TAPS * 2];

typedef void FT_DOWNSCALER(downscaler_ctxt_t *ps_scaler_state, buffer_container_t *ps_src,
                           buffer_container_t *ps_dst, FILTER_COEFF_ARRAY pai1_filters,
                           UWORD32 u4_blk_wd, UWORD32 u4_blk_ht, UWORD8 u1_is_chroma);

/* Structs */
typedef struct
{
    /**
     * pointer to scratch buf
     */
    void *pv_scratch_buf;

    /**
     * initial offset while calculating input pixel location
     */
    WORD32 i4_init_offset;

    /**
     * increment to the centre pixel in horizontal direction
     */
    UWORD32 u4_horz_increment;

    /**
     * increment to the centre pixel in vertical direction
     */
    UWORD32 u4_vert_increment;

    /**
     * pointer to the filter coefficients
     */
    FILTER_COEFF_ARRAY pai1_filters;

    /**
     * function pointer to the leaf level function for horizontal scaling
     */
    FT_DOWNSCALER *pf_downscaler;

    /**
     * width of the input (highest SVC layer)
     */
    UWORD32 u4_in_wd;

    /**
     * height of the input (highest SVC layer)
     */
    UWORD32 u4_in_ht;

} downscaler_state_t;

static FORCEINLINE UWORD32 get_filter_phase(UWORD32 u4_center_pixel_pos)
{
    UWORD32 au4_phase_binning_pos[NUM_SCALER_FILTER_PHASES + 1];
    UWORD32 i;

    ASSERT(NUM_SCALER_FILTER_PHASES == 8);

    for(i = 0; i < NUM_SCALER_FILTER_PHASES + 1; i++)
    {
        au4_phase_binning_pos[i] = (i << DOWNSCALER_Q) / NUM_SCALER_FILTER_PHASES;
    }

    u4_center_pixel_pos = u4_center_pixel_pos % (1 << DOWNSCALER_Q);

    for(i = 0; i < NUM_SCALER_FILTER_PHASES; i++)
    {
        if((u4_center_pixel_pos < au4_phase_binning_pos[i + 1]) &&
           (u4_center_pixel_pos >= au4_phase_binning_pos[i]))
        {
            return i;
        }
    }

    ASSERT(0);

    return 0;
}

/* SSE42 Declarations */
extern FT_DOWNSCALER isvce_horizontal_downscale_and_transpose_sse42;

/* NEON Declarations */
extern FT_DOWNSCALER isvce_horizontal_downscale_and_transpose_neon;

#endif
