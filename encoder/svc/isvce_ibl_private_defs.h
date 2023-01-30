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

/**
*******************************************************************************
* @file
*  isvce_intra_pred_private_defs.h
*
* @brief
*  Contains datatype and macro definitions used exclusively in
*  residual prediction
*
*******************************************************************************
*/

#ifndef _ISVCE_IBL_PRIVATE_DEFS_H_
#define _ISVCE_IBL_PRIVATE_DEFS_H_

#include "ih264_typedefs.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_structs.h"
#include "isvc_intra_resample.h"

/* Structs */
typedef struct intra_pred_mb_state_t
{
    coordinates_t s_offsets;

    coordinates_t s_ref_array_dims;

    WORD32 *pi4_ref_array_positions_x;

    WORD32 *pi4_ref_array_positions_y;

    coordinates_t *ps_ref_array_phases;

    coordinates_t s_min_pos;

    coordinates_t s_max_pos;

} intra_pred_mb_state_t;

typedef struct intra_pred_layer_state_t
{
    layer_resampler_props_t *ps_luma_props;

    layer_resampler_props_t *ps_chroma_props;

    intra_pred_mb_state_t *ps_luma_mb_states;

    intra_pred_mb_state_t *ps_chroma_mb_states;

    WORD8 *pi1_mb_mode;

    WORD32 i4_mb_mode_stride;

    /* buffer to store the reference
       layer data before intra sampling */
    UWORD8 *pu1_refarray_buffer;

    UWORD8 *pu1_refarray_cb;

    UWORD8 *pu1_refarray_cr;

    WORD32 *pi4_temp_interpolation_buffer;

} intra_pred_layer_state_t;

typedef struct intra_pred_state_t
{
    /* Array of size numSpatialLayers */
    intra_pred_layer_state_t *ps_layer_state;

} intra_pred_state_t;

#endif
