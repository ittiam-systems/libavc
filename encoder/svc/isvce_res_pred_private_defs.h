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
*  isvc_res_pred_private_defs.h
*
* @brief
*  Contains datatype and macro definitions used exclusively in
*  residual prediction
*
*******************************************************************************
*/

#ifndef _ISVCE_RES_PRED_PRIVATE_DEFS_H_
#define _ISVCE_RES_PRED_PRIVATE_DEFS_H_

#include "ih264_typedefs.h"
#include "isvc_defs.h"
#include "isvc_structs.h"

#define REF_ARRAY_MAX_WIDTH (MB_SIZE + 6)

#define REF_ARRAY_MAX_HEIGHT (MB_SIZE + 6)

typedef UWORD32 FT_GET_SAD_WITH_RES_PRED(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                         buffer_container_t *ps_res, UWORD32 u4_mb_wd,
                                         UWORD32 u4_mb_ht);

typedef void FT_RESIDUAL_SAMPLER(coordinates_t *ps_ref_array_positions,
                                 coordinates_t *ps_ref_array_phases, buffer_container_t *ps_inp,
                                 buffer_container_t *ps_out, buffer_container_t *ps_scratch,
                                 UWORD32 u4_ref_nnz, UWORD8 u1_ref_tx_size);

/* Structs */
/* Offsets, etc used for residual upsampling and interpolation */
/* Derived as per 'G.8.6.3.2', and 'G.8.6.3.3' for all MB's once during init */
typedef struct res_pred_mb_state_t
{
    coordinates_t s_offsets;

    coordinates_t s_ref_array_dims;

    coordinates_t *ps_ref_array_positions;

    coordinates_t *ps_ref_array_phases;
} res_pred_mb_state_t;

typedef struct res_pred_layer_state_t
{
    layer_resampler_props_t *ps_luma_props;

    layer_resampler_props_t *ps_chroma_props;

    res_pred_mb_state_t *ps_luma_mb_states;

    res_pred_mb_state_t *ps_chroma_mb_states;

    WORD8 *pi1_mb_mode;

    WORD32 i4_mb_mode_stride;

} res_pred_layer_state_t;

typedef struct res_pred_mem_store_t
{
    buffer_container_t s_scratch;

} res_pred_mem_store_t;

typedef struct res_pred_state_t
{
    /* Array of size numSpatialLayers */
    res_pred_layer_state_t *ps_layer_state;

    res_pred_mem_store_t s_mem_store;

    FT_RESIDUAL_SAMPLER *apf_residual_samplers[NUM_COMPONENTS];

    FT_GET_SAD_WITH_RES_PRED *pf_get_sad_with_residual_pred;

    UWORD8 *pu1_ref_x_ptr_incr; /*!< buffer to store the reference
                        array ptr increments for
                        operand 2 of interpolation
                    */
    UWORD8 *pu1_ref_y_ptr_incr; /*!< buffer to store the reference
                          array ptr increments for
                          operand 2 of interpolation
                          */

} res_pred_state_t;

/* C declarations */
extern FT_RESIDUAL_SAMPLER isvce_luma_residual_sampler_2x;
extern FT_RESIDUAL_SAMPLER isvce_chroma_residual_sampler_2x;
extern FT_GET_SAD_WITH_RES_PRED isvce_get_sad_with_residual_pred;

/* SSE42 declarations */
extern FT_RESIDUAL_SAMPLER isvce_luma_residual_sampler_2x_sse42;
extern FT_GET_SAD_WITH_RES_PRED isvce_get_sad_with_residual_pred_sse42;

/* NEON declarations */
extern FT_RESIDUAL_SAMPLER isvce_luma_residual_sampler_2x_neon;
extern FT_GET_SAD_WITH_RES_PRED isvce_get_sad_with_residual_pred_neon;

#endif
