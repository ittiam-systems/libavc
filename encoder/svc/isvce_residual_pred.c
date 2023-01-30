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
*  isvce_residual_pred.c
*
* @brief
*  Contains functions used for SVC residual prediction
*
*******************************************************************************
*/
#include <stdint.h>
#include <math.h>

#include "ih264_typedefs.h"
#include "iv2.h"
#include "isvc_macros.h"
#include "ih264_debug.h"
#include "isvc_defs.h"
#include "isvc_structs.h"
#include "isvce_defs.h"
#include "isvce_structs.h"
#include "isvce_res_pred_private_defs.h"
#include "isvce_residual_pred.h"
#include "isvce_utils.h"
#include "isvc_defs.h"

void isvce_chroma_residual_sampler_2x(coordinates_t *ps_ref_array_positions,
                                      coordinates_t *ps_ref_array_phases,
                                      buffer_container_t *ps_inp, buffer_container_t *ps_out,
                                      buffer_container_t *ps_scratch, UWORD32 u4_ref_nnz,
                                      UWORD8 u1_ref_tx_size)
{
    WORD32 i4_i;
    WORD16 *pi2_ref_data_byte;
    WORD32 *pi4_ref_array;
    WORD32 i4_phase1, i4_phase2;

    WORD16 *pi2_inp_data = ps_inp->pv_data;
    WORD16 *pi2_out_res = ps_out->pv_data;
    WORD32 i4_inp_data_stride = ps_inp->i4_data_stride;
    WORD32 i4_out_res_stride = ps_out->i4_data_stride;

    UNUSED(u4_ref_nnz);

    UNUSED(ps_ref_array_positions);
    UNUSED(u1_ref_tx_size);

    /* For 2x scaling, offsets always point to TL pixel outside MB */
    /* Hence, refTransBlkIdc will be different and since phase */
    /* for first refArray pos for horiz filtering samples > 8, */
    /* first row and first column from the refArray is never used */
    pi2_inp_data += 2 + i4_inp_data_stride;

    pi2_ref_data_byte = pi2_inp_data;

    i4_phase1 = ps_ref_array_phases[0].i4_abscissa;
    i4_phase2 = ps_ref_array_phases[1].i4_abscissa;

    ASSERT(i4_phase1 >= 8);

    pi4_ref_array = (WORD32 *) ps_scratch->pv_data;

    for(i4_i = 0; i4_i < BLK_SIZE; i4_i++)
    {
        WORD16 i2_coeff1, i2_coeff2;

        i2_coeff1 = (WORD16) (pi2_ref_data_byte[0]);

        /* populate the first inter sample */
        *pi4_ref_array++ = i2_coeff1 << 4;

        {
            /* unroll count 1 */
            i2_coeff2 = (WORD16) (pi2_ref_data_byte[2]);

            /* populate 2 samples based on current coeffs */
            *pi4_ref_array++ = ((16 - i4_phase2) * i2_coeff1 + i4_phase2 * i2_coeff2);

            /* unroll count 2 */
            *pi4_ref_array++ = ((16 - i4_phase1) * i2_coeff1 + i4_phase1 * i2_coeff2);

            /* unroll count 3 */
            i2_coeff1 = (WORD16) (pi2_ref_data_byte[4]);

            /* populate 2 samples based on current coeffs */
            *pi4_ref_array++ = ((16 - i4_phase2) * i2_coeff2 + i4_phase2 * i2_coeff1);

            /* unroll count 4 */
            *pi4_ref_array++ = ((16 - i4_phase1) * i2_coeff2 + i4_phase1 * i2_coeff1);

            /* unroll count 5 */
            i2_coeff2 = (WORD16) (pi2_ref_data_byte[6]);

            /* populate 2 samples based on current coeffs */
            *pi4_ref_array++ = ((16 - i4_phase2) * i2_coeff1 + i4_phase2 * i2_coeff2);

            /* unroll count 6 */
            *pi4_ref_array++ = ((16 - i4_phase1) * i2_coeff1 + i4_phase1 * i2_coeff2);
        }

        /* populate the last inter sample */
        *pi4_ref_array++ = i2_coeff2 << 4;

        /* vertical loop uopdates */
        pi2_ref_data_byte = pi2_inp_data + ((i4_i + 1) * i4_inp_data_stride);
    }

    /* ----------- Vertical Interpolation ---------------- */
    pi4_ref_array = (WORD32 *) ps_scratch->pv_data;

    i4_phase1 = ps_ref_array_phases[0].i4_ordinate;
    i4_phase2 = ps_ref_array_phases[2].i4_ordinate;

    for(i4_i = 0; i4_i < BLK8x8SIZE; i4_i++)
    {
        WORD16 *pi2_out;
        WORD32 *pi4_ref_array_temp;
        WORD32 i4_horz_samp_1, i4_horz_samp_2;
        pi2_out = pi2_out_res;
        pi4_ref_array_temp = pi4_ref_array;

        /* populate the first inter sample */
        i4_horz_samp_1 = *pi4_ref_array_temp;
        pi4_ref_array_temp += BLK8x8SIZE;
        *pi2_out = (i4_horz_samp_1 + 8) >> 4;
        pi2_out += i4_out_res_stride;

        {
            /* unroll count 1 */
            i4_horz_samp_2 = *pi4_ref_array_temp;
            pi4_ref_array_temp += BLK8x8SIZE;

            /* populate 2 samples based on current coeffs */
            *pi2_out = ((16 - i4_phase2) * i4_horz_samp_1 + i4_phase2 * i4_horz_samp_2 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 2 */
            *pi2_out = ((16 - i4_phase1) * i4_horz_samp_1 + i4_phase1 * i4_horz_samp_2 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 3 */
            i4_horz_samp_1 = *pi4_ref_array_temp;
            pi4_ref_array_temp += BLK8x8SIZE;

            /* populate 2 samples based on current coeffs */
            *pi2_out = ((16 - i4_phase2) * i4_horz_samp_2 + i4_phase2 * i4_horz_samp_1 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 4 */
            *pi2_out = ((16 - i4_phase1) * i4_horz_samp_2 + i4_phase1 * i4_horz_samp_1 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 5 */
            i4_horz_samp_2 = *pi4_ref_array_temp;

            /* populate 2 samples based on current coeffs */
            *pi2_out = ((16 - i4_phase2) * i4_horz_samp_1 + i4_phase2 * i4_horz_samp_2 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 6 */
            *pi2_out = ((16 - i4_phase1) * i4_horz_samp_1 + i4_phase1 * i4_horz_samp_2 + 128) >> 8;
            pi2_out += i4_out_res_stride;
        }

        /* populate the last inter sample */
        *pi2_out = (i4_horz_samp_2 + 8) >> 4;

        /* horizontal loop updates */
        pi4_ref_array++;
        pi2_out_res += 2;
    }
}

void isvce_luma_residual_sampler_2x(coordinates_t *ps_ref_array_positions,
                                    coordinates_t *ps_ref_array_phases, buffer_container_t *ps_inp,
                                    buffer_container_t *ps_out, buffer_container_t *ps_scratch,
                                    UWORD32 u4_ref_nnz, UWORD8 u1_ref_tx_size)
{
    WORD16 *pi2_inp_data = ps_inp->pv_data;
    WORD16 *pi2_out_res = ps_out->pv_data;
    WORD32 i4_inp_data_stride = ps_inp->i4_data_stride;
    WORD32 i4_out_res_stride = ps_out->i4_data_stride;
    WORD16 *pi2_refarray_buffer = ps_scratch->pv_data;
    WORD32 i4_blk_ctr;

    UNUSED(ps_ref_array_positions);
    UNUSED(ps_ref_array_phases);

    /* For 2x scaling, offsets always point to TL pixel outside MB */
    /* Hence, refTransBlkIdc will be different and since phase */
    /* for first refArray pos for horiz filtering samples > 8, */
    /* first row and first column from the refArray is never used */
    pi2_inp_data += 1 + i4_inp_data_stride;

    if((u1_ref_tx_size) && (0 != u4_ref_nnz))
    {
        WORD16 *pi2_ref_data_byte;
        WORD32 *pi4_ref_array;
        WORD32 i4_i, i4_j;

        pi2_ref_data_byte = pi2_inp_data;

        /* ----------- Horizontal Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) pi2_refarray_buffer;

        for(i4_i = 0; i4_i < BLK8x8SIZE; i4_i++)
        {
            WORD16 i2_coeff1, i2_coeff2;

            i2_coeff1 = (WORD16) (*pi2_ref_data_byte++);

            /* populate the first inter sample */
            *pi4_ref_array++ = i2_coeff1 << 2;

            for(i4_j = 0; i4_j < 14; i4_j += 2)
            {
                i2_coeff2 = (WORD16) (*pi2_ref_data_byte++);

                /* populate 2 samples based on current coeffs */
                *pi4_ref_array++ = ((i2_coeff1 << 1) + (i2_coeff1) + (i2_coeff2));

                *pi4_ref_array++ = ((i2_coeff2 << 1) + (i2_coeff2) + (i2_coeff1));

                /* store the coeff 2 to coeff 1 */
                /* (used in next iteration)     */
                i2_coeff1 = i2_coeff2;
            }

            /* populate the last inter sample */
            *pi4_ref_array++ = i2_coeff1 << 2;

            /* vertical loop uopdates */
            pi2_ref_data_byte = pi2_inp_data + ((i4_i + 1) * i4_inp_data_stride);
        }

        /* ----------- Vertical Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) pi2_refarray_buffer;

        for(i4_i = 0; i4_i < MB_SIZE; i4_i++)
        {
            WORD32 *pi4_ref_array_temp;
            WORD16 *pi2_out;
            WORD32 i4_horz_samp_1, i4_horz_samp_2;

            pi4_ref_array_temp = pi4_ref_array;
            pi2_out = pi2_out_res;
            i4_horz_samp_1 = *pi4_ref_array_temp;

            /* populate the first inter sample */
            *pi2_out = (i4_horz_samp_1 + 2) >> 2;
            pi2_out += i4_out_res_stride;

            for(i4_j = 0; i4_j < 14; i4_j += 2)
            {
                pi4_ref_array_temp += MB_SIZE;
                i4_horz_samp_2 = *pi4_ref_array_temp;

                /* populate 2 samples based on current coeffs */
                *pi2_out = ((i4_horz_samp_1 << 1) + (i4_horz_samp_1) + (i4_horz_samp_2) + 8) >> 4;
                pi2_out += i4_out_res_stride;

                *pi2_out = ((i4_horz_samp_2 << 1) + (i4_horz_samp_2) + (i4_horz_samp_1) + 8) >> 4;
                pi2_out += i4_out_res_stride;

                /* store the coeff 2 to coeff 1 */
                /* (used in next iteration)     */
                i4_horz_samp_1 = i4_horz_samp_2;
            }

            /* populate the first inter sample */
            *pi2_out = (i4_horz_samp_1 + 2) >> 2;

            /* horizontal loop updates */
            pi4_ref_array++;
            pi2_out_res++;
        }
    }
    else
    {
        /* ----------------------------------------------------------------- */
        /* LOOP over number of blocks                                        */
        /* ----------------------------------------------------------------- */
        for(i4_blk_ctr = 0; i4_blk_ctr < BLK_SIZE; i4_blk_ctr++)
        {
            WORD16 *pi2_ref_data_byte;
            WORD32 *pi4_ref_array;
            WORD32 i4_i;

            /* if reference layer is not coded then no processing */
            if(0 != (u4_ref_nnz & 0x1))
            {
                pi2_ref_data_byte = pi2_inp_data;

                /* ----------- Horizontal Interpolation ---------------- */
                pi4_ref_array = (WORD32 *) pi2_refarray_buffer;

                for(i4_i = 0; i4_i < BLK_SIZE; i4_i++)
                {
                    WORD16 i2_coeff1, i2_coeff2;

                    i2_coeff1 = (WORD16) (*pi2_ref_data_byte++);

                    /* populate the first inter sample */
                    *pi4_ref_array++ = i2_coeff1 << 2;

                    {
                        i2_coeff2 = (WORD16) (*pi2_ref_data_byte++);

                        /* populate 2 samples based on current coeffs */
                        *pi4_ref_array++ = ((i2_coeff1 << 1) + (i2_coeff1) + (i2_coeff2));

                        *pi4_ref_array++ = ((i2_coeff2 << 1) + (i2_coeff2) + (i2_coeff1));

                        i2_coeff1 = (WORD16) (*pi2_ref_data_byte++);

                        /* populate 2 samples based on current coeffs */
                        *pi4_ref_array++ = ((i2_coeff2 << 1) + (i2_coeff2) + (i2_coeff1));

                        *pi4_ref_array++ = ((i2_coeff1 << 1) + (i2_coeff1) + (i2_coeff2));

                        i2_coeff2 = (WORD16) (*pi2_ref_data_byte++);

                        /* populate 2 samples based on current coeffs */
                        *pi4_ref_array++ = ((i2_coeff1 << 1) + (i2_coeff1) + (i2_coeff2));

                        *pi4_ref_array++ = ((i2_coeff2 << 1) + (i2_coeff2) + (i2_coeff1));
                    }

                    /* populate the last inter sample */
                    *pi4_ref_array++ = i2_coeff2 << 2;

                    /* vertical loop uopdates */
                    pi2_ref_data_byte = pi2_inp_data + ((i4_i + 1) * i4_inp_data_stride);
                }

                /* ----------- Vertical Interpolation ---------------- */
                pi4_ref_array = (WORD32 *) pi2_refarray_buffer;

                for(i4_i = 0; i4_i < BLK8x8SIZE; i4_i++)
                {
                    WORD32 *pi4_ref_array_temp;
                    WORD16 *pi2_out;
                    WORD32 i4_horz_samp_1, i4_horz_samp_2;

                    pi4_ref_array_temp = pi4_ref_array;
                    pi2_out = pi2_out_res;
                    i4_horz_samp_1 = *pi4_ref_array_temp;

                    /* populate the first inter sample */
                    *pi2_out = (i4_horz_samp_1 + 2) >> 2;
                    pi2_out += i4_out_res_stride;

                    {
                        /* unroll loop count 1 */
                        pi4_ref_array_temp += BLK8x8SIZE;
                        i4_horz_samp_2 = *pi4_ref_array_temp;

                        /* populate 2 samples based on current coeffs */
                        *pi2_out =
                            ((i4_horz_samp_1 << 1) + (i4_horz_samp_1) + (i4_horz_samp_2) + 8) >> 4;
                        pi2_out += i4_out_res_stride;

                        *pi2_out =
                            ((i4_horz_samp_2 << 1) + (i4_horz_samp_2) + (i4_horz_samp_1) + 8) >> 4;
                        pi2_out += i4_out_res_stride;

                        /* unroll loop count 2 */
                        pi4_ref_array_temp += BLK8x8SIZE;
                        i4_horz_samp_1 = *pi4_ref_array_temp;

                        /* populate 2 samples based on current coeffs */
                        *pi2_out =
                            ((i4_horz_samp_2 << 1) + (i4_horz_samp_2) + (i4_horz_samp_1) + 8) >> 4;
                        pi2_out += i4_out_res_stride;

                        *pi2_out =
                            ((i4_horz_samp_1 << 1) + (i4_horz_samp_1) + (i4_horz_samp_2) + 8) >> 4;
                        pi2_out += i4_out_res_stride;

                        /* unroll loop count 3 */
                        pi4_ref_array_temp += BLK8x8SIZE;
                        i4_horz_samp_2 = *pi4_ref_array_temp;

                        /* populate 2 samples based on current coeffs */
                        *pi2_out =
                            ((i4_horz_samp_1 << 1) + (i4_horz_samp_1) + (i4_horz_samp_2) + 8) >> 4;
                        pi2_out += i4_out_res_stride;

                        *pi2_out =
                            ((i4_horz_samp_2 << 1) + (i4_horz_samp_2) + (i4_horz_samp_1) + 8) >> 4;
                        pi2_out += i4_out_res_stride;
                    }

                    /* populate the last inter sample */
                    *pi2_out = (i4_horz_samp_2 + 2) >> 2;

                    /* horizontal loop updates */
                    pi4_ref_array++;
                    pi2_out_res++;
                }
            }
            else
            {
                pi2_out_res += BLK8x8SIZE;
            }

            if(1 == i4_blk_ctr)
            {
                pi2_inp_data -= BLK_SIZE;
                pi2_inp_data += (i4_inp_data_stride * BLK_SIZE);
                pi2_out_res -= MB_SIZE;
                pi2_out_res += (i4_out_res_stride * BLK8x8SIZE);
                u4_ref_nnz >>= 2;
            }
            else
            {
                pi2_inp_data += BLK_SIZE;
            }

            u4_ref_nnz >>= 1;
        }
    }
}

/**
*******************************************************************************
*
* @brief
*  Returns size of buffers for storing residual pred ctxt
*
* @param[in] u1_num_spatial_layers
*  Num Spatial Layers
*
* @param[in] d_spatial_res_ratio
*  Resolution Ratio b/w spatial layers
*
* @param[in] u4_wd
*  Input Width
*
* @param[in] u4_ht
*  Input Height
*
* @returns  Size of buffers
*
*******************************************************************************
*/
UWORD32 isvce_get_svc_res_pred_ctxt_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                         UWORD32 u4_wd, UWORD32 u4_ht)
{
    UWORD32 u4_size = 0;

    if(u1_num_spatial_layers > 1)
    {
        WORD32 i;

        u4_size += MAX_PROCESS_CTXT * sizeof(svc_res_pred_ctxt_t);
        u4_size += MAX_PROCESS_CTXT * sizeof(res_pred_state_t);

        /* Mem for storing pred */
        u4_size += MAX_PROCESS_CTXT * MB_SIZE * MB_SIZE * sizeof(WORD16);
        u4_size += MAX_PROCESS_CTXT * MB_SIZE * (MB_SIZE / 2) * sizeof(WORD16);

        /* Mem for storing intermediates */
        u4_size += MAX_PROCESS_CTXT * REF_ARRAY_MAX_WIDTH * REF_ARRAY_MAX_HEIGHT * sizeof(WORD16);

        /* Mem for pu1_ref_x_ptr_incr and pu1_ref_y_ptr_incr*/
        u4_size +=
            2 * MAX_PROCESS_CTXT * REF_ARRAY_MAX_WIDTH * REF_ARRAY_MAX_HEIGHT * sizeof(UWORD8);

        u4_size += MAX_PROCESS_CTXT * u1_num_spatial_layers * sizeof(res_pred_layer_state_t);

        for(i = u1_num_spatial_layers - 1; i >= 1; i--)
        {
            WORD32 i4_layer_luma_wd =
                (WORD32) ((DOUBLE) u4_wd /
                          pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) +
                0.99;
            WORD32 i4_layer_luma_ht =
                ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
            WORD32 i4_layer_luma_mbs = (i4_layer_luma_wd / MB_SIZE) * (i4_layer_luma_ht / MB_SIZE);
            WORD32 i4_layer_u_wd = i4_layer_luma_wd / 2.0 + 0.99;
            WORD32 i4_layer_u_ht = i4_layer_luma_ht / 2.0 + 0.99;
            WORD32 i4_layer_u_mbs =
                (i4_layer_u_wd / (MB_SIZE / 2)) * (i4_layer_u_ht / (MB_SIZE / 2));

            /* ps_luma_mb_states */
            {
                u4_size += i4_layer_luma_mbs * sizeof(res_pred_mb_state_t);

                /* ps_ref_array_positions */
                u4_size +=
                    ((1.5 == d_spatial_res_ratio) ? (i4_layer_luma_mbs * MB_SIZE * MB_SIZE) : 0) *
                    sizeof(coordinates_t);

                /* ps_ref_array_phases */
                u4_size += ((1.5 == d_spatial_res_ratio) ? (i4_layer_luma_mbs * 5) : 0) *
                           sizeof(coordinates_t);
            }

            /* ps_chroma_mb_states */
            {
                u4_size += i4_layer_u_mbs * sizeof(res_pred_mb_state_t);

                /* ps_ref_array_positions */
                u4_size +=
                    ((1.5 == d_spatial_res_ratio) ? (i4_layer_u_mbs * (MB_SIZE / 2) * (MB_SIZE / 2))
                                                  : 0) *
                    sizeof(coordinates_t);

                /* ps_ref_array_phases */
                u4_size += ((1.5 == d_spatial_res_ratio) ? (i4_layer_u_mbs * 5) : 3) *
                           sizeof(coordinates_t);
            }
        }

        for(i = u1_num_spatial_layers - 1; i >= 0; i--)
        {
            WORD32 i4_layer_luma_wd =
                (WORD32) ((DOUBLE) u4_wd /
                          pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) +
                0.99;
            WORD32 i4_layer_luma_ht =
                ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
            WORD32 i4_layer_luma_mbs =
                ((i4_layer_luma_wd / MB_SIZE) + 2) * ((i4_layer_luma_ht / MB_SIZE) + 2);

            /* pi1_mb_mode */
            u4_size += i4_layer_luma_mbs * sizeof(WORD8);
        }
    }
    else
    {
        u4_size += MAX_PROCESS_CTXT * sizeof(yuv_buf_props_t);

        /* Mem for storing pred */
        u4_size += MAX_PROCESS_CTXT * MB_SIZE * MB_SIZE * sizeof(WORD16);
        u4_size += MAX_PROCESS_CTXT * MB_SIZE * (MB_SIZE / 2) * sizeof(WORD16);
    }

    return u4_size;
}

static FORCEINLINE WORD32 isvce_get_scaled_pixel_pos(layer_resampler_props_t *ps_layer_props,
                                                     WORD32 i4_pixel_pos, UWORD8 u1_dim_id)
{
    if(1 == u1_dim_id)
    {
        return (((i4_pixel_pos - ps_layer_props->i4_offset_y) *
                     ((WORD64) ps_layer_props->u4_scale_y) +
                 ps_layer_props->i4_add_y) >>
                (ps_layer_props->u4_shift_y - 4)) -
               ps_layer_props->i4_delta_y;
    }
    else
    {
        return (((i4_pixel_pos - ps_layer_props->i4_offset_x) *
                     ((WORD64) ps_layer_props->u4_scale_x) +
                 ps_layer_props->i4_add_x) >>
                (ps_layer_props->u4_shift_x - 4)) -
               ps_layer_props->i4_delta_x;
    }
}

static FORCEINLINE void isvce_ref_array_pos_and_phase_init_dyadic(
    layer_resampler_props_t *ps_layer_props, res_pred_mb_state_t *ps_mb_state,
    coordinates_t *ps_mb_pos, UWORD8 u1_frame_mbs_only_flag, UWORD8 u1_field_mb_flag,
    UWORD8 u1_ref_layer_frame_mbs_only_flag)
{
    UWORD32 i, j;

    coordinates_t *ps_ref_array_phases = ps_mb_state->ps_ref_array_phases;

    WORD32 i4_x_offset = ps_mb_state->s_offsets.i4_abscissa;
    WORD32 i4_y_offset = ps_mb_state->s_offsets.i4_ordinate;

    for(i = 0; i < 2; i++)
    {
        WORD32 i4_y_ref16;

        WORD32 i4_yc = ps_mb_pos->i4_ordinate * ps_layer_props->u4_mb_ht + i;

        if((0 == u1_frame_mbs_only_flag) || (0 == u1_ref_layer_frame_mbs_only_flag))
        {
            i4_yc = i4_yc >> (1 - u1_field_mb_flag);
        }

        i4_y_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_yc, 1);

        for(j = 0; j < ((0 == i) ? 2 : 1); j++)
        {
            WORD32 i4_x_ref16;

            WORD32 i4_xc = ps_mb_pos->i4_abscissa * ps_layer_props->u4_mb_wd + j;

            i4_x_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_xc, 0);

            ps_ref_array_phases[j + i * 2].i4_abscissa = (i4_x_ref16 - (16 * i4_x_offset)) & 15;
            ps_ref_array_phases[j + i * 2].i4_ordinate = (i4_y_ref16 - (16 * i4_y_offset)) & 15;
        }
    }
}

static FORCEINLINE void isvce_ref_array_pos_and_phase_init(layer_resampler_props_t *ps_layer_props,
                                                           res_pred_mb_state_t *ps_mb_state,
                                                           coordinates_t *ps_mb_pos,
                                                           UWORD8 u1_frame_mbs_only_flag,
                                                           UWORD8 u1_field_mb_flag,
                                                           UWORD8 u1_ref_layer_frame_mbs_only_flag)
{
    UWORD32 i, j;

    coordinates_t *ps_ref_array_positions = ps_mb_state->ps_ref_array_positions;
    coordinates_t *ps_ref_array_phases = ps_mb_state->ps_ref_array_phases;

    WORD32 i4_x_offset = ps_mb_state->s_offsets.i4_abscissa;
    WORD32 i4_y_offset = ps_mb_state->s_offsets.i4_ordinate;
    UWORD32 u4_phase_array_idx = 0;

    for(i = 0; i < ps_layer_props->u4_mb_ht; i++)
    {
        WORD32 i4_y_ref16;

        WORD32 i4_yc = ps_mb_pos->i4_ordinate * ps_layer_props->u4_mb_ht + i;

        if((0 == u1_frame_mbs_only_flag) || (0 == u1_ref_layer_frame_mbs_only_flag))
        {
            i4_yc = i4_yc >> (1 - u1_field_mb_flag);
        }

        i4_y_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_yc, 1);

        for(j = 0; j < ps_layer_props->u4_mb_wd; j++)
        {
            WORD32 i4_x_ref16;

            WORD32 i4_xc = ps_mb_pos->i4_abscissa * ps_layer_props->u4_mb_wd + j;

            i4_x_ref16 = isvce_get_scaled_pixel_pos(ps_layer_props, i4_xc, 0);

            ps_ref_array_positions[j + i * ps_layer_props->u4_mb_wd].i4_abscissa =
                (i4_x_ref16 >> 4) - i4_x_offset;
            ps_ref_array_positions[j + i * ps_layer_props->u4_mb_wd].i4_ordinate =
                (i4_y_ref16 >> 4) - i4_y_offset;

            if(((0 == i) && (j < 3)) || ((0 == j) && (i < 3)))
            {
                ps_ref_array_phases[u4_phase_array_idx].i4_abscissa =
                    (i4_x_ref16 - (16 * i4_x_offset)) & 15;
                ps_ref_array_phases[u4_phase_array_idx].i4_ordinate =
                    (i4_y_ref16 - (16 * i4_y_offset)) & 15;

                u4_phase_array_idx++;
            }
        }
    }
}

static void isvce_res_pred_layer_state_init(res_pred_layer_state_t *ps_layer_state,
                                            DOUBLE d_spatial_res_ratio, UWORD32 u4_wd,
                                            UWORD32 u4_ht, IV_COLOR_FORMAT_T e_color_format)
{
    UWORD32 i, j, k;

    const UWORD8 u1_ref_layer_field_pic_flag = 0;
    const UWORD8 u1_field_pic_flag = 0;
    const UWORD8 u1_frame_mbs_only_flag = 1;
    const UWORD8 u1_ref_layer_frame_mbs_only_flag = 1;
    const UWORD8 u1_field_mb_flag = 0;

    ASSERT((IV_YUV_420P == e_color_format) || (IV_YUV_420SP_UV == e_color_format));

    UNUSED(e_color_format);

    for(i = 0; i < 2; i++)
    {
        res_pred_mb_state_t *ps_mb_states;
        layer_resampler_props_t *ps_layer_props;

        UWORD32 u4_wd_in_mbs;
        UWORD32 u4_ht_in_mbs;

        UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) i));
        UWORD32 u4_ref_wd = (u4_wd / d_spatial_res_ratio);
        UWORD32 u4_ref_ht = (u4_ht / d_spatial_res_ratio) * (1 + u1_ref_layer_field_pic_flag);
        UWORD32 u4_scaled_wd = u4_wd;
        UWORD32 u4_scaled_ht = u4_ht * (1 + u1_field_pic_flag);

        ps_mb_states =
            u1_is_chroma ? ps_layer_state->ps_chroma_mb_states : ps_layer_state->ps_luma_mb_states;
        ps_layer_props =
            u1_is_chroma ? ps_layer_state->ps_chroma_props : ps_layer_state->ps_luma_props;

        u4_ref_wd = u4_ref_wd >> u1_is_chroma;
        u4_ref_ht = u4_ref_ht >> u1_is_chroma;
        u4_scaled_wd = u4_scaled_wd >> u1_is_chroma;
        u4_scaled_ht = u4_scaled_ht >> u1_is_chroma;

        u4_wd_in_mbs = u4_scaled_wd / ps_layer_props->u4_mb_wd;
        u4_ht_in_mbs = u4_scaled_ht / ps_layer_props->u4_mb_ht;

        for(j = 0; j < u4_ht_in_mbs; j++)
        {
            WORD32 i4_y_refmin16;
            WORD32 i4_y_refmax16;
            WORD32 i4_y_offset;

            i4_y_refmin16 =
                isvce_get_scaled_pixel_pos(ps_layer_props, j * ps_layer_props->u4_mb_ht, 1);
            i4_y_refmax16 = isvce_get_scaled_pixel_pos(
                ps_layer_props, j * ps_layer_props->u4_mb_ht + ps_layer_props->u4_mb_ht - 1, 1);
            i4_y_offset = i4_y_refmin16 >> 4;

            for(k = 0; k < u4_wd_in_mbs; k++)
            {
                WORD32 i4_x_refmin16;
                WORD32 i4_x_refmax16;
                WORD32 i4_x_offset;

                coordinates_t s_mb_pos = {k, j};

                i4_x_refmin16 =
                    isvce_get_scaled_pixel_pos(ps_layer_props, k * ps_layer_props->u4_mb_wd, 0);
                i4_x_refmax16 = isvce_get_scaled_pixel_pos(
                    ps_layer_props, k * ps_layer_props->u4_mb_wd + ps_layer_props->u4_mb_wd - 1, 0);
                i4_x_offset = i4_x_refmin16 >> 4;

                ps_mb_states[k + j * u4_wd_in_mbs].s_offsets.i4_abscissa = i4_x_offset;
                ps_mb_states[k + j * u4_wd_in_mbs].s_offsets.i4_ordinate = i4_y_offset;
                ps_mb_states[k + j * u4_wd_in_mbs].s_ref_array_dims.i4_abscissa =
                    (i4_x_refmax16 >> 4) - i4_x_offset + 2;
                ps_mb_states[k + j * u4_wd_in_mbs].s_ref_array_dims.i4_ordinate =
                    (i4_y_refmax16 >> 4) - i4_y_offset + 2;

                if((0 == k) && (0 == j) && (2 == d_spatial_res_ratio) && u1_is_chroma)
                {
                    isvce_ref_array_pos_and_phase_init_dyadic(
                        ps_layer_props, &ps_mb_states[k + j * u4_wd_in_mbs], &s_mb_pos,
                        u1_frame_mbs_only_flag, u1_field_mb_flag, u1_ref_layer_frame_mbs_only_flag);
                }
                else if(1.5 == d_spatial_res_ratio)
                {
                    isvce_ref_array_pos_and_phase_init(
                        ps_layer_props, &ps_mb_states[k + j * u4_wd_in_mbs], &s_mb_pos,
                        u1_frame_mbs_only_flag, u1_field_mb_flag, u1_ref_layer_frame_mbs_only_flag);
                }
            }
        }
    }
}

void isvce_svc_residual_sampling_function_selector(res_pred_state_t *ps_res_pred_state,
                                                   DOUBLE d_spatial_res_ratio, IV_ARCH_T e_arch)
{
    if(2. == d_spatial_res_ratio)
    {
        ps_res_pred_state->apf_residual_samplers[U] = isvce_chroma_residual_sampler_2x;
        ps_res_pred_state->apf_residual_samplers[V] = isvce_chroma_residual_sampler_2x;

        switch(e_arch)
        {
#if defined(X86)
            case ARCH_X86_SSE42:
            {
                ps_res_pred_state->apf_residual_samplers[Y] = isvce_luma_residual_sampler_2x_sse42;

                break;
            }
#elif defined(ARMV8)
            case ARCH_ARM_A53:
            case ARCH_ARM_A57:
            case ARCH_ARM_V8_NEON:
            {
                ps_res_pred_state->apf_residual_samplers[Y] = isvce_luma_residual_sampler_2x_neon;

                break;
            }
#elif !defined(DISABLE_NEON)
            case ARCH_ARM_A9Q:
            case ARCH_ARM_A9A:
            case ARCH_ARM_A9:
            case ARCH_ARM_A7:
            case ARCH_ARM_A5:
            case ARCH_ARM_A15:
            {
                ps_res_pred_state->apf_residual_samplers[Y] = isvce_luma_residual_sampler_2x_neon;

                break;
            }
#endif
            default:
            {
                ps_res_pred_state->apf_residual_samplers[Y] = isvce_luma_residual_sampler_2x;

                break;
            }
        }
    }

    switch(e_arch)
    {
#if defined(X86)
        case ARCH_X86_SSE42:
        {
            ps_res_pred_state->pf_get_sad_with_residual_pred =
                isvce_get_sad_with_residual_pred_sse42;

            break;
    }
#elif defined(ARMV8)
        case ARCH_ARM_A53:
        case ARCH_ARM_A57:
        case ARCH_ARM_V8_NEON:
        {
            ps_res_pred_state->pf_get_sad_with_residual_pred =
                isvce_get_sad_with_residual_pred_neon;

            break;
    }
#elif !defined(DISABLE_NEON)
        case ARCH_ARM_A9Q:
        case ARCH_ARM_A9A:
        case ARCH_ARM_A9:
        case ARCH_ARM_A7:
        case ARCH_ARM_A5:
        case ARCH_ARM_A15:
        {
            ps_res_pred_state->pf_get_sad_with_residual_pred =
                isvce_get_sad_with_residual_pred_neon;

            break;
    }
#endif
    default:
    {
            ps_res_pred_state->pf_get_sad_with_residual_pred = isvce_get_sad_with_residual_pred;

            break;
    }
    }
}

/**
*******************************************************************************
*
* @brief
*  Function to initialize svc ilp buffers
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @param[in] ps_mem_rec
*  Pointer to memory allocated for input buffers
*
*******************************************************************************
*/
void isvce_svc_res_pred_ctxt_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec)
{
    WORD32 i, j, k;

    const WORD32 i4_num_proc_ctxts = sizeof(ps_codec->as_process) / sizeof(ps_codec->as_process[0]);
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;
    DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;
    UWORD32 u4_wd = ps_codec->s_cfg.u4_wd;
    UWORD32 u4_ht = ps_codec->s_cfg.u4_ht;
    UWORD8 *pu1_buf = ps_mem_rec->pv_base;
    WORD64 i8_alloc_mem_size =
        isvce_get_svc_res_pred_ctxt_size(u1_num_spatial_layers, d_spatial_res_ratio, u4_wd, u4_ht);

    if(u1_num_spatial_layers > 1)
    {
        res_pred_mb_state_t *aps_luma_mb_states[MAX_NUM_SPATIAL_LAYERS];
        res_pred_mb_state_t *aps_chroma_mb_states[MAX_NUM_SPATIAL_LAYERS];

        WORD8 *api1_mb_mode[MAX_NUM_SPATIAL_LAYERS];
        WORD32 ai4_mb_mode_stride[MAX_NUM_SPATIAL_LAYERS];

        WORD32 i4_size;

        for(i = 0; i < i4_num_proc_ctxts; i++)
        {
            res_pred_state_t *ps_res_pred_state;
            svc_res_pred_ctxt_t *ps_res_pred_ctxt;
            yuv_buf_props_t *ps_mb_res_buf;
            res_pred_mem_store_t *ps_mem_store;

            isvce_process_ctxt_t *ps_proc = ps_codec->as_process + i;

            ps_res_pred_ctxt = ps_proc->ps_res_pred_ctxt = (svc_res_pred_ctxt_t *) pu1_buf;
            pu1_buf += sizeof(svc_res_pred_ctxt_t);
            i8_alloc_mem_size -= sizeof(svc_res_pred_ctxt_t);

            ps_res_pred_ctxt->s_res_pred_constants.pv_state = pu1_buf;
            ps_res_pred_state = (res_pred_state_t *) pu1_buf;
            pu1_buf += sizeof(res_pred_state_t);
            i8_alloc_mem_size -= sizeof(res_pred_state_t);

            ps_res_pred_state->ps_layer_state = (res_pred_layer_state_t *) pu1_buf;
            pu1_buf += u1_num_spatial_layers * sizeof(ps_res_pred_state->ps_layer_state[0]);
            i8_alloc_mem_size -=
                u1_num_spatial_layers * sizeof(ps_res_pred_state->ps_layer_state[0]);

            i4_size = REF_ARRAY_MAX_WIDTH * REF_ARRAY_MAX_HEIGHT * sizeof(UWORD8);
            ps_res_pred_state->pu1_ref_x_ptr_incr = (UWORD8 *) pu1_buf;
            pu1_buf += i4_size;
            ps_res_pred_state->pu1_ref_y_ptr_incr = (UWORD8 *) pu1_buf;
            pu1_buf += i4_size;

            ASSERT(i8_alloc_mem_size >= 0);

            if(0 == i)
            {
                UWORD32 au4_ref_pos_array_size[NUM_SP_COMPONENTS];
                UWORD32 au4_ref_phase_array_size[NUM_SP_COMPONENTS];

                if(1.5 == d_spatial_res_ratio)
                {
                    au4_ref_pos_array_size[Y] = MB_SIZE * MB_SIZE;
                    au4_ref_phase_array_size[Y] = 5;
                    au4_ref_pos_array_size[U] = (MB_SIZE / 2) * (MB_SIZE / 2);
                    au4_ref_phase_array_size[U] = 5;
                }
                else
                {
                    au4_ref_pos_array_size[Y] = au4_ref_pos_array_size[U] = 0;
                    au4_ref_phase_array_size[Y] = 0;
                    au4_ref_phase_array_size[U] = 3;
                }

                for(j = u1_num_spatial_layers - 1; j >= 1; j--)
                {
                    res_pred_layer_state_t *ps_layer = &ps_res_pred_state->ps_layer_state[j];

                    WORD32 i4_layer_luma_wd =
                        ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                        0.99;
                    WORD32 i4_layer_luma_ht =
                        ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                        0.99;
                    WORD32 i4_layer_luma_mbs =
                        (i4_layer_luma_wd / MB_SIZE) * (i4_layer_luma_ht / MB_SIZE);
                    WORD32 i4_layer_u_wd = i4_layer_luma_wd / 2.0 + 0.99;
                    WORD32 i4_layer_u_ht = i4_layer_luma_ht / 2.0 + 0.99;
                    WORD32 i4_layer_u_mbs =
                        (i4_layer_u_wd / (MB_SIZE / 2)) * (i4_layer_u_ht / (MB_SIZE / 2));

                    ps_layer->ps_luma_mb_states = (res_pred_mb_state_t *) pu1_buf;
                    aps_luma_mb_states[j] = ps_layer->ps_luma_mb_states;
                    pu1_buf += i4_layer_luma_mbs * sizeof(ps_layer->ps_luma_mb_states[0]);
                    i8_alloc_mem_size -=
                        u1_num_spatial_layers * sizeof(ps_layer->ps_luma_mb_states[0]);

                    ps_layer->ps_chroma_mb_states = (res_pred_mb_state_t *) pu1_buf;
                    aps_chroma_mb_states[j] = ps_layer->ps_chroma_mb_states;
                    pu1_buf += i4_layer_u_mbs * sizeof(ps_layer->ps_chroma_mb_states[0]);
                    i8_alloc_mem_size -= i4_layer_u_mbs * sizeof(ps_layer->ps_chroma_mb_states[0]);

                    if(1.5 == d_spatial_res_ratio)
                    {
                        coordinates_t *ps_ref_array_pos = (coordinates_t *) pu1_buf;
                        coordinates_t *ps_ref_array_phases =
                            ps_ref_array_pos + i4_layer_luma_mbs * au4_ref_pos_array_size[Y];

                        for(k = 0; k < i4_layer_luma_mbs; k++)
                        {
                            ps_layer->ps_luma_mb_states[k].ps_ref_array_positions =
                                ps_ref_array_pos + k * au4_ref_pos_array_size[Y];
                            ps_layer->ps_luma_mb_states[k].ps_ref_array_phases =
                                ps_ref_array_phases + k * au4_ref_phase_array_size[Y];
                            pu1_buf += au4_ref_pos_array_size[Y] * sizeof(ps_ref_array_pos[0]);
                            i8_alloc_mem_size -=
                                au4_ref_pos_array_size[Y] * sizeof(ps_ref_array_pos[0]);
                            pu1_buf += au4_ref_phase_array_size[Y] * sizeof(ps_ref_array_phases[0]);
                            i8_alloc_mem_size -=
                                au4_ref_phase_array_size[Y] * sizeof(ps_ref_array_phases[0]);
                        }

                        ps_ref_array_pos = (coordinates_t *) pu1_buf;
                        ps_ref_array_phases =
                            ps_ref_array_pos + i4_layer_u_mbs * au4_ref_pos_array_size[U];

                        for(k = 0; k < i4_layer_u_mbs; k++)
                        {
                            ps_layer->ps_chroma_mb_states[k].ps_ref_array_positions =
                                ps_ref_array_pos + k * au4_ref_pos_array_size[U];
                            ps_layer->ps_chroma_mb_states[k].ps_ref_array_phases =
                                ps_ref_array_phases + k * au4_ref_phase_array_size[U];
                            pu1_buf += au4_ref_pos_array_size[U] * sizeof(ps_ref_array_pos[0]);
                            i8_alloc_mem_size -=
                                au4_ref_pos_array_size[U] * sizeof(ps_ref_array_pos[0]);
                            pu1_buf += au4_ref_phase_array_size[U] * sizeof(ps_ref_array_phases[0]);
                            i8_alloc_mem_size -=
                                au4_ref_phase_array_size[U] * sizeof(ps_ref_array_phases[0]);
                        }
                    }
                    else
                    {
                        coordinates_t *ps_ref_array_pos = NULL;
                        coordinates_t *ps_ref_array_phases = NULL;

                        for(k = 0; k < i4_layer_luma_mbs; k++)
                        {
                            ps_layer->ps_luma_mb_states[k].ps_ref_array_positions =
                                ps_ref_array_pos;
                            ps_layer->ps_luma_mb_states[k].ps_ref_array_phases =
                                ps_ref_array_phases;
                        }

                        ps_ref_array_pos = NULL;
                        ps_ref_array_phases = (coordinates_t *) pu1_buf;

                        for(k = 0; k < i4_layer_u_mbs; k++)
                        {
                            ps_layer->ps_chroma_mb_states[k].ps_ref_array_positions =
                                ps_ref_array_pos;
                            ps_layer->ps_chroma_mb_states[k].ps_ref_array_phases =
                                ps_ref_array_phases;
                        }

                        pu1_buf += au4_ref_phase_array_size[U] * sizeof(ps_ref_array_pos[0]);
                        i8_alloc_mem_size -=
                            au4_ref_phase_array_size[U] * sizeof(ps_ref_array_phases[0]);
                    }

                    ASSERT(i8_alloc_mem_size >= 0);
                    /* Asserts below verify that
                     * 'ps_codec->s_svc_ilp_data.aps_layer_resampler_props' is initialised
                     */
                    ASSERT(ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][j].u4_mb_wd ==
                           MB_SIZE);
                    ASSERT(ps_codec->s_svc_ilp_data.aps_layer_resampler_props[UV][j].u4_mb_wd ==
                           (MB_SIZE / 2));

                    ps_layer->ps_luma_props =
                        &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][j];
                    ps_layer->ps_chroma_props =
                        &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[UV][j];

                    isvce_res_pred_layer_state_init(ps_layer, d_spatial_res_ratio, i4_layer_luma_wd,
                                                    i4_layer_luma_ht,
                                                    ps_codec->s_cfg.e_inp_color_fmt);
                }

                for(j = u1_num_spatial_layers - 1; j >= 0; j--)
                {
                    res_pred_layer_state_t *ps_layer = &ps_res_pred_state->ps_layer_state[j];

                    WORD32 i4_layer_luma_wd =
                        ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                        0.99;
                    WORD32 i4_layer_luma_ht =
                        ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                        0.99;
                    WORD32 i4_layer_luma_mbs =
                        ((i4_layer_luma_wd / MB_SIZE) + 2) * ((i4_layer_luma_ht / MB_SIZE) + 2);

                    ps_layer->pi1_mb_mode = (WORD8 *) pu1_buf;
                    pu1_buf += i4_layer_luma_mbs * sizeof(WORD8);
                    memset(ps_layer->pi1_mb_mode, -1, i4_layer_luma_mbs);

                    ps_layer->i4_mb_mode_stride = ai4_mb_mode_stride[j] =
                        (i4_layer_luma_wd / MB_SIZE) + 2;
                    ps_layer->pi1_mb_mode += 1 + ps_layer->i4_mb_mode_stride;
                    api1_mb_mode[j] = ps_layer->pi1_mb_mode;
                }
            }
            else
            {
                for(j = u1_num_spatial_layers - 1; j >= 1; j--)
                {
                    res_pred_layer_state_t *ps_layer = &ps_res_pred_state->ps_layer_state[j];

                    ps_layer->ps_luma_mb_states = aps_luma_mb_states[j];
                    ps_layer->ps_chroma_mb_states = aps_chroma_mb_states[j];

                    ps_layer->ps_luma_props =
                        &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[Y][j];
                    ps_layer->ps_chroma_props =
                        &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[UV][j];
                }
                for(j = u1_num_spatial_layers - 1; j >= 0; j--)
                {
                    res_pred_layer_state_t *ps_layer = &ps_res_pred_state->ps_layer_state[j];

                    ps_layer->pi1_mb_mode = api1_mb_mode[j];
                    ps_layer->i4_mb_mode_stride = ai4_mb_mode_stride[j];
                }
            }

            ps_mb_res_buf = &ps_res_pred_ctxt->s_res_pred_outputs.s_res_pred;
            ps_mem_store = &ps_res_pred_state->s_mem_store;
            ps_proc->ps_mb_res_buf = ps_mb_res_buf;

            for(j = 0; j < NUM_SP_COMPONENTS; j++)
            {
                buffer_container_t *ps_comp_buf = &ps_mb_res_buf->as_component_bufs[j];

                UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) j));

                ps_comp_buf->pv_data = pu1_buf;
                ps_comp_buf->i4_data_stride = MB_SIZE;
                pu1_buf += MB_SIZE * (MB_SIZE >> u1_is_chroma) * sizeof(WORD16);
                i8_alloc_mem_size -= MB_SIZE * (MB_SIZE >> u1_is_chroma) * sizeof(WORD16);
            }

            ps_mem_store->s_scratch.pv_data = pu1_buf;
            ps_mem_store->s_scratch.i4_data_stride = REF_ARRAY_MAX_WIDTH;
            pu1_buf += REF_ARRAY_MAX_WIDTH * REF_ARRAY_MAX_HEIGHT * sizeof(WORD16);
            i8_alloc_mem_size -= REF_ARRAY_MAX_WIDTH * REF_ARRAY_MAX_HEIGHT * sizeof(WORD16);

            ASSERT(i8_alloc_mem_size >= 0);

            ps_mb_res_buf->as_component_bufs[V].pv_data = NULL;
            ps_mb_res_buf->e_color_format = IV_YUV_420SP_UV;
            ps_mb_res_buf->u1_bit_depth = 10;
            ps_mb_res_buf->u4_width = MB_SIZE;
            ps_mb_res_buf->u4_height = MB_SIZE;

            isvce_svc_residual_sampling_function_selector(ps_res_pred_state, d_spatial_res_ratio,
                                                          ps_codec->s_cfg.e_arch);
        }
    }
    else
    {
        for(i = 0; i < i4_num_proc_ctxts; i++)
        {
            isvce_process_ctxt_t *ps_proc = ps_codec->as_process + i;

            ps_proc->ps_res_pred_ctxt = NULL;

            ps_proc->ps_mb_res_buf = (yuv_buf_props_t *) pu1_buf;
            pu1_buf += sizeof(yuv_buf_props_t);
            i8_alloc_mem_size -= sizeof(yuv_buf_props_t);

            for(j = 0; j < NUM_SP_COMPONENTS; j++)
            {
                buffer_container_t *ps_comp_buf = &ps_proc->ps_mb_res_buf->as_component_bufs[j];

                UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) j));

                ps_comp_buf->pv_data = pu1_buf;
                ps_comp_buf->i4_data_stride = MB_SIZE;
                pu1_buf += MB_SIZE * (MB_SIZE >> u1_is_chroma) * sizeof(WORD16);
                i8_alloc_mem_size -= MB_SIZE * (MB_SIZE >> u1_is_chroma) * sizeof(WORD16);
            }

            ASSERT(i8_alloc_mem_size >= 0);
        }
    }
}

void isvce_get_mb_residual_pred(svc_res_pred_ctxt_t *ps_res_pred_ctxt)
{
    buffer_container_t s_inp;
    buffer_container_t s_out;
    coordinates_t s_frame_dims;
    coordinates_t s_frame_dims_in_mbs;
    coordinates_t s_ref_array_offsets;
    svc_layer_data_t *ps_ref_layer_data;
    res_pred_layer_state_t *ps_layer_state;
    yuv_buf_props_t *ps_ref_residual_buf;
    res_pred_mb_state_t *ps_luma_mb_state;
    res_pred_mb_state_t *ps_chroma_mb_state;
    isvce_mb_info_t *ps_ref_mb;

    WORD32 i;

    res_pred_constants_t *ps_res_pred_constants = &ps_res_pred_ctxt->s_res_pred_constants;
    res_pred_variables_t *ps_res_pred_variables = &ps_res_pred_ctxt->s_res_pred_variables;
    res_pred_outputs_t *ps_res_pred_outputs = &ps_res_pred_ctxt->s_res_pred_outputs;
    res_pred_state_t *ps_res_pred_state = (res_pred_state_t *) ps_res_pred_constants->pv_state;
    res_pred_mem_store_t *ps_mem_store = &ps_res_pred_state->s_mem_store;
    svc_ilp_data_t *ps_svc_ilp_data = ps_res_pred_variables->ps_svc_ilp_data;
    coordinates_t *ps_mb_pos = &ps_res_pred_variables->s_mb_pos;

    UWORD8 u1_spatial_layer_id = ps_res_pred_variables->u1_spatial_layer_id;

    ASSERT(u1_spatial_layer_id > 0);

    s_frame_dims.i4_abscissa = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_width;
    s_frame_dims.i4_ordinate = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_height;
    s_frame_dims_in_mbs.i4_abscissa = s_frame_dims.i4_abscissa / MB_SIZE;
    s_frame_dims_in_mbs.i4_ordinate = s_frame_dims.i4_ordinate / MB_SIZE;

    ps_ref_layer_data =
        &ps_svc_ilp_data->ps_svc_au_data->ps_svc_layer_data[u1_spatial_layer_id - 1];
    ps_layer_state = &ps_res_pred_state->ps_layer_state[u1_spatial_layer_id];
    ps_ref_residual_buf = &ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id - 1];
    ps_luma_mb_state = ps_layer_state->ps_luma_mb_states + ps_mb_pos->i4_abscissa +
                       ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa;
    ps_chroma_mb_state = ps_layer_state->ps_chroma_mb_states + ps_mb_pos->i4_abscissa +
                         ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa;

    for(i = 0; i < NUM_COMPONENTS; i++)
    {
        res_pred_mb_state_t *ps_mb_state;
        layer_resampler_props_t *ps_layer_props;

        UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) i));

        ps_mb_state = u1_is_chroma ? ps_chroma_mb_state : ps_luma_mb_state;
        ps_layer_props =
            u1_is_chroma ? ps_layer_state->ps_chroma_props : ps_layer_state->ps_luma_props;

        /* Presence of appropriate padding is assumed */
        s_ref_array_offsets = ps_mb_state->s_offsets;

        s_inp = ps_ref_residual_buf->as_component_bufs[u1_is_chroma ? UV : Y];
        s_inp.pv_data = ((WORD16 *) s_inp.pv_data) + (V == ((COMPONENT_TYPE) i)) +
                        (s_ref_array_offsets.i4_abscissa << u1_is_chroma) +
                        s_ref_array_offsets.i4_ordinate * s_inp.i4_data_stride;

        s_out = ps_res_pred_outputs->s_res_pred.as_component_bufs[u1_is_chroma ? UV : Y];
        s_out.pv_data = ((WORD16 *) s_out.pv_data) + (V == ((COMPONENT_TYPE) i));

        ps_ref_mb =
            ps_ref_layer_data->ps_mb_info +
            ((s_ref_array_offsets.i4_abscissa + (ps_mb_state->s_ref_array_dims.i4_abscissa / 2)) /
             ps_layer_props->u4_mb_wd) +
            ((s_ref_array_offsets.i4_ordinate + (ps_mb_state->s_ref_array_dims.i4_ordinate / 2)) /
             ps_layer_props->u4_mb_ht) *
                (s_frame_dims_in_mbs.i4_abscissa / 2);

        ps_res_pred_state->apf_residual_samplers[i](
            ps_mb_state->ps_ref_array_positions, ps_mb_state->ps_ref_array_phases, &s_inp, &s_out,
            &ps_mem_store->s_scratch, UINT32_MAX, ps_ref_mb->u1_tx_size == 8);
    }
}

void isvce_get_ref_layer_mbtype_tx_size(WORD8 *pi1_ref_mb_modes, WORD32 i4_ref_mode_stride,
                                        WORD32 i4_element_size, WORD32 i4_x_ref, WORD32 i4_y_ref,
                                        WORD32 *pi4_mb_type, WORD32 *pi4_tx_size,
                                        WORD32 i4_chroma_flag)
{
    WORD32 i4_mb_wd_sft, i4_mb_ht_sft;
    WORD32 i4_mb_x, i4_mb_y;
    WORD8 i1_mb_mode;

    if(i4_x_ref < 0)
    {
        i4_x_ref = 0;
    }
    if(i4_y_ref < 0)
    {
        i4_y_ref = 0;
    }

    i4_mb_wd_sft = (MB_WIDTH_SHIFT - i4_chroma_flag);
    i4_mb_ht_sft = (MB_HEIGHT_SHIFT - i4_chroma_flag);
    i4_mb_x = (i4_x_ref >> i4_mb_wd_sft);
    i4_mb_y = (i4_y_ref >> i4_mb_ht_sft);

    pi1_ref_mb_modes += (i4_mb_y * i4_ref_mode_stride * i4_element_size);
    pi1_ref_mb_modes += (i4_mb_x * i4_element_size);
    i1_mb_mode = *pi1_ref_mb_modes;
    i1_mb_mode = (i1_mb_mode < 0) ? i1_mb_mode : SVC_EXTRACT_MB_MODE(*pi1_ref_mb_modes);

    if(i1_mb_mode <= SVC_INTER_MB)
    {
        *pi4_mb_type = SVC_INTER_MB;
        *pi4_tx_size = GET_BIT_TX_SIZE(*pi1_ref_mb_modes, 1);
    }
    else
    {
        *pi4_mb_type = SVC_INTRA_MB;
        *pi4_tx_size = 1;
    }
}

void isvce_ref_layer_ptr_incr(WORD8 *pi1_ref_mb_modes, WORD32 i4_ref_mode_stride,
                              WORD32 i4_element_size, WORD32 i4_x_offset, WORD32 i4_y_offset,
                              WORD32 i4_refary_wd, WORD32 i4_refary_ht, UWORD8 *pu1_ref_x_ptr_incr,
                              UWORD8 *pu1_ref_y_ptr_incr, WORD32 i4_chroma_flag)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_x_idx, i4_y_idx;
    WORD32 i4_prev_x, i4_prev_y;
    WORD32 i4_const_val;
    WORD32 i4_pos_x, i4_pos_y;
    WORD32 i4_trans_size;
    WORD32 i4_mb_type, i4_tx_size;
    WORD32 i4_act_ary_wd, i4_act_ary_ht;
    WORD32 i4_and_const;
    UWORD8 *pu1_incr_x, *pu1_incr_y;

    memset(pu1_ref_x_ptr_incr, 1, (i4_refary_wd * i4_refary_ht));
    memset(pu1_ref_y_ptr_incr, 1, (i4_refary_wd * i4_refary_ht));

    i4_act_ary_wd = i4_refary_wd;
    i4_act_ary_ht = i4_refary_ht;

    i4_x = 0;
    i4_y = 0;
    i4_prev_y = 0;

    if(0 == i4_chroma_flag)
    {
        do
        {
            WORD32 i4_x_ref, i4_y_ref;
            WORD32 i4_idx;
            WORD32 i4_wd, i4_ht;
            WORD32 i4_max_pos_x, i4_max_pos_y;

            i4_prev_x = i4_x;

            i4_x_ref = i4_x_offset + i4_x;
            i4_y_ref = i4_y_offset + i4_y;

            isvce_get_ref_layer_mbtype_tx_size(pi1_ref_mb_modes, i4_ref_mode_stride,
                                               i4_element_size, i4_x_ref, i4_y_ref, &i4_mb_type,
                                               &i4_tx_size, i4_chroma_flag);

            i4_trans_size = ((i4_tx_size + 1) << 2);
            i4_const_val = i4_trans_size - 1;
            i4_and_const = i4_const_val;

            /* Fill horizontal tx block edges of current reference mb with 0 */
            pu1_incr_x = pu1_ref_x_ptr_incr + i4_x;
            pu1_incr_x += (i4_y * i4_refary_wd);

            i4_ht = (16 - (i4_y_ref & 0xF));
            i4_ht = MIN((i4_act_ary_ht - i4_y), i4_ht);

            i4_x_idx = i4_x;

            i4_pos_x = i4_x_ref & 0xF;

            i4_max_pos_x = 16;
            i4_x += (16 - i4_pos_x);

            /* Get the transform block edge pos */
            i4_idx = (i4_const_val - (i4_pos_x & i4_and_const));

            i4_x_idx += i4_idx;

            while((i4_pos_x < i4_max_pos_x) && (i4_x_idx < i4_act_ary_wd))
            {
                WORD32 i4_i;
                UWORD8 *pu1_incr;

                pu1_incr = pu1_incr_x + i4_idx;

                for(i4_i = 0; i4_i < i4_ht; i4_i++)
                { /* Fill the block edge with 0s */
                    *pu1_incr = 0;
                    pu1_incr += i4_refary_wd;
                }

                i4_pos_x += i4_trans_size;
                pu1_incr_x += i4_trans_size;
                i4_x_idx += MIN(i4_trans_size, (i4_act_ary_wd - i4_x_idx));
            }

            /* Fill vertical tx block edges of current reference mb with 0 */
            pu1_incr_y = pu1_ref_y_ptr_incr + i4_prev_x;
            pu1_incr_y += (i4_y * i4_refary_wd);

            i4_wd = (16 - (i4_x_ref & 0xF));
            i4_wd = MIN((i4_act_ary_wd - i4_prev_x), i4_wd);

            i4_y_idx = i4_y;

            i4_pos_y = i4_y_ref & 0xF;

            i4_max_pos_y = 16;
            i4_y += (16 - i4_pos_y);

            /* Get the transform block edge pos */
            i4_idx = (i4_const_val - (i4_pos_y & i4_and_const));

            i4_y_idx += i4_idx;

            while((i4_pos_y < i4_max_pos_y) && (i4_y_idx < i4_act_ary_ht))
            {
                WORD32 i4_i;
                UWORD8 *pu1_incr;

                pu1_incr = pu1_incr_y + i4_idx * i4_refary_wd;

                for(i4_i = 0; i4_i < i4_wd; i4_i++)
                { /* Fill the block edge with 0s */
                    *pu1_incr = 0;
                    pu1_incr++;
                }

                i4_pos_y += i4_trans_size;
                pu1_incr_y += i4_trans_size * i4_refary_wd;
                i4_y_idx += MIN(i4_trans_size, (i4_act_ary_ht - i4_y_idx));
            }

            if(i4_x < i4_act_ary_wd)
            {
                i4_y = i4_prev_y;
            }
            else if(i4_y < i4_act_ary_ht)
            {
                i4_prev_y = i4_y;
                i4_x = 0;
            }
        } while((i4_y < i4_act_ary_ht) || (i4_x < i4_act_ary_wd));
    }
    else
    {
        i4_trans_size = 4;
        i4_const_val = 3;

        do
        {
            WORD32 i4_x_ref, i4_y_ref;
            WORD32 i4_idx;
            WORD32 i4_wd, i4_ht;
            WORD32 i4_max_pos_x, i4_max_pos_y;

            i4_prev_x = i4_x;

            i4_x_ref = i4_x_offset + i4_x;
            i4_y_ref = i4_y_offset + i4_y;

            /* Fill horizontal tx block edges of current reference mb with 0 */
            pu1_incr_x = pu1_ref_x_ptr_incr + i4_x;
            pu1_incr_x += (i4_y * i4_refary_wd);

            i4_ht = (8 - (i4_y_ref & 0x7));
            i4_ht = MIN((i4_act_ary_ht - i4_y), i4_ht);

            i4_x_idx = i4_x;

            i4_pos_x = i4_x_ref & 0x7;

            i4_max_pos_x = 8;
            i4_x += (8 - i4_pos_x);

            /* Get the transform block edge pos */
            i4_idx = (i4_const_val - (i4_pos_x & 0x3));

            i4_x_idx += i4_idx;

            while((i4_pos_x < i4_max_pos_x) && (i4_x_idx < i4_act_ary_wd))
            {
                WORD32 i4_i;
                UWORD8 *pu1_incr;

                pu1_incr = pu1_incr_x + i4_idx;

                for(i4_i = 0; i4_i < i4_ht; i4_i++)
                { /* Fill the block edge with 0s */
                    *pu1_incr = 0;
                    pu1_incr += i4_refary_wd;
                }

                i4_pos_x += i4_trans_size;
                pu1_incr_x += i4_trans_size;
                i4_x_idx += MIN(i4_trans_size, (i4_act_ary_wd - i4_x_idx));
            }

            /* Fill vertical tx block edges of current reference mb with 0 */
            pu1_incr_y = pu1_ref_y_ptr_incr + i4_prev_x;
            pu1_incr_y += (i4_y * i4_refary_wd);

            i4_wd = (8 - (i4_x_ref & 0x7));
            i4_wd = MIN((i4_act_ary_wd - i4_prev_x), i4_wd);

            i4_y_idx = i4_y;

            i4_pos_y = i4_y_ref & 0x7;

            i4_max_pos_y = 8;
            i4_y += (8 - i4_pos_y);

            /* Get the transform block edge pos */
            i4_idx = (i4_const_val - (i4_pos_y & 0x3));

            i4_y_idx += i4_idx;

            while((i4_pos_y < i4_max_pos_y) && (i4_y_idx < i4_act_ary_ht))
            {
                WORD32 i4_i;
                UWORD8 *pu1_incr;

                pu1_incr = pu1_incr_y + i4_idx * i4_refary_wd;

                for(i4_i = 0; i4_i < i4_wd; i4_i++)
                { /* Fill the block edge with 0s */
                    *pu1_incr = 0;
                    pu1_incr++;
                }

                i4_pos_y += i4_trans_size;
                pu1_incr_y += i4_trans_size * i4_refary_wd;
                i4_y_idx += MIN(i4_trans_size, (i4_act_ary_ht - i4_y_idx));
            }

            if(i4_x < i4_act_ary_wd)
            {
                i4_y = i4_prev_y;
            }
            else if(i4_y < i4_act_ary_ht)
            {
                i4_prev_y = i4_y;
                i4_x = 0;
            }
        } while((i4_y < i4_act_ary_ht) || (i4_x < i4_act_ary_wd));
    }
}

void isvce_residual_reflayer_const(svc_res_pred_ctxt_t *ps_res_pred_ctxt, WORD16 *pi2_inp_data,
                                   WORD32 i4_inp_data_stride, WORD8 *ps_ref_mb_mode,
                                   WORD32 i4_ref_mb_mode_stride, WORD32 *pi4_refarr_wd,
                                   WORD32 i4_chroma_flag)
{
    WORD8 *pi1_ref_mb_modes;
    WORD32 i4_ref_mode_stride;

    WORD32 i4_x, i4_y;
    WORD32 i4_ref_wd;
    WORD32 i4_ref_ht;
    WORD32 i4_x_offset;
    WORD32 i4_y_offset;
    WORD32 i4_refarray_wd;
    WORD32 i4_refarray_ht;

    WORD16 *pi2_ref_array;

    res_pred_mb_state_t *ps_mb_states;
    res_pred_layer_state_t *ps_layer_state;

    res_pred_constants_t *ps_res_pred_constants = &ps_res_pred_ctxt->s_res_pred_constants;
    res_pred_variables_t *ps_res_pred_variables = &ps_res_pred_ctxt->s_res_pred_variables;
    res_pred_state_t *ps_res_pred_state = (res_pred_state_t *) ps_res_pred_constants->pv_state;
    res_pred_mem_store_t *ps_mem_store = &ps_res_pred_state->s_mem_store;
    svc_ilp_data_t *ps_svc_ilp_data = ps_res_pred_variables->ps_svc_ilp_data;
    coordinates_t *ps_mb_pos = &ps_res_pred_variables->s_mb_pos;

    UWORD8 u1_spatial_layer_id = ps_res_pred_variables->u1_spatial_layer_id;

    ps_layer_state = &ps_res_pred_state->ps_layer_state[u1_spatial_layer_id];
    pi2_ref_array = (WORD16 *) ps_mem_store->s_scratch.pv_data;

    pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode;
    i4_ref_mode_stride = i4_ref_mb_mode_stride;

    ASSERT(NULL != pi1_ref_mb_modes);

    {
        WORD32 i4_base_width;
        WORD32 i4_base_height;

        coordinates_t s_frame_dims, s_frame_dims_in_mbs;

        s_frame_dims.i4_abscissa = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_width;
        s_frame_dims.i4_ordinate = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_height;
        s_frame_dims_in_mbs.i4_abscissa = s_frame_dims.i4_abscissa / MB_SIZE;
        s_frame_dims_in_mbs.i4_ordinate = s_frame_dims.i4_ordinate / MB_SIZE;

        ps_mb_states = i4_chroma_flag ? ps_layer_state->ps_chroma_mb_states
                                      : ps_layer_state->ps_luma_mb_states;

        ps_mb_states +=
            ps_mb_pos->i4_abscissa + ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa;

        i4_base_width = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id - 1].u4_width;
        i4_base_height = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id - 1].u4_height;

        i4_ref_wd = i4_base_width >> i4_chroma_flag;
        i4_ref_ht = i4_base_height >> i4_chroma_flag;

        i4_x_offset = ps_mb_states->s_offsets.i4_abscissa;
        i4_y_offset = ps_mb_states->s_offsets.i4_ordinate;
        i4_refarray_wd = ps_mb_states->s_ref_array_dims.i4_abscissa;
        i4_refarray_ht = ps_mb_states->s_ref_array_dims.i4_ordinate;
    }

    {
        isvce_ref_layer_ptr_incr(pi1_ref_mb_modes, i4_ref_mode_stride, 1, i4_x_offset, i4_y_offset,
                                 i4_refarray_wd, i4_refarray_ht,
                                 ps_res_pred_state->pu1_ref_x_ptr_incr,
                                 ps_res_pred_state->pu1_ref_y_ptr_incr, i4_chroma_flag);
    }

    for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
    {
        for(i4_x = 0; i4_x < i4_refarray_wd; i4_x++)
        {
            WORD32 i4_x_ref;
            WORD32 i4_y_ref;
            WORD32 i4_ref_mb_type, i4_ref_tx_size;
            WORD16 *pi2_ref_data_byte;
            WORD16 *pi2_ref_array_temp;

            i4_x_ref = MAX(0, MIN(i4_ref_wd - 1, i4_x + i4_x_offset));
            i4_y_ref = MAX(0, MIN(i4_ref_ht - 1, i4_y + i4_y_offset));

            isvce_get_ref_layer_mbtype_tx_size(pi1_ref_mb_modes, i4_ref_mode_stride, 1, i4_x_ref,
                                               i4_y_ref, &i4_ref_mb_type, &i4_ref_tx_size,
                                               i4_chroma_flag);

            if(0 <= i4_x_offset)
            {
                i4_x_ref = i4_x_ref - i4_x_offset;
            }

            if(0 <= i4_y_offset)
            {
                i4_y_ref = i4_y_ref - i4_y_offset;
            }

            pi2_ref_array_temp = pi2_ref_array + i4_x;
            pi2_ref_array_temp += i4_y * i4_refarray_wd;

            if(SVC_INTER_MB == i4_ref_mb_type)
            {
                pi2_ref_data_byte = pi2_inp_data + (i4_x_ref << i4_chroma_flag);
                pi2_ref_data_byte += i4_y_ref * i4_inp_data_stride;

                *pi2_ref_array_temp = (WORD16) (*pi2_ref_data_byte);
            }
            else
            {
                *pi2_ref_array_temp = 0;
            }
        }
    }
    *pi4_refarr_wd = i4_refarray_wd;
}

void isvce_interpolate_residual(svc_res_pred_ctxt_t *ps_res_pred_ctxt, WORD16 *pi2_out,
                                WORD32 i4_out_stride, WORD32 i4_refarray_wd, WORD32 i4_chroma_flag,
                                coordinates_t *ps_mb_pos)
{
    res_pred_constants_t *ps_res_pred_constants = &ps_res_pred_ctxt->s_res_pred_constants;
    res_pred_variables_t *ps_res_pred_variables = &ps_res_pred_ctxt->s_res_pred_variables;
    res_pred_state_t *ps_res_pred_state = (res_pred_state_t *) ps_res_pred_constants->pv_state;
    res_pred_mem_store_t *ps_mem_store = &ps_res_pred_state->s_mem_store;

    WORD32 i4_x, i4_y;
    WORD32 i4_temp_array_ht;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;
    WORD16 *pi2_ref_array;
    UWORD8 *pu1_ref_x_ptr_incr, *pu1_ref_y_ptr_incr;

    coordinates_t *ps_phase;
    coordinates_t *ps_pos;
    res_pred_mb_state_t *ps_mb_states;

    coordinates_t s_frame_dims;
    coordinates_t s_frame_dims_in_mbs;

    UWORD8 u1_spatial_layer_id = ps_res_pred_variables->u1_spatial_layer_id;

    svc_ilp_data_t *ps_svc_ilp_data = ps_res_pred_variables->ps_svc_ilp_data;

    res_pred_mb_state_t *ps_mb_state;

    s_frame_dims.i4_abscissa = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_width;
    s_frame_dims.i4_ordinate = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_height;
    s_frame_dims_in_mbs.i4_abscissa = s_frame_dims.i4_abscissa / MB_SIZE;
    s_frame_dims_in_mbs.i4_ordinate = s_frame_dims.i4_ordinate / MB_SIZE;

    pu1_ref_x_ptr_incr = ps_res_pred_state->pu1_ref_x_ptr_incr;
    pu1_ref_y_ptr_incr = ps_res_pred_state->pu1_ref_y_ptr_incr;

    ps_mb_states = i4_chroma_flag
                       ? ps_res_pred_state->ps_layer_state[u1_spatial_layer_id].ps_chroma_mb_states
                       : ps_res_pred_state->ps_layer_state[u1_spatial_layer_id].ps_luma_mb_states;

    i4_mb_wd = MB_SIZE >> i4_chroma_flag;
    i4_mb_ht = MB_SIZE >> i4_chroma_flag;

    ps_mb_state = &ps_mb_states[ps_mb_pos->i4_abscissa +
                                (ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa)];

    ps_phase = ps_mb_state->ps_ref_array_phases;
    ps_pos = ps_mb_state->ps_ref_array_positions;

    i4_temp_array_ht = i4_mb_ht;

    pi2_ref_array = (WORD16 *) ps_mem_store->s_scratch.pv_data;

    for(i4_y = 0; i4_y < i4_temp_array_ht; i4_y++)
    {
        for(i4_x = 0; i4_x < i4_mb_wd; i4_x++)
        {
            WORD32 i4_i;
            WORD32 i4_y_ref;
            WORD32 i4_y_phase;
            WORD32 i4_x_ref;
            WORD32 i4_x_phase;
            WORD32 i4_x_ref_round;
            WORD16 *pi2_out_curr;
            WORD32 ai4_temp_pred[2];
            UWORD8 *pu1_ref_y_ptr_incr_temp;
            WORD32 *pi4_temp_pred;
            UWORD8 u1_incr_y;
            WORD16 i2_res;

            pi2_out_curr = pi2_out + (i4_x << i4_chroma_flag) + (i4_y * i4_out_stride);

            i4_y_ref = ps_pos[(i4_mb_wd * i4_y) + i4_x].i4_ordinate;
            i4_y_phase = ps_phase[((i4_y % 3) > 0) * 2 + (i4_y % 3)].i4_ordinate;

            i4_x_ref = ps_pos[(i4_mb_wd * i4_y) + i4_x].i4_abscissa;
            i4_x_phase = ps_phase[i4_x % 3].i4_abscissa;

            /* horizontal processing*/
            for(i4_i = 0; i4_i < 2; i4_i++)
            {
                UWORD8 *pu1_ref_x_ptr_incr_temp;
                UWORD8 u1_incr;
                WORD16 *pi2_ref_array_1, *pi2_ref_array_2;

                pu1_ref_x_ptr_incr_temp = pu1_ref_x_ptr_incr + i4_x_ref;
                pu1_ref_x_ptr_incr_temp += ((i4_y_ref + i4_i) * i4_refarray_wd);
                u1_incr = *pu1_ref_x_ptr_incr_temp;

                pi2_ref_array_1 = pi2_ref_array + i4_x_ref;
                pi2_ref_array_1 += ((i4_y_ref + i4_i) * i4_refarray_wd);

                if(!u1_incr)
                {
                    pi2_ref_array_1 += (i4_x_phase >> 3);
                }

                pi2_ref_array_2 = pi2_ref_array_1 + u1_incr;

                ai4_temp_pred[i4_i] =
                    (16 - i4_x_phase) * (*pi2_ref_array_1) + i4_x_phase * (*pi2_ref_array_2);
            }

            /* vertical processing */
            i4_x_ref_round = (i4_x_ref + (i4_x_phase >> 3));

            pu1_ref_y_ptr_incr_temp =
                pu1_ref_y_ptr_incr + i4_x_ref_round + (i4_y_ref * i4_refarray_wd);
            u1_incr_y = *pu1_ref_y_ptr_incr_temp;

            pi4_temp_pred = &ai4_temp_pred[0];
            if(!u1_incr_y)
            {
                pi4_temp_pred += (i4_y_phase >> 3);
            }
            i2_res = (((16 - i4_y_phase) * pi4_temp_pred[0] +
                       i4_y_phase * pi4_temp_pred[u1_incr_y] + 128) >>
                      8);
            *pi2_out_curr = i2_res;
        }
    }
}

void isvce_get_mb_residual_pred_non_dyadic(svc_res_pred_ctxt_t *ps_res_pred_ctxt)
{
    buffer_container_t s_inp;
    buffer_container_t s_out;
    coordinates_t s_frame_dims;
    coordinates_t s_frame_dims_in_mbs;
    coordinates_t s_ref_array_offsets;
    res_pred_layer_state_t *ps_layer_state, *ps_ref_layer_state;
    yuv_buf_props_t *ps_ref_residual_buf;
    res_pred_mb_state_t *ps_luma_mb_state;
    res_pred_mb_state_t *ps_chroma_mb_state;

    WORD16 *pi2_inp, *pi2_out;
    WORD32 i4_inp_stride, i4_out_stride;

    res_pred_constants_t *ps_res_pred_constants = &ps_res_pred_ctxt->s_res_pred_constants;
    res_pred_variables_t *ps_res_pred_variables = &ps_res_pred_ctxt->s_res_pred_variables;
    res_pred_outputs_t *ps_res_pred_outputs = &ps_res_pred_ctxt->s_res_pred_outputs;
    res_pred_state_t *ps_res_pred_state = (res_pred_state_t *) ps_res_pred_constants->pv_state;
    svc_ilp_data_t *ps_svc_ilp_data = ps_res_pred_variables->ps_svc_ilp_data;
    coordinates_t *ps_mb_pos = &ps_res_pred_variables->s_mb_pos;

    UWORD8 u1_spatial_layer_id = ps_res_pred_variables->u1_spatial_layer_id;

    WORD32 i4_refarray_wd;

    WORD32 i;

    ASSERT(u1_spatial_layer_id > 0);

    s_frame_dims.i4_abscissa = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_width;
    s_frame_dims.i4_ordinate = ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id].u4_height;
    s_frame_dims_in_mbs.i4_abscissa = s_frame_dims.i4_abscissa / MB_SIZE;
    s_frame_dims_in_mbs.i4_ordinate = s_frame_dims.i4_ordinate / MB_SIZE;

    ps_layer_state = &ps_res_pred_state->ps_layer_state[u1_spatial_layer_id];
    ps_ref_layer_state = &ps_res_pred_state->ps_layer_state[u1_spatial_layer_id - 1];
    ps_ref_residual_buf = &ps_svc_ilp_data->ps_residual_bufs[u1_spatial_layer_id - 1];
    ps_luma_mb_state = ps_layer_state->ps_luma_mb_states + ps_mb_pos->i4_abscissa +
                       ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa;
    ps_chroma_mb_state = ps_layer_state->ps_chroma_mb_states + ps_mb_pos->i4_abscissa +
                         ps_mb_pos->i4_ordinate * s_frame_dims_in_mbs.i4_abscissa;

    for(i = 0; i < NUM_COMPONENTS; i++)
    {
        res_pred_mb_state_t *ps_mb_state;

        UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) i));

        ps_mb_state = u1_is_chroma ? ps_chroma_mb_state : ps_luma_mb_state;

        s_ref_array_offsets.i4_abscissa =
            MAX(0, MIN(ps_mb_state->s_offsets.i4_abscissa,
                       (s_frame_dims.i4_abscissa >> u1_is_chroma) - 1));
        s_ref_array_offsets.i4_ordinate =
            MAX(0, MIN(ps_mb_state->s_offsets.i4_ordinate,
                       (s_frame_dims.i4_ordinate >> u1_is_chroma) - 1));

        s_inp = ps_ref_residual_buf->as_component_bufs[u1_is_chroma ? UV : Y];
        s_inp.pv_data = ((WORD16 *) s_inp.pv_data) + (V == ((COMPONENT_TYPE) i)) +
                        (s_ref_array_offsets.i4_abscissa << u1_is_chroma) +
                        s_ref_array_offsets.i4_ordinate * s_inp.i4_data_stride;

        s_out = ps_res_pred_outputs->s_res_pred.as_component_bufs[u1_is_chroma ? UV : Y];
        s_out.pv_data = ((WORD16 *) s_out.pv_data) + (V == ((COMPONENT_TYPE) i));

        pi2_inp = (WORD16 *) s_inp.pv_data;
        pi2_out = (WORD16 *) s_out.pv_data;

        i4_inp_stride = s_inp.i4_data_stride;
        i4_out_stride = s_out.i4_data_stride;

        /* ------- Constructing refSampleArray ----------------------- */
        isvce_residual_reflayer_const(
            ps_res_pred_ctxt, pi2_inp, i4_inp_stride, ps_ref_layer_state->pi1_mb_mode,
            ps_ref_layer_state->i4_mb_mode_stride, &i4_refarray_wd, u1_is_chroma);

        /* ---- Interpolation process for Residual prediction	 ------ */
        isvce_interpolate_residual(ps_res_pred_ctxt, pi2_out, i4_out_stride, i4_refarray_wd,
                                   u1_is_chroma, ps_mb_pos);
    }
}

UWORD32 isvce_get_sad_with_residual_pred(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                         buffer_container_t *ps_res, UWORD32 u4_mb_wd,
                                         UWORD32 u4_mb_ht)
{
    UWORD32 i, j;

    UWORD32 u4_sad = 0;

    for(i = 0; i < u4_mb_ht; i++)
    {
        for(j = 0; j < u4_mb_wd; j++)
        {
            WORD16 i2_src = ((UWORD8 *) ps_src->pv_data)[j + i * ps_src->i4_data_stride];
            WORD16 i2_pred = ((UWORD8 *) ps_pred->pv_data)[j + i * ps_pred->i4_data_stride];
            WORD16 i2_res = ((WORD16 *) ps_res->pv_data)[j + i * ps_res->i4_data_stride];

            u4_sad += ABS(i2_src - i2_pred - i2_res);
        }
    }
    return u4_sad;
}

/**
*******************************************************************************
*
* @brief
*  Function to evaluate residual_prediction_flag
*
* @param[in] ps_src
*  Pointer to MB src buffers
*
* @param[in] ps_pred
*  Pointer to MB pred buffers
*
* @param[in] ps_res
*  Pointer to MB res buffers
*
* @param[out] pu4_res_pred_sad
*  Output variable for SAD
*
* @param[out] pu1_residual_prediction_flag
*  Output variable for residual_prediction_flag
*
* @param[in] u4_winning_sad
*  Winning mode's SAD
*
* @notes The algorithm currently uses only luma for evaluating
*        residual_prediction_flag.
*
*******************************************************************************
*/
void isvce_residual_pred_eval(svc_res_pred_ctxt_t *ps_res_pred_ctxt, yuv_buf_props_t *ps_src,
                              yuv_buf_props_t *ps_pred, yuv_buf_props_t *ps_res,
                              UWORD32 *pu4_res_pred_sad, UWORD8 *pu1_residual_prediction_flag,
                              UWORD32 u4_winning_sad)
{
    res_pred_constants_t *ps_res_pred_constants = &ps_res_pred_ctxt->s_res_pred_constants;
    res_pred_state_t *ps_res_pred_state = (res_pred_state_t *) ps_res_pred_constants->pv_state;
    pu4_res_pred_sad[0] = ps_res_pred_state->pf_get_sad_with_residual_pred(
        &ps_src->as_component_bufs[Y], &ps_pred->as_component_bufs[Y],
        &ps_res->as_component_bufs[Y], MB_SIZE, MB_SIZE);

    pu1_residual_prediction_flag[0] = pu4_res_pred_sad[0] < u4_winning_sad;
}

void isvce_update_res_pred_info(isvce_process_ctxt_t *ps_proc)
{
    if(ps_proc->s_svc_params.u1_num_spatial_layers > 1)
    {
        svc_res_pred_ctxt_t *ps_res_pred_ctxt = ps_proc->ps_res_pred_ctxt;
        res_pred_constants_t *ps_res_pred_constants = &ps_res_pred_ctxt->s_res_pred_constants;
        res_pred_state_t *ps_res_pred_state = (res_pred_state_t *) ps_res_pred_constants->pv_state;
        res_pred_layer_state_t *ps_layer_state =
            &ps_res_pred_state->ps_layer_state[ps_proc->u1_spatial_layer_id];

        WORD8 i1_is_intra = ps_proc->ps_mb_info->u1_is_intra;

        WORD8 *pi1_mb_mode =
            &ps_layer_state->pi1_mb_mode[ps_proc->i4_mb_x +
                                         (ps_proc->i4_mb_y * (ps_layer_state->i4_mb_mode_stride))];

        if(ps_proc->ps_mb_info->u1_base_mode_flag == 1 && i1_is_intra)
        {
            *pi1_mb_mode = SVC_IBL_MB;
        }
        else
        {
            if(i1_is_intra)
            {
                *pi1_mb_mode = SVC_INTRA_MB;
            }
            else
            {
                *pi1_mb_mode = SVC_INTER_MB;
            }
        }
    }
}
