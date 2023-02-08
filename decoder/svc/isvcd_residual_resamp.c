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
 *  isvcd_residual_resamp.c
 *
 * @brief
 *  Contains routines that resample for SVC resampling
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_residual_chroma_dyadic_alt()
 *  - isvcd_residual_chroma_dyadic()
 *  - isvcd_residual_luma_dyadic()
 *  - isvcd_ref_layer_ptr_incr()
 *  - isvcd_residual_reflayer_const_non_boundary_mb()
 *  - isvcd_residual_reflayer_const_boundary_mb()
 *  - isvcd_residual_reflayer_const()
 *  - isvcd_interpolate_residual()
 *  - isvcd_residual_samp_mb()
 *  - isvcd_residual_samp_mb_dyadic()
 *  - isvcd_residual_samp_populate_list()
 *  - isvcd_residual_samp_res_init()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <assert.h>
#include <string.h>
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_bitstrm.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "isvcd_structs.h"
#include "ih264d_defs.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_mb_utils.h"
#include "ih264d_deblocking.h"
#include "ih264d_dpb_manager.h"
#include "ih264d_mvpred.h"
#include "ih264d_inter_pred.h"
#include "ih264d_process_pslice.h"
#include "ih264d_error_handler.h"
#include "ih264d_cabac.h"
#include "ih264d_debug.h"
#include "ih264d_tables.h"
#include "ih264d_parse_slice.h"
#include "ih264d_utils.h"
#include "ih264d_parse_islice.h"
#include "ih264d_process_bslice.h"
#include "ih264d_process_intra_mb.h"
#include "ih264_debug.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_chroma_dyadic_alt                          */
/*                                                                           */
/*  Description   : this fucntion does the upsampling of chroma residuals for*/
/*                  Dyadic cases and specific chroma phase cases             */
/*                                                                           */
/*  Inputs        : pv_residual_samp_ctxt : Residual upsampling context      */
/*                  pu1_inp_data : input 8 bit data pointer                  */
/*                  i4_inp_data_stride : input buffer stride                 */
/*                  pi2_out_res : output 16 bit buffer pointer               */
/*                  i4_out_res_stride : Output buffer stride                 */
/*                  pu1_inp_bitmap : input packed sign bit data pointer      */
/*                  i4_inp_bitmap_stride : sign bit buffer stride            */
/*                  i4_start_bit_pos : bit position in the byte of packed    */
/*                                      sign values                          */
/*  Globals       : none                                                     */
/*  Processing    : it does the upsampling with intial phase values          */
/*                                                                           */
/*  Outputs       : Upsampled residuals for chroma                           */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_residual_chroma_dyadic_alt(void *pv_residual_samp_ctxt, UWORD16 u2_mb_x, UWORD16 u2_mb_y,
                                      mem_element_t *ps_ref_mb_mode, WORD16 *pi2_inp_data,
                                      WORD32 i4_inp_data_stride, WORD16 *pi2_out_res,
                                      WORD32 i4_out_res_stride, WORD32 i4_cr_flag)
{
    residual_sampling_ctxt_t *ps_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    /* ----------------- Processing ------------------------------- */
    {
        ref_pixel_map_t *ps_pos_phase;
        residual_samp_map_ctxt_t *ps_chroma_map;
        ref_mb_map_t *ps_x_off_len_chroma;
        ref_mb_map_t *ps_y_off_len_chroma;
        WORD32 i4_i;
        WORD16 *pi2_ref_data_byte;
        WORD32 *pi4_ref_array;
        WORD32 i4_phase1, i4_phase2;
        WORD32 i4_offset_x, i4_offset_y;
        WORD32 i4_chrm_horz_int_mode, i4_chrm_vert_int_mode;
        WORD32 i4_horz_intp_ctr = SUB_BLOCK_HEIGHT;

        ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;
        ps_x_off_len_chroma = ps_chroma_map->ps_x_offset_length;
        ps_y_off_len_chroma = ps_chroma_map->ps_y_offset_length;

        /* get the actual offset for the buffers */
        i4_offset_x = ps_x_off_len_chroma[u2_mb_x].i2_offset;
        i4_offset_y = ps_y_off_len_chroma[u2_mb_y].i2_offset;

        {
            UWORD8 u1_mask;
            WORD32 i4_mb_x, i4_mb_y;
            WORD32 i4_chrm_nnz;
            WORD32 i4_num_element_stride;
            inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms, *ps_inter_lyr_mb_prms_curr;

            u1_mask = (SVCD_TRUE == i4_cr_flag) ? 0xF0 : 0x0F;

            /* Top Left */
            i4_mb_x = i4_offset_x >> 3;
            i4_mb_y = i4_offset_y >> 3;

            /* get the location of the byte which has the current mb mode */
            ps_inter_lyr_mb_prms = ps_ref_mb_mode->pv_buffer;
            i4_num_element_stride = ps_ref_mb_mode->i4_num_element_stride;

            ps_inter_lyr_mb_prms_curr = ps_inter_lyr_mb_prms + i4_mb_x;
            ps_inter_lyr_mb_prms_curr += i4_mb_y * i4_num_element_stride;
            i4_chrm_nnz = ps_inter_lyr_mb_prms_curr->u1_chroma_nnz & u1_mask;

            /* Top Right */
            i4_mb_x = (i4_offset_x + 4) >> 3;
            ps_inter_lyr_mb_prms_curr = ps_inter_lyr_mb_prms + i4_mb_x;
            ps_inter_lyr_mb_prms_curr += i4_mb_y * i4_num_element_stride;
            i4_chrm_nnz |= ps_inter_lyr_mb_prms_curr->u1_chroma_nnz & u1_mask;

            /* Bottom Left */
            i4_mb_x = i4_offset_x >> 3;
            i4_mb_y = (i4_offset_y + 4) >> 3;
            ps_inter_lyr_mb_prms_curr = ps_inter_lyr_mb_prms + i4_mb_x;
            ps_inter_lyr_mb_prms_curr += i4_mb_y * i4_num_element_stride;
            i4_chrm_nnz |= ps_inter_lyr_mb_prms_curr->u1_chroma_nnz & u1_mask;

            /* Bottom Right */
            i4_mb_x = (i4_offset_x + 4) >> 3;
            ps_inter_lyr_mb_prms_curr = ps_inter_lyr_mb_prms + i4_mb_x;
            ps_inter_lyr_mb_prms_curr += i4_mb_y * i4_num_element_stride;

            i4_chrm_nnz |= ps_inter_lyr_mb_prms_curr->u1_chroma_nnz & u1_mask;
            if(0 == i4_chrm_nnz)
            {
                return;
            }
        }

        i4_chrm_horz_int_mode = ps_lyr_ctxt->i4_chrm_horz_int_mode;
        i4_chrm_vert_int_mode = ps_lyr_ctxt->i4_chrm_vert_int_mode;

        if(0 == i4_chrm_horz_int_mode)
        {
            if(i4_offset_x >= 0)
            {
                pi2_inp_data++;
            }
        }

        if(0 == i4_chrm_vert_int_mode)
        {
            if(i4_offset_y >= 0)
            {
                pi2_inp_data += i4_inp_data_stride;
            }
        }
        else
        {
            /* extra additional row of interpolation required for this case */
            i4_horz_intp_ctr++;
        }

        /* ----------- Horizontal Interpolation ---------------- */
        pi2_ref_data_byte = pi2_inp_data;
        ps_pos_phase = ps_lyr_ctxt->s_chroma_map_ctxt.ps_x_pos_phase;

        pi4_ref_array = (WORD32 *) ps_ctxt->pi2_refarray_buffer;
        i4_phase1 = ps_pos_phase[0].i2_phase;
        i4_phase2 = (i4_phase1 + 8) & 0x0F;

        /* interchange the phase values for corner case */
        if(1 == i4_chrm_horz_int_mode)
        {
            WORD32 i4_temp;
            i4_temp = i4_phase1;
            i4_phase1 = i4_phase2;
            i4_phase2 = i4_temp;
        }

        for(i4_i = 0; i4_i < i4_horz_intp_ctr; i4_i++)
        {
            WORD16 i2_coeff1, i2_coeff2;

            i2_coeff1 = (WORD16) (pi2_ref_data_byte[0]);

            if(0 == i4_chrm_horz_int_mode)
            {
                /* populate the first inter sample */
                *pi4_ref_array++ = i2_coeff1 << 4;
            }

            /* unroll count 1 */
            i2_coeff2 = (WORD16) (pi2_ref_data_byte[2]);

            /* populate 2 samples based on current coeffs */
            *pi4_ref_array++ = ((16 - i4_phase2) * i2_coeff1 + i4_phase2 * i2_coeff2);
            *pi4_ref_array++ = ((16 - i4_phase1) * i2_coeff1 + i4_phase1 * i2_coeff2);

            /* unroll count 2 */
            i2_coeff1 = (WORD16) (pi2_ref_data_byte[4]);

            /* populate 2 samples based on current coeffs */
            *pi4_ref_array++ = ((16 - i4_phase2) * i2_coeff2 + i4_phase2 * i2_coeff1);
            *pi4_ref_array++ = ((16 - i4_phase1) * i2_coeff2 + i4_phase1 * i2_coeff1);

            /* unroll count 3 */
            i2_coeff2 = (WORD16) (pi2_ref_data_byte[6]);

            /* populate 2 samples based on current coeffs */
            *pi4_ref_array++ = ((16 - i4_phase2) * i2_coeff1 + i4_phase2 * i2_coeff2);
            *pi4_ref_array++ = ((16 - i4_phase1) * i2_coeff1 + i4_phase1 * i2_coeff2);

            /* populate the last inter sample */
            *pi4_ref_array++ = i2_coeff2 << 4;

            if(1 == i4_chrm_horz_int_mode)
            {
                i2_coeff1 = (WORD16) (pi2_ref_data_byte[4]);

                /* populate the last inter sample */
                *pi4_ref_array++ = i2_coeff1 << 4;
            }

            /* vertical loop updates */
            pi2_ref_data_byte = pi2_inp_data + ((i4_i + 1) * i4_inp_data_stride);
        }

        /* ----------- Vertical Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) ps_ctxt->pi2_refarray_buffer;
        ps_pos_phase = ps_lyr_ctxt->s_chroma_map_ctxt.ps_y_pos_phase;
        i4_phase1 = ps_pos_phase[0].i2_phase;
        i4_phase2 = (i4_phase1 + 8) & 0x0F;

        /* interchange the phase values for corner case */
        if(0 != i4_chrm_vert_int_mode)
        {
            WORD32 i4_temp;
            i4_temp = i4_phase1;
            i4_phase1 = i4_phase2;
            i4_phase2 = i4_temp;
        }

        for(i4_i = 0; i4_i < BLOCK_WIDTH; i4_i++)
        {
            WORD16 *pi2_out;
            WORD32 *pi4_ref_array_temp;
            WORD32 i4_horz_samp_1, i4_horz_samp_2;
            pi2_out = pi2_out_res;
            pi4_ref_array_temp = pi4_ref_array;

            /* populate the first inter sample */
            i4_horz_samp_1 = *pi4_ref_array_temp;
            pi4_ref_array_temp += BLOCK_WIDTH;

            if(1 != i4_chrm_vert_int_mode)
            {
                *pi2_out = (i4_horz_samp_1 + 8) >> 4;
                pi2_out += i4_out_res_stride;
            }

            if(2 == i4_chrm_vert_int_mode)
            {
                i4_horz_samp_1 = *pi4_ref_array_temp;
                pi4_ref_array_temp += BLOCK_WIDTH;
                *pi2_out = (i4_horz_samp_1 + 8) >> 4;
                pi2_out += i4_out_res_stride;
            }

            /* unroll count 1 */
            i4_horz_samp_2 = *pi4_ref_array_temp;
            pi4_ref_array_temp += BLOCK_WIDTH;

            /* populate 2 samples based on current coeffs */
            *pi2_out = ((16 - i4_phase2) * i4_horz_samp_1 + i4_phase2 * i4_horz_samp_2 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 2 */
            *pi2_out = ((16 - i4_phase1) * i4_horz_samp_1 + i4_phase1 * i4_horz_samp_2 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 3 */
            i4_horz_samp_1 = *pi4_ref_array_temp;
            pi4_ref_array_temp += BLOCK_WIDTH;

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

            if(2 != i4_chrm_vert_int_mode)
            {
                /* populate the last inter sample */
                *pi2_out = (i4_horz_samp_2 + 8) >> 4;

                if(1 == i4_chrm_vert_int_mode)
                {
                    pi2_out += i4_out_res_stride;
                    pi4_ref_array_temp += BLOCK_WIDTH;
                    i4_horz_samp_1 = *pi4_ref_array_temp;

                    /* populate the last inter sample */
                    *pi2_out = (i4_horz_samp_1 + 8) >> 4;
                }
            }
            /* horizontal loop updates */
            pi4_ref_array++;
            pi2_out_res += 2;
        }
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_chroma_dyadic                              */
/*                                                                           */
/*  Description   : this fucntion does the upsampling of chroma residuals for*/
/*                  Dyadic cases                                             */
/*                                                                           */
/*  Inputs        : pv_residual_samp_ctxt : Residual upsampling context      */
/*                  pu1_inp_data : input 8 bit data pointer                  */
/*                  i4_inp_data_stride : input buffer stride                 */
/*                  pi2_out_res : output 16 bit buffer pointer               */
/*                  i4_out_res_stride : Output buffer stride                 */
/*                  pu1_inp_bitmap : input packed sign bit data pointer      */
/*                  i4_inp_bitmap_stride : sign bit buffer stride            */
/*                  i4_start_bit_pos : bit position in the byte of packed    */
/*                                      sign values                          */
/*  Globals       : none                                                     */
/*  Processing    : it does the upsampling with intial phase values          */
/*                                                                           */
/*  Outputs       : Upsampled residuals for chroma                           */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_residual_chroma_dyadic(void *pv_residual_samp_ctxt, WORD16 *pi2_inp_data,
                                  WORD32 i4_inp_data_stride, WORD16 *pi2_out_res,
                                  WORD32 i4_out_res_stride)
{
    residual_sampling_ctxt_t *ps_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    /* ----------------- Processing ------------------------------- */
    {
        WORD32 i4_i;
        WORD16 *pi2_ref_data_byte;
        WORD32 *pi4_ref_array;
        ref_pixel_map_t *ps_pos_phase;
        WORD32 i4_phase1, i4_phase2;

        pi2_ref_data_byte = pi2_inp_data;
        ps_pos_phase = ps_lyr_ctxt->s_chroma_map_ctxt.ps_x_pos_phase;

        /* ----------- Horizontal Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) ps_ctxt->pi2_refarray_buffer;
        i4_phase1 = ps_pos_phase[0].i2_phase;
        i4_phase2 = (i4_phase1 + 8) & 0x0F;

        for(i4_i = 0; i4_i < SUB_BLOCK_HEIGHT; i4_i++)
        {
            WORD16 i2_coeff1, i2_coeff2;
            i2_coeff1 = (WORD16) (pi2_ref_data_byte[0]);

            /* populate the first inter sample */
            *pi4_ref_array++ = i2_coeff1 << 4;

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

            /* populate the last inter sample */
            *pi4_ref_array++ = i2_coeff2 << 4;

            /* vertical loop uopdates */
            pi2_ref_data_byte = pi2_inp_data + ((i4_i + 1) * i4_inp_data_stride);
        }

        /* ----------- Vertical Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) ps_ctxt->pi2_refarray_buffer;
        ps_pos_phase = ps_lyr_ctxt->s_chroma_map_ctxt.ps_y_pos_phase;
        i4_phase1 = ps_pos_phase[0].i2_phase;
        i4_phase2 = (i4_phase1 + 8) & 0x0F;

        for(i4_i = 0; i4_i < BLOCK_WIDTH; i4_i++)
        {
            WORD16 *pi2_out;
            WORD32 *pi4_ref_array_temp;
            WORD32 i4_horz_samp_1, i4_horz_samp_2;
            pi2_out = pi2_out_res;
            pi4_ref_array_temp = pi4_ref_array;

            /* populate the first inter sample */
            i4_horz_samp_1 = *pi4_ref_array_temp;
            pi4_ref_array_temp += BLOCK_WIDTH;
            *pi2_out = (i4_horz_samp_1 + 8) >> 4;
            pi2_out += i4_out_res_stride;

            /* unroll count 1 */
            i4_horz_samp_2 = *pi4_ref_array_temp;
            pi4_ref_array_temp += BLOCK_WIDTH;

            /* populate 2 samples based on current coeffs */
            *pi2_out = ((16 - i4_phase2) * i4_horz_samp_1 + i4_phase2 * i4_horz_samp_2 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 2 */
            *pi2_out = ((16 - i4_phase1) * i4_horz_samp_1 + i4_phase1 * i4_horz_samp_2 + 128) >> 8;
            pi2_out += i4_out_res_stride;

            /* unroll count 3 */
            i4_horz_samp_1 = *pi4_ref_array_temp;
            pi4_ref_array_temp += BLOCK_WIDTH;

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

            /* populate the last inter sample */
            *pi2_out = (i4_horz_samp_2 + 8) >> 4;

            /* horizontal loop updates */
            pi4_ref_array++;
            pi2_out_res += 2;
        }
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_dyadic                                */
/*                                                                           */
/*  Description   : this fucntion does the upsampling of luma residuals for  */
/*                  Dyadic cases                                             */
/*                                                                           */
/*  Inputs        : pv_residual_samp_ctxt : Residual upsampling context      */
/*                  pu1_inp_data : input 8 bit data pointer                  */
/*                  i4_inp_data_stride : input buffer stride                 */
/*                  pi2_out_res : output 16 bit buffer pointer               */
/*                  i4_out_res_stride : Output buffer stride                 */
/*                  pu1_inp_bitmap : input packed sign bit data pointer      */
/*                  i4_inp_bitmap_stride : sign bit buffer stride            */
/*                  ps_ref_mb_mode : reference mb mode pointer of base layer */
/*                  ps_coord : mb co-ordinate pointer                        */
/*  Globals       : none                                                     */
/*  Processing    : it does the upsampling with fixed phase values and       */
/*                  reference layer transform size                           */
/*  Outputs       : Upsampled residuals for luma                             */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_residual_luma_dyadic(void *pv_residual_samp_ctxt, WORD16 *pi2_inp_data,
                                WORD32 i4_inp_data_stride, WORD16 *pi2_out_res,
                                WORD32 i4_out_res_stride, mem_element_t *ps_ref_mb_mode,
                                UWORD16 u2_mb_x, UWORD16 u2_mb_y, WORD32 i4_ref_nnz,
                                WORD32 i4_ref_tx_size)

{
    WORD16 *pi2_refarray_buffer;
    WORD32 i4_blk_ctr;
    residual_sampling_ctxt_t *ps_ctxt;

    UNUSED(ps_ref_mb_mode);
    UNUSED(u2_mb_x);
    UNUSED(u2_mb_y);

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    pi2_refarray_buffer = ps_ctxt->pi2_refarray_buffer;

    /* based on transform size the counter and interpolation width and */
    /* height are intialised as follows                                */
    if((i4_ref_tx_size) && (0 != i4_ref_nnz))
    {
        WORD16 *pi2_ref_data_byte;
        WORD32 *pi4_ref_array;
        WORD32 i4_i, i4_j;

        pi2_ref_data_byte = pi2_inp_data;

        /* ----------- Horizontal Interpolation ---------------- */
        pi4_ref_array = (WORD32 *) pi2_refarray_buffer;
        for(i4_i = 0; i4_i < BLOCK_HEIGHT; i4_i++)
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

        for(i4_i = 0; i4_i < MB_WIDTH; i4_i++)
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
                pi4_ref_array_temp += MB_WIDTH;
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
        for(i4_blk_ctr = 0; i4_blk_ctr < 4; i4_blk_ctr++)
        {
            WORD16 *pi2_ref_data_byte;
            WORD32 *pi4_ref_array;
            WORD32 i4_i;

            /* if reference layer is not coded then no processing */
            if(0 != (i4_ref_nnz & 0x1))
            {
                pi2_ref_data_byte = pi2_inp_data;

                /* ----------- Horizontal Interpolation ---------------- */
                pi4_ref_array = (WORD32 *) pi2_refarray_buffer;

                for(i4_i = 0; i4_i < SUB_BLOCK_HEIGHT; i4_i++)
                {
                    WORD16 i2_coeff1, i2_coeff2;
                    i2_coeff1 = (WORD16) (*pi2_ref_data_byte++);

                    /* populate the first inter sample */
                    *pi4_ref_array++ = i2_coeff1 << 2;

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

                    /* populate the last inter sample */
                    *pi4_ref_array++ = i2_coeff2 << 2;

                    /* vertical loop uopdates */
                    pi2_ref_data_byte = pi2_inp_data + ((i4_i + 1) * i4_inp_data_stride);
                }

                /* ----------- Vertical Interpolation ---------------- */
                pi4_ref_array = (WORD32 *) pi2_refarray_buffer;

                for(i4_i = 0; i4_i < BLOCK_WIDTH; i4_i++)
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

                    /* unroll loop count 1 */
                    pi4_ref_array_temp += BLOCK_WIDTH;
                    i4_horz_samp_2 = *pi4_ref_array_temp;

                    /* populate 2 samples based on current coeffs */
                    *pi2_out =
                        ((i4_horz_samp_1 << 1) + (i4_horz_samp_1) + (i4_horz_samp_2) + 8) >> 4;
                    pi2_out += i4_out_res_stride;
                    *pi2_out =
                        ((i4_horz_samp_2 << 1) + (i4_horz_samp_2) + (i4_horz_samp_1) + 8) >> 4;
                    pi2_out += i4_out_res_stride;

                    /* unroll loop count 2 */
                    pi4_ref_array_temp += BLOCK_WIDTH;
                    i4_horz_samp_1 = *pi4_ref_array_temp;

                    /* populate 2 samples based on current coeffs */
                    *pi2_out =
                        ((i4_horz_samp_2 << 1) + (i4_horz_samp_2) + (i4_horz_samp_1) + 8) >> 4;
                    pi2_out += i4_out_res_stride;
                    *pi2_out =
                        ((i4_horz_samp_1 << 1) + (i4_horz_samp_1) + (i4_horz_samp_2) + 8) >> 4;
                    pi2_out += i4_out_res_stride;

                    /* unroll loop count 3 */
                    pi4_ref_array_temp += BLOCK_WIDTH;
                    i4_horz_samp_2 = *pi4_ref_array_temp;

                    /* populate 2 samples based on current coeffs */
                    *pi2_out =
                        ((i4_horz_samp_1 << 1) + (i4_horz_samp_1) + (i4_horz_samp_2) + 8) >> 4;
                    pi2_out += i4_out_res_stride;
                    *pi2_out =
                        ((i4_horz_samp_2 << 1) + (i4_horz_samp_2) + (i4_horz_samp_1) + 8) >> 4;
                    pi2_out += i4_out_res_stride;

                    /* populate the last inter sample */
                    *pi2_out = (i4_horz_samp_2 + 2) >> 2;

                    /* horizontal loop updates */
                    pi4_ref_array++;
                    pi2_out_res++;
                }
            }
            else
            {
                pi2_out_res += BLOCK_WIDTH;
            }

            /* Block level loop updates */
            if(1 == i4_blk_ctr)
            {
                pi2_inp_data -= SUB_BLOCK_WIDTH;
                pi2_inp_data += (i4_inp_data_stride * SUB_BLOCK_HEIGHT);
                pi2_out_res -= MB_WIDTH;
                pi2_out_res += (i4_out_res_stride * BLOCK_HEIGHT);
                i4_ref_nnz >>= 2;
            }
            else
            {
                pi2_inp_data += SUB_BLOCK_WIDTH;
            }

            i4_ref_nnz >>= 1;

        } /* end of loop over all the blocks */
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_ref_layer_ptr_incr                                 */
/*                                                                           */
/*  Description   :  this function returns the pointer increments for        */
/*                   the operand2 of the bilinear interpolation              */
/*  Inputs        : pi1_ref_mb_modes : reference mb modes                    */
/*                  i4_ref_mode_stride : mb mode buffer stride               */
/*                    i4_element_size : size of reference mb mode            */
/*                  i4_x_offset : ref offset x                               */
/*                  i4_y_offset : ref offset y                               */
/*                    i4_refary_wd : reference array width                   */
/*                    i4_refary_ht : reference array height                  */
/*                  pu1_ref_x_ptr_incr : ptr increment buffer for x          */
/*                  pu1_ref_y_ptr_incr : ptr increment buffer for y          */
/*                    i4_chroma_flag : chroma processing flag                */
/*  Globals       : none                                                     */
/*  Processing    : it calculates the increment as per the transform size    */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         18 08 2021   Kishore              creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_ref_layer_ptr_incr(WORD8 *pi1_ref_mb_modes, WORD32 i4_ref_mode_stride,
                                WORD32 i4_element_size, WORD32 i4_x_offset, WORD32 i4_y_offset,
                                WORD32 i4_refary_wd, WORD32 i4_refary_ht,
                                UWORD8 *pu1_ref_x_ptr_incr, UWORD8 *pu1_ref_y_ptr_incr,
                                WORD32 i4_chroma_flag)
{
    WORD32 i4_x, i4_y;
    WORD32 i4_x_idx, i4_y_idx;
    WORD32 i4_prev_x, i4_prev_y;
    WORD32 i4_const_val;
    WORD32 i4_pos_x, i4_pos_y;
    WORD32 i4_trans_size;
    WORD32 i4_act_ary_wd, i4_act_ary_ht;
    WORD32 i4_and_const;
    UWORD8 *pu1_incr_x, *pu1_incr_y;
    WORD32 i4_mb_sft;
    WORD32 i4_mb_x, i4_mb_y;
    WORD8 *pi1_ref_mb_modes_incr;
    WORD8 *pi1_ref_mb_modes_incr_temp;
    inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms;
    WORD32 i4_mb_x_strt, i4_mb_y_strt;
    WORD32 i4_mb_quard1_part_x, i4_mb_quard1_part_y;
    WORD32 i4_x_ref, i4_y_ref;
    WORD32 i4_tx_size, i4_tx_size_q0, i4_tx_size_q1, i4_tx_size_q2, i4_tx_size_q3;
    WORD8 i1_mb_mode_q0, i1_mb_mode_q1, i1_mb_mode_q2, i1_mb_mode_q3;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;

    i4_mb_wd = MB_WIDTH >> i4_chroma_flag;
    i4_mb_ht = MB_HEIGHT >> i4_chroma_flag;

    /* Memset to 1 the increment buffers */
    memset(pu1_ref_x_ptr_incr, 1, (i4_refary_wd * i4_refary_ht));
    memset(pu1_ref_y_ptr_incr, 1, (i4_refary_wd * i4_refary_ht));

    /* Initialise actual width and height */
    i4_act_ary_wd = i4_refary_wd;
    i4_act_ary_ht = i4_refary_ht;

    /* Initialize x and y */
    i4_x = 0;
    i4_y = 0;
    i4_prev_y = 0;
    i4_mb_sft = (MB_WIDTH_SHIFT - i4_chroma_flag);

    /* Loop over all MBs in the reference array */
    if(0 == i4_chroma_flag)
    {
        i4_x_ref = i4_x_offset + 0;
        i4_y_ref = i4_y_offset + 0;
        i4_mb_x_strt = i4_x_ref % i4_mb_wd;
        i4_mb_y_strt = i4_y_ref % i4_mb_ht;
        i4_mb_quard1_part_x = i4_mb_wd - i4_mb_x_strt;
        i4_mb_quard1_part_y = i4_mb_ht - i4_mb_y_strt;

        if(!(i4_mb_quard1_part_x >= 0))
        {
            return NOT_OK;
        }
        if(!(i4_mb_quard1_part_y >= 0))
        {
            return NOT_OK;
        }

        /* Take care of negative offsets */
        if(i4_x_ref > 0)
        {
            i4_mb_x = (i4_x_ref >> i4_mb_sft);
        }
        else
        {
            i4_mb_x = 0;
        }
        if(i4_y_ref > 0)
        {
            i4_mb_y = (i4_y_ref >> i4_mb_sft);
        }
        else
        {
            i4_mb_y = 0;
        }

        /* get the location of the byte which has the current mb mode */
        pi1_ref_mb_modes_incr = pi1_ref_mb_modes + (i4_mb_y * i4_ref_mode_stride * i4_element_size);
        pi1_ref_mb_modes_incr += (i4_mb_x * i4_element_size);
        ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes_incr;
        i1_mb_mode_q0 = ps_inter_lyr_mb_prms->i1_mb_mode;
        i4_tx_size_q0 =
            (i1_mb_mode_q0 <= SVC_INTER_MB)
                ? ((ps_inter_lyr_mb_prms->i1_tx_size < 0) ? 1 : ps_inter_lyr_mb_prms->i1_tx_size)
                : 1;

        pi1_ref_mb_modes_incr_temp = pi1_ref_mb_modes_incr;
        if(i4_mb_quard1_part_x > 0)
        {
            pi1_ref_mb_modes_incr_temp = pi1_ref_mb_modes_incr + i4_element_size;
            ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes_incr_temp;
            i1_mb_mode_q1 = ps_inter_lyr_mb_prms->i1_mb_mode;
            i4_tx_size_q1 =
                (i1_mb_mode_q1 <= SVC_INTER_MB)
                    ? ((ps_inter_lyr_mb_prms->i1_tx_size < 0) ? 1
                                                              : ps_inter_lyr_mb_prms->i1_tx_size)
                    : 1;
        }

        if(i4_mb_quard1_part_y > 0)
        {
            pi1_ref_mb_modes_incr_temp =
                pi1_ref_mb_modes_incr + (i4_ref_mode_stride * i4_element_size);
            ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes_incr_temp;
            i1_mb_mode_q2 = ps_inter_lyr_mb_prms->i1_mb_mode;
            i4_tx_size_q2 =
                (i1_mb_mode_q2 <= SVC_INTER_MB)
                    ? ((ps_inter_lyr_mb_prms->i1_tx_size < 0) ? 1
                                                              : ps_inter_lyr_mb_prms->i1_tx_size)
                    : 1;
        }

        if((i4_mb_quard1_part_x > 0) && (i4_mb_quard1_part_y > 0))
        {
            pi1_ref_mb_modes_incr_temp =
                pi1_ref_mb_modes_incr + (i4_ref_mode_stride * i4_element_size) + i4_element_size;
            ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes_incr_temp;
            i1_mb_mode_q3 = ps_inter_lyr_mb_prms->i1_mb_mode;
            i4_tx_size_q3 =
                (i1_mb_mode_q3 <= SVC_INTER_MB)
                    ? ((ps_inter_lyr_mb_prms->i1_tx_size < 0) ? 1
                                                              : ps_inter_lyr_mb_prms->i1_tx_size)
                    : 1;
        }

        do
        {
            WORD32 i4_idx;
            WORD32 i4_wd, i4_ht;
            WORD32 i4_max_pos_x, i4_max_pos_y;

            i4_prev_x = i4_x;
            i4_x_ref = i4_x_offset + i4_x;
            i4_y_ref = i4_y_offset + i4_y;
            i4_tx_size = i4_tx_size_q0;
            if(i4_x >= i4_mb_quard1_part_x)
            {
                if(i4_y < i4_mb_quard1_part_y)
                {
                    i4_tx_size = i4_tx_size_q1;
                }
                else if(i4_y >= i4_mb_quard1_part_y)
                {
                    i4_tx_size = i4_tx_size_q3;
                }
            }
            else if(i4_x < i4_mb_quard1_part_x)
            {
                if(i4_y >= i4_mb_quard1_part_y)
                {
                    i4_tx_size = i4_tx_size_q2;
                }
            }

            /* Get the transform size as 4 or 8 */
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
                {
                    /* Fill the block edge with 0s */
                    *pu1_incr = 0;
                    pu1_incr += i4_refary_wd;
                }

                /* Updates */
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
                {
                    /* Fill the block edge with 0s */
                    *pu1_incr = 0;
                    pu1_incr++;
                }

                /* Updates */
                i4_pos_y += i4_trans_size;
                pu1_incr_y += i4_trans_size * i4_refary_wd;
                i4_y_idx += MIN(i4_trans_size, (i4_act_ary_ht - i4_y_idx));
            }

            /* Loop updates */
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

    } /* End of if 0 == i4_chroma_flag */
    else
    {
        /* Set the transform size as 4 */
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
                {
                    /* Fill the block edge with 0s */
                    *pu1_incr = 0;
                    pu1_incr += i4_refary_wd;
                }

                /* Updates */
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
                {
                    /* Fill the block edge with 0s */
                    *pu1_incr = 0;
                    pu1_incr++;
                }

                /* Updates */
                i4_pos_y += i4_trans_size;
                pu1_incr_y += i4_trans_size * i4_refary_wd;
                i4_y_idx += MIN(i4_trans_size, (i4_act_ary_ht - i4_y_idx));
            }

            /* Loop updates */
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

    } /* End of chroma */
    return OK;

} /* End of "isvcd_ref_layer_ptr_incr" */

void isvcd_residual_reflayer_const_non_boundary_mb(
    WORD16 *pi2_inp_data, WORD32 i4_inp_data_stride, WORD16 *pi2_ref_array, WORD32 i4_refarray_wd,
    WORD32 i4_refarray_ht, WORD32 i4_ref_mb_type_q0, WORD32 i4_ref_mb_type_q1,
    WORD32 i4_ref_mb_type_q2, WORD32 i4_ref_mb_type_q3, WORD32 i4_mb_quard1_part_x,
    WORD32 i4_mb_quard1_part_y, WORD32 i4_chroma_flag)
{
    WORD32 i4_x_ref, i4_y_ref;
    WORD32 i4_x, i4_y;
    WORD32 i4_ref_mb_type;
    WORD16 *pi2_ref_data_byte;
    WORD16 *pi2_ref_array_temp;

    for(i4_y = 0; i4_y < i4_refarray_ht; i4_y++)
    {
        for(i4_x = 0; i4_x < i4_refarray_wd; i4_x++)
        {
            i4_y_ref = i4_y;
            i4_x_ref = i4_x;

            i4_ref_mb_type = i4_ref_mb_type_q0;
            if(i4_x >= i4_mb_quard1_part_x)
            {
                if(i4_y < i4_mb_quard1_part_y)
                {
                    i4_ref_mb_type = i4_ref_mb_type_q1;
                }
                else if(i4_y >= i4_mb_quard1_part_y)
                {
                    i4_ref_mb_type = i4_ref_mb_type_q3;
                }
            }
            else if(i4_x < i4_mb_quard1_part_x)
            {
                if(i4_y >= i4_mb_quard1_part_y)
                {
                    i4_ref_mb_type = i4_ref_mb_type_q2;
                }
            }

            /****************************************************************/
            /* Reference layer Residual Buffer is maintained as 8-bit data  */
            /* Buffer and 1-bit sign bit packed buffer. Sign Byte is read   */
            /* and sign of the data sample is extracted depending upon bit  */
            /* postition.                                                   */
            /****************************************************************/

            /* update the buffer pointers to appropriate locations */
            pi2_ref_array_temp = pi2_ref_array + i4_x;
            pi2_ref_array_temp += i4_y * i4_refarray_wd;

            /* extract the residual value and fill the buffer */
            if(SVC_INTER_MB == i4_ref_mb_type)
            {
                /* derive the reference data pointers */
                pi2_ref_data_byte = pi2_inp_data + (i4_x_ref << i4_chroma_flag);
                pi2_ref_data_byte += i4_y_ref * i4_inp_data_stride;

                /* store the residual value */
                *pi2_ref_array_temp = (WORD16) (*pi2_ref_data_byte);
            }
            else
            {
                /* if non inter MB then store the 0 */
                *pi2_ref_array_temp = 0;
            }
        }
    }
}

void isvcd_residual_reflayer_const_boundary_mb(WORD16 *pi2_inp_data, WORD32 i4_inp_data_stride,
                                               WORD16 *pi2_ref_array, WORD32 i4_refarray_wd,
                                               WORD32 i4_refarray_ht, WORD32 i4_ref_wd,
                                               WORD32 i4_ref_ht, WORD32 i4_x_offset,
                                               WORD32 i4_y_offset, WORD32 i4_ref_mb_type_q0,
                                               WORD32 i4_ref_mb_type_q1, WORD32 i4_ref_mb_type_q2,
                                               WORD32 i4_ref_mb_type_q3, WORD32 i4_mb_quard1_part_x,
                                               WORD32 i4_mb_quard1_part_y, WORD32 i4_chroma_flag)
{
    WORD32 i4_x_ref, i4_y_ref;
    WORD32 i4_x, i4_y;
    WORD16 *pi2_ref_data_byte;
    WORD16 *pi2_ref_array_temp;

    /*Quard 0*/
    for(i4_y = 0; i4_y < i4_mb_quard1_part_y; i4_y++)
    {
        for(i4_x = 0; i4_x < i4_mb_quard1_part_x; i4_x++)
        {
            i4_y_ref = MAX(0, MIN(i4_ref_ht - 1, i4_y + i4_y_offset));
            i4_x_ref = MAX(0, MIN(i4_ref_wd - 1, i4_x + i4_x_offset));

            /****************************************************************/
            /* Reference layer Residual Buffer is maintained as 8-bit data  */
            /* Buffer and 1-bit sign bit packed buffer. Sign Byte is read   */
            /* and sign of the data sample is extracted depending upon bit  */
            /* postition.                                                   */
            /****************************************************************/

            /* update the buffer pointers to appropriate locations */
            pi2_ref_array_temp = pi2_ref_array + i4_x;
            pi2_ref_array_temp += i4_y * i4_refarray_wd;

            /* extract the residual value and fill the buffer */
            if(SVC_INTER_MB == i4_ref_mb_type_q0)
            {
                /* input pointer will be pointing to (xoffset,yoffset) */
                /* So subtract the correction to reference location    */
                if(0 <= i4_x_offset)
                {
                    /* if only inside frame dimension */
                    i4_x_ref = i4_x_ref - i4_x_offset;
                }

                if(0 <= i4_y_offset)
                {
                    /* if only inside frame dimension */
                    i4_y_ref = i4_y_ref - i4_y_offset;
                }
                /* derive the reference data pointers */

                pi2_ref_data_byte = pi2_inp_data + (i4_x_ref << i4_chroma_flag);
                pi2_ref_data_byte += i4_y_ref * i4_inp_data_stride;

                /* store the residual value */
                *pi2_ref_array_temp = (WORD16) (*pi2_ref_data_byte);
            }
            else
            {
                /* if non inter MB then store the 0 */
                *pi2_ref_array_temp = 0;
            }
        }
    }

    /*Quard 1*/
    for(i4_y = 0; i4_y < i4_mb_quard1_part_y; i4_y++)
    {
        for(i4_x = i4_mb_quard1_part_x; i4_x < i4_refarray_wd; i4_x++)
        {
            i4_y_ref = MAX(0, MIN(i4_ref_ht - 1, i4_y + i4_y_offset));
            i4_x_ref = MAX(0, MIN(i4_ref_wd - 1, i4_x + i4_x_offset));

            /****************************************************************/
            /* Reference layer Residual Buffer is maintained as 8-bit data  */
            /* Buffer and 1-bit sign bit packed buffer. Sign Byte is read   */
            /* and sign of the data sample is extracted depending upon bit  */
            /* postition.                                                   */
            /****************************************************************/

            /* update the buffer pointers to appropriate locations */
            pi2_ref_array_temp = pi2_ref_array + i4_x;
            pi2_ref_array_temp += i4_y * i4_refarray_wd;

            /* extract the residual value and fill the buffer */
            if(SVC_INTER_MB == i4_ref_mb_type_q1)
            {
                /* input pointer will be pointing to (xoffset,yoffset) */
                /* So subtract the correction to reference location    */
                if(0 <= i4_x_offset)
                {
                    /* if only inside frame dimension */
                    i4_x_ref = i4_x_ref - i4_x_offset;
                }
                if(0 <= i4_y_offset)
                {
                    /* if only inside frame dimension */
                    i4_y_ref = i4_y_ref - i4_y_offset;
                }
                /* derive the reference data pointers */

                pi2_ref_data_byte = pi2_inp_data + (i4_x_ref << i4_chroma_flag);
                pi2_ref_data_byte += i4_y_ref * i4_inp_data_stride;

                /* store the residual value */
                *pi2_ref_array_temp = (WORD16) (*pi2_ref_data_byte);
            }
            else
            {
                /* if non inter MB then store the 0 */
                *pi2_ref_array_temp = 0;
            }
        }
    }

    /*Quard 2*/
    for(i4_y = i4_mb_quard1_part_y; i4_y < i4_refarray_ht; i4_y++)
    {
        for(i4_x = 0; i4_x < i4_mb_quard1_part_x; i4_x++)
        {
            i4_y_ref = MAX(0, MIN(i4_ref_ht - 1, i4_y + i4_y_offset));
            i4_x_ref = MAX(0, MIN(i4_ref_wd - 1, i4_x + i4_x_offset));

            /****************************************************************/
            /* Reference layer Residual Buffer is maintained as 8-bit data  */
            /* Buffer and 1-bit sign bit packed buffer. Sign Byte is read   */
            /* and sign of the data sample is extracted depending upon bit  */
            /* postition.                                                   */
            /****************************************************************/

            /* update the buffer pointers to appropriate locations */
            pi2_ref_array_temp = pi2_ref_array + i4_x;
            pi2_ref_array_temp += i4_y * i4_refarray_wd;

            /* extract the residual value and fill the buffer */
            if(SVC_INTER_MB == i4_ref_mb_type_q2)
            {
                /* input pointer will be pointing to (xoffset,yoffset) */
                /* So subtract the correction to reference location    */
                if(0 <= i4_x_offset)
                {
                    /* if only inside frame dimension */
                    i4_x_ref = i4_x_ref - i4_x_offset;
                }
                if(0 <= i4_y_offset)
                {
                    /* if only inside frame dimension */
                    i4_y_ref = i4_y_ref - i4_y_offset;
                }
                /* derive the reference data pointers */
                pi2_ref_data_byte = pi2_inp_data + (i4_x_ref << i4_chroma_flag);
                pi2_ref_data_byte += i4_y_ref * i4_inp_data_stride;

                /* store the residual value */
                *pi2_ref_array_temp = (WORD16) (*pi2_ref_data_byte);
            }
            else
            {
                /* if non inter MB then store the 0 */
                *pi2_ref_array_temp = 0;
            }
        }
    }

    /*Quard 3*/
    for(i4_y = i4_mb_quard1_part_y; i4_y < i4_refarray_ht; i4_y++)
    {
        for(i4_x = i4_mb_quard1_part_x; i4_x < i4_refarray_wd; i4_x++)
        {
            i4_y_ref = MAX(0, MIN(i4_ref_ht - 1, i4_y + i4_y_offset));
            i4_x_ref = MAX(0, MIN(i4_ref_wd - 1, i4_x + i4_x_offset));

            /****************************************************************/
            /* Reference layer Residual Buffer is maintained as 8-bit data  */
            /* Buffer and 1-bit sign bit packed buffer. Sign Byte is read   */
            /* and sign of the data sample is extracted depending upon bit  */
            /* postition.                                                   */
            /****************************************************************/

            /* update the buffer pointers to appropriate locations */
            pi2_ref_array_temp = pi2_ref_array + i4_x;
            pi2_ref_array_temp += i4_y * i4_refarray_wd;

            /* extract the residual value and fill the buffer */
            if(SVC_INTER_MB == i4_ref_mb_type_q3)
            {
                /* input pointer will be pointing to (xoffset,yoffset) */
                /* So subtract the correction to reference location    */
                if(0 <= i4_x_offset)
                {
                    /* if only inside frame dimension */
                    i4_x_ref = i4_x_ref - i4_x_offset;
                }
                if(0 <= i4_y_offset)
                {
                    /* if only inside frame dimension */
                    i4_y_ref = i4_y_ref - i4_y_offset;
                }
                /* derive the reference data pointers */
                pi2_ref_data_byte = pi2_inp_data + (i4_x_ref << i4_chroma_flag);
                pi2_ref_data_byte += i4_y_ref * i4_inp_data_stride;

                /* store the residual value */
                *pi2_ref_array_temp = (WORD16) (*pi2_ref_data_byte);
            }
            else
            {
                /* if non inter MB then store the 0 */
                *pi2_ref_array_temp = 0;
            }
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_reflayer_const                             */
/*                                                                           */
/*  Description   :  This function constructs the reference array buffer     */
/*                    used for residual resampling of a component in an MB   */
/*                                                                           */
/*  Inputs        : pv_residual_samp_ctxt: intra sampling context            */
/*                  pu1_inp : input (reference layer data)                   */
/*                  i4_inp_stride : input buffer stride                      */
/*                  pu1_inp_bitmap : input (reference layer sign bits)       */
/*                  i4_inp_stride : input buffer stride                      */
/*                  ps_ref_mb_mode : ref layer mb mode buffer desc           */
/*                  pi4_refarr_wd : pointer to store the reference array WD  */
/*                  pi4_refarr_ht : pointer to store the reference array HT  */
/*                  pi4_x_offset : pointer to store the reference X offset   */
/*                  pi4_y_offset : pointer to store the reference Y offset   */
/*                  ps_coord     : mb co-ordinate structure                  */
/*                  i4_chroma_flag : chroma processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : it fills the reference layer data if they are falling in */
/*                  INTRA MB region. If all the pixels are not filled  it    */
/*                  calls the border extension algorithm to fill them        */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 07 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_residual_reflayer_const(void *pv_residual_samp_ctxt, WORD16 *pi2_inp_data,
                                     WORD32 i4_inp_data_stride, mem_element_t *ps_ref_mb_mode,
                                     WORD32 *pi4_refarr_wd, mb_coord_t *ps_coord,
                                     WORD32 i4_chroma_flag)
{
    residual_sampling_ctxt_t *ps_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    WORD8 *pi1_ref_mb_modes;
    WORD32 i4_ref_mode_stride;
    WORD32 i4_element_size;
    WORD32 i4_ref_wd;
    WORD32 i4_ref_ht;
    WORD32 i4_x_offset;
    WORD32 i4_y_offset;
    WORD32 i4_refarray_wd;
    WORD32 i4_refarray_ht;
    WORD8 i1_edge_mb;
    WORD16 *pi2_ref_array;
    WORD32 i4_mb_sft;
    WORD32 i4_mb_x, i4_mb_y;
    WORD32 i4_mb_x_strt, i4_mb_y_strt;
    WORD32 i4_mb_quard1_part_x, i4_mb_quard1_part_y;
    WORD8 *pi1_ref_mb_modes_incr;
    WORD8 *pi1_ref_mb_modes_incr_temp;
    inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms;

    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;
    WORD32 i4_x_ref, i4_y_ref;
    WORD32 i4_ref_mb_type_q0, i4_ref_mb_type_q1, i4_ref_mb_type_q2, i4_ref_mb_type_q3;
    WORD8 i1_mb_mode_q0, i1_mb_mode_q1, i1_mb_mode_q2, i1_mb_mode_q3;
    WORD32 ret;

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    pi2_ref_array = ps_ctxt->pi2_refarray_buffer;

    pi1_ref_mb_modes = (WORD8 *) ps_ref_mb_mode->pv_buffer;
    i4_ref_mode_stride = ps_ref_mb_mode->i4_num_element_stride;
    i4_element_size = ps_ref_mb_mode->i4_element_size;

    if(NULL == pi1_ref_mb_modes)
    {
        return NOT_OK;
    }

    i4_mb_wd = MB_WIDTH >> i4_chroma_flag;
    i4_mb_ht = MB_HEIGHT >> i4_chroma_flag;

    /* ----------------------------------------------------------------- */
    /* Deriving the parameters required for further processing           */
    /* ----------------------------------------------------------------- */
    {
        ref_mb_map_t *ps_x_off_len;
        ref_mb_map_t *ps_y_off_len;
        WORD32 i4_mbaddr_x;
        WORD32 i4_mbaddr_y;
        WORD32 i4_base_width;
        WORD32 i4_base_height;
        residual_samp_map_ctxt_t *ps_map_ctxt;

        if(1 == i4_chroma_flag)
            ps_map_ctxt = &ps_lyr_ctxt->s_chroma_map_ctxt;
        else
            ps_map_ctxt = &ps_lyr_ctxt->s_luma_map_ctxt;

        i4_mbaddr_y = ps_coord->u2_mb_y;
        i4_mbaddr_x = ps_coord->u2_mb_x;
        i4_base_width = ps_lyr_ctxt->i4_ref_width;
        i4_base_height = ps_lyr_ctxt->i4_ref_height;
        i4_ref_wd = i4_base_width >> i4_chroma_flag;
        i4_ref_ht = i4_base_height >> i4_chroma_flag;

        /* --------------------------------------------------------------------- */
        /* Extracting information from the mapping context                       */
        /* --------------------------------------------------------------------- */
        ps_x_off_len = ps_map_ctxt->ps_x_offset_length;
        ps_y_off_len = ps_map_ctxt->ps_y_offset_length;
        i4_x_offset = ps_x_off_len[i4_mbaddr_x].i2_offset;
        i4_y_offset = ps_y_off_len[i4_mbaddr_y].i2_offset;
        i4_refarray_wd = ps_x_off_len[i4_mbaddr_x].i2_length;
        i4_refarray_ht = ps_y_off_len[i4_mbaddr_y].i2_length;
    }

    /* Call the module to fill the increments based on transform blocks */
    ret = isvcd_ref_layer_ptr_incr(pi1_ref_mb_modes, i4_ref_mode_stride, i4_element_size,
                                   i4_x_offset, i4_y_offset, i4_refarray_wd, i4_refarray_ht,
                                   ps_ctxt->pu1_ref_x_ptr_incr, ps_ctxt->pu1_ref_y_ptr_incr,
                                   i4_chroma_flag);

    if(ret != OK)
    {
        return ret;
    }
    i4_mb_sft = (MB_WIDTH_SHIFT - i4_chroma_flag);

    /* --------------------------------------------------------------------- */
    /* MB Level Resampling for the MB - Pointers sent for MB in both layers  */
    /* This has been written according to the dyadic case                    */
    /* --------------------------------------------------------------------- */
    i4_y_ref = MAX(0, MIN(i4_ref_ht - 1, 0 + i4_y_offset));
    i4_x_ref = MAX(0, MIN(i4_ref_wd - 1, 0 + i4_x_offset));
    i4_mb_x_strt = i4_x_ref % i4_mb_wd;
    i4_mb_y_strt = i4_y_ref % i4_mb_ht;

    i4_mb_quard1_part_x = i4_mb_wd - i4_mb_x_strt;
    i4_mb_quard1_part_y = i4_mb_ht - i4_mb_y_strt;
    if(!(i4_mb_quard1_part_x >= 0))
    {
        return NOT_OK;
    }
    if(!(i4_mb_quard1_part_y >= 0))
    {
        return NOT_OK;
    }

    i4_mb_x = (i4_x_ref >> i4_mb_sft);
    i4_mb_y = (i4_y_ref >> i4_mb_sft);

    /* get the location of the byte which has the current mb mode */
    pi1_ref_mb_modes_incr = pi1_ref_mb_modes + (i4_mb_y * i4_ref_mode_stride * i4_element_size);
    pi1_ref_mb_modes_incr += (i4_mb_x * i4_element_size);
    ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes_incr;
    i1_mb_mode_q0 = ps_inter_lyr_mb_prms->i1_mb_mode;
    i1_mb_mode_q1 = i1_mb_mode_q0;
    i1_mb_mode_q2 = i1_mb_mode_q0;
    i1_mb_mode_q3 = i1_mb_mode_q0;

    pi1_ref_mb_modes_incr_temp = pi1_ref_mb_modes_incr;
    if(i4_mb_quard1_part_x > 0)
    {
        pi1_ref_mb_modes_incr_temp = pi1_ref_mb_modes_incr + i4_element_size;
        ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes_incr_temp;
        i1_mb_mode_q1 = ps_inter_lyr_mb_prms->i1_mb_mode;
    }

    if(i4_mb_quard1_part_y > 0)
    {
        pi1_ref_mb_modes_incr_temp = pi1_ref_mb_modes_incr + (i4_ref_mode_stride * i4_element_size);
        ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes_incr_temp;
        i1_mb_mode_q2 = ps_inter_lyr_mb_prms->i1_mb_mode;
    }

    if((i4_mb_quard1_part_x > 0) && (i4_mb_quard1_part_y > 0))
    {
        pi1_ref_mb_modes_incr_temp =
            pi1_ref_mb_modes_incr + (i4_ref_mode_stride * i4_element_size) + i4_element_size;
        ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes_incr_temp;
        i1_mb_mode_q3 = ps_inter_lyr_mb_prms->i1_mb_mode;
    }

    i4_ref_mb_type_q0 = (i1_mb_mode_q0 <= SVC_INTER_MB) ? SVC_INTER_MB : SVC_INTRA_MB;
    i4_ref_mb_type_q1 = (i1_mb_mode_q1 <= SVC_INTER_MB) ? SVC_INTER_MB : SVC_INTRA_MB;
    i4_ref_mb_type_q2 = (i1_mb_mode_q2 <= SVC_INTER_MB) ? SVC_INTER_MB : SVC_INTRA_MB;
    i4_ref_mb_type_q3 = (i1_mb_mode_q3 <= SVC_INTER_MB) ? SVC_INTER_MB : SVC_INTRA_MB;

    i1_edge_mb = (ps_coord->u2_mb_x == 0 || ps_coord->u2_mb_y == 0 ||
                  (ps_coord->u2_mb_x == ((ps_lyr_ctxt->i4_curr_width >> MB_WIDTH_SHIFT) - 1)) ||
                  (ps_coord->u2_mb_y == ((ps_lyr_ctxt->i4_curr_height >> MB_HEIGHT_SHIFT) - 1)));
    if(i1_edge_mb)
    {
        ps_ctxt->pf_residual_reflayer_const_boundary_mb(
            pi2_inp_data, i4_inp_data_stride, pi2_ref_array, i4_refarray_wd, i4_refarray_ht,
            i4_ref_wd, i4_ref_ht, i4_x_offset, i4_y_offset, i4_ref_mb_type_q0, i4_ref_mb_type_q1,
            i4_ref_mb_type_q2, i4_ref_mb_type_q3, i4_mb_quard1_part_x, i4_mb_quard1_part_y,
            i4_chroma_flag);
    }
    else
    {
        ps_ctxt->pf_residual_reflayer_const_non_boundary_mb(
            pi2_inp_data, i4_inp_data_stride, pi2_ref_array, i4_refarray_wd, i4_refarray_ht,
            i4_ref_mb_type_q0, i4_ref_mb_type_q1, i4_ref_mb_type_q2, i4_ref_mb_type_q3,
            i4_mb_quard1_part_x, i4_mb_quard1_part_y, i4_chroma_flag);
    }
    /* store the values into the place holders */
    *pi4_refarr_wd = i4_refarray_wd;

    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_interpolate_residual                                */
/*                                                                           */
/*  Description   : This function takes the refernce array buffer and perform*/
/*                    interpolation of a component to find the residual      */
/*                     resampled value                                       */
/*  Inputs        : pv_residual_samp_ctxt : residual sampling context        */
/*                  pi2_out : output buffer pointer                          */
/*                  i4_out_stride : output buffer stride                     */
/*                  i4_refarray_wd : reference array width                   */
/*                  i4_x_offset : offset in reference layer in horz direction*/
/*                  ps_coord : current mb co-ordinate                        */
/*                  i4_chroma_flag : chroma processing flag                  */
/*  Globals       : none                                                     */
/*  Processing    : it does the interpolation in vertical direction followed */
/*                  by horizontal direction                                  */
/*  Outputs       : resampled pixels                                         */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 07 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_interpolate_residual(void *pv_residual_samp_ctxt, WORD16 *pi2_out, WORD32 i4_out_stride,
                                WORD32 i4_refarray_wd, UWORD16 u2_mb_x, UWORD16 u2_mb_y,
                                WORD32 i4_chroma_flag)
{
    residual_sampling_ctxt_t *ps_ctxt;
    residual_samp_map_ctxt_t *ps_map_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    ref_pixel_map_t *ps_x_pos_phase;
    ref_pixel_map_t *ps_y_pos_phase;

    WORD32 i4_x, i4_y;
    WORD32 i4_frm_mb_x, i4_frm_mb_y;
    WORD32 i4_temp_array_ht;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;
    WORD16 *pi2_ref_array;
    UWORD8 *pu1_ref_x_ptr_incr, *pu1_ref_y_ptr_incr;

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];
    pi2_ref_array = ps_ctxt->pi2_refarray_buffer;
    pu1_ref_x_ptr_incr = ps_ctxt->pu1_ref_x_ptr_incr;
    pu1_ref_y_ptr_incr = ps_ctxt->pu1_ref_y_ptr_incr;

    /* --------------------------------------------------------------------- */
    /* Extracting information from the mapping context                       */
    /* --------------------------------------------------------------------- */
    if(1 == i4_chroma_flag)
        ps_map_ctxt = &ps_lyr_ctxt->s_chroma_map_ctxt;
    else
        ps_map_ctxt = &ps_lyr_ctxt->s_luma_map_ctxt;

    i4_mb_wd = MB_WIDTH >> i4_chroma_flag;
    i4_mb_ht = MB_HEIGHT >> i4_chroma_flag;

    ps_x_pos_phase = ps_map_ctxt->ps_x_pos_phase;
    ps_y_pos_phase = ps_map_ctxt->ps_y_pos_phase;
    i4_temp_array_ht = i4_mb_ht;
    i4_frm_mb_y = u2_mb_y * i4_mb_ht;
    i4_frm_mb_x = u2_mb_x * i4_mb_wd;

    /* --------------------------------------------------------------------- */
    /* Loop for interpolation                                                */
    /* --------------------------------------------------------------------- */
    for(i4_y = 0; i4_y < (i4_temp_array_ht); i4_y++)
    {
        for(i4_x = 0; i4_x < (i4_mb_wd); i4_x++)
        {
            WORD32 i4_i;
            WORD32 i4_y_ref;
            WORD32 i4_y_phase;
            WORD32 i4_x_ref;
            WORD32 i4_x_phase;
            WORD32 i4_x_ref_round;
            WORD16 *pi2_out_curr;
            WORD32 ai4_temp_pred[2] = {0};
            UWORD8 *pu1_ref_y_ptr_incr_temp;
            WORD32 *pi4_temp_pred;
            UWORD8 u1_incr_y;
            WORD16 i2_res;

            /* derive the current output pointer */
            pi2_out_curr = pi2_out + (i4_x << i4_chroma_flag) + (i4_y * i4_out_stride);

            /* -------------------------------------------------------------- */
            /* Finding the offset                                             */
            /* -------------------------------------------------------------- */
            i4_y_ref = ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_ref_pos;
            i4_y_phase = ps_y_pos_phase[i4_y + i4_frm_mb_y].i2_phase;
            i4_x_ref = ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_ref_pos;
            i4_x_phase = ps_x_pos_phase[i4_x + i4_frm_mb_x].i2_phase;

            /* horizontal processing*/
            for(i4_i = 0; i4_i < 2; i4_i++)
            {
                UWORD8 *pu1_ref_x_ptr_incr_temp;
                UWORD8 u1_incr;
                WORD16 *pi2_ref_array_1, *pi2_ref_array_2;

                /* derive appropriate pointers */
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

            /* store back the final residual */
            *pi2_out_curr = i2_res;
        } /* end of loop over width */
    }     /* end of loop over height */

    return;
} /* End of Interpolation Function */
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_samp_mb                                    */
/*                                                                           */
/*  Description   : MB level function whcih perform the residual resampling  */
/*                  of data of an MB (luma and chroma insclusive)            */
/*                                                                           */
/*  Inputs        : pv_residual_samp_ctxt : residual sampling context        */
/*                  ps_ref_luma : reference layer luma data buffer desc      */
/*                  ps_ref_chroma : reference layer chroma data buffer desc  */
/*                  ps_ref_luma_bitmap : ref layer luma bit map buffer desc  */
/*                  ps_ref_chroma_bitmap : ref layer chroma bit map buff des */
/*                  ps_ref_mb_mode : ref layer mb mode map buff desc         */
/*                  ps_out_luma : current layer out luma buffer desc         */
/*                  ps_out_chroma : current layer out chroma buffer desc     */
/*                  ps_mb_coord : current mb coorinate                       */
/*  Globals       : none                                                     */
/*  Processing    : it calls the reference layer construction followed by    */
/*                   interplaotion function for luma and cb and cr           */
/*  Outputs       : inter resampled data of current MB                       */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 06 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_residual_samp_mb(void *pv_residual_samp_ctxt, mem_element_t *ps_ref_luma,
                              mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode,
                              mem_element_t *ps_out_luma, mem_element_t *ps_out_chroma,
                              UWORD16 u2_mb_x, UWORD16 u2_mb_y)
{
    /* --------------------------------------------------------------------- */
    /* I/O buffer params                                                     */
    /* --------------------------------------------------------------------- */
    residual_sampling_ctxt_t *ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    WORD16 *pi2_inp;
    WORD16 *pi2_out;
    WORD32 i4_inp_stride;
    WORD32 i4_out_stride;
    WORD32 i4_refarray_wd;
    mb_coord_t s_mb_coord = {0};
    WORD32 ret;
    s_mb_coord.u2_mb_x = u2_mb_x;
    s_mb_coord.u2_mb_y = u2_mb_y;

    /* --------------------------------------------------------------------- */
    /* LUMA PROCESSING                                                        */
    /* --------------------------------------------------------------------- */
    pi2_inp = (WORD16 *) ps_ref_luma->pv_buffer;
    pi2_out = (WORD16 *) ps_out_luma->pv_buffer;
    i4_inp_stride = ps_ref_luma->i4_num_element_stride;
    i4_out_stride = ps_out_luma->i4_num_element_stride;

    /* ------- Constructing refSampleArray ----------------------- */
    ret = isvcd_residual_reflayer_const(pv_residual_samp_ctxt, pi2_inp, i4_inp_stride,
                                        ps_ref_mb_mode, &i4_refarray_wd, &s_mb_coord, 0);

    if(ret != OK) return ret;
    /* ---- Interpolation process for Residual prediction     ------ */
    ps_ctxt->pf_interpolate_residual(pv_residual_samp_ctxt, pi2_out, i4_out_stride, i4_refarray_wd,
                                     s_mb_coord.u2_mb_x, s_mb_coord.u2_mb_y, 0);

    /* --------------------------------------------------------------------- */
    /* CHROMA PROCESSING                                                       */
    /* --------------------------------------------------------------------- */
    /* CB */
    pi2_inp = (WORD16 *) ps_ref_chroma->pv_buffer;
    pi2_out = (WORD16 *) ps_out_chroma->pv_buffer;
    i4_inp_stride = ps_ref_chroma->i4_num_element_stride;
    i4_out_stride = ps_out_chroma->i4_num_element_stride;

    /* ------- Constructing refSampleArray ----------------------- */
    ret = isvcd_residual_reflayer_const(pv_residual_samp_ctxt, pi2_inp, i4_inp_stride,
                                        ps_ref_mb_mode, &i4_refarray_wd, &s_mb_coord, 1);

    if(ret != OK) return ret;
    /* ---- Interpolation process for Residual prediction     ------ */
    ps_ctxt->pf_interpolate_residual(pv_residual_samp_ctxt, pi2_out, i4_out_stride, i4_refarray_wd,
                                     s_mb_coord.u2_mb_x, s_mb_coord.u2_mb_y, 1);

    /* CR */
    pi2_inp += 1;
    pi2_out += 1;

    /* ------- Constructing refSampleArray ----------------------- */
    ret = isvcd_residual_reflayer_const(pv_residual_samp_ctxt, pi2_inp, i4_inp_stride,
                                        ps_ref_mb_mode, &i4_refarray_wd, &s_mb_coord, 1);

    if(ret != OK) return ret;
    /* ---- Interpolation process for Residual prediction --------- */
    ps_ctxt->pf_interpolate_residual(pv_residual_samp_ctxt, pi2_out, i4_out_stride, i4_refarray_wd,
                                     s_mb_coord.u2_mb_x, s_mb_coord.u2_mb_y, 1);
    return OK;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_samp_mb_dyadic                             */
/*                                                                           */
/*  Description   : MB level function whcih perform the residual resampling  */
/*                  of data of an MB (luma and chroma insclusive)            */
/*                  for Dyadic cases                                         */
/*  Inputs        : pv_residual_samp_ctxt : residual sampling context        */
/*                  ps_ref_luma : reference layer luma data buffer desc      */
/*                  ps_ref_chroma : reference layer chroma data buffer desc  */
/*                  ps_ref_luma_bitmap : ref layer luma bit map buffer desc  */
/*                  ps_ref_chroma_bitmap : ref layer chroma bit map buff des */
/*                  ps_ref_mb_mode : ref layer mb mode map buff desc         */
/*                  ps_out_luma : current layer out luma buffer desc         */
/*                  ps_out_chroma : current layer out chroma buffer desc     */
/*                  ps_mb_coord : current mb coorinate                       */
/*  Globals       : none                                                     */
/*  Processing    : it calls the reference layer construction followed by    */
/*                   interplaotion function for luma and cb and cr           */
/*  Outputs       : inter resampled data of current MB                       */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 06 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_residual_samp_mb_dyadic(void *pv_residual_samp_ctxt, mem_element_t *ps_ref_luma,
                                     mem_element_t *ps_ref_chroma, mem_element_t *ps_ref_mb_mode,
                                     mem_element_t *ps_out_luma, mem_element_t *ps_out_chroma,
                                     UWORD16 u2_mb_x, UWORD16 u2_mb_y)
{
    residual_sampling_ctxt_t *ps_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    /* --------------------------------------------------------------------- */
    /* I/O buffer params                                                     */
    /* --------------------------------------------------------------------- */
    WORD16 *pi2_inp;
    WORD16 *pi2_out;
    WORD32 i4_inp_stride;
    WORD32 i4_out_stride;
    WORD32 i4_luma_nnz;
    WORD32 i4_chroma_nnz;
    WORD32 i4_tx_size;

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_ctxt->i4_res_lyr_id];

    /* --------------------------------------------------------------------- */
    /* LUMA PROCESSING                                                        */
    /* --------------------------------------------------------------------- */
    pi2_inp = (WORD16 *) ps_ref_luma->pv_buffer;
    pi2_out = (WORD16 *) ps_out_luma->pv_buffer;
    i4_inp_stride = ps_ref_luma->i4_num_element_stride;
    i4_out_stride = ps_out_luma->i4_num_element_stride;

    {
        WORD32 i4_offset_x, i4_offset_y;
        residual_samp_map_ctxt_t *ps_luma_map;
        ref_mb_map_t *ps_x_off_len_luma;
        ref_mb_map_t *ps_y_off_len_luma;

        ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
        ps_x_off_len_luma = ps_luma_map->ps_x_offset_length;
        ps_y_off_len_luma = ps_luma_map->ps_y_offset_length;

        /* get the actual offset for the buffers */
        i4_offset_x = ps_x_off_len_luma[u2_mb_x].i2_offset;
        i4_offset_y = ps_y_off_len_luma[u2_mb_y].i2_offset;

        {
            inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms;
            WORD32 i4_mb_x, i4_mb_y;
            UWORD16 u2_luma_mask = 0x0033;
            UWORD8 u1_chrm_mask = 0x11;
            WORD32 i4_luma_rt_sft_amt = 0;
            WORD32 i4_chrm_rt_sft_amt = 0;

            i4_mb_x = ((i4_offset_x + 1) >> MB_WIDTH_SHIFT);
            i4_mb_y = ((i4_offset_y + 1) >> MB_HEIGHT_SHIFT);

            /* get the location of the byte which has the current mb mode */
            ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) ps_ref_mb_mode->pv_buffer;
            ps_inter_lyr_mb_prms += i4_mb_x;
            ps_inter_lyr_mb_prms += i4_mb_y * ps_ref_mb_mode->i4_num_element_stride;

            /* get the approp block in base layer in horz direction */
            if(0 != ((i4_offset_x + 1) & 15))
            {
                u2_luma_mask <<= 2;
                i4_luma_rt_sft_amt += 2;
                u1_chrm_mask <<= 1;
                i4_chrm_rt_sft_amt += 1;
            }
            /* get the approp block in base layer in vert direction */
            if(0 != ((i4_offset_y + 1) & 15))
            {
                u2_luma_mask <<= 8;
                i4_luma_rt_sft_amt += 8;

                u1_chrm_mask <<= 2;
                i4_chrm_rt_sft_amt += 2;
            }

            /* extract the nnz and store it */
            i4_luma_nnz = (ps_inter_lyr_mb_prms->u2_luma_nnz & u2_luma_mask) >> i4_luma_rt_sft_amt;
            i4_chroma_nnz =
                (ps_inter_lyr_mb_prms->u1_chroma_nnz & u1_chrm_mask) >> i4_chrm_rt_sft_amt;
            i4_tx_size =
                (ps_inter_lyr_mb_prms->i1_tx_size < 0) ? 1 : ps_inter_lyr_mb_prms->i1_tx_size;
        }

        /* since in dyadic case the window width and height will be 10x10   */
        /* and the window start offsets will be always 1 column left and    */
        /* 1 row above the block boundary. so the pointer and the required  */
        /* positions are appropriately modified                             */
        if(i4_offset_x >= 0)
        {
            pi2_inp++;
        }

        if(i4_offset_y >= 0)
        {
            pi2_inp += i4_inp_stride;
        }

        ps_ctxt->pf_residual_luma_dyadic(pv_residual_samp_ctxt, pi2_inp, i4_inp_stride, pi2_out,
                                         i4_out_stride, ps_ref_mb_mode, u2_mb_x, u2_mb_y,
                                         i4_luma_nnz, i4_tx_size);
    }

    /* --------------------------------------------------------------------- */
    /* CHROMA PROCESSING                                                       */
    /* --------------------------------------------------------------------- */
    /* CB */
    pi2_inp = (WORD16 *) ps_ref_chroma->pv_buffer;
    pi2_out = (WORD16 *) ps_out_chroma->pv_buffer;
    i4_inp_stride = ps_ref_chroma->i4_num_element_stride;
    i4_out_stride = ps_out_chroma->i4_num_element_stride;

    /* choose the appropriate chroma processing routine */
    if(SVCD_FALSE == ps_lyr_ctxt->i4_chrm_alt_proc)
    {
        WORD32 i4_offset_x, i4_offset_y;
        residual_samp_map_ctxt_t *ps_chroma_map;
        ref_mb_map_t *ps_x_off_len_chroma;
        ref_mb_map_t *ps_y_off_len_chroma;

        ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;
        ps_x_off_len_chroma = ps_chroma_map->ps_x_offset_length;
        ps_y_off_len_chroma = ps_chroma_map->ps_y_offset_length;

        /* get the actual offset for the buffers */
        i4_offset_x = ps_x_off_len_chroma[u2_mb_x].i2_offset;
        i4_offset_y = ps_y_off_len_chroma[u2_mb_y].i2_offset;

        /* since in dyadic case the window width and height will be 6x6     */
        /* and the window start offsets will be always 1 column left and    */
        /* 1 row above the block boundary. so the pointer and the required  */
        /* positions are appropriately modified                             */
        if(i4_offset_x >= 0)
        {
            pi2_inp += 2;
        }

        if(i4_offset_y >= 0)
        {
            pi2_inp += i4_inp_stride;
        }

        if(0 != (i4_chroma_nnz & 0x01))
        {
            ps_ctxt->pf_residual_chroma_dyadic(pv_residual_samp_ctxt, pi2_inp, i4_inp_stride,
                                               pi2_out, i4_out_stride);
        }
    }
    else
    {
        ps_ctxt->pf_residual_chroma_dyadic_alt(pv_residual_samp_ctxt, u2_mb_x, u2_mb_y,
                                               ps_ref_mb_mode, pi2_inp, i4_inp_stride, pi2_out,
                                               i4_out_stride, SVCD_FALSE);
    }

    /* CR */
    pi2_inp += 1;
    pi2_out += 1;

    if(SVCD_FALSE == ps_lyr_ctxt->i4_chrm_alt_proc)
    {
        if(0 != (i4_chroma_nnz & 0x10))
        {
            ps_ctxt->pf_residual_chroma_dyadic(pv_residual_samp_ctxt, pi2_inp, i4_inp_stride,
                                               pi2_out, i4_out_stride);
        }
    }
    else
    {
        ps_ctxt->pf_residual_chroma_dyadic_alt(pv_residual_samp_ctxt, u2_mb_x, u2_mb_y,
                                               ps_ref_mb_mode, pi2_inp, i4_inp_stride, pi2_out,
                                               i4_out_stride, SVCD_TRUE);
    }
    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_samp_populate_list                         */
/*                                                                           */
/*  Description   : This is a seq or frame level init function which fills   */
/*                  all offsets, projected locations arrays based on         */
/*                  the two resolutions  and cropping parameters             */
/*  Inputs        : refer to doxygen comments below                          */
/*  Globals       : none                                                     */
/*  Processing    : it projects the locations and computes the values        */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 07 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_residual_samp_populate_list(residual_samp_map_ctxt_t *ps_map_ctxt,
                                       dec_seq_params_t *ps_sps,
                                       dec_svc_seq_params_t *ps_subset_sps,
                                       res_prms_t *ps_curr_res_prms, res_prms_t *ps_ref_res_prms,
                                       WORD32 i4_chroma_flag)
{
    /* --------------------------------------------------------------------- */
    /* Local variables required for finding the mapping between the layers     */
    /* --------------------------------------------------------------------- */
    UWORD32 u4_shift_x;
    UWORD32 u4_shift_y;
    UWORD32 u4_scale_x;
    UWORD32 u4_scale_y;
    WORD32 i4_offset_x;
    WORD32 i4_offset_y;
    WORD32 i4_add_x;
    WORD32 i4_add_y;
    WORD32 i4_delta_x;
    WORD32 i4_delta_y;
    WORD32 i4_refphase_x;
    WORD32 i4_refphase_y;
    WORD32 i4_phase_x;
    WORD32 i4_phase_y;
    WORD32 i4_sub_wd;
    WORD32 i4_sub_ht;
    WORD32 i4_mb_wd;
    WORD32 i4_mb_ht;
    /* --------------------------------------------------------------------- */
    /* Local Pointer Declaration for arrays in Mapping context                 */
    /* --------------------------------------------------------------------- */
    ref_mb_map_t *ps_x_off_len;
    ref_mb_map_t *ps_y_off_len;
    UWORD32 i4_ref_wd;
    UWORD32 i4_ref_ht;
    UWORD32 i4_scaled_wd;
    UWORD32 i4_scaled_ht;
    WORD32 i4_curr_lyr_width;
    WORD32 i4_curr_lyr_height;

    /* --------------------------------------------------------------------- */
    /* Local Flag Declaration                                                 */
    /* --------------------------------------------------------------------- */
    WORD32 i4_ref_layer_field_pic_flag;
    WORD32 i4_field_pic_flag;
    WORD32 i4_frame_mbs_only_flag;
    WORD32 i4_ref_layer_frame_Mbs_only_flag;
    WORD32 i4_field_Mb_flag;
    WORD32 i4_bot_field_flag;

    /* --------------------------------------------------------------------- */
    /* Cropping Parameters Declaration                                         */
    /* --------------------------------------------------------------------- */
    WORD32 i4_scaled_ref_layer_left_offset;
    WORD32 i4_scaled_ref_layer_top_offset;

    /* --------------------------------------------------------------------- */
    /* Hardcoding flag information    (assuming no field support) */
    /* --------------------------------------------------------------------- */
    i4_ref_layer_field_pic_flag = SVCD_FALSE;
    i4_field_pic_flag = SVCD_FALSE;
    i4_frame_mbs_only_flag = SVCD_TRUE;
    i4_field_Mb_flag = SVCD_FALSE;
    i4_bot_field_flag = SVCD_FALSE;
    i4_ref_layer_frame_Mbs_only_flag = SVCD_TRUE;

    /* --------------------------------------------------------------------- */
    /* Pointer and Paramater are intialized    - Chroma and Luma */
    /* --------------------------------------------------------------------- */
    {
        WORD32 i4_base_width;
        WORD32 i4_base_height;
        WORD32 i4_ref_layer_chroma_phase_x_plus1_flag;
        WORD32 i4_ref_layer_chroma_phase_y_plus1;
        WORD32 i4_chroma_phase_x_plus1_flag;
        WORD32 i4_chroma_phase_y_plus1;

        /* ------------------------------------------------------------- */
        /* HARD CODED FOR 420                                             */
        /* ------------------------------------------------------------- */
        WORD32 i4_sub_wd_chroma = 2;
        WORD32 i4_sub_ht_chroma = 2;

        i4_base_width = ps_ref_res_prms->i4_res_width;
        i4_base_height = ps_ref_res_prms->i4_res_height;

        i4_ref_layer_chroma_phase_x_plus1_flag =
            ps_curr_res_prms->i1_ref_lyr_chroma_phase_x_plus1_flag;
        i4_ref_layer_chroma_phase_y_plus1 = ps_curr_res_prms->i1_ref_lyr_chroma_phase_y_plus1;
        i4_chroma_phase_x_plus1_flag =
            ps_subset_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_x_plus1_flag;
        i4_chroma_phase_y_plus1 =
            ps_subset_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_y_plus1;
        i4_scaled_ref_layer_left_offset = ps_curr_res_prms->s_ref_lyr_scaled_offset.i2_left;
        i4_scaled_ref_layer_top_offset = ps_curr_res_prms->s_ref_lyr_scaled_offset.i2_top;

        /* ----------------------------------------------------------------- */
        /* Computing Effective Frame Dimensions                                 */
        /* ------------------------------------------------------------------*/
        i4_ref_wd = (i4_base_width >> i4_chroma_flag);
        i4_ref_ht = (i4_base_height >> i4_chroma_flag) * (1 + i4_ref_layer_field_pic_flag);
        i4_scaled_wd = ps_curr_res_prms->u2_scaled_ref_width;
        i4_scaled_ht = ps_curr_res_prms->u2_scaled_ref_height;
        i4_scaled_wd = (i4_scaled_wd >> i4_chroma_flag);
        i4_scaled_ht = (i4_scaled_ht >> i4_chroma_flag) * (1 + i4_field_pic_flag);

        if(1 == i4_chroma_flag)
        {
            i4_refphase_x = i4_ref_layer_chroma_phase_x_plus1_flag - 1;
            i4_refphase_y = i4_ref_layer_chroma_phase_y_plus1 - 1;
            i4_phase_x = i4_chroma_phase_x_plus1_flag - 1;
            i4_phase_y = i4_chroma_phase_y_plus1 - 1;
            i4_sub_wd = i4_sub_wd_chroma;
            i4_sub_ht = i4_sub_ht_chroma;
            i4_mb_wd = MB_WIDTH >> 1;
            i4_mb_ht = MB_HEIGHT >> 1;
        }
        else
        {
            i4_refphase_x = 0;
            i4_refphase_y = 0;
            i4_phase_x = 0;
            i4_phase_y = 0;
            i4_sub_wd = 1;
            i4_sub_ht = 1;
            i4_mb_wd = MB_WIDTH;
            i4_mb_ht = MB_HEIGHT;
        }
    }

    /* --------------------------------------------------------------------- */
    /* Derive shift x and y based on level idd                               */
    /* --------------------------------------------------------------------- */
    if(ps_sps->u1_level_idc <= 30)
    {
        u4_shift_x = 16;
        u4_shift_y = 16;
    }
    else
    {
        u4_shift_x = 31 - isvcd_get_ceil_log2(i4_ref_wd);
        u4_shift_y = 31 - isvcd_get_ceil_log2(i4_ref_ht);
    }

    /* --------------------------------------------------------------------- */
    /* The following condition is not true in our case for time being         */
    /* --------------------------------------------------------------------- */
    if((SVCD_FALSE == i4_frame_mbs_only_flag) || (SVCD_FALSE == i4_ref_layer_frame_Mbs_only_flag))
    {
        i4_phase_y = i4_phase_y + 4 * i4_bot_field_flag;

        if(1 == i4_ref_layer_frame_Mbs_only_flag)
            i4_refphase_y = (2 * i4_refphase_y) + 2;
        else
            i4_refphase_y = i4_refphase_y + (4 * i4_bot_field_flag);
    }

    /* --------------------------------------------------------------------- */
    /* Dx and Dy Computation - Ratio of the base and enhance layer width     */
    /* --------------------------------------------------------------------- */
    u4_scale_x = ((i4_ref_wd << u4_shift_x) + (i4_scaled_wd >> 1)) / (i4_scaled_wd);

    u4_scale_y = ((i4_ref_ht << u4_shift_y) + (i4_scaled_ht >> 1)) / (i4_scaled_ht);

    i4_offset_x = i4_scaled_ref_layer_left_offset / i4_sub_wd;
    i4_add_x = (((i4_ref_wd * (2 + i4_phase_x)) << (u4_shift_x - 2)) + (i4_scaled_wd >> 1)) /
                   i4_scaled_wd +
               (1 << (u4_shift_x - 5));
    i4_delta_x = 4 * (2 + i4_refphase_x);

    if((SVCD_TRUE == i4_frame_mbs_only_flag) && (SVCD_TRUE == i4_ref_layer_frame_Mbs_only_flag))
    {
        i4_offset_y = i4_scaled_ref_layer_top_offset / i4_sub_ht;
        i4_add_y = (((i4_ref_ht * (2 + i4_phase_y)) << (u4_shift_y - 2)) + (i4_scaled_ht >> 1)) /
                       i4_scaled_ht +
                   (1 << (u4_shift_y - 5));
        i4_delta_y = 4 * (2 + i4_refphase_y);
    }
    else
    {
        i4_offset_y = i4_scaled_ref_layer_top_offset / (2 * i4_sub_ht);
        i4_add_y = (((i4_ref_ht * (2 + i4_phase_y)) << (u4_shift_y - 3)) + (i4_scaled_ht >> 1)) /
                       i4_scaled_ht +
                   (1 << (u4_shift_y - 5));
        i4_delta_y = 2 * (2 + i4_refphase_y);
    }

    /* --------------------------------------------------------------------- */
    /* Intializing Local Pointers    - Chroma and Luma                       */
    /* --------------------------------------------------------------------- */
    ps_x_off_len = ps_map_ctxt->ps_x_offset_length;
    ps_y_off_len = ps_map_ctxt->ps_y_offset_length;
    i4_curr_lyr_width = ps_curr_res_prms->i4_res_width >> i4_chroma_flag;
    i4_curr_lyr_height = ps_curr_res_prms->i4_res_height >> i4_chroma_flag;

    {
        WORD32 i4_i, i4_j;

        /* ----------------------------------------------------------------- */
        /* Computation of offsetX refArrayW Xmin and Xmax Lists               */
        /* ----------------------------------------------------------------- */
        for(i4_i = 0; i4_i < i4_curr_lyr_width; i4_i = i4_i + i4_mb_wd)
        {
            WORD32 i4_x_refmin16;
            WORD32 i4_x_refmax16;
            WORD32 i4_x_offset;

            i4_x_refmin16 = (WORD64) (((WORD64) ((i4_i - i4_offset_x) * u4_scale_x) + i4_add_x) >>
                                      ((WORD32) (u4_shift_x - 4))) -
                            i4_delta_x;

            i4_x_refmax16 =
                (WORD64) (((WORD64) (i4_i + i4_mb_wd - 1 - i4_offset_x) * u4_scale_x + i4_add_x) >>
                          ((WORD32) (u4_shift_x - 4))) -
                i4_delta_x;

            /* AC205 */
            i4_x_offset = i4_x_refmin16 >> 4;
            ps_x_off_len->i2_offset = i4_x_offset;
            ps_x_off_len->i2_length = (i4_x_refmax16 >> 4) - i4_x_offset + 2;

            /* increment the pointer */
            ps_x_off_len++;

        } /* end of loop over current layer width */

        /* ----------------------------------------------------------------- */
        /* Computation of offsetY refArrayH Ymin and Ymax Lists              */
        /* ----------------------------------------------------------------- */
        for(i4_j = 0; i4_j < i4_curr_lyr_height; i4_j = i4_j + i4_mb_ht)
        {
            WORD32 i4_y_refmin16;
            WORD32 i4_y_refmax16;
            WORD32 i4_y_offset;

            i4_y_refmin16 = (WORD64) (((WORD64) (i4_j - i4_offset_y) * u4_scale_y + i4_add_y) >>
                                      ((WORD32) (u4_shift_y - 4))) -
                            i4_delta_y;

            i4_y_refmax16 =
                (WORD64) (((WORD64) (i4_j + i4_mb_ht - 1 - i4_offset_y) * u4_scale_y + i4_add_y) >>
                          ((WORD32) (u4_shift_y - 4))) -
                i4_delta_y;

            /* AC205 */
            i4_y_offset = i4_y_refmin16 >> 4;
            ps_y_off_len->i2_offset = i4_y_offset;
            ps_y_off_len->i2_length = (i4_y_refmax16 >> 4) - i4_y_offset + 2;

            /* increment the pointer */
            ps_y_off_len++;

        } /* end of loop over current layer height */
    }

    /* --------------------------------------------------------------------- */
    /* Computation of Xref and Xphase List as per standard                     */
    /* --------------------------------------------------------------------- */
    ps_x_off_len = ps_map_ctxt->ps_x_offset_length;
    ps_y_off_len = ps_map_ctxt->ps_y_offset_length;

    {
        WORD32 i4_xc;
        WORD32 i4_offset_x_index;
        ref_pixel_map_t *ps_x_pos_phase;

        ps_x_pos_phase = ps_map_ctxt->ps_x_pos_phase;

        for(i4_xc = 0; i4_xc < i4_curr_lyr_width; i4_xc++)
        {
            WORD32 i4_x_offset;
            WORD32 i4_x_ref16;

            i4_offset_x_index = i4_xc / i4_mb_wd;
            i4_x_offset = ps_x_off_len[i4_offset_x_index].i2_offset;
            i4_x_ref16 = (WORD64) (((WORD64) (i4_xc - i4_offset_x) * u4_scale_x + i4_add_x) >>
                                   ((WORD32) (u4_shift_x - 4))) -
                         i4_delta_x;

            /* store the values */
            ps_x_pos_phase->i2_ref_pos = (i4_x_ref16 >> 4) - i4_x_offset;
            ps_x_pos_phase->i2_phase = (i4_x_ref16 - (16 * i4_x_offset)) & 15;

            /* increment the pointer */
            ps_x_pos_phase++;
        } /* end of loop over scaled width */
    }

    /* --------------------------------------------------------------------- */
    /* Computation of Yref and Yphase List as per standard                     */
    /* --------------------------------------------------------------------- */
    {
        WORD32 i4_yc;
        WORD32 i4_offset_y_index;
        ref_pixel_map_t *ps_y_pos_phase;

        ps_y_pos_phase = ps_map_ctxt->ps_y_pos_phase;

        for(i4_yc = 0; i4_yc < i4_curr_lyr_height; i4_yc++)
        {
            WORD32 i4_y_offset;
            WORD32 i4_y_ref16;

            i4_offset_y_index = i4_yc / i4_mb_ht;
            i4_y_offset = ps_y_off_len[i4_offset_y_index].i2_offset;

            if((SVCD_FALSE == i4_frame_mbs_only_flag) ||
               (SVCD_FALSE == i4_ref_layer_frame_Mbs_only_flag))
            {
                i4_yc = i4_yc >> (1 - i4_field_Mb_flag);
            }

            i4_y_ref16 = (WORD64) ((((WORD64) (i4_yc - i4_offset_y) * u4_scale_y + i4_add_y) >>
                                    ((WORD32) (u4_shift_y - 4))) -
                                   i4_delta_y);
            ps_y_pos_phase->i2_ref_pos = (i4_y_ref16 >> 4) - i4_y_offset;
            ps_y_pos_phase->i2_phase = (i4_y_ref16 - (16 * i4_y_offset)) & 15;

            /* increment the pointer */
            ps_y_pos_phase++;
        } /* end of loop over scaled height */
    }
    return;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_samp_res_init                              */
/*                                                                           */
/*  Description   : this function calculates the scale factors and initialise*/
/*                  the context structure                                    */
/*                                                                           */
/*  Inputs        : pv_residual_samp_ctxt: handle to private structure       */
/*                  ps_curr_lyr_res_prms: pointer to current resolution      */
/*                                               params                      */
/*  Globals       : none                                                     */
/*  Processing    : it stores the layer dimensions                           */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 06 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_residual_samp_res_init(void *pv_svc_dec)
{
    residual_sampling_ctxt_t *ps_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    dec_seq_params_t *ps_sps;
    dec_svc_seq_params_t *ps_subset_sps;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) pv_svc_dec;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;

    void *pv_residual_samp_ctxt = ps_svc_lyr_dec->pv_residual_sample_ctxt;
    res_prms_t *ps_curr_lyr_res_prms = &ps_svc_lyr_dec->s_res_prms;
    ref_mb_map_t **pps_luma_map_horz = &ps_svc_lyr_dec->ps_ressam_luma_map_horz;
    ref_mb_map_t **pps_chroma_map_horz = &ps_svc_lyr_dec->ps_ressam_chroma_map_horz;
    ref_mb_map_t **pps_luma_map_vert = &ps_svc_lyr_dec->ps_ressam_luma_map_vert;
    ref_mb_map_t **pps_chroma_map_vert = &ps_svc_lyr_dec->ps_ressam_chroma_map_vert;

    if((NULL == pv_residual_samp_ctxt) || (NULL == ps_curr_lyr_res_prms) ||
       (NULL == pps_luma_map_horz) || (NULL == pps_chroma_map_horz) ||
       (NULL == pps_luma_map_vert) || (NULL == pps_chroma_map_vert))
    {
        return NOT_OK;
    }

    ps_ctxt = (residual_sampling_ctxt_t *) pv_residual_samp_ctxt;

    /* if called for base resolution store deafult values */
    if(SVCD_TRUE == ps_svc_lyr_dec->u1_base_res_flag)
    {
        *pps_luma_map_horz = NULL;
        *pps_chroma_map_horz = NULL;
        *pps_luma_map_vert = NULL;
        *pps_chroma_map_vert = NULL;
        ps_ctxt->i4_res_lyr_id = -1;
        ps_ctxt->i4_ref_width = ps_curr_lyr_res_prms->i4_res_width;
        ps_ctxt->i4_ref_height = ps_curr_lyr_res_prms->i4_res_height;
        return OK;
    }

    /* derive the current sps */
    ps_sps = ps_dec->ps_cur_sps;
    ps_subset_sps = ps_svc_lyr_dec->ps_cur_subset_sps;

    /* store the res id appropriately */
    ps_ctxt->i4_res_lyr_id = ps_svc_lyr_dec->u1_layer_id - 1;

    /* get the current layer ctxt */
    ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[ps_svc_lyr_dec->u1_layer_id - 1];

    /* get the width and heights */
    ps_lyr_ctxt->i4_curr_width = ps_curr_lyr_res_prms->i4_res_width;
    ps_lyr_ctxt->i4_curr_height = ps_curr_lyr_res_prms->i4_res_height;
    ps_lyr_ctxt->i4_ref_width = ps_ctxt->i4_ref_width;
    ps_lyr_ctxt->i4_ref_height = ps_ctxt->i4_ref_height;

    /* store the strcuture pointer containing projected locations */
    *pps_luma_map_horz = ps_lyr_ctxt->s_luma_map_ctxt.ps_x_offset_length;
    *pps_chroma_map_horz = ps_lyr_ctxt->s_chroma_map_ctxt.ps_x_offset_length;
    *pps_luma_map_vert = ps_lyr_ctxt->s_luma_map_ctxt.ps_y_offset_length;
    *pps_chroma_map_vert = ps_lyr_ctxt->s_chroma_map_ctxt.ps_y_offset_length;

    /* check for recomputation of mapping required */
    if(SVCD_TRUE == ps_curr_lyr_res_prms->u1_remap_req_flag)
    {
        res_prms_t s_ref_res_prms = {0};

        /* store the reference layer resolution width and height */
        s_ref_res_prms.i4_res_width = ps_ctxt->i4_ref_width;
        s_ref_res_prms.i4_res_height = ps_ctxt->i4_ref_height;

        /* call the frame level projections calculation function */
        isvcd_residual_samp_populate_list(&ps_lyr_ctxt->s_luma_map_ctxt, ps_sps, ps_subset_sps,
                                          ps_curr_lyr_res_prms, &s_ref_res_prms, 0);
        isvcd_residual_samp_populate_list(&ps_lyr_ctxt->s_chroma_map_ctxt, ps_sps, ps_subset_sps,
                                          ps_curr_lyr_res_prms, &s_ref_res_prms, 1);

        /* default values for flags */
        ps_lyr_ctxt->pf_residual_samp_mb = &isvcd_residual_samp_mb;
        ps_lyr_ctxt->i4_chrm_horz_int_mode = 0;
        ps_lyr_ctxt->i4_chrm_vert_int_mode = 0;
        ps_lyr_ctxt->i4_chrm_alt_proc = SVCD_FALSE;

        /* Store the Dyadic flag */
        ps_lyr_ctxt->i4_dyadic_flag = ps_curr_lyr_res_prms->u1_dyadic_flag;

        /* set the appropriate chroma processing routine based on */
        /* phase values */
        if(SVCD_TRUE == ps_curr_lyr_res_prms->u1_dyadic_flag)
        {
            WORD32 i4_ref_layer_chroma_phase_x_plus1_flag;
            WORD32 i4_ref_layer_chroma_phase_y_plus1;
            WORD32 i4_chroma_phase_x_plus1_flag;
            WORD32 i4_chroma_phase_y_plus1;

            ps_lyr_ctxt->pf_residual_samp_mb = &isvcd_residual_samp_mb_dyadic;
            i4_ref_layer_chroma_phase_x_plus1_flag =
                ps_curr_lyr_res_prms->i1_ref_lyr_chroma_phase_x_plus1_flag;
            i4_ref_layer_chroma_phase_y_plus1 =
                ps_curr_lyr_res_prms->i1_ref_lyr_chroma_phase_y_plus1;
            i4_chroma_phase_x_plus1_flag =
                ps_subset_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_x_plus1_flag;
            i4_chroma_phase_y_plus1 =
                ps_subset_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_y_plus1;
            if((0 == i4_ref_layer_chroma_phase_x_plus1_flag) && (1 == i4_chroma_phase_x_plus1_flag))
            {
                ps_lyr_ctxt->i4_chrm_horz_int_mode = 1;
                ps_lyr_ctxt->i4_chrm_alt_proc = SVCD_TRUE;
            }

            if((0 == i4_ref_layer_chroma_phase_y_plus1) && (1 == i4_chroma_phase_y_plus1))
            {
                ps_lyr_ctxt->i4_chrm_vert_int_mode = 1;
                ps_lyr_ctxt->i4_chrm_alt_proc = SVCD_TRUE;
            }

            if((0 == i4_ref_layer_chroma_phase_y_plus1) && (2 == i4_chroma_phase_y_plus1))
            {
                ps_lyr_ctxt->i4_chrm_vert_int_mode = 1;
                ps_lyr_ctxt->i4_chrm_alt_proc = SVCD_TRUE;
            }

            if((2 == i4_ref_layer_chroma_phase_y_plus1) && (0 == i4_chroma_phase_y_plus1))
            {
                ps_lyr_ctxt->i4_chrm_vert_int_mode = 2;
                ps_lyr_ctxt->i4_chrm_alt_proc = SVCD_TRUE;
            }
        }
    }
    else
    {
        /* should take false value */
        if(SVCD_FALSE != ps_curr_lyr_res_prms->u1_remap_req_flag)
        {
            return NOT_OK;
        }
    }

    /* store the current layer width and height to context */
    ps_ctxt->i4_ref_width = ps_curr_lyr_res_prms->i4_res_width;
    ps_ctxt->i4_ref_height = ps_curr_lyr_res_prms->i4_res_height;

    /* assert on max ranges of width and shift values */
    if((ps_lyr_ctxt->i4_curr_width > H264_MAX_FRAME_WIDTH) ||
       (ps_lyr_ctxt->i4_ref_width > H264_MAX_FRAME_WIDTH) ||
       (ps_lyr_ctxt->i4_curr_height > H264_MAX_FRAME_HEIGHT) ||
       (ps_lyr_ctxt->i4_ref_height > H264_MAX_FRAME_HEIGHT))
    {
        return NOT_OK;
    }
    return OK;
}