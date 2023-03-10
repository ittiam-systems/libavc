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
 *  isvce_core_coding.c
 *
 * @brief
 *  This file contains routines that perform luma and chroma core coding for
 *  intra macroblocks
 *
 * @author
 *  ittiam
 *
 * @par List of Functions:
 *  - isvce_pack_l_mb_i16()
 *  - isvce_pack_c_mb_i8()
 *  - isvce_code_luma_intra_macroblock_16x16()
 *  - isvce_code_luma_intra_macroblock_4x4()
 *  - isvce_code_chroma_intra_macroblock_8x8()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* User include files */
#include "ih264_typedefs.h"
#include "ih264_debug.h"
#include "ih264_platform_macros.h"
#include "iv2.h"
#include "ive2.h"
#include "isvc_macros.h"
#include "isvc_defs.h"
#include "ih264e_config.h"
#include "isvce_defs.h"
#include "ih264_trans_data.h"
#include "ih264e_error.h"
#include "ih264e_bitstream.h"
#include "ime_distortion_metrics.h"
#include "ime_defs.h"
#include "ime_structs.h"
#include "isvc_structs.h"
#include "isvc_trans_quant_itrans_iquant.h"
#include "isvc_inter_pred_filters.h"
#include "isvc_mem_fns.h"
#include "ih264_padding.h"
#include "ih264_intra_pred_filters.h"
#include "ih264_deblk_edge_filters.h"
#include "isvc_cabac_tables.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "isvce_rate_control.h"
#include "isvce_cabac_structs.h"
#include "isvce_structs.h"
#include "isvce_globals.h"
#include "isvce_core_coding.h"
#include "isvce_mc.h"
#include "isvce_ibl_eval.h"

/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

/**
*******************************************************************************
*
* @brief
*  This function performs does the DCT transform then Hadamard transform
*  and quantization for a macroblock when the mb mode is intra 16x16 mode
*
* @par Description:
*  First  cf4 is done on all 16 4x4 blocks of the 16x16 input block.
*  Then hadamard transform is done on the DC coefficients
*  Quantization is then performed on the 16x16 block, 4x4 wise
*
* @param[in] pu1_src
*  Pointer to source sub-block
*
* @param[in] pu1_pred
*  Pointer to prediction sub-block
*
* @param[in] pi2_out
*  Pointer to residual sub-block
*  The output will be in linear format
*  The first 16 continuous locations will contain the values of Dc block
*  After DC block and a stride 1st AC block will follow
*  After one more stride next AC block will follow
*  The blocks will be in raster scan order
*
* @param[in] i4_src_stride
*  Source stride
*
* @param[in] i4_pred_stride
*  Prediction stride
*
* @param[in] dst_strd
*  Destination stride
*
* @param[in] pu2_scale_matrix
*  The quantization matrix for 4x4 transform
*
* @param[in] pu2_threshold_matrix
*  Threshold matrix
*
* @param[in] u4_qbits
*  15+QP/6
*
* @param[in] u4_round_factor
*  Round factor for quant
*
* @param[out] pu1_nnz
*  Memory to store the non-zeros after transform
*  The first byte will be the nnz of DC block
*  From the next byte the AC nnzs will be stored in raster scan order
*
* @param u4_dc_flag
*  Signals if Dc transform is to be done or not
*   1 -> Dc transform will be done
*   0 -> Dc transform will not be done
*
* @remarks
*
*******************************************************************************
*/
void isvce_luma_16x16_resi_trans_dctrans_quant(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_quant_coeffs,
    buffer_container_t *ps_upsampled_res, isa_dependent_fxns_t *ps_isa_dependent_fxns,
    const UWORD16 *pu2_scale_matrix, const UWORD16 *pu2_threshold_matrix, UWORD8 *pu1_nnz,
    UWORD32 u4_qbits, UWORD32 u4_round_factor, UWORD32 u4_dc_flag, UWORD8 u1_use_upsampled_res)
{
    WORD32 blk_cntr;
    WORD32 i4_offsetx, i4_offsety;

    enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;
    buffer_container_t s_src = ps_src[0];
    buffer_container_t s_pred = ps_pred[0];
    buffer_container_t s_quant_coeffs = ps_quant_coeffs[0];
    buffer_container_t s_upsampled_res = {0};
    resi_trans_quant_constants_t s_resi_trans_quant_constants = {
        .pu2_scale_matrix = pu2_scale_matrix,
        .pu2_threshold_matrix = pu2_threshold_matrix,
        .u4_qbits = u4_qbits,
        .u4_round_factor = u4_round_factor};

    UWORD8 u1_resi_trans_fxn_idx = isvc_get_resi_trans_quant_variant_idx(u1_use_upsampled_res);

    /* Move to the ac addresses */
    pu1_nnz++;
    s_quant_coeffs.pv_data = ((WORD16 *) s_quant_coeffs.pv_data) + s_quant_coeffs.i4_data_stride;

    if(u1_use_upsampled_res)
    {
        s_upsampled_res = ps_upsampled_res[0];
    }

    for(blk_cntr = 0; blk_cntr < NUM_LUMA4x4_BLOCKS_IN_MB; blk_cntr++)
    {
        IND2SUB_LUMA_MB(blk_cntr, i4_offsetx, i4_offsety);

        s_src.pv_data =
            ((UWORD8 *) ps_src[0].pv_data) + i4_offsetx + i4_offsety * ps_src[0].i4_data_stride;
        s_pred.pv_data =
            ((UWORD8 *) ps_pred[0].pv_data) + i4_offsetx + i4_offsety * ps_pred[0].i4_data_stride;
        s_quant_coeffs.pv_data =
            ((WORD16 *) ps_quant_coeffs[0].pv_data) + blk_cntr * ps_quant_coeffs[0].i4_data_stride;

        if(u1_use_upsampled_res)
        {
            s_upsampled_res.pv_data = ((WORD16 *) ps_upsampled_res[0].pv_data) + i4_offsetx +
                                      i4_offsety * ps_upsampled_res[0].i4_data_stride;
        }

        /* Move to the ac addresses */
        s_quant_coeffs.pv_data =
            ((WORD16 *) s_quant_coeffs.pv_data) + ps_quant_coeffs[0].i4_data_stride;

        s_quant_coeffs.i4_data_stride = 4;

        ps_enc_loop_fxns->apf_resi_trans_quant_4x4[u1_resi_trans_fxn_idx](
            &s_src, &s_pred, &s_quant_coeffs, &s_upsampled_res, &s_resi_trans_quant_constants,
            &pu1_nnz[blk_cntr], ((WORD16 *) ps_quant_coeffs->pv_data) + blk_cntr,
            u1_use_upsampled_res);
    }

    if(!u4_dc_flag)
    {
        return;
    }

    /*
     * In case of i16x16, we need to remove the contribution of dc coeffs into
     * nnz of each block. We are doing that in the packing function
     */

    /* Adjust pointers to point to dc values */
    s_quant_coeffs = ps_quant_coeffs[0];
    pu1_nnz--;

    u4_qbits++;
    u4_round_factor <<= 1;

    ps_enc_loop_fxns->pf_hadamard_quant_4x4(((WORD16 *) s_quant_coeffs.pv_data),
                                            ((WORD16 *) s_quant_coeffs.pv_data),
                                            &s_resi_trans_quant_constants, &pu1_nnz[0]);
}

/**
*******************************************************************************
*
* @brief
*  This function performs the intra 16x16 inverse transform process for H264
*  it includes inverse Dc transform, inverse quant and then inverse transform
*
* @par Description:
*
* @param[in] pi2_src
*  Input data, 16x16 size
*  First 16 mem locations will have the Dc coffs in rater scan order in linear
*fashion after a stride 1st AC clock will be present again in raster can order
*  Then each AC block of the 16x16 block will follow in raster scan order
*
* @param[in] pu1_pred
*  The predicted data, 16x16 size
*  Block by block form
*
* @param[in] pu1_out
*  Output 16x16
*  In block by block form
*
* @param[in] i4_src_stride
*  Source stride
*
* @param[in] i4_pred_stride
*  input stride for prediction buffer
*
* @param[in] i4_out_stride
*  input stride for output buffer
*
* @param[in] pu2_iscale_mat
*  Inverse quantization matrix for 4x4 transform
*
* @param[in] pu2_weigh_mat
*  weight matrix of 4x4 transform
*
* @param[in] u4_qp_div_6
*  QP/6
*
* @param[in] pi4_tmp
*  Input temporary buffer
*  needs to be at least 20 in size
*
* @param[in] pu4_cntrl
*  Controls the transform path
*  total Last 17 bits are used
*  the 16th th bit will correspond to DC block
*  and 32-17 will correspond to the ac blocks in raster scan order
*  bit equaling zero indicates that the entire 4x4 block is zero for DC
*  For AC blocks a bit equaling zero will mean that all 15 AC coffs of the block
*is nonzero
*
* @param[in] pi4_tmp
*  Input temporary buffer
*  needs to be at least COFF_CNT_SUB_BLK_4x4+COFF_CNT_SUB_BLK_4x4 size
*
* @returns
*  none
*
* @remarks
*  The all zero case must be taken care outside
*
*******************************************************************************
*/
void isvce_luma_16x16_idctrans_iquant_itrans_recon(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_recon,
    buffer_container_t *ps_res, buffer_container_t *ps_res_pred,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants,
    isa_dependent_fxns_t *ps_isa_dependent_fxns, WORD32 *pi4_tmp, UWORD32 u4_cntrl,
    UWORD32 u4_dc_trans_flag, UWORD8 u1_res_accumulate)
{
    /* Cntrl bits for 4x4 transforms
     * u4_blk_cntrl       : controls if a 4x4 block should be processed in ac path
     * u4_dc_cntrl        : controls is a 4x4 block is to be processed in dc path
     *                    : dc block must contain only single dc coefficient
     * u4_empty_blk_cntrl : control fot 4x4 block with no coeffs, ie no dc and ac
     *                    : ie not (ac or dc)
     */
    UWORD32 u4_blk_cntrl, u4_dc_cntrl, u4_empty_blk_cntrl;
    UWORD32 u4_blk_id;
    WORD32 i4_offset_x, i4_offset_y;
    UWORD32 u4_dc_inc;
    WORD16 *pi2_dc_src;

    enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;
    buffer_container_t s_src = ps_src[0];
    buffer_container_t s_pred = ps_pred[0];
    buffer_container_t s_recon = ps_recon[0];
    buffer_container_t s_res = ps_res[0];
    buffer_container_t s_res_pred = ps_res_pred[0];

    /* Start index for inverse quant in a 4x4 block */
    WORD32 i4_iq_start_idx = (u4_dc_trans_flag == 0) ? 0 : 1;
    const UWORD16 *pu2_iscale_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    UWORD8 u1_iq_it_recon_fxn_idx =
        isvc_get_iq_it_recon_variant_idx(!!u4_dc_trans_flag, u1_res_accumulate);

    /*
     * For intra blocks we need to do inverse dc transform
     * In case if intra blocks, its here that we populate the dc bits in cntrl
     * as they cannot be populated any earlier
     */
    if(u4_dc_trans_flag)
    {
        UWORD32 cntr, u4_dc_cntrl;

        /* Do inv hadamard and place the results at the start of each AC block */
        ps_enc_loop_fxns->pf_ihadamard_scaling_4x4(ps_src->pv_data, ps_src->pv_data, pu2_iscale_mat,
                                                   pu2_weigh_mat, u4_qp_div_6, pi4_tmp);

        /* Update the cntrl flag */
        u4_dc_cntrl = 0;

        for(cntr = 0; cntr < DC_COEFF_CNT_LUMA_MB; cntr++)
        {
            u4_dc_cntrl |= ((((WORD16 *) ps_src->pv_data)[cntr] != 0) << (15 - cntr));
        }

        /* Mark dc bits as 1 if corresponding ac bit is 0 */
        u4_dc_cntrl = (~(u4_cntrl >> 16) & u4_dc_cntrl);

        /* Combine both ac and dc bits */
        u4_cntrl = (u4_cntrl & CNTRL_FLAG_AC_MASK_LUMA) | (u4_dc_cntrl & CNTRL_FLAG_DC_MASK_LUMA);
    }

    /* Source for dc coeffs
     * If the block is intra, we have to read dc values from first row of src
     * then stride for each block is 1, other wise its src stride
     */
    pi2_dc_src = ((WORD16 *) ps_src->pv_data) + (i4_iq_start_idx == 0) * ps_src->i4_data_stride;
    u4_dc_inc = (i4_iq_start_idx == 0) ? ps_src->i4_data_stride : 1;

    /* Get the block bits */
    u4_blk_cntrl = (u4_cntrl & CNTRL_FLAG_AC_MASK_LUMA);
    u4_dc_cntrl = (u4_cntrl & CNTRL_FLAG_DC_MASK_LUMA) << 16;
    u4_empty_blk_cntrl = (~(u4_dc_cntrl | u4_blk_cntrl)) & 0xFFFF0000;

    /* Get first block to process */
    DEQUEUE_BLKID_FROM_CONTROL(u4_dc_cntrl, u4_blk_id);

    while(u4_blk_id < NUM_LUMA4x4_BLOCKS_IN_MB)
    {
        /* Compute address of src blocks */
        WORD32 i4_src_offset = u4_dc_inc * u4_blk_id;

        /* Tx blk coeffs are stored blk by blk */
        /* Hence, in order to access rows of each Tx blk, one needs to stride of
         * TxxSize */
        s_src.i4_data_stride = 4;
        s_src.pv_data = pi2_dc_src + i4_src_offset;

        IND2SUB_LUMA_MB(u4_blk_id, i4_offset_x, i4_offset_y);

        /* Compute address of out and pred blocks */
        s_pred.pv_data =
            ((UWORD8 *) ps_pred->pv_data) + i4_offset_x + i4_offset_y * ps_pred->i4_data_stride;
        s_recon.pv_data =
            ((UWORD8 *) ps_recon->pv_data) + i4_offset_x + i4_offset_y * ps_recon->i4_data_stride;
        s_res.pv_data =
            ((WORD16 *) ps_res->pv_data) + i4_offset_x + i4_offset_y * ps_res->i4_data_stride;
        s_res_pred.pv_data = ((WORD16 *) ps_res_pred->pv_data) + i4_offset_x +
                             i4_offset_y * ps_res_pred->i4_data_stride;

        ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4_dc[u1_iq_it_recon_fxn_idx](
            &s_src, &s_pred, &s_res_pred, &s_res, &s_recon, ps_iq_it_res_rec_constants, NULL,
            pi2_dc_src + i4_src_offset, i4_iq_start_idx, u1_res_accumulate);

        DEQUEUE_BLKID_FROM_CONTROL(u4_dc_cntrl, u4_blk_id);
    }

    /* now process ac/mixed blocks */
    DEQUEUE_BLKID_FROM_CONTROL(u4_blk_cntrl, u4_blk_id);
    while(u4_blk_id < NUM_LUMA4x4_BLOCKS_IN_MB)
    {
        IND2SUB_LUMA_MB(u4_blk_id, i4_offset_x, i4_offset_y);

        /* Tx blk coeffs are stored blk by blk */
        /* Hence, in order to access rows of each Tx blk, one needs to stride of
         * TxxSize */
        s_src.i4_data_stride = 4;
        /* The AC blocks starts from 2nd row */
        s_src.pv_data = ((WORD16 *) ps_src->pv_data) + (u4_blk_id + 1) * ps_src->i4_data_stride;

        s_pred.pv_data =
            ((UWORD8 *) ps_pred->pv_data) + i4_offset_x + i4_offset_y * ps_pred->i4_data_stride;
        s_recon.pv_data =
            ((UWORD8 *) ps_recon->pv_data) + i4_offset_x + i4_offset_y * ps_recon->i4_data_stride;
        s_res.pv_data =
            ((WORD16 *) ps_res->pv_data) + i4_offset_x + i4_offset_y * ps_res->i4_data_stride;
        s_res_pred.pv_data = ((WORD16 *) ps_res_pred->pv_data) + i4_offset_x +
                             i4_offset_y * ps_res_pred->i4_data_stride;

        ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4[u1_iq_it_recon_fxn_idx](
            &s_src, &s_pred, &s_res_pred, &s_res, &s_recon, ps_iq_it_res_rec_constants,
            (WORD16 *) pi4_tmp, pi2_dc_src + u4_blk_id, i4_iq_start_idx, u1_res_accumulate);

        DEQUEUE_BLKID_FROM_CONTROL(u4_blk_cntrl, u4_blk_id);
    }

    /* Now process empty blocks */
    DEQUEUE_BLKID_FROM_CONTROL(u4_empty_blk_cntrl, u4_blk_id);
    while(u4_blk_id < NUM_LUMA4x4_BLOCKS_IN_MB)
    {
        IND2SUB_LUMA_MB(u4_blk_id, i4_offset_x, i4_offset_y);

        /* Tx blk coeffs are stored blk by blk */
        /* Hence, in order to access rows of each Tx blk, one needs to stride of
         * TxxSize */
        s_src.i4_data_stride = 4;
        /* The AC blocks starts from 2nd row */
        s_src.pv_data = ((WORD16 *) ps_src->pv_data) + (u4_blk_id + 1) * ps_src->i4_data_stride;

        s_pred.pv_data =
            ((UWORD8 *) ps_pred->pv_data) + i4_offset_x + i4_offset_y * ps_pred->i4_data_stride;
        s_recon.pv_data =
            ((UWORD8 *) ps_recon->pv_data) + i4_offset_x + i4_offset_y * ps_recon->i4_data_stride;
        s_res.pv_data =
            ((WORD16 *) ps_res->pv_data) + i4_offset_x + i4_offset_y * ps_res->i4_data_stride;
        s_res_pred.pv_data = ((WORD16 *) ps_res_pred->pv_data) + i4_offset_x +
                             i4_offset_y * ps_res_pred->i4_data_stride;

        ps_enc_loop_fxns->pf_zcbf_iquant_itrans_recon_4x4(
            &s_src, &s_pred, &s_res_pred, &s_res, &s_recon, ps_iq_it_res_rec_constants, NULL,
            pi2_dc_src + u4_blk_id, i4_iq_start_idx, u1_res_accumulate);

        DEQUEUE_BLKID_FROM_CONTROL(u4_empty_blk_cntrl, u4_blk_id);
    }
}

/**
*******************************************************************************
*
* @brief
*  This function performs does the DCT transform then Hadamard transform
*  and quantization for a chroma macroblock
*
* @par Description:
*  First  cf4 is done on all 16 4x4 blocks of the 8x8input block
*  Then hadamard transform is done on the DC coefficients
*  Quantization is then performed on the 8x8 block, 4x4 wise
*
* @param[in] pu1_src
*  Pointer to source sub-block
*  The input is in interleaved format for two chroma planes
*
* @param[in] pu1_pred
*  Pointer to prediction sub-block
*  Prediction is in inter leaved format
*
* @param[in] pi2_out
*  Pointer to residual sub-block
*  The output will be in linear format
*  The first 4 continuous locations will contain the values of DC block for U
*  and then next 4 will contain for V.
*  After DC block and a stride 1st AC block of U plane will follow
*  After one more stride next AC block of V plane will follow
*  The blocks will be in raster scan order
*
*  After all the AC blocks of U plane AC blocks of V plane will follow in exact
*  same way
*
* @param[in] i4_src_stride
*  Source stride
*
* @param[in] i4_pred_stride
*  Prediction stride
*
* @param[in] dst_strd
*  Destination stride
*
* @param[in] pu2_scale_matrix
*  The quantization matrix for 4x4 transform
*
* @param[in] pu2_threshold_matrix
*  Threshold matrix
*
* @param[in] u4_qbits
*  15+QP/6
*
* @param[in] u4_round_factor
*  Round factor for quant
*
* @param[out] pu1_nnz
*  Memory to store the non-zeros after transform
*  The first byte will be the nnz od DC block for U plane
*  From the next byte the AC nnzs will be storerd in raster scan order
*  The fifth byte will be nnz of Dc block of V plane
*  Then Ac blocks will follow
*
* @param u4_dc_flag
*  Signals if Dc transform is to be done or not
*   1 -> Dc transform will be done
*   0 -> Dc transform will not be done
*
* @remarks
*
*******************************************************************************
*/
void isvce_chroma_8x8_resi_trans_dctrans_quant(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_quant_coeffs,
    buffer_container_t *ps_upsampled_res, isa_dependent_fxns_t *ps_isa_dependent_fxns,
    const UWORD16 *pu2_scale_matrix, const UWORD16 *pu2_threshold_matrix, UWORD8 *pu1_nnz,
    UWORD32 u4_qbits, UWORD32 u4_round_factor, UWORD8 u1_use_upsampled_res)
{
    WORD32 blk_cntr;
    WORD32 i4_offsetx, i4_offsety;
    UWORD8 au1_dcnnz[2];

    enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;
    buffer_container_t s_src = ps_src[0];
    buffer_container_t s_pred = ps_pred[0];
    buffer_container_t s_quant_coeffs = ps_quant_coeffs[0];
    buffer_container_t s_upsampled_res = {0};
    resi_trans_quant_constants_t s_resi_trans_quant_constants = {
        .pu2_scale_matrix = pu2_scale_matrix,
        .pu2_threshold_matrix = pu2_threshold_matrix,
        .u4_qbits = u4_qbits,
        .u4_round_factor = u4_round_factor};

    UWORD8 u1_resi_trans_fxn_idx = isvc_get_resi_trans_quant_variant_idx(u1_use_upsampled_res);

    if(u1_use_upsampled_res)
    {
        s_upsampled_res = ps_upsampled_res[0];
    }

    /* Move to the ac addresses */
    pu1_nnz++;

    for(blk_cntr = 0; blk_cntr < NUM_CHROMA4x4_BLOCKS_IN_MB; blk_cntr++)
    {
        IND2SUB_CHROMA_MB(blk_cntr, i4_offsetx, i4_offsety);

        s_src.pv_data =
            ((UWORD8 *) ps_src[0].pv_data) + i4_offsetx + i4_offsety * ps_src[0].i4_data_stride;
        s_pred.pv_data =
            ((UWORD8 *) ps_pred[0].pv_data) + i4_offsetx + i4_offsety * ps_pred[0].i4_data_stride;
        s_quant_coeffs.pv_data =
            ((WORD16 *) ps_quant_coeffs[0].pv_data) + blk_cntr * ps_quant_coeffs[0].i4_data_stride;

        if(u1_use_upsampled_res)
        {
            s_upsampled_res.pv_data = ((WORD16 *) ps_upsampled_res[0].pv_data) + i4_offsetx +
                                      i4_offsety * ps_upsampled_res[0].i4_data_stride;
        }

        /* Move to the ac addresses */
        s_quant_coeffs.pv_data =
            ((WORD16 *) s_quant_coeffs.pv_data) + ps_quant_coeffs[0].i4_data_stride;

        s_quant_coeffs.i4_data_stride = 4;

        /* For chroma, v plane nnz is populated from position 5 */
        ps_enc_loop_fxns->apf_resi_trans_quant_chroma_4x4[u1_resi_trans_fxn_idx](
            &s_src, &s_pred, &s_quant_coeffs, &s_upsampled_res, &s_resi_trans_quant_constants,
            &pu1_nnz[blk_cntr + (blk_cntr > 3)], ((WORD16 *) ps_quant_coeffs->pv_data) + blk_cntr,
            u1_use_upsampled_res);
    }

    /* Adjust pointers to point to dc values */
    s_quant_coeffs = ps_quant_coeffs[0];
    pu1_nnz--;

    s_resi_trans_quant_constants.u4_qbits++;
    s_resi_trans_quant_constants.u4_round_factor <<= 1;

    ps_enc_loop_fxns->pf_hadamard_quant_2x2_uv(((WORD16 *) ps_quant_coeffs->pv_data),
                                               ((WORD16 *) ps_quant_coeffs->pv_data),
                                               &s_resi_trans_quant_constants, au1_dcnnz);

    /* Copy the dc nnzs */
    pu1_nnz[0] = au1_dcnnz[0];
    pu1_nnz[5] = au1_dcnnz[1];
}

/**
*******************************************************************************
* @brief
*  This function performs the inverse transform with process for chroma MB of
*H264
*
* @par Description:
*  Does inverse DC transform ,inverse quantization inverse transform
*
* @param[in] pi2_src
*  Input data, 16x16 size
*  The input is in the form of, first 4 locations will contain DC coeffs of
*  U plane, next 4 will contain DC coeffs of V plane, then AC blocks of U plane
*  in raster scan order will follow, each block as linear array in raster scan
*order. After a stride next AC block will follow. After all AC blocks of U plane
*  V plane AC blocks will follow in exact same order.
*
* @param[in] pu1_pred
*  The predicted data, 8x16 size, U and V interleaved
*
* @param[in] pu1_out
*  Output 8x16, U and V interleaved
*
* @param[in] i4_src_stride
*  Source stride
*
* @param[in] i4_pred_stride
*  input stride for prediction buffer
*
* @param[in] i4_out_stride
*  input stride for output buffer
*
* @param[in] pu2_iscale_mat
*  Inverse quantization martix for 4x4 transform
*
* @param[in] pu2_weigh_mat
*  weight matrix of 4x4 transform
*
* @param[in] u4_qp_div_6
*  QP/6
*
* @param[in] pi4_tmp
*  Input temporary buffer
*  needs to be at least COFF_CNT_SUB_BLK_4x4 + Number of Dc cofss for chroma *
*number of planes in size
*
* @param[in] pu4_cntrl
*  Controls the transform path
*  the 15 th bit will correspond to DC block of U plane , 14th will indicate the
*V plane Dc block 32-28 bits will indicate AC blocks of U plane in raster scan
*order 27-23 bits will indicate AC blocks of V plane in rater scan order The bit
*1 implies that there is at least one non zero coeff in a block
*
* @returns
*  none
*
* @remarks
*******************************************************************************
*/
void isvce_chroma_8x8_idctrans_iquant_itrans_recon(
    buffer_container_t *ps_src, buffer_container_t *ps_pred, buffer_container_t *ps_recon,
    buffer_container_t *ps_res, buffer_container_t *ps_res_pred,
    iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants,
    isa_dependent_fxns_t *ps_isa_dependent_fxns, WORD32 *pi4_tmp, UWORD32 u4_cntrl,
    UWORD8 u1_res_accumulate)
{
    /* Cntrl bits for 4x4 transforms
     * u4_blk_cntrl       : controls if a 4x4 block should be processed in ac path
     * u4_dc_cntrl        : controls is a 4x4 block is to be processed in dc path
     *                    : dc block must contain only single dc coefficient
     * u4_empty_blk_cntrl : control fot 4x4 block with no coeffs, ie no dc and ac
     *                    : ie not (ac or dc)
     */
    UWORD32 u4_blk_cntrl, u4_dc_cntrl, u4_empty_blk_cntrl;
    WORD32 u4_blk_id;
    WORD32 i4_offset_x, i4_offset_y;
    WORD16 *pi2_dc_src;
    /* Increment for dc block */
    WORD32 i4_dc_inc;

    enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;
    buffer_container_t s_src = ps_src[0];
    buffer_container_t s_pred = ps_pred[0];
    buffer_container_t s_recon = ps_recon[0];
    buffer_container_t s_res = ps_res[0];
    buffer_container_t s_res_pred = ps_res_pred[0];

    WORD16 i2_zero = 0;
    const UWORD16 *pu2_iscale_mat = ps_iq_it_res_rec_constants->pu2_iscal_mat;
    const UWORD16 *pu2_weigh_mat = ps_iq_it_res_rec_constants->pu2_weigh_mat;
    UWORD32 u4_qp_div_6 = ps_iq_it_res_rec_constants->u4_qp_div_6;
    UWORD8 u1_iq_it_recon_fxn_idx = isvc_get_iq_it_recon_variant_idx(0, u1_res_accumulate);

    /*
     * Lets do the inverse transform for dc coeffs in chroma
     */
    if(u4_cntrl & CNTRL_FLAG_DCBLK_MASK_CHROMA)
    {
        UWORD32 cntr, u4_dc_cntrl;
        /* Do inv hadamard for u an v block */

        ps_enc_loop_fxns->pf_ihadamard_scaling_2x2_uv(s_src.pv_data, s_src.pv_data, pu2_iscale_mat,
                                                      pu2_weigh_mat, u4_qp_div_6, NULL);
        /*
         * Update the cntrl flag
         * Flag is updated as follows bits 15-11 -> u block dc bits
         */
        u4_dc_cntrl = 0;
        for(cntr = 0; cntr < 8; cntr++)
        {
            u4_dc_cntrl |= ((((WORD16 *) ps_src->pv_data)[cntr] != 0) << (15 - cntr));
        }

        /* Mark dc bits as 1 if corresponding ac bit is 0 */
        u4_dc_cntrl = (~(u4_cntrl >> 16) & u4_dc_cntrl);
        /* Combine both ac and dc bits */
        u4_cntrl =
            (u4_cntrl & CNTRL_FLAG_AC_MASK_CHROMA) | (u4_dc_cntrl & CNTRL_FLAG_DC_MASK_CHROMA);

        /* Since we populated the dc coffs, we have to read them from there */
        pi2_dc_src = ((WORD16 *) ps_src->pv_data);
        i4_dc_inc = 1;
    }
    else
    {
        u4_cntrl = u4_cntrl & CNTRL_FLAG_AC_MASK_CHROMA;
        pi2_dc_src = &i2_zero;
        i4_dc_inc = 0;
    }

    /* Get the block bits */
    u4_blk_cntrl = (u4_cntrl & CNTRL_FLAG_AC_MASK_CHROMA);
    u4_dc_cntrl = (u4_cntrl & CNTRL_FLAG_DC_MASK_CHROMA) << 16;
    u4_empty_blk_cntrl = (~(u4_dc_cntrl | u4_blk_cntrl)) & 0xFF000000;

    DEQUEUE_BLKID_FROM_CONTROL(u4_dc_cntrl, u4_blk_id);

    while(u4_blk_id < 8)
    {
        WORD32 dc_src_offset = u4_blk_id * i4_dc_inc;

        /* Tx blk coeffs are stored blk by blk */
        /* Hence, in order to access rows of each Tx blk, one needs to stride of
         * TxxSize */
        s_src.i4_data_stride = 4;
        s_src.pv_data = pi2_dc_src + dc_src_offset;

        IND2SUB_CHROMA_MB(u4_blk_id, i4_offset_x, i4_offset_y);

        s_pred.pv_data =
            ((UWORD8 *) ps_pred->pv_data) + i4_offset_x + i4_offset_y * ps_pred->i4_data_stride;
        s_recon.pv_data =
            ((UWORD8 *) ps_recon->pv_data) + i4_offset_x + i4_offset_y * ps_recon->i4_data_stride;
        s_res.pv_data =
            ((WORD16 *) ps_res->pv_data) + i4_offset_x + i4_offset_y * ps_res->i4_data_stride;
        s_res_pred.pv_data = ((WORD16 *) ps_res_pred->pv_data) + i4_offset_x +
                             i4_offset_y * ps_res_pred->i4_data_stride;

        ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4_dc[u1_iq_it_recon_fxn_idx](
            &s_src, &s_pred, &s_res_pred, &s_res, &s_recon, ps_iq_it_res_rec_constants, NULL,
            s_src.pv_data, 0, u1_res_accumulate);

        /* Get next DC block to process */
        DEQUEUE_BLKID_FROM_CONTROL(u4_dc_cntrl, u4_blk_id);
    }

    /* now process ac/mixed blocks */
    DEQUEUE_BLKID_FROM_CONTROL(u4_blk_cntrl, u4_blk_id);
    while(u4_blk_id < 8)
    {
        WORD32 dc_src_offset = i4_dc_inc * u4_blk_id;

        /* Tx blk coeffs are stored blk by blk */
        /* Hence, in order to access rows of each Tx blk, one needs to stride of
         * TxxSize */
        s_src.i4_data_stride = 4;
        /* The AC blocks starts from 2nd row */
        s_src.pv_data = ((WORD16 *) ps_src->pv_data) + (u4_blk_id + 1) * ps_src->i4_data_stride;

        IND2SUB_CHROMA_MB(u4_blk_id, i4_offset_x, i4_offset_y);

        s_pred.pv_data =
            ((UWORD8 *) ps_pred->pv_data) + i4_offset_x + i4_offset_y * ps_pred->i4_data_stride;
        s_recon.pv_data =
            ((UWORD8 *) ps_recon->pv_data) + i4_offset_x + i4_offset_y * ps_recon->i4_data_stride;
        s_res.pv_data =
            ((WORD16 *) ps_res->pv_data) + i4_offset_x + i4_offset_y * ps_res->i4_data_stride;
        s_res_pred.pv_data = ((WORD16 *) ps_res_pred->pv_data) + i4_offset_x +
                             i4_offset_y * ps_res_pred->i4_data_stride;

        ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4[u1_iq_it_recon_fxn_idx](
            &s_src, &s_pred, &s_res_pred, &s_res, &s_recon, ps_iq_it_res_rec_constants,
            (WORD16 *) pi4_tmp, pi2_dc_src + dc_src_offset, 0, u1_res_accumulate);

        DEQUEUE_BLKID_FROM_CONTROL(u4_blk_cntrl, u4_blk_id);
    }

    /* Now process empty blocks */
    DEQUEUE_BLKID_FROM_CONTROL(u4_empty_blk_cntrl, u4_blk_id);

    while(u4_blk_id < 8)
    {
        WORD32 dc_src_offset = i4_dc_inc * u4_blk_id;

        IND2SUB_CHROMA_MB(u4_blk_id, i4_offset_x, i4_offset_y);

        /* Tx blk coeffs are stored blk by blk */
        /* Hence, in order to access rows of each Tx blk, one needs to stride of
         * TxxSize */
        s_src.i4_data_stride = 4;
        /* The AC blocks starts from 2nd row */
        s_src.pv_data = ((WORD16 *) ps_src->pv_data) + (u4_blk_id + 1) * ps_src->i4_data_stride;

        s_pred.pv_data =
            ((UWORD8 *) ps_pred->pv_data) + i4_offset_x + i4_offset_y * ps_pred->i4_data_stride;
        s_recon.pv_data =
            ((UWORD8 *) ps_recon->pv_data) + i4_offset_x + i4_offset_y * ps_recon->i4_data_stride;
        s_res.pv_data =
            ((WORD16 *) ps_res->pv_data) + i4_offset_x + i4_offset_y * ps_res->i4_data_stride;
        s_res_pred.pv_data = ((WORD16 *) ps_res_pred->pv_data) + i4_offset_x +
                             i4_offset_y * ps_res_pred->i4_data_stride;

        ps_enc_loop_fxns->pf_chroma_zcbf_iquant_itrans_recon_4x4(
            &s_src, &s_pred, &s_res_pred, &s_res, &s_recon, ps_iq_it_res_rec_constants,
            (WORD16 *) pi4_tmp, pi2_dc_src + dc_src_offset, 0, u1_res_accumulate);

        DEQUEUE_BLKID_FROM_CONTROL(u4_empty_blk_cntrl, u4_blk_id);
    }
}

/**
******************************************************************************
*
* @brief  This function packs residue of an i16x16 luma mb for entropy coding
*
* @par   Description
*  An i16 macro block contains two classes of units, dc 4x4 block and
*  4x4 ac blocks. while packing the mb, the dc block is sent first, and
*  the 16 ac blocks are sent next in scan order. Each and every block is
*  represented by 3 parameters (nnz, significant coefficient map and the
*  residue coefficients itself). If a 4x4 unit does not have any coefficients
*  then only nnz is sent. Inside a 4x4 block the individual coefficients are
*  sent in scan order.
*
*  The first byte of each block will be nnz of the block, if it is non zero,
*  a 2 byte significance map is sent. This is followed by nonzero coefficients.
*  This is repeated for 1 dc + 16 ac blocks.
*
* @param[in]  pi2_res_mb
*  pointer to residue mb
*
* @param[in, out]  pv_mb_coeff_data
*  buffer pointing to packed residue coefficients
*
* @param[in]  u4_res_strd
*  residual block stride
*
* @param[out]  u1_cbp_l
*  coded block pattern luma
*
* @param[in]   pu1_nnz
*  number of non zero coefficients in each 4x4 unit
*
* @param[out]
*  Control signal for inverse transform of 16x16 blocks
*
* @return none
*
* @ remarks
*
******************************************************************************
*/
void isvce_pack_l_mb_i16(WORD16 *pi2_res_mb, void **pv_mb_coeff_data, WORD32 i4_res_strd,
                         UWORD8 *u1_cbp_l, UWORD8 *pu1_nnz, UWORD32 *pu4_cntrl)
{
    /* pointer to packed sub block buffer space */
    tu_sblk_coeff_data_t *ps_mb_coeff_data = (*pv_mb_coeff_data), *ps_mb_coeff_data_ac;

    /* no of non zero coefficients in the current sub block */
    UWORD32 u4_nnz_cnt;

    /* significant coefficient map */
    UWORD32 u4_s_map;

    /* pointer to scanning matrix */
    const UWORD8 *pu1_scan_order;

    /* number of non zeros in sub block */
    UWORD32 u4_nnz;

    /* coeff scan order */
    const UWORD8 u1_scan_order[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

    /* temp var */
    UWORD32 coeff_cnt, mask, b4, u4_cntrl = 0;

    /*DC and AC coeff pointers*/
    WORD16 *pi2_res_mb_ac, *pi2_res_mb_dc;

    /********************************************************/
    /*  pack dc coeff data for entropy coding               */
    /********************************************************/

    pi2_res_mb_dc = pi2_res_mb;
    pu1_scan_order = gu1_luma_scan_order_dc;

    u4_nnz = *pu1_nnz;
    u4_cntrl = 0;

    /* write number of non zero coefficients */
    ps_mb_coeff_data->i4_sig_map_nnz = u4_nnz;

    if(u4_nnz)
    {
        for(u4_nnz_cnt = 0, coeff_cnt = 0, mask = 1, u4_s_map = 0; u4_nnz_cnt < u4_nnz; coeff_cnt++)
        {
            if(pi2_res_mb_dc[pu1_scan_order[coeff_cnt]])
            {
                /* write residue */
                ps_mb_coeff_data->ai2_residue[u4_nnz_cnt++] =
                    pi2_res_mb_dc[pu1_scan_order[coeff_cnt]];
                u4_s_map |= mask;
            }
            mask <<= 1;
        }
        /* write significant coeff map */
        ps_mb_coeff_data->i4_sig_map_nnz |= (u4_s_map << 16);
        (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue + ALIGN2(u4_nnz_cnt);

        u4_cntrl = 0x00008000;  // Set DC bit in ctrl code
    }
    else
    {
        (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue;
    }

    /********************************************************/
    /*  pack ac coeff data for entropy coding               */
    /********************************************************/

    pu1_nnz++;
    pu1_scan_order = gu1_luma_scan_order;
    pi2_res_mb += i4_res_strd; /*Move to AC block*/

    ps_mb_coeff_data_ac = (*pv_mb_coeff_data);

    for(b4 = 0; b4 < 16; b4++)
    {
        ps_mb_coeff_data = (*pv_mb_coeff_data);

        u4_nnz = pu1_nnz[u1_scan_order[b4]];

        /* Jump according to the scan order */
        pi2_res_mb_ac = pi2_res_mb + (i4_res_strd * u1_scan_order[b4]);

        /*
         * Since this is a i16x16 block, we should not count dc coeff on indi
         * vidual 4x4 blocks to nnz. But due to the implementation of 16x16
         * trans function, we add dc's nnz to u4_nnz too. Hence we adjust that
         * here
         */
        u4_nnz -= (pi2_res_mb_ac[0] != 0);

        /* write number of non zero coefficients */
        ps_mb_coeff_data->i4_sig_map_nnz = u4_nnz;

        if(u4_nnz)
        {
            for(u4_nnz_cnt = 0, coeff_cnt = 1, mask = 1, u4_s_map = 0; u4_nnz_cnt < u4_nnz;
                coeff_cnt++)
            {
                if(pi2_res_mb_ac[pu1_scan_order[coeff_cnt]])
                {
                    /* write residue */
                    ps_mb_coeff_data->ai2_residue[u4_nnz_cnt++] =
                        pi2_res_mb_ac[pu1_scan_order[coeff_cnt]];
                    u4_s_map |= mask;
                }
                mask <<= 1;
            }
            /* write significant coeff map */
            ps_mb_coeff_data->i4_sig_map_nnz |= (u4_s_map << 16);
            (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue + ALIGN2(u4_nnz_cnt);
            *u1_cbp_l = 15;

            u4_cntrl |= (1 << (31 - u1_scan_order[b4]));
        }
        else
        {
            (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue;
        }
    }

    if(!(*u1_cbp_l))
    {
        (*pv_mb_coeff_data) = ps_mb_coeff_data_ac;
    }

    /* Store the cntrl signal */
    (*pu4_cntrl) = u4_cntrl;
    return;
}

/**
******************************************************************************
*
* @brief  This function packs residue of an p16x16 luma mb for entropy coding
*
* @par   Description
*  A p16x16 macro block contains two classes of units 16  4x4 ac blocks.
*  while packing the mb, the dc block is sent first, and
*  the 16 ac blocks are sent next in scan order. Each and every block is
*  represented by 3 parameters (nnz, significant coefficient map and the
*  residue coefficients itself). If a 4x4 unit does not have any coefficients
*  then only nnz is sent. Inside a 4x4 block the individual coefficients are
*  sent in scan order.
*
*  The first byte of each block will be nnz of the block, if it is non zero,
*  a 2 byte significance map is sent. This is followed by nonzero coefficients.
*  This is repeated for 1 dc + 16 ac blocks.
*
* @param[in]  pi2_res_mb
*  pointer to residue mb
*
* @param[in, out]  pv_mb_coeff_data
*  buffer pointing to packed residue coefficients
*
* @param[in]  i4_res_strd
*  residual block stride
*
* @param[out]  u1_cbp_l
*  coded block pattern luma
*
* @param[in]   pu1_nnz
*  number of non zero coefficients in each 4x4 unit
*
* @param[out] pu4_cntrl
*  Control signal for inverse transform
*
* @return none
*
* @remarks Killing coffs not yet coded
*
******************************************************************************
*/
void isvce_pack_l_mb(WORD16 *pi2_res_mb, void **pv_mb_coeff_data, WORD32 i4_res_strd,
                     UWORD8 *u1_cbp_l, UWORD8 *pu1_nnz, UWORD32 u4_thres_resi, UWORD32 *pu4_cntrl)
{
    /* pointer to packed sub block buffer space */
    tu_sblk_coeff_data_t *ps_mb_coeff_data, *ps_mb_coeff_data_b8, *ps_mb_coeff_data_mb;

    /* no of non zero coefficients in the current sub block */
    UWORD32 u4_nnz_cnt;

    /* significant coefficient map */
    UWORD32 u4_s_map;

    /* pointer to scanning matrix */
    const UWORD8 *pu1_scan_order = gu1_luma_scan_order;

    /* number of non zeros in sub block */
    UWORD32 u4_nnz;

    /* pointer to residual sub block */
    WORD16 *pi2_res_sb;

    /* coeff scan order */
    const UWORD8 u1_scan_order[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

    /* coeff cost */
    const UWORD8 *pu1_coeff_cost = gu1_coeff_cost;

    /* temp var */
    UWORD32 u4_mb_coeff_cost = 0, u4_b8_coeff_cost = 0, coeff_cnt, mask, u4_cntrl = 0, b4, b8;

    /* temp var */
    WORD32 i4_res_val, i4_run = -1, dcac_block;

    /* When Hadamard transform is disabled, first row values are dont care, ignore
     * them */
    pi2_res_mb += i4_res_strd;

    /* When Hadamard transform is disabled, first unit value is dont care, ignore
     * this */
    pu1_nnz++;

    ps_mb_coeff_data_mb = ps_mb_coeff_data_b8 = (*pv_mb_coeff_data);

    /********************************************************/
    /*  pack coeff data for entropy coding                  */
    /********************************************************/

    for(b4 = 0; b4 < 16; b4++)
    {
        ps_mb_coeff_data = (*pv_mb_coeff_data);

        b8 = b4 >> 2;

        u4_nnz = pu1_nnz[u1_scan_order[b4]];

        /* Jump according to the scan order */
        pi2_res_sb = pi2_res_mb + (i4_res_strd * u1_scan_order[b4]);

        /* write number of non zero coefficients */
        ps_mb_coeff_data->i4_sig_map_nnz = u4_nnz;

        if(u4_nnz)
        {
            for(u4_nnz_cnt = 0, coeff_cnt = 0, mask = 1, u4_s_map = 0; u4_nnz_cnt < u4_nnz;
                coeff_cnt++)
            {
                /* number of runs of zero before, this is used to compute coeff cost */
                i4_run++;

                i4_res_val = pi2_res_sb[pu1_scan_order[coeff_cnt]];

                if(i4_res_val)
                {
                    /* write residue */
                    ps_mb_coeff_data->ai2_residue[u4_nnz_cnt++] = i4_res_val;
                    u4_s_map |= mask;

                    if(u4_thres_resi)
                    {
                        /* compute coeff cost */
                        if(i4_res_val == 1 || i4_res_val == -1)
                        {
                            if(i4_run < 6) u4_b8_coeff_cost += pu1_coeff_cost[i4_run];
                        }
                        else
                            u4_b8_coeff_cost += 9;

                        i4_run = -1;
                    }
                }

                mask <<= 1;
            }

            /* write significant coeff map */
            ps_mb_coeff_data->i4_sig_map_nnz |= (u4_s_map << 16);
            (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue + ALIGN2(u4_nnz_cnt);

            /* cbp */
            *u1_cbp_l |= (1 << b8);

            /* Cntrl map for inverse transform computation
             *
             * If coeff_cnt is zero, it means that only nonzero was a dc coeff
             * Hence we have to set the 16 - u1_scan_order[b4]) position instead
             * of 31 - u1_scan_order[b4]
             */
            dcac_block = (coeff_cnt == 0) ? 16 : 31;
            u4_cntrl |= (1 << (dcac_block - u1_scan_order[b4]));
        }
        else
        {
            (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue;
        }

        /* Decide if the 8x8 unit has to be sent for entropy coding? */
        if((b4 + 1) % 4 == 0)
        {
            if(u4_thres_resi && (u4_b8_coeff_cost <= LUMA_SUB_BLOCK_SKIP_THRESHOLD) &&
               (*u1_cbp_l & (1 << b8)))
            {
                /*
                 * When we want to reset the full 8x8 block, we have to reset
                 * both the dc and ac coeff bits hence we have the symmetric
                 * arrangement of bits
                 */
                const UWORD32 cntrl_mask_map[4] = {0xcc00cc00, 0x33003300, 0x00cc00cc, 0x00330033};

                /* restore cbp */
                *u1_cbp_l = (*u1_cbp_l & (~(1 << b8)));

                /* correct cntrl flag */
                u4_cntrl = u4_cntrl & (~cntrl_mask_map[(b4 >> 2)]);

                /* correct nnz */
                pu1_nnz[u1_scan_order[b4 - 3]] = 0;
                pu1_nnz[u1_scan_order[b4 - 2]] = 0;
                pu1_nnz[u1_scan_order[b4 - 1]] = 0;
                pu1_nnz[u1_scan_order[b4]] = 0;

                /* reset blk cost */
                u4_b8_coeff_cost = 0;
            }

            if(!(*u1_cbp_l & (1 << b8)))
            {
                (*pv_mb_coeff_data) = ps_mb_coeff_data_b8;
            }

            u4_mb_coeff_cost += u4_b8_coeff_cost;

            u4_b8_coeff_cost = 0;
            i4_run = -1;
            ps_mb_coeff_data_b8 = (*pv_mb_coeff_data);
        }
    }

    if(u4_thres_resi && (u4_mb_coeff_cost <= LUMA_BLOCK_SKIP_THRESHOLD) && (*u1_cbp_l))
    {
        (*pv_mb_coeff_data) = ps_mb_coeff_data_mb;
        *u1_cbp_l = 0;
        u4_cntrl = 0;
        memset(pu1_nnz, 0, 16);
    }

    (*pu4_cntrl) = u4_cntrl;

    return;
}

/**
******************************************************************************
*
* @brief  This function packs residue of an i8x8 chroma mb for entropy coding
*
* @par   Description
*  An i8 chroma macro block contains two classes of units, dc 2x2 block and
*  4x4 ac blocks. while packing the mb, the dc block is sent first, and
*  the 4 ac blocks are sent next in scan order. Each and every block is
*  represented by 3 parameters (nnz, significant coefficient map and the
*  residue coefficients itself). If a 4x4 unit does not have any coefficients
*  then only nnz is sent. Inside a 4x4 block the individual coefficients are
*  sent in scan order.
*
*  The first byte of each block will be nnz of the block, if it is non zero,
*  a 2 byte significance map is sent. This is followed by nonzero coefficients.
*  This is repeated for 1 dc + 4 ac blocks.
*
* @param[in]  pi2_res_mb
*  pointer to residue mb
*
* @param[in, out]  pv_mb_coeff_data
*  buffer pointing to packed residue coefficients
*
* @param[in]  u4_res_strd
*  residual block stride
*
* @param[out]  u1_cbp_c
*  coded block pattern chroma
*
* @param[in]   pu1_nnz
*  number of non zero coefficients in each 4x4 unit
*
* @param[out]   pu1_nnz
*  Control signal for inverse transform
*
* @param[in]   u4_swap_uv
*  Swaps the order of U and V planes in entropy bitstream
*
* @return none
*
* @ remarks
*
******************************************************************************
*/
void isvce_pack_c_mb(WORD16 *pi2_res_mb, void **pv_mb_coeff_data, WORD32 i4_res_strd,
                     UWORD8 *u1_cbp_c, UWORD8 *pu1_nnz, UWORD32 u4_thres_resi, UWORD32 *pu4_cntrl,
                     UWORD32 u4_swap_uv)
{
    /* pointer to packed sub block buffer space */
    tu_sblk_coeff_data_t *ps_mb_coeff_data = (*pv_mb_coeff_data);
    tu_sblk_coeff_data_t *ps_mb_coeff_data_dc, *ps_mb_coeff_data_ac;

    /* nnz pointer */
    UWORD8 *pu1_nnz_ac, *pu1_nnz_dc;

    /* nnz counter */
    UWORD32 u4_nnz_cnt;

    /* significant coefficient map */
    UWORD32 u4_s_map;

    /* pointer to scanning matrix */
    const UWORD8 *pu1_scan_order;

    /* no of non zero coefficients in the current sub block */
    UWORD32 u4_nnz;

    /* pointer to residual sub block, res val */
    WORD16 *pi2_res_sb, i2_res_val;

    /* temp var */
    UWORD32 coeff_cnt, mask, b4, plane;

    /* temp var */
    UWORD32 u4_coeff_cost;
    WORD32 i4_run;

    /* coeff cost */
    const UWORD8 *pu1_coeff_cost = gu1_coeff_cost;

    /* pointer to packed buffer space */
    UWORD32 *pu4_mb_coeff_data = NULL;

    /* ac coded block pattern */
    UWORD8 u1_cbp_ac;

    /* Variable to store the current bit pos in cntrl variable*/
    UWORD32 cntrl_pos = 0;

    /********************************************************/
    /*  pack dc coeff data for entropy coding               */
    /********************************************************/
    pu1_scan_order = gu1_chroma_scan_order_dc;
    pi2_res_sb = pi2_res_mb;
    pu1_nnz_dc = pu1_nnz;
    (*pu4_cntrl) = 0;
    cntrl_pos = 15;
    ps_mb_coeff_data_dc = (*pv_mb_coeff_data);

    /* Color space conversion between SP_UV and SP_VU
     * We always assume SP_UV for all the processing
     * Hence to get proper stream output we need to swap U and V channels here
     *
     * For that there are two paths we need to look for
     * One is the path to bitstream , these variables should have the proper input
     * configured UV or VU
     * For the other path the inverse transform variables should have what ever
     * ordering the input had
     */

    if(u4_swap_uv)
    {
        pu1_nnz_dc += 5; /* Move to NNZ of V planve */
        pi2_res_sb += 4; /* Move to DC coff of V plane */

        cntrl_pos = 14; /* Control bit for V plane */
    }

    for(plane = 0; plane < 2; plane++)
    {
        ps_mb_coeff_data = (*pv_mb_coeff_data);

        u4_nnz = *pu1_nnz_dc;
        /* write number of non zero coefficients U/V */
        ps_mb_coeff_data->i4_sig_map_nnz = u4_nnz;

        if(u4_nnz)
        {
            for(u4_nnz_cnt = 0, coeff_cnt = 0, mask = 1, u4_s_map = 0; u4_nnz_cnt < u4_nnz;
                coeff_cnt++)
            {
                i2_res_val = pi2_res_sb[pu1_scan_order[coeff_cnt]];
                if(i2_res_val)
                {
                    /* write residue U/V */
                    ps_mb_coeff_data->ai2_residue[u4_nnz_cnt++] = i2_res_val;
                    u4_s_map |= mask;
                }
                mask <<= 1;
            }
            /* write significant coeff map U/V */
            ps_mb_coeff_data->i4_sig_map_nnz |= (u4_s_map << 16);
            (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue + ALIGN2(u4_nnz_cnt);
            *u1_cbp_c = 1;

            (*pu4_cntrl) |= (1 << cntrl_pos);
        }
        else
        {
            (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue;
        }

        if(u4_swap_uv)
        {
            cntrl_pos++;     /* Control bit for U plane */
            pu1_nnz_dc -= 5; /* Move to NNZ of U plane */
            pi2_res_sb -= 4; /* Move to DC coff of U plane */
        }
        else
        {
            cntrl_pos--;     /* Control bit for U plane */
            pu1_nnz_dc += 5; /* 4 for AC NNZ and 1 for DC */
            pi2_res_sb += 4; /* Move to DC coff of V plane */
        }
    }

    /********************************************************/
    /*  pack ac coeff data for entropy coding               */
    /********************************************************/

    pu1_scan_order = gu1_chroma_scan_order;
    ps_mb_coeff_data_ac = (*pv_mb_coeff_data);

    if(u4_swap_uv)
    {
        pi2_res_sb = pi2_res_mb + i4_res_strd * 5; /* Move to V plane ,ie 1dc row+ 4 ac row */
        cntrl_pos = 27;           /* The control bits are to be added for V bloc ie 31-4 th bit */
        pu1_nnz_ac = pu1_nnz + 6; /*Move the nnz to V block NNZ 1 dc + 1dc + 4 ac */
    }
    else
    {
        pi2_res_sb = pi2_res_mb + i4_res_strd; /* Move to U plane ,ie 1dc row */
        cntrl_pos = 31;
        pu1_nnz_ac = pu1_nnz + 1; /* Move the nnz to V block NNZ 1 dc */
    }

    for(plane = 0; plane < 2; plane++)
    {
        pu4_mb_coeff_data = (*pv_mb_coeff_data);

        u4_coeff_cost = 0;
        i4_run = -1;

        /* get the current cbp, so that it automatically
         * gets reverted in case of zero ac values */
        u1_cbp_ac = *u1_cbp_c;

        for(b4 = 0; b4 < 4; b4++)
        {
            ps_mb_coeff_data = (*pv_mb_coeff_data);

            u4_nnz = *pu1_nnz_ac;

            /*
             * We are scanning only ac coeffs, but the nnz is for the
             * complete 4x4 block. Hence we have to discount the nnz contributed
             * by the dc coefficient
             */
            u4_nnz -= (pi2_res_sb[0] != 0);

            /* write number of non zero coefficients U/V */
            ps_mb_coeff_data->i4_sig_map_nnz = u4_nnz;

            if(u4_nnz)
            {
                for(u4_nnz_cnt = 0, coeff_cnt = 0, mask = 1, u4_s_map = 0; u4_nnz_cnt < u4_nnz;
                    coeff_cnt++)
                {
                    i2_res_val = pi2_res_sb[pu1_scan_order[coeff_cnt]];

                    i4_run++;

                    if(i2_res_val)
                    {
                        /* write residue U/V */
                        ps_mb_coeff_data->ai2_residue[u4_nnz_cnt++] = i2_res_val;
                        u4_s_map |= mask;

                        if(u4_thres_resi && (u4_coeff_cost < CHROMA_BLOCK_SKIP_THRESHOLD))
                        {
                            /* compute coeff cost */
                            if(i2_res_val == 1 || i2_res_val == -1)
                            {
                                if(i4_run < 6) u4_coeff_cost += pu1_coeff_cost[i4_run];
                            }
                            else
                                u4_coeff_cost += 9;

                            i4_run = -1;
                        }
                    }
                    mask <<= 1;
                }

                /* write significant coeff map U/V */
                ps_mb_coeff_data->i4_sig_map_nnz |= (u4_s_map << 16);
                (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue + ALIGN2(u4_nnz_cnt);
                u1_cbp_ac = 2;

                (*pu4_cntrl) |= 1 << cntrl_pos;
            }
            else
            {
                (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue;
            }

            pu1_nnz_ac++;
            pi2_res_sb += i4_res_strd;
            cntrl_pos--;
        }

        /* reset block */
        if(u4_thres_resi && (u4_coeff_cost < CHROMA_BLOCK_SKIP_THRESHOLD))
        {
            pu4_mb_coeff_data[0] = 0;
            pu4_mb_coeff_data[1] = 0;
            pu4_mb_coeff_data[2] = 0;
            pu4_mb_coeff_data[3] = 0;
            (*pv_mb_coeff_data) = pu4_mb_coeff_data + 4;

            /* Generate the control signal */
            /* Zero out the current plane's AC coefficients */
            (*pu4_cntrl) &= ((plane == u4_swap_uv) ? 0x0FFFFFFF : 0xF0FFFFFF);

            /* Similarly do for the NNZ also */
            *(pu1_nnz_ac - 4) = 0;
            *(pu1_nnz_ac - 3) = 0;
            *(pu1_nnz_ac - 2) = 0;
            *(pu1_nnz_ac - 1) = 0;
        }
        else
        {
            *u1_cbp_c = u1_cbp_ac;
        }

        if(u4_swap_uv)
        {
            pi2_res_sb =
                pi2_res_mb + i4_res_strd; /* Move to V plane ,ie 1dc row+ 4 ac row + 1 dc row */
            cntrl_pos = 31; /* The control bits are to be added for V bloc ie 31-4 th bit */
            pu1_nnz_ac = pu1_nnz + 1; /* Move the nnz to V block NNZ 1 dc + 1dc + 4 ac */

            pu1_nnz_ac = pu1_nnz + 1;
        }
        else
            pu1_nnz_ac = pu1_nnz + 6; /* Go to nnz of V plane */
    }

    /* restore the ptr basing on cbp */
    if(*u1_cbp_c == 0)
    {
        (*pv_mb_coeff_data) = ps_mb_coeff_data_dc;
    }
    else if(*u1_cbp_c == 1)
    {
        (*pv_mb_coeff_data) = ps_mb_coeff_data_ac;
    }

    return;
}

/**
*******************************************************************************
*
* @brief performs luma core coding when intra mode is i16x16
*
* @par Description:
*  If the current mb is to be coded as intra of mb type i16x16, the mb is first
*  predicted using one of i16x16 prediction filters, basing on the intra mode
*  chosen. Then, error is computed between the input blk and the estimated blk.
*  This error is transformed (hierarchical transform i.e., dct followed by hada-
*  -mard), quantized. The quantized coefficients are packed in scan order for
*  entropy coding.
*
* @param[in] ps_proc_ctxt
*  pointer to the current macro block context
*
* @returns u1_cbp_l
*  coded block pattern luma
*
* @remarks none
*
*******************************************************************************
*/

UWORD8 isvce_code_luma_intra_macroblock_16x16(isvce_process_ctxt_t *ps_proc)
{
    buffer_container_t s_src;
    buffer_container_t s_pred;
    buffer_container_t s_recon;
    buffer_container_t s_res;
    buffer_container_t s_quant_coeffs;

    /*Cntrol signal for itrans*/
    UWORD32 u4_cntrl;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    quant_params_t *ps_qp_params = ps_proc->ps_qp_params[0];
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    inter_pred_fxns_t *ps_inter_pred_fxns = &ps_isa_dependent_fxns->s_inter_pred_fxns;
    iq_it_res_rec_constants_t s_iq_it_res_rec_constants = {
        .pu2_iscal_mat = ps_qp_params->pu2_iscale_mat,
        .pu2_weigh_mat = ps_qp_params->pu2_weigh_mat,
        .u4_qp_div_6 = ps_qp_params->u1_qp_div};

    UWORD8 *pu1_pred_mb = NULL;
    WORD16 *pi2_res_mb = ps_proc->pi2_res_buf_intra_4x4;
    WORD32 i4_pred_stride = ps_proc->i4_pred_strd;
    WORD32 i4_res_strd = ps_proc->i4_res_strd;
    UWORD8 u1_intra_mode = ps_proc->u1_l_i16_mode;
    UWORD32 au4_nnz[5] = {0};
    UWORD8 u1_cbp_l = 0;
    UWORD8 *pu1_nnz = (UWORD8 *) au4_nnz;
    /* pointer to packed mb coeff data */
    void **pv_mb_coeff_data = &(ps_proc->pv_mb_coeff_data);

    if(u1_intra_mode == PLANE_I16x16)
    {
        pu1_pred_mb = ps_proc->pu1_pred_mb_intra_16x16_plane;
    }
    else
    {
        pu1_pred_mb = ps_proc->pu1_pred_mb_intra_16x16;
    }

    s_src = ps_proc->s_src_buf_props.as_component_bufs[Y];
    s_recon = ps_proc->s_rec_buf_props.as_component_bufs[Y];
    s_pred.pv_data = pu1_pred_mb;
    s_pred.i4_data_stride = i4_pred_stride;
    s_quant_coeffs.pv_data = pi2_res_mb;
    s_quant_coeffs.i4_data_stride = i4_res_strd;

    s_res = ps_codec->s_svc_ilp_data.ps_residual_bufs[ps_proc->u1_spatial_layer_id]
                .as_component_bufs[Y];
    s_res.pv_data = ((WORD16 *) s_res.pv_data) + ps_proc->i4_mb_x * MB_SIZE +
                    ps_proc->i4_mb_y * MB_SIZE * s_res.i4_data_stride;

    /********************************************************/
    /*  error estimation,                                   */
    /*  transform                                           */
    /*  quantization                                        */
    /********************************************************/
    isvce_luma_16x16_resi_trans_dctrans_quant(
        &s_src, &s_pred, &s_quant_coeffs, &ps_proc->ps_mb_res_buf->as_component_bufs[Y],
        ps_isa_dependent_fxns, ps_qp_params->pu2_scale_mat, ps_qp_params->pu2_thres_mat, pu1_nnz,
        ps_qp_params->u1_qbits, ps_qp_params->u4_dead_zone, ENABLE_DC_TRANSFORM,
        ps_proc->ps_mb_info->u1_residual_prediction_flag);

    /********************************************************/
    /*  pack coeff data for entropy coding                  */
    /********************************************************/
    isvce_pack_l_mb_i16(pi2_res_mb, pv_mb_coeff_data, i4_res_strd, &u1_cbp_l, pu1_nnz, &u4_cntrl);

    /********************************************************/
    /*  ierror estimation,                                  */
    /*  itransform                                          */
    /*  iquantization                                       */
    /********************************************************/
    /*
     *if refernce frame is not to be computed
     *we only need the right and bottom border 4x4 blocks to predict next intra
     *blocks, hence only compute them
     */
    if(!ps_proc->u4_compute_recon)
    {
        u4_cntrl &= 0x111F8000;
    }

    if(u4_cntrl)
    {
        isvce_luma_16x16_idctrans_iquant_itrans_recon(
            &s_quant_coeffs, &s_pred, &s_recon, &s_res,
            &ps_proc->ps_mb_res_buf->as_component_bufs[Y], &s_iq_it_res_rec_constants,
            ps_isa_dependent_fxns, ps_proc->pv_scratch_buff, u4_cntrl, ENABLE_DC_TRANSFORM, 0);
    }
    else
    {
        ps_inter_pred_fxns->pf_inter_pred_luma_copy(pu1_pred_mb, (UWORD8 *) s_recon.pv_data,
                                                    i4_pred_stride, s_recon.i4_data_stride, MB_SIZE,
                                                    MB_SIZE, NULL, 0);
    }

    return (u1_cbp_l);
}

/**
*******************************************************************************
*
* @brief performs luma core coding when intra mode is i4x4
*
* @par Description:
*  If the current mb is to be coded as intra of mb type i4x4, the mb is first
*  predicted using one of i4x4 prediction filters, basing on the intra mode
*  chosen. Then, error is computed between the input blk and the estimated blk.
*  This error is dct transformed and quantized. The quantized coefficients are
*  packed in scan order for entropy coding.
*
* @param[in] ps_proc_ctxt
*  pointer to the current macro block context
*
* @returns u1_cbp_l
*  coded block pattern luma
*
* @remarks
*  The traversal of 4x4 subblocks in the 16x16 macroblock is as per the scan
*order mentioned in h.264 specification
*
*******************************************************************************
*/
UWORD8 isvce_code_luma_intra_macroblock_4x4(isvce_process_ctxt_t *ps_proc)
{
    buffer_container_t s_src;
    buffer_container_t s_pred;
    buffer_container_t s_recon;
    buffer_container_t s_res;
    buffer_container_t s_res_pred;
    buffer_container_t s_quant_coeffs;

    /* pointer to neighbors: left, top, top-left */
    UWORD8 *pu1_mb_a;
    UWORD8 *pu1_mb_b;
    UWORD8 *pu1_mb_c;
    UWORD8 *pu1_mb_d;
    WORD32 i4_ngbr_avbl;
    UWORD8 u1_nnz;
    UWORD32 u4_nnz_cnt;
    /* significant coefficient map */
    UWORD32 u4_s_map;
    /*Dummy variable for 4x4 trans fucntion*/
    WORD16 i2_dc_dummy;
    UWORD32 i, b8, b4, u1_blk_x, u1_blk_y, u1_pix_x, u1_pix_y, coeff_cnt, mask;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    quant_params_t *ps_qp_params = ps_proc->ps_qp_params[0];
    /* pointer to packed mb coeff data */
    tu_sblk_coeff_data_t *ps_mb_coeff_data, *ps_mb_coeff_data_b8;
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;
    inter_pred_fxns_t *ps_inter_pred_fxns = &ps_isa_dependent_fxns->s_inter_pred_fxns;
    resi_trans_quant_constants_t s_resi_trans_quant_constants = {
        .pu2_scale_matrix = ps_qp_params->pu2_scale_mat,
        .pu2_threshold_matrix = ps_qp_params->pu2_thres_mat,
        .u4_qbits = ps_qp_params->u1_qbits,
        .u4_round_factor = ps_qp_params->u4_dead_zone};
    iq_it_res_rec_constants_t s_iq_it_res_rec_constants = {
        .pu2_iscal_mat = ps_qp_params->pu2_iscale_mat,
        .pu2_weigh_mat = ps_qp_params->pu2_weigh_mat,
        .u4_qp_div_6 = ps_qp_params->u1_qp_div};

    UWORD8 *pu1_pred_mb = ps_proc->pu1_pred_mb;
    WORD16 *pi2_res_mb = ps_proc->pi2_res_buf_intra_4x4;
    WORD32 i4_pred_stride = ps_proc->i4_pred_strd;
    UWORD8 u1_intra_mode = ps_proc->u1_l_i16_mode;
    UWORD8 *pu1_ngbr_pels_i4 = ps_proc->au1_ngbr_pels;
    UWORD8 u1_cbp_l = 0;
    /* pointer to packed mb coeff data */
    void **pv_mb_coeff_data = &(ps_proc->pv_mb_coeff_data);
    const UWORD8 *pu1_scan_order = gu1_luma_scan_order;
    UWORD8 u1_resi_trans_fxn_idx = isvc_get_resi_trans_quant_variant_idx(0);
    UWORD8 u1_iq_it_recon_fxn_idx = isvc_get_iq_it_recon_variant_idx(1, 0);

    s_src = ps_proc->s_src_buf_props.as_component_bufs[Y];
    s_recon = ps_proc->s_rec_buf_props.as_component_bufs[Y];
    s_pred.pv_data = pu1_pred_mb;
    s_pred.i4_data_stride = i4_pred_stride;
    s_quant_coeffs.pv_data = pi2_res_mb;
    s_quant_coeffs.i4_data_stride = 4;

    /* Process 16 4x4 lum sub-blocks of the MB in scan order */
    for(b8 = 0; b8 < 4; b8++)
    {
        u1_blk_x = GET_BLK_RASTER_POS_X(b8) << 3;
        u1_blk_y = GET_BLK_RASTER_POS_Y(b8) << 3;

        /* if in case cbp for the 8x8 block is zero, send no residue */
        ps_mb_coeff_data_b8 = *pv_mb_coeff_data;

        for(b4 = 0; b4 < 4; b4++)
        {
            /* index of pel in MB */
            u1_pix_x = u1_blk_x + (GET_SUB_BLK_RASTER_POS_X(b4) << 2);
            u1_pix_y = u1_blk_y + (GET_SUB_BLK_RASTER_POS_Y(b4) << 2);

            /* Initialize source and reference pointers */
            s_src = ps_proc->s_src_buf_props.as_component_bufs[Y];
            s_recon = ps_proc->s_rec_buf_props.as_component_bufs[Y];
            s_src.pv_data = ((UWORD8 *) s_src.pv_data) + u1_pix_x + u1_pix_y * s_src.i4_data_stride;
            s_recon.pv_data =
                ((UWORD8 *) s_recon.pv_data) + u1_pix_x + u1_pix_y * s_recon.i4_data_stride;

            s_res = ps_codec->s_svc_ilp_data.ps_residual_bufs[ps_proc->u1_spatial_layer_id]
                        .as_component_bufs[Y];
            s_res.pv_data = ((WORD16 *) s_res.pv_data) + ps_proc->i4_mb_x * MB_SIZE +
                            ps_proc->i4_mb_y * MB_SIZE * s_res.i4_data_stride;
            s_res.pv_data = ((WORD16 *) s_res.pv_data) + u1_pix_x + u1_pix_y * s_res.i4_data_stride;

            s_res_pred = ps_proc->ps_mb_res_buf->as_component_bufs[Y];
            s_res_pred.pv_data =
                ((WORD16 *) s_res_pred.pv_data) + u1_pix_x + u1_pix_y * s_res_pred.i4_data_stride;

            /* pointer to left of ref macro block */
            pu1_mb_a = ((UWORD8 *) s_recon.pv_data) - 1;
            /* pointer to top of ref macro block */
            pu1_mb_b = ((UWORD8 *) s_recon.pv_data) - s_recon.i4_data_stride;
            /* pointer to topright of ref macro block */
            pu1_mb_c = pu1_mb_b + 4;
            /* pointer to topleft macro block */
            pu1_mb_d = pu1_mb_b - 1;

            /* compute neighbor availability */
            i4_ngbr_avbl = ps_proc->au1_ngbr_avbl_4x4_subblks[(b8 << 2) + b4];

            /* sub block intra mode */
            u1_intra_mode = ps_proc->au1_intra_luma_mb_4x4_modes[(b8 << 2) + b4];

            /********************************************************/
            /* gather prediction pels from neighbors for prediction */
            /********************************************************/
            /* left pels */
            if(i4_ngbr_avbl & LEFT_MB_AVAILABLE_MASK)
            {
                for(i = 0; i < 4; i++)
                    pu1_ngbr_pels_i4[4 - 1 - i] = pu1_mb_a[i * s_recon.i4_data_stride];
            }
            else
            {
                memset(pu1_ngbr_pels_i4, 0, 4);
            }

            /* top pels */
            if(i4_ngbr_avbl & TOP_MB_AVAILABLE_MASK)
            {
                memcpy(pu1_ngbr_pels_i4 + 4 + 1, pu1_mb_b, 4);
            }
            else
            {
                memset(pu1_ngbr_pels_i4 + 5, 0, 4);
            }
            /* top left pels */
            if(i4_ngbr_avbl & TOP_LEFT_MB_AVAILABLE_MASK)
            {
                pu1_ngbr_pels_i4[4] = *pu1_mb_d;
            }
            else
            {
                pu1_ngbr_pels_i4[4] = 0;
            }
            /* top right pels */
            if(i4_ngbr_avbl & TOP_RIGHT_MB_AVAILABLE_MASK)
            {
                memcpy(pu1_ngbr_pels_i4 + 8 + 1, pu1_mb_c, 4);
            }
            else if(i4_ngbr_avbl & TOP_MB_AVAILABLE_MASK)
            {
                memset(pu1_ngbr_pels_i4 + 8 + 1, pu1_ngbr_pels_i4[8], 4);
            }

            /********************************************************/
            /*  prediction                                          */
            /********************************************************/
            (ps_codec->apf_intra_pred_4_l)[u1_intra_mode](pu1_ngbr_pels_i4, pu1_pred_mb, 0,
                                                          i4_pred_stride, i4_ngbr_avbl);

            /********************************************************/
            /*  error estimation,                                   */
            /*  transform                                           */
            /*  quantization                                        */
            /********************************************************/
            ps_enc_loop_fxns->apf_resi_trans_quant_4x4[u1_resi_trans_fxn_idx](
                &s_src, &s_pred, &s_quant_coeffs, &s_res_pred, &s_resi_trans_quant_constants,
                &u1_nnz, &i2_dc_dummy, 0);

            /********************************************************/
            /*  pack coeff data for entropy coding                  */
            /********************************************************/
            ps_mb_coeff_data = *pv_mb_coeff_data;

            /* write number of non zero coefficients */
            ps_mb_coeff_data->i4_sig_map_nnz = u1_nnz;

            if(u1_nnz)
            {
                for(u4_nnz_cnt = 0, coeff_cnt = 0, mask = 1, u4_s_map = 0; u4_nnz_cnt < u1_nnz;
                    coeff_cnt++)
                {
                    if(pi2_res_mb[pu1_scan_order[coeff_cnt]])
                    {
                        /* write residue */
                        ps_mb_coeff_data->ai2_residue[u4_nnz_cnt++] =
                            pi2_res_mb[pu1_scan_order[coeff_cnt]];
                        u4_s_map |= mask;
                    }
                    mask <<= 1;
                }
                /* write significant coeff map */
                ps_mb_coeff_data->i4_sig_map_nnz |= (u4_s_map << 16);

                /* update ptr to coeff data */
                (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue + ALIGN2(u4_nnz_cnt);

                /* cbp */
                u1_cbp_l |= (1 << b8);
            }
            else
            {
                (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue;
            }

            /********************************************************/
            /*  ierror estimation,                                  */
            /*  itransform                                          */
            /*  iquantization                                       */
            /********************************************************/
            if(u1_nnz)
            {
                buffer_container_t s_src = s_quant_coeffs;

                /* Tx blk coeffs are stored blk by blk */
                /* Hence, in order to access rows of each Tx blk, one needs to stride of
                 * TxxSize */
                s_src.i4_data_stride = 4;

                ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4[u1_iq_it_recon_fxn_idx](
                    &s_src, &s_pred, &s_res_pred, &s_res, &s_recon, &s_iq_it_res_rec_constants,
                    (WORD16 *) ps_proc->pv_scratch_buff, s_src.pv_data, 0, 0);
            }
            else
            {
                ps_inter_pred_fxns->pf_inter_pred_luma_copy(
                    (UWORD8 *) s_pred.pv_data, (UWORD8 *) s_recon.pv_data, s_pred.i4_data_stride,
                    s_recon.i4_data_stride, BLK_SIZE, BLK_SIZE, NULL, 0);
            }
        }

        /* if the 8x8 block has no residue, nothing needs to be sent to entropy */
        if(!(u1_cbp_l & (1 << b8)))
        {
            *pv_mb_coeff_data = ps_mb_coeff_data_b8;
        }
    }

    return (u1_cbp_l);
}

/**
*******************************************************************************
*
* @brief performs luma core coding when intra mode is i4x4
*
* @par Description:
*  If the current mb is to be coded as intra of mb type i4x4, the mb is first
*  predicted using one of i4x4 prediction filters, basing on the intra mode
*  chosen. Then, error is computed between the input blk and the estimated blk.
*  This error is dct transformed and quantized. The quantized coefficients are
*  packed in scan order for entropy coding.
*
* @param[in] ps_proc_ctxt
*  pointer to the current macro block context
*
* @returns u1_cbp_l
*  coded block pattern luma
*
* @remarks
*  The traversal of 4x4 subblocks in the 16x16 macroblock is as per the scan
*order mentioned in h.264 specification
*
*******************************************************************************
*/
UWORD8 isvce_code_luma_intra_macroblock_4x4_rdopt_on(isvce_process_ctxt_t *ps_proc)
{
    /* pointer to packed mb coeff data */
    tu_sblk_coeff_data_t *ps_mb_coeff_data, *ps_mb_coeff_data_b8;

    UWORD32 u4_nnz_cnt;
    /* significant coefficient map */
    UWORD32 u4_s_map;
    UWORD32 b8, b4, coeff_cnt, mask;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    inter_pred_fxns_t *ps_inter_pred_fxns = &ps_isa_dependent_fxns->s_inter_pred_fxns;

    UWORD8 *pu1_ref_mb_intra_4x4 = ps_proc->pu1_ref_mb_intra_4x4;
    UWORD8 *pu1_rec_mb = ((UWORD8 *) ps_proc->s_rec_buf_props.as_component_bufs[0].pv_data);
    WORD16 *pi2_res_mb = ps_proc->pi2_res_buf_intra_4x4;
    WORD32 i4_rec_strd = ps_proc->s_rec_buf_props.as_component_bufs[0].i4_data_stride;
    UWORD8 *pu1_nnz = (UWORD8 *) ps_proc->au4_nnz_intra_4x4;
    UWORD8 u1_cbp_l = 0;
    void **pv_mb_coeff_data = &(ps_proc->pv_mb_coeff_data);
    const UWORD8 *pu1_scan_order = gu1_luma_scan_order;

    /* Process 16 4x4 lum sub-blocks of the MB in scan order */
    for(b8 = 0; b8 < 4; b8++)
    {
        /* if in case cbp for the 8x8 block is zero, send no residue */
        ps_mb_coeff_data_b8 = *pv_mb_coeff_data;

        for(b4 = 0; b4 < 4; b4++, pu1_nnz++, pi2_res_mb += MB_SIZE)
        {
            /********************************************************/
            /*  pack coeff data for entropy coding                  */
            /********************************************************/
            ps_mb_coeff_data = *pv_mb_coeff_data;

            /* write number of non zero coefficients */
            ps_mb_coeff_data->i4_sig_map_nnz = *pu1_nnz;

            if(*pu1_nnz)
            {
                for(u4_nnz_cnt = 0, coeff_cnt = 0, mask = 1, u4_s_map = 0; u4_nnz_cnt < *pu1_nnz;
                    coeff_cnt++)
                {
                    if(pi2_res_mb[pu1_scan_order[coeff_cnt]])
                    {
                        /* write residue */
                        ps_mb_coeff_data->ai2_residue[u4_nnz_cnt++] =
                            pi2_res_mb[pu1_scan_order[coeff_cnt]];
                        u4_s_map |= mask;
                    }
                    mask <<= 1;
                }
                /* write significant coeff map */
                ps_mb_coeff_data->i4_sig_map_nnz |= (u4_s_map << 16);

                /* update ptr to coeff data */
                (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue + ALIGN2(u4_nnz_cnt);

                /* cbp */
                u1_cbp_l |= (1 << b8);
            }
            else
            {
                (*pv_mb_coeff_data) = ps_mb_coeff_data->ai2_residue;
            }
        }

        /* if the 8x8 block has no residue, nothing needs to be sent to entropy */
        if(!(u1_cbp_l & (1 << b8)))
        {
            *pv_mb_coeff_data = ps_mb_coeff_data_b8;
        }
    }

    ps_inter_pred_fxns->pf_inter_pred_luma_copy(pu1_ref_mb_intra_4x4, pu1_rec_mb, MB_SIZE,
                                                i4_rec_strd, MB_SIZE, MB_SIZE, NULL, 0);

    return (u1_cbp_l);
}

/**
*******************************************************************************
*
* @brief performs chroma core coding for intra macro blocks
*
* @par Description:
*  If the current MB is to be intra coded with mb type chroma I8x8, the MB is
*  first predicted using intra 8x8 prediction filters. The predicted data is
*  compared with the input for error and the error is transformed. The DC
*  coefficients of each transformed sub blocks are further transformed using
*  Hadamard transform. The resulting coefficients are quantized, packed and sent
*  for entropy coding.
*
* @param[in] ps_proc_ctxt
*  pointer to the current macro block context
*
* @returns u1_cbp_c
*  coded block pattern chroma
*
* @remarks
*  The traversal of 4x4 subblocks in the 8x8 macroblock is as per the scan order
*  mentioned in h.264 specification
*
*******************************************************************************
*/
UWORD8 isvce_code_chroma_intra_macroblock_8x8(isvce_process_ctxt_t *ps_proc)
{
    buffer_container_t s_src;
    buffer_container_t s_pred;
    buffer_container_t s_recon;
    buffer_container_t s_res;
    buffer_container_t s_res_pred;
    buffer_container_t s_quant_coeffs;

    /* Control signal for inverse transform */
    UWORD32 u4_cntrl;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    quant_params_t *ps_qp_params = ps_proc->ps_qp_params[1];
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    iq_it_res_rec_constants_t s_iq_it_res_rec_constants = {
        .pu2_iscal_mat = ps_qp_params->pu2_iscale_mat,
        .pu2_weigh_mat = ps_qp_params->pu2_weigh_mat,
        .u4_qp_div_6 = ps_qp_params->u1_qp_div};

    UWORD8 *pu1_pred_mb = NULL;
    WORD16 *pi2_res_mb = ps_proc->pi2_res_buf;
    WORD32 i4_pred_stride = ps_proc->i4_pred_strd;
    WORD32 i4_res_strd = ps_proc->i4_res_strd;
    UWORD8 u1_intra_mode = ps_proc->u1_c_i8_mode;
    UWORD8 u1_cbp_c = 0;
    UWORD8 au1_nnz[2 * (NUM_4x4_IN_8x8 + 1)] = {0};
    /* pointer to packed mb coeff data */
    void **pv_mb_coeff_data = &(ps_proc->pv_mb_coeff_data);
    /* See if we need to swap U and V plances for entropy */
    UWORD32 u4_swap_uv = (ps_codec->s_cfg.e_inp_color_fmt == IV_YUV_420SP_VU);

    if(PLANE_CH_I8x8 == u1_intra_mode)
    {
        pu1_pred_mb = ps_proc->pu1_pred_mb_intra_chroma_plane;
    }
    else
    {
        pu1_pred_mb = ps_proc->pu1_pred_mb_intra_chroma;
    }

    s_src = ps_proc->s_src_buf_props.as_component_bufs[UV];
    s_recon = ps_proc->s_rec_buf_props.as_component_bufs[UV];
    s_pred.pv_data = pu1_pred_mb;
    s_pred.i4_data_stride = i4_pred_stride;
    s_quant_coeffs.pv_data = pi2_res_mb;
    s_quant_coeffs.i4_data_stride = i4_res_strd;

    s_res = ps_codec->s_svc_ilp_data.ps_residual_bufs[ps_proc->u1_spatial_layer_id]
                .as_component_bufs[UV];
    s_res.pv_data = ((WORD16 *) s_res.pv_data) + ps_proc->i4_mb_x * MB_SIZE +
                    ps_proc->i4_mb_y * (MB_SIZE / 2) * s_res.i4_data_stride;

    s_res_pred = ps_proc->ps_mb_res_buf->as_component_bufs[U];

    /********************************************************/
    /*  error estimation,                                   */
    /*  transform                                           */
    /*  quantization                                        */
    /********************************************************/
    isvce_chroma_8x8_resi_trans_dctrans_quant(
        &s_src, &s_pred, &s_quant_coeffs, &s_res_pred, ps_isa_dependent_fxns,
        ps_qp_params->pu2_scale_mat, ps_qp_params->pu2_thres_mat, au1_nnz, ps_qp_params->u1_qbits,
        ps_qp_params->u4_dead_zone, 0);

    /********************************************************/
    /*  pack coeff data for entropy coding                  */
    /********************************************************/
    isvce_pack_c_mb(pi2_res_mb, pv_mb_coeff_data, i4_res_strd, &u1_cbp_c, au1_nnz,
                    ps_codec->u4_thres_resi, &u4_cntrl, u4_swap_uv);

    /********************************************************/
    /*  ierror estimation,                                  */
    /*  itransform                                          */
    /*  iquantization                                       */
    /********************************************************/
    isvce_chroma_8x8_idctrans_iquant_itrans_recon(
        &s_quant_coeffs, &s_pred, &s_recon, &s_res, &s_res_pred, &s_iq_it_res_rec_constants,
        ps_isa_dependent_fxns, ps_proc->pv_scratch_buff, u4_cntrl, 0);

    memcpy(ps_proc->au1_chroma_nnz, au1_nnz, sizeof(ps_proc->au1_chroma_nnz));

    return (u1_cbp_c);
}

/**
*******************************************************************************
*
* @brief performs luma core coding when  mode is inter
*
* @par Description:
*  If the current mb is to be coded as inter the mb is predicted based on the
*  sub mb partitions and corresponding motion vectors generated by ME. Then,
*  error is computed between the input blk and the estimated blk. This error is
*  transformed, quantized. The quantized coefficients are packed in scan order
*  for entropy coding
*
* @param[in] ps_proc_ctxt
*  pointer to the current macro block context
*
* @returns u1_cbp_l
*  coded block pattern luma
*
* @remarks none
*
*******************************************************************************
*/

UWORD8 isvce_code_luma_inter_macroblock_16x16(isvce_process_ctxt_t *ps_proc)
{
    buffer_container_t s_src;
    buffer_container_t s_pred;
    buffer_container_t s_recon;
    buffer_container_t s_res;
    buffer_container_t s_res_pred;
    buffer_container_t s_quant_coeffs;

    /*Control signal of itrans*/
    UWORD32 u4_cntrl;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    quant_params_t *ps_qp_params = ps_proc->ps_qp_params[0];
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    iq_it_res_rec_constants_t s_iq_it_res_rec_constants = {
        .pu2_iscal_mat = ps_qp_params->pu2_iscale_mat,
        .pu2_weigh_mat = ps_qp_params->pu2_weigh_mat,
        .u4_qp_div_6 = ps_qp_params->u1_qp_div};

    WORD16 *pi2_res_mb = ps_proc->pi2_res_buf_intra_4x4;
    WORD32 i4_res_strd = ps_proc->i4_res_strd;
    UWORD8 u1_cbp_l = 0;
    UWORD8 *pu1_nnz = (UWORD8 *) ps_proc->au4_nnz;
    /* pointer to packed mb coeff data */
    void **pv_mb_coeff_data = &(ps_proc->pv_mb_coeff_data);

    ps_proc->au4_nnz[0] = 0;
    ps_proc->au4_nnz[1] = 0;
    ps_proc->au4_nnz[2] = 0;
    ps_proc->au4_nnz[3] = 0;
    ps_proc->au4_nnz[4] = 0;

    /********************************************************/
    /*  prediction                                          */
    /********************************************************/
    isvce_motion_comp_luma(ps_proc, &s_pred);

    s_src = ps_proc->s_src_buf_props.as_component_bufs[0];
    s_recon = ps_proc->s_rec_buf_props.as_component_bufs[0];
    s_quant_coeffs.pv_data = pi2_res_mb;
    s_quant_coeffs.i4_data_stride = i4_res_strd;

    s_res = ps_codec->s_svc_ilp_data.ps_residual_bufs[ps_proc->u1_spatial_layer_id]
                .as_component_bufs[Y];
    s_res.pv_data = ((WORD16 *) s_res.pv_data) + ps_proc->i4_mb_x * MB_SIZE +
                    ps_proc->i4_mb_y * MB_SIZE * s_res.i4_data_stride;

    s_res_pred = ps_proc->ps_mb_res_buf->as_component_bufs[Y];

    /********************************************************/
    /*  error estimation,                                   */
    /*  transform                                           */
    /*  quantization                                        */
    /********************************************************/
    if(ps_proc->u4_min_sad_reached == 0 || ps_proc->u4_min_sad != 0)
    {
        isvce_luma_16x16_resi_trans_dctrans_quant(
            &s_src, &s_pred, &s_quant_coeffs, &s_res_pred, ps_isa_dependent_fxns,
            ps_qp_params->pu2_scale_mat, ps_qp_params->pu2_thres_mat, pu1_nnz,
            ps_qp_params->u1_qbits, ps_qp_params->u4_dead_zone, DISABLE_DC_TRANSFORM,
            ps_proc->ps_mb_info->u1_residual_prediction_flag);

        /********************************************************/
        /*  pack coeff data for entropy coding                  */
        /********************************************************/
        isvce_pack_l_mb(pi2_res_mb, pv_mb_coeff_data, i4_res_strd, &u1_cbp_l, pu1_nnz,
                        ps_codec->u4_thres_resi, &u4_cntrl);
    }
    else
    {
        u1_cbp_l = 0;
        u4_cntrl = 0;
    }

    /********************************************************/
    /*  ierror estimation,                                  */
    /*  itransform                                          */
    /*  iquantization                                       */
    /********************************************************/

    /*If the frame is not to be used for P frame reference or dumping recon
     * we only will use the reocn for only predicting intra Mbs
     * THis will need only right and bottom edge 4x4 blocks recon
     * Hence we selectively enable them using control signal(including DC)
     */
    if(ps_proc->u4_compute_recon != 1)
    {
        u4_cntrl &= 0x111F0000;
    }

    isvce_luma_16x16_idctrans_iquant_itrans_recon(
        &s_quant_coeffs, &s_pred, &s_recon, &s_res, &s_res_pred, &s_iq_it_res_rec_constants,
        ps_isa_dependent_fxns, ps_proc->pv_scratch_buff, u4_cntrl, DISABLE_DC_TRANSFORM,
        ps_proc->ps_mb_info->u1_residual_prediction_flag);

    return (u1_cbp_l);
}

/**
*******************************************************************************
*
* @brief performs chroma core coding for inter macro blocks
*
* @par Description:
*  If the current mb is to be coded as inter predicted mb,based on the sub mb
*partitions and corresponding motion vectors generated by ME  ,prediction is
*done. Then, error is computed between the input blk and the estimated blk. This
*error is transformed , quantized. The quantized coefficients are packed in scan
*order for entropy coding.
*
* @param[in] ps_proc_ctxt
*  pointer to the current macro block context
*
* @returns u1_cbp_l
*  coded block pattern chroma
*
* @remarks none
*
*******************************************************************************
*/
UWORD8 isvce_code_chroma_inter_macroblock_8x8(isvce_process_ctxt_t *ps_proc)
{
    buffer_container_t s_src;
    buffer_container_t s_pred;
    buffer_container_t s_recon;
    buffer_container_t s_res;
    buffer_container_t s_res_pred;
    buffer_container_t s_quant_coeffs;

    /*Control signal for inverse transform*/
    UWORD32 u4_cntrl;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    quant_params_t *ps_qp_params = ps_proc->ps_qp_params[1];
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    iq_it_res_rec_constants_t s_iq_it_res_rec_constants = {
        .pu2_iscal_mat = ps_qp_params->pu2_iscale_mat,
        .pu2_weigh_mat = ps_qp_params->pu2_weigh_mat,
        .u4_qp_div_6 = ps_qp_params->u1_qp_div};

    WORD16 *pi2_res_mb = ps_proc->pi2_res_buf;
    WORD32 i4_res_strd = ps_proc->i4_res_strd;
    UWORD8 u1_cbp_c = 0;
    UWORD8 au1_nnz[2 * (NUM_4x4_IN_8x8 + 1)] = {0};
    /* pointer to packed mb coeff data */
    void **pv_mb_coeff_data = &(ps_proc->pv_mb_coeff_data);
    /*See if we need to swap U and V plances for entropy*/
    UWORD32 u4_swap_uv = ps_codec->s_cfg.e_inp_color_fmt == IV_YUV_420SP_VU;

    isvce_motion_comp_chroma(ps_proc, &s_pred);

    s_src = ps_proc->s_src_buf_props.as_component_bufs[UV];
    s_recon = ps_proc->s_rec_buf_props.as_component_bufs[UV];
    s_quant_coeffs.pv_data = pi2_res_mb;
    s_quant_coeffs.i4_data_stride = i4_res_strd;

    s_res = ps_codec->s_svc_ilp_data.ps_residual_bufs[ps_proc->u1_spatial_layer_id]
                .as_component_bufs[UV];
    s_res.pv_data = ((WORD16 *) s_res.pv_data) + ps_proc->i4_mb_x * MB_SIZE +
                    ps_proc->i4_mb_y * (MB_SIZE / 2) * s_res.i4_data_stride;

    s_res_pred = ps_proc->ps_mb_res_buf->as_component_bufs[UV];

    /********************************************************/
    /*  error estimation,                                   */
    /*  transform                                           */
    /*  quantization                                        */
    /********************************************************/
    isvce_chroma_8x8_resi_trans_dctrans_quant(
        &s_src, &s_pred, &s_quant_coeffs, &s_res_pred, ps_isa_dependent_fxns,
        ps_qp_params->pu2_scale_mat, ps_qp_params->pu2_thres_mat, au1_nnz, ps_qp_params->u1_qbits,
        ps_qp_params->u4_dead_zone, ps_proc->ps_mb_info->u1_residual_prediction_flag);

    /********************************************************/
    /*  pack coeff data for entropy coding                  */
    /********************************************************/
    isvce_pack_c_mb(pi2_res_mb, pv_mb_coeff_data, i4_res_strd, &u1_cbp_c, au1_nnz,
                    ps_codec->u4_thres_resi, &u4_cntrl, u4_swap_uv);

    /********************************************************/
    /*  ierror estimation,                                  */
    /*  itransform                                          */
    /*  iquantization                                       */
    /********************************************************/

    /* If the frame is not to be used for P frame reference or dumping recon
     * we only will use the reocn for only predicting intra Mbs
     * THis will need only right and bottom edge 4x4 blocks recon
     * Hence we selectively enable them using control signal(including DC)
     */
    if(!ps_proc->u4_compute_recon)
    {
        u4_cntrl &= 0x7700C000;
    }

    isvce_chroma_8x8_idctrans_iquant_itrans_recon(
        &s_quant_coeffs, &s_pred, &s_recon, &s_res, &s_res_pred, &s_iq_it_res_rec_constants,
        ps_isa_dependent_fxns, ps_proc->pv_scratch_buff, u4_cntrl,
        ps_proc->ps_mb_info->u1_residual_prediction_flag);

    memcpy(ps_proc->au1_chroma_nnz, au1_nnz, sizeof(ps_proc->au1_chroma_nnz));

    return (u1_cbp_c);
}
