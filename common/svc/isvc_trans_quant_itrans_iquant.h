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
 *  isvc_trans_quant.h
 *
 * @brief
 *  Contains declarations for forward and inverse transform paths for H264
 *
 * @author
 *  Ittiam
 *
 * @remarks
 *
 *******************************************************************************
 */

#ifndef _ISVC_TRANS_QUANT_ITRANS_IQUANT_H_
#define _ISVC_TRANS_QUANT_ITRANS_IQUANT_H_

#include <stdint.h>

#include "ih264_typedefs.h"
#include "ih264_debug.h"
#include "ih264_macros.h"
#include "isvc_macros.h"
#include "isvc_structs.h"

/* With and without residual_pred use */
#define NUM_RESI_TRANS_QUANT_VARIANTS 2

#define NUM_IQ_IT_RECON_VARIANTS 3

/* Structs */
typedef struct resi_trans_quant_constants_t
{
    const UWORD16 *pu2_scale_matrix;

    const UWORD16 *pu2_threshold_matrix;

    UWORD32 u4_qbits;

    UWORD32 u4_round_factor;
} resi_trans_quant_constants_t;

typedef struct iq_it_res_rec_constants_t
{
    const UWORD16 *pu2_iscal_mat;

    const UWORD16 *pu2_weigh_mat;

    UWORD32 u4_qp_div_6;
} iq_it_res_rec_constants_t;

/* Typedefs */
typedef void FT_RESI_TRANS_DCTRANS_QUANT(UWORD8 *pu1_src, UWORD8 *pu1_pred, WORD16 *pi2_out,
                                         WORD32 src_strd, WORD32 pred_strd, WORD32 dst_strd,
                                         const UWORD16 *pu2_scale_mat,
                                         const UWORD16 *pu2_thresh_mat, UWORD32 u4_qbit,
                                         UWORD32 u4_round_fact, UWORD8 *pu1_nnz);

typedef void FT_IDCTRANS_IQUANT_ITRANS_RECON(WORD16 *pi2_src, UWORD8 *pu1_pred, UWORD8 *pu1_out,
                                             WORD32 src_strd, WORD32 pred_strd, WORD32 out_strd,
                                             const UWORD16 *pu2_iscale_mat,
                                             const UWORD16 *pu2_weigh_mat, UWORD32 qp_div,
                                             UWORD32 pi4_cntrl, WORD32 *pi4_tmp);

typedef void FT_RESI_TRANS_QUANT(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                                 buffer_container_t *ps_out, buffer_container_t *ps_upsampled_res,
                                 resi_trans_quant_constants_t *ps_quant_constants, UWORD8 *pu1_nnz,
                                 WORD16 *pi2_dc_out, UWORD8 u1_use_upsampled_res);

typedef void FT_LUMA_16X16_RESI_TRANS_DCTRANS_QUANT(
    UWORD8 *pu1_src, UWORD8 *pu1_pred, WORD16 *pi2_out, WORD32 src_strd, WORD32 pred_strd,
    WORD32 dst_strd, const UWORD16 *pu2_scale_matrix, const UWORD16 *pu2_threshold_matrix,
    UWORD32 u4_qbits, UWORD32 u4_round_factor, UWORD8 *pu1_nnz, UWORD32 u4_dc_flag);

typedef void FT_CHROMA_8X8_RESI_TRANS_DCTRANS_QUANT(
    UWORD8 *pu1_src, UWORD8 *pu1_pred, WORD16 *pi2_out, WORD32 src_strd, WORD32 pred_strd,
    WORD32 dst_strd, const UWORD16 *pu2_scale_matrix, const UWORD16 *pu2_threshold_matrix,
    UWORD32 u4_qbits, UWORD32 u4_round_factor, UWORD8 *pu1_nnz);

typedef void FT_IQ_IT_RECON(buffer_container_t *ps_src, buffer_container_t *ps_pred,
                            buffer_container_t *ps_res_pred, buffer_container_t *ps_res,
                            buffer_container_t *ps_rec,
                            iq_it_res_rec_constants_t *ps_iq_it_res_rec_constants, WORD16 *pi2_tmp,
                            WORD16 *pi2_dc_src, WORD32 i4_iq_start_idx, UWORD8 u1_res_accumulate);

typedef void FT_LUMA_16X16_IDCTRANS_IQUANT_ITRANS_RECON(
    WORD16 *pi2_src, UWORD8 *pu1_pred, UWORD8 *pu1_out, WORD32 src_strd, WORD32 pred_strd,
    WORD32 out_strd, const UWORD16 *pu2_iscale_mat, const UWORD16 *pu2_weigh_mat, UWORD32 qp_div,
    UWORD32 pi4_cntrl, UWORD32 u4_dc_trans_flag, WORD32 *pi4_tmp);

typedef void FT_CHROMA_8X8_IDCTRANS_IQUANT_ITRANS_RECON(
    WORD16 *pi2_src, UWORD8 *pu1_pred, UWORD8 *pu1_out, WORD32 src_strd, WORD32 pred_strd,
    WORD32 out_strd, const UWORD16 *pu2_iscale_mat, const UWORD16 *pu2_weigh_mat, UWORD32 qp_div,
    UWORD32 pi4_cntrl, WORD32 *pi4_tmp);

typedef void FT_IHADAMARD_SCALING(WORD16 *pi2_src, WORD16 *pi2_out, const UWORD16 *pu2_iscal_mat,
                                  const UWORD16 *pu2_weigh_mat, UWORD32 u4_qp_div_6,
                                  WORD32 *pi4_tmp);

typedef void FT_HADAMARD_QUANT(WORD16 *pi2_src, WORD16 *pi2_dst,
                               resi_trans_quant_constants_t *ps_quant_constants, UWORD8 *pu1_nnz);

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/

extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_4x4;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_chroma_4x4;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_8x8;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_8x8;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_dc;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_dc;
extern FT_IQ_IT_RECON isvc_zcbf_iquant_itrans_recon_4x4;
extern FT_IQ_IT_RECON isvc_chroma_zcbf_iquant_itrans_recon_4x4;
extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_4x4;
extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_2x2_uv;
extern FT_HADAMARD_QUANT isvc_hadamard_quant_4x4;
extern FT_HADAMARD_QUANT isvc_hadamard_quant_2x2_uv;

/* A9 Declarations */
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_4x4_a9;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_chroma_4x4_a9;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_a9;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_8x8_a9;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_dc_a9;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_8x8_dc_a9;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_a9;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_dc_a9;
extern FT_LUMA_16X16_RESI_TRANS_DCTRANS_QUANT isvc_luma_16x16_resi_trans_dctrans_quant_a9;
extern FT_CHROMA_8X8_RESI_TRANS_DCTRANS_QUANT isvc_chroma_8x8_resi_trans_dctrans_quant_a9;
extern FT_LUMA_16X16_IDCTRANS_IQUANT_ITRANS_RECON isvc_luma_16x16_idctrans_iquant_itrans_recon_a9;
extern FT_CHROMA_8X8_IDCTRANS_IQUANT_ITRANS_RECON isvc_chroma_8x8_idctrans_iquant_itrans_recon_a9;
extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_4x4_a9;
extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_2x2_uv_a9;
extern FT_HADAMARD_QUANT isvc_hadamard_quant_4x4_a9;
extern FT_HADAMARD_QUANT isvc_hadamard_quant_2x2_uv_a9;

/* Av8 Declarations */
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_4x4_av8;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_chroma_4x4_av8;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_av8;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_8x8_av8;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_dc_av8;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_8x8_dc_av8;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_av8;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_dc_av8;
extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_4x4_av8;
extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_2x2_uv_av8;

/* NEON Declarations */
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_4x4_neon;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_4x4_with_residual_sub_neon;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_chroma_4x4_neon;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_chroma_4x4_with_residual_sub_neon;

/* SSSE3 Declarations */
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_ssse3;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_8x8_ssse3;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_dc_ssse3;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_8x8_dc_ssse3;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_dc_ssse3;
extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_4x4_ssse3;
extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_2x2_uv_ssse3;

/* SSSE42 Declarations */
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_4x4_sse42;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_4x4_with_res_pred_sse42;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_chroma_4x4_sse42;
extern FT_RESI_TRANS_QUANT isvc_resi_trans_quant_chroma_4x4_with_res_pred_sse42;

extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_sse42;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_res_4x4_sse42;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_res_4x4_with_res_acc_sse42;

extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_res_chroma_4x4_sse42;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_res_chroma_4x4_with_res_acc_sse42;

extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_dc_4x4_sse42;

extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_res_chroma_4x4_dc_sse42;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_res_chroma_4x4_dc_with_res_acc_sse42;

extern FT_IHADAMARD_SCALING ih264_ihadamard_scaling_4x4_sse42;

extern FT_HADAMARD_QUANT isvc_hadamard_quant_4x4_sse42;
extern FT_HADAMARD_QUANT isvc_hadamard_quant_2x2_uv_sse42;

/* NEON Declarations */
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_neon;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_with_res_output_neon;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_with_res_accumulate_neon;

extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_with_res_output_neon;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_with_res_accumulate_neon;

extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_4x4_dc_neon;

extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_dc_with_res_output_neon;
extern FT_IQ_IT_RECON isvc_iquant_itrans_recon_chroma_4x4_dc_with_res_accumulate_neon;

static FORCEINLINE UWORD8 isvc_get_resi_trans_quant_variant_idx(UWORD8 u1_use_upsampled_res)
{
    return u1_use_upsampled_res;
}

static FORCEINLINE UWORD8 isvc_get_iq_it_recon_variant_idx(UWORD8 u1_is_intra,
                                                           UWORD8 u1_res_accumulate)
{
    ASSERT(!((1 == u1_is_intra) && (1 == u1_res_accumulate)));

    return u1_is_intra * 2 + u1_res_accumulate;
}

static FORCEINLINE WORD16 isvc_get_residue(WORD16 i2_it_out, WORD16 i2_res_pred,
                                           UWORD8 u1_res_accumulate)
{
    return (u1_res_accumulate
                ? (CLIP3(-((WORD16) UINT8_MAX), ((WORD16) UINT8_MAX), i2_it_out + i2_res_pred))
                : (CLIP3(-((WORD16) UINT8_MAX), ((WORD16) UINT8_MAX), i2_it_out)));
}

#endif
