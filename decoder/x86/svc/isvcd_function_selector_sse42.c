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
*  isvcd_function_selector_sse42.c
*
* @brief
*  Contains functions to initialize function pointers of codec context
*
* @author
*  Kishore
*
* @par List of Functions:
*  - isvcd_init_function_ptr_sse42
*
* @remarks
*  None
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System Include files */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* User Include files */
#include "isvcd_structs.h"
#include "ih264d_function_selector.h"

/**
 *******************************************************************************
 *
 * @brief Initialize the intra/inter/transform/deblk function pointers of
 * codec context
 *
 * @par Description: the current routine initializes the function pointers of
 * codec context based on x86_sse42 arch in use for svc
 *
 * @param[in] ps_svc_lyr_dec
 *  svc dec layer context pointer
 *
 * @returns  none
 *
 * @remarks none
 *
 *******************************************************************************
 */
void isvcd_init_function_ptr_sse42(svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    dec_struct_t *ps_codec = &ps_svc_lyr_dec->s_dec;
    residual_sampling_ctxt_t *ps_resd_samp_ctx;
    intra_sampling_ctxt_t *ps_intra_samp_ctxt;
    /* call default init func prt sse42 */
    ih264d_init_function_ptr_sse42(ps_codec);

    ps_resd_samp_ctx = (residual_sampling_ctxt_t *) ps_svc_lyr_dec->pv_residual_sample_ctxt;
    ps_intra_samp_ctxt = (intra_sampling_ctxt_t *) ps_svc_lyr_dec->pv_intra_sample_ctxt;

    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_4x4 =
        isvcd_iquant_itrans_residual_recon_4x4_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_4x4_dc =
        isvcd_iquant_itrans_residual_recon_4x4_dc_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_8x8 =
        isvcd_iquant_itrans_residual_recon_8x8_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_8x8_dc =
        isvcd_iquant_itrans_residual_recon_8x8_dc_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_chroma_4x4 =
        isvcd_iquant_itrans_residual_recon_chroma_4x4_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_chroma_4x4_dc =
        isvcd_iquant_itrans_residual_recon_chroma_4x4_dc_sse42;

    ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_4x4 = isvcd_iquant_itrans_residual_4x4_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_4x4_dc =
        isvcd_iquant_itrans_residual_4x4_dc_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_8x8 = isvcd_iquant_itrans_residual_8x8_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_8x8_dc =
        isvcd_iquant_itrans_residual_8x8_dc_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_chroma_4x4 =
        isvcd_iquant_itrans_residual_chroma_4x4_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_chroma_4x4_dc =
        isvcd_iquant_itrans_residual_chroma_4x4_dc_sse42;

    ps_svc_lyr_dec->pf_pred_residual_recon_luma_4x4 = isvcd_pred_residual_recon_4x4_sse42;
    ps_svc_lyr_dec->pf_pred_residual_recon_luma_8x8 = isvcd_pred_residual_recon_8x8_sse42;
    ps_svc_lyr_dec->pf_pred_residual_recon_luma_16x16 = isvcd_pred_residual_recon_16x16_sse42;
    ps_svc_lyr_dec->pf_pred_residual_recon_chroma_4x4 = isvcd_pred_residual_recon_chroma_4x4_sse42;
    ps_svc_lyr_dec->pf_pred_residual_recon_chroma_8x8 = isvcd_pred_residual_recon_chroma_8x8_sse42;

    ps_svc_lyr_dec->pf_residual_luma_4x4 = isvcd_residual_luma_4x4_sse42;
    ps_svc_lyr_dec->pf_residual_luma_8x8 = isvcd_residual_luma_8x8_sse42;
    ps_svc_lyr_dec->pf_residual_luma_16x16 = isvcd_residual_luma_16x16_sse42;
    ps_svc_lyr_dec->pf_residual_chroma_cb_cr_8x8 = isvcd_residual_chroma_cb_cr_8x8_sse42;

    ps_svc_lyr_dec->pf_iquant_itrans_luma_4x4 = isvcd_iquant_itrans_4x4_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_luma_4x4_dc = isvcd_iquant_itrans_4x4_dc_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_luma_8x8 = isvcd_iquant_itrans_8x8_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_luma_8x8_dc = isvcd_iquant_itrans_8x8_dc_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_chroma_4x4 = isvcd_iquant_itrans_chroma_4x4_sse42;
    ps_svc_lyr_dec->pf_iquant_itrans_chroma_4x4_dc = isvcd_iquant_itrans_chroma_4x4_dc_sse42;

    ps_intra_samp_ctxt->pf_interpolate_base_luma_dyadic = isvcd_interpolate_base_luma_dyadic_sse42;
    ps_intra_samp_ctxt->pf_interpolate_intra_base = isvcd_interpolate_intra_base_sse42;

    ps_intra_samp_ctxt->pf_vert_chroma_interpol[0] = isvcd_vert_interpol_chroma_dyadic_1_sse42;
    ps_intra_samp_ctxt->pf_vert_chroma_interpol[1] = isvcd_vert_interpol_chroma_dyadic_2_sse42;
    ps_intra_samp_ctxt->pf_vert_chroma_interpol[2] = isvcd_vert_interpol_chroma_dyadic_3_sse42;

    ps_intra_samp_ctxt->pf_horz_chroma_interpol[0] = isvcd_horz_interpol_chroma_dyadic_1_sse42;
    ps_intra_samp_ctxt->pf_horz_chroma_interpol[1] = isvcd_horz_interpol_chroma_dyadic_2_sse42;

    ps_resd_samp_ctx->pf_residual_luma_dyadic = isvcd_residual_luma_dyadic_sse42;
    ps_resd_samp_ctx->pf_interpolate_residual = isvcd_interpolate_residual_sse42;
    ps_resd_samp_ctx->pf_residual_reflayer_const_non_boundary_mb =
        isvcd_residual_reflayer_const_non_boundary_mb_sse42;

    return;
}
