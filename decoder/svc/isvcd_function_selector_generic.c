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
 *  isvcd_function_selector_generic.c
 *
 * @brief
 *  Contains functions to initialize function pointers of codec context
 *
 * @author
 *  kishore
 *
 * @par List of Functions:
 *  - isvcd_init_function_ptr_generic()
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
#include "ih264_typedefs.h"
#include "iv.h"
#include "ivd.h"
#include "ih264_defs.h"
#include "ih264_size_defs.h"
#include "ih264_error.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include "ih264_inter_pred_filters.h"
#include "isvcd_iquant_itrans_residual_recon.h"
#include "isvcd_iquant_itrans_residual.h"
#include "isvcd_iquant_itrans.h"
#include "isvcd_pred_residual_recon.h"
#include "isvcd_structs.h"
#include "ih264d_function_selector.h"
#include "isvcd_function_selector.h"
/**
 *******************************************************************************
 *
 * @brief Initialize the intra/inter/transform/deblk function pointers of
 * codec context
 *
 * @par Description: the current routine initializes the function pointers of
 * codec context basing on the architecture in use for svc
 *
 * @param[in] ps_codec
 *  Codec context pointer
 *
 * @returns  none
 *
 * @remarks none
 *
 *******************************************************************************
 */
void isvcd_init_function_ptr_generic(svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    residual_sampling_ctxt_t *ps_resd_samp_ctx;
    intra_sampling_ctxt_t *ps_intra_samp_ctxt;
    dec_struct_t *ps_codec = &ps_svc_lyr_dec->s_dec;

    /* call the ih264 init ptr generic fn*/
    ih264d_init_function_ptr_generic(ps_codec);

    ps_resd_samp_ctx = (residual_sampling_ctxt_t *) ps_svc_lyr_dec->pv_residual_sample_ctxt;
    ps_intra_samp_ctxt = (intra_sampling_ctxt_t *) ps_svc_lyr_dec->pv_intra_sample_ctxt;

    ps_svc_lyr_dec->pf_pred_residual_recon_luma_4x4 = isvcd_pred_residual_recon_4x4;
    ps_svc_lyr_dec->pf_pred_residual_recon_luma_8x8 = isvcd_pred_residual_recon_8x8;
    ps_svc_lyr_dec->pf_pred_residual_recon_luma_16x16 = isvcd_pred_residual_recon_16x16;
    ps_svc_lyr_dec->pf_pred_residual_recon_chroma_4x4 = isvcd_pred_residual_recon_chroma_4x4;
    ps_svc_lyr_dec->pf_pred_residual_recon_chroma_8x8 = isvcd_pred_residual_recon_chroma_8x8;

    ps_svc_lyr_dec->pf_residual_luma_4x4 = isvcd_residual_luma_4x4;
    ps_svc_lyr_dec->pf_residual_luma_8x8 = isvcd_residual_luma_8x8;
    ps_svc_lyr_dec->pf_residual_luma_16x16 = isvcd_residual_luma_16x16;
    ps_svc_lyr_dec->pf_residual_chroma_cb_cr_8x8 = isvcd_residual_chroma_cb_cr_8x8;

    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_4x4 =
        isvcd_iquant_itrans_residual_recon_4x4;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_4x4_dc =
        isvcd_iquant_itrans_residual_recon_4x4_dc;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_8x8 =
        isvcd_iquant_itrans_residual_recon_8x8;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_luma_8x8_dc =
        isvcd_iquant_itrans_residual_recon_8x8_dc;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_chroma_4x4 =
        isvcd_iquant_itrans_residual_recon_chroma_4x4;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_recon_chroma_4x4_dc =
        isvcd_iquant_itrans_residual_recon_chroma_4x4_dc;

    ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_4x4 = isvcd_iquant_itrans_residual_4x4;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_4x4_dc = isvcd_iquant_itrans_residual_4x4_dc;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_8x8 = isvcd_iquant_itrans_residual_8x8;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_luma_8x8_dc = isvcd_iquant_itrans_residual_8x8_dc;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_chroma_4x4 = isvcd_iquant_itrans_residual_chroma_4x4;
    ps_svc_lyr_dec->pf_iquant_itrans_residual_chroma_4x4_dc =
        isvcd_iquant_itrans_residual_chroma_4x4_dc;

    ps_svc_lyr_dec->pf_iquant_itrans_luma_4x4 = isvcd_iquant_itrans_4x4;
    ps_svc_lyr_dec->pf_iquant_itrans_luma_4x4_dc = isvcd_iquant_itrans_4x4_dc;
    ps_svc_lyr_dec->pf_iquant_itrans_luma_8x8 = isvcd_iquant_itrans_8x8;
    ps_svc_lyr_dec->pf_iquant_itrans_luma_8x8_dc = isvcd_iquant_itrans_8x8_dc;
    ps_svc_lyr_dec->pf_iquant_itrans_chroma_4x4 = isvcd_iquant_itrans_chroma_4x4;
    ps_svc_lyr_dec->pf_iquant_itrans_chroma_4x4_dc = isvcd_iquant_itrans_chroma_4x4_dc;

    ps_intra_samp_ctxt->pf_interpolate_base_luma_dyadic = isvcd_interpolate_base_luma_dyadic;
    ps_intra_samp_ctxt->pf_interpolate_intra_base = isvcd_interpolate_intra_base;
    ps_intra_samp_ctxt->pf_vert_chroma_interpol[0] = isvcd_vert_interpol_chroma_dyadic_1;
    ps_intra_samp_ctxt->pf_vert_chroma_interpol[1] = isvcd_vert_interpol_chroma_dyadic_2;
    ps_intra_samp_ctxt->pf_vert_chroma_interpol[2] = isvcd_vert_interpol_chroma_dyadic_3;

    ps_intra_samp_ctxt->pf_horz_chroma_interpol[0] = isvcd_horz_interpol_chroma_dyadic_1;
    ps_intra_samp_ctxt->pf_horz_chroma_interpol[1] = isvcd_horz_interpol_chroma_dyadic_2;

    /*Dyadic Residual Resampling*/
    ps_resd_samp_ctx->pf_residual_luma_dyadic = isvcd_residual_luma_dyadic;
    ps_resd_samp_ctx->pf_residual_chroma_dyadic = isvcd_residual_chroma_dyadic;
    ps_resd_samp_ctx->pf_residual_chroma_dyadic_alt = isvcd_residual_chroma_dyadic_alt;

    /*Non-dyadic Residual Resampling*/
    ps_resd_samp_ctx->pf_interpolate_residual = isvcd_interpolate_residual;

    ps_resd_samp_ctx->pf_residual_reflayer_const_non_boundary_mb =
        isvcd_residual_reflayer_const_non_boundary_mb;
    ps_resd_samp_ctx->pf_residual_reflayer_const_boundary_mb =
        isvcd_residual_reflayer_const_boundary_mb;

    return;
}
