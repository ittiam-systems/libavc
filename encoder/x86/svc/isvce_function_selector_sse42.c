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
*  isvce_function_selector_sse42.c
*
* @brief
*  Contains functions to initialize function pointers of codec context
*
* @author
*  Ittiam
*
* @par List of Functions:
*  - isvce_init_function_ptr_sse42
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
#include "iv2.h"
#include "ive2.h"
#include "isvc_defs.h"
#include "ih264_size_defs.h"
#include "isvce_defs.h"
#include "ih264e_error.h"
#include "ih264e_bitstream.h"
#include "ime_distortion_metrics.h"
#include "ime_defs.h"
#include "ime_structs.h"
#include "ih264_error.h"
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
#include "isvce_cabac.h"
#include "ih264e_platform_macros.h"
#include "isvce_core_coding.h"
#include "ih264_cavlc_tables.h"
#include "isvce_cavlc.h"
#include "ih264e_intra_modes_eval.h"
#include "ih264e_fmt_conv.h"
#include "ih264e_half_pel.h"

/**
*******************************************************************************
*
* @brief Initialize the intra/inter/transform/deblk function pointers of
* codec context
*
* @par Description: the current routine initializes the function pointers of
* codec context basing on the architecture in use
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
void isvce_init_function_ptr_sse42(isvce_codec_t *ps_codec)
{
    WORD32 i;
    isvce_process_ctxt_t *ps_proc = NULL;
    isvce_me_ctxt_t *ps_me_ctxt = NULL;
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;
    mem_fxns_t *ps_mem_fxns = &ps_isa_dependent_fxns->s_mem_fxns;

    ps_enc_loop_fxns->pf_hadamard_quant_4x4 = isvc_hadamard_quant_4x4_sse42;
    ps_enc_loop_fxns->pf_hadamard_quant_2x2_uv = isvc_hadamard_quant_2x2_uv_sse42;

    ps_enc_loop_fxns->apf_resi_trans_quant_4x4[0] = isvc_resi_trans_quant_4x4_sse42;
    ps_enc_loop_fxns->apf_resi_trans_quant_4x4[1] = isvc_resi_trans_quant_4x4_with_res_pred_sse42;

    ps_enc_loop_fxns->apf_resi_trans_quant_chroma_4x4[0] = isvc_resi_trans_quant_chroma_4x4_sse42;
    ps_enc_loop_fxns->apf_resi_trans_quant_chroma_4x4[1] =
        isvc_resi_trans_quant_chroma_4x4_with_res_pred_sse42;

    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4[0] = isvc_iquant_itrans_recon_res_4x4_sse42;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4[1] =
        isvc_iquant_itrans_recon_res_4x4_with_res_acc_sse42;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4[2] = isvc_iquant_itrans_recon_4x4_sse42;

    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4[0] =
        isvc_iquant_itrans_recon_res_chroma_4x4_sse42;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4[1] =
        isvc_iquant_itrans_recon_res_chroma_4x4_with_res_acc_sse42;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4[2] =
        isvc_iquant_itrans_recon_chroma_4x4_sse42;

    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4_dc[0] = isvc_iquant_itrans_recon_res_dc_4x4_sse42;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4_dc[1] =
        isvc_iquant_itrans_recon_res_dc_with_res_acc_4x4_sse42;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4_dc[2] = isvc_iquant_itrans_recon_dc_4x4_sse42;

    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4_dc[0] =
        isvc_iquant_itrans_recon_res_chroma_4x4_dc_sse42;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4_dc[1] =
        isvc_iquant_itrans_recon_res_chroma_4x4_dc_with_res_acc_sse42;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4_dc[2] =
        isvc_iquant_itrans_recon_chroma_4x4_dc_sse42;

    ps_enc_loop_fxns->pf_ihadamard_scaling_4x4 = ih264_ihadamard_scaling_4x4_sse42;

    /* sad me level functions */
    ps_codec->apf_compute_sad_16x16[0] = ime_compute_sad_16x16_sse42;
    ps_codec->apf_compute_sad_16x16[1] = ime_compute_sad_16x16_fast_sse42;
    ps_codec->pf_compute_sad_16x8 = ime_compute_sad_16x8_sse42;

    ps_mem_fxns->pf_copy_2d = isvc_copy_2d_ssse3;
    ps_mem_fxns->pf_memset_2d = isvc_memset_2d_sse42;

    /* sad me level functions */
    for(i = 0; i < (MAX_PROCESS_CTXT); i++)
    {
        ps_proc = &ps_codec->as_process[i];

        ps_me_ctxt = &ps_proc->s_me_ctxt;
        ps_me_ctxt->pf_ime_compute_sad_16x16[0] = ime_compute_sad_16x16_sse42;
        ps_me_ctxt->pf_ime_compute_sad_16x16[1] = ime_compute_sad_16x16_fast_sse42;
        ps_me_ctxt->pf_ime_compute_sad_16x8 = ime_compute_sad_16x8_sse42;
        ps_me_ctxt->pf_ime_compute_sad4_diamond = ime_calculate_sad4_prog_sse42;
        ps_me_ctxt->pf_ime_sub_pel_compute_sad_16x16 = ime_sub_pel_compute_sad_16x16_sse42;
        ps_me_ctxt->pf_ime_compute_sad_stat_luma_16x16 = ime_compute_satqd_16x16_lumainter_sse42;
    }
}
