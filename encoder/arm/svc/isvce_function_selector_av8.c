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
*  isvce_function_selector_av8.c
*
* @brief
*  Contains functions to initialize function pointers of codec context
*
* @author
*  Ittiam
*
* @par List of Functions:
*  - isvce_init_function_ptr_generic
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
#include "ih264_inter_pred_filters.h"
#include "ih264_mem_fns.h"
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
#include "ih264e_platform_macros.h"
#include "isvce_cabac.h"
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
void isvce_init_function_ptr_neon_av8(isvce_codec_t *ps_codec)
{
    WORD32 i = 0;

    /* curr proc ctxt */
    isvce_process_ctxt_t *ps_proc = NULL;
    isvce_me_ctxt_t *ps_me_ctxt = NULL;
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;
    inter_pred_fxns_t *ps_inter_pred_fxns = &ps_isa_dependent_fxns->s_inter_pred_fxns;
    mem_fxns_t *ps_mem_fxns = &ps_isa_dependent_fxns->s_mem_fxns;

    /* Init function pointers for intra pred leaf level functions luma
     * Intra 16x16 */
    ps_codec->apf_intra_pred_16_l[0] = ih264_intra_pred_luma_16x16_mode_vert_av8;
    ps_codec->apf_intra_pred_16_l[1] = ih264_intra_pred_luma_16x16_mode_horz_av8;
    ps_codec->apf_intra_pred_16_l[2] = ih264_intra_pred_luma_16x16_mode_dc_av8;
    ps_codec->apf_intra_pred_16_l[3] = ih264_intra_pred_luma_16x16_mode_plane_av8;

    /* Init function pointers for intra pred leaf level functions luma
     * Intra 4x4 */
    ps_codec->apf_intra_pred_4_l[0] = ih264_intra_pred_luma_4x4_mode_vert_av8;
    ps_codec->apf_intra_pred_4_l[1] = ih264_intra_pred_luma_4x4_mode_horz_av8;
    ps_codec->apf_intra_pred_4_l[2] = ih264_intra_pred_luma_4x4_mode_dc_av8;
    ps_codec->apf_intra_pred_4_l[3] = ih264_intra_pred_luma_4x4_mode_diag_dl_av8;
    ps_codec->apf_intra_pred_4_l[4] = ih264_intra_pred_luma_4x4_mode_diag_dr_av8;
    ps_codec->apf_intra_pred_4_l[5] = ih264_intra_pred_luma_4x4_mode_vert_r_av8;
    ps_codec->apf_intra_pred_4_l[6] = ih264_intra_pred_luma_4x4_mode_horz_d_av8;
    ps_codec->apf_intra_pred_4_l[7] = ih264_intra_pred_luma_4x4_mode_vert_l_av8;
    ps_codec->apf_intra_pred_4_l[8] = ih264_intra_pred_luma_4x4_mode_horz_u_av8;

    /* Init function pointers for intra pred leaf level functions luma
     * Intra 8x8 */
    ps_codec->apf_intra_pred_8_l[0] = ih264_intra_pred_luma_8x8_mode_vert_av8;
    ps_codec->apf_intra_pred_8_l[2] = ih264_intra_pred_luma_8x8_mode_dc_av8;
    ps_codec->apf_intra_pred_8_l[3] = ih264_intra_pred_luma_8x8_mode_diag_dl_av8;
    ps_codec->apf_intra_pred_8_l[4] = ih264_intra_pred_luma_8x8_mode_diag_dr_av8;
    ps_codec->apf_intra_pred_8_l[5] = ih264_intra_pred_luma_8x8_mode_vert_r_av8;
    ps_codec->apf_intra_pred_8_l[6] = ih264_intra_pred_luma_8x8_mode_horz_d_av8;
    ps_codec->apf_intra_pred_8_l[7] = ih264_intra_pred_luma_8x8_mode_vert_l_av8;
    ps_codec->apf_intra_pred_8_l[8] = ih264_intra_pred_luma_8x8_mode_horz_u_av8;

    /* Init function pointers for intra pred leaf level functions chroma
     * Intra 8x8 */
    ps_codec->apf_intra_pred_c[0] = ih264_intra_pred_chroma_8x8_mode_dc_av8;
    ps_codec->apf_intra_pred_c[1] = ih264_intra_pred_chroma_8x8_mode_horz_av8;
    ps_codec->apf_intra_pred_c[2] = ih264_intra_pred_chroma_8x8_mode_vert_av8;
    ps_codec->apf_intra_pred_c[3] = ih264_intra_pred_chroma_8x8_mode_plane_av8;

    /* Init forward transform fn ptr */
    ps_enc_loop_fxns->apf_resi_trans_quant_8x8[0] = isvc_resi_trans_quant_8x8;
    ps_enc_loop_fxns->apf_resi_trans_quant_8x8[1] = isvc_resi_trans_quant_8x8;
    ps_enc_loop_fxns->apf_resi_trans_quant_4x4[0] = isvc_resi_trans_quant_4x4_neon;
    ps_enc_loop_fxns->apf_resi_trans_quant_4x4[1] =
        isvc_resi_trans_quant_4x4_with_residual_sub_neon;
    ps_enc_loop_fxns->apf_resi_trans_quant_chroma_4x4[0] = isvc_resi_trans_quant_chroma_4x4_neon;
    ps_enc_loop_fxns->apf_resi_trans_quant_chroma_4x4[1] =
        isvc_resi_trans_quant_chroma_4x4_with_residual_sub_neon;

    /* Init inverse transform fn ptr */
    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4[0] =
        isvc_iquant_itrans_recon_4x4_with_res_output_neon;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4[1] =
        isvc_iquant_itrans_recon_4x4_with_res_accumulate_neon;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4[2] = isvc_iquant_itrans_recon_4x4_neon;
    ;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_4x4_dc[2] = isvc_iquant_itrans_recon_4x4_dc_neon;

    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4[0] =
        isvc_iquant_itrans_recon_chroma_4x4_with_res_output_neon;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4[1] =
        isvc_iquant_itrans_recon_chroma_4x4_with_res_accumulate_neon;

    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4_dc[0] =
        isvc_iquant_itrans_recon_chroma_4x4_dc_with_res_output_neon;
    ps_enc_loop_fxns->apf_iquant_itrans_recon_chroma_4x4_dc[1] =
        isvc_iquant_itrans_recon_chroma_4x4_dc_with_res_accumulate_neon;

    ps_enc_loop_fxns->pf_ihadamard_scaling_4x4 = ih264_ihadamard_scaling_4x4_av8;
    ps_enc_loop_fxns->pf_ihadamard_scaling_2x2_uv = ih264_ihadamard_scaling_2x2_uv_av8;

    /* Init fn ptr luma core coding */
    ps_enc_loop_fxns->apf_luma_energy_compaction[0] = isvce_code_luma_intra_macroblock_16x16;
    ps_enc_loop_fxns->apf_luma_energy_compaction[1] = isvce_code_luma_intra_macroblock_4x4;
    ps_enc_loop_fxns->apf_luma_energy_compaction[3] = isvce_code_luma_inter_macroblock_16x16;

    /* Init fn ptr chroma core coding */
    ps_enc_loop_fxns->apf_chroma_energy_compaction[0] = isvce_code_chroma_intra_macroblock_8x8;
    ps_enc_loop_fxns->apf_chroma_energy_compaction[1] = isvce_code_chroma_inter_macroblock_8x8;

    /* Init fn ptr luma deblocking */
    ps_codec->pf_deblk_luma_vert_bs4 = ih264_deblk_luma_vert_bs4_av8;
    ps_codec->pf_deblk_luma_vert_bslt4 = ih264_deblk_luma_vert_bslt4_av8;
    ps_codec->pf_deblk_luma_horz_bs4 = ih264_deblk_luma_horz_bs4_av8;
    ps_codec->pf_deblk_luma_horz_bslt4 = ih264_deblk_luma_horz_bslt4_av8;

    /* Init fn ptr chroma deblocking */
    ps_codec->pf_deblk_chroma_vert_bs4 = ih264_deblk_chroma_vert_bs4_av8;
    ps_codec->pf_deblk_chroma_vert_bslt4 = ih264_deblk_chroma_vert_bslt4_av8;
    ps_codec->pf_deblk_chroma_horz_bs4 = ih264_deblk_chroma_horz_bs4_av8;
    ps_codec->pf_deblk_chroma_horz_bslt4 = ih264_deblk_chroma_horz_bslt4_av8;

    /* write mb syntax layer */
    ps_codec->pf_write_mb_syntax_layer[CAVLC][ISLICE] = isvce_write_islice_mb_cavlc;
    ps_codec->pf_write_mb_syntax_layer[CAVLC][PSLICE] = isvce_write_pslice_mb_cavlc;
    ps_codec->pf_write_mb_syntax_layer[CAVLC][BSLICE] = isvce_write_bslice_mb_cavlc;
    ps_codec->pf_write_mb_syntax_layer[CABAC][ISLICE] = isvce_write_islice_mb_cabac;
    ps_codec->pf_write_mb_syntax_layer[CABAC][PSLICE] = isvce_write_pslice_mb_cabac;

    /* Padding Functions */
    ps_codec->pf_pad_top = ih264_pad_top_av8;
    ps_codec->pf_pad_bottom = ih264_pad_bottom;
    ps_codec->pf_pad_left_luma = ih264_pad_left_luma_av8;
    ps_codec->pf_pad_left_chroma = ih264_pad_left_chroma_av8;
    ps_codec->pf_pad_right_luma = ih264_pad_right_luma_av8;
    ps_codec->pf_pad_right_chroma = ih264_pad_right_chroma_av8;

    /* Inter pred leaf level functions */
    ps_inter_pred_fxns->pf_inter_pred_luma_copy = ih264_inter_pred_luma_copy_av8;
    ps_inter_pred_fxns->pf_inter_pred_luma_horz = ih264_inter_pred_luma_horz_av8;
    ps_inter_pred_fxns->pf_inter_pred_luma_vert = ih264_inter_pred_luma_vert_av8;
    ps_inter_pred_fxns->pf_inter_pred_luma_bilinear = ih264_inter_pred_luma_bilinear;
    ps_inter_pred_fxns->pf_inter_pred_chroma = ih264_inter_pred_chroma_av8;

    /* sad me level functions */
    ps_codec->apf_compute_sad_16x16[0] = ime_compute_sad_16x16_av8;
    ps_codec->apf_compute_sad_16x16[1] = ime_compute_sad_16x16_fast_av8;
    ps_codec->pf_compute_sad_16x8 = ime_compute_sad_16x8_av8;

    /* memor handling operations */
    ps_mem_fxns->pf_mem_cpy_mul8 = ih264_memcpy_mul_8_av8;
    ps_mem_fxns->pf_mem_set_mul8 = ih264_memset_mul_8_av8;

    /* sad me level functions */
    for(i = 0; i < (MAX_PROCESS_CTXT); i++)
    {
        ps_proc = &ps_codec->as_process[i];
        ps_me_ctxt = &ps_proc->s_me_ctxt;
        ps_me_ctxt->pf_ime_compute_sad_16x16[0] = ime_compute_sad_16x16_av8;
        ps_me_ctxt->pf_ime_compute_sad_16x16[1] = ime_compute_sad_16x16_fast_av8;
        ps_me_ctxt->pf_ime_compute_sad_16x8 = ime_compute_sad_16x8_av8;
        ps_me_ctxt->pf_ime_compute_sad4_diamond = ime_calculate_sad4_prog_av8;
        ps_me_ctxt->pf_ime_compute_sad3_diamond = ime_calculate_sad3_prog_av8;
        ps_me_ctxt->pf_ime_compute_sad2_diamond = ime_calculate_sad2_prog_av8;
        ps_me_ctxt->pf_ime_sub_pel_compute_sad_16x16 = ime_sub_pel_compute_sad_16x16_av8;
        ps_me_ctxt->pf_ime_compute_sad_stat_luma_16x16 = ime_compute_satqd_16x16_lumainter_av8;
    }

    /* intra mode eval -encoder level function */
    ps_codec->pf_ih264e_evaluate_intra16x16_modes = ih264e_evaluate_intra16x16_modes_av8;
    ps_codec->pf_ih264e_evaluate_intra_chroma_modes = ih264e_evaluate_intra_chroma_modes_av8;
    ps_codec->pf_ih264e_evaluate_intra_4x4_modes = ih264e_evaluate_intra_4x4_modes;

    /* csc */
    ps_codec->pf_ih264e_conv_420p_to_420sp = ih264e_fmt_conv_420p_to_420sp;
    ps_codec->pf_ih264e_fmt_conv_422i_to_420sp = ih264e_fmt_conv_422i_to_420sp;

    /* Halp pel generation function - encoder level*/
    ps_codec->pf_ih264e_sixtapfilter_horz = ih264e_sixtapfilter_horz_av8;
    ps_codec->pf_ih264e_sixtap_filter_2dvh_vert = ih264e_sixtap_filter_2dvh_vert_av8;
}
