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
* isvce_cabac_init.c
*
* @brief
*  Contains all initialization functions for cabac contexts
*
* @author
*  Doney Alex
*
* @par List of Functions:
*
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
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

/* User include files */
#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "isvc_defs.h"
#include "ih264_debug.h"
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
#include "ih264_platform_macros.h"
#include "isvc_macros.h"
#include "ih264_buf_mgr.h"
#include "ih264e_error.h"
#include "ih264e_bitstream.h"
#include "isvc_common_tables.h"
#include "isvc_cabac_tables.h"
#include "ih264_list.h"
#include "isvce_defs.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "isvce_rate_control.h"
#include "isvce_cabac_structs.h"
#include "isvce_structs.h"
#include "isvce_cabac.h"
#include "isvce_process.h"
#include "ithread.h"
#include "isvce_encode_header.h"
#include "isvce_globals.h"
#include "ih264e_config.h"
#include "ih264e_trace.h"
#include "ih264e_statistics.h"
#include "ih264_cavlc_tables.h"
#include "isvce_deblk.h"
#include "isvce_me.h"
#include "ih264e_debug.h"
#include "ih264e_master.h"
#include "isvce_utils.h"
#include "irc_mem_req_and_acq.h"
#include "irc_rate_control_api.h"
#include "ih264e_platform_macros.h"
#include "ime_statistics.h"

/*****************************************************************************/
/*  Function definitions .                                                   */
/*****************************************************************************/

/**
 *******************************************************************************
 *
 * @brief
 * Initialize cabac encoding environment
 *
 * @param[in] ps_cab_enc_env
 *  Pointer to encoding_envirnoment_t structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_init_cabac_enc_envirnoment(encoding_envirnoment_t *ps_cab_enc_env)
{
    ps_cab_enc_env->u4_code_int_low = 0;
    ps_cab_enc_env->u4_code_int_range = 0x1fe;
    ps_cab_enc_env->u4_out_standing_bytes = 0;
    ps_cab_enc_env->u4_bits_gen = 0;
}

/**
 *******************************************************************************
 *
 * @brief
 * Initialize default context values and pointers (Called once at the beginning
 *of encoding).
 *
 * @param[in] ps_ent_ctxt
 *  Pointer to entropy context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
void isvce_init_cabac_table(isvce_entropy_ctxt_t *ps_ent_ctxt)
{
    /* CABAC context */
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;
    ps_cabac_ctxt->ps_mb_map_ctxt_inc = ps_cabac_ctxt->ps_mb_map_ctxt_inc_base + 1;
    ps_cabac_ctxt->ps_lft_csbp = &ps_cabac_ctxt->s_lft_csbp;
    ps_cabac_ctxt->ps_bitstrm = ps_ent_ctxt->ps_bitstrm;

    {
        /* 0th entry of mb_map_ctxt_inc will be always be containing default values
         */
        /* for CABAC context representing MB not available                       */
        isvce_mb_info_ctxt_t *ps_def_ctxt = ps_cabac_ctxt->ps_mb_map_ctxt_inc - 1;

        ps_def_ctxt->u1_mb_type = CAB_SKIP;
        ps_def_ctxt->u1_cbp = 0x0f;
        ps_def_ctxt->u1_intrapred_chroma_mode = 0;
        ps_def_ctxt->u1_base_mode_flag = 0;

        memset(ps_def_ctxt->i1_ref_idx, 0, sizeof(ps_def_ctxt->i1_ref_idx));
        memset(ps_def_ctxt->u1_mv, 0, sizeof(ps_def_ctxt->u1_mv));
        ps_cabac_ctxt->ps_def_ctxt_mb_info = ps_def_ctxt;
    }
}

/**
 *******************************************************************************
 *
 * @brief
 * Initialize cabac context: Initialize all contest with init values given in
 *the spec. Called at the beginning of entropy coding of each slice for CABAC
 *encoding.
 *
 * @param[in] ps_ent_ctxt
 *  Pointer to entropy context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
void isvce_init_cabac_ctxt(isvce_entropy_ctxt_t *ps_ent_ctxt, slice_header_t *ps_slice_hdr)
{
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;

    const UWORD8 u1_slice_type = ps_slice_hdr->u1_slice_type;
    WORD8 i1_cabac_init_idc = 0;
    bin_ctxt_model *au1_cabac_ctxt_table = ps_cabac_ctxt->au1_cabac_ctxt_table;
    UWORD8 u1_qp_y = ps_slice_hdr->i1_slice_qp;

    isvce_init_cabac_enc_envirnoment(&ps_cabac_ctxt->s_cab_enc_env);

    ps_cabac_ctxt->i1_prevps_mb_qp_delta_ctxt = 0;

    if(ISLICE != u1_slice_type)
    {
        i1_cabac_init_idc = ps_slice_hdr->i1_cabac_init_idc;
    }
    else
    {
        i1_cabac_init_idc = 3;
    }

    memcpy(au1_cabac_ctxt_table, gau1_isvc_cabac_ctxt_init_table[i1_cabac_init_idc][u1_qp_y],
           NUM_SVC_CABAC_CTXTS * sizeof(bin_ctxt_model));
}
