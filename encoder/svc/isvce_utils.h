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
*  isvce_utils.h
*
* @brief
*  Contains function declarations for function declared in ih264e_svc_utils.c
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_UTILS_H_
#define _ISVCE_UTILS_H_

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "ih264_typedefs.h"
#include "ih264e_bitstream.h"
/* Dependencies of 'irc_picture_type.h' */
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "irc_mem_req_and_acq.h"
/* Dependencies of 'irc_rate_control_api_structs' */
#include "irc_picture_type.h"
#include "irc_rd_model.h"
#include "irc_vbr_storage_vbv.h"
#include "irc_est_sad.h"
#include "irc_bit_allocation.h"
#include "irc_mb_model_based.h"
#include "irc_cbr_buffer_control.h"
#include "irc_vbr_str_prms.h"
#include "irc_common.h"
#include "irc_rate_control_api_structs.h"
#include "isvc_macros.h"
#include "isvc_structs.h"
#include "isvce_defs.h"
#include "isvce_globals.h"
#include "isvce_interface_structs.h"
#include "isvce_structs.h"

static FORCEINLINE void isvce_svc_au_buf_init(svc_au_buf_t *ps_svc_pic_buf,
                                              svc_params_t *ps_svc_params)
{
    ps_svc_pic_buf->i1_temporal_id = -1;
    ps_svc_pic_buf->u1_num_spatial_layers = ps_svc_params->u1_num_spatial_layers;
    ps_svc_pic_buf->d_spatial_res_ratio = ps_svc_params->d_spatial_res_ratio;
}

static FORCEINLINE WORD8 isvce_svc_temporal_id_compute(WORD32 i4_poc, UWORD8 u1_num_temporal_layers,
                                                       PIC_TYPE_T e_pic_type)
{
    if(e_pic_type == PIC_IDR)
    {
        return 0;
    }
    else
    {
        return i4_poc % u1_num_temporal_layers;
    }
}

static FORCEINLINE WORD32 isvcd_get_ceil_log2(WORD32 i4_input)
{
    WORD32 i4_bits = 0;

    /* check for negative number */
    ASSERT(i4_input >= 0);

    i4_input--;

    while(i4_input > 0)
    {
        i4_bits++;
        i4_input >>= 1;
    }

    return (i4_bits);
}
/**
*******************************************************************************
*
* @brief calculate coded subblock pattern from nnz
*
* @par Description:
*  calculate coded subblock pattern from nnz
*
* @param[in] ps_proc
*  process context
*
* @returns  csbp
*
* @remarks  none
*
*******************************************************************************
*/
static FORCEINLINE UWORD32 isvce_calculate_csbp(isvce_process_ctxt_t *ps_proc)
{
    WORD32 i;

    UWORD8 *pu1_curr_nnz = ((UWORD8 *) ps_proc->au4_nnz) + 1;
    UWORD32 u4_csbp = 0;

    for(i = 0; i < 16; i++)
    {
        UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[i];

        u4_csbp |= ((!!pu1_curr_nnz[i]) << u1_zscan_idx);
    }

    return u4_csbp;
}

static FORCEINLINE UWORD8 isvce_check_identical_mv(isvce_enc_pu_mv_t *ps_mv1,
                                                   isvce_enc_pu_mv_t *ps_mv2,
                                                   PRED_MODE_T e_pred_mode)
{
    if(e_pred_mode != L0)
    {
        if(!((ps_mv1[L1].i1_ref_idx == ps_mv2[L1].i1_ref_idx) &&
             (ps_mv1[L1].s_mv.i2_mvx == ps_mv2[L1].s_mv.i2_mvx) &&
             (ps_mv1[L1].s_mv.i2_mvy == ps_mv2[L1].s_mv.i2_mvy)))
        {
            return 0;
        }
    }

    if(e_pred_mode != L1)
    {
        if(!((ps_mv1[L0].i1_ref_idx == ps_mv2[L0].i1_ref_idx) &&
             (ps_mv1[L0].s_mv.i2_mvx == ps_mv2[L0].s_mv.i2_mvx) &&
             (ps_mv1[L0].s_mv.i2_mvy == ps_mv2[L0].s_mv.i2_mvy)))
        {
            return 0;
        }
    }

    return 1;
}

static FORCEINLINE WORD32 isvce_get_num_bits(bitstrm_t *ps_bitstream)
{
    return GET_NUM_BITS(ps_bitstream);
}

extern WORD32 ih264e_get_min_level(WORD32 wd, WORD32 ht);

extern WORD32 isvce_svc_au_props_validate(svc_inp_params_t *ps_svc_inp_params, UWORD32 u4_inp_wd,
                                          UWORD32 u4_inp_ht, UWORD32 u4_svc_comp_wd,
                                          UWORD32 u4_svc_comp_ht);

extern WORD32 isvce_svc_inp_params_validate(isvce_init_ip_t *ps_ip, isvce_cfg_params_t *ps_cfg);

extern WORD32 isvce_svc_rc_params_validate(isvce_cfg_params_t *ps_cfg);

extern WORD32 isvce_svc_frame_params_validate(
    rate_control_api_t *aps_rate_control_api[MAX_NUM_SPATIAL_LAYERS], UWORD8 u1_num_spatial_layers);

extern WORD32 isvce_get_total_svc_au_buf_size(svc_inp_params_t *ps_svc_inp_params,
                                              WORD32 i4_pic_size, WORD32 i4_level,
                                              WORD32 i4_horz_pad, WORD32 i4_vert_pad,
                                              WORD32 i4_num_ref_frames,
                                              WORD32 i4_num_reorder_frames);

extern UWORD32 isvce_get_total_svc_au_data_size(WORD32 i4_num_luma_samples,
                                                UWORD8 u1_num_spatial_layers,
                                                DOUBLE d_spatial_res_ratio);

extern IH264E_ERROR_T isvce_svc_au_data_mgr_add_bufs(isvce_codec_t *ps_codec);

extern IH264E_ERROR_T isvce_svc_au_buf_mgr_add_bufs(isvce_codec_t *ps_codec);

extern UWORD32 isvce_get_svc_inp_buf_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                          UWORD32 u4_wd, UWORD32 u4_ht);

extern void isvce_svc_inp_buf_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec);

extern void isvce_init_svc_dimension(isvce_inp_buf_t *ps_inp);

extern void isvce_svc_inp_buf_populate(isvce_codec_t *ps_codec, isvce_inp_buf_t *ps_inp);

extern void isvce_get_svc_compliant_dimensions(UWORD8 u1_num_spatial_layers,
                                               DOUBLE d_scaling_factor, UWORD32 u4_wd,
                                               UWORD32 u4_ht, UWORD32 *pu4_svc_comp_wd,
                                               UWORD32 *pu4_svc_comp_ht);

extern UWORD32 isvce_get_svc_nbr_info_buf_size(UWORD8 u1_num_spatial_layers,
                                               DOUBLE d_spatial_res_ratio, UWORD32 u4_wd,
                                               UWORD32 u4_ht);

extern void isvce_svc_nbr_info_buf_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec);

extern IH264E_ERROR_T isvce_svc_au_init(isvce_codec_t *ps_codec, isvce_inp_buf_t *ps_inp_buf);

extern IH264E_ERROR_T isvce_svc_layer_pic_init(isvce_codec_t *ps_codec, isvce_inp_buf_t *ps_inp_buf,
                                               UWORD8 u1_spatial_layer_id);

extern IH264E_ERROR_T isvce_init_layer_proc_ctxt(isvce_process_ctxt_t *ps_proc);

extern UWORD32 isvce_get_svc_ilp_buf_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                          UWORD32 u4_wd, UWORD32 u4_ht);

extern void isvce_svc_ilp_buf_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec);

extern void isvce_svc_ilp_buf_update(isvce_process_ctxt_t *ps_proc);

extern void isvce_svc_pad_frame(isvce_process_ctxt_t *ps_proc);

extern IH264E_ERROR_T isvce_init_air_map(isvce_codec_t *ps_codec);

extern void isvce_derive_nghbr_avbl_of_mbs(isvce_process_ctxt_t *ps_proc);

extern void isvce_init_quant_params(isvce_process_ctxt_t *ps_proc, WORD32 qp);

extern IH264E_ERROR_T isvce_codec_init(isvce_codec_t *ps_codec);

extern IH264E_ERROR_T isvce_codec_update_config(isvce_codec_t *ps_codec,
                                                isvce_cfg_params_t *ps_cfg);

extern WORD32 isvce_input_queue_update(isvce_codec_t *ps_codec, ive_video_encode_ip_t *ps_ive_ip,
                                       isvce_inp_buf_t *ps_enc_buff, WORD8 i1_layer_id);

extern void isvce_join_threads(isvce_codec_t *ps_codec);

extern UWORD32 isvce_get_min_outbuf_size(UWORD32 u4_wd, UWORD32 u4_ht,
                                         UWORD8 u1_num_spatial_layers);

#endif
