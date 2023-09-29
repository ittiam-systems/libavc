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
*  ih264e_svc_utils.c
*
* @brief
*  Contains utility functions used for SVC encoding
*
* @author
*  ittiam
*
* @par List of Functions:
*  - ih264e_svc_ref_list_refresh()
*  - ih264e_svc_inp_params_validate()
*
* @remarks
*  None
*
*******************************************************************************
*/
#include <math.h>
#include <limits.h>

#include "ih264_typedefs.h"

/* Dependencies of ih264_buf_mgr.h */
/* Dependencies of ih264_list.h */
#include "ih264_error.h"

#include "ih264_buf_mgr.h"
#include "ih264_list.h"
#include "ih264_trans_data.h"
#include "ih264_size_defs.h"

/* Dependencies of ih264_common_tables.h */
#include "ih264_defs.h"
#include "ih264_structs.h"

#include "ih264_common_tables.h"

/* Dependencies of ih264e_bitstream.h */
#include "ih264e_error.h"

/* Dependencies of ih264e_cabac_structs.h */
#include "ih264_cabac_tables.h"

/* Dependencies of ime_structs.h */
#include "ime_defs.h"
#include "ime_distortion_metrics.h"

/* Dependencies of ih264e_structs.h */
#include "iv2.h"
#include "ive2.h"
#include "ih264_defs.h"
#include "ih264_deblk_edge_filters.h"
#include "ih264_inter_pred_filters.h"
#include "ih264_structs.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include "ih264e_bitstream.h"
#include "ih264e_cabac_structs.h"
#include "ime_statistics.h"
#include "ime_structs.h"
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
#include "irc_rate_control_api.h"
#include "irc_svc_rate_control_api.h"
/* Dependencies of 'ih264e_utils.h' */
#include "ih264e_defs.h"
#include "ih264e_structs.h"
/* Dependencies of 'ih264e_utils.h' */
#include "irc_mem_req_and_acq.h"
#include "ih264e_rc_mem_interface.h"
#include "ih264e_time_stamp.h"
#include "ih264e_utils.h"
#include "ime.h"
#include "isvc_macros.h"
#include "isvce_cabac.h"
#include "isvce_core_coding.h"
#include "isvce_defs.h"
#include "isvce_error.h"
#include "isvce_me.h"
#include "isvce_utils.h"
#include "isvce_downscaler.h"
#include "isvce_encode_header.h"
#include "isvce_rate_control.h"
#include "isvce_sub_pic_rc.h"

static const UWORD32 gu4_downscaler_blk_size = 96;

static FORCEINLINE UWORD32 isvce_get_downscaler_blk_dims(UWORD32 u4_frame_dim, UWORD32 u4_blk_pos,
                                                         UWORD32 u4_default_blk_size)
{
    return ((u4_frame_dim - u4_blk_pos * u4_default_blk_size) < u4_default_blk_size)
               ? (u4_frame_dim - u4_blk_pos * u4_default_blk_size)
               : u4_default_blk_size;
}

/**
*******************************************************************************
*
* @brief
*  Reference and MV bank Buffer Manager for SVC
*
* @par Description:
*  Here we will
*      1) Find the correct ref pics for the current frame
*      2) Free the ref pics that are not going to be used anymore
*
*  1) Finding correct ref pic
*      All pics needed for future are arranged in a picture list called
*      ps_codec->as_ref_set. Each picture in this will have a pic buffer and
*      MV buffer that is marked appropriately as BUF_MGR_REF, BUF_MGR_IO or
*      BUF_MGR_CODEC. pic_cnt, poc, and temporal_id will also be present.
*      The strategy is to pick the closest references that belongs to the
*      same temporal_id or lesser. The closeness is measured via the
*      smallest absolute difference between ref and cur pocs.
*
*      Note that i4_pic_cnt == -1 is used to filter uninitialised ref pics.
*      Now since we only have max two ref pics, we will always find max 2
*      ref pics.
*
*  2) Self explanatory
*
* @param[in] ps_codec
*  Pointer to codeec context
*
* @param[in] pps_ref_pics
*  Array of pointers to refPicBufs
*
* @param[in] pps_mv_bufs
*  Array of pointers to refMVBufs
*
* @param[in] e_pic_type
*  Picture type
*
* @returns Nothing
*
*******************************************************************************
*/
static WORD32 isvce_ref_list_refresh(isvce_codec_t *ps_codec, svc_au_buf_t **pps_ref_pics,
                                     svc_au_data_t **pps_mv_bufs, WORD32 *pi4_ref_set_id,
                                     PIC_TYPE_T e_pic_type)
{
    typedef struct
    {
        WORD32 i4_buf_id;

        WORD32 i4_abs_poc_diff;

        WORD8 i1_temporal_id;
    } ref_pic_props_t;

    ref_pic_props_t s_ref_pic_props = {0, 0, -1};

    WORD32 i, buf_status;

    WORD32 i4_cur_pic_poc = ps_codec->i4_poc;
    WORD32 i4_cur_pic_temporal_id = isvce_svc_temporal_id_compute(
        ps_codec->i4_poc, ps_codec->s_cfg.s_svc_params.u1_num_temporal_layers, e_pic_type);

    if(e_pic_type == PIC_B)
    {
        return IH264E_FAIL;
    }

    ASSERT(1 == MAX_LAYER_REFERENCE_PICS);

    /* Pick a ref_pic for the current picture */
    if(e_pic_type != PIC_IDR)
    {
        for(i = 0; i < ps_codec->i4_ref_buf_cnt; i++)
        {
            WORD32 i4_abs_poc_diff;
            WORD8 i1_temporal_id;

            if(ps_codec->as_ref_set[i].i4_pic_cnt == -1)
            {
                continue;
            }

            buf_status = ih264_buf_mgr_get_status(ps_codec->pv_ref_buf_mgr,
                                                  ps_codec->as_ref_set[i].ps_pic_buf->i4_buf_id);

            if(buf_status & BUF_MGR_REF)
            {
                i4_abs_poc_diff = ABS(ps_codec->as_ref_set[i].i4_poc - i4_cur_pic_poc);
                i1_temporal_id = ps_codec->as_ref_set[i].ps_pic_buf->i1_temporal_id;

                if(s_ref_pic_props.i1_temporal_id > -1)
                {
                    if((i1_temporal_id <= i4_cur_pic_temporal_id) &&
                       (s_ref_pic_props.i4_abs_poc_diff > i4_abs_poc_diff))
                    {
                        s_ref_pic_props.i4_abs_poc_diff = i4_abs_poc_diff;
                        s_ref_pic_props.i1_temporal_id = i1_temporal_id;
                        s_ref_pic_props.i4_buf_id = i;
                    }
                }
                else if(i1_temporal_id <= i4_cur_pic_temporal_id)
                {
                    s_ref_pic_props.i4_abs_poc_diff = i4_abs_poc_diff;
                    s_ref_pic_props.i1_temporal_id = i1_temporal_id;
                    s_ref_pic_props.i4_buf_id = i;
                }
            }
        }

        if(s_ref_pic_props.i1_temporal_id < 0)
        {
            return IH264E_FAIL;
        }

        pps_ref_pics[0] = pps_ref_pics[1] =
            ps_codec->as_ref_set[s_ref_pic_props.i4_buf_id].ps_pic_buf;
        pps_mv_bufs[0] = pps_mv_bufs[1] =
            ps_codec->as_ref_set[s_ref_pic_props.i4_buf_id].ps_svc_au_data;

        /* Pick all ref pic_bufs to be freed. */
        for(i = 0; i < ps_codec->i4_ref_buf_cnt; i++)
        {
            if(ps_codec->as_ref_set[i].i4_pic_cnt == -1)
            {
                continue;
            }

            buf_status = ih264_buf_mgr_get_status(ps_codec->pv_ref_buf_mgr,
                                                  ps_codec->as_ref_set[i].ps_pic_buf->i4_buf_id);

            if((buf_status & (BUF_MGR_REF | BUF_MGR_CODEC | BUF_MGR_IO)) == 0)
            {
                ps_codec->as_ref_set[i].i4_pic_cnt = -1;
                ps_codec->as_ref_set[i].i4_poc = 32768;

                continue;
            }

            if(buf_status & BUF_MGR_REF)
            {
                if((i4_cur_pic_temporal_id <= ps_codec->as_ref_set[i].ps_pic_buf->i1_temporal_id) &&
                   (pps_ref_pics[0]->i4_frame_num !=
                    ps_codec->as_ref_set[i].ps_pic_buf->i4_frame_num) &&
                   (pps_ref_pics[0]->i4_frame_num !=
                    ps_codec->as_ref_set[i].ps_pic_buf->i4_frame_num))
                {
                    ih264_buf_mgr_release(ps_codec->pv_svc_au_data_store_mgr,
                                          ps_codec->as_ref_set[i].ps_pic_buf->i4_buf_id,
                                          BUF_MGR_REF);

                    ih264_buf_mgr_release(ps_codec->pv_ref_buf_mgr,
                                          ps_codec->as_ref_set[i].ps_pic_buf->i4_buf_id,
                                          BUF_MGR_REF);
                }
            }
        }
    }
    else
    {
        for(i = 0; i < ps_codec->i4_ref_buf_cnt; i++)
        {
            if(ps_codec->as_ref_set[i].i4_pic_cnt == -1)
            {
                continue;
            }

            buf_status = ih264_buf_mgr_get_status(ps_codec->pv_ref_buf_mgr,
                                                  ps_codec->as_ref_set[i].ps_pic_buf->i4_buf_id);

            if((buf_status & (BUF_MGR_REF | BUF_MGR_CODEC | BUF_MGR_IO)) == 0)
            {
                ps_codec->as_ref_set[i].i4_pic_cnt = -1;
                ps_codec->as_ref_set[i].i4_poc = 32768;

                continue;
            }

            if(buf_status & BUF_MGR_REF)
            {
                ih264_buf_mgr_release(ps_codec->pv_svc_au_data_store_mgr,
                                      ps_codec->as_ref_set[i].ps_pic_buf->i4_buf_id, BUF_MGR_REF);

                ih264_buf_mgr_release(ps_codec->pv_ref_buf_mgr,
                                      ps_codec->as_ref_set[i].ps_pic_buf->i4_buf_id, BUF_MGR_REF);
            }
        }
    }

    /*
     * Mark all reference pic with unused buffers to be free
     * We need this step since each one, ie ref, recon io etc only unset their
     * respective flags. Hence we need to combine togather and mark the ref set
     * accordingly
     */
    pi4_ref_set_id[0] = -1;

    for(i = 0; i < ps_codec->i4_ref_buf_cnt; i++)
    {
        if(ps_codec->as_ref_set[i].i4_pic_cnt == -1)
        {
            pi4_ref_set_id[0] = i;
            continue;
        }

        buf_status = ih264_buf_mgr_get_status(ps_codec->pv_ref_buf_mgr,
                                              ps_codec->as_ref_set[i].ps_pic_buf->i4_buf_id);

        if((buf_status & (BUF_MGR_REF | BUF_MGR_CODEC | BUF_MGR_IO)) == 0)
        {
            ps_codec->as_ref_set[i].i4_pic_cnt = -1;
            ps_codec->as_ref_set[i].i4_poc = 32768;

            pi4_ref_set_id[0] = i;
        }
    }

    /* An asssert failure here means we donot have any free buffs */
    if(pi4_ref_set_id[0] < 0)
    {
        return IH264E_FAIL;
    }

    return IH264E_SUCCESS;
}

/**
*******************************************************************************
*
* @brief
*  Validates SVC AU properties
*
* @param[in] ps_cfg
*  Cfg parameters
*
* @returns  error code in conformance with 'IH264E_ERROR_T'
*
*******************************************************************************
*/
WORD32 isvce_svc_au_props_validate(svc_inp_params_t *ps_svc_inp_params, UWORD32 u4_inp_wd,
                                   UWORD32 u4_inp_ht, UWORD32 u4_svc_comp_wd,
                                   UWORD32 u4_svc_comp_ht)
{
    typedef struct
    {
        DOUBLE d_spatial_res_ratio;

        UWORD8 u1_max_num_spatial_layers;
    } spatial_layer_props_t;

    UWORD8 i;
    UWORD32 au4_svc_wd[MAX_NUM_SPATIAL_LAYERS];
    UWORD32 au4_svc_ht[MAX_NUM_SPATIAL_LAYERS];

    DOUBLE d_scaling_factor = ps_svc_inp_params->d_spatial_res_ratio;
    UWORD8 u1_num_spatial_layers = ps_svc_inp_params->u1_num_spatial_layers;
    const spatial_layer_props_t gas_valid_spatial_layer_props[] = {{1.5, 2}, {2, 3}};
    UWORD32 u4_error_code = IV_SUCCESS;
    const UWORD8 u1_min_num_temporal_layers = 1;
    const UWORD8 u1_min_num_spatial_layers = 1;
    const UWORD8 u1_max_num_temporal_layers = MAX_NUM_TEMPORAL_LAYERS;
    const UWORD8 u1_max_num_spatial_layers = MAX_NUM_SPATIAL_LAYERS;
    const UWORD8 u1_num_valid_spatial_layer_props =
        sizeof(gas_valid_spatial_layer_props) / sizeof(gas_valid_spatial_layer_props[0]);

    if((ps_svc_inp_params->u1_num_temporal_layers < u1_min_num_temporal_layers) ||
       (ps_svc_inp_params->u1_num_temporal_layers > u1_max_num_temporal_layers))
    {
        u4_error_code |= IH264E_INVALID_SVC_PARAMS | IH264E_INVALID_NUM_TEMPORAL_LAYERS;
    }

    if((ps_svc_inp_params->u1_num_spatial_layers < u1_min_num_spatial_layers) ||
       (ps_svc_inp_params->u1_num_spatial_layers > u1_max_num_spatial_layers))
    {
        u4_error_code |= IH264E_INVALID_SVC_PARAMS | IH264E_INVALID_NUM_SPATIAL_LAYERS;
    }

    {
        UWORD8 u1_is_input_ratio_valid = 0;

        for(i = 0; i < u1_num_valid_spatial_layer_props; i++)
        {
            if(ps_svc_inp_params->d_spatial_res_ratio ==
               gas_valid_spatial_layer_props[i].d_spatial_res_ratio)
            {
                u1_is_input_ratio_valid = 1;

                if(ps_svc_inp_params->u1_num_spatial_layers >
                   gas_valid_spatial_layer_props[i].u1_max_num_spatial_layers)
                {
                    u4_error_code |= IH264E_INVALID_SVC_PARAMS | IH264E_INVALID_NUM_SPATIAL_LAYERS;
                }

                break;
            }
        }

        if(!u1_is_input_ratio_valid)
        {
            u4_error_code |= IH264E_INVALID_SVC_PARAMS | IH264E_INVALID_SPATIAL_RES_RATIO;
        }
    }

    if((u4_svc_comp_wd > SVCE_MAX_INP_DIM) || (u4_svc_comp_ht > SVCE_MAX_INP_DIM) ||
       ((u4_svc_comp_wd * u4_svc_comp_ht) > SVCE_MAX_INP_FRAME_SIZE) ||
       (u4_svc_comp_wd % 16 != 0) || (u4_svc_comp_ht % 16 != 0))
    {
        u4_error_code |= IH264E_INVALID_SVC_INPUT_DIMENSIONS;
    }

    /* Constraint from padding intrinsics */
    if((u4_svc_comp_wd - u4_inp_wd) % 16)
    {
        u4_error_code |= IH264E_INVALID_SVC_INPUT_DIMENSIONS;
    }

    /* Constraint from 420p to 420sp conversion */
    if((u4_svc_comp_ht - u4_inp_ht) % 4)
    {
        u4_error_code |= IH264E_INVALID_SVC_INPUT_DIMENSIONS;
    }

    au4_svc_wd[u1_num_spatial_layers - 1] = u4_svc_comp_wd;
    au4_svc_ht[u1_num_spatial_layers - 1] = u4_svc_comp_ht;

    for(i = (u1_num_spatial_layers - 1); i > 0; i--)
    {
        au4_svc_wd[i - 1] = au4_svc_wd[i] / d_scaling_factor;
        au4_svc_ht[i - 1] = au4_svc_ht[i] / d_scaling_factor;

        if((au4_svc_wd[i - 1] * d_scaling_factor != au4_svc_wd[i]) ||
           (au4_svc_ht[i - 1] * d_scaling_factor != au4_svc_ht[i]) ||
           (au4_svc_ht[i - 1] % 16 != 0) || (au4_svc_ht[i - 1] % 16 != 0))
        {
            u4_error_code |= IH264E_INVALID_SVC_INPUT_DIMENSIONS;
        }
    }

    return u4_error_code;
}

/**
*******************************************************************************
*
* @brief
*  Validates SVC input params
*
* @param[in] ps_cfg
*  Cfg parameters
*
* @returns  error code in conformance with 'IH264E_ERROR_T'
*
*******************************************************************************
*/
WORD32 isvce_svc_inp_params_validate(isvce_init_ip_t *ps_ip, isvce_cfg_params_t *ps_cfg)
{
    UWORD32 u4_error_code = isvce_svc_au_props_validate(&ps_ip->s_svc_inp_params, ps_ip->u4_wd,
                                                        ps_ip->u4_ht, ps_cfg->u4_wd, ps_cfg->u4_ht);

    if(ps_cfg->u4_enable_alt_ref)
    {
        u4_error_code |= IH264E_INVALID_ALT_REF_OPTION;
    }

    if(ps_cfg->u4_num_bframes)
    {
        u4_error_code |= IH264E_BFRAMES_NOT_SUPPORTED;
    }

    if(ps_cfg->e_slice_mode != IVE_SLICE_MODE_NONE)
    {
        u4_error_code |= IH264E_SLICE_TYPE_INPUT_INVALID;
    }

    if(ps_cfg->e_content_type != IV_PROGRESSIVE)
    {
        u4_error_code |= IH264E_CONTENT_TYPE_NOT_SUPPORTED;
    }

    if(ps_cfg->u4_weighted_prediction)
    {
        u4_error_code |= IH264E_WEIGHTED_PRED_NOT_SUPPORTED;
    }

    return u4_error_code;
}

/**
*******************************************************************************
*
* @brief
*  Validates SVC frame-level input params
*
* @param[in] ps_cfg
*  Cfg parameters
*
* @returns  error code in conformance with 'IH264E_ERROR_T'
*
*******************************************************************************
*/
WORD32 isvce_svc_frame_params_validate(
    rate_control_api_t *aps_rate_control_api[MAX_NUM_SPATIAL_LAYERS], UWORD8 u1_num_spatial_layers)
{
    WORD32 i;

    /* RC requires total bits in a second to fit int32_t */
    for(i = 0; i < u1_num_spatial_layers; i++)
    {
        if((((UWORD64) irc_get_bits_per_frame(aps_rate_control_api[i])) *
            irc_get_intra_frame_interval(aps_rate_control_api[i])) > ((UWORD64) INT32_MAX))
        {
            return IH264E_BITRATE_NOT_SUPPORTED;
        }
    }

    return IV_SUCCESS;
}

/**
*******************************************************************************
*
* @brief
*  Used to get reference picture buffer size for a given level and
*  and padding used
*
* @param[in] ps_svc_inp_params
*  Struct containing SVC specific input params
*
* @param[in] i4_pic_size
*  Number of luma samples (Width * Height)
*
* @param[in] i4_level
*  Level
*
* @param[in] i4_horz_pad
*  Total padding used in horizontal direction
*
* @param[in] i4_vert_pad
*  Total padding used in vertical direction
*
* @param[in] i4_num_ref_frames
*  Num Reference Frames
*
* @param[in] i4_num_reorder_frames
*  Num Reorder Frames
*
* @returns  Total picture buffer size
*
*******************************************************************************
*/
WORD32 isvce_get_total_svc_au_buf_size(svc_inp_params_t *ps_svc_inp_params, WORD32 i4_pic_size,
                                       WORD32 i4_level, WORD32 i4_horz_pad, WORD32 i4_vert_pad,
                                       WORD32 i4_num_ref_frames, WORD32 i4_num_reorder_frames)
{
    WORD32 i;
    WORD32 size;
    WORD32 num_luma_samples;
    WORD32 lvl_idx;
    WORD32 max_wd, min_ht;
    WORD32 num_samples;
    WORD32 max_num_bufs;

    WORD32 pad = MAX(i4_horz_pad, i4_vert_pad);
    DOUBLE d_svc_size_multiplier = 1;

    for(i = 1; i < ps_svc_inp_params->u1_num_spatial_layers; i++)
    {
        d_svc_size_multiplier += 1. / pow(ps_svc_inp_params->d_spatial_res_ratio, i);
    }

    /*
     * If i4_num_ref_frames and num_reorder_frmaes is specified
     * Use minimum value
     */
    max_num_bufs = (i4_num_ref_frames + i4_num_reorder_frames + MAX_CTXT_SETS +
                    ps_svc_inp_params->u1_num_temporal_layers);

    /* Get i4_level index */
    lvl_idx = ih264e_get_lvl_idx(i4_level);

    /* Maximum number of luma samples in a picture at given i4_level */
    num_luma_samples = gai4_ih264_max_luma_pic_size[lvl_idx];
    num_luma_samples = MAX(num_luma_samples, i4_pic_size);

    /* Account for chroma */
    num_samples = num_luma_samples * 3 / 2;

    /* Maximum width of luma samples in a picture at given i4_level */
    max_wd = gai4_ih264_max_wd_ht[lvl_idx];

    /* Minimum height of luma samples in a picture at given i4_level */
    min_ht = gai4_ih264_min_wd_ht[lvl_idx];

    /* Allocation is required for
     * (Wd + i4_horz_pad) * (Ht + i4_vert_pad) * (2 * max_dpb_size + 1)
     *
     * Above expanded as
     * ((Wd * Ht) + (i4_horz_pad * i4_vert_pad) + Wd * i4_vert_pad + Ht *
     * i4_horz_pad) * (2 * max_dpb_size + 1) (Wd * Ht) * (2 * max_dpb_size + 1) +
     * ((i4_horz_pad * i4_vert_pad) + Wd * i4_vert_pad + Ht * i4_horz_pad) * (2 *
     * max_dpb_size + 1) Now max_dpb_size increases with smaller Wd and Ht, but Wd
     * * ht * max_dpb_size will still be lesser or equal to max_wd * max_ht *
     * dpb_size
     *
     * In the above equation (Wd * Ht) * (2 * max_dpb_size + 1) is accounted by
     * using num_samples * (2 * max_dpb_size + 1) below
     *
     * For the padded area use MAX(i4_horz_pad, i4_vert_pad) as pad
     * ((pad * pad) + pad * (Wd + Ht)) * (2 * max_dpb_size + 1) has to accounted
     * from the above for padding
     *
     * Since Width and Height can change worst Wd + Ht is when One of the
     * dimensions is max and other is min So use max_wd and min_ht
     */

    /* Number of bytes in reference pictures */
    size = num_samples * max_num_bufs;

    /* Account for Spatial Layers */
    size = (WORD32) (size * d_svc_size_multiplier + 0.99);

    /* Account for padding area */
    size += ((pad * pad) + pad * (max_wd + min_ht)) * 3 / 2 * max_num_bufs *
            ps_svc_inp_params->u1_num_spatial_layers;

    size += ps_svc_inp_params->u1_num_spatial_layers * sizeof(yuv_buf_props_t);

    return size;
}

/**
*******************************************************************************
*
* @brief
*  Used to get size of buffers used for storing prediction data
*
* @param[in] ps_svc_inp_params
*  Struct containing SVC specific input params
*
* @param[in] i4_num_luma_samples
*  Number of luma samples (Width * Height)
*
* @returns  Size of buffers used for storing prediction data
*
*******************************************************************************
*/
UWORD32 isvce_get_total_svc_au_data_size(WORD32 i4_num_luma_samples, UWORD8 u1_num_spatial_layers,
                                         DOUBLE d_spatial_res_ratio)
{
    WORD32 i;

    UWORD32 u4_svc_au_data_size = 0;

    u4_svc_au_data_size += u1_num_spatial_layers * sizeof(svc_layer_data_t);

    for(i = 0; i < u1_num_spatial_layers; i++)
    {
        WORD32 i4_layer_luma_samples =
            ((DOUBLE) i4_num_luma_samples) / pow(pow(d_spatial_res_ratio, i), 2) + 0.99;
        WORD32 i4_num_mbs = i4_layer_luma_samples / (MB_SIZE * MB_SIZE);

        /* isvce_mb_info_t */
        u4_svc_au_data_size += i4_num_mbs * sizeof(isvce_mb_info_t);

        /* pu4_num_pus_in_mb */
        u4_svc_au_data_size += i4_num_mbs * sizeof(UWORD32);
    }

    return u4_svc_au_data_size;
}

/**
*******************************************************************************
*
* @brief Function to add buffers to SVC AU Data Store Manager
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @returns  error status
*
*******************************************************************************
*/
IH264E_ERROR_T isvce_svc_au_data_mgr_add_bufs(isvce_codec_t *ps_codec)
{
    IH264_ERROR_T ret;

    WORD32 i, j;
    UWORD8 *pu1_buf;

    svc_au_data_t *ps_svc_au_data = ps_codec->ps_svc_au_data_base;

    WORD32 i4_max_dpb_size = ps_codec->i4_ref_buf_cnt;
    WORD64 i8_alloc_mem_size = ps_codec->i4_svc_au_data_size;
    WORD32 i4_num_luma_samples = ALIGN16(ps_codec->s_cfg.u4_wd) * ALIGN16(ps_codec->s_cfg.u4_ht);
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;
    DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;

    ps_codec->ps_svc_au_data = ps_svc_au_data;
    pu1_buf = (UWORD8 *) ps_svc_au_data;
    pu1_buf += BUF_MGR_MAX_CNT * sizeof(ps_svc_au_data[0]);

    i8_alloc_mem_size -= (BUF_MGR_MAX_CNT * sizeof(ps_svc_au_data[0]));

    i = 0;

    while(i < i4_max_dpb_size)
    {
        ps_svc_au_data->ps_svc_layer_data = (svc_layer_data_t *) pu1_buf;
        pu1_buf += u1_num_spatial_layers * sizeof(ps_svc_au_data->ps_svc_layer_data[0]);
        i8_alloc_mem_size -= u1_num_spatial_layers * sizeof(ps_svc_au_data->ps_svc_layer_data[0]);

        for(j = u1_num_spatial_layers - 1; j >= 0; j--)
        {
            WORD32 i4_layer_luma_samples =
                ((DOUBLE) i4_num_luma_samples) /
                    pow(pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j), 2) +
                0.99;
            WORD32 i4_num_mbs = i4_layer_luma_samples / (MB_SIZE * MB_SIZE);

            ps_svc_au_data->ps_svc_layer_data[j].pu4_num_pus_in_mb = (UWORD32 *) pu1_buf;
            pu1_buf +=
                i4_num_mbs * sizeof(ps_svc_au_data->ps_svc_layer_data[j].pu4_num_pus_in_mb[0]);
            i8_alloc_mem_size -=
                i4_num_mbs * sizeof(ps_svc_au_data->ps_svc_layer_data[j].pu4_num_pus_in_mb[0]);

            ps_svc_au_data->ps_svc_layer_data[j].ps_mb_info = (isvce_mb_info_t *) pu1_buf;
            pu1_buf += i4_num_mbs * sizeof(ps_svc_au_data->ps_svc_layer_data[j].ps_mb_info[0]);
            i8_alloc_mem_size -=
                i4_num_mbs * sizeof(ps_svc_au_data->ps_svc_layer_data[j].ps_mb_info[0]);

            ASSERT(i8_alloc_mem_size >= 0);
        }

        if(i8_alloc_mem_size < 0)
        {
            ps_codec->i4_error_code = IH264E_INSUFFICIENT_MEM_MVBANK;

            return IH264E_INSUFFICIENT_MEM_MVBANK;
        }

        ret =
            ih264_buf_mgr_add((buf_mgr_t *) ps_codec->pv_svc_au_data_store_mgr, ps_svc_au_data, i);

        if(IH264_SUCCESS != ret)
        {
            ps_codec->i4_error_code = IH264E_BUF_MGR_ERROR;

            return IH264E_BUF_MGR_ERROR;
        }

        ps_svc_au_data++;
        i++;
    }

    return IH264E_SUCCESS;
}

/**
*******************************************************************************
*
* @brief
*  Function to initialize svc_au_buf_t structs add au buffers to
*  buffer manager in case of non-shared mode
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @returns  error status
*
*******************************************************************************
*/
IH264E_ERROR_T isvce_svc_au_buf_mgr_add_bufs(isvce_codec_t *ps_codec)
{
    WORD32 i, j;
    WORD32 buf_ret;

    svc_au_buf_t *ps_pic_buf = ps_codec->ps_pic_buf;

    IH264E_ERROR_T ret = IH264E_SUCCESS;

    WORD32 i4_max_dpb_size = ps_codec->i4_ref_buf_cnt;
    WORD64 i8_alloc_mem_size =
        ps_codec->i4_total_pic_buf_size - BUF_MGR_MAX_CNT * sizeof(ps_pic_buf[0]);
    UWORD8 *pu1_buf = (UWORD8 *) ps_codec->ps_pic_buf;
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;
    DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;

    pu1_buf += BUF_MGR_MAX_CNT * sizeof(svc_au_buf_t);

    for(i = 0; i < i4_max_dpb_size; i++)
    {
        WORD32 i4_total_fpel_mem_size = 0;

        ps_pic_buf->ps_layer_yuv_buf_props = (yuv_buf_props_t *) pu1_buf;
        pu1_buf += u1_num_spatial_layers * sizeof(ps_pic_buf->ps_layer_yuv_buf_props[0]);
        i8_alloc_mem_size -= u1_num_spatial_layers * sizeof(ps_pic_buf->ps_layer_yuv_buf_props[0]);

        if(i8_alloc_mem_size < 0)
        {
            ps_codec->i4_error_code = IH264E_INSUFFICIENT_MEM_PICBUF;
            return IH264E_INSUFFICIENT_MEM_PICBUF;
        }

        for(j = u1_num_spatial_layers - 1; j >= 0; j--)
        {
            WORD32 i4_layer_luma_wd = ((DOUBLE) ps_codec->s_cfg.u4_wd /
                                       pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                                      0.99;
            WORD32 i4_layer_luma_ht = ((DOUBLE) ps_codec->s_cfg.u4_ht /
                                       pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) +
                                      0.99;
            WORD32 i4_layer_luma_samples =
                (ALIGN16(i4_layer_luma_wd) + PAD_WD) * (i4_layer_luma_ht + PAD_HT);
            WORD32 i4_layer_uv_wd = i4_layer_luma_wd;
            WORD32 i4_layer_uv_ht = i4_layer_luma_ht / 2.0 + 0.99;
            WORD32 i4_layer_uv_samples =
                (ALIGN16(i4_layer_uv_wd) + PAD_WD) * (i4_layer_uv_ht + PAD_HT);

            ps_pic_buf->ps_layer_yuv_buf_props[j].as_component_bufs[0].i4_data_stride =
                ALIGN16(i4_layer_luma_wd) + PAD_WD;
            ps_pic_buf->ps_layer_yuv_buf_props[j].as_component_bufs[0].pv_data =
                pu1_buf +
                ps_pic_buf->ps_layer_yuv_buf_props[j].as_component_bufs[0].i4_data_stride *
                    PAD_TOP +
                PAD_LEFT;

            pu1_buf += i4_layer_luma_samples;

            ps_pic_buf->ps_layer_yuv_buf_props[j].as_component_bufs[1].i4_data_stride =
                ALIGN16(i4_layer_uv_wd) + PAD_WD;
            ps_pic_buf->ps_layer_yuv_buf_props[j].as_component_bufs[1].pv_data =
                pu1_buf +
                ps_pic_buf->ps_layer_yuv_buf_props[j].as_component_bufs[1].i4_data_stride *
                    (PAD_TOP / 2) +
                PAD_LEFT;

            pu1_buf += i4_layer_uv_samples;

            ps_pic_buf->ps_layer_yuv_buf_props[j].u4_width = i4_layer_luma_wd;
            ps_pic_buf->ps_layer_yuv_buf_props[j].u4_height = i4_layer_luma_ht;
            ps_pic_buf->ps_layer_yuv_buf_props[j].u1_bit_depth = 8;
            ps_pic_buf->ps_layer_yuv_buf_props[j].e_color_format = IV_YUV_420SP_UV;

            i8_alloc_mem_size -= i4_layer_luma_samples + i4_layer_uv_samples;
            i4_total_fpel_mem_size += i4_layer_luma_samples + i4_layer_uv_samples;

            if(i8_alloc_mem_size < 0)
            {
                ps_codec->i4_error_code = IH264E_INSUFFICIENT_MEM_PICBUF;
                return IH264E_INSUFFICIENT_MEM_PICBUF;
            }
        }

        buf_ret = ih264_buf_mgr_add((buf_mgr_t *) ps_codec->pv_ref_buf_mgr, ps_pic_buf, i);

        if(0 != buf_ret)
        {
            ps_codec->i4_error_code = IH264E_BUF_MGR_ERROR;
            return IH264E_BUF_MGR_ERROR;
        }

        pu1_buf += (HPEL_PLANES_CNT - 1) * i4_total_fpel_mem_size;
        ps_pic_buf++;
    }

    return ret;
}

/**
*******************************************************************************
*
* @brief
*  Returns size of buffers for storing SVC input data
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
UWORD32 isvce_get_svc_inp_buf_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                   UWORD32 u4_wd, UWORD32 u4_ht)
{
    padding_dims_t s_pad_dims;

    UWORD32 i;
    UWORD8 u1_filter_padding_size_x, u1_filter_padding_size_y;

    UWORD32 u4_size = 0;

    isvce_get_downscaler_padding_dims(&s_pad_dims);

    u1_filter_padding_size_x = s_pad_dims.u1_left_pad_size + s_pad_dims.u1_right_pad_size;

    u1_filter_padding_size_y = s_pad_dims.u1_top_pad_size + s_pad_dims.u1_bottom_pad_size;

    for(i = 0; i < u1_num_spatial_layers; i++)
    {
        WORD32 i4_layer_luma_wd = ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, i)) + 0.99;
        WORD32 i4_layer_luma_ht = ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, i)) + 0.99;
        WORD32 i4_layer_luma_samples =
            (ALIGN16(i4_layer_luma_wd) + PAD_WD + u1_filter_padding_size_x) *
            (i4_layer_luma_ht + PAD_HT + u1_filter_padding_size_y);
        WORD32 i4_layer_uv_wd = i4_layer_luma_wd;
        WORD32 i4_layer_uv_ht = i4_layer_luma_ht / 2.0 + 0.99;
        /* u1_filter_padding_size_x * 2 because U and V
        both need same amount of padding */
        WORD32 i4_layer_uv_samples =
            (ALIGN16(i4_layer_uv_wd) + PAD_WD + u1_filter_padding_size_x * 2) *
            (i4_layer_uv_ht + PAD_HT + u1_filter_padding_size_y);

        u4_size += (i4_layer_luma_samples + i4_layer_uv_samples) * sizeof(UWORD8);
    }

    return SVC_MAX_NUM_INP_FRAMES * u4_size;
}

/**
*******************************************************************************
*
* @brief
*  Function to initialize svc input buffers
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @param[in] ps_mem_rec
*  Pointer to memory allocated for input buffers
*
*******************************************************************************
*/
void isvce_svc_inp_buf_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec)
{
    padding_dims_t s_pad_dims;

    WORD32 i, j;
    UWORD8 u1_filter_padding_size_x, u1_filter_padding_size_y;

    DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;
    UWORD32 u4_wd = ps_codec->s_cfg.u4_wd;
    UWORD32 u4_ht = ps_codec->s_cfg.u4_ht;
    UWORD8 *pu1_buf = ps_mem_rec->pv_base;
    WORD64 i8_alloc_mem_size =
        isvce_get_svc_inp_buf_size(u1_num_spatial_layers, d_spatial_res_ratio, u4_wd, u4_ht);

    isvce_get_downscaler_padding_dims(&s_pad_dims);

    u1_filter_padding_size_x = s_pad_dims.u1_left_pad_size + s_pad_dims.u1_right_pad_size;

    u1_filter_padding_size_y = s_pad_dims.u1_top_pad_size + s_pad_dims.u1_bottom_pad_size;

    for(i = 0; i < SVC_MAX_NUM_INP_FRAMES; i++)
    {
        ps_codec->as_inp_list[i].s_svc_params = ps_codec->s_cfg.s_svc_params;

        for(j = u1_num_spatial_layers - 1; j >= 0; j--)
        {
            WORD32 i4_layer_luma_wd =
                ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) + 0.99;
            WORD32 i4_layer_luma_ht =
                ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - j)) + 0.99;
            WORD32 i4_layer_luma_samples =
                (ALIGN16(i4_layer_luma_wd) + PAD_WD + u1_filter_padding_size_x) *
                (i4_layer_luma_ht + PAD_HT + u1_filter_padding_size_y);
            WORD32 i4_layer_uv_wd = i4_layer_luma_wd;
            WORD32 i4_layer_uv_ht = i4_layer_luma_ht / 2.0 + 0.99;
            /* u1_filter_padding_size_x * 2 because U and V
            both need same amount of padding */
            WORD32 i4_layer_uv_samples =
                (ALIGN16(i4_layer_uv_wd) + PAD_WD + u1_filter_padding_size_x * 2) *
                (i4_layer_uv_ht + PAD_HT + u1_filter_padding_size_y);

            ps_codec->as_inp_list[i].as_layer_yuv_buf_props[j].as_component_bufs[Y].i4_data_stride =
                ALIGN16(i4_layer_luma_wd) + PAD_WD + u1_filter_padding_size_x;
            ps_codec->as_inp_list[i].as_layer_yuv_buf_props[j].as_component_bufs[Y].pv_data =
                pu1_buf +
                ps_codec->as_inp_list[i]
                        .as_layer_yuv_buf_props[j]
                        .as_component_bufs[Y]
                        .i4_data_stride *
                    (PAD_TOP + s_pad_dims.u1_top_pad_size) +
                (PAD_LEFT + s_pad_dims.u1_left_pad_size);
            pu1_buf += i4_layer_luma_samples * sizeof(UWORD8);
            i8_alloc_mem_size -= i4_layer_luma_samples * sizeof(UWORD8);

            ps_codec->as_inp_list[i]
                .as_layer_yuv_buf_props[j]
                .as_component_bufs[UV]
                .i4_data_stride = ALIGN16(i4_layer_uv_wd) + PAD_WD + u1_filter_padding_size_x * 2;
            ps_codec->as_inp_list[i].as_layer_yuv_buf_props[j].as_component_bufs[UV].pv_data =
                pu1_buf +
                ps_codec->as_inp_list[i]
                        .as_layer_yuv_buf_props[j]
                        .as_component_bufs[UV]
                        .i4_data_stride *
                    (PAD_TOP + s_pad_dims.u1_top_pad_size) +
                (PAD_LEFT + s_pad_dims.u1_left_pad_size * 2);
            pu1_buf += i4_layer_uv_samples * sizeof(UWORD8);
            i8_alloc_mem_size -= i4_layer_uv_samples * sizeof(UWORD8);

            /* Chroma is always stored interleaved */
            ps_codec->as_inp_list[i].as_layer_yuv_buf_props[j].as_component_bufs[V].pv_data = NULL;

            ps_codec->as_inp_list[i].as_layer_yuv_buf_props[j].u1_bit_depth = 8;
            ps_codec->as_inp_list[i].as_layer_yuv_buf_props[j].e_color_format = IV_YUV_420SP_UV;
            ps_codec->as_inp_list[i].as_layer_yuv_buf_props[j].u4_width = i4_layer_luma_wd;
            ps_codec->as_inp_list[i].as_layer_yuv_buf_props[j].u4_height = i4_layer_luma_ht;

            ASSERT(i8_alloc_mem_size >= 0);
        }
    }
}

void isvce_init_svc_dimension(isvce_inp_buf_t *ps_inp)
{
    WORD32 i;

    UWORD8 u1_num_spatial_layers = ps_inp->s_svc_params.u1_num_spatial_layers;
    DOUBLE d_spatial_res_ratio = ps_inp->s_svc_params.d_spatial_res_ratio;
    UWORD32 u4_wd = ps_inp->s_inp_props.s_raw_buf.au4_wd[Y];
    UWORD32 u4_ht = ps_inp->s_inp_props.s_raw_buf.au4_ht[Y];

    for(i = 0; i < u1_num_spatial_layers; i++)
    {
        ps_inp->as_layer_yuv_buf_props[i].u4_width =
            ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
        ps_inp->as_layer_yuv_buf_props[i].u4_height =
            ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
    }
}

/**
*******************************************************************************
*
* @brief
*  Pads input buf as assumed by the downscaler filter
*
* @param[in] ps_codec
*  Pointer to codec ctxt
*
* @param[in] ps_inp
*  Pointer to svc input buffer
*
* @param[in] u1_svc_layer_index
*  SVC layer index of the buffer
*
*******************************************************************************
*/

static void isvce_pad_buf_for_filtering(isvce_codec_t *ps_codec, isvce_inp_buf_t *ps_inp,
                                        UWORD8 u1_svc_layer_index)
{
    padding_dims_t s_pad_dims;

    UWORD8 *pu1_buf;
    UWORD32 u4_buf_width, u4_buf_height;

    UWORD8 u1_pad_left_size;
    UWORD8 u1_pad_right_size;
    UWORD8 u1_pad_top_size;
    UWORD8 u1_pad_bottom_size;
    UWORD8 u1_filter_padding_size_x;
    UWORD8 u1_filter_padding_size_chroma_x;

    ASSERT(ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].e_color_format == IV_YUV_420SP_UV);

    isvce_get_downscaler_padding_dims(&s_pad_dims);

    u1_pad_left_size = s_pad_dims.u1_left_pad_size;
    u1_pad_right_size = s_pad_dims.u1_right_pad_size;
    u1_pad_top_size = s_pad_dims.u1_top_pad_size;
    u1_pad_bottom_size = s_pad_dims.u1_bottom_pad_size;
    u1_filter_padding_size_x = u1_pad_left_size + u1_pad_right_size;
    u1_filter_padding_size_chroma_x = u1_filter_padding_size_x * 2;

    u4_buf_width = ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].u4_width;

    u4_buf_height = ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].u4_height;

    pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index]
                              .as_component_bufs[0]
                              .pv_data);

    ps_codec->pf_pad_left_luma(
        pu1_buf,
        ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[0].i4_data_stride,
        u4_buf_height, u1_pad_left_size);

    pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index]
                              .as_component_bufs[0]
                              .pv_data);

    pu1_buf += u4_buf_width;

    ps_codec->pf_pad_right_luma(
        pu1_buf,
        ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[0].i4_data_stride,
        u4_buf_height, u1_pad_right_size);

    pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index]
                              .as_component_bufs[1]
                              .pv_data);

    ps_codec->pf_pad_left_chroma(
        pu1_buf,
        ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[1].i4_data_stride,
        u4_buf_height / 2, u1_pad_left_size * 2);

    pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index]
                              .as_component_bufs[1]
                              .pv_data);

    pu1_buf += u4_buf_width;

    ps_codec->pf_pad_right_chroma(
        pu1_buf,
        ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[1].i4_data_stride,
        u4_buf_height / 2, u1_pad_right_size * 2);

    pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index]
                              .as_component_bufs[0]
                              .pv_data) -
              u1_pad_left_size;

    ps_codec->pf_pad_top(
        pu1_buf,
        ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[0].i4_data_stride,
        (u4_buf_width + u1_filter_padding_size_x), u1_pad_top_size);

    pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index]
                              .as_component_bufs[0]
                              .pv_data) -
              u1_pad_left_size;

    pu1_buf +=
        (u4_buf_height *
         ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[0].i4_data_stride);

    ps_codec->pf_pad_bottom(
        pu1_buf,
        ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[0].i4_data_stride,
        (u4_buf_width + u1_filter_padding_size_x), u1_pad_bottom_size);

    pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index]
                              .as_component_bufs[1]
                              .pv_data) -
              u1_pad_left_size * 2;

    ps_codec->pf_pad_top(
        pu1_buf,
        ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[1].i4_data_stride,
        (u4_buf_width + u1_filter_padding_size_chroma_x), u1_pad_top_size);

    pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index]
                              .as_component_bufs[1]
                              .pv_data) -
              u1_pad_left_size * 2;

    pu1_buf +=
        ((u4_buf_height / 2) *
         ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[1].i4_data_stride);

    ps_codec->pf_pad_bottom(
        pu1_buf,
        ps_inp->as_layer_yuv_buf_props[u1_svc_layer_index].as_component_bufs[1].i4_data_stride,
        (u4_buf_width + u1_filter_padding_size_chroma_x), u1_pad_bottom_size);
}

/**
*******************************************************************************
*
* @brief
*  Pads raw input to satisfy SVC compliant input dimensions
*
* @param[in] ps_codec
*  Pointer to codec ctxt
*
* @param[in] ps_inp
*  Pointer to svc input buffer
*
*******************************************************************************
*/

static void isvce_pad_input_to_svc_compliant_dims(isvce_codec_t *ps_codec, isvce_inp_buf_t *ps_inp)
{
    UWORD8 *pu1_buf;
    UWORD32 u4_raw_input_wd, u4_raw_input_ht, u4_padded_width, u4_padded_height, u4_width_delta,
        u4_height_delta;
    UWORD8 u1_num_layers = ps_inp->s_svc_params.u1_num_spatial_layers;

    ASSERT(ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].e_color_format == IV_YUV_420SP_UV);

    u4_padded_width = ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].u4_width;
    u4_padded_height = ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].u4_height;
    u4_raw_input_wd = ps_inp->s_inp_props.s_raw_buf.au4_wd[0];
    u4_raw_input_ht = ps_inp->s_inp_props.s_raw_buf.au4_ht[0];
    u4_width_delta = u4_padded_width - u4_raw_input_wd;
    u4_height_delta = u4_padded_height - u4_raw_input_ht;

    ASSERT(!(u4_width_delta & 1));
    ASSERT(!(u4_height_delta & 1));

    if(u4_width_delta)
    {
        pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                  .as_component_bufs[0]
                                  .pv_data);

        pu1_buf += ((u4_width_delta / 2) + (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                                .as_component_bufs[0]
                                                .i4_data_stride) *
                                               (u4_height_delta / 2));

        ps_codec->pf_pad_left_luma(
            pu1_buf,
            ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[0].i4_data_stride,
            u4_padded_height, u4_width_delta / 2);

        pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                  .as_component_bufs[0]
                                  .pv_data);

        pu1_buf += ((u4_width_delta / 2) + (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                                .as_component_bufs[0]
                                                .i4_data_stride) *
                                               (u4_height_delta / 2));

        pu1_buf += u4_raw_input_wd;

        ps_codec->pf_pad_right_luma(
            pu1_buf,
            ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[0].i4_data_stride,
            u4_padded_height, u4_width_delta / 2);

        pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                  .as_component_bufs[1]
                                  .pv_data);

        pu1_buf += ((u4_width_delta / 2) + (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                                .as_component_bufs[1]
                                                .i4_data_stride) *
                                               (u4_height_delta / 4));

        ps_codec->pf_pad_left_chroma(
            pu1_buf,
            ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[1].i4_data_stride,
            u4_padded_height / 2, u4_width_delta / 2);

        pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                  .as_component_bufs[1]
                                  .pv_data);

        pu1_buf += ((u4_width_delta / 2) + (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                                .as_component_bufs[1]
                                                .i4_data_stride) *
                                               (u4_height_delta / 4));

        pu1_buf += u4_raw_input_wd;

        ps_codec->pf_pad_right_chroma(
            pu1_buf,
            ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[1].i4_data_stride,
            u4_padded_height / 2, u4_width_delta / 2);
    }

    if(u4_height_delta)
    {
        pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                  .as_component_bufs[0]
                                  .pv_data);

        pu1_buf += ((ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                         .as_component_bufs[0]
                         .i4_data_stride) *
                    (u4_height_delta / 2));

        ps_codec->pf_pad_top(
            pu1_buf,
            ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[0].i4_data_stride,
            u4_padded_width, u4_height_delta / 2);

        pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                  .as_component_bufs[0]
                                  .pv_data);

        pu1_buf += ((ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                         .as_component_bufs[0]
                         .i4_data_stride) *
                    (u4_height_delta / 2));

        pu1_buf +=
            (u4_raw_input_ht *
             ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[0].i4_data_stride);

        ps_codec->pf_pad_bottom(
            pu1_buf,
            ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[0].i4_data_stride,
            u4_padded_width, u4_height_delta / 2);

        pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                  .as_component_bufs[1]
                                  .pv_data);

        pu1_buf += ((ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                         .as_component_bufs[1]
                         .i4_data_stride) *
                    (u4_height_delta / 4));

        ps_codec->pf_pad_top(
            pu1_buf,
            ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[1].i4_data_stride,
            u4_padded_width, u4_height_delta / 4);

        pu1_buf = (UWORD8 *) (ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                                  .as_component_bufs[1]
                                  .pv_data);

        pu1_buf += ((ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1]
                         .as_component_bufs[1]
                         .i4_data_stride) *
                    (u4_height_delta / 4));

        pu1_buf +=
            ((u4_raw_input_ht / 2) *
             ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[1].i4_data_stride);

        ps_codec->pf_pad_bottom(
            pu1_buf,
            ps_inp->as_layer_yuv_buf_props[u1_num_layers - 1].as_component_bufs[1].i4_data_stride,
            u4_padded_width, u4_height_delta / 4);
    }
}

/**
*******************************************************************************
*
* @brief
*  Format conversion and downsampling for deriving spatial layer inputs
*
* @param[in] ps_inp
*  Pointer to input buffer
*
*******************************************************************************
*/
void isvce_svc_inp_buf_populate(isvce_codec_t *ps_codec, isvce_inp_buf_t *ps_inp)
{
    yuv_buf_props_t s_src_buf_props, s_dst_buf_props;

    UWORD32 i;
    UWORD32 u4_blk_x, u4_blk_y;
    UWORD8 *pu1_planar_y, *pu1_planar_u, *pu1_planar_v, *pu1_semi_planar_y, *pu1_semi_planar_uv;
    UWORD8 *pu1_src_luma, *pu1_src_chroma, *pu1_dst_luma, *pu1_dst_chroma;
    UWORD32 u4_num_blocks_x, u4_num_blocks_y;
    UWORD32 u4_scaled_block_wd, u4_scaled_block_ht;
    UWORD32 u4_blk_wd_luma, u4_blk_ht_luma;

    downscaler_ctxt_t *ps_scaler = &ps_codec->s_scaler;
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    mem_fxns_t *ps_mem_fxns = &ps_isa_dependent_fxns->s_mem_fxns;

    const UWORD8 u1_num_yuv_components_420sp = NUM_SP_COMPONENTS;
    UWORD8 u1_num_spatial_layers = ps_inp->s_svc_params.u1_num_spatial_layers;
    UWORD32 u4_padded_width = ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1].u4_width;
    UWORD32 u4_padded_height = ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1].u4_height;
    UWORD32 u4_raw_input_wd = ps_inp->s_inp_props.s_raw_buf.au4_wd[0];
    UWORD32 u4_raw_input_ht = ps_inp->s_inp_props.s_raw_buf.au4_ht[0];
    UWORD32 u4_width_delta = u4_padded_width - u4_raw_input_wd;
    UWORD32 u4_height_delta = u4_padded_height - u4_raw_input_ht;

    ASSERT(!(u4_width_delta & 1));
    ASSERT(!(u4_height_delta & 1));

    ASSERT((ps_inp->s_inp_props.s_raw_buf.e_color_fmt == IV_YUV_420P) ||
           (ps_inp->s_inp_props.s_raw_buf.e_color_fmt == IV_YUV_420SP_UV));

    /* Check is input is valid */
    if(!(ps_inp->s_inp_props.s_raw_buf.apv_bufs[0]))
    {
        ASSERT(0);

        return;
    }

    /* Convert the input into semi-planar in case of other formats */
    if(ps_inp->s_inp_props.s_raw_buf.e_color_fmt == IV_YUV_420P)
    {
        pu1_planar_y = (UWORD8 *) ps_inp->s_inp_props.s_raw_buf.apv_bufs[0];
        pu1_planar_u = (UWORD8 *) ps_inp->s_inp_props.s_raw_buf.apv_bufs[1];
        pu1_planar_v = (UWORD8 *) ps_inp->s_inp_props.s_raw_buf.apv_bufs[2];

        pu1_semi_planar_y = (UWORD8 *) ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                                .as_component_bufs[0]
                                .pv_data;

        pu1_semi_planar_uv = (UWORD8 *) ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                                 .as_component_bufs[1]
                                 .pv_data;

        pu1_semi_planar_y +=
            ((u4_width_delta / 2) + (ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                                         .as_component_bufs[0]
                                         .i4_data_stride) *
                                        (u4_height_delta / 2));

        pu1_semi_planar_uv +=
            ((u4_width_delta / 2) + (ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                                         .as_component_bufs[1]
                                         .i4_data_stride) *
                                        (u4_height_delta / 4));

        ps_codec->pf_ih264e_conv_420p_to_420sp(
            pu1_planar_y, pu1_planar_u, pu1_planar_v, pu1_semi_planar_y, pu1_semi_planar_uv,
            ps_inp->s_inp_props.s_raw_buf.au4_ht[0], ps_inp->s_inp_props.s_raw_buf.au4_wd[0],
            ps_inp->s_inp_props.s_raw_buf.au4_strd[0], ps_inp->s_inp_props.s_raw_buf.au4_strd[1],
            ps_inp->s_inp_props.s_raw_buf.au4_strd[2],
            ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                .as_component_bufs[0]
                .i4_data_stride,
            ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                .as_component_bufs[1]
                .i4_data_stride,
            0);
    }
    else
    {
        UWORD32 u4_wd, u4_ht;
        UWORD8 u1_comp;
        UWORD32 au4_arr_dims[4];
        UWORD8 *pu1_src, *pu1_dst;

        au4_arr_dims[0] = ps_inp->s_inp_props.s_raw_buf.au4_wd[0];
        au4_arr_dims[1] = ps_inp->s_inp_props.s_raw_buf.au4_ht[0];
        au4_arr_dims[2] = ps_inp->s_inp_props.s_raw_buf.au4_wd[1];
        au4_arr_dims[3] = ps_inp->s_inp_props.s_raw_buf.au4_ht[1];

        for(u1_comp = 0; u1_comp < u1_num_yuv_components_420sp; u1_comp++)
        {
            u4_wd = au4_arr_dims[u1_comp * 2];
            u4_ht = au4_arr_dims[(u1_comp * 2) + 1];

            pu1_dst = (UWORD8 *) ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                          .as_component_bufs[u1_comp]
                          .pv_data;

            pu1_dst +=
                ((u4_width_delta / 2) + (ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                                             .as_component_bufs[u1_comp]
                                             .i4_data_stride) *
                                            ((u4_height_delta / 2) / (u1_comp + 1)));

            pu1_src = ps_inp->s_inp_props.s_raw_buf.apv_bufs[u1_comp];

            ps_mem_fxns->pf_copy_2d(pu1_dst,
                                    ps_inp->as_layer_yuv_buf_props[u1_num_spatial_layers - 1]
                                        .as_component_bufs[u1_comp]
                                        .i4_data_stride,
                                    pu1_src, ps_inp->s_inp_props.s_raw_buf.au4_strd[u1_comp], u4_wd,
                                    u4_ht);
        }
    }

    /* Padding input to satisfy SVC constraints */
    isvce_pad_input_to_svc_compliant_dims(ps_codec, ps_inp);

    /* Downscaling */
    for(i = u1_num_spatial_layers - 1; i > 0; i--)
    {
        const UWORD32 u4_default_scaled_blk_wd =
            gu4_downscaler_blk_size / ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio + 0.5;
        const UWORD32 u4_default_scaled_blk_ht =
            gu4_downscaler_blk_size / ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio + 0.5;

        isvce_pad_buf_for_filtering(ps_codec, ps_inp, i);

        s_src_buf_props = ps_inp->as_layer_yuv_buf_props[i];
        s_dst_buf_props = ps_inp->as_layer_yuv_buf_props[i - 1];

        u4_num_blocks_x =
            (s_src_buf_props.u4_width + (gu4_downscaler_blk_size - 1)) / gu4_downscaler_blk_size;

        u4_num_blocks_y =
            (s_src_buf_props.u4_height + (gu4_downscaler_blk_size - 1)) / gu4_downscaler_blk_size;

        pu1_src_luma = (UWORD8 *) s_src_buf_props.as_component_bufs[Y].pv_data;
        pu1_src_chroma = (UWORD8 *) s_src_buf_props.as_component_bufs[U].pv_data;
        pu1_dst_luma = (UWORD8 *) s_dst_buf_props.as_component_bufs[Y].pv_data;
        pu1_dst_chroma = (UWORD8 *) s_dst_buf_props.as_component_bufs[U].pv_data;

        for(u4_blk_x = 0; u4_blk_x < u4_num_blocks_x; u4_blk_x++)
        {
            for(u4_blk_y = 0; u4_blk_y < u4_num_blocks_y; u4_blk_y++)
            {
                u4_blk_wd_luma = isvce_get_downscaler_blk_dims(s_src_buf_props.u4_width, u4_blk_x,
                                                               gu4_downscaler_blk_size);

                u4_blk_ht_luma = isvce_get_downscaler_blk_dims(s_src_buf_props.u4_height, u4_blk_y,
                                                               gu4_downscaler_blk_size);

                u4_scaled_block_wd = isvce_get_downscaler_blk_dims(
                    s_dst_buf_props.u4_width, u4_blk_x, u4_default_scaled_blk_wd);

                u4_scaled_block_ht = isvce_get_downscaler_blk_dims(
                    s_dst_buf_props.u4_height, u4_blk_y, u4_default_scaled_blk_ht);

                s_src_buf_props.as_component_bufs[Y].pv_data =
                    pu1_src_luma + (u4_blk_x * gu4_downscaler_blk_size +
                                    u4_blk_y * gu4_downscaler_blk_size *
                                        s_src_buf_props.as_component_bufs[Y].i4_data_stride);

                s_src_buf_props.as_component_bufs[U].pv_data =
                    pu1_src_chroma + (u4_blk_x * gu4_downscaler_blk_size +
                                      u4_blk_y * (gu4_downscaler_blk_size / 2) *
                                          s_src_buf_props.as_component_bufs[U].i4_data_stride);

                s_dst_buf_props.as_component_bufs[Y].pv_data =
                    pu1_dst_luma + (u4_blk_x * u4_default_scaled_blk_wd +
                                    u4_blk_y * u4_default_scaled_blk_ht *
                                        s_dst_buf_props.as_component_bufs[Y].i4_data_stride);

                s_dst_buf_props.as_component_bufs[U].pv_data =
                    pu1_dst_chroma + (u4_blk_x * u4_default_scaled_blk_wd +
                                      u4_blk_y * (u4_default_scaled_blk_ht / 2) *
                                          s_dst_buf_props.as_component_bufs[U].i4_data_stride);

                ASSERT(!(u4_scaled_block_wd & 1));
                ASSERT(!(u4_scaled_block_ht & 1));

                isvce_process_downscaler(ps_scaler, &s_src_buf_props, &s_dst_buf_props,
                                         u4_blk_wd_luma, u4_blk_ht_luma);
            }
        }
    }

    UNUSED(u4_scaled_block_wd);
    UNUSED(u4_scaled_block_ht);
}

/**
*******************************************************************************
*
* @brief
*  calculates the greatest common divisor between the two parameters.
*
*******************************************************************************
*/

static DOUBLE isvce_get_GCD(DOUBLE a, DOUBLE b)
{
    if(b == 0)
    {
        return a;
    }

    return isvce_get_GCD(b, fmod(a, b));
}

/**
*******************************************************************************
*
* @brief
*  calculates the least common multiple between the two parameters
*
*******************************************************************************
*/

static DOUBLE isvce_get_LCM(DOUBLE a, DOUBLE b) { return (a / isvce_get_GCD(a, b)) * b; }

/**
*******************************************************************************
*
* @brief
*  sets the width and height in config structure to SVC compliant width and
*   height
*
* @param[in] ps_cfg
*  Pointer to config struct
*
* @param[in] u4_app_wd
*  width of the YUV as read by the app
*
* @param[in] u4_app_ht
*  height of the YUV as read by the app
*
*******************************************************************************
*/

void isvce_get_svc_compliant_dimensions(UWORD8 u1_num_spatial_layers, DOUBLE d_scaling_factor,
                                        UWORD32 u4_wd, UWORD32 u4_ht, UWORD32 *pu4_svc_comp_wd,
                                        UWORD32 *pu4_svc_comp_ht)
{
    DOUBLE d_scaling_factor_power_num_layers_minus1 = 0;
    UWORD32 u4_constraint_offset = 0;

    d_scaling_factor_power_num_layers_minus1 = pow(d_scaling_factor, u1_num_spatial_layers - 1);

    if(fmod(16, d_scaling_factor_power_num_layers_minus1))
    {
        u4_constraint_offset =
            (UWORD32) isvce_get_LCM(16, d_scaling_factor_power_num_layers_minus1);
    }
    else
    {
        u4_constraint_offset = (UWORD32) (16 * d_scaling_factor_power_num_layers_minus1);
    }

    if(u4_wd % u4_constraint_offset)
    {
        *pu4_svc_comp_wd = u4_wd - ((u4_wd) % u4_constraint_offset) + u4_constraint_offset;
    }
    else
    {
        *pu4_svc_comp_wd = u4_wd;
    }

    if(u4_ht % u4_constraint_offset)
    {
        *pu4_svc_comp_ht = u4_ht - ((u4_ht) % u4_constraint_offset) + u4_constraint_offset;
    }
    else
    {
        *pu4_svc_comp_ht = u4_ht;
    }
}

/**
*******************************************************************************
*
* @brief
*  Returns size of buffers for storing SVC layer nbr info
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
* @returns  Size of buffers
*
*******************************************************************************
*/
UWORD32 isvce_get_svc_nbr_info_buf_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                        UWORD32 u4_wd, UWORD32 u4_ht)
{
    UWORD32 i;

    UWORD32 u4_size = 0;

    ASSERT(1 == MAX_CTXT_SETS);

    u4_size += MAX_PROCESS_CTXT * u1_num_spatial_layers * sizeof(nbr_info_t);

    for(i = 0; i < u1_num_spatial_layers; i++)
    {
        WORD32 i4_layer_luma_wd = ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, i)) + 0.99;
        WORD32 i4_layer_luma_ht = ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, i)) + 0.99;
        WORD32 i4_num_mbs_in_row = i4_layer_luma_wd / MB_SIZE;
        WORD32 i4_num_mbs_in_col = i4_layer_luma_ht / MB_SIZE;

        /* ps_top_row_mb_info */
        u4_size += (i4_num_mbs_in_row + 1) * i4_num_mbs_in_col * sizeof(isvce_mb_info_t);

        /* ps_left_mb_info */
        u4_size += MAX_PROCESS_CTXT * sizeof(isvce_mb_info_t);

        /* ps_top_mb_intra_modes */
        u4_size += (i4_num_mbs_in_row + 1) * i4_num_mbs_in_col * sizeof(mb_intra_modes_t);

        /* ps_left_mb_intra_modes */
        u4_size += MAX_PROCESS_CTXT * sizeof(mb_intra_modes_t);
    }

    return u4_size;
}

/**
*******************************************************************************
*
* @brief
*  Function to initialize svc nbr info buffers
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @param[in] ps_mem_rec
*  Pointer to memory allocated for input buffers
*
*******************************************************************************
*/
void isvce_svc_nbr_info_buf_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec)
{
    WORD32 i, j;

    DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;
    UWORD32 u4_wd = ps_codec->s_cfg.u4_wd;
    UWORD32 u4_ht = ps_codec->s_cfg.u4_ht;

    UWORD8 *pu1_buf = ps_mem_rec->pv_base;
    WORD64 i8_alloc_mem_size =
        isvce_get_svc_nbr_info_buf_size(u1_num_spatial_layers, d_spatial_res_ratio, u4_wd, u4_ht);

    ASSERT(1 == MAX_CTXT_SETS);

    for(i = 0; i < MAX_PROCESS_CTXT; i++)
    {
        ps_codec->as_process[i].s_nbr_info_base.ps_layer_nbr_info = (nbr_info_t *) pu1_buf;
        pu1_buf += u1_num_spatial_layers *
                   sizeof(ps_codec->as_process[i].s_nbr_info_base.ps_layer_nbr_info[0]);
        i8_alloc_mem_size -= u1_num_spatial_layers *
                             sizeof(ps_codec->as_process[i].s_nbr_info_base.ps_layer_nbr_info[0]);

        for(j = u1_num_spatial_layers - 1; j >= 0; j--)
        {
            ps_codec->as_process[i].s_nbr_info_base.ps_layer_nbr_info[j].ps_left_mb_info =
                (isvce_mb_info_t *) pu1_buf;
            ps_codec->as_process[i].s_nbr_info.ps_left_mb_info = (isvce_mb_info_t *) pu1_buf;
            pu1_buf += sizeof(ps_codec->as_process[i].s_nbr_info.ps_left_mb_info[0]);
            i8_alloc_mem_size -= sizeof(ps_codec->as_process[i].s_nbr_info.ps_left_mb_info[0]);

            ps_codec->as_process[i].s_nbr_info_base.ps_layer_nbr_info[j].ps_left_mb_intra_modes =
                (mb_intra_modes_t *) pu1_buf;
            ps_codec->as_process[i].s_nbr_info.ps_left_mb_intra_modes =
                (mb_intra_modes_t *) pu1_buf;
            pu1_buf += sizeof(ps_codec->as_process[i].s_nbr_info.ps_left_mb_intra_modes[0]);
            i8_alloc_mem_size -=
                sizeof(ps_codec->as_process[i].s_nbr_info.ps_left_mb_intra_modes[0]);
        }

        ASSERT(i8_alloc_mem_size >= 0);
    }

    for(i = u1_num_spatial_layers - 1; i >= 0; i--)
    {
        isvce_mb_info_t *ps_top_mb_info;
        mb_intra_modes_t *ps_top_intra_modes;

        WORD32 i4_layer_luma_wd =
            ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
        WORD32 i4_layer_luma_ht =
            ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
        WORD32 i4_num_mbs_in_row = i4_layer_luma_wd / MB_SIZE;
        WORD32 i4_num_mbs_in_col = i4_layer_luma_ht / MB_SIZE;

        ps_top_mb_info = (isvce_mb_info_t *) pu1_buf;
        pu1_buf += (i4_num_mbs_in_row + 1) * i4_num_mbs_in_col * sizeof(ps_top_mb_info[0]);
        i8_alloc_mem_size -=
            (i4_num_mbs_in_row + 1) * i4_num_mbs_in_col * sizeof(ps_top_mb_info[0]);

        ps_top_intra_modes = (mb_intra_modes_t *) pu1_buf;
        pu1_buf += (i4_num_mbs_in_row + 1) * i4_num_mbs_in_col * sizeof(ps_top_intra_modes[0]);
        i8_alloc_mem_size -=
            (i4_num_mbs_in_row + 1) * i4_num_mbs_in_col * sizeof(ps_top_intra_modes[0]);

        for(j = 0; j < MAX_PROCESS_CTXT; j++)
        {
            ps_codec->as_process[j].s_nbr_info_base.ps_layer_nbr_info[i].ps_top_row_mb_info =
                ps_top_mb_info;
            ps_codec->as_process[j].s_nbr_info.ps_top_row_mb_info = NULL;

            ps_codec->as_process[j].s_nbr_info_base.ps_layer_nbr_info[i].ps_top_mb_intra_modes =
                ps_top_intra_modes;
            ps_codec->as_process[j].s_nbr_info.ps_top_mb_intra_modes = NULL;
        }

        ASSERT(i8_alloc_mem_size >= 0);
    }
}

/**
*******************************************************************************
*
* @brief
*  isvce_codec_t and proc_t initialisations for an Access Unit
*
* @par Description:
*  Before beginning to encode the frame, the current function initializes all
*  the ctxts (proc, entropy, me, ...) basing on the input configured params.
*  It locates space for storing recon in the encoder picture buffer set, fetches
*  reference frame from encoder picture buffer set. Calls RC pre-enc to get
*  qp and pic type for the current frame. Queues proc jobs so that
*  the other threads can begin encoding. In brief, this function sets up the
*  tone for the entire encoder.
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @param[in] ps_inp_buf
*  Pointer to input buffer context
*
* @returns  error_status
*
* @remarks
*
*
*******************************************************************************
*/
IH264E_ERROR_T isvce_svc_au_init(isvce_codec_t *ps_codec, isvce_inp_buf_t *ps_inp_buf)
{
    svc_au_buf_t *ps_cur_pic;

    WORD32 cur_mv_bank_buf_id;
    WORD32 cur_pic_buf_id;
    WORD32 ref_set_id;
    WORD32 i, j;

    svc_au_data_t *ps_mv_buf = NULL;
    svc_au_buf_t *aps_ref_pic[MAX_REF_PIC_CNT] = {NULL, NULL};
    svc_au_data_t *aps_mv_buf[MAX_REF_PIC_CNT] = {NULL, NULL};

    IH264E_ERROR_T error_status = IH264E_SUCCESS;
    PIC_TYPE_T *pic_type = &ps_codec->pic_type;

    UWORD32 u4_timestamp_high = ps_inp_buf->s_inp_props.u4_timestamp_high;
    UWORD32 u4_timestamp_low = ps_inp_buf->s_inp_props.u4_timestamp_low;
    WORD32 ctxt_sel = ps_codec->i4_encode_api_call_cnt % MAX_CTXT_SETS;
    /* Diamond search Iteration Max Cnt */
    UWORD32 u4_num_layers =
        (ps_codec->s_cfg.u4_enc_speed_preset == IVE_FASTEST) ? (NUM_LAYERS >> 2) : NUM_LAYERS;
    UWORD32 u4_enable_fast_sad = ps_codec->s_cfg.u4_enable_fast_sad;

    if((PIC_I == *pic_type) || (PIC_IDR == *pic_type))
    {
        ps_codec->i4_slice_type = ISLICE;
    }
    else if(PIC_P == *pic_type)
    {
        ps_codec->i4_slice_type = PSLICE;
    }
    else if(PIC_B == *pic_type)
    {
        ps_codec->i4_slice_type = BSLICE;
    }

    ps_codec->u4_is_curr_frm_ref = 0;
    ps_codec->u4_is_curr_frm_ref = (*pic_type != PIC_B);

    if(ps_codec->s_cfg.u4_enable_alt_ref && (*pic_type == PIC_P) &&
       (ps_codec->i4_pic_cnt % (ps_codec->s_cfg.u4_enable_alt_ref + 1)))
    {
        ps_codec->u4_is_curr_frm_ref = 0;
    }

    ps_codec->u4_is_idr = 0;

    if(PIC_IDR == *pic_type)
    {
        ps_codec->u4_is_idr = 1;

        ps_codec->i4_frame_num = 0;

        ps_codec->i4_idr_pic_id++;
    }

    ps_codec->u4_disable_deblock_level = 1;

    if(ps_codec->s_cfg.u4_disable_deblock_level == DISABLE_DEBLK_LEVEL_0)
    {
        ps_codec->u4_disable_deblock_level = 0;
    }
    else if(ps_codec->s_cfg.u4_disable_deblock_level == DISABLE_DEBLK_LEVEL_2)
    {
        if(ps_codec->u4_disable_deblock_level_cnt == DISABLE_DEBLOCK_INTERVAL ||
           ps_codec->i4_slice_type == ISLICE)
        {
            ps_codec->u4_disable_deblock_level = 0;
        }
    }
    else if(ps_codec->s_cfg.u4_disable_deblock_level == DISABLE_DEBLK_LEVEL_3)
    {
        if(ps_codec->i4_slice_type == ISLICE)
        {
            ps_codec->u4_disable_deblock_level = 0;
        }
    }

    if(ps_codec->u4_disable_deblock_level)
    {
        ps_codec->u4_disable_deblock_level_cnt++;
    }
    else
    {
        ps_codec->u4_disable_deblock_level_cnt = 0;
    }

    if(ps_codec->u4_disable_deblock_level == 0)
    {
        if(ps_codec->s_cfg.e_slice_mode != IVE_SLICE_MODE_NONE)
        {
            ps_codec->i4_error_code = IH264E_SLICE_TYPE_INPUT_INVALID;

            return IH264E_SLICE_TYPE_INPUT_INVALID;
        }
    }

    ps_codec->i4_error_code = IH264E_SUCCESS;

    if(ps_codec->i4_gen_header)
    {
        sps_t *ps_sps = NULL;
        pps_t *ps_pps = NULL;
        subset_sps_t *ps_subset_sps = NULL;
        UWORD8 u1_profile_idc = IH264_PROFILE_BASELINE;

        if(ps_codec->as_process[ctxt_sel * MAX_PROCESS_THREADS].u1_spatial_layer_id > 0)
        {
            u1_profile_idc = IH264_SCALABLE_BASELINE;
        }

        ps_sps = ps_codec->ps_sps_base;
        isvce_populate_sps(ps_codec, ps_sps, 0, u1_profile_idc, ps_inp_buf, 0);

        ps_pps = ps_codec->ps_pps_base;
        isvce_populate_pps(ps_codec, ps_pps, 0, 0, 0);

        for(i = 1; i < ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers; i++)
        {
            ps_subset_sps = ps_codec->ps_subset_sps_base + i;
            isvce_populate_subset_sps(ps_codec, ps_subset_sps, i, ps_inp_buf, i);

            /* populate pps header */
            ps_pps = ps_codec->ps_pps_base + i;
            isvce_populate_pps(ps_codec, ps_pps, i, i, i);
        }
    }

    if(IH264E_SUCCESS !=
       isvce_ref_list_refresh(ps_codec, aps_ref_pic, aps_mv_buf, &ref_set_id, pic_type[0]))
    {
        ps_codec->i4_error_code = IH264E_NO_FREE_PICBUF;

        return IH264E_NO_FREE_PICBUF;
    }

    {
        ps_mv_buf = (svc_au_data_t *) ih264_buf_mgr_get_next_free(
            (buf_mgr_t *) ps_codec->pv_svc_au_data_store_mgr, &cur_mv_bank_buf_id);

        if(NULL == ps_mv_buf)
        {
            ps_codec->i4_error_code = IH264E_NO_FREE_MVBANK;

            return IH264E_NO_FREE_MVBANK;
        }

        if(ps_codec->u4_is_curr_frm_ref)
        {
            ih264_buf_mgr_set_status(ps_codec->pv_svc_au_data_store_mgr, cur_mv_bank_buf_id,
                                     BUF_MGR_REF);
        }

        ps_mv_buf->i4_abs_poc = ps_codec->i4_abs_pic_order_cnt;
        ps_mv_buf->i4_buf_id = cur_mv_bank_buf_id;
    }

    {
        ps_cur_pic = (svc_au_buf_t *) ih264_buf_mgr_get_next_free(
            (buf_mgr_t *) ps_codec->pv_ref_buf_mgr, &cur_pic_buf_id);

        if(NULL == ps_cur_pic)
        {
            ps_codec->i4_error_code = IH264E_NO_FREE_PICBUF;

            return IH264E_NO_FREE_PICBUF;
        }

        if(ps_codec->u4_is_curr_frm_ref)
        {
            ih264_buf_mgr_set_status(ps_codec->pv_ref_buf_mgr, cur_pic_buf_id, BUF_MGR_REF);
        }

        if(1 == ps_codec->s_cfg.u4_enable_recon)
        {
            ih264_buf_mgr_set_status(ps_codec->pv_ref_buf_mgr, cur_pic_buf_id, BUF_MGR_IO);
        }

        ps_cur_pic->u4_timestamp_high = ps_inp_buf->s_inp_props.u4_timestamp_high;
        ps_cur_pic->u4_timestamp_low = ps_inp_buf->s_inp_props.u4_timestamp_low;

        ps_cur_pic->i4_abs_poc = ps_codec->i4_poc;
        ps_cur_pic->i4_poc_lsb = ps_codec->i4_pic_order_cnt_lsb;
        ps_cur_pic->i4_frame_num = ps_codec->i4_frame_num;

        ps_cur_pic->i4_buf_id = cur_pic_buf_id;

        ps_cur_pic->i1_temporal_id = isvce_svc_temporal_id_compute(
            ps_codec->i4_poc, ps_codec->s_cfg.s_svc_params.u1_num_temporal_layers, pic_type[0]);
    }

    /*
     * Add the current picture to ref list independent of the fact that it is used
     * as reference or not. This is because, now recon is not in sync with output
     * hence we may need the current recon after some delay. By adding it to ref
     * list we can retrieve the recon any time we want. The information that it is
     * used for ref can still be found by checking the buffer status of pic buf.
     */
    ps_codec->as_ref_set[ref_set_id].i4_pic_cnt = ps_codec->i4_pic_cnt;
    ps_codec->as_ref_set[ref_set_id].i4_poc = ps_codec->i4_poc;
    ps_codec->as_ref_set[ref_set_id].ps_svc_au_data = ps_mv_buf;
    ps_codec->as_ref_set[ref_set_id].ps_pic_buf = ps_cur_pic;

    ps_codec->s_svc_ilp_data.ps_svc_au_data = ps_mv_buf;

    {
        isvce_process_ctxt_t *ps_proc = NULL;

        j = ctxt_sel * MAX_PROCESS_THREADS;

        for(i = j; i < (j + MAX_PROCESS_THREADS); i++)
        {
            ps_proc = &ps_codec->as_process[i];

            ps_proc->s_svc_params = ps_codec->s_cfg.s_svc_params;

            ps_proc->i4_frame_num = ps_codec->i4_frame_num;
            ps_proc->u4_is_idr = ps_codec->u4_is_idr;
            ps_proc->u4_idr_pic_id = ps_codec->i4_idr_pic_id;
            ps_proc->i4_slice_type = ps_codec->i4_slice_type;

            ps_proc->u4_half_x_offset = 0;
            ps_proc->u4_half_y_offset = 0;
            ps_proc->u4_half_xy_offset = 0;

            ps_proc->u4_disable_deblock_level = ps_codec->u4_disable_deblock_level;

            ps_proc->i4_cur_mv_bank_buf_id = cur_mv_bank_buf_id;
            ps_proc->ps_cur_pic = ps_cur_pic;
            ps_proc->ps_cur_mv_buf = ps_mv_buf;

            /*
             * pointer to ref picture
             * 0    : Temporal back reference
             * 1    : Temporal forward reference
             */
            ps_proc->aps_ref_pic[L0] = aps_ref_pic[L0];
            ps_proc->aps_ref_pic[L1] = aps_ref_pic[L1];
            if(ps_codec->pic_type == PIC_B)
            {
                ps_proc->aps_mv_buf[L0] = aps_mv_buf[L0];
                ps_proc->aps_mv_buf[L1] = aps_mv_buf[L1];
            }
            else
            {
                /*
                 * Else is dummy since for non B pic we does not need this
                 * But an assignment here will help in not having a segfault
                 * when we calcualte colpic in P slices
                 */
                ps_proc->aps_mv_buf[L0] = ps_mv_buf;
                ps_proc->aps_mv_buf[L1] = ps_mv_buf;
            }

            ps_proc->s_inp_buf = ps_inp_buf[0];

            ps_proc->i4_encode_api_call_cnt = ps_codec->i4_encode_api_call_cnt;

            ps_proc->i4_pic_cnt = ps_codec->i4_pic_cnt;

            ps_proc->i4_error_code = 0;

            {
                isvce_entropy_ctxt_t *ps_entropy = &ps_proc->s_entropy;

                ps_entropy->i4_sof = 0;
                ps_entropy->i4_eof = 0;
                ps_entropy->ps_sps_base = ps_codec->ps_sps_base;
                ps_entropy->ps_pps_base = ps_codec->ps_pps_base;
                ps_entropy->pu1_slice_idx = ps_proc->pu1_slice_idx;
                ps_entropy->ps_svc_nalu_ext_base = ps_proc->ps_svc_nalu_ext_base;
                ps_entropy->ps_subset_sps_base = ps_proc->ps_subset_sps_base;
                ps_entropy->ps_slice_hdr_base = ps_proc->ps_slice_hdr_base;
                ps_entropy->ps_svc_slice_hdr_base = ps_proc->ps_svc_slice_hdr_base;
                ps_entropy->i4_abs_pic_order_cnt = ps_codec->i4_poc;

                ps_entropy->i1_transform_8x8_mode_flag = 0;

                ps_entropy->i4_error_code = IH264E_SUCCESS;
                ps_proc->s_entropy.u4_is_last = ps_inp_buf->s_inp_props.u4_is_last;
                ps_proc->s_entropy.i4_pic_cnt = ps_codec->i4_pic_cnt;

                ps_entropy->u4_timestamp_low = u4_timestamp_low;
                ps_entropy->u4_timestamp_high = u4_timestamp_high;
            }

            {
                isvce_me_ctxt_t *ps_me_ctxt = &ps_proc->s_me_ctxt;

                ps_me_ctxt->ai2_srch_boundaries[0] = ps_codec->s_cfg.u4_srch_rng_x;
                ps_me_ctxt->ai2_srch_boundaries[1] = ps_codec->s_cfg.u4_srch_rng_y;

                ps_me_ctxt->u4_half_x_offset = ps_proc->u4_half_x_offset;
                ps_me_ctxt->u4_half_y_offset = ps_proc->u4_half_y_offset;
                ps_me_ctxt->u4_half_xy_offset = ps_proc->u4_half_xy_offset;

                ps_me_ctxt->u4_enable_fast_sad = u4_enable_fast_sad;
                ps_me_ctxt->u4_enable_hpel = ps_codec->s_cfg.u4_enable_hpel;
                ps_me_ctxt->u4_num_layers = u4_num_layers;
                ps_me_ctxt->u4_me_speed_preset = ps_codec->s_cfg.u4_me_speed_preset;

                if((i == j) && (0 == ps_codec->i4_poc))
                {
                    isvce_init_mv_bits(ps_me_ctxt);
                }
            }

            ps_proc->ps_ngbr_avbl = &(ps_proc->s_ngbr_avbl);
        }
    }

    return error_status;
}

void isvce_init_quant_params(isvce_process_ctxt_t *ps_proc, WORD32 qp)
{
    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    /* quant params */
    quant_params_t *ps_qp_params;

    /* ptr to forward quant threshold matrix */
    const UWORD16 *pu2_thres_mat = NULL;

    /* ptr to forward scale matrix */
    const UWORD16 *pu2_scale_mat = gu2_quant_scale_matrix_4x4;

    /* ptr to inverse scale matrix */
    const UWORD16 *pu2_iscale_mat = gau2_ih264_iquant_scale_matrix_4x4;

    /* temp var */
    UWORD32 u4_qp[3], u4_qp_div6, u4_qp_mod6;
    COMPONENT_TYPE plane;
    WORD32 i;
    UWORD32 u4_satdq_t;
    const UWORD16 *pu2_smat;

    /********************************************************************/
    /* init quant params for all planes Y, U and V                      */
    /********************************************************************/
    /* luma qp */
    u4_qp[Y] = qp;

    /* chroma qp
     * TODO_LATER : just in case if the chroma planes use different qp's this
     * needs to be corrected accordingly.
     */
    u4_qp[U] = gu1_qpc_fqpi[qp];
    u4_qp[V] = gu1_qpc_fqpi[qp];

    plane = Y;
    while(plane <= V)
    {
        u4_qp_div6 = (u4_qp[plane] / 6);
        u4_qp_mod6 = (u4_qp[plane] % 6);

        ps_qp_params = ps_proc->ps_qp_params[plane];

        /* mb qp */
        ps_qp_params->u1_mb_qp = u4_qp[plane];

        /* mb qp / 6 */
        ps_qp_params->u1_qp_div = u4_qp_div6;

        /* mb qp % 6 */
        ps_qp_params->u1_qp_rem = u4_qp_mod6;

        /* QP bits */
        ps_qp_params->u1_qbits = QP_BITS_h264_4x4 + u4_qp_div6;

        /* forward scale matrix */
        ps_qp_params->pu2_scale_mat = pu2_scale_mat + (u4_qp_mod6 * 16);

        /* threshold matrix & weight for quantization */
        pu2_thres_mat = gu2_forward_quant_threshold_4x4 + (u4_qp_mod6 * 16);
        for(i = 0; i < 16; i++)
        {
            ps_qp_params->pu2_thres_mat[i] = pu2_thres_mat[i] >> (8 - u4_qp_div6);
            ps_qp_params->pu2_weigh_mat[i] = 16;
        }

        /* qp dependent rounding constant */
        ps_qp_params->u4_dead_zone = gu4_forward_quant_round_factor_4x4[u4_qp_div6];

        /* slice dependent rounding constant */
        if(ps_proc->i4_slice_type != ISLICE && ps_proc->i4_slice_type != SISLICE)
        {
            ps_qp_params->u4_dead_zone >>= 1;
        }

        /* SATQD threshold for zero block prediction */
        if(ps_codec->s_cfg.u4_enable_satqd)
        {
            pu2_smat = ps_qp_params->pu2_scale_mat;

            u4_satdq_t = ((1 << (ps_qp_params->u1_qbits)) - ps_qp_params->u4_dead_zone);

            ps_qp_params->pu2_sad_thrsh[0] = u4_satdq_t / MAX(pu2_smat[3], pu2_smat[11]);
            ps_qp_params->pu2_sad_thrsh[1] = u4_satdq_t / MAX(pu2_smat[1], pu2_smat[9]);
            ps_qp_params->pu2_sad_thrsh[2] = u4_satdq_t / pu2_smat[15];
            ps_qp_params->pu2_sad_thrsh[3] = u4_satdq_t / pu2_smat[7];
            ps_qp_params->pu2_sad_thrsh[4] = u4_satdq_t / MAX(pu2_smat[12], pu2_smat[14]);
            ps_qp_params->pu2_sad_thrsh[5] = u4_satdq_t / MAX(pu2_smat[4], pu2_smat[6]);
            ps_qp_params->pu2_sad_thrsh[6] = u4_satdq_t / pu2_smat[13];
            ps_qp_params->pu2_sad_thrsh[7] = u4_satdq_t / pu2_smat[5];
            ps_qp_params->pu2_sad_thrsh[8] =
                u4_satdq_t / MAX(MAX3(pu2_smat[0], pu2_smat[2], pu2_smat[8]), pu2_smat[10]);
        }

        /* inverse scale matrix */
        ps_qp_params->pu2_iscale_mat = pu2_iscale_mat + (u4_qp_mod6 * 16);

        plane += 1;
    }
}

/**
*******************************************************************************
*
* @brief
*  isvce_codec_t and proc_t initialisations for an Access Unit
*
* @par Description:
*  Before beginning to encode the frame, the current function initializes all
*  the ctxts (proc, entropy, me, ...) basing on the input configured params.
*  It locates space for storing recon in the encoder picture buffer set, fetches
*  reference frame from encoder picture buffer set. Calls RC pre-enc to get
*  qp and pic type for the current frame. Queues proc jobs so that
*  the other threads can begin encoding. In brief, this function sets up the
*  tone for the entire encoder.
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @param[in] ps_inp_buf
*  Pointer to input buffer context
*
* @param[in] u1_spatial_layer_id
*  Spatial Layer IDl 0 => Base layer
*
* @returns  error_status
*
* @remarks
*
*
*******************************************************************************
*/
IH264E_ERROR_T isvce_svc_layer_pic_init(isvce_codec_t *ps_codec, isvce_inp_buf_t *ps_inp_buf,
                                        UWORD8 u1_spatial_layer_id)
{
    WORD32 i;

    IH264E_ERROR_T error_status = IH264E_SUCCESS;
    IH264_ERROR_T ret = IH264_SUCCESS;
    PIC_TYPE_T e_pic_type = ps_codec->pic_type;

    ASSERT(MAX_CTXT_SETS == 1);

    for(i = 0; i < MAX_PROCESS_THREADS; i++)
    {
        isvce_process_ctxt_t *ps_proc = &ps_codec->as_process[i];
        isvce_entropy_ctxt_t *ps_entropy = &ps_proc->s_entropy;
        isvce_deblk_ctxt_t *ps_deblk = &ps_proc->s_deblk_ctxt;
        isvce_me_ctxt_t *ps_me_ctxt = &ps_proc->s_me_ctxt;
        svc_au_buf_t *ps_cur_pic = ps_proc->ps_cur_pic;
        svc_au_buf_t *aps_ref_pic[MAX_REF_PIC_CNT] = {ps_proc->aps_ref_pic[L0],
                                                      ps_proc->aps_ref_pic[L1]};

        ps_proc->u1_spatial_layer_id = u1_spatial_layer_id;

        ps_proc->s_src_pic_buf_props = ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id];

        ps_proc->s_rec_pic_buf_props = ps_cur_pic->ps_layer_yuv_buf_props[u1_spatial_layer_id];

        ASSERT(0 == (ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id].u4_width % MB_SIZE));
        ASSERT(0 == (ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id].u4_height % MB_SIZE));

        ps_proc->i4_wd_mbs =
            ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id].u4_width / MB_SIZE;
        ps_proc->i4_ht_mbs =
            ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id].u4_height / MB_SIZE;

        ps_proc->u1_frame_qp = ps_codec->au4_frame_qp[u1_spatial_layer_id];

        ps_proc->u1_mb_qp = ps_proc->u1_frame_qp;
        ps_entropy->ps_mb_qp_ctxt->u1_cur_mb_qp = ps_proc->u1_frame_qp;

        isvce_init_quant_params(ps_proc, ps_proc->u1_frame_qp);

        memset(&ps_proc->s_frame_info, 0, sizeof(frame_info_t));

        /* row '-1' */
        memset(ps_proc->pu1_proc_map - ps_proc->i4_wd_mbs, 1,
               ps_proc->i4_wd_mbs * sizeof(ps_proc->pu1_proc_map[0]));

        /* row 0 to ht in mbs */
        memset(ps_proc->pu1_proc_map, 0,
               ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs * sizeof(ps_proc->pu1_proc_map[0]));

        /* row '-1' */
        memset(ps_proc->pu1_deblk_map - ps_proc->i4_wd_mbs, 1,
               ps_proc->i4_wd_mbs * sizeof(ps_proc->pu1_deblk_map[0]));

        /* row 0 to ht in mbs */
        memset(ps_proc->pu1_deblk_map, 0,
               ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs * sizeof(ps_proc->pu1_deblk_map[0]));

        /* row '-1' */
        memset(ps_proc->pu1_me_map - ps_proc->i4_wd_mbs, 1,
               ps_proc->i4_wd_mbs * sizeof(ps_proc->pu1_me_map[0]));

        /* row 0 to ht in mbs */
        memset(ps_proc->pu1_me_map, 0,
               ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs * sizeof(ps_proc->pu1_me_map[0]));

        if(IVE_AIR_MODE_NONE != ps_codec->s_cfg.e_air_mode)
        {
            ps_codec->i4_air_pic_cnt =
                (ps_codec->i4_air_pic_cnt + 1) % ps_codec->s_cfg.u4_air_refresh_period;

            if(!ps_codec->i4_air_pic_cnt)
            {
                memset(ps_proc->pu1_is_intra_coded, 0,
                       ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs *
                           sizeof(ps_proc->pu1_is_intra_coded[0]));
            }
        }

        if(ps_codec->s_cfg.e_slice_mode == IVE_SLICE_MODE_NONE)
        {
            memset(ps_proc->pu1_slice_idx, 0,
                   ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs * sizeof(ps_proc->pu1_slice_idx[0]));
        }
        else if(ps_codec->s_cfg.e_slice_mode == IVE_SLICE_MODE_BLOCKS)
        {
            UWORD8 *pu1_slice_idx = ps_proc->pu1_slice_idx;
            WORD32 i4_mb_y = 0, slice_idx = 0, cnt;

            while(i4_mb_y < ps_proc->i4_ht_mbs)
            {
                if(i4_mb_y + (WORD32) ps_codec->s_cfg.u4_slice_param < ps_proc->i4_ht_mbs)
                {
                    cnt = ps_codec->s_cfg.u4_slice_param * ps_proc->i4_wd_mbs;
                    i4_mb_y += ps_codec->s_cfg.u4_slice_param;
                }
                else
                {
                    cnt = (ps_proc->i4_ht_mbs - i4_mb_y) * ps_proc->i4_wd_mbs;
                    i4_mb_y += (ps_proc->i4_ht_mbs - i4_mb_y);
                }

                memset(pu1_slice_idx, slice_idx, cnt * sizeof(pu1_slice_idx[0]));

                slice_idx++;
                pu1_slice_idx += cnt;
            }
        }

        if((e_pic_type != PIC_IDR) && (e_pic_type != PIC_I))
        {
            ps_proc->as_ref_pic_buf_props[L0] =
                aps_ref_pic[L0]->ps_layer_yuv_buf_props[u1_spatial_layer_id];
            ps_proc->as_ref_pic_buf_props[L1] =
                aps_ref_pic[L1]->ps_layer_yuv_buf_props[u1_spatial_layer_id];
        }

        ps_entropy->i4_gen_header = ps_codec->i4_gen_header && (0 == u1_spatial_layer_id);
        ps_entropy->i4_gen_subset_sps =
            (ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers > 1) && ps_codec->i4_gen_header;

        /* row '-1' */
        memset(ps_entropy->pu1_entropy_map - ps_proc->i4_wd_mbs, 1,
               ps_proc->i4_wd_mbs * sizeof(ps_entropy->pu1_entropy_map[0]));

        /* row 0 to ht in mbs */
        memset(ps_entropy->pu1_entropy_map, 0,
               ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs * sizeof(ps_entropy->pu1_entropy_map[0]));

        isvce_init_cabac_table(ps_entropy);

        ps_entropy->i4_wd_mbs = ps_proc->i4_wd_mbs;
        ps_entropy->i4_ht_mbs = ps_proc->i4_ht_mbs;

        ps_entropy->u1_entropy_coding_mode_flag =
            ((ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers > 1) && (0 == u1_spatial_layer_id))
                ? CAVLC
                : ps_codec->s_cfg.u4_entropy_coding_mode;

        ps_proc->s_entropy.pi4_mb_skip_run[0] = 0;

        ps_entropy->u4_header_bits[MB_TYPE_INTRA] = 0;
        ps_entropy->u4_header_bits[MB_TYPE_INTER] = 0;
        ps_entropy->u4_residue_bits[MB_TYPE_INTRA] = 0;
        ps_entropy->u4_residue_bits[MB_TYPE_INTER] = 0;

        ps_entropy->u1_spatial_layer_id = ps_proc->u1_spatial_layer_id;

        ps_deblk->pu1_slice_idx = ps_proc->pu1_slice_idx;

        ps_me_ctxt->u1_mb_qp = ps_codec->au4_frame_qp[u1_spatial_layer_id];

        {
            UWORD8 u1_min_qp;
            UWORD8 u1_max_qp;

            svc_sub_pic_rc_ctxt_t *ps_sub_pic_rc_ctxt = ps_proc->ps_sub_pic_rc_ctxt;
            svc_sub_pic_rc_layer_variables_t *ps_layer_variables =
                &ps_sub_pic_rc_ctxt->s_sub_pic_rc_variables.s_layer_variables;

            switch(ps_proc->i4_slice_type)
            {
                case ISLICE:
                {
                    u1_min_qp = ps_codec->s_cfg.au4_i_qp_min[u1_spatial_layer_id];
                    u1_max_qp = ps_codec->s_cfg.au4_i_qp_max[u1_spatial_layer_id];

                    break;
                }
                case PSLICE:
                {
                    u1_min_qp = ps_codec->s_cfg.au4_p_qp_min[u1_spatial_layer_id];
                    u1_max_qp = ps_codec->s_cfg.au4_p_qp_max[u1_spatial_layer_id];

                    break;
                }
                default:
                {
                    u1_min_qp = ps_codec->s_cfg.au4_b_qp_min[u1_spatial_layer_id];
                    u1_max_qp = ps_codec->s_cfg.au4_b_qp_max[u1_spatial_layer_id];

                    break;
                }
            }

            ps_layer_variables->i4_max_num_reference_frames = ps_codec->i4_max_num_reference_frames;
            ps_layer_variables->i4_slice_type = ps_proc->i4_slice_type;
            ps_layer_variables->i4_frame_num = ps_proc->i4_frame_num;
            ps_layer_variables->u1_frame_qp = ps_proc->u1_frame_qp;
            ps_layer_variables->u1_spatial_layer_id = u1_spatial_layer_id;
            ps_layer_variables->u1_min_qp = u1_min_qp;
            ps_layer_variables->u1_max_qp = u1_max_qp;

            isvce_sub_pic_rc_ctxt_layer_init(ps_proc->ps_sub_pic_rc_ctxt);
        }
    }

    {
        job_t s_job;

        s_job.i4_cmd = CMD_PROCESS;
        s_job.i2_mb_cnt =
            ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id].u4_width / MB_SIZE;
        s_job.i2_mb_x = 0;

        for(i = 0; i < (WORD32) (ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id].u4_height /
                                 MB_SIZE);
            i++)
        {
            s_job.i2_mb_y = i;

            ret = ih264_list_queue(ps_codec->pv_proc_jobq, &s_job, 1);

            if(ret != IH264_SUCCESS)
            {
                ps_codec->i4_error_code = ret;

                return IH264E_FAIL;
            }
        }

        /* Once all the jobs are queued, terminate the queue */
        /* Since the threads are created and deleted in each call, terminating
        here is not an issue */
        ih264_list_terminate(ps_codec->pv_proc_jobq);
    }

    ps_codec->i4_gen_header = 0;

    return error_status;
}

/**
*******************************************************************************
*
* @brief   initialize process context.
*
* @par Description:
*  Before dispatching the current job to process thread, the process context
*  associated with the job is initialized. Usually every job aims to encode one
*  row of mb's. Basing on the row indices provided by the job, the process
*  context's buffer ptrs, slice indices and other elements that are necessary
*  during core-coding are initialized.
*
* @param[in] ps_proc
*  Pointer to the current process context
*
* @returns error status
*
* @remarks none
*
*******************************************************************************
*/
IH264E_ERROR_T isvce_init_layer_proc_ctxt(isvce_process_ctxt_t *ps_proc)
{
    WORD32 i4_mb_x, i4_mb_y;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    n_mb_process_ctxt_t *ps_n_mb_ctxt = &ps_proc->s_n_mb_ctxt;
    quant_params_t *ps_qp_params = ps_proc->ps_qp_params[0];
    isvce_deblk_ctxt_t *ps_deblk = &ps_proc->s_deblk_ctxt;
    isvce_bs_ctxt_t *ps_bs = &(ps_deblk->s_bs_ctxt);
    svc_au_data_t *ps_cur_mv_buf = ps_proc->ps_cur_mv_buf;

    i4_mb_x = ps_proc->i4_mb_x;
    i4_mb_y = ps_proc->i4_mb_y;

    ASSERT((ps_codec->s_cfg.u4_wd - ps_codec->s_cfg.u4_disp_wd) == 0);
    ASSERT((ps_codec->s_cfg.u4_ht - ps_codec->s_cfg.u4_disp_ht) == 0);

    ps_proc->i4_nmb_ntrpy = ps_proc->i4_wd_mbs;
    ps_proc->u4_nmb_me = 1;

    ps_proc->s_src_buf_props = ps_proc->s_src_pic_buf_props;
    ps_proc->s_rec_buf_props = ps_proc->s_rec_pic_buf_props;
    ps_proc->as_ref_buf_props[0] = ps_proc->as_ref_pic_buf_props[0];
    ps_proc->as_ref_buf_props[1] = ps_proc->as_ref_pic_buf_props[1];

    ps_proc->s_src_buf_props.as_component_bufs[0].pv_data =
        ((UWORD8 *) ps_proc->s_src_buf_props.as_component_bufs[0].pv_data) + (i4_mb_x * MB_SIZE) +
        ps_proc->s_src_buf_props.as_component_bufs[0].i4_data_stride * (i4_mb_y * MB_SIZE);
    ps_proc->s_src_buf_props.as_component_bufs[1].pv_data =
        ((UWORD8 *) ps_proc->s_src_pic_buf_props.as_component_bufs[1].pv_data) +
        (i4_mb_x * MB_SIZE) +
        ps_proc->s_src_buf_props.as_component_bufs[1].i4_data_stride * (i4_mb_y * BLK8x8SIZE);

    ps_proc->s_rec_buf_props.as_component_bufs[0].pv_data =
        ((UWORD8 *) ps_proc->s_rec_buf_props.as_component_bufs[0].pv_data) + (i4_mb_x * MB_SIZE) +
        ps_proc->s_rec_buf_props.as_component_bufs[0].i4_data_stride * (i4_mb_y * MB_SIZE);
    ps_proc->s_rec_buf_props.as_component_bufs[1].pv_data =
        ((UWORD8 *) ps_proc->s_rec_buf_props.as_component_bufs[1].pv_data) + (i4_mb_x * MB_SIZE) +
        ps_proc->s_rec_buf_props.as_component_bufs[1].i4_data_stride * (i4_mb_y * BLK8x8SIZE);

    ps_proc->as_ref_buf_props[0].as_component_bufs[0].pv_data =
        ((UWORD8 *) ps_proc->as_ref_buf_props[0].as_component_bufs[0].pv_data) +
        (i4_mb_x * MB_SIZE) +
        ps_proc->as_ref_buf_props[0].as_component_bufs[0].i4_data_stride * (i4_mb_y * MB_SIZE);
    ps_proc->as_ref_buf_props[0].as_component_bufs[1].pv_data =
        ((UWORD8 *) ps_proc->as_ref_buf_props[0].as_component_bufs[1].pv_data) +
        (i4_mb_x * MB_SIZE) +
        ps_proc->as_ref_buf_props[0].as_component_bufs[1].i4_data_stride * (i4_mb_y * BLK8x8SIZE);

    ps_proc->as_ref_buf_props[1].as_component_bufs[0].pv_data =
        ((UWORD8 *) ps_proc->as_ref_buf_props[1].as_component_bufs[0].pv_data) +
        (i4_mb_x * MB_SIZE) +
        ps_proc->as_ref_buf_props[1].as_component_bufs[0].i4_data_stride * (i4_mb_y * MB_SIZE);
    ps_proc->as_ref_buf_props[1].as_component_bufs[1].pv_data =
        ((UWORD8 *) ps_proc->as_ref_buf_props[1].as_component_bufs[1].pv_data) +
        (i4_mb_x * MB_SIZE) +
        ps_proc->as_ref_buf_props[1].as_component_bufs[1].i4_data_stride * (i4_mb_y * BLK8x8SIZE);

    ps_proc->pv_mb_coeff_data =
        ((UWORD8 *) ps_proc->pv_pic_mb_coeff_data) + i4_mb_y * ps_codec->u4_size_coeff_data;

    ps_proc->pv_mb_header_data =
        ((UWORD8 *) ps_proc->pv_pic_mb_header_data) + i4_mb_y * ps_codec->u4_size_header_data;

    ps_proc->i4_cur_slice_idx = ps_proc->pu1_slice_idx[i4_mb_y * ps_proc->i4_wd_mbs + i4_mb_x];

    ps_proc->ps_mb_info =
        ps_cur_mv_buf->ps_svc_layer_data[ps_proc->u1_spatial_layer_id].ps_mb_info +
        i4_mb_y * ps_proc->i4_wd_mbs;

    ps_proc->ps_col_mb =
        ps_proc->aps_mv_buf[1]->ps_svc_layer_data[ps_proc->u1_spatial_layer_id].ps_mb_info +
        i4_mb_y * ps_proc->i4_wd_mbs;

    {
        ps_proc->s_nbr_info.ps_top_row_mb_info =
            ps_proc->s_nbr_info_base.ps_layer_nbr_info[ps_proc->u1_spatial_layer_id]
                .ps_top_row_mb_info +
            (i4_mb_x + (i4_mb_y - 1) * ps_proc->i4_wd_mbs);

        ps_proc->s_nbr_info.ps_top_mb_intra_modes =
            ps_proc->s_nbr_info_base.ps_layer_nbr_info[ps_proc->u1_spatial_layer_id]
                .ps_top_mb_intra_modes +
            (i4_mb_x + (i4_mb_y - 1) * ps_proc->i4_wd_mbs);
    }

    ps_proc->pu4_mb_pu_cnt =
        ps_cur_mv_buf->ps_svc_layer_data[ps_proc->u1_spatial_layer_id].pu4_num_pus_in_mb +
        (i4_mb_y * ps_proc->i4_wd_mbs);

    ps_proc->ps_mb_info->u2_mb_type = I16x16;

    ps_proc->u4_lambda = gu1_qp0[ps_qp_params->u1_mb_qp];

    ps_proc->i4_mb_distortion = SHRT_MAX;

    if(i4_mb_x == 0)
    {
        ps_proc->s_nbr_info.ps_left_mb_info[0].i4_mb_distortion = 0;
    }

    ps_proc->i4_mb_cost = INT_MAX;

    ps_deblk->i4_mb_x = ps_proc->i4_mb_x;
    /* deblk lags the current mb proc by 1 row */
    /* NOTE: Intra prediction has to happen with non deblocked samples used as
     * reference */
    /* Hence to deblk MB 0 of row 0, you have wait till MB 0 of row 1 is encoded.
     */
    /* For simplicity, we chose to lag deblking by 1 Row wrt to proc */
    ps_deblk->i4_mb_y = ps_proc->i4_mb_y - 1;

    ps_deblk->s_rec_pic_buf_props = ps_proc->s_rec_pic_buf_props;

    ps_bs->i4_mb_x = ps_proc->i4_mb_x;
    ps_bs->i4_mb_y = ps_proc->i4_mb_y;

    ps_n_mb_ctxt->i4_mb_x = 0;
    ps_n_mb_ctxt->i4_mb_y = ps_deblk->i4_mb_y;
    ps_n_mb_ctxt->i4_n_mbs = ps_proc->i4_nmb_ntrpy;

    return IH264E_SUCCESS;
}

/**
*******************************************************************************
*
* @brief
*  Returns size of buffers for storing SVC ILP data
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
UWORD32 isvce_get_svc_ilp_buf_size(UWORD8 u1_num_spatial_layers, DOUBLE d_spatial_res_ratio,
                                   UWORD32 u4_wd, UWORD32 u4_ht)
{
    WORD32 i;

    UWORD32 u4_size = 0;

    if(u1_num_spatial_layers > 1)
    {
        /* ps_intra_recon_bufs */
        u4_size += u1_num_spatial_layers * sizeof(yuv_buf_props_t);

        /* ps_residual_bufs */
        u4_size += u1_num_spatial_layers * sizeof(yuv_buf_props_t);

        /* aps_layer_resampler_props[Y] */
        u4_size += u1_num_spatial_layers * sizeof(layer_resampler_props_t);

        /* aps_layer_resampler_props[UV] */
        u4_size += u1_num_spatial_layers * sizeof(layer_resampler_props_t);

        for(i = u1_num_spatial_layers - 1; i >= 0; i--)
        {
            WORD32 i4_layer_luma_wd =
                ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
            WORD32 i4_layer_luma_ht =
                ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
            WORD32 i4_layer_luma_samples =
                (ALIGN16(i4_layer_luma_wd) + PAD_WD) * (i4_layer_luma_ht + PAD_HT);
            WORD32 i4_layer_uv_wd = i4_layer_luma_wd;
            WORD32 i4_layer_uv_ht = i4_layer_luma_ht / 2.0 + 0.99;
            WORD32 i4_layer_uv_samples =
                (ALIGN16(i4_layer_uv_wd) + PAD_WD) * (i4_layer_uv_ht + PAD_HT);

            /* ps_intra_recon_bufs */
            u4_size += (i4_layer_luma_samples + i4_layer_uv_samples) * sizeof(UWORD8);

            /* ps_residual_bufs */
            u4_size += (i4_layer_luma_samples + i4_layer_uv_samples) * sizeof(WORD16);
        }
    }
    else
    {
        WORD32 i4_layer_luma_wd = u4_wd;
        WORD32 i4_layer_luma_ht = u4_ht;
        WORD32 i4_layer_luma_samples =
            (ALIGN16(i4_layer_luma_wd) + PAD_WD) * (i4_layer_luma_ht + PAD_HT);
        WORD32 i4_layer_uv_wd = i4_layer_luma_wd;
        WORD32 i4_layer_uv_ht = i4_layer_luma_ht / 2.0 + 0.99;
        WORD32 i4_layer_uv_samples = (ALIGN16(i4_layer_uv_wd) + PAD_WD) * (i4_layer_uv_ht + PAD_HT);

        /* ps_residual_bufs */
        u4_size += sizeof(yuv_buf_props_t);

        /* ps_residual_bufs */
        u4_size += (i4_layer_luma_samples + i4_layer_uv_samples) * sizeof(WORD16);
    }

    return u4_size;
}

static void isvce_layer_resampler_props_init(layer_resampler_props_t *ps_layer_props,
                                             DOUBLE d_spatial_res_ratio, UWORD32 u4_wd,
                                             UWORD32 u4_ht, UWORD8 u1_level_idc,
                                             UWORD8 u1_is_chroma)
{
    const UWORD8 u1_ref_layer_field_pic_flag = 0;
    const UWORD8 u1_field_pic_flag = 0;
    const UWORD8 u1_frame_mbs_only_flag = 1;
    const UWORD8 u1_ref_layer_frame_mbs_only_flag = 1;
    const UWORD8 u1_bot_field_flag = 0;
    const WORD32 i4_scaled_ref_layer_left_offset = 0;
    const WORD32 i4_scaled_ref_layer_top_offset = 0;
    const WORD32 i4_ref_layer_chroma_phase_x_plus1 = 1;
    const WORD32 i4_ref_layer_chroma_phase_y_plus1 = 1;
    const WORD32 i4_chroma_phase_x_plus1 = 1;
    const WORD32 i4_chroma_phase_y_plus1 = 1;
    const WORD32 i4_sub_wd_chroma = 2;
    const WORD32 i4_sub_ht_chroma = 2;
    UWORD32 u4_ref_wd = (u4_wd / d_spatial_res_ratio);
    UWORD32 u4_ref_ht = (u4_ht / d_spatial_res_ratio) * (1 + u1_ref_layer_field_pic_flag);
    UWORD32 u4_scaled_wd = u4_wd;
    UWORD32 u4_scaled_ht = u4_ht * (1 + u1_field_pic_flag);

    u4_ref_wd = u4_ref_wd >> u1_is_chroma;
    u4_ref_ht = u4_ref_ht >> u1_is_chroma;
    u4_scaled_wd = u4_scaled_wd >> u1_is_chroma;
    u4_scaled_ht = u4_scaled_ht >> u1_is_chroma;

    if(u1_is_chroma)
    {
        ps_layer_props->i4_refphase_x = i4_ref_layer_chroma_phase_x_plus1 - 1;
        ps_layer_props->i4_refphase_y = i4_ref_layer_chroma_phase_y_plus1 - 1;
        ps_layer_props->i4_phase_x = i4_chroma_phase_x_plus1 - 1;
        ps_layer_props->i4_phase_y = i4_chroma_phase_y_plus1 - 1;
        ps_layer_props->u4_sub_wd = i4_sub_wd_chroma;
        ps_layer_props->u4_sub_ht = i4_sub_ht_chroma;
        ps_layer_props->u4_mb_wd = MB_SIZE >> 1;
        ps_layer_props->u4_mb_ht = MB_SIZE >> 1;
    }
    else
    {
        ps_layer_props->i4_refphase_x = 0;
        ps_layer_props->i4_refphase_y = 0;
        ps_layer_props->i4_phase_x = 0;
        ps_layer_props->i4_phase_y = 0;
        ps_layer_props->u4_sub_wd = 1;
        ps_layer_props->u4_sub_ht = 1;
        ps_layer_props->u4_mb_wd = MB_SIZE;
        ps_layer_props->u4_mb_ht = MB_SIZE;
    }

    if(u1_level_idc <= 30)
    {
        ps_layer_props->u4_shift_x = 16;
        ps_layer_props->u4_shift_y = 16;
    }
    else
    {
        ps_layer_props->u4_shift_x = 31 - isvcd_get_ceil_log2(u4_ref_wd);
        ps_layer_props->u4_shift_y = 31 - isvcd_get_ceil_log2(u4_ref_ht);
    }

    if((0 == u1_frame_mbs_only_flag) || (0 == u1_ref_layer_frame_mbs_only_flag))
    {
        ps_layer_props->i4_phase_y = ps_layer_props->i4_phase_y + 4 * u1_bot_field_flag;

        if(1 == u1_ref_layer_frame_mbs_only_flag)
        {
            ps_layer_props->i4_refphase_y = (2 * ps_layer_props->i4_refphase_y) + 2;
        }
        else
        {
            ps_layer_props->i4_refphase_y = ps_layer_props->i4_refphase_y + (4 * u1_bot_field_flag);
        }
    }

    ps_layer_props->u4_scale_x =
        ((u4_ref_wd << ps_layer_props->u4_shift_x) + (u4_scaled_wd >> 1)) / (u4_scaled_wd);
    ps_layer_props->u4_scale_y =
        ((u4_ref_ht << ps_layer_props->u4_shift_y) + (u4_scaled_ht >> 1)) / (u4_scaled_ht);

    ps_layer_props->i4_offset_x = i4_scaled_ref_layer_left_offset / ps_layer_props->u4_sub_wd;
    ps_layer_props->i4_add_x =
        (((u4_ref_wd * (2 + ps_layer_props->i4_phase_x)) << (ps_layer_props->u4_shift_x - 2)) +
         (u4_scaled_wd >> 1)) /
            u4_scaled_wd +
        (1 << (ps_layer_props->u4_shift_x - 5));
    ps_layer_props->i4_delta_x = 4 * (2 + ps_layer_props->i4_refphase_x);

    if((1 == u1_frame_mbs_only_flag) && (1 == u1_ref_layer_frame_mbs_only_flag))
    {
        ps_layer_props->i4_offset_y = i4_scaled_ref_layer_top_offset / ps_layer_props->u4_sub_ht;
        ps_layer_props->i4_add_y =
            (((u4_ref_ht * (2 + ps_layer_props->i4_phase_y)) << (ps_layer_props->u4_shift_y - 2)) +
             (u4_scaled_ht >> 1)) /
                u4_scaled_ht +
            (1 << (ps_layer_props->u4_shift_y - 5));
        ps_layer_props->i4_delta_y = 4 * (2 + ps_layer_props->i4_refphase_y);
    }
    else
    {
        ps_layer_props->i4_offset_y =
            i4_scaled_ref_layer_top_offset / (2 * ps_layer_props->u4_sub_ht);
        ps_layer_props->i4_add_y =
            (((u4_ref_ht * (2 + ps_layer_props->i4_phase_y)) << (ps_layer_props->u4_shift_y - 3)) +
             (u4_scaled_ht >> 1)) /
                u4_scaled_ht +
            (1 << (ps_layer_props->u4_shift_y - 5));
        ps_layer_props->i4_delta_y = 2 * (2 + ps_layer_props->i4_refphase_y);
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
void isvce_svc_ilp_buf_init(isvce_codec_t *ps_codec, iv_mem_rec_t *ps_mem_rec)
{
    UWORD8 u1_num_spatial_layers = ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers;
    DOUBLE d_spatial_res_ratio = ps_codec->s_cfg.s_svc_params.d_spatial_res_ratio;
    UWORD32 u4_wd = ps_codec->s_cfg.u4_wd;
    UWORD32 u4_ht = ps_codec->s_cfg.u4_ht;
    UWORD8 *pu1_buf = ps_mem_rec->pv_base;
    WORD64 i8_alloc_mem_size =
        isvce_get_svc_ilp_buf_size(u1_num_spatial_layers, d_spatial_res_ratio, u4_wd, u4_ht);

    if(u1_num_spatial_layers > 1)
    {
        WORD32 i, j;

        ps_codec->s_svc_ilp_data.ps_intra_recon_bufs = (yuv_buf_props_t *) pu1_buf;
        pu1_buf += u1_num_spatial_layers * sizeof(ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[0]);
        i8_alloc_mem_size -=
            u1_num_spatial_layers * sizeof(ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[0]);

        ps_codec->s_svc_ilp_data.ps_residual_bufs = (yuv_buf_props_t *) pu1_buf;
        pu1_buf += u1_num_spatial_layers * sizeof(ps_codec->s_svc_ilp_data.ps_residual_bufs[0]);
        i8_alloc_mem_size -=
            u1_num_spatial_layers * sizeof(ps_codec->s_svc_ilp_data.ps_residual_bufs[0]);

        for(i = 0; i < NUM_SP_COMPONENTS; i++)
        {
            ps_codec->s_svc_ilp_data.aps_layer_resampler_props[i] =
                (layer_resampler_props_t *) pu1_buf;
            pu1_buf += u1_num_spatial_layers *
                       sizeof(ps_codec->s_svc_ilp_data.aps_layer_resampler_props[i][0]);
            i8_alloc_mem_size -= u1_num_spatial_layers *
                                 sizeof(ps_codec->s_svc_ilp_data.aps_layer_resampler_props[i][0]);
        }

        ASSERT(i8_alloc_mem_size >= 0);

        for(i = u1_num_spatial_layers - 1; i >= 0; i--)
        {
            WORD32 i4_stride;

            WORD32 i4_layer_luma_wd =
                ((DOUBLE) u4_wd / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
            WORD32 i4_layer_luma_ht =
                ((DOUBLE) u4_ht / pow(d_spatial_res_ratio, u1_num_spatial_layers - 1 - i)) + 0.99;
            WORD32 i4_layer_luma_samples =
                (ALIGN16(i4_layer_luma_wd) + PAD_WD) * (i4_layer_luma_ht + PAD_HT);
            WORD32 i4_layer_uv_wd = i4_layer_luma_wd;
            WORD32 i4_layer_uv_ht = i4_layer_luma_ht / 2.0 + 0.99;
            WORD32 i4_layer_uv_samples =
                (ALIGN16(i4_layer_uv_wd) + PAD_WD) * (i4_layer_uv_ht + PAD_HT);

            ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[i].u4_width = i4_layer_luma_wd;
            ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[i].u4_height = i4_layer_luma_ht;
            ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[i].e_color_format = IV_YUV_420SP_UV;
            ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[i].u1_bit_depth = 8;

            i4_stride = ALIGN16(i4_layer_luma_wd) + PAD_WD;
            ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[i].as_component_bufs[Y].pv_data =
                pu1_buf + PAD_LEFT + PAD_TOP * i4_stride;
            ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[i].as_component_bufs[Y].i4_data_stride =
                ALIGN16(i4_layer_luma_wd) + PAD_WD;
            pu1_buf += i4_layer_luma_samples * sizeof(UWORD8);
            i8_alloc_mem_size -= i4_layer_luma_samples * sizeof(UWORD8);

            i4_stride = ALIGN16(i4_layer_uv_wd) + PAD_WD;
            ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[i].as_component_bufs[UV].pv_data =
                pu1_buf + PAD_LEFT + PAD_TOP * i4_stride;
            ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[i].as_component_bufs[UV].i4_data_stride =
                ALIGN16(i4_layer_uv_wd) + PAD_WD;
            pu1_buf += i4_layer_uv_samples * sizeof(UWORD8);
            i8_alloc_mem_size -= i4_layer_uv_samples * sizeof(UWORD8);

            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].u4_width = i4_layer_luma_wd;
            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].u4_height = i4_layer_luma_ht;
            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].e_color_format = IV_YUV_420SP_UV;
            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].u1_bit_depth = 10;

            i4_stride = ALIGN16(i4_layer_luma_wd) + PAD_WD;
            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].as_component_bufs[Y].pv_data =
                pu1_buf + (PAD_LEFT + PAD_TOP * i4_stride) * (sizeof(WORD16) / sizeof(pu1_buf[0]));
            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].as_component_bufs[Y].i4_data_stride =
                i4_stride;
            pu1_buf += i4_layer_luma_samples * sizeof(WORD16);
            i8_alloc_mem_size -= i4_layer_luma_samples * sizeof(WORD16);

            i4_stride = ALIGN16(i4_layer_uv_wd) + PAD_WD;
            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].as_component_bufs[UV].pv_data =
                pu1_buf + (PAD_LEFT + PAD_TOP * i4_stride) * (sizeof(WORD16) / sizeof(pu1_buf[0]));
            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].as_component_bufs[UV].i4_data_stride =
                i4_stride;
            pu1_buf += i4_layer_uv_samples * sizeof(WORD16);
            i8_alloc_mem_size -= i4_layer_uv_samples * sizeof(WORD16);

            ps_codec->s_svc_ilp_data.ps_residual_bufs[i].as_component_bufs[V].pv_data = NULL;

            ASSERT(i8_alloc_mem_size >= 0);

            if(i >= 1)
            {
                for(j = 0; j < NUM_SP_COMPONENTS; j++)
                {
                    isvce_layer_resampler_props_init(
                        &ps_codec->s_svc_ilp_data.aps_layer_resampler_props[j][i],
                        d_spatial_res_ratio, i4_layer_luma_wd, i4_layer_luma_ht,
                        ps_codec->s_cfg.u4_max_level, ((COMPONENT_TYPE) j) == UV);
                }
            }
        }
    }
    else
    {
        WORD32 i4_stride;

        WORD32 i4_layer_luma_wd = u4_wd;
        WORD32 i4_layer_luma_ht = u4_ht;
        WORD32 i4_layer_luma_samples =
            (ALIGN16(i4_layer_luma_wd) + PAD_WD) * (i4_layer_luma_ht + PAD_HT);
        WORD32 i4_layer_uv_wd = i4_layer_luma_wd;
        WORD32 i4_layer_uv_ht = i4_layer_luma_ht / 2.0 + 0.99;
        WORD32 i4_layer_uv_samples = (ALIGN16(i4_layer_uv_wd) + PAD_WD) * (i4_layer_uv_ht + PAD_HT);

        ps_codec->s_svc_ilp_data.ps_residual_bufs = (yuv_buf_props_t *) pu1_buf;
        pu1_buf += sizeof(ps_codec->s_svc_ilp_data.ps_residual_bufs[0]);
        i8_alloc_mem_size -= sizeof(ps_codec->s_svc_ilp_data.ps_residual_bufs[0]);

        ASSERT(i8_alloc_mem_size >= 0);

        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].u4_width = i4_layer_luma_wd;
        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].u4_height = i4_layer_luma_ht;
        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].e_color_format = IV_YUV_420SP_UV;
        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].u1_bit_depth = 10;

        i4_stride = ALIGN16(i4_layer_luma_wd) + PAD_WD;
        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].as_component_bufs[Y].pv_data =
            pu1_buf + (PAD_LEFT + PAD_TOP * i4_stride) * (sizeof(WORD16) / sizeof(pu1_buf[0]));
        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].as_component_bufs[Y].i4_data_stride =
            i4_stride;
        pu1_buf += i4_layer_luma_samples * sizeof(WORD16);
        i8_alloc_mem_size -= i4_layer_luma_samples * sizeof(WORD16);

        i4_stride = ALIGN16(i4_layer_uv_wd) + PAD_WD;
        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].as_component_bufs[UV].pv_data =
            pu1_buf + (PAD_LEFT + PAD_TOP * i4_stride) * (sizeof(WORD16) / sizeof(pu1_buf[0]));
        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].as_component_bufs[UV].i4_data_stride =
            i4_stride;
        pu1_buf += i4_layer_uv_samples * sizeof(WORD16);
        i8_alloc_mem_size -= i4_layer_uv_samples * sizeof(WORD16);

        ps_codec->s_svc_ilp_data.ps_residual_bufs[0].as_component_bufs[V].pv_data = NULL;

        ASSERT(i8_alloc_mem_size >= 0);
    }
}

static FORCEINLINE UWORD32 isvce_get_residual_csbf(mem_fxns_t *ps_mem_fxns,
                                                   buffer_container_t *ps_comp_buf)
{
    WORD32 i;

    UWORD32 u4_csbf = 0;

    for(i = 0; i < MAX_TU_IN_MB; i++)
    {
        UWORD8 u1_zscan_idx = gau1_raster_to_zscan_map[i];
        UWORD8 u1_offset_x = (i % MAX_TU_IN_MB_ROW) * MIN_TU_SIZE;
        UWORD8 u1_offset_y = (i / MAX_TU_IN_MB_ROW) * MIN_TU_SIZE;
        WORD16 *pi2_res = ((WORD16 *) ps_comp_buf->pv_data) + u1_offset_x +
                          u1_offset_y * ps_comp_buf->i4_data_stride;
        UWORD8 u1_cbf = ps_mem_fxns->pf_nonzero_checker(
            (UWORD8 *) pi2_res, ps_comp_buf->i4_data_stride * (sizeof(WORD16) / sizeof(UWORD8)),
            MIN_TU_SIZE * (sizeof(WORD16) / sizeof(UWORD8)), MIN_TU_SIZE);

        u4_csbf |= (u1_cbf << u1_zscan_idx);
    }

    return u4_csbf;
}

/**
*******************************************************************************
*
* @brief
*  Function to update svc ilp buffers after every MB
*
* @param[in] ps_proc
*  Pointer to process context
*
*******************************************************************************
*/
void isvce_svc_ilp_buf_update(isvce_process_ctxt_t *ps_proc)
{
    isvce_codec_t *ps_codec = ps_proc->ps_codec;
    svc_params_t *ps_svc_params = &ps_codec->s_cfg.s_svc_params;

    UWORD8 u1_spatial_layer_id = ps_proc->u1_spatial_layer_id;

    if(ps_svc_params->u1_num_spatial_layers > 1)
    {
        buffer_container_t s_src;
        buffer_container_t s_dst;

        WORD32 i;

        svc_ilp_data_t *ps_svc_ilp_data = &ps_codec->s_svc_ilp_data;
        isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
        mem_fxns_t *ps_mem_fxns = &ps_isa_dependent_fxns->s_mem_fxns;
        yuv_buf_props_t *ps_residual_buf =
            &ps_codec->s_svc_ilp_data.ps_residual_bufs[u1_spatial_layer_id];

        WORD32 i4_mb_x = ps_proc->i4_mb_x;
        WORD32 i4_mb_y = ps_proc->i4_mb_y;

        ASSERT(ps_proc->s_rec_buf_props.e_color_format == IV_YUV_420SP_UV);

        if(u1_spatial_layer_id < (ps_svc_params->u1_num_spatial_layers - 1))
        {
            if(ps_proc->ps_mb_info->u1_is_intra)
            {
                for(i = 0; i < NUM_SP_COMPONENTS; i++)
                {
                    UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) i));

                    s_src = ps_proc->s_rec_buf_props.as_component_bufs[i];

                    s_dst.i4_data_stride = ps_svc_ilp_data->ps_intra_recon_bufs[u1_spatial_layer_id]
                                               .as_component_bufs[i]
                                               .i4_data_stride;
                    s_dst.pv_data =
                        ((UWORD8 *) ps_svc_ilp_data->ps_intra_recon_bufs[u1_spatial_layer_id]
                             .as_component_bufs[i]
                             .pv_data) +
                        i4_mb_x * MB_SIZE +
                        i4_mb_y * (MB_SIZE >> u1_is_chroma) * s_dst.i4_data_stride;

                    ps_mem_fxns->pf_copy_2d((UWORD8 *) s_dst.pv_data, s_dst.i4_data_stride,
                                            (UWORD8 *) s_src.pv_data, s_src.i4_data_stride, MB_SIZE,
                                            (MB_SIZE >> u1_is_chroma));
                }
            }
            else
            {
                for(i = 0; i < NUM_SP_COMPONENTS; i++)
                {
                    UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) i));

                    s_dst.i4_data_stride = ps_svc_ilp_data->ps_intra_recon_bufs[u1_spatial_layer_id]
                                               .as_component_bufs[i]
                                               .i4_data_stride;
                    s_dst.pv_data =
                        ((UWORD8 *) ps_svc_ilp_data->ps_intra_recon_bufs[u1_spatial_layer_id]
                             .as_component_bufs[i]
                             .pv_data) +
                        i4_mb_x * MB_SIZE +
                        i4_mb_y * (MB_SIZE >> u1_is_chroma) * s_dst.i4_data_stride;

                    ps_mem_fxns->pf_memset_2d((UWORD8 *) s_dst.pv_data, s_dst.i4_data_stride, 0,
                                              MB_SIZE, (MB_SIZE >> u1_is_chroma));
                }
            }
        }

        if(ENABLE_RESIDUAL_PREDICTION && (ps_proc->i4_slice_type != ISLICE) &&
           (u1_spatial_layer_id < (ps_svc_params->u1_num_spatial_layers - 1)))
        {
            if(ps_proc->ps_mb_info->u1_is_intra || (ps_proc->ps_mb_info->u2_mb_type == PSKIP) ||
               (ps_proc->ps_mb_info->u2_mb_type == BSKIP))
            {
                for(i = 0; i < NUM_SP_COMPONENTS; i++)
                {
                    buffer_container_t *ps_comp_buf;

                    WORD16 *pi2_res;

                    UWORD8 u1_is_chroma = (Y != ((COMPONENT_TYPE) i));

                    ps_comp_buf = &ps_residual_buf->as_component_bufs[u1_is_chroma ? UV : Y];
                    pi2_res =
                        ((WORD16 *) ps_comp_buf->pv_data) + ps_proc->i4_mb_x * MB_SIZE +
                        ps_proc->i4_mb_y * (MB_SIZE >> u1_is_chroma) * ps_comp_buf->i4_data_stride;

                    ps_mem_fxns->pf_memset_2d(
                        (UWORD8 *) pi2_res,
                        ps_comp_buf->i4_data_stride * (sizeof(WORD16) / sizeof(UWORD8)), 0,
                        MB_SIZE * (sizeof(WORD16) / sizeof(UWORD8)), MB_SIZE >> u1_is_chroma);
                }
            }
        }

        if(ENABLE_RESIDUAL_PREDICTION && (u1_spatial_layer_id > 0) &&
           !(ps_proc->ps_mb_info->u1_is_intra || (ps_proc->ps_mb_info->u2_mb_type == PSKIP) ||
             (ps_proc->ps_mb_info->u2_mb_type == BSKIP)))
        {
            s_src = ps_residual_buf->as_component_bufs[Y];
            s_src.pv_data = ((WORD16 *) s_src.pv_data) + ps_proc->i4_mb_x * MB_SIZE +
                            ps_proc->i4_mb_y * MB_SIZE * s_src.i4_data_stride;

            ps_proc->ps_mb_info->u4_res_csbp = isvce_get_residual_csbf(ps_mem_fxns, &s_src);
        }
        else
        {
            ps_proc->ps_mb_info->u4_res_csbp = 0;
        }
    }
    else
    {
        ps_proc->ps_mb_info->u4_res_csbp = 0;
    }
}

/*
 * Padding has a one MB row dependency on deblock  which
 * in turn has a one MB row dependency on encode
 */
static IH264E_ERROR_T isvce_pad_frame(isvce_process_ctxt_t *ps_proc, yuv_buf_props_t *ps_pad_buf)
{
    /* codec context */
    isvce_codec_t *ps_codec = ps_proc->ps_codec;

    WORD32 i4_element_size = (ps_pad_buf->u1_bit_depth > 8) ? 2 : 1;

    /* src buffers luma */
    WORD32 i4_luma_stride = ps_pad_buf->as_component_bufs[0].i4_data_stride * i4_element_size;
    UWORD8 *pu1_curr_pic_luma = (UWORD8 *) (ps_pad_buf->as_component_bufs[0].pv_data);

    /* src buffers chroma */
    WORD32 i4_chroma_stride = ps_pad_buf->as_component_bufs[1].i4_data_stride * i4_element_size;
    UWORD8 *pu1_curr_pic_chroma = (UWORD8 *) (ps_pad_buf->as_component_bufs[1].pv_data);

    WORD32 i4_bottom_offset_luma = ps_pad_buf->u4_height * i4_luma_stride;
    WORD32 i4_bottom_offset_chroma = (ps_pad_buf->u4_height >> 1) * i4_chroma_stride;

    /* Pad left */
    ps_codec->pf_pad_left_luma(pu1_curr_pic_luma, i4_luma_stride, ps_pad_buf->u4_height,
                               PAD_LEFT * i4_element_size);
    ps_codec->pf_pad_left_chroma(pu1_curr_pic_chroma, i4_chroma_stride, ps_pad_buf->u4_height >> 1,
                                 PAD_LEFT * i4_element_size);

    /* Pad right */
    ps_codec->pf_pad_right_luma(pu1_curr_pic_luma + ps_pad_buf->u4_width * i4_element_size,
                                i4_luma_stride, ps_pad_buf->u4_height, PAD_RIGHT * i4_element_size);
    ps_codec->pf_pad_right_chroma(pu1_curr_pic_chroma + ps_pad_buf->u4_width * i4_element_size,
                                  i4_chroma_stride, ps_pad_buf->u4_height >> 1,
                                  PAD_RIGHT * i4_element_size);

    /* Pad top */
    ps_codec->pf_pad_top(pu1_curr_pic_luma - (PAD_LEFT * i4_element_size), i4_luma_stride,
                         (ps_pad_buf->u4_width + PAD_WD) * i4_element_size, PAD_TOP);
    ps_codec->pf_pad_top(pu1_curr_pic_chroma - (PAD_LEFT * i4_element_size), i4_chroma_stride,
                         (ps_pad_buf->u4_width + PAD_WD) * i4_element_size, PAD_TOP >> 1);

    /* Pad bottom */
    ps_codec->pf_pad_bottom(
        pu1_curr_pic_luma + i4_bottom_offset_luma - (PAD_LEFT * i4_element_size), i4_luma_stride,
        (ps_pad_buf->u4_width + PAD_WD) * i4_element_size, PAD_BOT);
    ps_codec->pf_pad_bottom(
        pu1_curr_pic_chroma + i4_bottom_offset_chroma - (PAD_LEFT * i4_element_size),
        i4_chroma_stride, (ps_pad_buf->u4_width + PAD_WD) * i4_element_size, PAD_BOT >> 1);

    return IH264E_SUCCESS;
}

void isvce_svc_pad_frame(isvce_process_ctxt_t *ps_proc)
{
    isvce_codec_t *ps_codec = ps_proc->ps_codec;

    isvce_pad_frame(ps_proc, &(ps_proc->s_rec_pic_buf_props));

    if(ps_proc->s_svc_params.u1_num_spatial_layers > 1)
    {
        isvce_pad_frame(
            ps_proc, &(ps_codec->s_svc_ilp_data.ps_intra_recon_bufs[ps_proc->u1_spatial_layer_id]));
        isvce_pad_frame(ps_proc,
                        &(ps_codec->s_svc_ilp_data.ps_residual_bufs[ps_proc->u1_spatial_layer_id]));
    }
}

/**
*******************************************************************************
*
* @brief
*  Initialize AIR mb frame Map
*
* @par Description:
*  Initialize AIR mb frame map
*  MB frame map indicates which frame an Mb should be coded as intra according
*to AIR
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @returns  error_status
*
* @remarks
*
*
*******************************************************************************
*/
IH264E_ERROR_T isvce_init_air_map(isvce_codec_t *ps_codec)
{
    /* intra refresh map */
    UWORD16 *pu2_intr_rfrsh_map = ps_codec->pu2_intr_rfrsh_map;

    /* air mode */
    IVE_AIR_MODE_T air_mode = ps_codec->s_cfg.e_air_mode;

    /* refresh period */
    UWORD32 air_period = ps_codec->s_cfg.u4_air_refresh_period;

    /* mb cnt */
    UWORD32 u4_mb_cnt = ps_codec->s_cfg.i4_wd_mbs * ps_codec->s_cfg.i4_ht_mbs;

    /* temp var */
    UWORD32 curr_mb, seed_rand = 1;

    switch(air_mode)
    {
        case IVE_AIR_MODE_CYCLIC:

            for(curr_mb = 0; curr_mb < u4_mb_cnt; curr_mb++)
            {
                pu2_intr_rfrsh_map[curr_mb] = curr_mb % air_period;
            }
            break;

        case IVE_AIR_MODE_RANDOM:

            for(curr_mb = 0; curr_mb < u4_mb_cnt; curr_mb++)
            {
                seed_rand = (seed_rand * 32719 + 3) % 32749;
                pu2_intr_rfrsh_map[curr_mb] = seed_rand % air_period;
            }
            break;

        default:

            break;
    }

    return IH264E_SUCCESS;
}

/**
******************************************************************************
*
* @brief
*  derivation process for macroblock availability
*
* @par   Description
*  Calculates the availability of the left, top, topright and topleft macroblocks.
*
* @param[in] ps_proc_ctxt
*  pointer to proc context (handle)
*
* @remarks Based on section 6.4.5 in H264 spec
*
* @return  none
*
******************************************************************************
*/
void isvce_derive_nghbr_avbl_of_mbs(isvce_process_ctxt_t *ps_proc)
{
    UWORD8 *pu1_slice_idx_curr = ps_proc->pu1_slice_idx;
    UWORD8 *pu1_slice_idx_b;
    UWORD8 *pu1_slice_idx_a;
    UWORD8 *pu1_slice_idx_c;
    UWORD8 *pu1_slice_idx_d;
    block_neighbors_t *ps_ngbr_avbl;
    WORD32 i4_mb_x, i4_mb_y;
    WORD32 i4_wd_mbs;

    i4_mb_x = ps_proc->i4_mb_x;
    i4_mb_y = ps_proc->i4_mb_y;

    i4_wd_mbs = ps_proc->i4_wd_mbs;

    pu1_slice_idx_curr += (i4_mb_y * i4_wd_mbs) + i4_mb_x;
    pu1_slice_idx_a = pu1_slice_idx_curr - 1;
    pu1_slice_idx_b = pu1_slice_idx_curr - i4_wd_mbs;
    pu1_slice_idx_c = pu1_slice_idx_b + 1;
    pu1_slice_idx_d = pu1_slice_idx_b - 1;
    ps_ngbr_avbl = ps_proc->ps_ngbr_avbl;

    /**********************************************************************/
    /* The macroblock is marked as available, unless one of the following */
    /* conditions is true in which case the macroblock shall be marked as */
    /* not available.                                                     */
    /* 1. mbAddr < 0                                                      */
    /* 2  mbAddr > CurrMbAddr                                             */
    /* 3. the macroblock with address mbAddr belongs to a different slice */
    /* than the macroblock with address CurrMbAddr                        */
    /**********************************************************************/

    /* left macroblock availability */
    if(i4_mb_x == 0)
    { /* macroblocks along first column */
        ps_ngbr_avbl->u1_mb_a = 0;
    }
    else
    { /* macroblocks belong to same slice? */
        if(*pu1_slice_idx_a != *pu1_slice_idx_curr)
            ps_ngbr_avbl->u1_mb_a = 0;
        else
            ps_ngbr_avbl->u1_mb_a = 1;
    }

    /* top macroblock availability */
    if(i4_mb_y == 0)
    { /* macroblocks along first row */
        ps_ngbr_avbl->u1_mb_b = 0;
    }
    else
    { /* macroblocks belong to same slice? */
        if(*pu1_slice_idx_b != *pu1_slice_idx_curr)
            ps_ngbr_avbl->u1_mb_b = 0;
        else
            ps_ngbr_avbl->u1_mb_b = 1;
    }

    /* top right macroblock availability */
    if(i4_mb_x == i4_wd_mbs - 1 || i4_mb_y == 0)
    { /* macroblocks along last column */
        ps_ngbr_avbl->u1_mb_c = 0;
    }
    else
    { /* macroblocks belong to same slice? */
        if(*pu1_slice_idx_c != *pu1_slice_idx_curr)
            ps_ngbr_avbl->u1_mb_c = 0;
        else
            ps_ngbr_avbl->u1_mb_c = 1;
    }

    /* top left macroblock availability */
    if(i4_mb_x == 0 || i4_mb_y == 0)
    { /* macroblocks along first column */
        ps_ngbr_avbl->u1_mb_d = 0;
    }
    else
    { /* macroblocks belong to same slice? */
        if(*pu1_slice_idx_d != *pu1_slice_idx_curr)
            ps_ngbr_avbl->u1_mb_d = 0;
        else
            ps_ngbr_avbl->u1_mb_d = 1;
    }
}

/**
*******************************************************************************
*
* @brief
*  Codec level initializations
*
* @par Description:
*  Initializes the codec with parameters that needs to be set before encoding
*  first frame
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @param[in] ps_inp_buf
*  Pointer to input buffer context
*
* @returns  error_status
*
* @remarks
*
*
*******************************************************************************
*/
IH264E_ERROR_T isvce_codec_init(isvce_codec_t *ps_codec)
{
    isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
    enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;
    WORD8 i;

    /********************************************************************
     *                     INITIALIZE CODEC CONTEXT                     *
     ********************************************************************/
    /* encoder presets */
    if(ps_codec->s_cfg.u4_enc_speed_preset != IVE_CONFIG)
    {
        if(ps_codec->s_cfg.u4_enc_speed_preset == IVE_SLOWEST)
        { /* high quality */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;
            ps_codec->s_cfg.u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 1;
            if(!FORCE_FAST_INTRA4X4)
            {
                ps_enc_loop_fxns->apf_luma_energy_compaction[1] =
                    isvce_code_luma_intra_macroblock_4x4_rdopt_on;
            }

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 1;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if(ps_codec->s_cfg.u4_enc_speed_preset == IVE_NORMAL)
        { /* normal */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;
            ps_codec->s_cfg.u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 1;

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 1;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if(ps_codec->s_cfg.u4_enc_speed_preset == IVE_FAST)
        { /* normal */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;
            ps_codec->s_cfg.u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 0;

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 1;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 1;
        }
        else if(ps_codec->s_cfg.u4_enc_speed_preset == IVE_HIGH_SPEED)
        { /* fast */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;
            ps_codec->s_cfg.u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 0;

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 0;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if(ps_codec->s_cfg.u4_enc_speed_preset == IVE_FASTEST)
        { /* fastest */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 0;

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 0;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 1;
        }
    }

    /*****************************************************************
     * Initialize AIR inside codec
     *****************************************************************/
    if(IVE_AIR_MODE_NONE != ps_codec->s_cfg.e_air_mode)
    {
        isvce_init_air_map(ps_codec);

        ps_codec->i4_air_pic_cnt = -1;
    }

    /****************************************************/
    /*           INITIALIZE RATE CONTROL                */
    /****************************************************/
    {
        for(i = 0; i < MAX_NUM_SPATIAL_LAYERS; i++)
        {
            UWORD8 au1_init_qp[MAX_PIC_TYPE];
            UWORD8 au1_min_max_qp[2 * MAX_PIC_TYPE];
            UWORD8 au1_min_max_avc_qp[2 * MAX_PIC_TYPE];

            /* update rc lib with modified qp */
            au1_init_qp[0] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp[i]];
            au1_init_qp[1] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp[i]];
            au1_init_qp[2] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp[i]];

            au1_min_max_qp[2 * I_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp_min[i]];
            au1_min_max_qp[2 * I_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp_max[i]];

            au1_min_max_qp[2 * P_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp_min[i]];
            au1_min_max_qp[2 * P_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp_max[i]];

            au1_min_max_qp[2 * B_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp_min[i]];
            au1_min_max_qp[2 * B_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp_max[i]];

            /* get rc mode */
            switch(ps_codec->s_cfg.e_rc_mode)
            {
                case IVE_RC_STORAGE:
                    ps_codec->s_rate_control.e_rc_type = VBR_STORAGE;
                    break;
                case IVE_RC_CBR_NON_LOW_DELAY:
                    ps_codec->s_rate_control.e_rc_type = CBR_NLDRC;
                    break;
                case IVE_RC_CBR_LOW_DELAY:
                    ps_codec->s_rate_control.e_rc_type = CBR_LDRC;
                    break;
                case IVE_RC_NONE:
                    ps_codec->s_rate_control.e_rc_type = CONST_QP;
                    break;
                default:
                    break;
            }

            ps_codec->u1_enable_init_qp = DEFAULT_INIT_QP;

            /* init rate control */
            isvce_rc_init(
                ps_codec->s_rate_control.apps_rate_control_api[i],
                ps_codec->s_rate_control.pps_frame_time, ps_codec->s_rate_control.pps_time_stamp,
                ps_codec->s_rate_control.pps_pd_frm_rate, ps_codec->s_cfg.u4_max_framerate,
                ps_codec->s_cfg.u4_src_frame_rate, ps_codec->s_cfg.u4_tgt_frame_rate,
                ps_codec->s_rate_control.e_rc_type, ps_codec->s_cfg.au4_target_bitrate[i],
                ps_codec->s_cfg.au4_max_bitrate[i], ps_codec->s_cfg.au4_vbv_buffer_delay[i],
                ps_codec->s_cfg.u4_i_frm_interval, ps_codec->s_cfg.u4_num_bframes + 1, au1_init_qp,
                ps_codec->s_cfg.u4_num_bframes + 2, au1_min_max_qp,
                MAX(ps_codec->s_cfg.u4_max_level,
                    (UWORD32) ih264e_get_min_level(ps_codec->s_cfg.u4_max_wd,
                                                   ps_codec->s_cfg.u4_max_ht)));

            au1_min_max_avc_qp[2 * I_PIC] = ps_codec->s_cfg.au4_i_qp_min[i];
            au1_min_max_avc_qp[2 * I_PIC + 1] = ps_codec->s_cfg.au4_i_qp_max[i];

            au1_min_max_avc_qp[2 * P_PIC] = ps_codec->s_cfg.au4_p_qp_min[i];
            au1_min_max_avc_qp[2 * P_PIC + 1] = ps_codec->s_cfg.au4_p_qp_max[i];

            au1_min_max_avc_qp[2 * B_PIC] = ps_codec->s_cfg.au4_b_qp_min[i];
            au1_min_max_avc_qp[2 * B_PIC + 1] = ps_codec->s_cfg.au4_b_qp_max[i];

            irc_change_qp_constraints(ps_codec->s_rate_control.apps_rate_control_api[i],
                                      au1_min_max_qp, au1_min_max_avc_qp);
        }
    }

    /* recon stride */
    ps_codec->i4_rec_strd = ALIGN16(ps_codec->s_cfg.u4_max_wd) + PAD_WD;

    /* max ref and reorder cnt */
    ps_codec->i4_ref_buf_cnt = ps_codec->s_cfg.u4_max_ref_cnt + ps_codec->s_cfg.u4_max_reorder_cnt;
    ps_codec->i4_ref_buf_cnt += MAX_CTXT_SETS;
    ps_codec->i4_ref_buf_cnt += ps_codec->s_cfg.s_svc_params.u1_num_temporal_layers;

    DEBUG_HISTOGRAM_INIT();

    /* Init dependecy vars */
    ps_codec->i4_last_inp_buff_received = 0;

    /* At codec start no IDR is pending */
    ps_codec->i4_pending_idr_flag = 0;

    for(i = 0; i < ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers - 1; i++)
    {
        ps_codec->au4_constrained_intra_pred[i] = 1;
    }

    ps_codec->au4_constrained_intra_pred[ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers - 1] =
        0;

    return IH264E_SUCCESS;
}

/**
*******************************************************************************
*
* @brief update encoder configuration parameters
*
* @par Description:
*  updates encoder configuration parameters from the given config set.
*  Initialize/reinitialize codec parameters according to new configurations.
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @param[in] ps_cfg
*  Pointer to config param set
*
* @remarks none
*
*******************************************************************************
*/
IH264E_ERROR_T isvce_codec_update_config(isvce_codec_t *ps_codec, isvce_cfg_params_t *ps_cfg)
{
    /* config params */
    isvce_cfg_params_t *ps_curr_cfg = &ps_codec->s_cfg;

    /* error status */
    IH264E_ERROR_T err = IH264E_SUCCESS;

    /* temp var */
    UWORD32 u4_init_rc = 0;

    WORD8 i;

    /***********************/
    /* UPDATE CODEC CONFIG */
    /***********************/
    if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_DIMENSIONS)
    {
        UWORD32 wd_aln = ALIGN16(ps_cfg->u4_wd);
        UWORD32 ht_aln = ALIGN16(ps_cfg->u4_ht);

        if(ps_curr_cfg->u4_wd != wd_aln || ps_curr_cfg->u4_ht != ht_aln ||
           ps_curr_cfg->u4_disp_wd != ps_cfg->u4_disp_wd ||
           ps_curr_cfg->u4_disp_ht != ps_cfg->u4_disp_ht)
        {
            ps_curr_cfg->u4_wd = wd_aln;
            ps_curr_cfg->u4_ht = ht_aln;

            ps_curr_cfg->u4_disp_wd = ps_cfg->u4_disp_wd;
            ps_curr_cfg->u4_disp_ht = ps_cfg->u4_disp_ht;

            ps_curr_cfg->i4_wd_mbs = ps_curr_cfg->u4_wd >> 4;
            ps_curr_cfg->i4_ht_mbs = ps_curr_cfg->u4_ht >> 4;

            ps_codec->i4_rec_strd = ALIGN16(ps_cfg->u4_wd) + PAD_WD;

            /* If number of MBs in a frame changes the air map also changes.
             * Hence recompute air map also reset air pic cnt */
            if(ps_codec->s_cfg.e_air_mode != IVE_AIR_MODE_NONE)
            {
                /* re-init the air map */
                isvce_init_air_map(ps_codec);

                /* reset air counter */
                ps_codec->i4_air_pic_cnt = -1;
            }

            /* initialize mv bank buffer manager */
            err = isvce_svc_au_data_mgr_add_bufs(ps_codec);
            if(err != IH264E_SUCCESS) return err;

            /* initialize ref bank buffer manager */
            err = isvce_svc_au_buf_mgr_add_bufs(ps_codec);
            if(err != IH264E_SUCCESS) return err;

            /* since dimension changed, start new sequence by forcing IDR */
            ps_codec->force_curr_frame_type = IV_IDR_FRAME;

            /* in case dimension changes, we need to reinitialize RC as the
             * old model shall not fit further */
            u4_init_rc = 1;

            /* when the dimension changes, the header needs to be regenerated */
            ps_codec->i4_gen_header = 1;
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_FRAMERATE)
    {
        /* temp var */
        UWORD32 u4_src_ticks, u4_tgt_ticks;

        u4_src_ticks = ih264e_frame_time_get_src_ticks(ps_codec->s_rate_control.pps_frame_time);

        u4_tgt_ticks = ih264e_frame_time_get_tgt_ticks(ps_codec->s_rate_control.pps_frame_time);

        /* Change frame rate */
        if(ps_codec->s_cfg.u4_src_frame_rate != ps_cfg->u4_src_frame_rate * 1000)
        {
            ps_codec->s_cfg.u4_src_frame_rate = ps_cfg->u4_src_frame_rate * 1000;

            ih264e_frame_time_update_src_frame_rate(ps_codec->s_rate_control.pps_frame_time,
                                                    ps_codec->s_cfg.u4_src_frame_rate);

            ih264_time_stamp_update_frame_rate(ps_codec->s_rate_control.pps_time_stamp,
                                               ps_codec->s_cfg.u4_src_frame_rate);

            for(i = 0; i < ps_cfg->s_svc_params.u1_num_spatial_layers; i++)
            {
                irc_change_frame_rate(ps_codec->s_rate_control.apps_rate_control_api[i],
                                      ps_codec->s_cfg.u4_src_frame_rate, u4_src_ticks,
                                      u4_tgt_ticks);
            }
        }

        if(ps_codec->s_cfg.u4_tgt_frame_rate != ps_cfg->u4_tgt_frame_rate * 1000)
        {
            ps_codec->s_cfg.u4_tgt_frame_rate = ps_cfg->u4_tgt_frame_rate * 1000;

            ih264e_frame_time_update_tgt_frame_rate(ps_codec->s_rate_control.pps_frame_time,
                                                    ps_codec->s_cfg.u4_tgt_frame_rate);

            for(i = 0; i < ps_cfg->s_svc_params.u1_num_spatial_layers; i++)
            {
                irc_change_frame_rate(ps_codec->s_rate_control.apps_rate_control_api[i],
                                      ps_codec->s_cfg.u4_src_frame_rate, u4_src_ticks,
                                      u4_tgt_ticks);

                irc_change_frm_rate_for_bit_alloc(ps_codec->s_rate_control.apps_rate_control_api[i],
                                                  ps_codec->s_cfg.u4_tgt_frame_rate);
            }
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_BITRATE)
    {
        for(i = 0; i < MAX_NUM_SPATIAL_LAYERS; i++)
        {
            if(ps_curr_cfg->au4_target_bitrate[i] != ps_cfg->au4_target_bitrate[i])
            {
                if(IVE_RC_NONE != ps_curr_cfg->e_rc_mode)
                    irc_change_avg_bit_rate(ps_codec->s_rate_control.apps_rate_control_api[i],
                                            ps_cfg->au4_target_bitrate[i]);

                ps_curr_cfg->au4_target_bitrate[i] = ps_cfg->au4_target_bitrate[i];
            }
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_FRAMETYPE)
    {
        switch(ps_cfg->e_frame_type)
        {
            case IV_I_FRAME:
                ps_codec->force_curr_frame_type = IV_I_FRAME;
                break;

            case IV_IDR_FRAME:
                ps_codec->force_curr_frame_type = IV_IDR_FRAME;
                break;

            case IV_P_FRAME:
            default:
                break;
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_ME_PARAMS)
    {
        if(ps_curr_cfg->u4_enc_speed_preset == IVE_CONFIG)
        {
            ps_codec->s_cfg.u4_enable_hpel = ps_cfg->u4_enable_hpel;
            ps_codec->s_cfg.u4_enable_fast_sad = ps_cfg->u4_enable_fast_sad;
            ps_codec->s_cfg.u4_me_speed_preset = ps_cfg->u4_me_speed_preset;
            ps_codec->s_cfg.u4_enable_qpel = ps_cfg->u4_enable_qpel;
        }
        else if(ps_curr_cfg->u4_enc_speed_preset == IVE_FASTEST)
        {
            ps_codec->s_cfg.u4_enable_fast_sad = ps_cfg->u4_enable_fast_sad;
        }
        ps_codec->s_cfg.u4_srch_rng_x = ps_cfg->u4_srch_rng_x;
        ps_codec->s_cfg.u4_srch_rng_y = ps_cfg->u4_srch_rng_y;

        if(ps_codec->s_cfg.u4_enable_alt_ref != ps_cfg->u4_enable_alt_ref)
        {
            ps_codec->s_cfg.u4_enable_alt_ref = ps_cfg->u4_enable_alt_ref;
            ps_codec->u4_is_curr_frm_ref = 1;
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_IPE_PARAMS)
    {
        ps_curr_cfg->u4_enc_speed_preset = ps_cfg->u4_enc_speed_preset;

        if(ps_curr_cfg->u4_enc_speed_preset == IVE_SLOWEST)
        {
            isa_dependent_fxns_t *ps_isa_dependent_fxns = &ps_codec->s_isa_dependent_fxns;
            enc_loop_fxns_t *ps_enc_loop_fxns = &ps_isa_dependent_fxns->s_enc_loop_fxns;

            /* enable diamond search */
            ps_curr_cfg->u4_me_speed_preset = DMND_SRCH;
            ps_curr_cfg->u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_curr_cfg->u4_enable_intra_4x4 = 1;
            ps_enc_loop_fxns->apf_luma_energy_compaction[1] =
                isvce_code_luma_intra_macroblock_4x4_rdopt_on;

            /* sub pel off */
            ps_curr_cfg->u4_enable_hpel = 1;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if(ps_curr_cfg->u4_enc_speed_preset == IVE_NORMAL)
        { /* normal */
            /* enable diamond search */
            ps_curr_cfg->u4_me_speed_preset = DMND_SRCH;
            ps_curr_cfg->u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_curr_cfg->u4_enable_intra_4x4 = 1;

            /* sub pel off */
            ps_curr_cfg->u4_enable_hpel = 1;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if(ps_curr_cfg->u4_enc_speed_preset == IVE_FAST)
        { /* normal */
            /* enable diamond search */
            ps_curr_cfg->u4_me_speed_preset = DMND_SRCH;
            ps_curr_cfg->u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_curr_cfg->u4_enable_intra_4x4 = 0;

            /* sub pel off */
            ps_curr_cfg->u4_enable_hpel = 1;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 1;
        }
        else if(ps_curr_cfg->u4_enc_speed_preset == IVE_HIGH_SPEED)
        { /* fast */
            /* enable diamond search */
            ps_curr_cfg->u4_me_speed_preset = DMND_SRCH;
            ps_curr_cfg->u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_curr_cfg->u4_enable_intra_4x4 = 0;

            /* sub pel off */
            ps_curr_cfg->u4_enable_hpel = 0;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if(ps_curr_cfg->u4_enc_speed_preset == IVE_FASTEST)
        { /* fastest */
            /* enable diamond search */
            ps_curr_cfg->u4_me_speed_preset = DMND_SRCH;
            // u4_num_layers = 4;

            /* disable intra 4x4 */
            ps_curr_cfg->u4_enable_intra_4x4 = 0;

            /* sub pel off */
            ps_curr_cfg->u4_enable_hpel = 0;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 1;
        }
        else if(ps_curr_cfg->u4_enc_speed_preset == IVE_CONFIG)
        {
            ps_curr_cfg->u4_enable_intra_4x4 = ps_cfg->u4_enable_intra_4x4;
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_GOP_PARAMS)
    {
        if(ps_curr_cfg->u4_i_frm_interval != ps_cfg->u4_i_frm_interval)
        {
            ps_curr_cfg->u4_i_frm_interval = ps_cfg->u4_i_frm_interval;

            /* reset air counter */
            ps_codec->i4_air_pic_cnt = -1;

            /* re-init air map */
            isvce_init_air_map(ps_codec);

            /*Effect intra frame interval change*/
            for(i = 0; i < ps_cfg->s_svc_params.u1_num_spatial_layers; i++)
            {
                irc_change_intra_frm_int_call(ps_codec->s_rate_control.apps_rate_control_api[i],
                                              ps_curr_cfg->u4_i_frm_interval);
            }
        }

        ps_curr_cfg->u4_idr_frm_interval = ps_cfg->u4_idr_frm_interval;
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_DEBLOCK_PARAMS)
    {
        ps_curr_cfg->u4_disable_deblock_level = ps_cfg->u4_disable_deblock_level;
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_QP)
    {
        for(i = 0; i < ps_cfg->s_svc_params.u1_num_spatial_layers; i++)
        {
            UWORD8 au1_init_qp[MAX_PIC_TYPE];
            UWORD8 au1_min_max_qp[2 * MAX_PIC_TYPE];
            UWORD8 au1_min_max_avc_qp[2 * MAX_PIC_TYPE];

            ps_codec->s_cfg.au4_i_qp_max[i] = ps_cfg->au4_i_qp_max[i];
            ps_codec->s_cfg.au4_i_qp_min[i] = ps_cfg->au4_i_qp_min[i];
            ps_codec->s_cfg.au4_i_qp[i] = ps_cfg->au4_i_qp[i];

            ps_codec->s_cfg.au4_p_qp_max[i] = ps_cfg->au4_p_qp_max[i];
            ps_codec->s_cfg.au4_p_qp_min[i] = ps_cfg->au4_p_qp_min[i];
            ps_codec->s_cfg.au4_p_qp[i] = ps_cfg->au4_p_qp[i];

            ps_codec->s_cfg.au4_b_qp_max[i] = ps_cfg->au4_b_qp_max[i];
            ps_codec->s_cfg.au4_b_qp_min[i] = ps_cfg->au4_b_qp_min[i];
            ps_codec->s_cfg.au4_b_qp[i] = ps_cfg->au4_b_qp[i];

            /* update rc lib with modified qp */
            au1_init_qp[0] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp[i]];
            au1_init_qp[1] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp[i]];
            au1_init_qp[2] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp[i]];

            irc_change_init_qp(ps_codec->s_rate_control.apps_rate_control_api[i], au1_init_qp);

            au1_min_max_qp[2 * I_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp_min[i]];
            au1_min_max_qp[2 * I_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp_max[i]];

            au1_min_max_qp[2 * P_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp_min[i]];
            au1_min_max_qp[2 * P_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp_max[i]];

            au1_min_max_qp[2 * B_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp_min[i]];
            au1_min_max_qp[2 * B_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp_max[i]];

            au1_min_max_avc_qp[2 * I_PIC] = ps_codec->s_cfg.au4_i_qp_min[i];
            au1_min_max_avc_qp[2 * I_PIC + 1] = ps_codec->s_cfg.au4_i_qp_max[i];

            au1_min_max_avc_qp[2 * P_PIC] = ps_codec->s_cfg.au4_p_qp_min[i];
            au1_min_max_avc_qp[2 * P_PIC + 1] = ps_codec->s_cfg.au4_p_qp_max[i];

            au1_min_max_avc_qp[2 * B_PIC] = ps_codec->s_cfg.au4_b_qp_min[i];
            au1_min_max_avc_qp[2 * B_PIC + 1] = ps_codec->s_cfg.au4_b_qp_max[i];

            irc_change_qp_constraints(ps_codec->s_rate_control.apps_rate_control_api[i],
                                      au1_min_max_qp, au1_min_max_avc_qp);
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_ENC_MODE)
    {
        ps_codec->s_cfg.e_enc_mode = ps_cfg->e_enc_mode;

        if(ps_codec->s_cfg.e_enc_mode == IVE_ENC_MODE_HEADER)
        {
            ps_codec->i4_header_mode = 1;
            ps_codec->s_cfg.e_enc_mode = IVE_ENC_MODE_PICTURE;
        }
        else
        {
            ps_codec->i4_header_mode = 0;
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_VBV_PARAMS &&
            IVE_RC_NONE != ps_codec->s_cfg.e_rc_mode)
    {
        for(i = 0; i < ps_cfg->s_svc_params.u1_num_spatial_layers; i++)
        {
            ps_codec->s_cfg.au4_vbv_buffer_delay[i] = ps_cfg->au4_vbv_buffer_delay[i];
        }
        // irc_change_buffer_delay(ps_codec->s_rate_control.pps_rate_control_api,
        // ps_codec->s_cfg.u4_vbv_buffer_delay);

        // TODO: remove this when the support for changing buffer dynamically
        // is yet to be added.
        u4_init_rc = 1;
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_AIR_PARAMS)
    {
        if(ps_curr_cfg->e_air_mode != ps_cfg->e_air_mode ||
           ps_curr_cfg->u4_air_refresh_period != ps_cfg->u4_air_refresh_period)
        {
            ps_curr_cfg->e_air_mode = ps_cfg->e_air_mode;
            ps_curr_cfg->u4_air_refresh_period = ps_cfg->u4_air_refresh_period;

            isvce_init_air_map(ps_codec);

            /* reset air counter */
            ps_codec->i4_air_pic_cnt = -1;
        }
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_PROFILE_PARAMS)
    {
        ps_codec->s_cfg.e_profile = ps_cfg->e_profile;
        ps_codec->s_cfg.u4_entropy_coding_mode = ps_cfg->u4_entropy_coding_mode;
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_NUM_CORES)
    {
        ps_codec->s_cfg.u4_num_cores = ps_cfg->u4_num_cores;
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_VUI_PARAMS)
    {
        ps_codec->s_cfg.s_vui = ps_cfg->s_vui;
    }

    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_SEI_MDCV_PARAMS)
    {
        ps_codec->s_cfg.s_sei.u1_sei_mdcv_params_present_flag =
            ps_cfg->s_sei.u1_sei_mdcv_params_present_flag;
        ps_codec->s_cfg.s_sei.s_sei_mdcv_params = ps_cfg->s_sei.s_sei_mdcv_params;
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_SEI_CLL_PARAMS)
    {
        ps_codec->s_cfg.s_sei.u1_sei_cll_params_present_flag =
            ps_cfg->s_sei.u1_sei_cll_params_present_flag;
        ps_codec->s_cfg.s_sei.s_sei_cll_params = ps_cfg->s_sei.s_sei_cll_params;
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_SEI_AVE_PARAMS)
    {
        ps_codec->s_cfg.s_sei.u1_sei_ave_params_present_flag =
            ps_cfg->s_sei.u1_sei_ave_params_present_flag;
        ps_codec->s_cfg.s_sei.s_sei_ave_params = ps_cfg->s_sei.s_sei_ave_params;
    }
    else if(ps_cfg->e_cmd == ISVCE_CMD_CTL_SET_SEI_CCV_PARAMS)
    {
        ps_codec->s_cfg.s_sei.u1_sei_ccv_params_present_flag =
            ps_cfg->s_sei.u1_sei_ccv_params_present_flag;
        ps_codec->s_cfg.s_sei.s_sei_ccv_params = ps_cfg->s_sei.s_sei_ccv_params;
    }

    /* reset RC model */
    if(u4_init_rc)
    {
        for(i = 0; i < ps_cfg->s_svc_params.u1_num_spatial_layers; i++)
        {
            /* init qp */
            UWORD8 au1_init_qp[MAX_PIC_TYPE];

            /* min max qp */
            UWORD8 au1_min_max_qp[2 * MAX_PIC_TYPE];

            /* init i,p,b qp */
            au1_init_qp[0] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp[i]];
            au1_init_qp[1] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp[i]];
            au1_init_qp[2] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp[i]];

            /* init min max qp */
            au1_min_max_qp[2 * I_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp_min[i]];
            au1_min_max_qp[2 * I_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_i_qp_max[i]];

            au1_min_max_qp[2 * P_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp_min[i]];
            au1_min_max_qp[2 * P_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_p_qp_max[i]];

            au1_min_max_qp[2 * B_PIC] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp_min[i]];
            au1_min_max_qp[2 * B_PIC + 1] =
                gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.au4_b_qp_max[i]];

            /* get rc mode */
            switch(ps_codec->s_cfg.e_rc_mode)
            {
                case IVE_RC_STORAGE:
                    ps_codec->s_rate_control.e_rc_type = VBR_STORAGE;
                    break;

                case IVE_RC_CBR_NON_LOW_DELAY:
                    ps_codec->s_rate_control.e_rc_type = CBR_NLDRC;
                    break;

                case IVE_RC_CBR_LOW_DELAY:
                    ps_codec->s_rate_control.e_rc_type = CBR_LDRC;
                    break;

                case IVE_RC_NONE:
                    ps_codec->s_rate_control.e_rc_type = CONST_QP;
                    break;

                default:
                    break;
            }

            /* init rate control */
            for(i = 0; i < MAX_NUM_SPATIAL_LAYERS; i++)
            {
                isvce_rc_init(
                    ps_codec->s_rate_control.apps_rate_control_api[i],
                    ps_codec->s_rate_control.pps_frame_time,
                    ps_codec->s_rate_control.pps_time_stamp,
                    ps_codec->s_rate_control.pps_pd_frm_rate, ps_codec->s_cfg.u4_max_framerate,
                    ps_codec->s_cfg.u4_src_frame_rate, ps_codec->s_cfg.u4_tgt_frame_rate,
                    ps_codec->s_rate_control.e_rc_type, ps_codec->s_cfg.au4_target_bitrate[i],
                    ps_codec->s_cfg.au4_max_bitrate[i], ps_codec->s_cfg.au4_vbv_buffer_delay[i],
                    ps_codec->s_cfg.u4_i_frm_interval, ps_codec->s_cfg.u4_num_bframes + 1,
                    au1_init_qp, ps_codec->s_cfg.u4_num_bframes + 2, au1_min_max_qp,
                    ps_codec->s_cfg.u4_max_level);
            }
        }
    }

    return err;
}

static FORCEINLINE void isvce_change_rc_init_qp(void *pv_rate_control_api, UWORD8 u1_qp)
{
    UWORD8 au1_pic_qps[MAX_PIC_TYPE];
    WORD32 i;

    for(i = 0; i < MAX_PIC_TYPE; i++)
    {
        au1_pic_qps[i] = gau1_h264_to_mpeg2_qmap[CLIP3(MIN_H264_QP, MAX_H264_QP, u1_qp + i)];
    }

    irc_change_init_qp(pv_rate_control_api, au1_pic_qps);
}

/**
 *******************************************************************************
 *
 * @brief
 *  Queues the current buffer, gets back a another buffer for encoding with
 *corrent picture type
 *
 * @par Description:
 *      This function performs 3 distinct but related functions.
 *      1) Maintains an input queue [Note the the term queue donot imply a
 *         first-in first-out logic here] that queues input and dequeues them so
 *         that input frames can be encoded at any predetermined encoding order
 *      2) Uses RC library to decide which frame must be encoded in current pass
 *         and which picture type it must be encoded to.
 *      3) Uses RC library to decide the QP at which current frame has to be
 *         encoded
 *      4) Determines if the current picture must be encoded or not based on
 *         PRE-ENC skip
 *
 *     Input queue is used for storing input buffers till they are used for
 *     encoding. This queue is maintained at ps_codec->as_inp_list. Whenever a
 *     valid input comes, it is added to the end of queue. This same input is
 *     added to RC queue using the identifier as ps_codec->i4_pic_cnt. Hence any
 *     pic from RC can be located in the input queue easily.
 *
 *     The dequeue operation does not start till we have
 *ps_codec->s_cfg.u4_max_num_bframes frames in the queue. THis is done in order
 *to ensure that once output starts we will have a constant stream of output
 *with no gaps.
 *
 *     THe output frame order is governed by RC library. When ever we dequeue a
 *     buffer from RC library, it ensures that we will get them in encoding
 *order With the output of RC library, we can use the picture id to dequeue the
 *     corresponding buffer from input queue and encode it.
 *
 *     Condition at the end of stream.
 *     -------------------------------
 *      At the last valid buffer from the app, we will get ps_ive_ip->u4_is_last
 *      to be set. This will the given to lib when appropriate input buffer is
 *      given to encoding.
 *
 *      Since we have to output is not in sync with input, we will have frames
 *to encode even after we recive the last vaild input buffer. Hence we have to
 *      make sure that we donot queue any new buffers once we get the flag [It
 *may mess up GOP ?]. This is acheived by setting
 *ps_codec->i4_last_inp_buff_received to act as a permenent marker for last
 *frame recived [This may not be needed, because in our current app, all buffers
 *after the last are marked as last. But can we rely on that?] . Hence after
 *this flgag is set no new buffers are queued.
 *
 * @param[in] ps_codec
 *   Pointer to codec descriptor
 *
 * @param[in] ps_ive_ip
 *   Current input buffer to the encoder
 *
 * @param[out] ps_inp
 *   Buffer to be encoded in the current pass
 *
 * @returns
 *   Flag indicating if we have a pre-enc skip or not
 *
 * @remarks
 * TODO (bpic)
 *  The check for null ans is last is redudent.
 *  Need to see if we can remove it
 *
 *******************************************************************************
 */
WORD32 isvce_input_queue_update(isvce_codec_t *ps_codec, ive_video_encode_ip_t *ps_ive_ip,
                                isvce_inp_buf_t *ps_enc_buff, WORD8 i1_layer_id)
{
    isvce_inp_buf_t *ps_inp_buf;
    picture_type_e e_pictype;
    WORD32 i4_skip;
    UWORD32 ctxt_sel, u4_pic_id, u4_pic_disp_id;
    UWORD8 u1_frame_qp = MAX_H264_QP;
    UWORD32 max_frame_bits = 0x7FFFFFFF;

    WORD32 i;

    /*  Mark that the last input frame has been received */
    if(ps_ive_ip->u4_is_last == 1)
    {
        ps_codec->i4_last_inp_buff_received = 1;
    }

    if(ps_ive_ip->s_inp_buf.apv_bufs[0] == NULL && !ps_codec->i4_last_inp_buff_received)
    {
        ps_enc_buff->s_inp_props.s_raw_buf.apv_bufs[0] = NULL;
        ps_enc_buff->s_inp_props.u4_is_last = ps_ive_ip->u4_is_last;
        return 0;
    }

    /***************************************************************************
     * Check for pre enc skip
     *   When src and target frame rates donot match, we skip some frames to
     *   maintain the relation ship between them
     **************************************************************************/
    {
        WORD32 skip_src;

        skip_src = isvce_update_rc_framerates(
            ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id],
            ps_codec->s_rate_control.pps_pd_frm_rate, ps_codec->s_rate_control.pps_time_stamp,
            ps_codec->s_rate_control.pps_frame_time);

        if(skip_src)
        {
            ps_enc_buff->s_inp_props.u4_is_last = ps_ive_ip->u4_is_last;
            return 1;
        }
    }

    /***************************************************************************
     *Queue the input to the queue
     **************************************************************************/
    ps_inp_buf = &(ps_codec->as_inp_list[ps_codec->i4_pic_cnt % SVC_MAX_NUM_INP_FRAMES]);

    /* copy input info. to internal structure */
    ps_inp_buf->s_inp_props.s_raw_buf = ps_ive_ip->s_inp_buf;
    ps_inp_buf->s_inp_props.u4_timestamp_low = ps_ive_ip->u4_timestamp_low;
    ps_inp_buf->s_inp_props.u4_timestamp_high = ps_ive_ip->u4_timestamp_high;
    ps_inp_buf->s_inp_props.u4_is_last = ps_ive_ip->u4_is_last;
    ps_inp_buf->s_inp_props.pv_mb_info = ps_ive_ip->pv_mb_info;
    ps_inp_buf->s_inp_props.u4_mb_info_type = ps_ive_ip->u4_mb_info_type;
    ps_inp_buf->s_inp_props.pv_pic_info = ps_ive_ip->pv_pic_info;
    ps_inp_buf->s_inp_props.u4_pic_info_type = ps_ive_ip->u4_pic_info_type;

    ps_inp_buf->s_inp_props.u1_sei_ccv_params_present_flag =
        ps_codec->s_cfg.s_sei.u1_sei_ccv_params_present_flag;
    ps_inp_buf->s_inp_props.s_sei_ccv = ps_codec->s_cfg.s_sei.s_sei_ccv_params;

    if(ps_inp_buf->s_inp_props.s_raw_buf.apv_bufs[0])
        isvce_svc_inp_buf_populate(ps_codec, ps_inp_buf);

    /***************************************************************************
     * Now we should add the picture to RC stack here
     **************************************************************************/
    /*
     * If an I frame has been requested, ask  RC to force it
     * For IDR requests, we have to ask RC to force I and set IDR by our selves
     * since RC Donot know about IDR. For forcing an IDR at dequeue stage we
     * should record that an IDR has been requested some where. Hence we will
     * store it in the u4_idr_inp_list at a position same as that of input frame
     */
    {
        WORD32 i4_force_idr, i4_force_i;

        i4_force_idr = (ps_codec->force_curr_frame_type == IV_IDR_FRAME);
        i4_force_idr |= !(ps_codec->i4_pic_cnt % ps_codec->s_cfg.u4_idr_frm_interval);

        i4_force_i = (ps_codec->force_curr_frame_type == IV_I_FRAME);

        ps_codec->i4_pending_idr_flag |= i4_force_idr;

        if((ps_codec->i4_pic_cnt > 0) && (i4_force_idr || i4_force_i))
        {
            irc_force_I_frame(ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id]);
        }

        if(i1_layer_id == (ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers - 1))
        {
            ps_codec->force_curr_frame_type = IV_NA_FRAME;
        }
    }

    irc_add_picture_to_stack(ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id],
                             ps_codec->i4_pic_cnt);

    /* Delay */
    if(ps_codec->i4_encode_api_call_cnt < (WORD32) (ps_codec->s_cfg.u4_num_bframes))
    {
        ps_enc_buff->s_inp_props.s_raw_buf.apv_bufs[0] = NULL;
        ps_enc_buff->s_inp_props.u4_is_last = 0;
        return 0;
    }

    /***************************************************************************
     * Get a new pic to encode
     **************************************************************************/
    /* Query the picture_type */
    e_pictype =
        isvce_rc_get_picture_details(ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id],
                                     (WORD32 *) (&u4_pic_id), (WORD32 *) (&u4_pic_disp_id));

    switch(e_pictype)
    {
        case I_PIC:
            ps_codec->pic_type = PIC_I;
            break;
        case P_PIC:
            ps_codec->pic_type = PIC_P;
            break;
        case B_PIC:
            ps_codec->pic_type = PIC_B;
            break;
        default:
            ps_codec->pic_type = PIC_NA;
            ps_enc_buff->s_inp_props.s_raw_buf.apv_bufs[0] = NULL;
            return 0;
    }

    /* Set IDR if it has been requested */
    if(ps_codec->pic_type == PIC_I)
    {
        ps_codec->pic_type = ps_codec->i4_pending_idr_flag ? PIC_IDR : ps_codec->pic_type;
        ps_codec->i4_pending_idr_flag = 0;
    }

    if(ps_codec->s_rate_control.e_rc_type != CONST_QP && ps_codec->u1_enable_init_qp &&
       (u4_pic_id == 0 ||
        irc_is_scenecut(ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id])))
    {
        DOUBLE d_bpp;

        svc_rc_utils_ctxt_t *ps_svc_rc_utils = &ps_codec->s_rate_control.s_rc_utils;

        UWORD32 u4_src_fps = ps_codec->s_cfg.u4_src_frame_rate / 1000;
        UWORD32 u4_wd = ps_inp_buf->as_layer_yuv_buf_props[i1_layer_id].u4_width;
        UWORD32 u4_ht = ps_inp_buf->as_layer_yuv_buf_props[i1_layer_id].u4_height;
        DOUBLE d_gpp =
            isvce_compute_gpp(ps_svc_rc_utils, &ps_inp_buf->as_layer_yuv_buf_props[i1_layer_id]);

        d_bpp = ((DOUBLE) irc_get_vbv_buf_size(
                     ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id]) /
                 10.) /
                ((DOUBLE) (u4_src_fps * u4_wd * u4_ht));

        u1_frame_qp = (UWORD8) irc_get_frame_level_init_qp(
            ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id],
            ps_codec->s_rate_control.e_rc_type, e_pictype, d_bpp, d_gpp);

        isvce_change_rc_init_qp(ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id],
                                u1_frame_qp);

        ps_codec->au4_frame_qp[i1_layer_id] = u1_frame_qp;
    }
    else
    {
        /* Get current frame Qp */
        u1_frame_qp = (UWORD8) irc_get_frame_level_qp(
            ps_codec->s_rate_control.apps_rate_control_api[i1_layer_id], e_pictype, max_frame_bits);
        ps_codec->au4_frame_qp[i1_layer_id] = gau1_mpeg2_to_h264_qmap[u1_frame_qp];
    }

    /*
     * copy the pic id to poc because the display order is assumed to be same
     * as input order
     */
    ps_codec->i4_poc = u4_pic_id;

    /***************************************************************************
     * Now retrieve the correct picture from the queue
     **************************************************************************/

    /* Mark the skip flag   */
    i4_skip = 0;
    ctxt_sel = ps_codec->i4_encode_api_call_cnt % MAX_CTXT_SETS;
    ps_codec->s_rate_control.pre_encode_skip[ctxt_sel] = i4_skip;

    /* Get a buffer to encode */
    ps_inp_buf = &(ps_codec->as_inp_list[u4_pic_id % SVC_MAX_NUM_INP_FRAMES]);

    /* copy dequeued input to output */
    ps_enc_buff[0] = ps_inp_buf[0];

    /* Special case for encoding trailing B frames
     *
     * In encoding streams with B frames it may happen that we have a B frame
     * at the end without a P/I frame after it. Hence when we are dequeing from
     * the RC, it will return the P frame [next in display order but before in
     * encoding order] first. Since the dequeue happens for an invalid frame we
     * will get a frame with null buff and set u4_is_last. Hence lib with return
     * last frame flag at this point and will stop encoding.
     *
     * Since for the last B frame, we does not have the forward ref frame
     * it makes sense to force it into P.
     *
     * To solve this, in case the current frame is P and if the last frame flag
     * is set, we need to see if there is and pending B frames. If there are any,
     * we should just encode that picture as the current P frame and set
     * that B frame as the last frame. Hence the encoder will terminate naturally
     * once that B-frame is encoded after all the in between frames.
     *
     * Since we cannot touch RC stack directly, the option of actually swapping
     * frames in RC is ruled out. We have to modify the as_inp_list to simulate
     * such a behavior by RC. We can do that by
     *  1) Search through as_inp_list to locate the largest u4_timestamp_low less
     *     than current u4_timestamp_low. This will give us the last B frame
     * before the current P frame. Note that this will handle pre encode skip too
     * since queue happens after pre enc skip. 2) Swap the position in
     * as_inp_list. Hence now the last B frame is encoded as P frame. And the new
     * last B frame will have u4_is_last set so that encoder will end naturally
     * once we reached that B frame or any subsequent frame. Also the current GOP
     * will have 1 less B frame Since we are swapping, the poc will also be
     * in-order. 3) In case we have an IPP stream, the result of our search will
     * be an I/P frame which is already encoded. Thus swap and encode will result
     *     in encoding of duplicate frames. Hence to avoid this we will only
     *     have this work around in case of u4_num_bframes > 0.
     *
     *     In case we have forced an I/IDR frame In between this P frame and
     *     the last B frame -> This cannot happen as the current P frame is
     *     supposed to have u4_is_last set. Thus forcing an I/ IDR after this
     *     is illogical.
     *
     *     In cae if we have forced an I such that the frame just before last
     * frame in is I/P -> This case will never arise. Since we have a closed GOP
     * now, once we force an I, the gop gets reset, hence there will be a B
     * between I/P and I/P.
     */
    if(ps_enc_buff->s_inp_props.u4_is_last && (ps_codec->pic_type == PIC_P) &&
       ps_codec->s_cfg.u4_num_bframes)
    {
        WORD32 cntr;
        WORD32 lst_bframe = -1;
        UWORD32 u4_timestamp_low = 0;
        UWORD32 u4_timestamp_high = 0;
        isvce_inp_buf_t *ps_swap_buff, *ps_inp_list;

        ps_inp_list = &ps_codec->as_inp_list[0];

        /* Now search the inp list for highest timestamp */
        for(cntr = 0; cntr < SVC_MAX_NUM_INP_FRAMES; cntr++)
        {
            if(ps_inp_list[cntr].s_inp_props.s_raw_buf.apv_bufs[0] != NULL)
            {
                if((ps_inp_list[cntr].s_inp_props.u4_timestamp_high > u4_timestamp_high) ||
                   (ps_inp_list[cntr].s_inp_props.u4_timestamp_high == u4_timestamp_high &&
                    ps_inp_list[cntr].s_inp_props.u4_timestamp_low > u4_timestamp_low))
                {
                    u4_timestamp_low = ps_inp_list[cntr].s_inp_props.u4_timestamp_low;
                    u4_timestamp_high = ps_inp_list[cntr].s_inp_props.u4_timestamp_high;
                    lst_bframe = cntr;
                }
            }
        }

        if(lst_bframe != -1)
        {
            ps_swap_buff = &(ps_codec->as_inp_list[lst_bframe]);

            /* copy the last B buffer to output */
            *ps_enc_buff = *ps_swap_buff;

            /* Store the current buf into the queue in place of last B buf */
            *ps_swap_buff = *ps_inp_buf;
        }
    }

    if(ps_enc_buff->s_inp_props.u4_is_last)
    {
        ps_codec->pic_type = PIC_NA;
    }

    /* The buffer in the queue is set to NULL to specify that encoding is done for
     * that frame */
    for(i = 0; i < 3; i++)
    {
        ps_inp_buf->s_inp_props.s_raw_buf.apv_bufs[i] = NULL;
    }

    /* Return the buffer status */
    return (0);
}

/**
******************************************************************************
*
* @brief
*  This function joins all the spawned threads after successful completion of
*  their tasks
*
* @par   Description
*
* @param[in] ps_codec
*  pointer to codec context
*
* @returns  none
*
******************************************************************************
*/
void isvce_join_threads(isvce_codec_t *ps_codec)
{
    WORD32 i = 0;
    WORD32 ret = 0;

    /* join spawned threads */
    while(i < ps_codec->i4_proc_thread_cnt)
    {
        if(ps_codec->ai4_process_thread_created[i])
        {
            ret = ithread_join(ps_codec->apv_proc_thread_handle[i], NULL);

            if(ret != 0)
            {
                ASSERT(0);
            }

            ps_codec->ai4_process_thread_created[i] = 0;
            i++;
        }
    }

    ps_codec->i4_proc_thread_cnt = 0;
}

UWORD32 isvce_get_min_outbuf_size(UWORD32 u4_wd, UWORD32 u4_ht, UWORD8 u1_num_spatial_layers)
{
    return MAX((u4_wd * u4_ht * 3), MIN_STREAM_SIZE) * u1_num_spatial_layers;
}