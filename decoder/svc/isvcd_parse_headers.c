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
 *  isvcd_parse_headers.c
 *
 * @brief
 *  Contains High level syntax[above slice] parsing routines
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_set_default_seq_svc_ext()
 *  - isvcd_parse_subset_sps()
 *  - isvcd_dec_ref_base_pic_marking()
 *  - isvcd_parse_nal_unit()
 *  - isvcd_parse_sps()
 *  - isvcd_parse_pps()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <string.h>
#include <assert.h>

#include "ih264_typedefs.h"
#include "ih264_defs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_bitstrm.h"
#include "isvcd_structs.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_defs.h"
#include "ih264d_parse_slice.h"
#include "ih264d_tables.h"
#include "ih264d_utils.h"
#include "ih264d_nal.h"
#include "ih264d_deblocking.h"
#include "ih264d_mem_request.h"
#include "ih264d_debug.h"
#include "ih264_debug.h"
#include "ih264d_error_handler.h"
#include "ih264d_mb_utils.h"
#include "ih264d_sei.h"
#include "ih264d_vui.h"
#include "ih264d_thread_parse_decode.h"
#include "ih264d_thread_compute_bs.h"
#include "ih264d_quant_scaling.h"
#include "ih264d_defs.h"
#include "ivd.h"
#include "ih264d_parse_islice.h"
#include "isvcd_parse_slice.h"
#include "ih264d_process_bslice.h"
#include "ih264d_process_pslice.h"
#include "isvcd_vui.h"

WORD32 ih264d_access_unit_delimiter_rbsp(dec_struct_t *ps_dec);
void ih264d_get_pre_sei_params(dec_struct_t *ps_dec, UWORD8 u1_nal_unit_type);
UWORD32 ih264d_correct_level_idc(UWORD32 u4_level_idc, UWORD32 u4_total_mbs);
WORD32 ih264d_parse_filler_data(dec_struct_t *ps_dec, dec_bit_stream_t *ps_bitstrm);
void ih264d_parse_end_of_stream(dec_struct_t *ps_dec);
WORD32 ih264d_parse_slice_partition(dec_struct_t *ps_dec, dec_bit_stream_t *ps_bitstrm);
/*!
**************************************************************************
* \if Function name : isvcd_set_default_seq_svc_ext \en
dif
*
* \brief
*    Sets the default values for the svc params in the SVC bitstream
*
* \return
**************************************************************************
*/
void isvcd_set_default_seq_svc_ext(dec_subset_seq_params_t *ps_seq_svc_ext)
{
    ps_seq_svc_ext->u1_inter_layer_deblocking_filter_control_present_flag = 0;
    ps_seq_svc_ext->u1_extended_spatial_scalability_idc = 0;
    ps_seq_svc_ext->u1_chroma_phase_x_plus1_flag = 1;
    ps_seq_svc_ext->u1_chroma_phase_y_plus1 = 1;
    ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_x_plus1_flag =
        ps_seq_svc_ext->u1_chroma_phase_x_plus1_flag;
    ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_y_plus1 = ps_seq_svc_ext->u1_chroma_phase_y_plus1;
    ps_seq_svc_ext->i4_seq_scaled_ref_layer_left_offset = 0;
    ps_seq_svc_ext->i4_seq_scaled_ref_layer_top_offset = 0;
    ps_seq_svc_ext->i4_seq_scaled_ref_layer_right_offset = 0;
    ps_seq_svc_ext->i4_seq_scaled_ref_layer_bottom_offset = 0;
    ps_seq_svc_ext->u1_seq_tcoeff_level_prediction_flag =
        ps_seq_svc_ext->u1_adaptive_tcoeff_level_prediction_flag = 0;
    ps_seq_svc_ext->u1_slice_header_restriction_flag = 0;
    ps_seq_svc_ext->u1_svc_vui_parameters_present_flag = 0;
}
/*!
**************************************************************************
* \if Function name : isvcd_parse_subset_sps \en
dif
*
* \brief
*    Decodes Sequence parameter set from the SVC bitstream
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_subset_sps(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_bit_stream_t *ps_bitstrm)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD8 i;
    dec_seq_params_t *ps_seq = NULL;
    dec_svc_seq_params_t *ps_subset_seq = NULL;
    dec_subset_seq_params_t *ps_seq_svc_ext;
    UWORD8 u1_profile_idc, u1_level_idc, u1_seq_parameter_set_id, u1_mb_aff_flag = 0;
    UWORD16 i2_max_frm_num;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD8 u1_frm, uc_constraint_set0_flag, uc_constraint_set1_flag, uc_constraint_set2_flag;
    WORD32 i4_cropped_ht, i4_cropped_wd;
    UWORD32 u4_temp;
    UWORD64 u8_temp;
    UWORD32 u4_pic_height_in_map_units, u4_pic_width_in_mbs;
    UWORD32 u2_pic_wd = 0;
    UWORD32 u2_pic_ht = 0;
    UWORD32 u2_frm_wd_y = 0;
    UWORD32 u2_frm_ht_y = 0;
    UWORD32 u2_frm_wd_uv = 0;
    UWORD32 u2_frm_ht_uv = 0;
    UWORD32 u2_crop_offset_y = 0;
    UWORD32 u2_crop_offset_uv = 0;
    WORD32 ret;
    /* High profile related syntax element */
    WORD32 i4_i;
    /* G050 */
    UWORD8 u1_frame_cropping_flag,
        u1_frame_cropping_rect_left_ofst = 0, u1_frame_cropping_rect_right_ofst = 0,
        u1_frame_cropping_rect_top_ofst = 0, u1_frame_cropping_rect_bottom_ofst = 0;
    /* G050 */
    /*--------------------------------------------------------------------*/
    /* Decode seq_parameter_set_id and profile and level values           */
    /*--------------------------------------------------------------------*/
    SWITCHONTRACE;
    u1_profile_idc = ih264d_get_bits_h264(ps_bitstrm, 8);
    COPYTHECONTEXT("SPS: profile_idc", u1_profile_idc);

    /* G050 */
    uc_constraint_set0_flag = ih264d_get_bit_h264(ps_bitstrm);
    uc_constraint_set1_flag = ih264d_get_bit_h264(ps_bitstrm);
    uc_constraint_set2_flag = ih264d_get_bit_h264(ps_bitstrm);
    UNUSED(uc_constraint_set1_flag);
    UNUSED(uc_constraint_set2_flag);

    /*****************************************************/
    /* Read 5 bits for uc_constraint_set3_flag (1 bit)   */
    /* and reserved_zero_4bits (4 bits) - Sushant        */
    /*****************************************************/
    ih264d_get_bits_h264(ps_bitstrm, 5);
    /* G050 */
    u1_level_idc = (UWORD8) ih264d_get_bits_h264(ps_bitstrm, 8);
    COPYTHECONTEXT("SPS: u4_level_idc", u1_level_idc);

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u4_temp & MASK_ERR_SEQ_SET_ID) return ERROR_INV_SPS_PPS_T;
    u1_seq_parameter_set_id = u4_temp;
    COPYTHECONTEXT("SPS: seq_parameter_set_id", u1_seq_parameter_set_id);

    if(u1_seq_parameter_set_id >= MAX_NUM_SEQ_PARAMS) return ERROR_INV_SPS_PPS_T;

    /*--------------------------------------------------------------------*/
    /* Find an seq param entry in seqparam array of decStruct             */
    /*--------------------------------------------------------------------*/
    ps_subset_seq = ps_svc_lyr_dec->pv_scratch_subset_sps;
    memset(ps_subset_seq, 0, sizeof(dec_svc_seq_params_t));
    ps_seq = ps_dec->pv_scratch_sps_pps;
    memset(ps_seq, 0, sizeof(dec_seq_params_t));

    ps_seq->u1_profile_idc = u1_profile_idc;
    ps_seq->u1_level_idc = u1_level_idc;
    ps_seq->u1_seq_parameter_set_id = u1_seq_parameter_set_id;

    /* subset_seq_sps_will be stored from location 32 : MAX_NUM_SEQ_PARAMS*/
    u1_seq_parameter_set_id += MAX_NUM_SEQ_PARAMS;
    ps_subset_seq->ps_seq = &ps_dec->ps_sps[u1_seq_parameter_set_id];

    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_profile_idc != u1_profile_idc))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }

    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_level_idc != u1_level_idc))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }
    /*******************************************************************/
    /* Initializations for high profile - Sushant                      */
    /*******************************************************************/
    ps_seq->i4_chroma_format_idc = 1;
    ps_seq->i4_bit_depth_luma_minus8 = 0;
    ps_seq->i4_bit_depth_chroma_minus8 = 0;
    ps_seq->i4_qpprime_y_zero_transform_bypass_flag = 0;
    ps_seq->i4_seq_scaling_matrix_present_flag = 0;
    if(u1_profile_idc == HIGH_PROFILE_IDC || u1_profile_idc == SCALABLE_BASELINE_PROFILE_IDC ||
       u1_profile_idc == SCALABLE_HIGH_PROFILE_IDC)
    {
        /* reading chroma_format_idc   */
        ps_seq->i4_chroma_format_idc = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        /* Monochrome is not supported */
        if(ps_seq->i4_chroma_format_idc != 1)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        /* reading bit_depth_luma_minus8   */
        ps_seq->i4_bit_depth_luma_minus8 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_seq->i4_bit_depth_luma_minus8 != 0)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        /* reading bit_depth_chroma_minus8   */
        ps_seq->i4_bit_depth_chroma_minus8 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_seq->i4_bit_depth_chroma_minus8 != 0)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        /* reading qpprime_y_zero_transform_bypass_flag   */
        ps_seq->i4_qpprime_y_zero_transform_bypass_flag = (WORD32) ih264d_get_bit_h264(ps_bitstrm);

        if(ps_seq->i4_qpprime_y_zero_transform_bypass_flag != 0)
        {
            return ERROR_INV_SPS_PPS_T;
        }

        /* reading seq_scaling_matrix_present_flag   */
        ps_seq->i4_seq_scaling_matrix_present_flag = (WORD32) ih264d_get_bit_h264(ps_bitstrm);

        if(ps_seq->i4_seq_scaling_matrix_present_flag)
        {
            for(i4_i = 0; i4_i < 8; i4_i++)
            {
                ps_seq->u1_seq_scaling_list_present_flag[i4_i] = ih264d_get_bit_h264(ps_bitstrm);

                /* initialize u1_use_default_scaling_matrix_flag[i4_i] to zero */
                /* before calling scaling list                             */
                ps_seq->u1_use_default_scaling_matrix_flag[i4_i] = 0;

                if(ps_seq->u1_seq_scaling_list_present_flag[i4_i])
                {
                    if(i4_i < 6)
                    {
                        ret = ih264d_scaling_list(ps_seq->i2_scalinglist4x4[i4_i], 16,
                                                  &ps_seq->u1_use_default_scaling_matrix_flag[i4_i],
                                                  ps_bitstrm);
                    }
                    else
                    {
                        ret = ih264d_scaling_list(ps_seq->i2_scalinglist8x8[i4_i - 6], 64,
                                                  &ps_seq->u1_use_default_scaling_matrix_flag[i4_i],
                                                  ps_bitstrm);
                    }
                    if(ret != OK)
                    {
                        return ret;
                    }
                }
            }
        }
    }
    /*--------------------------------------------------------------------*/
    /* Decode MaxFrameNum                                                 */
    /*--------------------------------------------------------------------*/
    u8_temp = (UWORD64) 4 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u8_temp > MAX_BITS_IN_FRAME_NUM)
    {
        return ERROR_INV_SPS_PPS_T;
    }
    ps_seq->u1_bits_in_frm_num = (UWORD8) u8_temp;
    COPYTHECONTEXT("SPS: log2_max_frame_num_minus4", (ps_seq->u1_bits_in_frm_num - 4));

    i2_max_frm_num = (1 << (ps_seq->u1_bits_in_frm_num));
    ps_seq->u2_u4_max_pic_num_minus1 = i2_max_frm_num - 1;
    /*--------------------------------------------------------------------*/
    /* Decode picture order count and related values                      */
    /*--------------------------------------------------------------------*/
    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(u4_temp > MAX_PIC_ORDER_CNT_TYPE)
    {
        return ERROR_INV_POC_TYPE_T;
    }
    ps_seq->u1_pic_order_cnt_type = u4_temp;
    COPYTHECONTEXT("SPS: pic_order_cnt_type", ps_seq->u1_pic_order_cnt_type);

    ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle = 1;
    if(ps_seq->u1_pic_order_cnt_type == 0)
    {
        u8_temp = (UWORD64) 4 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u8_temp > MAX_BITS_IN_POC_LSB)
        {
            return ERROR_INV_SPS_PPS_T;
        }
        ps_seq->u1_log2_max_pic_order_cnt_lsb_minus = (UWORD8) u8_temp;
        ps_seq->i4_max_pic_order_cntLsb = (1 << u8_temp);
        COPYTHECONTEXT("SPS: log2_max_pic_order_cnt_lsb_minus4", (u8_temp - 4));
    }
    else if(ps_seq->u1_pic_order_cnt_type == 1)
    {
        ps_seq->u1_delta_pic_order_always_zero_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SPS: delta_pic_order_always_zero_flag",
                       ps_seq->u1_delta_pic_order_always_zero_flag);

        ps_seq->i4_ofst_for_non_ref_pic = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: offset_for_non_ref_pic", ps_seq->i4_ofst_for_non_ref_pic);

        ps_seq->i4_ofst_for_top_to_bottom_field = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: offset_for_top_to_bottom_field",
                       ps_seq->i4_ofst_for_top_to_bottom_field);

        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > 255) return ERROR_INV_SPS_PPS_T;
        ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle = u4_temp;
        COPYTHECONTEXT("SPS: num_ref_frames_in_pic_order_cnt_cycle",
                       ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle);

        for(i = 0; i < ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle; i++)
        {
            ps_seq->i4_ofst_for_ref_frame[i] = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("SPS: offset_for_ref_frame", ps_seq->i4_ofst_for_ref_frame[i]);
        }
    }

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if((u4_temp > H264_MAX_REF_PICS))
    {
        return ERROR_NUM_REF;
    }

    /* Compare with older num_ref_frames is header is already once */
    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_num_ref_frames != u4_temp))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }
    ps_seq->u1_num_ref_frames = u4_temp;
    COPYTHECONTEXT("SPS: num_ref_frames", ps_seq->u1_num_ref_frames);

    ps_seq->u1_gaps_in_frame_num_value_allowed_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SPS: gaps_in_frame_num_value_allowed_flag",
                   ps_seq->u1_gaps_in_frame_num_value_allowed_flag);
    /* SVC_DEC_REVIEW */
    ps_seq->u1_gaps_in_frame_num_value_allowed_flag = 0;

    /*--------------------------------------------------------------------*/
    /* Decode FrameWidth and FrameHeight and related values               */
    /*--------------------------------------------------------------------*/
    u8_temp = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    /* Check  for unsupported resolutions*/
    if(u8_temp > (H264_MAX_FRAME_WIDTH >> 4))
    {
        return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
    }
    u4_pic_width_in_mbs = (UWORD32) u8_temp;
    COPYTHECONTEXT("SPS: pic_width_in_mbs_minus1", u4_pic_width_in_mbs - 1);

    u8_temp = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u8_temp > (H264_MAX_FRAME_HEIGHT >> 4))
    {
        return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
    }
    u4_pic_height_in_map_units = (UWORD32) u8_temp;

    ps_seq->u2_frm_wd_in_mbs = u4_pic_width_in_mbs;
    ps_seq->u2_frm_ht_in_mbs = u4_pic_height_in_map_units;

    u2_pic_wd = (u4_pic_width_in_mbs << 4);
    u2_pic_ht = (u4_pic_height_in_map_units << 4);
    if(ps_svc_lyr_dec->pic_width < u2_pic_wd)
    {
        ps_svc_lyr_dec->pic_width = u2_pic_wd;
    }
    if(ps_svc_lyr_dec->pic_height < u2_pic_ht)
    {
        ps_svc_lyr_dec->pic_height = u2_pic_ht;
    }

    /*--------------------------------------------------------------------*/
    /* Get the value of MaxMbAddress and Number of bits needed for it     */
    /*--------------------------------------------------------------------*/
    ps_seq->u4_max_mb_addr = (ps_seq->u2_frm_wd_in_mbs * ps_seq->u2_frm_ht_in_mbs) - 1;

    ps_seq->u4_total_num_of_mbs = ps_seq->u4_max_mb_addr + 1;

    ps_seq->u1_level_idc = ih264d_correct_level_idc(u1_level_idc, ps_seq->u4_total_num_of_mbs);

    u1_frm = ih264d_get_bit_h264(ps_bitstrm);

    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_frame_mbs_only_flag != u1_frm))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }
    ps_seq->u1_frame_mbs_only_flag = u1_frm;

    COPYTHECONTEXT("SPS: frame_mbs_only_flag", u1_frm);

    if(!u1_frm) u1_mb_aff_flag = ih264d_get_bit_h264(ps_bitstrm);
    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_mb_aff_flag != u1_mb_aff_flag))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }
    if(!u1_frm)
    {
        u2_pic_ht <<= 1;
        ps_seq->u1_mb_aff_flag = u1_mb_aff_flag;
        COPYTHECONTEXT("SPS: mb_adaptive_frame_field_flag", ps_seq->u1_mb_aff_flag);
    }
    else
        ps_seq->u1_mb_aff_flag = 0;

    ps_seq->u1_direct_8x8_inference_flag = ih264d_get_bit_h264(ps_bitstrm);

    COPYTHECONTEXT("SPS: direct_8x8_inference_flag", ps_seq->u1_direct_8x8_inference_flag);

    /* G050 */
    u1_frame_cropping_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SPS: frame_cropping_flag", u1_frame_cropping_flag);

    if(u1_frame_cropping_flag)
    {
        u1_frame_cropping_rect_left_ofst = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: frame_cropping_rect_left_offset", u1_frame_cropping_rect_left_ofst);
        u1_frame_cropping_rect_right_ofst = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: frame_cropping_rect_right_offset", u1_frame_cropping_rect_right_ofst);
        u1_frame_cropping_rect_top_ofst = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: frame_cropping_rect_top_offset", u1_frame_cropping_rect_top_ofst);
        u1_frame_cropping_rect_bottom_ofst = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: frame_cropping_rect_bottom_offset",
                       u1_frame_cropping_rect_bottom_ofst);
    }
    /* G050 */
    ps_seq->u1_vui_parameters_present_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SPS: vui_parameters_present_flag", ps_seq->u1_vui_parameters_present_flag);

    u2_frm_wd_y = u2_pic_wd + (UWORD8) (PAD_LEN_Y_H << 1);
    if(1 == ps_dec->u4_share_disp_buf)
    {
        if(ps_dec->u4_app_disp_width > u2_frm_wd_y) u2_frm_wd_y = ps_dec->u4_app_disp_width;
    }

    u2_frm_ht_y = u2_pic_ht + (UWORD8) (PAD_LEN_Y_V << 2);
    u2_frm_wd_uv = u2_pic_wd + (UWORD8) (PAD_LEN_UV_H << 2);
    u2_frm_wd_uv = MAX(u2_frm_wd_uv, u2_frm_wd_y);

    u2_frm_ht_uv = (u2_pic_ht >> 1) + (UWORD8) (PAD_LEN_UV_V << 2);
    u2_frm_ht_uv = MAX(u2_frm_ht_uv, (u2_frm_ht_y >> 1));

    /* Calculate display picture width, height and start u4_ofst from YUV420 */
    /* pictute buffers as per cropping information parsed above             */
    {
        UWORD16 u2_rgt_ofst = 0;
        UWORD16 u2_lft_ofst = 0;
        UWORD16 u2_top_ofst = 0;
        UWORD16 u2_btm_ofst = 0;
        UWORD8 u1_frm_mbs_flag;
        UWORD8 u1_vert_mult_factor;

        if(u1_frame_cropping_flag)
        {
            /* Calculate right and left u4_ofst for cropped picture           */
            u2_rgt_ofst = u1_frame_cropping_rect_right_ofst << 1;
            u2_lft_ofst = u1_frame_cropping_rect_left_ofst << 1;

            /* Know frame MBs only u4_flag                                      */
            u1_frm_mbs_flag = (1 == ps_seq->u1_frame_mbs_only_flag);

            /* Simplify the vertical u4_ofst calculation from field/frame     */
            u1_vert_mult_factor = (2 - u1_frm_mbs_flag);

            /* Calculate bottom and top u4_ofst for cropped  picture          */
            u2_btm_ofst = (u1_frame_cropping_rect_bottom_ofst << u1_vert_mult_factor);
            u2_top_ofst = (u1_frame_cropping_rect_top_ofst << u1_vert_mult_factor);
        }

        /* Calculate u4_ofst from start of YUV 420 picture buffer to start of*/
        /* cropped picture buffer                                           */
        u2_crop_offset_y = (u2_frm_wd_y * u2_top_ofst) + (u2_lft_ofst);
        u2_crop_offset_uv =
            (u2_frm_wd_uv * (u2_top_ofst >> 1)) + (u2_lft_ofst >> 1) * YUV420SP_FACTOR;
        /* Calculate the display picture width and height based on crop      */
        /* information                                                       */
        i4_cropped_ht = (WORD32) u2_pic_ht - (WORD32) (u2_btm_ofst + u2_top_ofst);
        i4_cropped_wd = (WORD32) u2_pic_wd - (WORD32) (u2_rgt_ofst + u2_lft_ofst);

        if((i4_cropped_ht < MB_SIZE) || (i4_cropped_wd < MB_SIZE))
        {
            return ERROR_INV_SPS_PPS_T;
        }

        if((ps_dec->i4_header_decoded & 1) &&
           (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
           (ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].u2_pic_wd != u2_pic_wd))
        {
            ps_dec->u1_res_changed = 1;
            return IVD_RES_CHANGED;
        }

        if((ps_dec->i4_header_decoded & 1) &&
           (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
           (ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].u2_disp_width != i4_cropped_wd))
        {
            ps_dec->u1_res_changed = 1;
            return IVD_RES_CHANGED;
        }

        if((ps_dec->i4_header_decoded & 1) &&
           (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
           (ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].u2_pic_ht != u2_pic_ht))
        {
            ps_dec->u1_res_changed = 1;
            return IVD_RES_CHANGED;
        }

        if((ps_dec->i4_header_decoded & 1) &&
           (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
           (ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].u2_disp_height != i4_cropped_ht))
        {
            ps_dec->u1_res_changed = 1;
            return IVD_RES_CHANGED;
        }
        /* Check again for unsupported resolutions with updated values*/
        if((u2_pic_wd > SVCD_MAX_FRAME_WIDTH) || (u2_pic_ht > SVCD_MAX_FRAME_HEIGHT) ||
           (u2_pic_wd < SVCD_MIN_FRAME_WIDTH) || (u2_pic_ht < SVCD_MIN_FRAME_HEIGHT) ||
           (u2_pic_wd * (UWORD32) u2_pic_ht > SVCD_MAX_FRAME_SIZE))
        {
            return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
        }

        /* If MBAff is enabled, decoder support is limited to streams with
         * width less than half of H264_MAX_FRAME_WIDTH.
         * In case of MBAff decoder processes two rows at a time
         */
        if((u2_pic_wd << ps_seq->u1_mb_aff_flag) > H264_MAX_FRAME_WIDTH)
        {
            return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
        }
    }

    if(1 == ps_seq->u1_vui_parameters_present_flag)
    {
        ret = ih264d_parse_vui_parametres(&ps_seq->s_vui, ps_bitstrm);
        if(ret != OK) return ret;
    }
    ps_seq_svc_ext = &ps_subset_seq->s_sps_svc_ext;

    isvcd_set_default_seq_svc_ext(ps_seq_svc_ext);

    if(SCALABLE_BASELINE_PROFILE_IDC == ps_seq->u1_profile_idc ||
       SCALABLE_HIGH_PROFILE_IDC == ps_seq->u1_profile_idc)
    {
        SWITCHONTRACE;
        ps_seq_svc_ext->u1_inter_layer_deblocking_filter_control_present_flag =
            ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SPS_EXt: u1_inter_layer_deblocking_filter_control_present_flag",
                       ps_seq_svc_ext->u1_inter_layer_deblocking_filter_control_present_flag);

        ps_seq_svc_ext->u1_extended_spatial_scalability_idc = ih264d_get_bits_h264(ps_bitstrm, 2);
        COPYTHECONTEXT("SPS_EXt: u1_extended_spatial_scalability_idc",
                       ps_seq_svc_ext->u1_extended_spatial_scalability_idc);

        /* u1_extended_spatial_scalability_idc value 0, 1 and 2 are supported */
        if(ps_seq_svc_ext->u1_extended_spatial_scalability_idc > 2)
        {
            return ERROR_SVC_INV_SUBSET_SPS;
        }

        /* ChromaArrayType = i4_chroma_format_idc  if  separate_colour_plane_flag =
         * 0 for all chroma format except 4:4:4 */
        if(1 == ps_seq->i4_chroma_format_idc || 2 == ps_seq->i4_chroma_format_idc)
        {
            ps_seq_svc_ext->u1_chroma_phase_x_plus1_flag = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("SPS_EXt: u1_chroma_phase_x_plus1_flag",
                           ps_seq_svc_ext->u1_chroma_phase_x_plus1_flag);
        }

        if(1 == ps_seq->i4_chroma_format_idc)
        {
            ps_seq_svc_ext->u1_chroma_phase_y_plus1 = ih264d_get_bits_h264(ps_bitstrm, 2);
            COPYTHECONTEXT("SPS_EXt: u1_chroma_phase_y_plus1",
                           ps_seq_svc_ext->u1_chroma_phase_y_plus1);

            if(ps_seq_svc_ext->u1_chroma_phase_y_plus1 >= 3)
            {
                return ERROR_SVC_INV_SUBSET_SPS;
            }
        }

        /* inferred values not covered in isvcd_set_default_seq_svc_ext*/
        ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_x_plus1_flag =
            ps_seq_svc_ext->u1_chroma_phase_x_plus1_flag;
        ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_y_plus1 =
            ps_seq_svc_ext->u1_chroma_phase_y_plus1;

        if(1 == ps_seq_svc_ext->u1_extended_spatial_scalability_idc)
        {
            if(ps_seq->i4_chroma_format_idc > 0)
            {
                ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_x_plus1_flag =
                    ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("SPS_EXt: u1_seq_ref_layer_chroma_phase_x_plus1_flag",
                               ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_x_plus1_flag);

                ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_y_plus1 =
                    ih264d_get_bits_h264(ps_bitstrm, 2);
                COPYTHECONTEXT("SPS_EXt: u1_seq_ref_layer_chroma_phase_y_plus1",
                               ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_y_plus1);

                if(ps_seq_svc_ext->u1_seq_ref_layer_chroma_phase_y_plus1 >= 3)
                {
                    return ERROR_SVC_INV_SUBSET_SPS;
                }
            }

            ps_seq_svc_ext->i4_seq_scaled_ref_layer_left_offset =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("SPS_EXt: i4_seq_scaled_ref_layer_left_offset",
                           ps_seq_svc_ext->i4_seq_scaled_ref_layer_left_offset);

            if(ps_seq_svc_ext->i4_seq_scaled_ref_layer_left_offset != 0)
            {
                return ERROR_SVC_INV_SUBSET_SPS;
            }

            if(ps_seq_svc_ext->i4_seq_scaled_ref_layer_left_offset >= MAX_SCLD_REF_LAYER_OFFSET ||
               ps_seq_svc_ext->i4_seq_scaled_ref_layer_left_offset < MIN_SCLD_REF_LAYER_OFFSET)
            {
                return ERROR_SVC_INV_SUBSET_SPS;
            }

            ps_seq_svc_ext->i4_seq_scaled_ref_layer_top_offset =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("SPS_EXt: i4_seq_scaled_ref_layer_top_offset",
                           ps_seq_svc_ext->i4_seq_scaled_ref_layer_top_offset);

            if(ps_seq_svc_ext->i4_seq_scaled_ref_layer_top_offset != 0)
            {
                return ERROR_SVC_INV_SUBSET_SPS;
            }

            if(ps_seq_svc_ext->i4_seq_scaled_ref_layer_top_offset >= MAX_SCLD_REF_LAYER_OFFSET ||
               ps_seq_svc_ext->i4_seq_scaled_ref_layer_top_offset < MIN_SCLD_REF_LAYER_OFFSET)
            {
                return ERROR_SVC_INV_SUBSET_SPS;
            }

            ps_seq_svc_ext->i4_seq_scaled_ref_layer_right_offset =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("SPS_EXt: i4_seq_scaled_ref_layer_right_offset",
                           ps_seq_svc_ext->i4_seq_scaled_ref_layer_right_offset);

            if(ps_seq_svc_ext->i4_seq_scaled_ref_layer_right_offset >= MAX_SCLD_REF_LAYER_OFFSET ||
               ps_seq_svc_ext->i4_seq_scaled_ref_layer_right_offset < MIN_SCLD_REF_LAYER_OFFSET)
            {
                return ERROR_SVC_INV_SUBSET_SPS;
            }

            ps_seq_svc_ext->i4_seq_scaled_ref_layer_bottom_offset =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("SPS_EXt: i4_seq_scaled_ref_layer_bottom_offset",
                           ps_seq_svc_ext->i4_seq_scaled_ref_layer_bottom_offset);

            if(ps_seq_svc_ext->i4_seq_scaled_ref_layer_bottom_offset >= MAX_SCLD_REF_LAYER_OFFSET ||
               ps_seq_svc_ext->i4_seq_scaled_ref_layer_bottom_offset < MIN_SCLD_REF_LAYER_OFFSET)
            {
                return ERROR_INV_SLICE_HDR_T;
            }
        }

        ps_seq_svc_ext->u1_seq_tcoeff_level_prediction_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SPS_EXt: u1_seq_tcoeff_level_prediction_flag",
                       ps_seq_svc_ext->u1_seq_tcoeff_level_prediction_flag);

        if(1 == ps_seq_svc_ext->u1_seq_tcoeff_level_prediction_flag)
        {
            ps_seq_svc_ext->u1_adaptive_tcoeff_level_prediction_flag =
                ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("SPS_EXt: u1_adaptive_tcoeff_level_prediction_flag",
                           ps_seq_svc_ext->u1_adaptive_tcoeff_level_prediction_flag);
        }

        ps_seq_svc_ext->u1_slice_header_restriction_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SPS_EXt: u1_slice_header_restriction_flag",
                       ps_seq_svc_ext->u1_slice_header_restriction_flag);

        ps_seq_svc_ext->u1_svc_vui_parameters_present_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SPS_EXt: u1_svc_vui_parameters_present_flag",
                       ps_seq_svc_ext->u1_svc_vui_parameters_present_flag);

        if(1 == ps_seq_svc_ext->u1_svc_vui_parameters_present_flag)
        {
            if(NULL ==
               ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].s_sps_svc_ext.ps_svc_vui_ext)
            {
                void *pv_buf;
                UWORD32 size;
                /* Memory allocation only if VUI is enabled in a particular subset SPS*/
                size = sizeof(svc_vui_ext_t);
                pv_buf = ps_dec->pf_aligned_alloc(ps_dec->pv_mem_ctxt, 128, size);
                RETURN_IF((NULL == pv_buf), IV_FAIL);
                memset(pv_buf, 0, size);
                ps_seq_svc_ext->ps_svc_vui_ext = pv_buf;
                ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id]
                    .s_sps_svc_ext.ps_svc_vui_ext = pv_buf;
            }
            else
            {
                ps_seq_svc_ext->ps_svc_vui_ext =
                    ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id]
                        .s_sps_svc_ext.ps_svc_vui_ext;
            }
            ret = isvcd_parse_vui_ext_parametres(ps_seq_svc_ext->ps_svc_vui_ext, ps_bitstrm);
            if(ret != OK) return ret;
        }
    }
    /* Add conditions for SCALABLE BASELINE PROFILE */
    if(SCALABLE_BASELINE_PROFILE_IDC == ps_seq->u1_profile_idc ||
       ((SCALABLE_HIGH_PROFILE_IDC == ps_seq->u1_profile_idc) && (1 == uc_constraint_set0_flag)))
    {
        if(ps_seq->i4_chroma_format_idc != 1)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        if(ps_seq->i4_bit_depth_luma_minus8 != 0)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        if(ps_seq->i4_bit_depth_chroma_minus8 != 0)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        if(ps_seq->i4_qpprime_y_zero_transform_bypass_flag != 0)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        if(ps_seq->u1_frame_mbs_only_flag != 1)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        if((0 != ps_seq_svc_ext->i4_seq_scaled_ref_layer_left_offset % 16) &&
           (0 != ps_seq_svc_ext->i4_seq_scaled_ref_layer_top_offset % 16))
        {
            return ERROR_FEATURE_UNAVAIL;
        }
    }
    /* Compare older num_reorder_frames with the new one if header is already
     * decoded */
    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_vui_parameters_present_flag) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].s_vui.u1_bitstream_restriction_flag))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }
    /* In case bitstream read has exceeded the filled size, then return an error */
    if(EXCEED_OFFSET(ps_bitstrm))
    {
        return ERROR_INV_SPS_PPS_T;
    }

    /*--------------------------------------------------------------------*/
    /* All initializations to ps_dec are beyond this point                */
    /*--------------------------------------------------------------------*/
    {
        WORD32 reorder_depth = ih264d_get_dpb_size(ps_seq);
        if((1 == ps_seq->u1_vui_parameters_present_flag) &&
           (1 == ps_seq->s_vui.u1_bitstream_restriction_flag))
        {
            reorder_depth = ps_seq->s_vui.u4_num_reorder_frames + 1;
        }

        if(reorder_depth > H264_MAX_REF_PICS)
        {
            return ERROR_INV_SPS_PPS_T;
        }

        if(ps_seq->u1_frame_mbs_only_flag != 1) reorder_depth *= 2;
        ps_subset_seq->i4_reorder_depth = reorder_depth + DISPLAY_LATENCY;
    }
    ps_subset_seq->u2_disp_height = i4_cropped_ht;
    ps_subset_seq->u2_disp_width = i4_cropped_wd;
    ps_subset_seq->u2_pic_wd = u2_pic_wd;
    ps_subset_seq->u2_pic_ht = u2_pic_ht;

    /* Assuming 8k is the maximum resolution svc dec supports*/
    if(u2_frm_wd_y > H264_MAX_FRAME_WIDTH) return (NOT_OK);
    if(u2_frm_ht_y > H264_MAX_FRAME_HEIGHT) return (NOT_OK);
    if(u2_frm_wd_uv > H264_MAX_FRAME_WIDTH) return (NOT_OK);
    if(u2_frm_ht_uv > H264_MAX_FRAME_HEIGHT) return (NOT_OK);

    /* Determining the Width and Height of Frame from that of Picture */
    ps_subset_seq->u2_frm_wd_y = u2_frm_wd_y;
    ps_subset_seq->u2_frm_ht_y = u2_frm_ht_y;
    ps_subset_seq->u2_frm_wd_uv = u2_frm_wd_uv;
    ps_subset_seq->u2_frm_ht_uv = u2_frm_ht_uv;

    ps_subset_seq->u1_pad_len_y_v = (UWORD8) (PAD_LEN_Y_V << (1 - u1_frm));
    ps_subset_seq->u1_pad_len_cr_v = (UWORD8) (PAD_LEN_UV_V << (1 - u1_frm));

    ps_subset_seq->u2_crop_offset_y = u2_crop_offset_y;
    ps_subset_seq->u2_crop_offset_uv = u2_crop_offset_uv;

    ps_seq->u1_is_valid = TRUE;
    ps_dec->ps_sps[u1_seq_parameter_set_id] = *ps_seq;
    if(NULL != ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].s_sps_svc_ext.ps_svc_vui_ext)
    {
        ps_seq_svc_ext->ps_svc_vui_ext =
            ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].s_sps_svc_ext.ps_svc_vui_ext;
    }
    ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id] = *ps_subset_seq;
    ps_svc_lyr_dec->ps_cur_subset_sps = &ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id];

    return OK;
}
/*!
 **************************************************************************
 * \if Function name : isvcd_dec_ref_base_pic_marking \endif
 *
 * \brief
 *    Decodes reference base pic marking params
 *
 * \return
 *    0 on Success and error code otherwise
 **************************************************************************
 */

WORD32 isvcd_dec_ref_base_pic_marking(
    dec_ref_base_pic_marking_params_t *ps_ref_base_pic_marking_svc_ext,
    dec_bit_stream_t *ps_bitstrm)
{
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;

    SWITCHONTRACE;

    ps_ref_base_pic_marking_svc_ext->u1_adaptive_ref_base_pic_marking_mode_flag =
        ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT(
        "Dec ref base pic marking params : "
        "u1_adaptive_ref_base_pic_marking_mode_flag",
        ps_ref_base_pic_marking_svc_ext->u1_adaptive_ref_base_pic_marking_mode_flag);

    if(1 == ps_ref_base_pic_marking_svc_ext->u1_adaptive_ref_base_pic_marking_mode_flag)
    {
        do
        {
            ps_ref_base_pic_marking_svc_ext->u4_memory_management_base_control_operation =
                ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT(
                "Dec ref base pic marking params : "
                "u4_memory_management_base_control_operation",
                ps_ref_base_pic_marking_svc_ext->u4_memory_management_base_control_operation);

            if(1 == ps_ref_base_pic_marking_svc_ext->u4_memory_management_base_control_operation)
            {
                ps_ref_base_pic_marking_svc_ext->u4_difference_of_base_pic_nums_minus1 =
                    ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                COPYTHECONTEXT(
                    "Dec ref base pic marking params : "
                    "u4_difference_of_base_pic_nums_minus1",
                    ps_ref_base_pic_marking_svc_ext->u4_difference_of_base_pic_nums_minus1);
            }

            if(2 == ps_ref_base_pic_marking_svc_ext->u4_memory_management_base_control_operation)
            {
                ps_ref_base_pic_marking_svc_ext->u4_long_term_base_pic_num =
                    ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                COPYTHECONTEXT("Dec ref base pic marking params : u4_long_term_base_pic_num",
                               ps_ref_base_pic_marking_svc_ext->u4_long_term_base_pic_num);
            }

        } while(0 != ps_ref_base_pic_marking_svc_ext->u4_memory_management_base_control_operation);
    }
    SWITCHOFFTRACE;

    return OK;
}

/*!
 **************************************************************************
 * \if Function name : isvcd_parse_nal_unit \endif
 *
 * \brief
 *    Decodes NAL unit
 *
 * \return
 *    0 on Success and error code otherwise
 **************************************************************************
 */

WORD32 isvcd_parse_nal_unit(svc_dec_lyr_struct_t *dec_svc_hdl, UWORD8 u1_nal_ref_idc)
{
    dec_bit_stream_t *ps_bitstrm;

    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    UWORD8 u1_nal_unit_type;
    WORD32 i_status = OK;

    ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) dec_svc_hdl;
    ps_dec = &ps_svc_lyr_dec->s_dec;

    {
        SWITCHOFFTRACE;
        u1_nal_unit_type = ps_dec->u1_nal_unit_type;

        ps_bitstrm = ps_dec->ps_bitstrm;

        // Skip all NALUs if SPS and PPS are not decoded
        switch(u1_nal_unit_type)
        {
            case SLICE_DATA_PARTITION_A_NAL:
            case SLICE_DATA_PARTITION_B_NAL:
            case SLICE_DATA_PARTITION_C_NAL:
                if(!ps_dec->i4_decode_header) ih264d_parse_slice_partition(ps_dec, ps_bitstrm);
                break;

            case IDR_SLICE_NAL:
            case SLICE_NAL:

                if(ps_svc_lyr_dec->u1_base_res_flag != 1)
                {
                    return NOT_OK;
                }
                if(!ps_dec->i4_decode_header)
                {
                    if(ps_dec->i4_header_decoded == 3)
                    {
                        /* ! */
                        DEBUG_THREADS_PRINTF("Decoding  a slice NAL\n");
                        {
                            ih264d_get_pre_sei_params(ps_dec, u1_nal_unit_type);
                            /* ! */
                            ps_dec->u4_slice_start_code_found = 1;

                            i_status = isvcd_parse_decode_slice(
                                (UWORD8) (u1_nal_unit_type == IDR_SLICE_NAL), u1_nal_ref_idc,
                                ps_svc_lyr_dec);

                            if(i_status != OK)
                            {
                                return i_status;
                            }
                        }
                    }
                }
                break;

            case SEI_NAL:
            case PREFIX_UNIT_NAL:
            case SEQ_PARAM_NAL:
            case PIC_PARAM_NAL:
            case SUBSET_SPS_NAL:
                H264_DEC_DEBUG_PRINT("\nUnknown NAL type %d\n", u1_nal_unit_type);
                break;

            case ACCESS_UNIT_DELIMITER_RBSP:
                if(!ps_dec->i4_decode_header)
                {
                    ih264d_access_unit_delimiter_rbsp(ps_dec);
                }
                break;
                // ignore the END_OF_SEQ_RBSP NAL and decode even after this NAL
            case END_OF_STREAM_RBSP:
                if(!ps_dec->i4_decode_header)
                {
                    ih264d_parse_end_of_stream(ps_dec);
                }
                break;
            case FILLER_DATA_NAL:
                if(!ps_dec->i4_decode_header)
                {
                    ih264d_parse_filler_data(ps_dec, ps_bitstrm);
                }
                break;
            case CODED_SLICE_EXTENSION_NAL:

                if(ps_svc_lyr_dec->u1_base_res_flag == 1)
                {
                    return NOT_OK;
                }
                if(!ps_dec->i4_decode_header)
                {
                    if(ps_dec->i4_header_decoded == 3)
                    {
                        /* ! */
                        DEBUG_THREADS_PRINTF("Decoding  an SVC slice NAL\n");
                        {
                            {
                                ih264d_get_pre_sei_params(ps_dec, u1_nal_unit_type);
                                /* ! */
                                ps_dec->u4_slice_start_code_found = 1;

                                i_status = isvcd_parse_decode_slice_ext_nal(
                                    (UWORD8) (ps_svc_lyr_dec->ps_nal_svc_ext->u1_idr_flag),
                                    u1_nal_ref_idc, ps_svc_lyr_dec);

                                if(i_status != OK)
                                {
                                    return i_status;
                                }
                            }
                        }
                    }
                }
                break;

            default:
                H264_DEC_DEBUG_PRINT("\nUnknown NAL type %d\n", u1_nal_unit_type);
                break;
        }
    }
    return i_status;
}

/*!
**************************************************************************
* \if Function name : isvcd_parse_sps \endif
*
* \brief
*    Decodes Picture Parameter set
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_sps(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_bit_stream_t *ps_bitstrm)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD8 i;
    dec_seq_params_t *ps_seq = NULL;
    dec_svc_seq_params_t *ps_subset_seq = NULL;
    UWORD8 u1_profile_idc, u1_level_idc, u1_seq_parameter_set_id, u1_mb_aff_flag = 0;
    UWORD16 i2_max_frm_num;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD8 u1_frm, uc_constraint_set0_flag, uc_constraint_set1_flag, uc_constraint_set2_flag;
    WORD32 i4_cropped_ht, i4_cropped_wd;
    UWORD32 u4_temp;
    UWORD64 u8_temp;
    UWORD32 u4_pic_height_in_map_units, u4_pic_width_in_mbs;
    UWORD32 u2_pic_wd = 0;
    UWORD32 u2_pic_ht = 0;
    UWORD32 u2_frm_wd_y = 0;
    UWORD32 u2_frm_ht_y = 0;
    UWORD32 u2_frm_wd_uv = 0;
    UWORD32 u2_frm_ht_uv = 0;
    UWORD32 u2_crop_offset_y = 0;
    UWORD32 u2_crop_offset_uv = 0;
    WORD32 ret;
    WORD32 num_reorder_frames;
    /* High profile related syntax element */
    WORD32 i4_i;
    /* G050 */
    UWORD8 u1_frame_cropping_flag,
        u1_frame_cropping_rect_left_ofst = 0, u1_frame_cropping_rect_right_ofst = 0,
        u1_frame_cropping_rect_top_ofst = 0, u1_frame_cropping_rect_bottom_ofst = 0;
    /* G050 */
    /*--------------------------------------------------------------------*/
    /* Decode seq_parameter_set_id and profile and level values           */
    /*--------------------------------------------------------------------*/
    SWITCHONTRACE;
    u1_profile_idc = ih264d_get_bits_h264(ps_bitstrm, 8);
    COPYTHECONTEXT("SPS: profile_idc", u1_profile_idc);

    /* G050 */
    uc_constraint_set0_flag = ih264d_get_bit_h264(ps_bitstrm);
    uc_constraint_set1_flag = ih264d_get_bit_h264(ps_bitstrm);
    uc_constraint_set2_flag = ih264d_get_bit_h264(ps_bitstrm);
    UNUSED(uc_constraint_set2_flag);
    /*****************************************************/
    /* Read 5 bits for uc_constraint_set3_flag (1 bit)   */
    /* and reserved_zero_4bits (4 bits) - Sushant        */
    /*****************************************************/
    ih264d_get_bits_h264(ps_bitstrm, 5);
    /* G050 */
    /* Check whether particular profile is suported or not */
    /* Check whether particular profile is suported or not */
    if((u1_profile_idc != MAIN_PROFILE_IDC) && (u1_profile_idc != BASE_PROFILE_IDC) &&
       (u1_profile_idc != HIGH_PROFILE_IDC))
    {
        /* Apart from Baseline, main and high profile,
         * only extended profile is supported provided
         * uc_constraint_set0_flag or uc_constraint_set1_flag are set to 1
         */
        if((u1_profile_idc != EXTENDED_PROFILE_IDC) ||
           ((uc_constraint_set1_flag != 1) && (uc_constraint_set0_flag != 1)))
        {
            return (ERROR_FEATURE_UNAVAIL);
        }
    }

    u1_level_idc = ih264d_get_bits_h264(ps_bitstrm, 8);
    COPYTHECONTEXT("SPS: u4_level_idc", u1_level_idc);

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u4_temp & MASK_ERR_SEQ_SET_ID) return ERROR_INV_SPS_PPS_T;
    u1_seq_parameter_set_id = u4_temp;
    COPYTHECONTEXT("SPS: seq_parameter_set_id", u1_seq_parameter_set_id);

    /*--------------------------------------------------------------------*/
    /* Find an seq param entry in seqparam array of decStruct             */
    /*--------------------------------------------------------------------*/
    ps_subset_seq = ps_svc_lyr_dec->pv_scratch_subset_sps;
    memset(ps_subset_seq, 0, sizeof(dec_svc_seq_params_t));
    ps_seq = ps_dec->pv_scratch_sps_pps;
    memset(ps_seq, 0, sizeof(dec_seq_params_t));

    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_profile_idc != u1_profile_idc))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }

    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_level_idc != u1_level_idc))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }

    ps_seq->u1_profile_idc = u1_profile_idc;
    ps_seq->u1_level_idc = u1_level_idc;
    ps_seq->u1_seq_parameter_set_id = u1_seq_parameter_set_id;
    ps_subset_seq->ps_seq = &ps_dec->ps_sps[u1_seq_parameter_set_id];

    /*******************************************************************/
    /* Initializations for high profile - Sushant                      */
    /*******************************************************************/
    ps_seq->i4_chroma_format_idc = 1;
    ps_seq->i4_bit_depth_luma_minus8 = 0;
    ps_seq->i4_bit_depth_chroma_minus8 = 0;
    ps_seq->i4_qpprime_y_zero_transform_bypass_flag = 0;
    ps_seq->i4_seq_scaling_matrix_present_flag = 0;
    if(u1_profile_idc == HIGH_PROFILE_IDC || u1_profile_idc == SCALABLE_BASELINE_PROFILE_IDC ||
       u1_profile_idc == SCALABLE_HIGH_PROFILE_IDC)
    {
        /* reading chroma_format_idc   */
        ps_seq->i4_chroma_format_idc = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        /* Monochrome is not supported */
        if(ps_seq->i4_chroma_format_idc != 1)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        /* reading bit_depth_luma_minus8   */
        ps_seq->i4_bit_depth_luma_minus8 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_seq->i4_bit_depth_luma_minus8 != 0)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        /* reading bit_depth_chroma_minus8   */
        ps_seq->i4_bit_depth_chroma_minus8 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_seq->i4_bit_depth_chroma_minus8 != 0)
        {
            return ERROR_FEATURE_UNAVAIL;
        }

        /* reading qpprime_y_zero_transform_bypass_flag   */
        ps_seq->i4_qpprime_y_zero_transform_bypass_flag = (WORD32) ih264d_get_bit_h264(ps_bitstrm);

        if(ps_seq->i4_qpprime_y_zero_transform_bypass_flag != 0)
        {
            return ERROR_INV_SPS_PPS_T;
        }

        /* reading seq_scaling_matrix_present_flag   */
        ps_seq->i4_seq_scaling_matrix_present_flag = (WORD32) ih264d_get_bit_h264(ps_bitstrm);

        if(ps_seq->i4_seq_scaling_matrix_present_flag)
        {
            for(i4_i = 0; i4_i < 8; i4_i++)
            {
                ps_seq->u1_seq_scaling_list_present_flag[i4_i] = ih264d_get_bit_h264(ps_bitstrm);

                /* initialize u1_use_default_scaling_matrix_flag[i4_i] to zero */
                /* before calling scaling list                             */
                ps_seq->u1_use_default_scaling_matrix_flag[i4_i] = 0;

                if(ps_seq->u1_seq_scaling_list_present_flag[i4_i])
                {
                    if(i4_i < 6)
                    {
                        ret = ih264d_scaling_list(ps_seq->i2_scalinglist4x4[i4_i], 16,
                                                  &ps_seq->u1_use_default_scaling_matrix_flag[i4_i],
                                                  ps_bitstrm);
                    }
                    else
                    {
                        ret = ih264d_scaling_list(ps_seq->i2_scalinglist8x8[i4_i - 6], 64,
                                                  &ps_seq->u1_use_default_scaling_matrix_flag[i4_i],
                                                  ps_bitstrm);
                    }
                    if(ret != OK)
                    {
                        return ret;
                    }
                }
            }
        }
    }
    /*--------------------------------------------------------------------*/
    /* Decode MaxFrameNum                                                 */
    /*--------------------------------------------------------------------*/
    u8_temp = (UWORD64) 4 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u8_temp > MAX_BITS_IN_FRAME_NUM)
    {
        return ERROR_INV_SPS_PPS_T;
    }
    ps_seq->u1_bits_in_frm_num = (UWORD8) u8_temp;
    COPYTHECONTEXT("SPS: log2_max_frame_num_minus4", (ps_seq->u1_bits_in_frm_num - 4));

    i2_max_frm_num = (1 << (ps_seq->u1_bits_in_frm_num));
    ps_seq->u2_u4_max_pic_num_minus1 = i2_max_frm_num - 1;
    /*--------------------------------------------------------------------*/
    /* Decode picture order count and related values                      */
    /*--------------------------------------------------------------------*/
    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u4_temp > MAX_PIC_ORDER_CNT_TYPE)
    {
        return ERROR_INV_POC_TYPE_T;
    }
    ps_seq->u1_pic_order_cnt_type = u4_temp;
    COPYTHECONTEXT("SPS: pic_order_cnt_type", ps_seq->u1_pic_order_cnt_type);

    ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle = 1;
    if(ps_seq->u1_pic_order_cnt_type == 0)
    {
        u8_temp = (UWORD64) 4 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u8_temp > MAX_BITS_IN_POC_LSB)
        {
            return ERROR_INV_SPS_PPS_T;
        }
        ps_seq->u1_log2_max_pic_order_cnt_lsb_minus = (UWORD8) u8_temp;
        ps_seq->i4_max_pic_order_cntLsb = (1 << u8_temp);
        COPYTHECONTEXT("SPS: log2_max_pic_order_cnt_lsb_minus4", (u8_temp - 4));
    }
    else if(ps_seq->u1_pic_order_cnt_type == 1)
    {
        ps_seq->u1_delta_pic_order_always_zero_flag = ih264d_get_bit_h264(ps_bitstrm);
        COPYTHECONTEXT("SPS: delta_pic_order_always_zero_flag",
                       ps_seq->u1_delta_pic_order_always_zero_flag);

        ps_seq->i4_ofst_for_non_ref_pic = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: offset_for_non_ref_pic", ps_seq->i4_ofst_for_non_ref_pic);

        ps_seq->i4_ofst_for_top_to_bottom_field = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: offset_for_top_to_bottom_field",
                       ps_seq->i4_ofst_for_top_to_bottom_field);

        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > 255) return ERROR_INV_SPS_PPS_T;
        ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle = u4_temp;
        COPYTHECONTEXT("SPS: num_ref_frames_in_pic_order_cnt_cycle",
                       ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle);

        for(i = 0; i < ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle; i++)
        {
            ps_seq->i4_ofst_for_ref_frame[i] = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            COPYTHECONTEXT("SPS: offset_for_ref_frame", ps_seq->i4_ofst_for_ref_frame[i]);
        }
    }

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if((u4_temp > H264_MAX_REF_PICS))
    {
        return ERROR_NUM_REF;
    }

    /* Compare with older num_ref_frames is header is already once */
    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_num_ref_frames != u4_temp))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }

    ps_seq->u1_num_ref_frames = u4_temp;
    COPYTHECONTEXT("SPS: num_ref_frames", ps_seq->u1_num_ref_frames);

    ps_seq->u1_gaps_in_frame_num_value_allowed_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SPS: gaps_in_frame_num_value_allowed_flag",
                   ps_seq->u1_gaps_in_frame_num_value_allowed_flag);

    ps_seq->u1_gaps_in_frame_num_value_allowed_flag = 0;

    /*--------------------------------------------------------------------*/
    /* Decode FrameWidth and FrameHeight and related values               */
    /*--------------------------------------------------------------------*/
    u8_temp = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    /* Check  for unsupported resolutions*/
    if(u8_temp > (H264_MAX_FRAME_WIDTH >> 4))
    {
        return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
    }
    u4_pic_width_in_mbs = (UWORD32) u8_temp;
    COPYTHECONTEXT("SPS: pic_width_in_mbs_minus1", u4_pic_width_in_mbs - 1);

    u8_temp = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u8_temp > (H264_MAX_FRAME_HEIGHT >> 4))
    {
        return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
    }
    u4_pic_height_in_map_units = (UWORD32) u8_temp;

    ps_seq->u2_frm_wd_in_mbs = u4_pic_width_in_mbs;
    ps_seq->u2_frm_ht_in_mbs = u4_pic_height_in_map_units;
    u2_pic_wd = (u4_pic_width_in_mbs << 4);
    u2_pic_ht = (u4_pic_height_in_map_units << 4);
    if(ps_svc_lyr_dec->pic_width < u2_pic_wd)
    {
        ps_svc_lyr_dec->pic_width = u2_pic_wd;
    }
    if(ps_svc_lyr_dec->pic_height < u2_pic_ht)
    {
        ps_svc_lyr_dec->pic_height = u2_pic_ht;
    }

    /*--------------------------------------------------------------------*/
    /* Get the value of MaxMbAddress and Number of bits needed for it     */
    /*--------------------------------------------------------------------*/
    ps_seq->u4_max_mb_addr = ((UWORD32)ps_seq->u2_frm_wd_in_mbs * (UWORD32)ps_seq->u2_frm_ht_in_mbs) - 1;
    ps_seq->u4_total_num_of_mbs = ps_seq->u4_max_mb_addr + 1;
    ps_seq->u1_level_idc = ih264d_correct_level_idc(u1_level_idc, ps_seq->u4_total_num_of_mbs);

    u1_frm = ih264d_get_bit_h264(ps_bitstrm);
    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_frame_mbs_only_flag != u1_frm))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }
    ps_seq->u1_frame_mbs_only_flag = u1_frm;
    COPYTHECONTEXT("SPS: frame_mbs_only_flag", u1_frm);

    if(!u1_frm) u1_mb_aff_flag = ih264d_get_bit_h264(ps_bitstrm);

    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
       (ps_dec->ps_sps[u1_seq_parameter_set_id].u1_mb_aff_flag != u1_mb_aff_flag))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }

    if(!u1_frm)
    {
        u2_pic_ht <<= 1;
        ps_seq->u1_mb_aff_flag = u1_mb_aff_flag;
        COPYTHECONTEXT("SPS: mb_adaptive_frame_field_flag", ps_seq->u1_mb_aff_flag);
    }
    else
        ps_seq->u1_mb_aff_flag = 0;

    ps_seq->u1_direct_8x8_inference_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SPS: direct_8x8_inference_flag", ps_seq->u1_direct_8x8_inference_flag);

    /* G050 */
    u1_frame_cropping_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SPS: frame_cropping_flag", u1_frame_cropping_flag);

    if(u1_frame_cropping_flag)
    {
        u1_frame_cropping_rect_left_ofst = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: frame_cropping_rect_left_offset", u1_frame_cropping_rect_left_ofst);
        u1_frame_cropping_rect_right_ofst = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: frame_cropping_rect_right_offset", u1_frame_cropping_rect_right_ofst);
        u1_frame_cropping_rect_top_ofst = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: frame_cropping_rect_top_offset", u1_frame_cropping_rect_top_ofst);
        u1_frame_cropping_rect_bottom_ofst = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SPS: frame_cropping_rect_bottom_offset",
                       u1_frame_cropping_rect_bottom_ofst);
    }
    /* G050 */
    ps_seq->u1_vui_parameters_present_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("SPS: vui_parameters_present_flag", ps_seq->u1_vui_parameters_present_flag);

    u2_frm_wd_y = u2_pic_wd + (UWORD8) (PAD_LEN_Y_H << 1);

    if(1 == ps_dec->u4_share_disp_buf)
    {
        if(ps_dec->u4_app_disp_width > u2_frm_wd_y) u2_frm_wd_y = ps_dec->u4_app_disp_width;
    }

    u2_frm_ht_y = u2_pic_ht + (UWORD8) (PAD_LEN_Y_V << 2);
    u2_frm_wd_uv = u2_pic_wd + (UWORD8) (PAD_LEN_UV_H << 2);
    u2_frm_wd_uv = MAX(u2_frm_wd_uv, u2_frm_wd_y);
    u2_frm_ht_uv = (u2_pic_ht >> 1) + (UWORD8) (PAD_LEN_UV_V << 2);
    u2_frm_ht_uv = MAX(u2_frm_ht_uv, (u2_frm_ht_y >> 1));

    /* Calculate display picture width, height and start u4_ofst from YUV420 */
    /* pictute buffers as per cropping information parsed above             */
    {
        UWORD16 u2_rgt_ofst = 0;
        UWORD16 u2_lft_ofst = 0;
        UWORD16 u2_top_ofst = 0;
        UWORD16 u2_btm_ofst = 0;
        UWORD8 u1_frm_mbs_flag;
        UWORD8 u1_vert_mult_factor;

        if(u1_frame_cropping_flag)
        {
            /* Calculate right and left u4_ofst for cropped picture           */
            u2_rgt_ofst = u1_frame_cropping_rect_right_ofst << 1;
            u2_lft_ofst = u1_frame_cropping_rect_left_ofst << 1;

            /* Know frame MBs only u4_flag                                      */
            u1_frm_mbs_flag = (1 == ps_seq->u1_frame_mbs_only_flag);

            /* Simplify the vertical u4_ofst calculation from field/frame     */
            u1_vert_mult_factor = (2 - u1_frm_mbs_flag);

            /* Calculate bottom and top u4_ofst for cropped  picture          */
            u2_btm_ofst = (u1_frame_cropping_rect_bottom_ofst << u1_vert_mult_factor);
            u2_top_ofst = (u1_frame_cropping_rect_top_ofst << u1_vert_mult_factor);
        }

        /* Calculate u4_ofst from start of YUV 420 picture buffer to start of*/
        /* cropped picture buffer                                           */
        u2_crop_offset_y = (u2_frm_wd_y * u2_top_ofst) + (u2_lft_ofst);
        u2_crop_offset_uv =
            (u2_frm_wd_uv * (u2_top_ofst >> 1)) + (u2_lft_ofst >> 1) * YUV420SP_FACTOR;
        /* Calculate the display picture width and height based on crop      */
        /* information                                                       */
        i4_cropped_ht = (WORD32) u2_pic_ht - (WORD32) (u2_btm_ofst + u2_top_ofst);
        i4_cropped_wd = (WORD32) u2_pic_wd - (WORD32) (u2_rgt_ofst + u2_lft_ofst);

        if((i4_cropped_ht < MB_SIZE) || (i4_cropped_wd < MB_SIZE))
        {
            return ERROR_INV_SPS_PPS_T;
        }

        if((ps_dec->i4_header_decoded & 1) &&
           (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
           (ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].u2_pic_wd != u2_pic_wd))
        {
            ps_dec->u1_res_changed = 1;
            return IVD_RES_CHANGED;
        }

        if((ps_dec->i4_header_decoded & 1) &&
           (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
           (ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].u2_disp_width != i4_cropped_wd))
        {
            ps_dec->u1_res_changed = 1;
            return IVD_RES_CHANGED;
        }

        if((ps_dec->i4_header_decoded & 1) &&
           (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
           (ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].u2_pic_ht != u2_pic_ht))
        {
            ps_dec->u1_res_changed = 1;
            return IVD_RES_CHANGED;
        }

        if((ps_dec->i4_header_decoded & 1) &&
           (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) &&
           (ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id].u2_disp_height != i4_cropped_ht))
        {
            ps_dec->u1_res_changed = 1;
            return IVD_RES_CHANGED;
        }
        /* Check again for unsupported resolutions with updated values*/
        if((u2_pic_wd > SVCD_MAX_FRAME_WIDTH) || (u2_pic_ht > SVCD_MAX_FRAME_HEIGHT) ||
           (u2_pic_wd < SVCD_MIN_FRAME_WIDTH) || (u2_pic_ht < SVCD_MIN_FRAME_HEIGHT) ||
           (u2_pic_wd * (UWORD32) u2_pic_ht > SVCD_MAX_FRAME_SIZE))
        {
            return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
        }

        /* If MBAff is enabled, decoder support is limited to streams with
         * width less than half of H264_MAX_FRAME_WIDTH.
         * In case of MBAff decoder processes two rows at a time
         */
        if((u2_pic_wd << ps_seq->u1_mb_aff_flag) > H264_MAX_FRAME_WIDTH)
        {
            return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
        }
    }

    /* Backup num_reorder_frames if header is already decoded */
    if((ps_dec->i4_header_decoded & 1) && (1 == ps_seq->u1_vui_parameters_present_flag) &&
       (1 == ps_seq->s_vui.u1_bitstream_restriction_flag))
    {
        num_reorder_frames = (WORD32) ps_seq->s_vui.u4_num_reorder_frames;
    }
    else
    {
        num_reorder_frames = -1;
    }
    if(1 == ps_seq->u1_vui_parameters_present_flag)
    {
        ret = ih264d_parse_vui_parametres(&ps_seq->s_vui, ps_bitstrm);
        if(ret != OK) return ret;
    }

    /* Compare older num_reorder_frames with the new one if header is already
     * decoded */
    if((ps_dec->i4_header_decoded & 1) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_is_valid) && (-1 != num_reorder_frames) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].u1_vui_parameters_present_flag) &&
       (1 == ps_dec->ps_sps[u1_seq_parameter_set_id].s_vui.u1_bitstream_restriction_flag) &&
       ((WORD32) ps_dec->ps_sps[u1_seq_parameter_set_id].s_vui.u4_num_reorder_frames !=
        num_reorder_frames))
    {
        ps_dec->u1_res_changed = 1;
        return IVD_RES_CHANGED;
    }

    /* In case bitstream read has exceeded the filled size, then return an error */
    if(EXCEED_OFFSET(ps_bitstrm))
    {
        return ERROR_INV_SPS_PPS_T;
    }

    /*--------------------------------------------------------------------*/
    /* All initializations to ps_dec are beyond this point                */
    /*--------------------------------------------------------------------*/
    {
        WORD32 reorder_depth = ih264d_get_dpb_size(ps_seq);
        if((1 == ps_seq->u1_vui_parameters_present_flag) &&
           (1 == ps_seq->s_vui.u1_bitstream_restriction_flag))
        {
            reorder_depth = ps_seq->s_vui.u4_num_reorder_frames + 1;
        }

        if(reorder_depth > H264_MAX_REF_PICS)
        {
            return ERROR_INV_SPS_PPS_T;
        }

        if(ps_seq->u1_frame_mbs_only_flag != 1) reorder_depth *= 2;
        ps_subset_seq->i4_reorder_depth = reorder_depth + DISPLAY_LATENCY;
    }
    ps_subset_seq->u2_disp_height = i4_cropped_ht;
    ps_subset_seq->u2_disp_width = i4_cropped_wd;
    ps_subset_seq->u2_pic_wd = u2_pic_wd;
    ps_subset_seq->u2_pic_ht = u2_pic_ht;

    /* Determining the Width and Height of Frame from that of Picture */
    ps_subset_seq->u2_frm_wd_y = u2_frm_wd_y;
    ps_subset_seq->u2_frm_ht_y = u2_frm_ht_y;
    ps_subset_seq->u2_frm_wd_uv = u2_frm_wd_uv;
    ps_subset_seq->u2_frm_ht_uv = u2_frm_ht_uv;

    ps_subset_seq->u1_pad_len_y_v = (UWORD8) (PAD_LEN_Y_V << (1 - u1_frm));
    ps_subset_seq->u1_pad_len_cr_v = (UWORD8) (PAD_LEN_UV_V << (1 - u1_frm));

    ps_subset_seq->u2_crop_offset_y = u2_crop_offset_y;
    ps_subset_seq->u2_crop_offset_uv = u2_crop_offset_uv;

    ps_seq->u1_is_valid = TRUE;
    ps_dec->ps_sps[u1_seq_parameter_set_id] = *ps_seq;
    ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id] = *ps_subset_seq;
    ps_svc_lyr_dec->ps_cur_subset_sps = &ps_svc_lyr_dec->ps_subset_sps[u1_seq_parameter_set_id];

    return OK;
}

/*!
**************************************************************************
* \if Function name : isvcd_parse_pps \endif
*
* \brief
*    Decodes Picture Parameter set
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_pps(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_bit_stream_t *ps_bitstrm)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD8 uc_temp;
    dec_seq_params_t *ps_sps = NULL;
    dec_pic_params_t *ps_pps = NULL;
    UWORD32 *pu4_bitstrm_buf = ps_dec->ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_dec->ps_bitstrm->u4_ofst;

    /* Variables used for error resilience checks */
    UWORD64 u8_temp;
    UWORD32 u4_temp;
    WORD32 i_temp;

    /* For High profile related syntax elements */
    UWORD8 u1_more_data_flag;
    WORD32 i4_i;

    /*--------------------------------------------------------------------*/
    /* Decode pic_parameter_set_id and find corresponding pic params      */
    /*--------------------------------------------------------------------*/
    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u4_temp & MASK_ERR_PIC_SET_ID) return ERROR_INV_SPS_PPS_T;
    ps_pps = ps_dec->pv_scratch_sps_pps;
    *ps_pps = ps_dec->ps_pps[u4_temp];
    ps_pps->u1_pic_parameter_set_id = (UWORD8) u4_temp;
    COPYTHECONTEXT("PPS: pic_parameter_set_id", ps_pps->u1_pic_parameter_set_id);

    /************************************************/
    /* initilization of High profile syntax element */
    /************************************************/
    ps_pps->i4_transform_8x8_mode_flag = 0;
    ps_pps->i4_pic_scaling_matrix_present_flag = 0;

    /*--------------------------------------------------------------------*/
    /* Decode seq_parameter_set_id and map it to a seq_parameter_set      */
    /*--------------------------------------------------------------------*/
    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u4_temp & MASK_ERR_SEQ_SET_ID) return ERROR_INV_SPS_PPS_T;
    COPYTHECONTEXT("PPS: seq_parameter_set_id", u4_temp);
    ps_sps = &ps_dec->ps_sps[u4_temp];
    ps_pps->ps_sps = ps_sps;

    /*--------------------------------------------------------------------*/
    /* Decode entropy_coding_mode                                         */
    /*--------------------------------------------------------------------*/
    ps_pps->u1_entropy_coding_mode = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("PPS: entropy_coding_mode_flag", ps_pps->u1_entropy_coding_mode);

    ps_pps->u1_pic_order_present_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("PPS: pic_order_present_flag", ps_pps->u1_pic_order_present_flag);

    /*--------------------------------------------------------------------*/
    /* Decode num_slice_groups_minus1                                     */
    /*--------------------------------------------------------------------*/
    u8_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf) + (UWORD64) 1;
    if(u8_temp != 1)
    {
        return ERROR_FEATURE_UNAVAIL;
    }
    ps_pps->u1_num_slice_groups = (UWORD8) u8_temp;
    COPYTHECONTEXT("PPS: num_slice_groups_minus1", ps_pps->u1_num_slice_groups - 1);

    /*--------------------------------------------------------------------*/
    /* Other parameter set values                                         */
    /*--------------------------------------------------------------------*/
    u8_temp = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u8_temp >= H264_MAX_REF_IDX) return ERROR_REF_IDX;
    ps_pps->u1_num_ref_idx_lx_active[0] = (UWORD8) u8_temp;
    COPYTHECONTEXT("PPS: num_ref_idx_l0_active_minus1", ps_pps->u1_num_ref_idx_lx_active[0] - 1);

    u8_temp = (UWORD64) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(u8_temp >= H264_MAX_REF_IDX) return ERROR_REF_IDX;
    ps_pps->u1_num_ref_idx_lx_active[1] = (UWORD8) u8_temp;
    COPYTHECONTEXT("PPS: num_ref_idx_l1_active_minus1", ps_pps->u1_num_ref_idx_lx_active[1] - 1);

    ps_pps->u1_wted_pred_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("PPS: weighted prediction u4_flag", ps_pps->u1_wted_pred_flag);
    uc_temp = (UWORD8) ih264d_get_bits_h264(ps_bitstrm, 2);
    COPYTHECONTEXT("PPS: weighted_bipred_idc", uc_temp);
    ps_pps->u1_wted_bipred_idc = uc_temp;

    if(ps_pps->u1_wted_bipred_idc > MAX_WEIGHT_BIPRED_IDC) return ERROR_INV_SPS_PPS_T;

    {
        WORD64 i8_temp = (WORD64) 26 + ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if((i8_temp < MIN_H264_QP) || (i8_temp > MAX_H264_QP)) return ERROR_INV_RANGE_QP_T;

        ps_pps->u1_pic_init_qp = (UWORD8) i8_temp;
        COPYTHECONTEXT("PPS: pic_init_qp_minus26", ps_pps->u1_pic_init_qp - 26);

        i8_temp = (WORD64) 26 + ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if((i8_temp < MIN_H264_QP) || (i8_temp > MAX_H264_QP)) return ERROR_INV_RANGE_QP_T;

        ps_pps->u1_pic_init_qs = (UWORD8) i8_temp;
        COPYTHECONTEXT("PPS: pic_init_qs_minus26", ps_pps->u1_pic_init_qs - 26);
    }

    i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if((i_temp < -12) || (i_temp > 12)) return ERROR_INV_RANGE_QP_T;
    ps_pps->i1_chroma_qp_index_offset = i_temp;
    COPYTHECONTEXT("PPS: chroma_qp_index_offset", ps_pps->i1_chroma_qp_index_offset);

    /***************************************************************************/
    /* initialize second_chroma_qp_index_offset to i1_chroma_qp_index_offset if */
    /* second_chroma_qp_index_offset is not present in bit-ps_bitstrm */
    /***************************************************************************/
    ps_pps->i1_second_chroma_qp_index_offset = ps_pps->i1_chroma_qp_index_offset;

    ps_pps->u1_deblocking_filter_parameters_present_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("PPS: deblocking_filter_control_present_flag",
                   ps_pps->u1_deblocking_filter_parameters_present_flag);
    ps_pps->u1_constrained_intra_pred_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("PPS: constrained_intra_pred_flag", ps_pps->u1_constrained_intra_pred_flag);
    ps_pps->u1_redundant_pic_cnt_present_flag = ih264d_get_bit_h264(ps_bitstrm);
    COPYTHECONTEXT("PPS: redundant_pic_cnt_present_flag",
                   ps_pps->u1_redundant_pic_cnt_present_flag);

    /* High profile related syntax elements */
    u1_more_data_flag = MORE_RBSP_DATA(ps_bitstrm);

    if(u1_more_data_flag)
    {
        /* read transform_8x8_mode_flag  */
        ps_pps->i4_transform_8x8_mode_flag = (WORD32) ih264d_get_bit_h264(ps_bitstrm);

        /* read pic_scaling_matrix_present_flag */
        ps_pps->i4_pic_scaling_matrix_present_flag = (WORD32) ih264d_get_bit_h264(ps_bitstrm);

        if(ps_pps->i4_pic_scaling_matrix_present_flag)
        {
            /* read the scaling matrices */
            for(i4_i = 0; i4_i < (6 + (ps_pps->i4_transform_8x8_mode_flag << 1)); i4_i++)
            {
                ps_pps->u1_pic_scaling_list_present_flag[i4_i] = ih264d_get_bit_h264(ps_bitstrm);

                if(ps_pps->u1_pic_scaling_list_present_flag[i4_i])
                {
                    WORD32 ret;
                    if(i4_i < 6)
                    {
                        ret = ih264d_scaling_list(
                            ps_pps->i2_pic_scalinglist4x4[i4_i], 16,
                            &ps_pps->u1_pic_use_default_scaling_matrix_flag[i4_i], ps_bitstrm);
                    }
                    else
                    {
                        ret = ih264d_scaling_list(
                            ps_pps->i2_pic_scalinglist8x8[i4_i - 6], 64,
                            &ps_pps->u1_pic_use_default_scaling_matrix_flag[i4_i], ps_bitstrm);
                    }

                    if(ret != OK)
                    {
                        return ret;
                    }
                }
            }
        }

        /* read second_chroma_qp_index_offset syntax element */
        i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if((i_temp < -12) || (i_temp > 12)) return ERROR_INV_RANGE_QP_T;

        ps_pps->i1_second_chroma_qp_index_offset = i_temp;
    }

    if(SCALABLE_BASELINE_PROFILE_IDC == ps_sps->u1_profile_idc)

    {
        if(ps_pps->u1_num_slice_groups > 7)
        {
            return ERROR_INV_SPS_PPS_T;
        }
    }

    /* In case bitstream read has exceeded the filled size, then return an error */
    if(EXCEED_OFFSET(ps_bitstrm))
    {
        return ERROR_INV_SPS_PPS_T;
    }
    ps_pps->u1_is_valid = TRUE;
    ps_dec->ps_pps[ps_pps->u1_pic_parameter_set_id] = *ps_pps;
    return OK;
}
