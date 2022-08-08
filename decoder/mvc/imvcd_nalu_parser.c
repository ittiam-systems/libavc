/******************************************************************************
 *
 * Copyright (C) 2021 The Android Open Source Project
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

/*****************************************************************************/
/*                                                                           */
/*  File Name         : imvcd_nalu_parser.h                                  */
/*                                                                           */
/*  Description       : Functions for MVC NALU parsing                       */
/*                                                                           */
/*****************************************************************************/
#include <string.h>

#include "ih264_typedefs.h"
#include "iv.h"
#include "imvcd.h"
#include "ih264d_error_handler.h"
#include "ih264d_bitstrm.h"
#include "ih264d_defs.h"
#include "ih264d_nal.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_structs.h"
#include "ih264d_vui.h"
#include "imvcd_defs.h"
#include "imvcd_slice_functions.h"
#include "imvcd_structs.h"
#include "imvcd_utils.h"

static WORD32 imvcd_nalu_mvc_ext_parser(mvc_dec_ctxt_t *ps_mvcd_ctxt, dec_bit_stream_t *ps_bitstrm)
{
    nalu_mvc_ext_t *ps_nalu_mvc_ext =
        &ps_mvcd_ctxt->as_nalu_mvc_ext[ps_mvcd_ctxt->u2_num_views_decoded];

    ps_nalu_mvc_ext->u1_non_idr_flag = ih264d_get_bit_h264(ps_bitstrm);
    ps_nalu_mvc_ext->u1_priority_id = ih264d_get_bits_h264(ps_bitstrm, 6);
    ps_nalu_mvc_ext->u2_view_id = ih264d_get_bits_h264(ps_bitstrm, 10);
    ps_nalu_mvc_ext->u1_temporal_id = ih264d_get_bits_h264(ps_bitstrm, 3);
    ps_nalu_mvc_ext->u1_anchor_pic_flag = ih264d_get_bit_h264(ps_bitstrm);
    ps_nalu_mvc_ext->u1_inter_view_flag = ih264d_get_bit_h264(ps_bitstrm);

    if(0 == ih264d_get_bit_h264(ps_bitstrm))
    {
        return IVD_INVALID_BITSTREAM;
    }

    if(ps_nalu_mvc_ext->u2_view_id >= MAX_NUM_VIEWS)
    {
        return IVD_INVALID_BITSTREAM;
    }

    return OK;
}

static WORD32 imvcd_parse_subset_sps(mvc_dec_ctxt_t *ps_mvcd_ctxt, dec_bit_stream_t *ps_bitstrm)
{
    subset_sps_t *ps_subset_sps;

    WORD32 i, j, k;
    UWORD8 u1_profile_idc;
    UWORD8 au1_constraint_set_flags[6];
    UWORD8 u1_level_idc;
    UWORD8 u1_seq_parameter_set_id;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD32 u4_temp;
    UWORD16 u2_num_views_m1;
    UWORD8 u1_num_refs;
    UWORD8 u1_num_level_values_signalled_m1;
    UWORD8 u2_num_ops_m1;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    WORD32 i4_error_code = OK;

    u1_profile_idc = ih264d_get_bits_h264(ps_bitstrm, 8);

    for(i = 0; i < 6; i++)
    {
        au1_constraint_set_flags[i] = ih264d_get_bit_h264(ps_bitstrm);
    }

    if((u1_profile_idc != MULTIVIEW_HIGH_PROFILE_IDC) || (au1_constraint_set_flags[1] == 1))
    {
        return ERROR_FEATURE_UNAVAIL;
    }

    /*****************************************************/
    /* Read reserved_zero_2bits (2 bits)                 */
    /*****************************************************/
    ih264d_get_bits_h264(ps_bitstrm, 2);

    u1_level_idc = ih264d_get_bits_h264(ps_bitstrm, 8);

    u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    u1_seq_parameter_set_id = u4_temp;

    if(u4_temp & MASK_ERR_SEQ_SET_ID)
    {
        return ERROR_INV_SPS_PPS_T;
    }

    ps_subset_sps = &ps_mvcd_ctxt->as_subset_sps[u1_seq_parameter_set_id];

    /* Accounting for the idiocy in 'ih264d_parse_pps' */
    ps_subset_sps->s_sps_data.u1_profile_idc = HIGH_PROFILE_IDC;
    ps_subset_sps->s_sps_data.u1_level_idc = u1_level_idc;
    ps_subset_sps->s_sps_data.u1_seq_parameter_set_id = u1_seq_parameter_set_id;

    ps_subset_sps->s_sps_data.i4_chroma_format_idc = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(ps_subset_sps->s_sps_data.i4_chroma_format_idc != 1)
    {
        return ERROR_FEATURE_UNAVAIL;
    }

    ps_subset_sps->s_sps_data.i4_bit_depth_luma_minus8 =
        ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(ps_subset_sps->s_sps_data.i4_bit_depth_luma_minus8 != 0)
    {
        return ERROR_FEATURE_UNAVAIL;
    }

    ps_subset_sps->s_sps_data.i4_bit_depth_chroma_minus8 =
        ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(ps_subset_sps->s_sps_data.i4_bit_depth_chroma_minus8 != 0)
    {
        return ERROR_FEATURE_UNAVAIL;
    }

    ps_subset_sps->s_sps_data.i4_qpprime_y_zero_transform_bypass_flag =
        ih264d_get_bit_h264(ps_bitstrm);

    if(ps_subset_sps->s_sps_data.i4_qpprime_y_zero_transform_bypass_flag != 0)
    {
        return ERROR_INV_SPS_PPS_T;
    }

    ps_subset_sps->s_sps_data.i4_seq_scaling_matrix_present_flag = ih264d_get_bit_h264(ps_bitstrm);

    if(ps_subset_sps->s_sps_data.i4_seq_scaling_matrix_present_flag)
    {
        for(i = 0; i < 8; i++)
        {
            ps_subset_sps->s_sps_data.u1_seq_scaling_list_present_flag[i] =
                ih264d_get_bit_h264(ps_bitstrm);
        }
    }

    ps_subset_sps->s_sps_data.u1_bits_in_frm_num =
        4 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(ps_subset_sps->s_sps_data.u1_bits_in_frm_num > MAX_BITS_IN_FRAME_NUM)
    {
        return ERROR_INV_SPS_PPS_T;
    }

    ps_subset_sps->s_sps_data.u2_u4_max_pic_num_minus1 =
        (1 << (ps_subset_sps->s_sps_data.u1_bits_in_frm_num)) - 1;

    ps_subset_sps->s_sps_data.u1_pic_order_cnt_type = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(ps_subset_sps->s_sps_data.u1_pic_order_cnt_type > MAX_PIC_ORDER_CNT_TYPE)
    {
        return ERROR_INV_POC_TYPE_T;
    }

    ps_subset_sps->s_sps_data.u1_num_ref_frames_in_pic_order_cnt_cycle = 1;

    if(ps_subset_sps->s_sps_data.u1_pic_order_cnt_type == 0)
    {
        ps_subset_sps->s_sps_data.u1_log2_max_pic_order_cnt_lsb_minus =
            4 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_subset_sps->s_sps_data.u1_log2_max_pic_order_cnt_lsb_minus > MAX_BITS_IN_POC_LSB)
        {
            return ERROR_INV_SPS_PPS_T;
        }

        ps_subset_sps->s_sps_data.i4_max_pic_order_cntLsb =
            (1 << ps_subset_sps->s_sps_data.u1_log2_max_pic_order_cnt_lsb_minus);
    }
    else if(ps_subset_sps->s_sps_data.u1_pic_order_cnt_type == 1)
    {
        ps_subset_sps->s_sps_data.u1_delta_pic_order_always_zero_flag =
            ih264d_get_bit_h264(ps_bitstrm);

        ps_subset_sps->s_sps_data.i4_ofst_for_non_ref_pic =
            ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        ps_subset_sps->s_sps_data.i4_ofst_for_top_to_bottom_field =
            ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        ps_subset_sps->s_sps_data.u1_num_ref_frames_in_pic_order_cnt_cycle =
            ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_subset_sps->s_sps_data.u1_num_ref_frames_in_pic_order_cnt_cycle > MVC_MAX_REF_PICS)
        {
            return ERROR_INV_SPS_PPS_T;
        }

        for(i = 0; i < ps_subset_sps->s_sps_data.u1_num_ref_frames_in_pic_order_cnt_cycle; i++)
        {
            ps_subset_sps->s_sps_data.i4_ofst_for_ref_frame[i] =
                ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }
    }

    ps_subset_sps->s_sps_data.u1_num_ref_frames = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if((ps_subset_sps->s_sps_data.u1_num_ref_frames > MVC_MAX_REF_PICS))
    {
        return ERROR_NUM_REF;
    }

    ps_subset_sps->s_sps_data.u1_gaps_in_frame_num_value_allowed_flag =
        ih264d_get_bit_h264(ps_bitstrm);

    ps_subset_sps->s_sps_data.u2_frm_wd_in_mbs = 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(ps_subset_sps->s_sps_data.u2_frm_wd_in_mbs > (H264_MAX_FRAME_WIDTH >> 4))
    {
        return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
    }

    ps_subset_sps->s_sps_data.u2_frm_ht_in_mbs = 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(ps_subset_sps->s_sps_data.u2_frm_ht_in_mbs > (H264_MAX_FRAME_HEIGHT >> 4))
    {
        return IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED;
    }

    ps_subset_sps->s_sps_data.u2_max_mb_addr =
        ps_subset_sps->s_sps_data.u2_frm_wd_in_mbs * ps_subset_sps->s_sps_data.u2_frm_ht_in_mbs - 1;

    ps_subset_sps->s_sps_data.u2_total_num_of_mbs = ps_subset_sps->s_sps_data.u2_max_mb_addr + 1;

    ps_subset_sps->s_sps_data.u1_frame_mbs_only_flag = ih264d_get_bit_h264(ps_bitstrm);

    if(!ps_subset_sps->s_sps_data.u1_frame_mbs_only_flag)
    {
        return ERROR_FEATURE_UNAVAIL;
    }

    ps_subset_sps->s_sps_data.u1_mb_aff_flag = 0;

    ps_subset_sps->s_sps_data.u1_direct_8x8_inference_flag = ih264d_get_bit_h264(ps_bitstrm);

    /* Frame cropping flag */
    u4_temp = ih264d_get_bit_h264(ps_bitstrm);

    if(u4_temp)
    {
        ps_subset_sps->s_disp_offsets.u2_left_offset =
            ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        ps_subset_sps->s_disp_offsets.u2_right_offset =
            ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        ps_subset_sps->s_disp_offsets.u2_top_offset = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        ps_subset_sps->s_disp_offsets.u2_bottom_offset =
            ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    }
    else
    {
        memset(&ps_subset_sps->s_disp_offsets, 0, sizeof(ps_subset_sps->s_disp_offsets));
    }

    ps_subset_sps->s_sps_data.u1_vui_parameters_present_flag = ih264d_get_bit_h264(ps_bitstrm);

    if(ps_subset_sps->s_sps_data.u1_vui_parameters_present_flag)
    {
        i4_error_code = ih264d_parse_vui_parametres(&ps_subset_sps->s_sps_data.s_vui, ps_bitstrm);

        if(i4_error_code != OK)
        {
            return i4_error_code;
        }
    }

    if(ih264d_get_bit_h264(ps_bitstrm) != 1)
    {
        return ERROR_INV_SPS_PPS_T;
    }

    u2_num_views_m1 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    ps_subset_sps->s_sps_mvc_ext.u2_num_views = 1 + u2_num_views_m1;

    if(u2_num_views_m1 > MAX_NUM_VIEWS)
    {
        return ERROR_INVALID_SEQ_PARAM;
    }

    if(ps_view_ctxt->i4_decode_header)
    {
        ps_mvcd_ctxt->u2_num_views = MAX(ps_mvcd_ctxt->u2_num_views, 1 + u2_num_views_m1);
    }
    else if(ps_mvcd_ctxt->u2_num_views != (1 + u2_num_views_m1))
    {
        return ERROR_INVALID_SEQ_PARAM;
    }

    for(i = 0; i <= u2_num_views_m1; i++)
    {
        ps_subset_sps->s_sps_mvc_ext.au2_view_ids[i] =
            ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    }

    for(i = 0; i < 2; i++)
    {
        ps_subset_sps->s_sps_mvc_ext.as_anchor_ref_data[i][0].u1_num_refs = 0;
        ps_subset_sps->s_sps_mvc_ext.as_non_anchor_ref_data[i][0].u1_num_refs = 0;
    }

    for(i = 1; i <= u2_num_views_m1; i++)
    {
        u1_num_refs = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        ps_subset_sps->s_sps_mvc_ext.as_anchor_ref_data[0][i].u1_num_refs = u1_num_refs;

        if(u1_num_refs > MAX_NUM_IVP_REFS)
        {
            return ERROR_INVALID_SEQ_PARAM;
        }

        for(j = 0; j < u1_num_refs; j++)
        {
            ps_subset_sps->s_sps_mvc_ext.as_anchor_ref_data[0][i].au2_ref_view_ids[j] =
                ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }

        u1_num_refs = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        ps_subset_sps->s_sps_mvc_ext.as_anchor_ref_data[1][i].u1_num_refs = u1_num_refs;

        if(u1_num_refs > MAX_NUM_IVP_REFS)
        {
            return ERROR_INVALID_SEQ_PARAM;
        }

        for(j = 0; j < u1_num_refs; j++)
        {
            ps_subset_sps->s_sps_mvc_ext.as_anchor_ref_data[1][i].au2_ref_view_ids[j] =
                ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }
    }

    for(i = 1; i <= u2_num_views_m1; i++)
    {
        u1_num_refs = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        ps_subset_sps->s_sps_mvc_ext.as_non_anchor_ref_data[0][i].u1_num_refs = u1_num_refs;

        if(u1_num_refs > MAX_NUM_IVP_REFS)
        {
            return ERROR_INVALID_SEQ_PARAM;
        }

        for(j = 0; j < u1_num_refs; j++)
        {
            ps_subset_sps->s_sps_mvc_ext.as_non_anchor_ref_data[0][i].au2_ref_view_ids[j] =
                ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }

        u1_num_refs = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        ps_subset_sps->s_sps_mvc_ext.as_non_anchor_ref_data[1][i].u1_num_refs = u1_num_refs;

        if(u1_num_refs > MAX_NUM_IVP_REFS)
        {
            return ERROR_INVALID_SEQ_PARAM;
        }

        for(j = 0; j < u1_num_refs; j++)
        {
            ps_subset_sps->s_sps_mvc_ext.as_non_anchor_ref_data[1][i].au2_ref_view_ids[j] =
                ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }
    }

    u1_num_level_values_signalled_m1 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    ps_subset_sps->s_sps_mvc_ext.u1_num_level_values_signalled =
        u1_num_level_values_signalled_m1 + 1;

    if(u1_num_level_values_signalled_m1 >= MAX_NUM_LEVEL_VALUES_SIGNALLED)
    {
        return ERROR_INVALID_SEQ_PARAM;
    }

    for(i = 0; i <= u1_num_level_values_signalled_m1; i++)
    {
        ps_subset_sps->s_sps_mvc_ext.as_mvc_level_info[i].u4_level_idc =
            ih264d_get_bits_h264(ps_bitstrm, 8);

        u2_num_ops_m1 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        ps_subset_sps->s_sps_mvc_ext.as_mvc_level_info->as_mvc_op_data[i].u2_num_ops =
            1 + u2_num_ops_m1;

        if(u2_num_ops_m1 >= MAX_NUM_OPERATING_POINTS)
        {
            return ERROR_INVALID_SEQ_PARAM;
        }

        for(j = 0; j <= u2_num_ops_m1; j++)
        {
            UWORD16 u2_num_target_views_m1;

            ps_subset_sps->s_sps_mvc_ext.as_mvc_level_info->as_mvc_op_data[j].u1_temporal_id =
                ih264d_get_bits_h264(ps_bitstrm, 3);

            u2_num_target_views_m1 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

            ps_subset_sps->s_sps_mvc_ext.as_mvc_level_info->as_mvc_op_data[j].u2_num_target_views =
                1 + u2_num_target_views_m1;

            if(u2_num_target_views_m1 >= MAX_NUM_VIEWS)
            {
                return ERROR_INVALID_SEQ_PARAM;
            }

            for(k = 0; k <= u2_num_target_views_m1; k++)
            {
                ps_subset_sps->s_sps_mvc_ext.as_mvc_level_info->as_mvc_op_data[j]
                    .au2_target_view_ids[k] = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            }

            ps_subset_sps->s_sps_mvc_ext.as_mvc_level_info->as_mvc_op_data[j].u2_num_views =
                (UWORD16) 1 + ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }
    }

    ps_subset_sps->u1_mvc_vui_parameters_present_flag = ih264d_get_bit_h264(ps_bitstrm);

    if(ps_subset_sps->u1_mvc_vui_parameters_present_flag)
    {
        return ERROR_INV_SPS_PPS_T;
    }

    /* In case bitstream read has exceeded the filled size, then
     return an error */
    if(EXCEED_OFFSET(ps_bitstrm))
    {
        return ERROR_INV_SPS_PPS_T;
    }

    ps_subset_sps->s_sps_data.u1_is_valid = 1;

    /* This ensures PPS has valid data in SPS array for reference */
    ps_view_ctxt->ps_sps[ps_subset_sps->s_sps_data.u1_seq_parameter_set_id] =
        ps_subset_sps->s_sps_data;

    ps_mvcd_ctxt->u1_num_subset_sps++;

    return OK;
}

/* This function removes emulation byte "0x03" from bitstream(EBSP to RBSP).
   It also converts bytestream format into 32 bit little - endian format. */
static WORD32 imvcd_transform_nalu(dec_bit_stream_t *ps_bitstrm, UWORD8 *pu1_nal_unit,
                                   UWORD32 u4_numbytes_in_nal_unit)
{
    UWORD32 ui4_word;
    UWORD8 u1_cur_byte;

    static const UWORD32 u4_num_bytes_in_word = sizeof(ui4_word) / sizeof(u1_cur_byte);
    UWORD32 u4_num_bytes_in_rbsp = 0;
    WORD32 i = 0, j;
    WORD8 c_count = 0;
    UWORD32 *puc_bitstream_buffer = (UWORD32 *) pu1_nal_unit;
    UWORD8 u1_nal_header_size = 1;
    UWORD8 u1_num_bytes_copied = 0;

    ps_bitstrm->pu4_buffer = puc_bitstream_buffer;

    ui4_word = *pu1_nal_unit++;
    u1_num_bytes_copied++;

    if((NAL_UNIT_TYPE(ui4_word) == PREFIX_NAL) ||
       (NAL_UNIT_TYPE(ui4_word) == CODED_SLICE_EXTENSION))
    {
        u1_nal_header_size += 3;
    }

    for(j = 0; j < 2; j++)
    {
        u1_cur_byte = *pu1_nal_unit++;

        ui4_word = ((ui4_word << 8) | u1_cur_byte);
        u1_num_bytes_copied++;

        c_count++;
        u4_num_bytes_in_rbsp++;

        if(u1_cur_byte != 0x00)
        {
            c_count = 0;
        }
    }

    if(u4_numbytes_in_nal_unit > 2)
    {
        i = ((u4_numbytes_in_nal_unit - 3));
    }

    for(; i > 8; i -= 4)
    {
        // loop 0
        u1_cur_byte = *pu1_nal_unit++;

        if(c_count == NUM_OF_ZERO_BYTES_BEFORE_START_CODE &&
           u1_cur_byte == EMULATION_PREVENTION_BYTE)
        {
            c_count = 0;
            u1_cur_byte = *pu1_nal_unit++;
            i--;
        }

        ui4_word = ((ui4_word << 8) | u1_cur_byte);
        u1_num_bytes_copied++;
        if(u4_num_bytes_in_word == u1_num_bytes_copied)
        {
            *puc_bitstream_buffer = ui4_word;
            puc_bitstream_buffer++;
            u1_num_bytes_copied = 0;
        }

        c_count++;
        if(u1_cur_byte != 0x00) c_count = 0;

        // loop 1
        u1_cur_byte = *pu1_nal_unit++;

        if(c_count == NUM_OF_ZERO_BYTES_BEFORE_START_CODE &&
           u1_cur_byte == EMULATION_PREVENTION_BYTE)
        {
            c_count = 0;
            u1_cur_byte = *pu1_nal_unit++;
            i--;
        }
        ui4_word = ((ui4_word << 8) | u1_cur_byte);
        u1_num_bytes_copied++;
        if(u4_num_bytes_in_word == u1_num_bytes_copied)
        {
            *puc_bitstream_buffer = ui4_word;
            puc_bitstream_buffer++;
            u1_num_bytes_copied = 0;
        }

        c_count++;
        if(u1_cur_byte != 0x00) c_count = 0;

        // loop 2
        u1_cur_byte = *pu1_nal_unit++;

        if(c_count == NUM_OF_ZERO_BYTES_BEFORE_START_CODE &&
           u1_cur_byte == EMULATION_PREVENTION_BYTE)
        {
            c_count = 0;
            u1_cur_byte = *pu1_nal_unit++;
            i--;
        }

        ui4_word = ((ui4_word << 8) | u1_cur_byte);
        u1_num_bytes_copied++;
        if(u4_num_bytes_in_word == u1_num_bytes_copied)
        {
            *puc_bitstream_buffer = ui4_word;
            puc_bitstream_buffer++;
            u1_num_bytes_copied = 0;
        }

        c_count++;
        if(u1_cur_byte != 0x00) c_count = 0;

        // loop 3
        u1_cur_byte = *pu1_nal_unit++;

        if(c_count == NUM_OF_ZERO_BYTES_BEFORE_START_CODE &&
           u1_cur_byte == EMULATION_PREVENTION_BYTE)
        {
            c_count = 0;
            u1_cur_byte = *pu1_nal_unit++;
            i--;
        }

        ui4_word = ((ui4_word << 8) | u1_cur_byte);
        u1_num_bytes_copied++;
        if(u4_num_bytes_in_word == u1_num_bytes_copied)
        {
            *puc_bitstream_buffer = ui4_word;
            puc_bitstream_buffer++;
            u1_num_bytes_copied = 0;
        }

        c_count++;
        if(u1_cur_byte != 0x00) c_count = 0;

        u4_num_bytes_in_rbsp += 4;
    }

    for(; i > 0; i--)
    {
        u1_cur_byte = *pu1_nal_unit++;

        if(c_count == NUM_OF_ZERO_BYTES_BEFORE_START_CODE &&
           u1_cur_byte == EMULATION_PREVENTION_BYTE)
        {
            c_count = 0;
            i--;
            u1_cur_byte = *pu1_nal_unit++;
        }

        ui4_word = ((ui4_word << 8) | u1_cur_byte);
        u4_num_bytes_in_rbsp++;

        if((u4_num_bytes_in_rbsp & 0x03) == 0x03)
        {
            *puc_bitstream_buffer = ui4_word;
            puc_bitstream_buffer++;
        }
        c_count++;
        if(u1_cur_byte != 0x00) c_count = 0;
    }

    *puc_bitstream_buffer = (ui4_word << ((3 - (((u4_num_bytes_in_rbsp << 30) >> 30))) << 3));
    ps_bitstrm->u4_ofst = 0;
    ps_bitstrm->u4_max_ofst = ((u4_num_bytes_in_rbsp + u1_nal_header_size) << 3);

    return (u4_num_bytes_in_rbsp);
}

WORD32 imvcd_nalu_parser(mvc_dec_ctxt_t *ps_mvcd_ctxt, UWORD8 *pu1_bitstream_buf,
                         UWORD32 i4_nalu_length)
{
    AVC_EXT_NALU_ID_T e_nalu_id;

    UWORD8 u1_first_byte;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    dec_bit_stream_t *ps_bitstrm = ps_view_ctxt->ps_bitstrm;

    WORD32 i4_error_code = NOT_OK;

    if((NULL != pu1_bitstream_buf) && (i4_nalu_length > 0))
    {
        imvcd_transform_nalu(ps_bitstrm, pu1_bitstream_buf, i4_nalu_length);

        u1_first_byte = ih264d_get_bits_h264(ps_bitstrm, 8);

        if(NAL_FORBIDDEN_BIT(u1_first_byte))
        {
            return NOT_OK;
        }

        e_nalu_id = NAL_UNIT_TYPE(u1_first_byte);
        ps_view_ctxt->u1_nal_unit_type = e_nalu_id;

        // if any other nal unit other than slice nal is encountered in between a
        // frame break out of loop without consuming header
        if((ps_view_ctxt->u4_slice_start_code_found == 1) &&
           (ps_view_ctxt->u1_pic_decode_done != 1) && is_slice_nalu_type(e_nalu_id))
        {
            return ERROR_INCOMPLETE_FRAME;
        }

        switch(e_nalu_id)
        {
            case PREFIX_NAL:
            {
                if(!ps_view_ctxt->i4_decode_header)
                {
                    if(1 == ih264d_get_bit_h264(ps_bitstrm))
                    {
                        return IVD_INVALID_BITSTREAM;
                    }

                    i4_error_code = imvcd_nalu_mvc_ext_parser(ps_mvcd_ctxt, ps_bitstrm);

                    if(i4_error_code != OK)
                    {
                        return i4_error_code;
                    }
                }

                break;
            }
            case SUBSET_SPS:
            {
                ih264d_rbsp_to_sodb(ps_view_ctxt->ps_bitstrm);

                i4_error_code = imvcd_parse_subset_sps(ps_mvcd_ctxt, ps_bitstrm);

                if(OK != i4_error_code)
                {
                    return i4_error_code;
                }

                ps_view_ctxt->i4_header_decoded |= 1 << SUBSET_SPS;

                break;
            }
            case SLICE_NON_IDR:
            case SLICE_IDR:
            {
                if(!ps_view_ctxt->i4_decode_header)
                {
                    if(is_header_decoded(ps_view_ctxt->i4_header_decoded, SPS) &&
                       is_header_decoded(ps_view_ctxt->i4_header_decoded, PPS))
                    {
                        nalu_mvc_ext_t *ps_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);

                        ps_view_ctxt->u4_slice_start_code_found = 1;

                        if((0 == ps_mvcd_ctxt->u2_num_views_decoded) &&
                           !ps_nalu_mvc_ext->u1_inter_view_flag)
                        {
                            ps_nalu_mvc_ext->u1_inter_view_flag = 1;
                        }

                        ih264d_rbsp_to_sodb(ps_view_ctxt->ps_bitstrm);

                        i4_error_code = imvcd_parse_decode_slice(ps_mvcd_ctxt);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }
                    else
                    {
                        return IVD_INVALID_BITSTREAM;
                    }
                }

                break;
            }
            case CODED_SLICE_EXTENSION:
            {
                if(!ps_view_ctxt->i4_decode_header)
                {
                    if(is_header_decoded(ps_view_ctxt->i4_header_decoded, SPS) &&
                       is_header_decoded(ps_view_ctxt->i4_header_decoded, PPS) &&
                       is_header_decoded(ps_view_ctxt->i4_header_decoded, SUBSET_SPS))
                    {
                        ps_view_ctxt->u4_slice_start_code_found = 1;

                        if(1 == ih264d_get_bit_h264(ps_bitstrm))
                        {
                            return IVD_INVALID_BITSTREAM;
                        }

                        i4_error_code = imvcd_nalu_mvc_ext_parser(ps_mvcd_ctxt, ps_bitstrm);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }

                        ih264d_rbsp_to_sodb(ps_view_ctxt->ps_bitstrm);

                        i4_error_code = imvcd_parse_decode_slice(ps_mvcd_ctxt);

                        if(i4_error_code != OK)
                        {
                            return i4_error_code;
                        }
                    }
                    else
                    {
                        return IVD_INVALID_BITSTREAM;
                    }
                }

                break;
            }
            default:
            {
                i4_error_code = ERROR_UNKNOWN_NAL;

                break;
            }
        }
    }

    return i4_error_code;
}
