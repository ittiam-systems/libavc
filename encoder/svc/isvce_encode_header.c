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
*  isvce_encode_header.c
*
* @brief
*  This file contains function definitions related to header encoding.
*
* @author
*  ittiam
*
* @par List of Functions:
*  - isvce_generate_sps()
*  - isvce_generate_pps()
*  - isvce_generate_slice_header()
*  - isvce_populate_sps()
*  - isvce_populate_pps()
*  - isvce_populate_slice_header()
*
*******************************************************************************
*/

#include "ih264_typedefs.h"
#include "ih264_debug.h"

/* Dependencies of ih264e_bitstream.h */
#include "ih264e_error.h"

#include "ih264e_bitstream.h"

#include "isvce_encode_header.h"
#include "isvce_utils.h"

static FORCEINLINE IH264E_ERROR_T isvce_generate_nal_unit_header(bitstrm_t *ps_bitstrm,
                                                                 WORD32 nal_unit_type,
                                                                 WORD32 nal_ref_idc)
{
    WORD32 return_status = IH264E_SUCCESS;

    if(!((nal_unit_type > 0) && (nal_unit_type < 32)))
    {
        return IH264E_FAIL;
    }

    /* forbidden_zero_bit + nal_ref_idc + nal_unit_type */
    PUT_BITS(ps_bitstrm, ((nal_ref_idc << 5) + nal_unit_type),
             (1 + 2 + 5), /*1 forbidden zero bit + 2 nal_ref_idc + 5 nal_unit_type */
             return_status, "nal_unit_header");

    return return_status;
}

/**
******************************************************************************
*
* @brief Generates SPS (Sequence Parameter Set)
*
* @par   Description
*  This function generates Sequence Parameter Set header as per the spec
*
* @param[in]   ps_bitstrm
*  pointer to bitstream context (handle)
*
* @param[in]   ps_sps
*  pointer to structure containing SPS data
*
* @param[in]   ps_vui
*  pointer to structure containing VUI data
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_generate_sps(bitstrm_t *ps_bitstrm, sps_t *ps_sps, NAL_UNIT_TYPE_T nal_type)
{
    WORD32 return_status = IH264E_SUCCESS;
    WORD32 i;
    WORD8 i1_nal_ref_idc = 3;
    vui_t *ps_vui = &ps_sps->s_vui_parameters;

    /* Insert Start Code */
    return_status = ih264e_put_nal_start_code_prefix(ps_bitstrm, 1);
    if(return_status != IH264E_SUCCESS)
    {
        return return_status;
    }
    /* Insert Nal Unit Header */
    return_status = isvce_generate_nal_unit_header(ps_bitstrm, nal_type, i1_nal_ref_idc);
    if(return_status != IH264E_SUCCESS)
    {
        return return_status;
    }

    /* profile_idc */
    PUT_BITS(ps_bitstrm, ps_sps->u1_profile_idc, 8, return_status, "profile_idc");

    /* constrained_set_flags */
    PUT_BITS(ps_bitstrm, ps_sps->u1_constraint_set0_flag, 1, return_status,
             "constrained_set0_flag");
    PUT_BITS(ps_bitstrm, ps_sps->u1_constraint_set1_flag, 1, return_status,
             "constrained_set1_flag");
    PUT_BITS(ps_bitstrm, ps_sps->u1_constraint_set2_flag, 1, return_status,
             "constrained_set2_flag");
    PUT_BITS(ps_bitstrm, ps_sps->u1_constraint_set3_flag, 1, return_status,
             "constrained_set3_flag");

    /* reserved_zero_four_bits */
    PUT_BITS(ps_bitstrm, 0, 4, return_status, "reserved_zero_four_bits");

    /* level_idc */
    PUT_BITS(ps_bitstrm, ps_sps->u1_level_idc, 8, return_status, "level_idc");

    /* seq_parameter_set_id */
    PUT_BITS_UEV(ps_bitstrm, ps_sps->u1_sps_id, return_status, "seq_parameter_set_id");

    if((ps_sps->u1_profile_idc == IH264_SCALABLE_BASELINE) ||
       (ps_sps->u1_profile_idc >= IH264_PROFILE_HIGH))
    {
        /* chroma_format_idc */
        PUT_BITS_UEV(ps_bitstrm, ps_sps->u1_chroma_format_idc, return_status, "chroma_format_idc");

        if(ps_sps->u1_chroma_format_idc == CHROMA_FMT_IDC_YUV444)
        {
            /* i1_residual_colour_transform_flag */
            PUT_BITS(ps_bitstrm, ps_sps->i1_residual_colour_transform_flag, 1, return_status,
                     "i1_residual_colour_transform_flag");
        }

        /* bit_depth_luma_minus8 */
        PUT_BITS_UEV(ps_bitstrm, (ps_sps->i1_bit_depth_luma - 8), return_status,
                     "bit_depth_luma_minus8");

        /* bit_depth_chroma_minus8 */
        PUT_BITS_UEV(ps_bitstrm, (ps_sps->i1_bit_depth_chroma - 8), return_status,
                     "bit_depth_chroma_minus8");

        /* qpprime_y_zero_transform_bypass_flag */
        PUT_BITS(ps_bitstrm, ps_sps->i1_qpprime_y_zero_transform_bypass_flag, 1, return_status,
                 "qpprime_y_zero_transform_bypass_flag");

        /* seq_scaling_matrix_present_flag */
        PUT_BITS(ps_bitstrm, ps_sps->i1_seq_scaling_matrix_present_flag, 1, return_status,
                 "seq_scaling_matrix_present_flag");

        /* seq_scaling_list */
        if(ps_sps->i1_seq_scaling_matrix_present_flag)
        {
            /* TODO_LATER: Will be enabled once scaling list support is added */
        }
    }

    /* log2_max_frame_num_minus4 */
    PUT_BITS_UEV(ps_bitstrm, (ps_sps->i1_log2_max_frame_num - 4), return_status,
                 "log2_max_frame_num_minus4");

    /* pic_order_cnt_type */
    PUT_BITS_UEV(ps_bitstrm, ps_sps->i1_pic_order_cnt_type, return_status, "pic_order_cnt_type");

    if(ps_sps->i1_pic_order_cnt_type == 0)
    {
        /* log2_max_pic_order_cnt_lsb_minus4 */
        PUT_BITS_UEV(ps_bitstrm, (ps_sps->i1_log2_max_pic_order_cnt_lsb - 4), return_status,
                     "log2_max_pic_order_cnt_lsb_minus4");
    }
    else if(ps_sps->i1_pic_order_cnt_type == 1)
    {
        /* delta_pic_order_always_zero_flag */
        PUT_BITS(ps_bitstrm, ps_sps->i1_delta_pic_order_always_zero_flag, 1, return_status,
                 "delta_pic_order_always_zero_flag");

        /* offset_for_non_ref_pic */
        PUT_BITS_SEV(ps_bitstrm, ps_sps->i4_offset_for_non_ref_pic, return_status,
                     "offset_for_non_ref_pic");

        /* offset_for_top_to_bottom_field */
        PUT_BITS_SEV(ps_bitstrm, ps_sps->i4_offset_for_top_to_bottom_field, return_status,
                     "offset_for_top_to_bottom_field");

        /* num_ref_frames_in_pic_order_cnt_cycle */
        PUT_BITS_UEV(ps_bitstrm, ps_sps->u1_num_ref_frames_in_pic_order_cnt_cycle, return_status,
                     "num_ref_frames_in_pic_order_cnt_cycle");

        /* Offset for ref frame */
        for(i = 0; i < ps_sps->u1_num_ref_frames_in_pic_order_cnt_cycle; i++)
        {
            /* offset_for_ref_frame */
            PUT_BITS_SEV(ps_bitstrm, ps_sps->ai4_offset_for_ref_frame[i], return_status,
                         "offset_for_ref_frame");
        }
    }

    /* num_ref_frames */
    PUT_BITS_UEV(ps_bitstrm, ps_sps->u1_max_num_ref_frames, return_status, "num_ref_frames");

    /* gaps_in_frame_num_value_allowed_flag */
    PUT_BITS(ps_bitstrm, ps_sps->i1_gaps_in_frame_num_value_allowed_flag, 1, return_status,
             "gaps_in_frame_num_value_allowed_flag");

    /* pic_width_in_mbs_minus1 */
    PUT_BITS_UEV(ps_bitstrm, ps_sps->i2_pic_width_in_mbs_minus1, return_status,
                 "pic_width_in_mbs_minus1");

    /* pic_height_in_map_units_minus1 */
    PUT_BITS_UEV(ps_bitstrm, ps_sps->i2_pic_height_in_map_units_minus1, return_status,
                 "pic_height_in_map_units_minus1");

    /* frame_mbs_only_flag */
    PUT_BITS(ps_bitstrm, ps_sps->i1_frame_mbs_only_flag, 1, return_status, "frame_mbs_only_flag");

    if(!ps_sps->i1_frame_mbs_only_flag)
    {
        /* mb_adaptive_frame_field_flag */
        PUT_BITS(ps_bitstrm, ps_sps->i1_mb_adaptive_frame_field_flag, 1, return_status,
                 "mb_adaptive_frame_field_flag");
    }

    /* direct_8x8_inference_flag */
    PUT_BITS(ps_bitstrm, ps_sps->i1_direct_8x8_inference_flag, 1, return_status,
             "direct_8x8_inference_flag");

    /* frame_cropping_flag */
    PUT_BITS(ps_bitstrm, ps_sps->i1_frame_cropping_flag, 1, return_status, "frame_cropping_flag");

    if(ps_sps->i1_frame_cropping_flag)
    {
        /* frame_crop_left_offset */
        PUT_BITS_UEV(ps_bitstrm, ps_sps->i2_frame_crop_left_offset, return_status,
                     "frame_crop_left_offset");

        /* frame_crop_right_offset */
        PUT_BITS_UEV(ps_bitstrm, ps_sps->i2_frame_crop_right_offset, return_status,
                     "frame_crop_right_offset");

        /* frame_crop_top_offset */
        PUT_BITS_UEV(ps_bitstrm, ps_sps->i2_frame_crop_top_offset, return_status,
                     "frame_crop_top_offset");

        /* frame_crop_bottom_offset */
        PUT_BITS_UEV(ps_bitstrm, ps_sps->i2_frame_crop_bottom_offset, return_status,
                     "frame_crop_bottom_offset");
    }

    /* vui_parameters_present_flag */
    PUT_BITS(ps_bitstrm, ps_sps->i1_vui_parameters_present_flag, 1, return_status,
             "vui_parameters_present_flag");

    if(ps_sps->i1_vui_parameters_present_flag)
    {
        /* Add vui parameters to the bitstream */;
        return_status = ih264e_generate_vui(ps_bitstrm, ps_vui);
        if(return_status != IH264E_SUCCESS)
        {
            return return_status;
        }
    }

    if(nal_type != NAL_SUBSET_SPS)
    {
        /* rbsp trailing bits */
        return_status = ih264e_put_rbsp_trailing_bits(ps_bitstrm);
    }

    return return_status;
}

/**
******************************************************************************
*
* @brief Generates PPS (Picture Parameter Set)
*
* @par   Description
*  Generate Picture Parameter Set as per Section 7.3.2.2
*
* @param[in]   ps_bitstrm
*  pointer to bitstream context (handle)
*
* @param[in]   ps_pps
*  pointer to structure containing PPS data
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_generate_pps(bitstrm_t *ps_bitstrm, pps_t *ps_pps, sps_t *ps_sps)
{
    WORD32 return_status = IH264E_SUCCESS;

    /* Insert the NAL start code */
    return_status = ih264e_put_nal_start_code_prefix(ps_bitstrm, 1);
    if(return_status != IH264E_SUCCESS)
    {
        return return_status;
    }

    /* Insert Nal Unit Header */
    PUT_BITS(ps_bitstrm, NAL_PPS_FIRST_BYTE, 8, return_status, "pps_header");

    /* pic_parameter_set_id */
    PUT_BITS_UEV(ps_bitstrm, ps_pps->u1_pps_id, return_status, "pic_parameter_set_id");

    /* seq_parameter_set_id */
    PUT_BITS_UEV(ps_bitstrm, ps_pps->u1_sps_id, return_status, "seq_parameter_set_id");

    /* Entropy coding : 0-VLC; 1 - CABAC */
    PUT_BITS(ps_bitstrm, ps_pps->u1_entropy_coding_mode_flag, 1, return_status,
             "Entropy coding : 0-VLC; 1 - CABAC");

    /* Pic order present flag */
    PUT_BITS(ps_bitstrm, ps_pps->u1_pic_order_present_flag, 1, return_status,
             "Pic order present flag");

    /* Number of slice groups */
    PUT_BITS_UEV(ps_bitstrm, ps_pps->u1_num_slice_groups - 1, return_status,
                 "Number of slice groups");

    if(ps_pps->u1_num_slice_groups > 1)
    {
        /* TODO_LATER: Currently the number of slice groups minus 1 is 0.
         * If this is not the case, we have to add Slice group map type to the bit
         * stream*/
    }

    /* num_ref_idx_l0_default_active_minus1 */
    PUT_BITS_UEV(ps_bitstrm, ps_pps->i1_num_ref_idx_l0_default_active - 1, return_status,
                 "num_ref_idx_l0_default_active_minus1");

    /* num_ref_idx_l1_default_active_minus1 */
    PUT_BITS_UEV(ps_bitstrm, ps_pps->i1_num_ref_idx_l1_default_active - 1, return_status,
                 "num_ref_idx_l1_default_active_minus1");

    /* weighted_pred_flag */
    PUT_BITS(ps_bitstrm, ps_pps->i1_weighted_pred_flag, 1, return_status, "weighted_pred_flag");

    /* weighted_bipred_flag */
    PUT_BITS(ps_bitstrm, ps_pps->i1_weighted_bipred_idc, 2, return_status, "weighted_bipred_idc");

    /* pic_init_qp_minus26 */
    PUT_BITS_SEV(ps_bitstrm, ps_pps->i1_pic_init_qp - 26, return_status, "pic_init_qp_minus26");

    /* pic_init_qs_minus26 */
    PUT_BITS_SEV(ps_bitstrm, ps_pps->i1_pic_init_qs - 26, return_status, "pic_init_qs_minus26");

    /* chroma_qp_index_offset */
    PUT_BITS_SEV(ps_bitstrm, ps_pps->i1_chroma_qp_index_offset, return_status,
                 "chroma_qp_index_offset");

    /* deblocking_filter_control_present_flag */
    PUT_BITS(ps_bitstrm, ps_pps->i1_deblocking_filter_control_present_flag, 1, return_status,
             "deblocking_filter_control_present_flag");

    /* constrained_intra_pred_flag */
    PUT_BITS(ps_bitstrm, ps_pps->i1_constrained_intra_pred_flag, 1, return_status,
             "constrained_intra_pred_flag");

    /*redundant_pic_cnt_present_flag */
    PUT_BITS(ps_bitstrm, ps_pps->i1_redundant_pic_cnt_present_flag, 1, return_status,
             "redundant_pic_cnt_present_flag");

    if(ps_sps->u1_profile_idc >= IH264_PROFILE_HIGH)
    {
        /* transform_8x8_mode_flag */
        PUT_BITS(ps_bitstrm, ps_pps->i1_transform_8x8_mode_flag, 1, return_status,
                 "transform_8x8_mode_flag");

        /* pic_scaling_matrix_present_flag */
        PUT_BITS(ps_bitstrm, ps_pps->i1_pic_scaling_matrix_present_flag, 1, return_status,
                 "pic_scaling_matrix_present_flag");

        if(ps_pps->i1_pic_scaling_matrix_present_flag)
        {
            /* TODO_LATER: Will be enabled once scaling list support is added */
        }

        /* Second chroma QP offset */
        PUT_BITS_SEV(ps_bitstrm, ps_pps->i1_second_chroma_qp_index_offset, return_status,
                     "Second chroma QP offset");
    }

    return_status = ih264e_put_rbsp_trailing_bits(ps_bitstrm);

    return return_status;
}

/**
******************************************************************************
*
* @brief Generates Slice Header
*
* @par   Description
*  Generate Slice Header as per Section 7.3.5.1
*
* @param[inout]   ps_bitstrm
*  pointer to bitstream context for generating slice header
*
* @param[in]   ps_slice_hdr
*  pointer to slice header params
*
* @param[in]   ps_pps
*  pointer to pps params referred by slice
*
* @param[in]   ps_sps
*  pointer to sps params referred by slice
*
* @param[out]   ps_dup_bit_strm_ent_offset
*  Bitstream struct to store bitstream state
*
* @param[out]   pu4_first_slice_start_offset
*  first slice offset is returned
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_generate_slice_header(bitstrm_t *ps_bitstrm, slice_header_t *ps_slice_hdr,
                                   pps_t *ps_pps, sps_t *ps_sps, UWORD8 u1_idr_flag)
{
    WORD32 return_status = IH264E_SUCCESS;
    UWORD8 u1_slice_type;

    /* Insert start code */
    return_status = ih264e_put_nal_start_code_prefix(ps_bitstrm, 1);
    if(return_status != IH264E_SUCCESS)
    {
        return return_status;
    }
    /* Insert Nal Unit Header */
    return_status = isvce_generate_nal_unit_header(ps_bitstrm, ps_slice_hdr->i1_nal_unit_type,
                                                   ps_slice_hdr->i1_nal_unit_idc);
    if(return_status != IH264E_SUCCESS)
    {
        return return_status;
    }
    /* first_mb_in_slice */
    PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u2_first_mb_in_slice, return_status,
                 "first_mb_in_slice");

    /* slice_type */
    PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u1_slice_type, return_status, "slice_type");

    /* pic_parameter_set_id */
    PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u1_pps_id, return_status, "pic_parameter_set_id");

    /* frame_num */
    PUT_BITS(ps_bitstrm, ps_slice_hdr->i4_frame_num, ps_sps->i1_log2_max_frame_num, return_status,
             "frame_num");

    if(!ps_sps->i1_frame_mbs_only_flag)
    {
        /* field_pic_flag */
        PUT_BITS(ps_bitstrm, ps_slice_hdr->i1_field_pic_flag, 1, return_status, "field_pic_flag");

        if(ps_slice_hdr->i1_field_pic_flag)
        {
            /* bottom_field_flag */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->i1_bottom_field_flag, 1, return_status,
                     "bottom_field_flag");
        }
    }

    if(u1_idr_flag == 1)
    {
        /* u2_idr_pic_id */
        PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u2_idr_pic_id, return_status, "idr_pic_id");
    }

    if(ps_sps->i1_pic_order_cnt_type == 0)
    {
        /* pic_order_cnt_lsb */
        PUT_BITS(ps_bitstrm, ps_slice_hdr->i4_pic_order_cnt_lsb,
                 ps_sps->i1_log2_max_pic_order_cnt_lsb, return_status, "pic_order_cnt_lsb");

        if(ps_pps->u1_pic_order_present_flag && !ps_slice_hdr->i1_field_pic_flag)
        {
            /* delta_pic_order_cnt_bottom */
            PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->i4_delta_pic_order_cnt_bottom, return_status,
                         "delta_pic_order_cnt_bottom");
        }
    }

    if(ps_sps->i1_pic_order_cnt_type == 1 && !ps_sps->i1_delta_pic_order_always_zero_flag)
    {
        /* delta_pic_order_cnt[0] */
        PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->ai4_delta_pic_order_cnt[0], return_status,
                     "delta_pic_order_cnt[0]");

        if(ps_pps->u1_pic_order_present_flag && !ps_slice_hdr->i1_field_pic_flag)
        {
            /* delta_pic_order_cnt[1] */
            PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->ai4_delta_pic_order_cnt[1], return_status,
                         "delta_pic_order_cnt[1]");
        }
    }

    if(ps_pps->i1_redundant_pic_cnt_present_flag)
    {
        /* redundant_pic_cnt */
        PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u1_redundant_pic_cnt, return_status,
                     "redundant_pic_cnt");
    }

    u1_slice_type = ps_slice_hdr->u1_slice_type % EPSLICE;

    if(u1_slice_type == BSLICE)
    {
        /* direct_spatial_mv_pred_flag */
        PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_direct_spatial_mv_pred_flag, 1, return_status,
                 "direct_spatial_mv_pred_flag");
    }

    if(u1_slice_type == PSLICE || u1_slice_type == SPSLICE || u1_slice_type == BSLICE)
    {
        /* num_ref_idx_active_override_flag */
        PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_num_ref_idx_active_override_flag, 1, return_status,
                 "num_ref_idx_active_override_flag");

        if(ps_slice_hdr->u1_num_ref_idx_active_override_flag)
        {
            /* num_ref_idx_l0_active_minus1 */
            PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->i1_num_ref_idx_l0_active - 1, return_status,
                         "num_ref_idx_l0_active_minus1");

            if(u1_slice_type == BSLICE)
            {
                /* num_ref_idx_l1_active_minus1 */
                PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->i1_num_ref_idx_l1_active - 1, return_status,
                             "num_ref_idx_l1_active_minus1");
            }
        }
    }

    /* ref_pic_list_modification */
    if((u1_slice_type != ISLICE) && (u1_slice_type != SISLICE))
    {
        /* ref_pic_list_modification_flag_l0 */
        PUT_BITS(ps_bitstrm, ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l0, 1,
                 return_status, "ref_pic_list_modification_flag_l0");

        if(ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l0)
        {
            UWORD8 i = 0;

            WORD8 *pi1_modification_of_pic_nums_idc_l0 =
                ps_slice_hdr->s_rplm.i1_modification_of_pic_nums_idc_l0;
            UWORD32 *pu4_abs_diff_pic_num_minus1_l0 =
                ps_slice_hdr->s_rplm.u4_abs_diff_pic_num_minus1_l0;
            UWORD8 *pu1_long_term_pic_num_l0 = ps_slice_hdr->s_rplm.u1_long_term_pic_num_l0;

            do
            {
                /* modification_of_pic_nums_idc */
                PUT_BITS_UEV(ps_bitstrm, pi1_modification_of_pic_nums_idc_l0[i], return_status,
                             "modification_of_pic_nums_idc");

                if((0 == pi1_modification_of_pic_nums_idc_l0[i]) ||
                   (1 == pi1_modification_of_pic_nums_idc_l0[i]))
                {
                    /* abs_diff_pic_num_minus1 */
                    PUT_BITS_UEV(ps_bitstrm, pu4_abs_diff_pic_num_minus1_l0[0], return_status,
                                 "abs_diff_pic_num_minus1");

                    pu4_abs_diff_pic_num_minus1_l0++;
                }
                else if(2 == pi1_modification_of_pic_nums_idc_l0[i])
                {
                    /* long_term_pic_num */
                    PUT_BITS_UEV(ps_bitstrm, pu1_long_term_pic_num_l0[0], return_status,
                                 "abs_diff_pic_num_minus1");

                    pu1_long_term_pic_num_l0++;
                }
            } while(pi1_modification_of_pic_nums_idc_l0[i++] != 3);
        }
    }

    if(u1_slice_type == BSLICE)
    {
        /* ref_pic_list_modification_flag_l1 */
        PUT_BITS(ps_bitstrm, ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l1, 1,
                 return_status, "ref_pic_list_modification_flag_l1");

        if(ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l1)
        {
            UWORD8 i = 0;

            WORD8 *pi1_modification_of_pic_nums_idc_l1 =
                ps_slice_hdr->s_rplm.i1_modification_of_pic_nums_idc_l1;
            UWORD32 *pu4_abs_diff_pic_num_minus1_l1 =
                ps_slice_hdr->s_rplm.u4_abs_diff_pic_num_minus1_l1;
            UWORD8 *pu1_long_term_pic_num_l1 = ps_slice_hdr->s_rplm.u1_long_term_pic_num_l1;

            do
            {
                /* modification_of_pic_nums_idc */
                PUT_BITS_UEV(ps_bitstrm, pi1_modification_of_pic_nums_idc_l1[i], return_status,
                             "modification_of_pic_nums_idc");

                if((0 == pi1_modification_of_pic_nums_idc_l1[i]) ||
                   (1 == pi1_modification_of_pic_nums_idc_l1[i]))
                {
                    /* abs_diff_pic_num_minus1 */
                    PUT_BITS_UEV(ps_bitstrm, pu4_abs_diff_pic_num_minus1_l1[0], return_status,
                                 "abs_diff_pic_num_minus1");

                    pu4_abs_diff_pic_num_minus1_l1++;
                }
                else if(2 == pi1_modification_of_pic_nums_idc_l1[i])
                {
                    /* long_term_pic_num */
                    PUT_BITS_UEV(ps_bitstrm, pu1_long_term_pic_num_l1[0], return_status,
                                 "abs_diff_pic_num_minus1");

                    pu1_long_term_pic_num_l1++;
                }
            } while(pi1_modification_of_pic_nums_idc_l1[i++] != 3);
        }
    }

    if((ps_pps->i1_weighted_pred_flag && u1_slice_type == PSLICE) ||
       (u1_slice_type == BSLICE && ps_pps->i1_weighted_bipred_idc == 1))
    {
        /* TODO_LATER: Currently there is no support for weighted prediction.
         This needs to be updated when the support is added */
    }

    if(ps_slice_hdr->i1_nal_unit_idc != 0)
    {
        if(u1_idr_flag == 1)
        {
            /* no_output_of_prior_pics_flag  */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_no_output_of_prior_pics_flag, 1, return_status,
                     "no_output_of_prior_pics_flag ");

            /* long_term_reference_flag  */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_long_term_reference_flag, 1, return_status,
                     "long_term_reference_flag ");
        }
        else
        {
            /* adaptive_ref_pic_marking_mode_flag  */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_adaptive_ref_pic_marking_mode_flag, 1,
                     return_status, "adaptive_ref_pic_marking_mode_flag ");

            if(ps_slice_hdr->u1_adaptive_ref_pic_marking_mode_flag)
            {
                /* TODO: if the reference picture marking mode is adaptive
                 add these fields in the bit-stream */
            }
        }
    }

    if(ps_slice_hdr->u1_entropy_coding_mode_flag && u1_slice_type != ISLICE &&
       u1_slice_type != SISLICE)
    {
        /* cabac_init_idc */
        PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->i1_cabac_init_idc, return_status, "cabac_init_idc");
    }

    /* slice_qp_delta */
    PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->i1_slice_qp - ps_pps->i1_pic_init_qp, return_status,
                 "slice_qp_delta");

    if(ps_slice_hdr->u1_slice_type == SPSLICE || ps_slice_hdr->u1_slice_type == SISLICE)
    {
        if(ps_slice_hdr->u1_slice_type == SPSLICE)
        {
            /* sp_for_switch_flag */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_sp_for_switch_flag, 1, return_status,
                     "sp_for_switch_flag");
        }
        /* slice_qs_delta */
        PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->u1_slice_qs - ps_pps->i1_pic_init_qs, return_status,
                     "slice_qs_delta");
    }

    if(ps_pps->i1_deblocking_filter_control_present_flag)
    {
        /* disable_deblocking_filter_idc */
        PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u1_disable_deblocking_filter_idc, return_status,
                     "disable_deblocking_filter_idc");

        if(ps_slice_hdr->u1_disable_deblocking_filter_idc != 1)
        {
            /* slice_alpha_c0_offset_div2 */
            PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->i1_slice_alpha_c0_offset_div2, return_status,
                         "slice_alpha_c0_offset_div2");

            /* slice_beta_offset_div2 */
            PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->i1_slice_beta_offset_div2, return_status,
                         "slice_beta_offset_div2");
        }
    }

    if(ps_slice_hdr->u1_num_slice_groups_minus1 > 0 && ps_pps->u1_slice_group_map_type >= 3 &&
       ps_pps->u1_slice_group_map_type <= 5)
    {
        /* slice_group_change_cycle */
        /* TODO_LATER: Currently the number of slice groups minus 1 is 0.
         * If this is not the case, we have to add Slice group map type to the bit
         * stream */
    }

    return return_status;
}

/**
******************************************************************************
*
* @brief Populates VUI structure
*
* @par   Description
*  Populates VUI structure for its use in header generation
*
* @param[in]   ps_codec
*  pointer to encoder context
*
* @return      success or failure error code
*
******************************************************************************
*/
static IH264E_ERROR_T isvce_populate_vui(isvce_codec_t *ps_codec, sps_t *ps_sps)
{
    vui_t *ps_vui = &ps_sps->s_vui_parameters;

    ps_vui->u1_nal_hrd_parameters_present_flag = 0;
    ps_vui->u1_vcl_hrd_parameters_present_flag = 0;

    ps_vui->u1_bitstream_restriction_flag = 1;
    ps_vui->u1_motion_vectors_over_pic_boundaries_flag = 1;
    ps_vui->u1_max_bytes_per_pic_denom = 0;
    ps_vui->u1_max_bits_per_mb_denom = 0;
    ps_vui->u1_log2_max_mv_length_horizontal = 16;
    ps_vui->u1_log2_max_mv_length_vertical = 16;

    if(ps_codec->s_cfg.u4_num_bframes == 0)
    {
        ps_vui->u1_num_reorder_frames = 0;
    }
    else
    {
        ps_vui->u1_num_reorder_frames = 1;
    }

    ps_vui->u1_max_dec_frame_buffering = ps_sps->u1_max_num_ref_frames;

    return 0;
}

/**
******************************************************************************
*
* @brief Populates sps structure
*
* @par   Description
*  Populates sps structure for its use in header generation
*
* @param[in]   ps_codec
*  pointer to encoder context
*
* @param[out]  ps_sps
*  pointer to sps params that needs to be populated
*
* @return      success or failure error code
*
******************************************************************************
*/
IH264E_ERROR_T isvce_populate_sps(isvce_codec_t *ps_codec, sps_t *ps_sps, UWORD8 u1_sps_id,
                                  UWORD8 u1_profile_idc, isvce_inp_buf_t *ps_inp_buf,
                                  UWORD8 u1_spatial_layer_id)
{
    /* active config parameters */
    isvce_cfg_params_t *ps_cfg = &(ps_codec->s_cfg);

    //    /* level */
    //    IH264_LEVEL_T   level_idc;

    /* error_status */
    IH264E_ERROR_T i4_err_code = IH264E_FAIL;

    /* profile */
    /*
     * Baseline profile supports, 8 bits per sample, 4:2:0 format, CAVLC.
     * B frames are not allowed. Further, Flexible mb ordering, Redundant slices,
     * Arbitrary slice ordering are supported. The constrained baseline profile is
     * baseline profile minus ASO, FMO and redundant slices. To the constrained
     * baseline profile if we add support for B slices, support for encoding
     * interlaced frames, support for weighted prediction and introduce CABAC
     * entropy coding then we have Main Profile.
     */
    ps_sps->u1_profile_idc = u1_profile_idc;

    /* level */
    ps_sps->u1_level_idc = MAX(
        ps_cfg->u4_max_level, (UWORD32) ih264e_get_min_level(ps_cfg->u4_max_wd, ps_cfg->u4_max_ht));

    /* constrained flags */
    /*
     * baseline profile automatically implies set 0 flag
     */
    ps_sps->u1_constraint_set0_flag = (ps_sps->u1_profile_idc == IH264_PROFILE_BASELINE);
    /*
     * main profile automatically implies set 1 flag
     * Although the encoder says it supports Baseline profile it actually supports
     * constrained baseline profile as ASO, FMO and redundant slices are not
     * supported
     */
    ps_sps->u1_constraint_set1_flag = (ps_sps->u1_profile_idc <= IH264_PROFILE_MAIN);
    /*
     * extended profile is not supported
     */
    ps_sps->u1_constraint_set2_flag = 0x00;
    /*
     * level 1b or level 11
     */
    if(ps_sps->u1_level_idc == IH264_LEVEL_1B)
    {
        ps_sps->u1_constraint_set3_flag = 0;
        ps_sps->u1_level_idc = IH264_LEVEL_11;
    }
    else
    {
        ps_sps->u1_constraint_set3_flag = 0;
    }

    /* active sps id */
    ps_sps->u1_sps_id = u1_sps_id;

    if((ps_sps->u1_profile_idc == IH264_SCALABLE_BASELINE) ||
       (ps_sps->u1_profile_idc >= IH264_PROFILE_HIGH))
    {
        /* chroma format idc */
        ps_sps->u1_chroma_format_idc = CHROMA_FMT_IDC_YUV420;

        /* residual_colour_transform_flag */
        ps_sps->i1_residual_colour_transform_flag = 0;

        /* luma bit depth 8 */
        ps_sps->i1_bit_depth_luma = 8;

        /* chroma bit depth 8 */
        ps_sps->i1_bit_depth_chroma = 8;

        /* qpprime_y_zero_transform_bypass_flag */
        ps_sps->i1_qpprime_y_zero_transform_bypass_flag = 0;

        /* seq_scaling_matrix_present_flag */
        ps_sps->i1_seq_scaling_matrix_present_flag = 0;

        if(ps_sps->i1_seq_scaling_matrix_present_flag)
        {
            /* TODO_LATER: Will be enabled once scaling list support is added */
        }
    }

    /* log2_max_frame_num_minus4 */
    ps_sps->i1_log2_max_frame_num = LOG2_MAX_FRAME_NUM_MINUS4 + 4;

    /* pic_order_cnt_type */
    ps_sps->i1_pic_order_cnt_type = 2;

    if(ps_codec->i4_non_ref_frames_in_stream)
    {
        ps_sps->i1_pic_order_cnt_type = 0;
    }

    /* log2_max_pic_order_cnt_lsb_minus4 */
    ps_sps->i1_log2_max_pic_order_cnt_lsb = 8;

    /* TODO : add support for other poc types */
    if(ps_sps->i1_pic_order_cnt_type == 0)
    {
    }
    else if(ps_sps->i1_pic_order_cnt_type == 1)
    {
    }

    ps_sps->u1_max_num_ref_frames = ps_codec->i4_max_num_reference_frames;

    /* gaps_in_frame_num_value_allowed_flag */
    ps_sps->i1_gaps_in_frame_num_value_allowed_flag = 0;

    /* pic width in mb - 1 */
    ps_sps->i2_pic_width_in_mbs_minus1 =
        (ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id].u4_width >> 4) - 1;

    /* pic height in mb - 1 */
    ps_sps->i2_pic_height_in_map_units_minus1 =
        (ps_inp_buf->as_layer_yuv_buf_props[u1_spatial_layer_id].u4_height >> 4) - 1;

    /* frame_mbs_only_flag, no support for interlace encoding */
    ps_sps->i1_frame_mbs_only_flag = 1;

    /* mb_adaptive_frame_field_flag */
    if(ps_sps->i1_frame_mbs_only_flag == 0)
    {
        ps_sps->i1_mb_adaptive_frame_field_flag = 0;
    }

    /* direct_8x8_inference_flag */
    if(ps_sps->u1_level_idc < IH264_LEVEL_30)
    {
        ps_sps->i1_direct_8x8_inference_flag = 0;
    }
    else
    {
        ps_sps->i1_direct_8x8_inference_flag = 1;
    }

    /* cropping params */
    /*NOTE : Cropping values depend on the chroma format
     * For our case ,decoder interprets the cropping values as 2*num pixels
     * Hence the difference in the disp width and width must be halved before
     * sending to get the expected results
     */
    ps_sps->i1_frame_cropping_flag = 0;
    ps_sps->i2_frame_crop_left_offset = 0;
    ps_sps->i2_frame_crop_right_offset = (ps_codec->s_cfg.u4_wd - ps_codec->s_cfg.u4_disp_wd) >> 1;
    ps_sps->i2_frame_crop_top_offset = 0;
    ps_sps->i2_frame_crop_bottom_offset = (ps_codec->s_cfg.u4_ht - ps_codec->s_cfg.u4_disp_ht) >> 1;

    if(ps_sps->i2_frame_crop_left_offset || ps_sps->i2_frame_crop_right_offset ||
       ps_sps->i2_frame_crop_top_offset || ps_sps->i2_frame_crop_bottom_offset)
    {
        ps_sps->i1_frame_cropping_flag = 1;
    }

    /* vui params */
    ps_sps->i1_vui_parameters_present_flag = !(ps_cfg->u4_disable_vui);

    if(!ps_sps->i1_vui_parameters_present_flag)
    {
        /* populate vui params */
        isvce_populate_vui(ps_codec, ps_sps);
    }
    else
    {
        ps_sps->s_vui_parameters = ps_cfg->s_vui;
    }

    return i4_err_code;
}

/**
******************************************************************************
*
* @brief Populates pps structure
*
* @par   Description
*  Populates pps structure for its use in header generation
*
* @param[in]   ps_codec
*  pointer to encoder context
*
* @param[out]  ps_pps
*  pointer to pps params that needs to be populated
*
* @return      success or failure error code
*
******************************************************************************
*/
IH264E_ERROR_T isvce_populate_pps(isvce_codec_t *ps_codec, pps_t *ps_pps, UWORD8 u1_sps_id,
                                  UWORD8 u1_pps_id, UWORD8 u1_spatial_layer_id)
{
    /* seq_parameter_set_id */
    ps_pps->u1_sps_id = u1_sps_id;

    /* pic_parameter_set_id */
    ps_pps->u1_pps_id = u1_pps_id;

    /* entropy_coding_mode */
    ps_pps->u1_entropy_coding_mode_flag =
        ((ps_codec->s_cfg.s_svc_params.u1_num_spatial_layers > 1) && (0 == u1_spatial_layer_id))
            ? CAVLC
            : ps_codec->s_cfg.u4_entropy_coding_mode;

    /* pic_order_present_flag is unset if we don't have feilds */
    ps_pps->u1_pic_order_present_flag = 0;

    /* Currently number of slice groups supported are 1 */
    ps_pps->u1_num_slice_groups = 1;

    if(ps_pps->u1_num_slice_groups - 1)
    {
        /* TODO_LATER: Currently the number of slice groups minus 1 is 0.
         * If this is not the case, we have to add Slice group map type to the bit
         * stream*/
    }

    /* number of reference frames for list 0 */
    /* FIXME : fix this hard coded value */
    ps_pps->i1_num_ref_idx_l0_default_active = 1;

    /* number of reference frames for list 1 */
    ps_pps->i1_num_ref_idx_l1_default_active = 1;

    /* weighted prediction for now is disabled */
    ps_pps->i1_weighted_pred_flag = 0;
    ps_pps->i1_weighted_bipred_idc = 0;

    /* The intent is to not signal qp from pps. Rather send the same in slice
     * headers */
    ps_pps->i1_pic_init_qp = 0;

    /* The intent is to not signal qp from pps. Rather send the same in slice
     * headers */
    ps_pps->i1_pic_init_qs = 0;

    /* The intent is to not signal qp from pps. Rather send the same in slice
     * headers */
    ps_pps->i1_chroma_qp_index_offset = 0;

    /* deblocking filter flags present in slice header */
    ps_pps->i1_deblocking_filter_control_present_flag = 1;

    /* constrained intra prediction */
    ps_pps->i1_constrained_intra_pred_flag =
        ps_codec->au4_constrained_intra_pred[u1_spatial_layer_id];

    /* sending redundant slices is not supported for now */
    ps_pps->i1_redundant_pic_cnt_present_flag = 0;

    ps_pps->u1_slice_group_map_type = 0;

    return IH264E_SUCCESS;
}

/**
******************************************************************************
*
* @brief Populates slice header structure
*
* @par   Description
*  Populates slice header structure for its use in header generation
*
* @param[in]  ps_proc
*  pointer to proc context
*
* @param[out]  ps_slice_hdr
*  pointer to slice header structure that needs to be populated
*
* @param[in]  ps_pps
*  pointer to pps params structure referred by the slice
*
* @param[in]   ps_sps
*  pointer to sps params referred by the pps
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_populate_slice_header(isvce_process_ctxt_t *ps_proc, slice_header_t *ps_slice_hdr,
                                   pps_t *ps_pps, sps_t *ps_sps, UWORD8 u1_is_idr)
{
    /* entropy context */
    isvce_entropy_ctxt_t *ps_entropy = &ps_proc->s_entropy;

    isvce_codec_t *ps_codec = ps_proc->ps_codec;

    if(ps_codec->u4_is_curr_frm_ref)
    {
        ps_slice_hdr->i1_nal_unit_idc = 3;
    }
    else
    {
        ps_slice_hdr->i1_nal_unit_idc = 0;
    }

    /* start mb address */
    ps_slice_hdr->u2_first_mb_in_slice = ps_entropy->i4_mb_start_add;

    /* slice type */
    ps_slice_hdr->u1_slice_type = ps_proc->i4_slice_type;

    /* pic_parameter_set_id */
    ps_slice_hdr->u1_pps_id = ps_pps->u1_pps_id;

    /* Separate color plane flag is 0,
     * hence the syntax element color_plane_id not included */

    /* frame num */
    ps_slice_hdr->i4_frame_num = ps_proc->i4_frame_num;

    /* frame_mbs_only_flag, no support for interlace encoding */
    if(!ps_sps->i1_frame_mbs_only_flag)
    {
        ps_slice_hdr->i1_field_pic_flag = 0;

        if(ps_slice_hdr->i1_field_pic_flag)
        {
            ps_slice_hdr->i1_bottom_field_flag = 0;
        }
    }

    /* idr pic id */
    if(u1_is_idr)
    {
        ps_slice_hdr->u2_idr_pic_id = ps_proc->u4_idr_pic_id;
        ps_slice_hdr->i1_nal_unit_type = NAL_SLICE_IDR;
    }
    else
    {
        ps_slice_hdr->i1_nal_unit_type = NAL_SLICE_NON_IDR;
    }

    if(ps_proc->u1_spatial_layer_id > 0)
    {
        ps_slice_hdr->i1_nal_unit_type = NAL_CODED_SLICE_EXTENSION;
    }

    if(ps_sps->i1_pic_order_cnt_type == 0)
    {
        WORD32 i4_poc;
        i4_poc = ps_codec->i4_poc;
        i4_poc %= (1 << ps_sps->i1_log2_max_pic_order_cnt_lsb);
        ps_slice_hdr->i4_pic_order_cnt_lsb = i4_poc;
    }
    /* TODO add support for poc type 1 */
    else if(ps_sps->i1_pic_order_cnt_type == 1)
    {
    }

    /*
     * redundant slices are not currently supported.
     * Hence the syntax element redundant slice cnt is not initialized
     */
    if(ps_pps->i1_redundant_pic_cnt_present_flag)
    {
    }

    /* direct spatial mv pred flag */
    if(ps_proc->i4_slice_type == BSLICE)
    {
        ps_slice_hdr->u1_direct_spatial_mv_pred_flag = 1;
    }

    if(ps_proc->i4_slice_type == PSLICE || ps_proc->i4_slice_type == SPSLICE ||
       ps_proc->i4_slice_type == BSLICE)
    {
        /* num_ref_idx_active_override_flag */
        ps_slice_hdr->u1_num_ref_idx_active_override_flag = 0;

        if(ps_slice_hdr->u1_num_ref_idx_active_override_flag)
        {
            /* num_ref_idx_l0_active_minus1 */

            if(ps_proc->i4_slice_type == BSLICE)
            {
                /* num_ref_idx_l1_active_minus1 */
            }
        }
    }

    /* ref_pic_list_modification */
    if((ps_proc->i4_slice_type != ISLICE) && (ps_proc->i4_slice_type != SISLICE))
    {
        /* ref_pic_list_modification_flag_l0 */
        ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l0 = 1;

        if(ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l0)
        {
            ps_slice_hdr->s_rplm.i1_modification_of_pic_nums_idc_l0[0] = 0;
            ps_slice_hdr->s_rplm.i1_modification_of_pic_nums_idc_l0[1] = 3;

            if((ps_codec->i4_frame_num - ps_proc->aps_ref_pic[L0]->i4_frame_num) < 1)
            {
                return IH264E_FAIL;
            }

            ps_slice_hdr->s_rplm.u4_abs_diff_pic_num_minus1_l0[0] =
                ps_codec->i4_frame_num - ps_proc->aps_ref_pic[L0]->i4_frame_num - 1;
        }
    }

    if(ps_proc->i4_slice_type == BSLICE)
    {
        /* ref_pic_list_modification_flag_l1 */
        ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l1 = 0;

        if(ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l1)
        {
        }
    }

    /* Currently we do not support weighted pred */
    /* ps_slice_hdr->u1_weighted_bipred_idc = 0; */

    if((ps_pps->i1_weighted_pred_flag &&
        (ps_proc->i4_slice_type == PSLICE || ps_proc->i4_slice_type == SPSLICE)) ||
       (ps_proc->i4_slice_type == BSLICE && ps_pps->i1_weighted_bipred_idc == 1))
    {
        /* TODO_LATER: Currently there is no support for weighted prediction.
             This needs to be updated when the support is added */
    }

    if(ps_slice_hdr->i1_nal_unit_idc != 0)
    {
        if(u1_is_idr)
        {
            /* no_output_of_prior_pics_flag  */
            ps_slice_hdr->u1_no_output_of_prior_pics_flag = 0;

            /* long_term_reference_flag  */
            ps_slice_hdr->u1_long_term_reference_flag = 0;
        }
        else
        {
            /* adaptive_ref_pic_marking_mode_flag  */
            ps_slice_hdr->u1_adaptive_ref_pic_marking_mode_flag = 0;

            if(ps_slice_hdr->u1_adaptive_ref_pic_marking_mode_flag)
            {
                /* TODO: if the reference picture marking mode is adaptive
                     add these fields in the bit-stream */
            }
        }
    }

    /* entropy coding mode flag */
    ps_slice_hdr->u1_entropy_coding_mode_flag = ps_entropy->u1_entropy_coding_mode_flag;

    if(ps_slice_hdr->u1_entropy_coding_mode_flag && ps_proc->i4_slice_type != ISLICE &&
       ps_proc->i4_slice_type != SISLICE)
    {
        /* cabac_init_idc */
        ps_slice_hdr->i1_cabac_init_idc = CABAC_INIT_IDC;
    }

    /* slice qp */
    ps_slice_hdr->i1_slice_qp = ps_proc->u1_frame_qp;

    if(ps_proc->i4_slice_type == SPSLICE || ps_proc->i4_slice_type == SISLICE)
    {
        if(ps_proc->i4_slice_type == SPSLICE)
        {
            /* sp_for_switch_flag */
        }
        /* slice_qs_delta */
    }

    if(ps_pps->i1_deblocking_filter_control_present_flag)
    {
        /* disable_deblocking_filter_idc */
        ps_slice_hdr->u1_disable_deblocking_filter_idc = ps_proc->u4_disable_deblock_level;

        if(ps_slice_hdr->u1_disable_deblocking_filter_idc != 1)
        {
            /* slice_alpha_c0_offset_div2 */
            ps_slice_hdr->i1_slice_alpha_c0_offset_div2 = 0;

            /* slice_beta_offset_div2 */
            ps_slice_hdr->i1_slice_beta_offset_div2 = 0;
        }
    }
    ps_slice_hdr->u1_num_slice_groups_minus1 = 0;
    if(ps_slice_hdr->u1_num_slice_groups_minus1 > 0 && ps_pps->u1_slice_group_map_type >= 3 &&
       ps_pps->u1_slice_group_map_type <= 5)
    {
        /* slice_group_change_cycle */
        /* TODO_LATER: Currently the number of slice groups minus 1 is 0.
         * If this is not the case, we have to add Slice group map type to the bit
         * stream */
    }

    ps_slice_hdr->i1_cabac_init_idc = CABAC_INIT_IDC;

    return IH264E_SUCCESS;
}

/**
******************************************************************************
*
* @brief Populates svc_nalu_ext structure
*
* @par   Description
*  Populates svc_nalu_ext structure for its use in header generation
*
* @param[in]  ps_proc
*  pointer to proc context
*
* @param[out]  ps_slice_hdr
*  pointer to slice header structure that needs to be populated
*
* @param[in]  ps_pps
*  pointer to pps params structure referred by the slice
*
* @param[in]   ps_sps
*  pointer to sps params referred by the pps
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_populate_svc_nalu_extension(isvce_process_ctxt_t *ps_proc,
                                         svc_nalu_ext_t *ps_svc_nalu_ext, NAL_UNIT_TYPE_T nalu_type,
                                         UWORD8 u1_idr_flag)
{
    isvce_codec_t *ps_codec = ps_proc->ps_codec;

    ps_svc_nalu_ext->u1_idr_flag = u1_idr_flag;

    ps_svc_nalu_ext->u1_priority_id = 0;

    ps_svc_nalu_ext->u1_no_inter_layer_pred_flag = ((nalu_type == NAL_PREFIX) ? 1 : 0);

    ps_svc_nalu_ext->u1_dependency_id =
        ((nalu_type == NAL_PREFIX) ? 0 : ps_proc->u1_spatial_layer_id);

    ps_svc_nalu_ext->u1_temporal_id = ps_proc->ps_cur_pic->i1_temporal_id;

    ps_svc_nalu_ext->u1_quality_id = 0;

    ps_svc_nalu_ext->u1_use_ref_base_pic_flag = 0;

    ps_svc_nalu_ext->u1_discardable_flag = 0;

    ps_svc_nalu_ext->u1_output_flag = 1;

    ps_svc_nalu_ext->u1_reserved_three_2bits = 3;

    ps_svc_nalu_ext->s_nalu_header.u1_nal_ref_idc = ps_codec->u4_is_curr_frm_ref ? 3 : 0;

    ps_svc_nalu_ext->s_nalu_header.u1_nal_unit_type = nalu_type;

    return IH264E_SUCCESS;
}

WORD32 isvce_populate_subset_sps(isvce_codec_t *ps_codec, subset_sps_t *ps_subset_sps,
                                 UWORD8 u1_sps_id, isvce_inp_buf_t *ps_inp_buf,
                                 UWORD8 u1_spatial_layer_id)
{
    sps_t *ps_sps = &ps_subset_sps->s_sps;

    isvce_populate_sps(ps_codec, ps_sps, u1_sps_id, IH264_SCALABLE_BASELINE, ps_inp_buf,
                       u1_spatial_layer_id);

    ps_subset_sps->s_sps_svc_ext.u1_inter_layer_deblocking_filter_control_present_flag = 1;

    ps_subset_sps->s_sps_svc_ext.i1_slice_header_restriction_flag = 1;

    ps_subset_sps->s_sps_svc_ext.u1_extended_spatial_scalability_idc = 0;

    ps_subset_sps->s_sps_svc_ext.i1_seq_tcoeff_level_prediction_flag = 0;

    ps_subset_sps->i1_svc_vui_parameters_present_flag = 0;

    ps_subset_sps->i1_additional_extension2_flag = 0;

    ps_subset_sps->s_sps_svc_ext.u1_chroma_phase_x_plus1 = 1;

    ps_subset_sps->s_sps_svc_ext.u1_chroma_phase_y_plus1 = 1;

    ps_subset_sps->s_sps_svc_ext.i1_adaptive_tcoeff_level_prediction_flag = 0;

    return IH264E_SUCCESS;
}

WORD32 isvce_populate_svc_slice(isvce_process_ctxt_t *ps_proc, svc_slice_header_t *ps_svc_slice_hdr,
                                pps_t *ps_pps, subset_sps_t *ps_subset_sps,
                                svc_nalu_ext_t *ps_svc_nalu_ext)
{
    WORD32 i4_return_status;

    i4_return_status =
        isvce_populate_slice_header(ps_proc, &ps_svc_slice_hdr->s_slice_header, ps_pps,
                                    &ps_subset_sps->s_sps, ps_svc_nalu_ext->u1_idr_flag);

    if(IH264E_SUCCESS != i4_return_status)
    {
        return IH264E_FAIL;
    }

    ps_svc_slice_hdr->i1_slice_skip_flag = 0;
    ps_svc_slice_hdr->i1_adaptive_residual_prediction_flag = ENABLE_RESIDUAL_PREDICTION;
    ps_svc_slice_hdr->i1_default_residual_prediction_flag = 0;
    ps_svc_slice_hdr->i1_adaptive_base_mode_flag = ENABLE_ILP_MV || ENABLE_IBL_MODE;
    ps_svc_slice_hdr->i1_default_base_mode_flag = 0;
    ps_svc_slice_hdr->i1_tcoeff_level_prediction_flag = 0;
    ps_svc_slice_hdr->i1_constrained_intra_resampling_flag = 0;
    ps_svc_slice_hdr->i1_adaptive_motion_prediction_flag = USE_ILP_MV_AS_MVP;
    ps_svc_slice_hdr->i1_default_motion_prediction_flag = 0;
    ps_svc_slice_hdr->u4_disable_inter_layer_deblocking_filter_idc =
        !ENABLE_INTRA_BASE_DEBLOCK || ps_proc->u4_disable_deblock_level;

    if(ps_svc_slice_hdr->u4_disable_inter_layer_deblocking_filter_idc != 1)
    {
        /* slice_alpha_c0_offset_div2 */
        ps_svc_slice_hdr->i4_inter_layer_slice_alpha_c0_offset_div2 = 0;

        /* slice_beta_offset_div2 */
        ps_svc_slice_hdr->i4_inter_layer_slice_beta_offset_div2 = 0;
    }

    if((ps_svc_nalu_ext->u1_quality_id == 0) && (ps_svc_nalu_ext->u1_no_inter_layer_pred_flag == 0))
    {
        ps_svc_slice_hdr->u4_ref_layer_dq_id = (ps_proc->u1_spatial_layer_id - 1) << 4;
    }

    return IH264E_SUCCESS;
}

/**
******************************************************************************
*
* @brief Signals prefix_nal_unit_rbsp
*
* @par   Description
*  prefix_nal_unit_rbsp as per Section G.7.3.2.12
*
* @param[inout]   ps_bitstrm
*  pointer to bitstream context for generating slice header
*
* @param[in]    svc_nalu_ext
*  pointer to svc NAL unit structure
*
* @param[in]    ps_slice_header
*  pointer to slice header
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_generate_prefix_nal(bitstrm_t *ps_bitstrm, svc_nalu_ext_t *ps_svc_nalu_ext,
                                 slice_header_t *ps_slice_header, UWORD8 u1_max_num_ref_frames,
                                 UWORD8 u1_num_spatial_layers)
{
    WORD32 return_status = IH264E_SUCCESS;

    WORD32 i4_store_ref_base_pic_flag = 0;
    WORD32 i4_additional_prefix_nal_unit_extension_flag = 0;

    if(ps_svc_nalu_ext->u1_dependency_id == (u1_num_spatial_layers - 1))
    {
        i4_store_ref_base_pic_flag = 1;
        if(u1_max_num_ref_frames < 2)
        {
            i4_store_ref_base_pic_flag = 0;
        }
    }

    /* store_ref_base_pic_flag */
    if(ps_svc_nalu_ext->s_nalu_header.u1_nal_ref_idc != 0)
    {
        PUT_BITS(ps_bitstrm, i4_store_ref_base_pic_flag, 1, return_status,
                 "store_ref_base_pic_flag");

        if((ps_svc_nalu_ext->u1_use_ref_base_pic_flag || i4_store_ref_base_pic_flag) &&
           !ps_svc_nalu_ext->u1_idr_flag)
        {
            PUT_BITS(ps_bitstrm, ps_slice_header->u1_adaptive_ref_pic_marking_mode_flag, 1,
                     return_status, "DPRM: adaptive_ref_base_pic_marking_mode_flag");
        }

        PUT_BITS(ps_bitstrm, i4_additional_prefix_nal_unit_extension_flag, 1, return_status,
                 "additional_prefix_nal_unit_extension_flag");
    }

    /* rbsp trailing bits */
    return_status = ih264e_put_rbsp_trailing_bits(ps_bitstrm);

    return return_status;
}

WORD32 isvce_generate_slice_header_svc(bitstrm_t *ps_bitstrm, pps_t *ps_pps,
                                       svc_nalu_ext_t *ps_svc_nalu_ext,
                                       svc_slice_header_t *ps_svc_slice_hdr,
                                       subset_sps_t *ps_subset_sps)
{
    WORD32 return_status = IH264E_SUCCESS;

    UWORD8 u1_slice_type;
    sps_t *ps_sps = &ps_subset_sps->s_sps;
    slice_header_t *ps_slice_hdr = &ps_svc_slice_hdr->s_slice_header;

    /* first_mb_in_slice */
    PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u2_first_mb_in_slice, return_status,
                 "SH: first_mb_in_slice");

    /* slice_type */
    PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u1_slice_type, return_status, "SH: slice_type");

    /* pic_parameter_set_id */
    PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u1_pps_id, return_status, "SH: pic_parameter_set_id");

    /* frame_num */
    PUT_BITS(ps_bitstrm, ps_slice_hdr->i4_frame_num, ps_sps->i1_log2_max_frame_num, return_status,
             "SH: frame_num");

    if(!ps_sps->i1_frame_mbs_only_flag)
    {
        /* field_pic_flag */
        PUT_BITS(ps_bitstrm, ps_slice_hdr->i1_field_pic_flag, 1, return_status,
                 "SH: field_pic_flag");

        if(ps_slice_hdr->i1_field_pic_flag)
        {
            /* bottom_field_flag */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->i1_bottom_field_flag, 1, return_status,
                     "SH: bottom_field_flag");
        }
    }

    if(ps_svc_nalu_ext->u1_idr_flag == 1)
    {
        /* u2_idr_pic_id */
        PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u2_idr_pic_id, return_status, "SH: idr_pic_id");
    }

    if(ps_sps->i1_pic_order_cnt_type == 0)
    {
        /* pic_order_cnt_lsb */
        PUT_BITS(ps_bitstrm, ps_slice_hdr->i4_pic_order_cnt_lsb,
                 ps_sps->i1_log2_max_pic_order_cnt_lsb, return_status, "SH: pic_order_cnt_lsb");

        if(ps_pps->u1_pic_order_present_flag && !ps_slice_hdr->i1_field_pic_flag)
        {
            /* delta_pic_order_cnt_bottom */
            PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->i4_delta_pic_order_cnt_bottom, return_status,
                         "SH: delta_pic_order_cnt_bottom");
        }
    }

    if(ps_sps->i1_pic_order_cnt_type == 1 && !ps_sps->i1_delta_pic_order_always_zero_flag)
    {
        /* delta_pic_order_cnt[0] */
        PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->ai4_delta_pic_order_cnt[0], return_status,
                     "SH: delta_pic_order_cnt[0]");

        if(ps_pps->u1_pic_order_present_flag && !ps_slice_hdr->i1_field_pic_flag)
        {
            /* delta_pic_order_cnt[1] */
            PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->ai4_delta_pic_order_cnt[1], return_status,
                         "SH: delta_pic_order_cnt[1]");
        }
    }

    if(ps_pps->i1_redundant_pic_cnt_present_flag)
    {
        /* redundant_pic_cnt */
        PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u1_redundant_pic_cnt, return_status,
                     "SH: redundant_pic_cnt");
    }

    u1_slice_type = ps_slice_hdr->u1_slice_type % EPSLICE;

    if(ps_svc_nalu_ext->u1_quality_id == 0)
    {
        if(u1_slice_type == BSLICE)
        {
            /* direct_spatial_mv_pred_flag */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_direct_spatial_mv_pred_flag, 1, return_status,
                     "SH: direct_spatial_mv_pred_flag");
        }

        if(u1_slice_type == PSLICE || u1_slice_type == BSLICE)
        {
            /* num_ref_idx_active_override_flag */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_num_ref_idx_active_override_flag, 1,
                     return_status, "SH: num_ref_idx_active_override_flag");

            if(ps_slice_hdr->u1_num_ref_idx_active_override_flag)
            {
                /* num_ref_idx_l0_active_minus1 */
                PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->i1_num_ref_idx_l0_active - 1, return_status,
                             "SH: num_ref_idx_l0_active_minus1");

                if(u1_slice_type == BSLICE)
                {
                    /* num_ref_idx_l1_active_minus1 */
                    PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->i1_num_ref_idx_l1_active - 1,
                                 return_status, "SH: num_ref_idx_l1_active_minus1");
                }
            }
        }

        /* ref_pic_list_modification */
        if((u1_slice_type != ISLICE) && (u1_slice_type != SISLICE))
        {
            /* ref_pic_list_modification_flag_l0 */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l0, 1,
                     return_status, "RPLR: ref_pic_list_reordering_flag");

            if(ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l0)
            {
                UWORD8 i = 0;

                WORD8 *pi1_modification_of_pic_nums_idc_l0 =
                    ps_slice_hdr->s_rplm.i1_modification_of_pic_nums_idc_l0;
                UWORD32 *pu4_abs_diff_pic_num_minus1_l0 =
                    ps_slice_hdr->s_rplm.u4_abs_diff_pic_num_minus1_l0;
                UWORD8 *pu1_long_term_pic_num_l0 = ps_slice_hdr->s_rplm.u1_long_term_pic_num_l0;

                do
                {
                    /* modification_of_pic_nums_idc */
                    PUT_BITS_UEV(ps_bitstrm, pi1_modification_of_pic_nums_idc_l0[i], return_status,
                                 "RPLR: reordering_of_pic_nums_idc");

                    if((0 == pi1_modification_of_pic_nums_idc_l0[i]) ||
                       (1 == pi1_modification_of_pic_nums_idc_l0[i]))
                    {
                        /* abs_diff_pic_num_minus1 */
                        PUT_BITS_UEV(ps_bitstrm, pu4_abs_diff_pic_num_minus1_l0[0], return_status,
                                     "RPLR: abs_diff_pic_num_minus1");

                        pu4_abs_diff_pic_num_minus1_l0++;
                    }
                    else if(2 == pi1_modification_of_pic_nums_idc_l0[i])
                    {
                        /* long_term_pic_num */
                        PUT_BITS_UEV(ps_bitstrm, pu1_long_term_pic_num_l0[0], return_status,
                                     "RPLR: long_term_pic_num");

                        pu1_long_term_pic_num_l0++;
                    }
                } while(pi1_modification_of_pic_nums_idc_l0[i++] != 3);
            }
        }

        if(u1_slice_type == BSLICE)
        {
            /* ref_pic_list_modification_flag_l1 */
            PUT_BITS(ps_bitstrm, ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l1, 1,
                     return_status, "SH: ref_pic_list_modification_flag_l1");

            if(ps_slice_hdr->s_rplm.i1_ref_pic_list_modification_flag_l1)
            {
                UWORD8 i = 0;

                WORD8 *pi1_modification_of_pic_nums_idc_l1 =
                    ps_slice_hdr->s_rplm.i1_modification_of_pic_nums_idc_l1;
                UWORD32 *pu4_abs_diff_pic_num_minus1_l1 =
                    ps_slice_hdr->s_rplm.u4_abs_diff_pic_num_minus1_l1;
                UWORD8 *pu1_long_term_pic_num_l1 = ps_slice_hdr->s_rplm.u1_long_term_pic_num_l1;

                do
                {
                    /* modification_of_pic_nums_idc */
                    PUT_BITS_UEV(ps_bitstrm, pi1_modification_of_pic_nums_idc_l1[i], return_status,
                                 "SH: modification_of_pic_nums_idc");

                    if((0 == pi1_modification_of_pic_nums_idc_l1[i]) ||
                       (1 == pi1_modification_of_pic_nums_idc_l1[i]))
                    {
                        /* abs_diff_pic_num_minus1 */
                        PUT_BITS_UEV(ps_bitstrm, pu4_abs_diff_pic_num_minus1_l1[0], return_status,
                                     "SH: abs_diff_pic_num_minus1");

                        pu4_abs_diff_pic_num_minus1_l1++;
                    }
                    else if(2 == pi1_modification_of_pic_nums_idc_l1[i])
                    {
                        /* long_term_pic_num */
                        PUT_BITS_UEV(ps_bitstrm, pu1_long_term_pic_num_l1[0], return_status,
                                     "SH: abs_diff_pic_num_minus1");

                        pu1_long_term_pic_num_l1++;
                    }
                } while(pi1_modification_of_pic_nums_idc_l1[i++] != 3);
            }
        }

        if((ps_pps->i1_weighted_pred_flag && u1_slice_type == PSLICE) ||
           (u1_slice_type == BSLICE && ps_pps->i1_weighted_bipred_idc == 1))
        {
            /* TODO_LATER: Currently there is no support for weighted prediction.
             This needs to be updated when the support is added */
        }

        if(ps_slice_hdr->i1_nal_unit_idc != 0)
        {
            if(ps_svc_nalu_ext->u1_idr_flag)
            {
                /* no_output_of_prior_pics_flag  */
                PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_no_output_of_prior_pics_flag, 1,
                         return_status, "DRPM: no_output_of_prior_pics_flag ");

                /* long_term_reference_flag  */
                PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_long_term_reference_flag, 1, return_status,
                         "DRPM: long_term_reference_flag ");
            }
            else
            {
                /* adaptive_ref_pic_marking_mode_flag  */
                PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_adaptive_ref_pic_marking_mode_flag, 1,
                         return_status, "DPRM: adaptive_ref_pic_marking_mode_flag ");

                if(ps_slice_hdr->u1_adaptive_ref_pic_marking_mode_flag)
                {
                    /* TODO: if the reference picture marking mode is adaptive
                     add these fields in the bit-stream */
                }
            }
            if(ps_subset_sps->s_sps_svc_ext.i1_slice_header_restriction_flag == 0)
            {
                WORD32 i4_store_ref_base_pic_flag = 0;

                if(ps_sps->u1_max_num_ref_frames >= 2)
                {
                    i4_store_ref_base_pic_flag = 1;
                }

                /* store_ref_base_pic_flag */
                PUT_BITS(ps_bitstrm, i4_store_ref_base_pic_flag, 1, return_status,
                         "SH: store_ref_base_pic_flag");

                if((ps_svc_nalu_ext->u1_use_ref_base_pic_flag || i4_store_ref_base_pic_flag) &&
                   (!ps_svc_nalu_ext->u1_idr_flag))
                {
                    PUT_BITS(ps_bitstrm, ps_slice_hdr->u1_adaptive_ref_pic_marking_mode_flag, 1,
                             return_status, "SH: adaptive_ref_base_pic_marking_mode_flag");
                }
            }
        }
    }

    if(ps_slice_hdr->u1_entropy_coding_mode_flag && u1_slice_type != ISLICE)
    {
        /* cabac_init_idc */
        PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->i1_cabac_init_idc, return_status,
                     "SH: cabac_init_idc");
    }

    /* slice_qp_delta */
    PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->i1_slice_qp - ps_pps->i1_pic_init_qp, return_status,
                 "SH: slice_qp_delta");

    if(ps_pps->i1_deblocking_filter_control_present_flag)
    {
        /* disable_deblocking_filter_idc */
        PUT_BITS_UEV(ps_bitstrm, ps_slice_hdr->u1_disable_deblocking_filter_idc, return_status,
                     "SH: disable_deblocking_filter_idc");

        if(ps_slice_hdr->u1_disable_deblocking_filter_idc != 1)
        {
            /* slice_alpha_c0_offset_div2 */
            PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->i1_slice_alpha_c0_offset_div2, return_status,
                         "SH: slice_alpha_c0_offset_div2");

            /* slice_beta_offset_div2 */
            PUT_BITS_SEV(ps_bitstrm, ps_slice_hdr->i1_slice_beta_offset_div2, return_status,
                         "SH: slice_beta_offset_div2");
        }
    }

    if(ps_slice_hdr->u1_num_slice_groups_minus1 > 0 && ps_pps->u1_slice_group_map_type >= 3 &&
       ps_pps->u1_slice_group_map_type <= 5)
    {
        /* slice_group_change_cycle */
        /* TODO_LATER: Currently the number of slice groups minus 1 is 0.
         * If this is not the case, we have to add Slice group map type to the bit
         * stream */
    }

    if((ps_svc_nalu_ext->u1_no_inter_layer_pred_flag == 0) && (ps_svc_nalu_ext->u1_quality_id == 0))
    {
        PUT_BITS_UEV(ps_bitstrm, ps_svc_slice_hdr->u4_ref_layer_dq_id, return_status,
                     "SH: ref_layer_dq_id");
        if(ps_subset_sps->s_sps_svc_ext.u1_inter_layer_deblocking_filter_control_present_flag)
        {
            PUT_BITS_UEV(ps_bitstrm, ps_svc_slice_hdr->u4_disable_inter_layer_deblocking_filter_idc,
                         return_status, "SH: disable_inter_layer_deblocking_filter_idc");
            if(ps_svc_slice_hdr->u4_disable_inter_layer_deblocking_filter_idc != 1)
            {
                PUT_BITS_SEV(ps_bitstrm,
                             ps_svc_slice_hdr->i4_inter_layer_slice_alpha_c0_offset_div2,
                             return_status, "SH: inter_layer_slice_alpha_c0_offset_div2");
                PUT_BITS_SEV(ps_bitstrm, ps_svc_slice_hdr->i4_inter_layer_slice_beta_offset_div2,
                             return_status, "SH: inter_layer_slice_beta_offset_div2");
            }
        }
        PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_constrained_intra_resampling_flag, 1,
                 return_status, "SH: constrained_intra_resampling_flag");
        if(ps_subset_sps->s_sps_svc_ext.u1_extended_spatial_scalability_idc == 2)
        {
            if(ps_sps->u1_chroma_format_idc > 0)
            {
                PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_ref_layer_chroma_phase_x_plus1_flag, 1,
                         return_status, "SH: ref_layer_chroma_phase_x_plus1_flag");
                PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_ref_layer_chroma_phase_y_plus1, 2,
                         return_status, "SH: ref_layer_chroma_phase_y_plus1");
            }
            PUT_BITS_SEV(ps_bitstrm, ps_svc_slice_hdr->i4_scaled_ref_layer_left, return_status,
                         "SH: scaled_ref_layer_left_offset");
            PUT_BITS_SEV(ps_bitstrm, ps_svc_slice_hdr->i4_scaled_ref_layer_top, return_status,
                         "SH: scaled_ref_layer_top_offset");
            PUT_BITS_SEV(ps_bitstrm, ps_svc_slice_hdr->i4_scaled_ref_layer_right, return_status,
                         "SH: scaled_ref_layer_right_offset");
            PUT_BITS_SEV(ps_bitstrm, ps_svc_slice_hdr->i4_scaled_ref_layer_bottom, return_status,
                         "SH: scaled_ref_layer_bottom_offset");
        }
    }

    if(!ps_svc_nalu_ext->u1_no_inter_layer_pred_flag)
    {
        PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_slice_skip_flag, 1, return_status,
                 "SH: slice_skip_flag");
        if(ps_svc_slice_hdr->i1_slice_skip_flag)
        {
            PUT_BITS_UEV(ps_bitstrm, ps_svc_slice_hdr->u4_num_mbs_in_slice_minus1, return_status,
                         "SH: num_mbs_in_slice_minus1");
        }
        else
        {
            PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_adaptive_base_mode_flag, 1, return_status,
                     "SH: adaptive_base_mode_flag");
            if(!ps_svc_slice_hdr->i1_adaptive_base_mode_flag)
                PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_default_base_mode_flag, 1, return_status,
                         "SH: default_base_mode_flag");

            if(!ps_svc_slice_hdr->i1_default_base_mode_flag)
            {
                PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_adaptive_motion_prediction_flag, 1,
                         return_status, "SH: adaptive_motion_prediction_flag");
                if(!ps_svc_slice_hdr->i1_adaptive_motion_prediction_flag)
                    PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_default_motion_prediction_flag, 1,
                             return_status, "SH: default_motion_prediction_flag");
            }
            PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_adaptive_residual_prediction_flag, 1,
                     return_status, "SH: adaptive_residual_prediction_flag");
            if(!ps_svc_slice_hdr->i1_adaptive_residual_prediction_flag)
                PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_default_residual_prediction_flag, 1,
                         return_status, "SH: default_residual_prediction_flag");
        }

        if(ps_subset_sps->s_sps_svc_ext.i1_adaptive_tcoeff_level_prediction_flag)
        {
            PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->i1_tcoeff_level_prediction_flag, 1,
                     return_status, "SH: tcoeff_level_prediction_flag");
        }
    }

    if(!ps_subset_sps->s_sps_svc_ext.i1_slice_header_restriction_flag &&
       !ps_svc_slice_hdr->i1_slice_skip_flag)
    {
        PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->u4_scan_idx_start, 4, return_status,
                 "SH: scan_idx_start");
        PUT_BITS(ps_bitstrm, ps_svc_slice_hdr->u4_scan_idx_end, 4, return_status,
                 "SH: scan_idx_end");
    }

    return return_status;
}

WORD32 isvce_seq_parameter_set_svc_extension(bitstrm_t *ps_bitstrm, subset_sps_t *ps_sub_sps,
                                             UWORD8 u1_chroma_format_idc)
{
    WORD32 return_status = IH264E_SUCCESS;

    /* inter_layer_deblocking_filter_control_present_flag */
    PUT_BITS(ps_bitstrm,
             ps_sub_sps->s_sps_svc_ext.u1_inter_layer_deblocking_filter_control_present_flag, 1,
             return_status, "SPS: inter_layer_deblocking_filter_control_present_flag");

    /* extended_spatial_scalability */
    PUT_BITS(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.u1_extended_spatial_scalability_idc, 2,
             return_status, "SPS: extended_spatial_scalability");

    if(u1_chroma_format_idc == 1 || u1_chroma_format_idc == 2)
    {
        /* chroma_phase_x_plus1_flag */
        PUT_BITS(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.u1_chroma_phase_x_plus1, 1, return_status,
                 "SPS: chroma_phase_x_plus1_flag");
    }

    if(u1_chroma_format_idc == 1)
    {
        /* chroma_phase_y_plus1 */
        PUT_BITS(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.u1_chroma_phase_y_plus1, 2, return_status,
                 "SPS: chroma_phase_y_plus1");
    }

    if(ps_sub_sps->s_sps_svc_ext.u1_extended_spatial_scalability_idc == 1)
    {
        if(u1_chroma_format_idc > 0)
        {
            /* seq_ref_layer_chroma_phase_x_plus1_flag */
            PUT_BITS(ps_bitstrm,
                     ps_sub_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_x_plus1_flag, 1,
                     return_status, "SPS: seq_ref_layer_chroma_phase_x_plus1_flag");

            /* seq_ref_layer_chroma_phase_y_plus1 */
            PUT_BITS(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_y_plus1, 2,
                     return_status, "SPS: seq_ref_layer_chroma_phase_y_plus1");
        }
        /* seq_scaled_ref_layer_left_offset */
        PUT_BITS_SEV(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_left_offset,
                     return_status, "SPS: seq_scaled_ref_layer_left_offset");

        /* seq_scaled_ref_layer_top_offset */
        PUT_BITS_SEV(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_top_offset,
                     return_status, "SPS: seq_scaled_ref_layer_top_offset");

        /* seq_scaled_ref_layer_right_offset */
        PUT_BITS_SEV(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_right_offset,
                     return_status, "SPS: seq_scaled_ref_layer_right_offset");

        /* seq_scaled_ref_layer_bottom_offset */
        PUT_BITS_SEV(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_bottom_offset,
                     return_status, "SPS: seq_scaled_ref_layer_bottom_offset");
    }

    /* seq_tcoeff_level_prediction_flag */
    PUT_BITS(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.i1_seq_tcoeff_level_prediction_flag, 1,
             return_status, "SPS: seq_tcoeff_level_prediction_flag");

    if(ps_sub_sps->s_sps_svc_ext.i1_seq_tcoeff_level_prediction_flag)
    {
        /* adaptive_tcoeff_level_prediction_flag */
        PUT_BITS(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.i1_adaptive_tcoeff_level_prediction_flag, 1,
                 return_status, "SPS: adaptive_tcoeff_level_prediction_flag");
    }

    /* slice_header_restriction_flag */
    PUT_BITS(ps_bitstrm, ps_sub_sps->s_sps_svc_ext.i1_slice_header_restriction_flag, 1,
             return_status, "SPS: slice_header_restriction_flag");

    return return_status;
}

WORD32 isvce_svc_vui_parameters_extension(bitstrm_t *ps_bitstrm, svc_vui_ext_t *ps_svc_vui)
{
    WORD32 return_status = IH264E_SUCCESS;
    UWORD32 i;

    PUT_BITS_UEV(ps_bitstrm, ps_svc_vui->u4_vui_ext_num_entries_minus1, return_status,
                 "num_layers_minus1");

    for(i = 0; i < ps_svc_vui->u4_vui_ext_num_entries_minus1; i++)
    {
        /* dependency_id */
        PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_dependency_id[i], 3, return_status,
                 "dependency_id");

        /* quality_id */
        PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_quality_id[i], 4, return_status, "quality_id");

        /* temporal_id */
        PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_temporal_id[i], 3, return_status,
                 "temporal_id");

        /* timing_info_present_flag */
        PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_timing_info_present_flag[i], 1, return_status,
                 "timing_info_present_flag");

        if(ps_svc_vui->u1_vui_ext_timing_info_present_flag[i])
        {
            /* num_units_in_tick */
            PUT_BITS(ps_bitstrm, ps_svc_vui->u4_vui_ext_num_units_in_tick[i], 32, return_status,
                     "num_units_in_tick");

            /* time_scale */
            PUT_BITS(ps_bitstrm, ps_svc_vui->u4_vui_ext_time_scale[i], 32, return_status,
                     "time_scale");

            /* fixed_frame_rate_flag */
            PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_fixed_frame_rate_flag[i], 1, return_status,
                     "fixed_frame_rate_flag");
        }

        /* nal_hrd_parameters_present_flag */
        PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_nal_hrd_params_present_flag[i], 1,
                 return_status, "nal_hrd_parameters_present_flag");

        if(ps_svc_vui->u1_vui_ext_nal_hrd_params_present_flag[i])
        {
        }

        /* nal_hrd_parameters_present_flag */
        PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_vcl_hrd_params_present_flag[i], 1,
                 return_status, "vcl_hrd_parameters_present_flag");

        if(ps_svc_vui->u1_vui_ext_vcl_hrd_params_present_flag[i])
        {
        }

        if(ps_svc_vui->u1_vui_ext_nal_hrd_params_present_flag[i] ||
           ps_svc_vui->u1_vui_ext_vcl_hrd_params_present_flag[i])
        {
            /* low_delay_hrd_flag */
            PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_low_delay_hrd_flag[i], 1, return_status,
                     "low_delay_hrd_flag");
        }

        /* pic_struct_present_flag */
        PUT_BITS(ps_bitstrm, ps_svc_vui->u1_vui_ext_pic_struct_present_flag[i], 1, return_status,
                 "pic_struct_present_flag");
    }

    return return_status;
}

WORD32 isvce_generate_subset_sps(bitstrm_t *ps_bitstrm, subset_sps_t *ps_subset_sps)
{
    WORD32 return_status = IH264E_SUCCESS;
    sps_t *ps_sps = &ps_subset_sps->s_sps;
    return_status = isvce_generate_sps(ps_bitstrm, &ps_subset_sps->s_sps, NAL_SUBSET_SPS);

    /* generate subset sps */
    if(ps_sps->u1_profile_idc == IH264_SCALABLE_BASELINE ||
       ps_sps->u1_profile_idc == IH264_SCALABLE_HIGH_PROFILE)
    {
        isvce_seq_parameter_set_svc_extension(ps_bitstrm, ps_subset_sps,
                                              ps_sps->u1_chroma_format_idc);

        /* svc_vui_parameters_present_flag */
        PUT_BITS(ps_bitstrm, ps_subset_sps->i1_svc_vui_parameters_present_flag, 1, return_status,
                 "SPS: svc_vui_parameters_present_flag");

        if(ps_subset_sps->i1_svc_vui_parameters_present_flag == 1)
        {
            svc_vui_ext_t *ps_svc_vui = NULL;
            isvce_svc_vui_parameters_extension(ps_bitstrm, ps_svc_vui);
        }

        /* additional_extension2_flag */
        PUT_BITS(ps_bitstrm, ps_subset_sps->i1_additional_extension2_flag, 1, return_status,
                 "SPS: additional_extension2_flag");
    }

    /* rbsp trailing bits */
    return_status = ih264e_put_rbsp_trailing_bits(ps_bitstrm);

    return return_status;
}
/**
******************************************************************************
*
* @brief Generates svc_nalu_ext
*
* @par   Description
*  Generate svc_nalu_ext as per Section G.7.3.1.1
*
* @param[inout]   ps_bitstrm
*  pointer to bitstream context for generating slice header
*
* @param[in]   ps_svc_nalu_ext
*  pointer to svc_nalu_ext struct
*
* @return    success or failure error code
*
******************************************************************************
*/
WORD32 isvce_generate_svc_nalu_extension(bitstrm_t *ps_bitstrm, svc_nalu_ext_t *ps_svc_nalu_ext,
                                         UWORD8 u1_nalu_id)
{
    WORD32 return_status = IH264E_SUCCESS;

    /* Insert start code */
    return_status = ih264e_put_nal_start_code_prefix(ps_bitstrm, 1);

    if(return_status != IH264E_SUCCESS)
    {
        return return_status;
    }

    /* Insert Nal Unit Header */
    return_status = isvce_generate_nal_unit_header(ps_bitstrm, u1_nalu_id, 3);

    if(return_status != IH264E_SUCCESS)
    {
        return return_status;
    }

    /* reserved_one_bit */
    PUT_BITS(ps_bitstrm, 1, 1, return_status, "NAL unit header: reserved_one_bit");

    /* idr_flag */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_idr_flag, 1, return_status,
             "NAL unit header: idr_flag");

    /* priority_id */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_priority_id, 6, return_status,
             "NAL unit header: priority_id");

    /* no_inter_layer_pred_flag */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_no_inter_layer_pred_flag, 1, return_status,
             "NAL unit header: no_inter_layer_pred_flag");

    /* dependency_id */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_dependency_id, 3, return_status,
             "NAL unit header: dependency_id");

    /* quality_id */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_quality_id, 4, return_status,
             "NAL unit header: quality_id");

    /* temporal_id */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_temporal_id, 3, return_status,
             "NAL unit header: temporal_id");

    /* use_ref_base_pic_flag */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_use_ref_base_pic_flag, 1, return_status,
             "NAL unit header: use_ref_base_pic_flag");

    /* discardable_flag */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_discardable_flag, 1, return_status,
             "NAL unit header: discardable_flag");

    /* output_flag */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_output_flag, 1, return_status,
             "NAL unit header: output_flag");

    /* reserved_three_2bits */
    PUT_BITS(ps_bitstrm, ps_svc_nalu_ext->u1_reserved_three_2bits, 2, return_status,
             "NAL unit header: reserved_three_2bits");

    return return_status;
}
