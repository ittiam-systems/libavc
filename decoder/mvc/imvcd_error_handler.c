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
/*  File Name         : imvcd_error_handler.c                                */
/*                                                                           */
/*  Description       : Functions for error handling                         */
/*                                                                           */
/*****************************************************************************/

#include "ih264_typedefs.h"
#include "iv.h"
#include "imvcd.h"
#include "ih264_macros.h"
#include "imvc_defs.h"
#include "ih264d_defs.h"
#include "ih264d_error_handler.h"
#include "ih264d_nal.h"
#include "ih264d_structs.h"
#include "imvcd_structs.h"
#include "imvcd_utils.h"

static IV_API_CALL_STATUS_T imvcd_check_invalid_numViews(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i;

    UWORD8 u1_num_valid_subset_sps_found = 0;
    UWORD16 u2_max_views = 1;

    if((ps_mvcd_ctxt->u1_num_subset_sps == 0) && (ps_mvcd_ctxt->u2_num_views_decoded > 0))
    {
        return IV_FAIL;
    }

    for(i = 0; i < MAX_NUM_SEQ_PARAMS; i++)
    {
        if(ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u1_is_valid)
        {
            if(ps_mvcd_ctxt->as_subset_sps[i].s_sps_mvc_ext.u2_num_views > MAX_NUM_VIEWS)
            {
                return IV_FAIL;
            }

            u2_max_views =
                MAX(u2_max_views, ps_mvcd_ctxt->as_subset_sps[i].s_sps_mvc_ext.u2_num_views);
            u1_num_valid_subset_sps_found++;
        }
    }

    if(u1_num_valid_subset_sps_found > ps_mvcd_ctxt->u1_num_subset_sps)
    {
        return IV_FAIL;
    }

    if(ps_mvcd_ctxt->u2_num_views > u2_max_views)
    {
        return IV_FAIL;
    }

    if(ps_mvcd_ctxt->u2_num_views_decoded >= u2_max_views)
    {
        return IV_FAIL;
    }

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_check_sps_and_subset_sps(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    UWORD32 i;
    UWORD32 u4_cnt;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    UWORD32 u4_num_sps = ps_mvcd_ctxt->u1_num_sps;
    UWORD32 u4_num_subset_sps = ps_mvcd_ctxt->u1_num_subset_sps;
    WORD32 i4_max_mb_addr = INT32_MIN;

    i = 0;
    u4_cnt = 0;

    while((u4_cnt < u4_num_sps) && (i < MAX_NUM_SEQ_PARAMS))
    {
        if(ps_view_ctxt->ps_sps[i].u1_is_valid)
        {
            u4_cnt++;

            if(i4_max_mb_addr == INT32_MIN)
            {
                i4_max_mb_addr = ps_view_ctxt->ps_sps[i].u2_max_mb_addr;
            }
            else if(i4_max_mb_addr != ps_view_ctxt->ps_sps[i].u2_max_mb_addr)
            {
                return IV_FAIL;
            }

            if(ps_view_ctxt->ps_sps[i].u2_max_mb_addr >
               imvcd_get_num_mbs_in_level(ps_view_ctxt->ps_sps[i].u1_level_idc))
            {
                return IV_FAIL;
            }

            if(ps_view_ctxt->ps_sps[i].u1_mb_aff_flag)
            {
                return IV_FAIL;
            }

            if(!ps_view_ctxt->ps_sps[i].u1_frame_mbs_only_flag)
            {
                return IV_FAIL;
            }
        }

        i++;
    }

    if(u4_cnt != u4_num_sps)
    {
        return IV_FAIL;
    }

    i = 0;
    u4_cnt = 0;

    while((u4_cnt < u4_num_subset_sps) && (i < MAX_NUM_SEQ_PARAMS))
    {
        if(ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u1_is_valid)
        {
            u4_cnt++;

            if(i4_max_mb_addr == INT32_MIN)
            {
                i4_max_mb_addr = ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u2_max_mb_addr;
            }
            else if(i4_max_mb_addr != ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u2_max_mb_addr)
            {
                return IV_FAIL;
            }

            if(ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u2_max_mb_addr >
               imvcd_get_num_mbs_in_level(ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u1_level_idc))
            {
                return IV_FAIL;
            }

            if(ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u1_mb_aff_flag)
            {
                return IV_FAIL;
            }

            if(!ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u1_frame_mbs_only_flag)
            {
                return IV_FAIL;
            }
        }

        i++;
    }

    if(u4_cnt != u4_num_subset_sps)
    {
        return IV_FAIL;
    }

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T imvcd_check_pps(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    bool b_is_valid_pps_found = false;

    for(i = 0; i < MAX_NUM_PIC_PARAMS; i++)
    {
        if(ps_view_ctxt->ps_pps[i].u1_is_valid)
        {
            b_is_valid_pps_found = true;

            if(ps_view_ctxt->ps_pps[i].u1_frame_cropping_flag)
            {
                return IV_FAIL;
            }

            if(ps_view_ctxt->ps_pps[i].u1_num_slice_groups != 1)
            {
                return IV_FAIL;
            }
        }
    }

    return b_is_valid_pps_found ? IV_SUCCESS : IV_FAIL;
}

static IV_API_CALL_STATUS_T imvcd_check_num_view_slices_in_au(mvc_dec_ctxt_t *ps_mvcd_ctxt,
                                                              imvcd_video_decode_ip_t *ps_ip)
{
    AVC_EXT_NALU_ID_T e_nalu_id;

    UWORD16 u2_num_view_slices_in_au = 0;
    UWORD8 *pu1_input_buffer = (UWORD8 *) ps_ip->s_ivd_ip.pv_stream_buffer;
    UWORD32 u4_num_bytes_remaining = ps_ip->s_ivd_ip.u4_num_Bytes;

    while(true)
    {
        UWORD32 u4_length_of_start_code = 0;
        UWORD32 u4_next_is_aud = 0;
        WORD32 i4_nalu_length = ih264d_find_start_code(pu1_input_buffer, 0, u4_num_bytes_remaining,
                                                       &u4_length_of_start_code, &u4_next_is_aud);

        if(i4_nalu_length <= 0)
        {
            break;
        }

        if((0 != u4_next_is_aud) && (1 != u4_next_is_aud))
        {
            break;
        }

        if(u4_length_of_start_code < (NUM_OF_ZERO_BYTES_BEFORE_START_CODE + 1))
        {
            break;
        }

        e_nalu_id = NAL_UNIT_TYPE(pu1_input_buffer[u4_length_of_start_code]);
        u2_num_view_slices_in_au += (SLICE_NON_IDR == e_nalu_id) || (SLICE_IDR == e_nalu_id) ||
                                    (CODED_SLICE_EXTENSION == e_nalu_id);

        if(((WORD64) u4_num_bytes_remaining) <=
           ((WORD64) (((WORD64) u4_length_of_start_code) + ((WORD64) i4_nalu_length))))
        {
            break;
        }
        else
        {
            pu1_input_buffer += u4_length_of_start_code + i4_nalu_length;
            u4_num_bytes_remaining -= u4_length_of_start_code + i4_nalu_length;
        }

        if(u2_num_view_slices_in_au == ps_mvcd_ctxt->u2_num_views)
        {
            break;
        }
    }

    return (u2_num_view_slices_in_au != ps_mvcd_ctxt->u2_num_views) ? IV_FAIL : IV_SUCCESS;
}

IV_API_CALL_STATUS_T imvcd_au_error_checks(mvc_dec_ctxt_t *ps_mvcd_ctxt,
                                           imvcd_video_decode_ip_t *ps_ip)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    if(ps_mvcd_ctxt->b_header_only_decode)
    {
        return IV_FAIL;
    }

    if(!ps_view_ctxt->init_done)
    {
        return IV_FAIL;
    }

    if(IV_FAIL == imvcd_check_invalid_numViews(ps_mvcd_ctxt))
    {
        return IV_FAIL;
    }

    if(IV_FAIL == imvcd_check_sps_and_subset_sps(ps_mvcd_ctxt))
    {
        return IV_FAIL;
    }

    if(IV_FAIL == imvcd_check_pps(ps_mvcd_ctxt))
    {
        return IV_FAIL;
    }

    if(!ps_mvcd_ctxt->b_flush_enabled)
    {
        if(IV_FAIL == imvcd_check_num_view_slices_in_au(ps_mvcd_ctxt, ps_ip))
        {
            return IV_FAIL;
        }
    }

    return IV_SUCCESS;
}

IV_API_CALL_STATUS_T imvcd_view_error_checks(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i;

    nalu_mvc_ext_t *ps_nalu_mvc_ext = imvcd_get_cur_nalu_mvc_ext(ps_mvcd_ctxt);
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    bool b_is_idr_slice = imvcd_is_idr_au(ps_mvcd_ctxt);

    if(b_is_idr_slice && (ps_view_ctxt->ps_cur_slice->u1_slice_type != ISLICE))
    {
        return IV_FAIL;
    }

    if(ps_view_ctxt->u1_first_slice_in_stream && !b_is_idr_slice)
    {
        return IV_FAIL;
    }

    if(ps_nalu_mvc_ext->u2_view_id != 0)
    {
        subset_sps_t *ps_subset_sps =
            ps_mvcd_ctxt
                ->aps_pps_id_to_subset_sps_map[ps_view_ctxt->ps_cur_pps->u1_pic_parameter_set_id];

        if((NULL == ps_subset_sps) || !ps_subset_sps->s_sps_data.u1_is_valid)
        {
            return IV_FAIL;
        }

        if(0 == ps_mvcd_ctxt->u1_num_subset_sps)
        {
            return IV_FAIL;
        }

        if(ps_nalu_mvc_ext->u2_view_id >= ps_subset_sps->s_sps_mvc_ext.u2_num_views)
        {
            return IV_FAIL;
        }
    }
    else
    {
        if(ps_mvcd_ctxt->u2_num_views_decoded > 0)
        {
            return IV_FAIL;
        }
    }

    if(!ps_view_ctxt->u4_first_slice_in_pic || (ps_view_ctxt->u2_cur_slice_num > 0))
    {
        return IV_FAIL;
    }

    if(ps_view_ctxt->u4_first_slice_in_pic &&
       (ps_view_ctxt->ps_cur_slice->u2_first_mb_in_slice != 0))
    {
        return IV_FAIL;
    }

    if(ps_view_ctxt->ps_cur_slice->u1_mmco_equalto5)
    {
        return IV_FAIL;
    }

    for(i = 0; i < ps_mvcd_ctxt->u2_num_views_decoded; i++)
    {
        if(ps_mvcd_ctxt->as_slices[i].i4_poc != ps_view_ctxt->ps_cur_slice->i4_poc)
        {
            return IV_FAIL;
        }

        if(ps_mvcd_ctxt->as_slices[i].u2_frame_num != ps_view_ctxt->ps_cur_slice->u2_frame_num)
        {
            return IV_FAIL;
        }
    }

    if(SKIP_NONE != ps_view_ctxt->u4_skip_frm_mask)
    {
        return IV_FAIL;
    }

    return IV_SUCCESS;
}
