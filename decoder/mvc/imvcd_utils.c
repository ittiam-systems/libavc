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
/*  File Name         : imvcd_utils.c                                        */
/*                                                                           */
/*  Description       : MVCD Utility functions used by 'imvcd_api.c'         */
/*                                                                           */
/*****************************************************************************/
#include <string.h>

#include "ih264_typedefs.h"
#include "iv.h"
#include "ih264_debug.h"
#include "ih264_disp_mgr.h"
#include "ih264_macros.h"
#include "ih264d_error_handler.h"
#include "ih264d_format_conv.h"
#include "ih264d_utils.h"
#include "imvcd_structs.h"
#include "imvcd_utils.h"

void imvcd_free_ref_bufs(mvc_au_buf_mgr_t *ps_mvc_au_buf_mgr,
                         mvc_au_mv_pred_buf_mgr_t *ps_mvc_au_mv_pred_buf_mgr, WORD32 i4_pic_buf_id)
{
    ih264_buf_mgr_release(ps_mvc_au_buf_mgr->ps_buf_mgr_ctxt, i4_pic_buf_id, BUF_MGR_REF);

    ih264_buf_mgr_release(ps_mvc_au_mv_pred_buf_mgr->ps_buf_mgr_ctxt,
                          ps_mvc_au_buf_mgr->au1_au_buf_id_to_mv_buf_id_map[i4_pic_buf_id],
                          BUF_MGR_REF);
}

void imvcd_release_all_ref_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt, WORD32 i4_num_bufs)
{
    WORD32 i;

    for(i = 0; i < i4_num_bufs; i++)
    {
        ih264_buf_mgr_release(ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt, i, BUF_MGR_REF);

        ih264_buf_mgr_release(ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.ps_buf_mgr_ctxt,
                              ps_mvcd_ctxt->s_mvc_au_buf_mgr.au1_au_buf_id_to_mv_buf_id_map[i],
                              BUF_MGR_REF);
    }
}

void imvcd_free_ref_and_io_bufs(mvc_au_buf_mgr_t *ps_mvc_au_buf_mgr,
                                mvc_au_mv_pred_buf_mgr_t *ps_mvc_au_mv_pred_buf_mgr,
                                WORD32 i4_pic_buf_id)
{
    ih264_buf_mgr_release(ps_mvc_au_buf_mgr->ps_buf_mgr_ctxt, i4_pic_buf_id,
                          BUF_MGR_REF | BUF_MGR_IO);

    ih264_buf_mgr_release(ps_mvc_au_mv_pred_buf_mgr->ps_buf_mgr_ctxt,
                          ps_mvc_au_buf_mgr->au1_au_buf_id_to_mv_buf_id_map[i4_pic_buf_id],
                          BUF_MGR_REF | BUF_MGR_IO);
}

void imvcd_release_all_ref_and_io_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt, WORD32 i4_num_bufs)
{
    WORD32 i;

    for(i = 0; i < i4_num_bufs; i++)
    {
        ih264_buf_mgr_release(ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt, i,
                              BUF_MGR_REF | BUF_MGR_IO);

        ih264_buf_mgr_release(ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.ps_buf_mgr_ctxt,
                              ps_mvcd_ctxt->s_mvc_au_buf_mgr.au1_au_buf_id_to_mv_buf_id_map[i],
                              BUF_MGR_REF | BUF_MGR_IO);
    }
}

bool is_header_decoded(WORD32 i4_header_decoded, AVC_EXT_NALU_ID_T e_nalu_id)
{
    /* Accounting for idiocy in 'ih264d_parse_nal_unit' */
    e_nalu_id = (SPS == e_nalu_id) ? UNSPEC_0 : ((PPS == e_nalu_id) ? SLICE_NON_IDR : e_nalu_id);
    return !!(i4_header_decoded & (1 << e_nalu_id));
}

bool is_mvc_nalu(AVC_EXT_NALU_ID_T e_nalu_id)
{
    switch(e_nalu_id)
    {
        case SLICE_NON_IDR:
        case SLICE_DPA:
        case SLICE_DPB:
        case SLICE_DPC:
        case SLICE_IDR:
        case PREFIX_NAL:
        case SUBSET_SPS:
        case CODED_SLICE_EXTENSION:
        {
            return true;
        }
        default:
        {
            return false;
        }
    }
}

bool is_slice_nalu_type(AVC_EXT_NALU_ID_T e_nalu_id)
{
    switch(e_nalu_id)
    {
        case SLICE_NON_IDR:
        case SLICE_DPA:
        case SLICE_DPB:
        case SLICE_DPC:
        case SLICE_IDR:
        case CODED_SLICE_EXTENSION:
        case PREFIX_NAL:
        {
            return true;
        }
        default:
        {
            return false;
        }
    }
}

nalu_mvc_ext_t *imvcd_get_cur_nalu_mvc_ext(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    return &ps_mvcd_ctxt->as_nalu_mvc_ext[ps_mvcd_ctxt->u2_num_views_decoded];
}

nalu_mvc_ext_t *imvcd_get_nalu_mvc_ext(nalu_mvc_ext_t *ps_nalu_mvc_exts,
                                       UWORD16 u2_num_views_decoded, UWORD16 u2_view_id)
{
    WORD32 i;

    for(i = 0; i < u2_num_views_decoded; i++)
    {
        if(ps_nalu_mvc_exts[i].u2_view_id == u2_view_id)
        {
            return &ps_nalu_mvc_exts[i];
        }
    }

    return NULL;
}

ref_pic_list_mod_data_t *imvcd_get_cur_ref_pic_list_mod_data(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    return &ps_mvcd_ctxt->as_ref_pic_list_mod_data[ps_mvcd_ctxt->u2_num_views_decoded];
}

subset_sps_t *imvcd_get_valid_subset_sps(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    if(0 != ps_mvcd_ctxt->u2_num_views_decoded)
    {
        dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

        return ps_mvcd_ctxt
            ->aps_pps_id_to_subset_sps_map[ps_view_ctxt->ps_cur_pps->u1_pic_parameter_set_id];
    }
    else
    {
        WORD32 i;

        for(i = 0; i < MAX_NUM_SEQ_PARAMS; i++)
        {
            if(ps_mvcd_ctxt->as_subset_sps[i].s_sps_data.u1_is_valid)
            {
                return &ps_mvcd_ctxt->as_subset_sps[i];
            }
        }

        return NULL;
    }
}

void imvcd_modulate_max_disp_seq(dec_struct_t *ps_view_ctxt)
{
    WORD64 i8_temp;

    i8_temp = ps_view_ctxt->i4_prev_max_display_seq + ps_view_ctxt->i4_max_poc +
              ps_view_ctxt->u1_max_dec_frame_buffering + 1;

    ps_view_ctxt->i4_prev_max_display_seq = IS_OUT_OF_RANGE_S32(i8_temp) ? 0 : ((WORD32) i8_temp);
    ps_view_ctxt->i4_max_poc = 0;
}

mv_pred_t imvcd_get_default_mv_pred(void)
{
    mv_pred_t s_mv_pred = {.i2_mv = {0},
                           .i1_ref_frame = {OUT_OF_RANGE_REF, OUT_OF_RANGE_REF},
                           .u1_col_ref_pic_idx = UINT8_MAX,
                           .u1_pic_type = UINT8_MAX};

    return s_mv_pred;
}

UWORD32 imvcd_get_max_num_ivp_refs(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    WORD32 i;

    subset_sps_t *ps_subset_sps = imvcd_get_valid_subset_sps(ps_mvcd_ctxt);

    UWORD32 u4_max_ivp_refs = 0;

    if(!ps_subset_sps)
    {
        return u4_max_ivp_refs;
    }

    for(i = 0; i < ps_subset_sps->s_sps_mvc_ext.u2_num_views; i++)
    {
        u4_max_ivp_refs = MAX(
            u4_max_ivp_refs, ps_subset_sps->s_sps_mvc_ext.as_anchor_ref_data[0][i].u1_num_refs +
                                 ps_subset_sps->s_sps_mvc_ext.as_anchor_ref_data[1][i].u1_num_refs);
        u4_max_ivp_refs =
            MAX(u4_max_ivp_refs,
                ps_subset_sps->s_sps_mvc_ext.as_non_anchor_ref_data[0][i].u1_num_refs +
                    ps_subset_sps->s_sps_mvc_ext.as_non_anchor_ref_data[1][i].u1_num_refs);
    }

    return u4_max_ivp_refs;
}

bool imvcd_is_idr_au(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    return (ps_mvcd_ctxt->u2_num_views_decoded > 1)
               ? !ps_mvcd_ctxt->as_nalu_mvc_ext->u1_non_idr_flag
               : (ps_mvcd_ctxt->ae_nalu_id[0] == SLICE_IDR);
}

coordinates_t imvcd_get_buf_pad_dims(bool b_is_chroma)
{
    coordinates_t s_dims;

    /* Vert pad is '4 * PAD_LEN_UV_V' to account for field Pics */
    if(b_is_chroma)
    {
        s_dims.i4_abscissa = (PAD_LEN_UV_H * 4);
        s_dims.i4_ordinate = (PAD_LEN_UV_V * 4);
    }
    else
    {
        s_dims.i4_abscissa = (PAD_LEN_Y_H * 2);
        s_dims.i4_ordinate = (PAD_LEN_Y_V * 4);
    }

    return s_dims;
}

WORD32 imvcd_get_ref_pic_pad_offset(WORD32 i4_stride, bool b_is_chroma)
{
    return !b_is_chroma ? (i4_stride * PAD_LEN_Y_V * 2 + PAD_LEN_Y_H)
                        : (i4_stride * PAD_LEN_UV_V * 2 + PAD_LEN_UV_H * 2);
}

UWORD32 imvcd_get_next_bits(dec_bit_stream_t *ps_bitstream)
{
    UWORD32 u4_next_word;

    NEXTBITS(u4_next_word, ps_bitstream->u4_ofst, ps_bitstream->pu4_buffer, 32);

    return u4_next_word;
}

void imvcd_set_view_buf_id_to_buf_map(dec_struct_t *ps_view_ctxt)
{
    WORD32 i, j;

    memset(ps_view_ctxt->apv_buf_id_pic_buf_map, 0, sizeof(ps_view_ctxt->apv_buf_id_pic_buf_map));

    for(i = 0; i < 2; i++)
    {
        for(j = 0; j < ps_view_ctxt->ps_cur_slice->u1_num_ref_idx_lx_active[i]; j++)
        {
            ps_view_ctxt
                ->apv_buf_id_pic_buf_map[ps_view_ctxt->ps_ref_pic_buf_lx[i][j]->u1_pic_buf_id] =
                (void *) ps_view_ctxt->ps_ref_pic_buf_lx[i][j];
        }
    }
}

IV_API_CALL_STATUS_T imvcd_get_next_display_au_buf(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    mvc_au_buffer_t *ps_au_buf;

    IV_API_CALL_STATUS_T e_retval = IV_FAIL;

    UWORD32 i;
    WORD32 i4_buf_id;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_au_buf =
        (mvc_au_buffer_t *) ih264_disp_mgr_get(&ps_mvcd_ctxt->s_mvc_disp_buf_mgr, &i4_buf_id);

    ps_view_ctxt->i4_display_index = DEFAULT_POC;

    if(ps_au_buf != NULL)
    {
        ps_view_ctxt->pv_disp_sei_params = &ps_au_buf->s_sei_pic;
        ps_view_ctxt->i4_display_index = ps_au_buf->i4_poc;
        ps_view_ctxt->u4_num_fld_in_frm += 2;
        ps_view_ctxt->s_disp_op.u4_ts = ps_au_buf->u4_time_stamp;

        e_retval = IV_SUCCESS;
    }

    if(ps_au_buf)
    {
        for(i = 0; i < ps_mvcd_ctxt->u2_num_views; i++)
        {
            yuv_buf_props_t *ps_src = &ps_au_buf->as_view_buffers[i];
            yuv_buf_props_t *ps_dst = &ps_mvcd_ctxt->s_out_buffer.as_view_buf_props[i];

            WORD32 i4_y_src_stride = ps_src->as_component_bufs[Y].i4_data_stride;
            WORD32 i4_uv_src_stride = ps_src->as_component_bufs[UV].i4_data_stride;
            UWORD8 *pu1_y_src = (UWORD8 *) ps_src->as_component_bufs[Y].pv_data;
            UWORD8 *pu1_uv_src = (UWORD8 *) ps_src->as_component_bufs[UV].pv_data;

            pu1_y_src += (0 == i) ? ps_view_ctxt->u2_crop_offset_y
                                  : (ps_au_buf->as_disp_offsets[i].u2_left_offset +
                                     ps_au_buf->as_disp_offsets[i].u2_top_offset * i4_y_src_stride);
            pu1_uv_src +=
                (0 == i) ? ps_view_ctxt->u2_crop_offset_uv
                         : (ps_au_buf->as_disp_offsets[i].u2_left_offset +
                            (ps_au_buf->as_disp_offsets[i].u2_top_offset * i4_uv_src_stride) / 2);

            ps_dst->u2_width = ps_au_buf->u2_disp_width;
            ps_dst->u2_height = ps_au_buf->u2_disp_height;

            ASSERT(ps_dst->as_component_bufs[U].i4_data_stride ==
                   ps_dst->as_component_bufs[V].i4_data_stride);

            ih264d_fmt_conv_420sp_to_420p(
                pu1_y_src, pu1_uv_src, (UWORD8 *) ps_dst->as_component_bufs[Y].pv_data,
                (UWORD8 *) ps_dst->as_component_bufs[U].pv_data,
                (UWORD8 *) ps_dst->as_component_bufs[V].pv_data, ps_dst->u2_width,
                ps_dst->u2_height, i4_y_src_stride, i4_uv_src_stride,
                ps_dst->as_component_bufs[Y].i4_data_stride,
                ps_dst->as_component_bufs[U].i4_data_stride, 1, 0);
        }

        ih264_buf_mgr_release(ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt,
                              ps_au_buf->i4_pic_buf_id, BUF_MGR_IO);
    }

    return e_retval;
}

UWORD32 imvcd_get_num_mbs_in_level(UWORD8 u1_level_idc)
{
    switch(u1_level_idc)
    {
        case H264_LEVEL_1_0:
        {
            return MAX_MBS_LEVEL_10;
        }
        case H264_LEVEL_1_1:
        {
            return MAX_MBS_LEVEL_11;
        }
        case H264_LEVEL_1_2:
        {
            return MAX_MBS_LEVEL_12;
        }
        case H264_LEVEL_1_3:
        {
            return MAX_MBS_LEVEL_13;
        }
        case H264_LEVEL_2_0:
        {
            return MAX_MBS_LEVEL_20;
        }
        case H264_LEVEL_2_1:
        {
            return MAX_MBS_LEVEL_21;
        }
        case H264_LEVEL_2_2:
        {
            return MAX_MBS_LEVEL_22;
        }
        case H264_LEVEL_3_0:
        {
            return MAX_MBS_LEVEL_30;
        }
        case H264_LEVEL_3_1:
        {
            return MAX_MBS_LEVEL_31;
        }
        case H264_LEVEL_3_2:
        {
            return MAX_MBS_LEVEL_32;
        }
        case H264_LEVEL_4_0:
        {
            return MAX_MBS_LEVEL_40;
        }
        case H264_LEVEL_4_1:
        {
            return MAX_MBS_LEVEL_41;
        }
        case H264_LEVEL_4_2:
        {
            return MAX_MBS_LEVEL_42;
        }
        case H264_LEVEL_5_0:
        {
            return MAX_MBS_LEVEL_50;
        }
        case H264_LEVEL_5_1:
        default:
        {
            return MAX_MBS_LEVEL_51;
        }
    }
}

WORD16 imvcd_free_dynamic_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_deblk_pic);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_dec_mb_map);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_recon_mb_map);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu2_slice_num_map);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_dec_slice_buf);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_frm_mb_info);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pi2_coeff_data);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_parse_mb_data);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_parse_part_params);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_deblk_top_mb);

    if(ps_view_ctxt->p_ctxt_inc_mb_map)
    {
        ps_view_ctxt->p_ctxt_inc_mb_map -= 1;
        PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->p_ctxt_inc_mb_map);
    }

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_mv_p[0]);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_mv_p[1]);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_pred_pkd);
    {
        UWORD8 i;
        for(i = 0; i < MV_SCRATCH_BUFS; i++)
        {
            PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_mv_top_p[i]);
        }
    }

    if(ps_view_ctxt->pu1_y_intra_pred_line)
    {
        ps_view_ctxt->pu1_y_intra_pred_line -= MB_SIZE;
    }
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_y_intra_pred_line);

    if(ps_view_ctxt->pu1_u_intra_pred_line)
    {
        ps_view_ctxt->pu1_u_intra_pred_line -= MB_SIZE;
    }
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_u_intra_pred_line);

    if(ps_view_ctxt->pu1_v_intra_pred_line)
    {
        ps_view_ctxt->pu1_v_intra_pred_line -= MB_SIZE;
    }
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_v_intra_pred_line);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->ps_nbr_mb_row);

    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_mvcd_ctxt->s_mvc_au_buf_mgr.pv_au_buf_base);
    PS_DEC_ALIGNED_FREE(ps_view_ctxt,
                        ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.pv_au_mv_pred_buf_base);

    return OK;
}

static UWORD32 imvcd_get_num_au_data_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    return ps_view_ctxt->u1_pic_bufs;
}

static UWORD32 imvcd_get_num_elements_in_mv_pred_buf(UWORD32 u4_view_wd, UWORD32 u4_view_ht)
{
    return (u4_view_wd * (u4_view_ht + PAD_MV_BANK_ROW)) / MB_SIZE;
}

static UWORD32 imvcd_get_mv_pred_buf_padding_length(UWORD32 u4_view_wd)
{
    return (u4_view_wd * OFFSET_MV_BANK_ROW) / MB_SIZE;
}

static UWORD32 imvcd_get_au_mv_pred_buf_size(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    UWORD32 u4_num_bufs = imvcd_get_num_au_data_bufs(ps_mvcd_ctxt);

    UWORD32 u4_size = 0;

    u4_size += sizeof(mvc_au_mv_pred_t);

    u4_size +=
        imvcd_get_num_elements_in_mv_pred_buf(ps_view_ctxt->u2_pic_wd, ps_view_ctxt->u2_pic_ht) *
        (sizeof(mv_pred_t) + sizeof(UWORD8));

    u4_size *= u4_num_bufs;
    u4_size *= ps_mvcd_ctxt->u2_num_views;

    return u4_size;
}

static UWORD32 imvcd_get_au_buf_size(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    UWORD32 u4_size = 0;
    UWORD32 u4_num_bufs = imvcd_get_num_au_data_bufs(ps_mvcd_ctxt);

    u4_size += sizeof(mvc_au_buffer_t);

    /* All rvalues below incorporate both padding and pic dimensions */
    u4_size += ALIGN64(ps_view_ctxt->u2_frm_wd_y * ps_view_ctxt->u2_frm_ht_y) * sizeof(UWORD8);
    u4_size += ALIGN64(ps_view_ctxt->u2_frm_wd_uv * ps_view_ctxt->u2_frm_ht_uv) * sizeof(UWORD8);

    u4_size *= ps_mvcd_ctxt->u2_num_views;
    u4_size *= u4_num_bufs;

    return u4_size;
}

WORD32 imvcd_init_au_buffers(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    UWORD32 i, j;
    UWORD32 u4_luma_size, u4_chroma_size;
    WORD32 i4_error_code;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    WORD64 i8_alloc_mem_size = imvcd_get_au_buf_size(ps_mvcd_ctxt);
    UWORD8 *pu1_buf = (UWORD8 *) ps_mvcd_ctxt->s_mvc_au_buf_mgr.pv_au_buf_base;
    UWORD32 u4_num_bufs = imvcd_get_num_au_data_bufs(ps_mvcd_ctxt);

    if(ps_mvcd_ctxt->u2_num_views > MAX_NUM_VIEWS)
    {
        ps_view_ctxt->i4_error_code = ERROR_BUF_MGR;
        return ERROR_BUF_MGR;
    }

    u4_luma_size = ps_view_ctxt->u2_frm_wd_y * ps_view_ctxt->u2_frm_ht_y;
    u4_chroma_size = ps_view_ctxt->u2_frm_wd_uv * ps_view_ctxt->u2_frm_ht_uv;

    for(i = 0; i < u4_num_bufs; i++)
    {
        WORD32 i4_stride;

        mvc_au_buffer_t *ps_au_buf = (mvc_au_buffer_t *) pu1_buf;

        pu1_buf += sizeof(ps_au_buf[0]);

        for(j = 0; j < ps_mvcd_ctxt->u2_num_views; j++)
        {
            i4_stride = ps_view_ctxt->u2_frm_wd_y;
            ps_au_buf->as_view_buffers[j].as_component_bufs[Y].i4_data_stride = i4_stride;
            ps_au_buf->as_view_buffers[j].as_component_bufs[Y].pv_data =
                pu1_buf + imvcd_get_ref_pic_pad_offset(i4_stride, false);
            pu1_buf += ALIGN64(u4_luma_size) * sizeof(pu1_buf[0]);
            i8_alloc_mem_size -= ALIGN64(u4_luma_size) * sizeof(pu1_buf[0]);

            i4_stride = ps_view_ctxt->u2_frm_wd_uv;
            ps_au_buf->as_view_buffers[j].as_component_bufs[UV].i4_data_stride = i4_stride;
            ps_au_buf->as_view_buffers[j].as_component_bufs[UV].pv_data =
                pu1_buf + imvcd_get_ref_pic_pad_offset(i4_stride, true);
            pu1_buf += ALIGN64(u4_chroma_size) * sizeof(pu1_buf[0]);
            i8_alloc_mem_size -= ALIGN64(u4_chroma_size) * sizeof(pu1_buf[0]);

            ps_au_buf->as_view_buffers[j].as_component_bufs[V].pv_data = NULL;

            ps_au_buf->as_view_buffers[j].u2_height =
                ps_view_ctxt->ps_cur_sps->u2_frm_ht_in_mbs * MB_SIZE;
            ps_au_buf->as_view_buffers[j].u2_width =
                ps_view_ctxt->ps_cur_sps->u2_frm_wd_in_mbs * MB_SIZE;
            ps_au_buf->as_view_buffers[j].u1_bit_depth = 8;

            ASSERT(i8_alloc_mem_size >= 0);
        }

        ps_au_buf->i4_pic_buf_id = i;

        i4_error_code =
            ih264_buf_mgr_add(ps_mvcd_ctxt->s_mvc_au_buf_mgr.ps_buf_mgr_ctxt, ps_au_buf, i);

        if(0 != i4_error_code)
        {
            ps_view_ctxt->i4_error_code = ERROR_BUF_MGR;

            return ERROR_BUF_MGR;
        }

        ps_mvcd_ctxt->s_mvc_au_buf_mgr.aps_buf_id_to_au_buf_map[i] = ps_au_buf;
    }

    return OK;
}

WORD32 imvcd_init_au_mv_pred_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    UWORD32 i, j;
    WORD32 buf_ret;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    UWORD32 u4_width = ps_view_ctxt->u2_pic_wd;
    UWORD32 u4_height = ps_view_ctxt->u2_pic_ht;
    UWORD32 u4_mode_info_buf_size = imvcd_get_num_elements_in_mv_pred_buf(u4_width, u4_height);
    UWORD8 *pu1_buf = ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.pv_au_mv_pred_buf_base;
    WORD64 i8_alloc_mem_size = imvcd_get_au_mv_pred_buf_size(ps_mvcd_ctxt);
    UWORD32 u4_num_bufs = imvcd_get_num_au_data_bufs(ps_mvcd_ctxt);

    if(ps_mvcd_ctxt->u2_num_views > MAX_NUM_VIEWS)
    {
        return ERROR_BUF_MGR;
    }

    for(i = 0; i < u4_num_bufs; i++)
    {
        mvc_au_mv_pred_t *ps_au_mv_data = (mvc_au_mv_pred_t *) pu1_buf;

        pu1_buf += sizeof(ps_au_mv_data[0]);

        buf_ret = ih264_buf_mgr_add(ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.ps_buf_mgr_ctxt,
                                    ps_au_mv_data, i);

        if(0 != buf_ret)
        {
            return ERROR_BUF_MGR;
        }

        for(j = 0; j < ps_mvcd_ctxt->u2_num_views; j++)
        {
            UWORD32 u4_mv_buf_size = u4_mode_info_buf_size * sizeof(ps_au_mv_data->aps_mvs[j][0]);
            UWORD32 u4_mode_desc_buf_size =
                u4_mode_info_buf_size * sizeof(ps_au_mv_data->apu1_mode_descriptors[j][0]);

            ps_au_mv_data->aps_mvs[j] = (mv_pred_t *) pu1_buf;
            ps_au_mv_data->aps_mvs[j] += imvcd_get_mv_pred_buf_padding_length(u4_width);
            pu1_buf += u4_mv_buf_size;
            i8_alloc_mem_size -= u4_mv_buf_size;

            ps_au_mv_data->apu1_mode_descriptors[j] = pu1_buf;
            pu1_buf += u4_mode_desc_buf_size;
            i8_alloc_mem_size -= u4_mode_desc_buf_size;

            memset(ps_au_mv_data->aps_mvs[j] - imvcd_get_mv_pred_buf_padding_length(u4_width), 0,
                   u4_mv_buf_size);

            memset(ps_au_mv_data->apu1_mode_descriptors[j], 0, u4_mode_desc_buf_size);

            ASSERT(i8_alloc_mem_size >= 0);
        }
    }

    return OK;
}

WORD32 imvcd_allocate_dynamic_bufs(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;
    dec_seq_params_t *ps_sps = ps_view_ctxt->ps_cur_sps;

    UWORD32 u4_total_mbs = ps_sps->u2_total_num_of_mbs;
    UWORD32 u4_wd_mbs = ps_view_ctxt->u2_frm_wd_in_mbs;
    UWORD32 u4_ht_mbs = ps_view_ctxt->u2_frm_ht_in_mbs;
    const WORD32 i4_default_alignment = 128;
    void *pv_mem_ctxt = ps_view_ctxt->pv_mem_ctxt;

    UWORD8 *pu1_buf;
    WORD32 i4_mem_size;
    WORD32 i;
    void *pv_buf;
    WORD32 i4_num_entries;

    if(ps_mvcd_ctxt->u2_num_views > MAX_NUM_VIEWS)
    {
        return IV_FAIL;
    }

    i4_mem_size = u4_total_mbs * sizeof(ps_view_ctxt->pu1_dec_mb_map[0]);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pu1_dec_mb_map = pv_buf;

    i4_mem_size = u4_total_mbs * sizeof(ps_view_ctxt->pu1_recon_mb_map[0]);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pu1_recon_mb_map = pv_buf;

    i4_mem_size = u4_total_mbs * sizeof(ps_view_ctxt->pu2_slice_num_map[0]);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->pu2_slice_num_map = pv_buf;

    ps_view_ctxt->ps_parse_cur_slice = ps_view_ctxt->ps_dec_slice_buf;
    ps_view_ctxt->ps_decode_cur_slice = ps_view_ctxt->ps_dec_slice_buf;
    ps_view_ctxt->ps_computebs_cur_slice = ps_view_ctxt->ps_dec_slice_buf;
    ps_view_ctxt->ps_pred_start = ps_view_ctxt->ps_pred;

    i4_mem_size = sizeof(parse_pmbarams_t) * (ps_view_ctxt->u1_recon_mb_grp);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_parse_mb_data = pv_buf;

    i4_mem_size = sizeof(parse_part_params_t) * ((ps_view_ctxt->u1_recon_mb_grp) << 4);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_parse_part_params = pv_buf;

    i4_mem_size = (u4_wd_mbs * sizeof(deblkmb_neighbour_t));
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_deblk_top_mb = pv_buf;

    i4_mem_size = sizeof(ctxt_inc_mb_info_t) * (u4_wd_mbs + 2);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->p_ctxt_inc_mb_map = pv_buf;
    /* 0th entry of CtxtIncMbMap will be always be containing default values
     for CABAC context representing MB not available */
    ps_view_ctxt->p_ctxt_inc_mb_map += 1;

    i4_mem_size = sizeof(mv_pred_t) * ps_view_ctxt->u1_recon_mb_grp * 16;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_mv_p[0] = pv_buf;

    i4_mem_size = sizeof(mv_pred_t) * ps_view_ctxt->u1_recon_mb_grp * 16;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_mv_p[1] = pv_buf;

    for(i = 0; i < MV_SCRATCH_BUFS; i++)
    {
        i4_mem_size = (sizeof(mv_pred_t) * ps_view_ctxt->u1_recon_mb_grp * 4);
        pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, i4_mem_size);
        ps_view_ctxt->ps_mv_top_p[i] = pv_buf;
    }

    i4_mem_size = sizeof(UWORD8) * ((u4_wd_mbs + 2) * MB_SIZE) * 2;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_view_ctxt->pu1_y_intra_pred_line = pv_buf;
    memset(ps_view_ctxt->pu1_y_intra_pred_line, 0, i4_mem_size);
    ps_view_ctxt->pu1_y_intra_pred_line += MB_SIZE;

    i4_mem_size = sizeof(UWORD8) * ((u4_wd_mbs + 2) * MB_SIZE) * 2;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_view_ctxt->pu1_u_intra_pred_line = pv_buf;
    memset(ps_view_ctxt->pu1_u_intra_pred_line, 0, i4_mem_size);
    ps_view_ctxt->pu1_u_intra_pred_line += MB_SIZE;

    i4_mem_size = sizeof(UWORD8) * ((u4_wd_mbs + 2) * MB_SIZE) * 2;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_view_ctxt->pu1_v_intra_pred_line = pv_buf;
    memset(ps_view_ctxt->pu1_v_intra_pred_line, 0, i4_mem_size);
    ps_view_ctxt->pu1_v_intra_pred_line += MB_SIZE;

    if(ps_view_ctxt->u1_separate_parse)
    {
        /* Needs one extra row of info, to hold top row data */
        i4_mem_size = sizeof(mb_neigbour_params_t) * 2 * ((u4_wd_mbs + 2) * (u4_ht_mbs + 1));
    }
    else
    {
        i4_mem_size = sizeof(mb_neigbour_params_t) * 2 * (u4_wd_mbs + 2);
    }

    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);

    ps_view_ctxt->ps_nbr_mb_row = pv_buf;
    memset(ps_view_ctxt->ps_nbr_mb_row, 0, i4_mem_size);

    i4_mem_size = (u4_total_mbs + u4_wd_mbs) * sizeof(deblk_mb_t);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_view_ctxt->ps_deblk_pic = pv_buf;
    memset(ps_view_ctxt->ps_deblk_pic, 0, i4_mem_size);

    i4_mem_size = sizeof(dec_mb_info_t) * u4_total_mbs;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_view_ctxt->ps_frm_mb_info = pv_buf;
    memset(ps_view_ctxt->ps_frm_mb_info, 0, i4_mem_size);

    if((1 >= ps_view_ctxt->ps_cur_sps->u1_num_ref_frames) && (0 == ps_view_ctxt->i4_display_delay))
    {
        i4_num_entries = 1;
    }
    else
    {
        i4_num_entries = MAX_FRAMES;
    }

    i4_num_entries = (2 * i4_num_entries) + 1;
    i4_num_entries *= 2;

    i4_mem_size = i4_num_entries * sizeof(void *);
    i4_mem_size += PAD_MAP_IDX_POC * sizeof(void *);
    i4_mem_size *= u4_total_mbs;
    i4_mem_size += sizeof(dec_slice_struct_t) * u4_total_mbs;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);

    ps_view_ctxt->ps_dec_slice_buf = pv_buf;
    memset(ps_view_ctxt->ps_dec_slice_buf, 0, i4_mem_size);
    pu1_buf = (UWORD8 *) ps_view_ctxt->ps_dec_slice_buf;
    pu1_buf += sizeof(dec_slice_struct_t) * u4_total_mbs;
    ps_view_ctxt->pv_map_ref_idx_to_poc_buf = (void *) pu1_buf;

    /* Allocate memory for packed pred info */
    i4_num_entries = u4_total_mbs;
    i4_num_entries *= 16 * 2;

    i4_mem_size = sizeof(pred_info_pkd_t) * i4_num_entries;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_view_ctxt->ps_pred_pkd = pv_buf;

    /* Allocate memory for coeff data */
    i4_mem_size = MB_LUM_SIZE * sizeof(WORD16);
    /*For I16x16 MBs, 16 4x4 AC coeffs and 1 4x4 DC coeff TU blocks will be sent
    For all MBs along with 8 4x4 AC coeffs 2 2x2 DC coeff TU blocks will be sent
    So use 17 4x4 TU blocks for luma and 9 4x4 TU blocks for chroma */
    i4_mem_size += u4_total_mbs *
                   (MAX(17 * sizeof(tu_sblk4x4_coeff_data_t), 4 * sizeof(tu_blk8x8_coeff_data_t)) +
                    9 * sizeof(tu_sblk4x4_coeff_data_t));
    // 32 bytes for each mb to store u1_prev_intra4x4_pred_mode and
    // u1_rem_intra4x4_pred_mode data
    i4_mem_size += u4_total_mbs * 32;
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);

    ps_view_ctxt->pi2_coeff_data = pv_buf;

    ps_view_ctxt->pv_pic_tu_coeff_data = (void *) (ps_view_ctxt->pi2_coeff_data + MB_LUM_SIZE);

    i4_mem_size = imvcd_get_au_mv_pred_buf_size(ps_mvcd_ctxt);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_mvcd_ctxt->s_mvc_au_mv_pred_buf_mgr.pv_au_mv_pred_buf_base = pv_buf;

    i4_mem_size = imvcd_get_au_buf_size(ps_mvcd_ctxt);
    pv_buf = ps_view_ctxt->pf_aligned_alloc(pv_mem_ctxt, i4_default_alignment, i4_mem_size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, i4_mem_size);
    ps_mvcd_ctxt->s_mvc_au_buf_mgr.pv_au_buf_base = pv_buf;

    /***************************************************************************/
    /*Initialize cabac context pointers for every SE that has fixed contextIdx */
    /***************************************************************************/
    {
        bin_ctxt_model_t *const p_cabac_ctxt_table_t = ps_view_ctxt->p_cabac_ctxt_table_t;
        bin_ctxt_model_t **p_coeff_abs_level_minus1_t = ps_view_ctxt->p_coeff_abs_level_minus1_t;
        bin_ctxt_model_t **p_cbf_t = ps_view_ctxt->p_cbf_t;

        ps_view_ctxt->p_mb_field_dec_flag_t = p_cabac_ctxt_table_t + MB_FIELD_DECODING_FLAG;
        ps_view_ctxt->p_prev_intra4x4_pred_mode_flag_t =
            p_cabac_ctxt_table_t + PREV_INTRA4X4_PRED_MODE_FLAG;
        ps_view_ctxt->p_rem_intra4x4_pred_mode_t = p_cabac_ctxt_table_t + REM_INTRA4X4_PRED_MODE;
        ps_view_ctxt->p_intra_chroma_pred_mode_t = p_cabac_ctxt_table_t + INTRA_CHROMA_PRED_MODE;
        ps_view_ctxt->p_mb_qp_delta_t = p_cabac_ctxt_table_t + MB_QP_DELTA;
        ps_view_ctxt->p_ref_idx_t = p_cabac_ctxt_table_t + REF_IDX;
        ps_view_ctxt->p_mvd_x_t = p_cabac_ctxt_table_t + MVD_X;
        ps_view_ctxt->p_mvd_y_t = p_cabac_ctxt_table_t + MVD_Y;
        p_cbf_t[0] = p_cabac_ctxt_table_t + CBF + 0;
        p_cbf_t[1] = p_cabac_ctxt_table_t + CBF + 4;
        p_cbf_t[2] = p_cabac_ctxt_table_t + CBF + 8;
        p_cbf_t[3] = p_cabac_ctxt_table_t + CBF + 12;
        p_cbf_t[4] = p_cabac_ctxt_table_t + CBF + 16;
        ps_view_ctxt->p_cbp_luma_t = p_cabac_ctxt_table_t + CBP_LUMA;
        ps_view_ctxt->p_cbp_chroma_t = p_cabac_ctxt_table_t + CBP_CHROMA;

        p_coeff_abs_level_minus1_t[LUMA_DC_CTXCAT] =
            p_cabac_ctxt_table_t + COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_0_OFFSET;

        p_coeff_abs_level_minus1_t[LUMA_AC_CTXCAT] =
            p_cabac_ctxt_table_t + COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_1_OFFSET;

        p_coeff_abs_level_minus1_t[LUMA_4X4_CTXCAT] =
            p_cabac_ctxt_table_t + COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_2_OFFSET;

        p_coeff_abs_level_minus1_t[CHROMA_DC_CTXCAT] =
            p_cabac_ctxt_table_t + COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_3_OFFSET;

        p_coeff_abs_level_minus1_t[CHROMA_AC_CTXCAT] =
            p_cabac_ctxt_table_t + COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_4_OFFSET;

        p_coeff_abs_level_minus1_t[LUMA_8X8_CTXCAT] =
            p_cabac_ctxt_table_t + COEFF_ABS_LEVEL_MINUS1_8X8 + COEFF_ABS_LEVEL_CAT_5_OFFSET;

        /********************************************************/
        /* context for the high profile related syntax elements */
        /* This is maintained seperately in s_high_profile     */
        /********************************************************/
        {
            ps_view_ctxt->s_high_profile.ps_transform8x8_flag =
                p_cabac_ctxt_table_t + TRANSFORM_SIZE_8X8_FLAG;

            ps_view_ctxt->s_high_profile.ps_sigcoeff_8x8_frame =
                p_cabac_ctxt_table_t + SIGNIFICANT_COEFF_FLAG_8X8_FRAME;

            ps_view_ctxt->s_high_profile.ps_last_sigcoeff_8x8_frame =
                p_cabac_ctxt_table_t + LAST_SIGNIFICANT_COEFF_FLAG_8X8_FRAME;

            ps_view_ctxt->s_high_profile.ps_coeff_abs_levelminus1 =
                p_cabac_ctxt_table_t + COEFF_ABS_LEVEL_MINUS1_8X8;

            ps_view_ctxt->s_high_profile.ps_sigcoeff_8x8_field =
                p_cabac_ctxt_table_t + SIGNIFICANT_COEFF_FLAG_8X8_FIELD;

            ps_view_ctxt->s_high_profile.ps_last_sigcoeff_8x8_field =
                p_cabac_ctxt_table_t + LAST_SIGNIFICANT_COEFF_FLAG_8X8_FIELD;
        }
    }

    return OK;
}

void imvcd_convert_au_buf_to_view_buf(mvc_au_buffer_t *ps_au_buf, pic_buffer_t *ps_view_buf,
                                      UWORD16 u2_view_order_id, UWORD16 u2_view_id)
{
    yuv_buf_props_t *ps_view_buffer = &ps_au_buf->as_view_buffers[u2_view_id];
    offsets_t *ps_disp_offsets = &ps_au_buf->as_disp_offsets[u2_view_id];
    mvc_au_mv_pred_t *ps_au_mv_data = ps_au_buf->ps_au_mv_data;

    ps_view_buf->pu1_buf1 = ps_view_buffer->as_component_bufs[Y].pv_data;
    ps_view_buf->pu1_buf2 = ps_view_buffer->as_component_bufs[UV].pv_data;
    ps_view_buf->pu1_buf3 = NULL;
    ps_view_buf->u2_frm_wd_y = ps_view_buffer->as_component_bufs[Y].i4_data_stride;
    ps_view_buf->u2_frm_wd_uv = ps_view_buffer->as_component_bufs[UV].i4_data_stride;
    ps_view_buf->pu1_buf3 = NULL;

    ps_view_buf->u2_disp_width = ps_au_buf->u2_disp_width;
    ps_view_buf->u2_disp_height = ps_au_buf->u2_disp_height;
    ps_view_buf->u2_frm_ht_y = ps_view_buffer->u2_height;
    ps_view_buf->u2_frm_ht_uv = ps_view_buffer->u2_height / 2;

    ps_view_buf->u4_time_stamp = ps_au_buf->u4_time_stamp;
    ps_view_buf->u4_ts = ps_au_buf->u4_time_stamp;

    ps_view_buf->u2_crop_offset_y =
        ps_disp_offsets->u2_left_offset + ps_disp_offsets->u2_top_offset * ps_view_buf->u2_frm_wd_y;
    ps_view_buf->u2_crop_offset_uv =
        ps_disp_offsets->u2_left_offset +
        (ps_disp_offsets->u2_top_offset / 2) * ps_view_buf->u2_frm_wd_uv;

    ps_view_buf->i4_poc = ps_au_buf->i4_poc;
    ps_view_buf->i4_pic_num = ps_au_buf->i4_frame_num;
    ps_view_buf->i4_frame_num = ps_au_buf->i4_frame_num;
    ps_view_buf->i4_avg_poc = ps_au_buf->i4_poc;
    ps_view_buf->u1_is_short = ps_au_buf->b_is_short_term_ref;
    ps_view_buf->u1_pic_type = ps_au_buf->u1_pic_type;
    ps_view_buf->i4_top_field_order_cnt = ps_au_buf->i4_poc;
    ps_view_buf->i4_bottom_field_order_cnt = ps_au_buf->i4_poc;
    ps_view_buf->u1_picturetype = FRM_PIC;
    ps_view_buf->u1_long_term_frm_idx = ps_au_buf->u1_long_term_frm_idx;
    ps_view_buf->u1_long_term_pic_num = ps_au_buf->u1_long_term_pic_num;
    ps_view_buf->u4_pack_slc_typ = ps_au_buf->au4_pack_slc_typ[u2_view_order_id];
    ps_view_buf->u1_pic_struct = ps_au_buf->u1_pic_struct;
    ps_view_buf->s_sei_pic = ps_au_buf->s_sei_pic;

    ps_view_buf->u1_pic_buf_id = ps_au_buf->i4_pic_buf_id;
    ps_view_buf->u1_mv_buf_id = ps_au_buf->i4_mv_buf_id;

    ps_view_buf->pu1_col_zero_flag = ps_au_mv_data->apu1_mode_descriptors[u2_view_id];
    ps_view_buf->ps_mv = ps_au_mv_data->aps_mvs[u2_view_id];
}

void imvcd_init_ref_idx_to_ref_buf_map(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    pic_buffer_t *ps_pic;

    void **ppv_map_ref_idx_to_poc_lx;
    WORD8 i, j;

    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    bool b_is_b_pic = !!(
        ps_mvcd_ctxt->ps_cur_au->au4_pack_slc_typ[ps_mvcd_ctxt->u2_num_views_decoded] & B_SLC_BIT);

    for(i = 0; i < 1 + ((WORD32) b_is_b_pic); i++)
    {
        ppv_map_ref_idx_to_poc_lx =
            ps_view_ctxt->ppv_map_ref_idx_to_poc + ((0 == i) ? FRM_LIST_L0 : FRM_LIST_L1);
        ppv_map_ref_idx_to_poc_lx[0] = NULL;
        ppv_map_ref_idx_to_poc_lx++;

        for(j = 0; j < ps_view_ctxt->ps_cur_slice->u1_num_ref_idx_lx_active[i]; j++)
        {
            ps_pic = ps_view_ctxt->ps_ref_pic_buf_lx[i][j];

            ppv_map_ref_idx_to_poc_lx[j] = ps_pic->pu1_buf1;
        }
    }

    if(!b_is_b_pic)
    {
        ppv_map_ref_idx_to_poc_lx = ps_view_ctxt->ppv_map_ref_idx_to_poc + FRM_LIST_L1;
        ppv_map_ref_idx_to_poc_lx[0] = NULL;
    }

    if(ps_view_ctxt->u4_num_cores >= 3)
    {
        WORD32 i4_size;

        WORD32 i4_num_entries = MAX_FRAMES;

        if((1 >= ps_view_ctxt->ps_cur_sps->u1_num_ref_frames) &&
           (0 == ps_view_ctxt->i4_display_delay))
        {
            i4_num_entries = 1;
        }

        i4_num_entries = 2 * i4_num_entries + 1;
        i4_num_entries *= 2;

        i4_size = i4_num_entries * sizeof(void *);
        i4_size += PAD_MAP_IDX_POC * sizeof(void *);

        memcpy(ps_view_ctxt->ps_parse_cur_slice->ppv_map_ref_idx_to_poc,
               ps_view_ctxt->ppv_map_ref_idx_to_poc, i4_size);
    }
}

void imvcd_ivp_buf_copier(mvc_au_buffer_t *ps_au_buf_src, mvc_au_buffer_t *ps_au_buf_dst,
                          mvc_au_mv_pred_t *ps_au_mv_data_src, mvc_au_mv_pred_t *ps_au_mv_data_dst,
                          UWORD16 u2_src_view_id, UWORD16 u2_dst_view_id)
{
    UWORD32 i, j;

    mv_pred_t *ps_mode_info_src = ps_au_mv_data_src->aps_mvs[u2_src_view_id];
    mv_pred_t *ps_mode_info_dst = ps_au_mv_data_dst->aps_mvs[u2_dst_view_id];

    UWORD32 u4_view_wd = ps_au_buf_src->as_view_buffers[u2_src_view_id].u2_width;
    UWORD32 u4_view_ht = ps_au_buf_src->as_view_buffers[u2_src_view_id].u2_height;
    UWORD32 u4_mode_info_buf_size = imvcd_get_num_elements_in_mv_pred_buf(u4_view_wd, u4_view_ht);
    UWORD32 u4_mode_info_pad_size = imvcd_get_mv_pred_buf_padding_length(u4_view_wd);

    ps_mode_info_src -= u4_mode_info_pad_size;
    ps_mode_info_dst -= u4_mode_info_pad_size;

    ps_au_buf_dst->ps_au_mv_data = ps_au_mv_data_dst;

    ps_au_buf_dst->as_disp_offsets[u2_dst_view_id] = ps_au_buf_src->as_disp_offsets[u2_src_view_id];

    for(i = 0; i < NUM_SP_COMPONENTS; i++)
    {
        bool b_is_chroma = ((COMPONENT_TYPES_T) i) != Y;

        coordinates_t s_pad_dims = imvcd_get_buf_pad_dims(b_is_chroma);
        buffer_container_t *ps_src =
            &ps_au_buf_src->as_view_buffers[u2_src_view_id].as_component_bufs[i];
        buffer_container_t *ps_dst =
            &ps_au_buf_dst->as_view_buffers[u2_dst_view_id].as_component_bufs[i];

        WORD32 i4_src_pad_offset =
            imvcd_get_ref_pic_pad_offset(ps_src->i4_data_stride, b_is_chroma);
        WORD32 i4_dst_pad_offset =
            imvcd_get_ref_pic_pad_offset(ps_dst->i4_data_stride, b_is_chroma);

        for(j = 0; j < ((u4_view_ht >> b_is_chroma) + s_pad_dims.i4_ordinate); j++)
        {
            UWORD8 *pu1_src =
                ((UWORD8 *) ps_src->pv_data) + j * ps_src->i4_data_stride - i4_src_pad_offset;
            UWORD8 *pu1_dst =
                ((UWORD8 *) ps_dst->pv_data) + j * ps_dst->i4_data_stride - i4_dst_pad_offset;

            memcpy(pu1_dst, pu1_src, (u4_view_wd + s_pad_dims.i4_abscissa) * sizeof(pu1_dst[0]));
        }
    }

    memcpy(ps_mode_info_dst, ps_mode_info_src, u4_mode_info_buf_size * sizeof(ps_mode_info_dst[0]));

    for(i = 0; i < u4_mode_info_buf_size; i++)
    {
        /* In accordance with 'H.8.4' */
        ps_au_mv_data_dst->apu1_mode_descriptors[u2_dst_view_id][i] =
            ps_au_mv_data_src->apu1_mode_descriptors[u2_src_view_id][i] & 0xFE;
    }

    ps_au_buf_dst->au4_pack_slc_typ[u2_dst_view_id] =
        ps_au_buf_src->au4_pack_slc_typ[u2_src_view_id];
    ps_au_buf_dst->b_is_short_term_ref = ps_au_buf_src->b_is_short_term_ref;
    ps_au_buf_dst->i4_avg_poc = ps_au_buf_src->i4_avg_poc;
    ps_au_buf_dst->i4_frame_num = ps_au_buf_src->i4_frame_num;
    ps_au_buf_dst->i4_pic_num = ps_au_buf_src->i4_pic_num;
    ps_au_buf_dst->i4_poc = ps_au_buf_src->i4_poc;
    ps_au_buf_dst->s_sei_pic = ps_au_buf_src->s_sei_pic;
    ps_au_buf_dst->u1_long_term_frm_idx = ps_au_buf_src->u1_long_term_frm_idx;
    ps_au_buf_dst->u1_long_term_pic_num = ps_au_buf_src->u1_long_term_pic_num;
    ps_au_buf_dst->u1_picturetype = ps_au_buf_src->u1_picturetype;
    ps_au_buf_dst->u1_pic_struct = ps_au_buf_src->u1_pic_struct;
    ps_au_buf_dst->u2_disp_height = ps_au_buf_src->u2_disp_height;
    ps_au_buf_dst->u2_disp_width = ps_au_buf_src->u2_disp_width;
    ps_au_buf_dst->u4_time_stamp = ps_au_buf_src->u4_time_stamp;
}
