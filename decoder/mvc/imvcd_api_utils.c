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
/*  File Name         : imvcd_api_utils.c                                    */
/*                                                                           */
/*  Description       : Utility functions used by 'imvcd_api.c'              */
/*                                                                           */
/*****************************************************************************/
#include <string.h>

#include "ih264_typedefs.h"
#include "iv.h"
#include "imvcd.h"
#include "ih264_disp_mgr.h"
#include "ih264d_structs.h"
#include "ih264d_tables.h"
#include "ih264d_utils.h"
#include "imvcd_structs.h"
#include "imvcd_utils.h"

IV_API_CALL_STATUS_T imvcd_check_dec_handle(iv_obj_t *ps_handle)
{
    if(ps_handle == NULL)
    {
        return IV_FAIL;
    }

    if(ps_handle->pv_codec_handle == NULL)
    {
        return IV_FAIL;
    }

    return IV_SUCCESS;
}

IV_API_CALL_STATUS_T imvcd_check_create_structs(imvcd_create_ip_t *ps_ip, imvcd_create_op_t *ps_op)
{
    if(NULL == ps_ip)
    {
        return IV_FAIL;
    }

    if(NULL == ps_op)
    {
        return IV_FAIL;
    }

    if(ps_ip->s_ivd_ip.e_output_format != IV_YUV_420P)
    {
        return IV_FAIL;
    }

    if(NULL == ps_ip->s_ivd_ip.pf_aligned_alloc)
    {
        return IV_FAIL;
    }

    if(NULL == ps_ip->s_ivd_ip.pf_aligned_free)
    {
        return IV_FAIL;
    }

    if(0 != ps_ip->s_ivd_ip.u4_share_disp_buf)
    {
        return IV_FAIL;
    }

    return IV_SUCCESS;
}

IV_API_CALL_STATUS_T imvcd_check_ctl_structs(void *pv_ip, void *pv_op)
{
    ivd_ctl_set_config_ip_t *ps_ip = (ivd_ctl_set_config_ip_t *) pv_ip;
    ivd_ctl_set_config_op_t *ps_op = (ivd_ctl_set_config_op_t *) pv_op;

    WORD32 i4_sub_cmd = ps_ip->e_sub_cmd;

    if(NULL == ps_ip)
    {
        return IV_FAIL;
    }

    if(NULL == ps_op)
    {
        return IV_FAIL;
    }

    switch(i4_sub_cmd)
    {
        case IVD_CMD_CTL_SETPARAMS:
        {
            if((ps_ip->u4_size != sizeof(ivd_ctl_set_config_ip_t)) ||
               (ps_op->u4_size != sizeof(ivd_ctl_set_config_op_t)))
            {
                return IV_FAIL;
            }
            else
            {
                return IV_SUCCESS;
            }
        }
        case IMVCD_CTL_SET_NUM_CORES:
        {
            if((ps_ip->u4_size != sizeof(imvcd_set_num_cores_ip_t)) ||
               (ps_op->u4_size != sizeof(imvcd_set_num_cores_op_t)))
            {
                return IV_FAIL;
            }
            else
            {
                return IV_SUCCESS;
            }
        }
        case IMVCD_CTL_SET_PROCESSOR:
        {
            if((ps_ip->u4_size != sizeof(imvcd_set_arch_ip_t)) ||
               (ps_op->u4_size != sizeof(imvcd_set_arch_op_t)))
            {
                return IV_FAIL;
            }
            else
            {
                return IV_SUCCESS;
            }
        }
        case IMVCD_CTL_DEGRADE:
        {
            if((ps_ip->u4_size != sizeof(imvcd_set_degrade_mode_ip_t)) ||
               (ps_op->u4_size != sizeof(imvcd_set_degrade_mode_op_t)))
            {
                return IV_FAIL;
            }
            else
            {
                return IV_SUCCESS;
            }
        }
        case IVD_CMD_CTL_FLUSH:
        {
            if((ps_ip->u4_size != sizeof(ivd_ctl_flush_ip_t)) ||
               (ps_op->u4_size != sizeof(ivd_ctl_flush_op_t)))
            {
                return IV_FAIL;
            }
            else
            {
                return IV_SUCCESS;
            }
        }
        case IVD_CMD_CTL_GETBUFINFO:
        {
            if((ps_ip->u4_size != sizeof(ivd_ctl_getbufinfo_ip_t)) ||
               (ps_op->u4_size != sizeof(ivd_ctl_getbufinfo_op_t)))
            {
                return IV_FAIL;
            }
            else
            {
                return IV_SUCCESS;
            }
        }
        case IMVCD_CTL_GET_VUI_PARAMS:
        {
            if((ps_ip->u4_size != sizeof(imvcd_get_vui_ip_t)) ||
               (ps_op->u4_size != sizeof(imvcd_get_vui_op_t)))
            {
                return IV_FAIL;
            }
            else
            {
                return IV_SUCCESS;
            }
        }
        default:
        {
            return IV_FAIL;
        }
    }
}

IV_API_CALL_STATUS_T imvcd_check_decode_structs(iv_obj_t *ps_dec_hdl,
                                                imvcd_video_decode_ip_t *ps_ip,
                                                imvcd_video_decode_op_t *ps_op)
{
    WORD32 i, j;

    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;

    UWORD16 u2_num_views = 1;

    if(NULL == ps_ip)
    {
        return IV_FAIL;
    }

    if(NULL == ps_op)
    {
        return IV_FAIL;
    }

    if(!ps_mvcd_ctxt->b_flush_enabled && !ps_mvcd_ctxt->b_header_only_decode)
    {
        if(NULL == ps_ip->s_ivd_ip.pv_stream_buffer)
        {
            return IV_FAIL;
        }

        for(i = 0; i < u2_num_views; i++)
        {
            for(j = 0; j < NUM_COMPONENTS; j++)
            {
                if(NULL == ps_ip->s_ivd_ip.s_out_buffer.pu1_bufs[i * NUM_COMPONENTS + j])
                {
                    return IV_FAIL;
                }
            }
        }

        if(0 == ps_ip->s_ivd_ip.u4_num_Bytes)
        {
            return IV_FAIL;
        }
    }

    return IV_SUCCESS;
}

static void imvcd_convert_app_out_buf(mvc_dec_ctxt_t *ps_mvcd_ctxt,
                                      ivd_out_bufdesc_t *ps_app_buffer)
{
    if(!ps_mvcd_ctxt->b_header_only_decode)
    {
        WORD32 i, j;

        subset_sps_t *ps_subset_sps = imvcd_get_valid_subset_sps(ps_mvcd_ctxt);
        dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

        UWORD16 u2_num_views =
            (NULL == ps_subset_sps) ? 1 : ps_subset_sps->s_sps_mvc_ext.u2_num_views;

        for(i = 0; i < u2_num_views; i++)
        {
            yuv_buf_props_t *ps_view_buf = &ps_mvcd_ctxt->s_out_buffer.as_view_buf_props[i];

            ps_view_buf->u1_bit_depth = 8;
            ps_view_buf->u2_height = ps_view_ctxt->u2_disp_height;
            ps_view_buf->u2_width = ps_view_ctxt->u2_disp_width;

            for(j = 0; j < NUM_COMPONENTS; j++)
            {
                buffer_container_t *ps_component_buf = &ps_view_buf->as_component_bufs[j];
                bool b_is_chroma = (((COMPONENT_TYPES_T) j) != Y);

                ps_component_buf->pv_data = ps_app_buffer->pu1_bufs[i * NUM_COMPONENTS + j];
                ps_component_buf->i4_data_stride = ps_view_buf->u2_width >> b_is_chroma;
            }
        }
    }
}

void imvcd_convert_to_app_disp_buf(mvc_dec_ctxt_t *ps_mvcd_ctxt, iv_yuv_buf_t *ps_view_disp_bufs)
{
    WORD32 i;

    UWORD16 u2_num_views = ps_mvcd_ctxt->u2_num_views;

    for(i = 0; i < u2_num_views; i++)
    {
        yuv_buf_props_t *ps_view_buf = &ps_mvcd_ctxt->s_out_buffer.as_view_buf_props[i];

        ps_view_disp_bufs[i].u4_size = sizeof(ps_view_disp_bufs[i]);

        ps_view_disp_bufs[i].pv_y_buf = ps_view_buf->as_component_bufs[Y].pv_data;
        ps_view_disp_bufs[i].u4_y_strd = ps_view_buf->as_component_bufs[Y].i4_data_stride;
        ps_view_disp_bufs[i].u4_y_wd = ps_view_buf->u2_width;
        ps_view_disp_bufs[i].u4_y_ht = ps_view_buf->u2_height;

        ps_view_disp_bufs[i].pv_u_buf = ps_view_buf->as_component_bufs[U].pv_data;
        ps_view_disp_bufs[i].u4_u_strd = ps_view_buf->as_component_bufs[U].i4_data_stride;
        ps_view_disp_bufs[i].u4_u_wd = ps_view_buf->u2_width / 2;
        ps_view_disp_bufs[i].u4_u_ht = ps_view_buf->u2_height / 2;

        ps_view_disp_bufs[i].pv_v_buf = ps_view_buf->as_component_bufs[V].pv_data;
        ps_view_disp_bufs[i].u4_v_strd = ps_view_buf->as_component_bufs[V].i4_data_stride;
        ps_view_disp_bufs[i].u4_v_wd = ps_view_buf->u2_width / 2;
        ps_view_disp_bufs[i].u4_v_ht = ps_view_buf->u2_height / 2;
    }
}

void imvcd_au_init(iv_obj_t *ps_dec_hdl, imvcd_video_decode_ip_t *ps_ip,
                   imvcd_video_decode_op_t *ps_op)
{
    subset_sps_t *ps_subset_sps;

    mvc_dec_ctxt_t *ps_mvcd_ctxt = (mvc_dec_ctxt_t *) ps_dec_hdl->pv_codec_handle;
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_mvcd_ctxt->u2_num_views_decoded = 0;
    ps_subset_sps = imvcd_get_valid_subset_sps(ps_mvcd_ctxt);
    ps_mvcd_ctxt->u2_num_views =
        (NULL == ps_subset_sps) ? 1 : ps_subset_sps->s_sps_mvc_ext.u2_num_views;
    imvcd_convert_app_out_buf(ps_mvcd_ctxt, &ps_ip->s_ivd_ip.s_out_buffer);

    ps_op->s_ivd_op.u4_num_bytes_consumed = 0;
    ps_op->s_ivd_op.u4_output_present = 0;
    ps_op->s_ivd_op.u4_error_code = 0;
    ps_op->s_ivd_op.e_pic_type = IV_FRAMETYPE_DEFAULT;
    ps_op->s_ivd_op.u4_frame_decoded_flag = 0;
    ps_op->s_ivd_op.u4_new_seq = 0;
    ps_op->s_ivd_op.u4_progressive_frame_flag = 1;
    ps_op->s_ivd_op.u4_is_ref_flag = 1;
    ps_op->s_ivd_op.e4_fld_type = IV_NA_FLD;

    ps_view_ctxt->u4_fmt_conv_cur_row = 0;
    ps_view_ctxt->u4_output_present = 0;
    ps_view_ctxt->u4_fmt_conv_num_rows = FMT_CONV_NUM_ROWS;
    ps_view_ctxt->u4_ts = ps_ip->s_ivd_ip.u4_ts;
    ps_view_ctxt->i4_frametype = IV_NA_FRAME;
    ps_view_ctxt->i4_content_type = IV_CONTENTTYPE_NA;

    /* Mimicking the hack in lines '2005' in 'ih264d_api.c' and '1323' in
     * 'ih264d_parse_headers.c */
    ps_view_ctxt->u4_sps_cnt_in_process = 0;

    memset(ps_mvcd_ctxt->as_nalu_mvc_ext, 0, sizeof(ps_mvcd_ctxt->as_nalu_mvc_ext));
}

void imvcd_view_init(mvc_dec_ctxt_t *ps_mvcd_ctxt)
{
    dec_struct_t *ps_view_ctxt = &ps_mvcd_ctxt->s_view_dec_ctxt;

    ps_view_ctxt->u4_num_fld_in_frm = 0;
    ps_view_ctxt->u4_slice_start_code_found = 0;
    ps_view_ctxt->u4_cur_mb_addr = 0;
    ps_view_ctxt->u4_total_mbs_coded = 0;
    ps_view_ctxt->u2_cur_slice_num = 0;
    ps_view_ctxt->cur_dec_mb_num = 0;
    ps_view_ctxt->cur_recon_mb_num = 0;
    ps_view_ctxt->u4_first_slice_in_pic = 1;
    ps_view_ctxt->u1_slice_header_done = 0;
    ps_view_ctxt->u1_dangling_field = 0;
    ps_view_ctxt->u4_dec_thread_created = 0;
    ps_view_ctxt->u4_bs_deblk_thread_created = 0;
    ps_view_ctxt->u4_cur_bs_mb_num = 0;
    ps_view_ctxt->u4_start_recon_deblk = 0;
    ps_view_ctxt->u4_pic_buf_got = 0;
    ps_view_ctxt->pu1_inv_scan = (UWORD8 *) gau1_ih264d_inv_scan;
    ps_view_ctxt->u1_pic_decode_done = 0;

    ps_view_ctxt->pu1_init_dpb_base = NULL;
    ps_view_ctxt->ps_dpb_mgr = NULL;
}

IV_API_CALL_STATUS_T imvcd_bitstream_buf_alloc(dec_struct_t *ps_view_ctxt, UWORD16 u2_num_views)
{
    UWORD32 u4_size;

    u4_size = MAX(MIN_BITSTREAMS_BUF_SIZE,
                  ((ps_view_ctxt->u2_pic_wd * ps_view_ctxt->u2_pic_ht * 3 / 2) + EXTRA_BS_OFFSET) *
                      u2_num_views * sizeof(ps_view_ctxt->pu1_bits_buf_dynamic[0]));
    ps_view_ctxt->pu1_bits_buf_dynamic =
        ps_view_ctxt->pf_aligned_alloc(ps_view_ctxt->pv_mem_ctxt, 128, u4_size);
    RETURN_IF((NULL == ps_view_ctxt->pu1_bits_buf_dynamic), IV_FAIL);

    memset(ps_view_ctxt->pu1_bits_buf_dynamic, 0, u4_size);
    ps_view_ctxt->u4_dynamic_bits_buf_size = u4_size;

    return IV_SUCCESS;
}

void imvcd_bitsteam_buf_free(dec_struct_t *ps_view_ctxt)
{
    PS_DEC_ALIGNED_FREE(ps_view_ctxt, ps_view_ctxt->pu1_bits_buf_dynamic);
}

IV_API_CALL_STATUS_T imvcd_bitstream_buf_realloc(dec_struct_t *ps_view_ctxt, UWORD32 u4_size)
{
    imvcd_bitsteam_buf_free(ps_view_ctxt);

    u4_size = MAX(MIN_BITSTREAMS_BUF_SIZE, u4_size);

    ps_view_ctxt->pu1_bits_buf_dynamic =
        ps_view_ctxt->pf_aligned_alloc(ps_view_ctxt->pv_mem_ctxt, 128, u4_size);
    RETURN_IF((NULL == ps_view_ctxt->pu1_bits_buf_dynamic), IV_FAIL);

    memset(ps_view_ctxt->pu1_bits_buf_dynamic, 0, u4_size);
    ps_view_ctxt->u4_dynamic_bits_buf_size = u4_size;

    return IV_SUCCESS;
}
