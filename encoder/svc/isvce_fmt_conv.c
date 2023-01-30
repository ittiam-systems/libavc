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
*  isvce_fmt_conv.c
*
* @brief
*  Contains functions for format conversion or frame copy of output buffer
*
* @author
*  ittiam
*
* @par List of Functions:
*  - isvce_fmt_conv()
*
* @remarks
*  None
*
*******************************************************************************
*/
#include "ih264_typedefs.h"
#include "ih264_macros.h"
/* Dependencies of ih264_buf_mgr.h */
/* Dependencies of ih264_list.h */
#include "ih264_error.h"
/* Dependencies of ih264_common_tables.h */
#include "ih264_defs.h"
#include "ih264_structs.h"
#include "ih264_buf_mgr.h"
#include "ih264_common_tables.h"
#include "ih264_list.h"
#include "ih264_platform_macros.h"
#include "ih264_trans_data.h"
#include "ih264_size_defs.h"
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
/* Dependencies of ih264e_bitstream.h */
#include "ih264e_error.h"
#include "ih264e_bitstream.h"
#include "ih264e_cabac_structs.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "ime_statistics.h"
#include "ime_structs.h"
/* Dependencies of 'ih264e_utils.h' */
#include "ih264e_defs.h"
#include "ih264e_structs.h"
#include "ih264e_fmt_conv.h"
#include "isvce_structs.h"

IH264E_ERROR_T isvce_fmt_conv(isvce_codec_t *ps_codec, svc_au_buf_t *ps_pic, UWORD8 *pu1_y_dst,
                              UWORD8 *pu1_u_dst, UWORD8 *pu1_v_dst, UWORD32 u4_dst_y_strd,
                              UWORD32 u4_dst_uv_strd, WORD32 cur_row, WORD32 num_rows)
{
    IH264E_ERROR_T ret = IH264E_SUCCESS;
    UWORD8 *pu1_y_src, *pu1_uv_src;
    UWORD8 *pu1_y_dst_tmp, *pu1_uv_dst_tmp;
    UWORD8 *pu1_u_dst_tmp, *pu1_v_dst_tmp;
    WORD32 is_u_first;
    UWORD8 *pu1_luma;
    UWORD8 *pu1_chroma;
    WORD32 wd;

    WORD32 src_y_strd;
    WORD32 src_uv_strd;

    WORD32 layer_id = ps_pic->u1_num_spatial_layers - 1;

    if(0 == num_rows)
    {
        return ret;
    }

    pu1_luma = ps_pic->ps_layer_yuv_buf_props[layer_id].as_component_bufs[0].pv_data;
    pu1_chroma = ps_pic->ps_layer_yuv_buf_props[layer_id].as_component_bufs[1].pv_data;

    src_y_strd = ps_pic->ps_layer_yuv_buf_props[layer_id].as_component_bufs[0].i4_data_stride;
    src_uv_strd = ps_pic->ps_layer_yuv_buf_props[layer_id].as_component_bufs[1].i4_data_stride;

    wd = ps_codec->s_cfg.u4_disp_wd;
    is_u_first = (IV_YUV_420SP_UV == ps_codec->e_codec_color_format) ? 1 : 0;

    /* In case of 420P output luma copy is disabled for shared mode */
    {
        pu1_y_src = pu1_luma + cur_row * src_y_strd;
        pu1_uv_src = pu1_chroma + (cur_row / 2) * src_uv_strd;

        pu1_y_dst_tmp = pu1_y_dst + cur_row * u4_dst_y_strd;
        pu1_uv_dst_tmp = pu1_u_dst + (cur_row / 2) * u4_dst_uv_strd;
        pu1_u_dst_tmp = pu1_u_dst + (cur_row / 2) * u4_dst_uv_strd;
        pu1_v_dst_tmp = pu1_v_dst + (cur_row / 2) * u4_dst_uv_strd;

        /* If the call is non-blocking and there are no rows to be copied then
         * return */
        /* In non-shared mode, reference buffers are in 420SP UV format,
         * if output also is in 420SP_UV, then just copy
         * if output is in 420SP_VU then swap UV values
         */
        if((IV_YUV_420SP_UV == ps_codec->s_cfg.e_recon_color_fmt) ||
           (IV_YUV_420SP_VU == ps_codec->s_cfg.e_recon_color_fmt))
        {
            ih264e_fmt_conv_420sp_to_420sp(pu1_y_src, pu1_uv_src, pu1_y_dst_tmp, pu1_uv_dst_tmp, wd,
                                           num_rows, ps_codec->i4_rec_strd, ps_codec->i4_rec_strd,
                                           u4_dst_y_strd, u4_dst_uv_strd);
        }
        else if(IV_YUV_420P == ps_codec->s_cfg.e_recon_color_fmt)
        {
            ih264e_fmt_conv_420sp_to_420p(pu1_y_src, pu1_uv_src, pu1_y_dst_tmp, pu1_u_dst_tmp,
                                          pu1_v_dst_tmp, wd, num_rows, ps_codec->i4_rec_strd,
                                          ps_codec->i4_rec_strd, u4_dst_y_strd, u4_dst_uv_strd,
                                          is_u_first, 0);
        }
    }
    return (ret);
}
