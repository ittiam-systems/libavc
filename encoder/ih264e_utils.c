/******************************************************************************
 *
 * Copyright (C) 2015 The Android Open Source Project
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
*  ih264e_utils.c
*
* @brief
*  Contains miscellaneous utility functions used by the encoder
*
* @author
*  ittiam
*
* @par List of Functions:
*  - ih264e_get_min_level()
*  - ih264e_get_lvl_idx()
*  - ih264e_get_dpb_size()
*  - ih264e_get_total_pic_buf_size()
*  - ih264e_get_pic_mv_bank_size()
*  - ih264e_pic_buf_mgr_add_bufs()
*  - ih264e_mv_buf_mgr_add_bufs()
*  - ih264e_init_quant_params()
*  - ih264e_init_air_map()
*  - ih264e_codec_init()
*  - ih264e_pic_init()
*
* @remarks
*  None
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* system include files */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* user include files */
#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "ih264e.h"
#include "ithread.h"
#include "ih264_defs.h"
#include "ih264_size_defs.h"
#include "ime_distortion_metrics.h"
#include "ime_structs.h"
#include "ih264_defs.h"
#include "ih264_error.h"
#include "ih264_structs.h"
#include "ih264_trans_quant_itrans_iquant.h"
#include "ih264_inter_pred_filters.h"
#include "ih264_mem_fns.h"
#include "ih264_padding.h"
#include "ih264_intra_pred_filters.h"
#include "ih264_deblk_edge_filters.h"
#include "ih264_macros.h"
#include "ih264_common_tables.h"
#include "ih264_debug.h"
#include "ih264_trans_data.h"
#include "ih264e_defs.h"
#include "ih264e_globals.h"
#include "ih264_buf_mgr.h"
#include "ih264_dpb_mgr.h"
#include "ih264e_error.h"
#include "ih264e_bitstream.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "ih264e_rate_control.h"
#include "ih264e_structs.h"
#include "ih264e_utils.h"
#include "ih264e_config.h"
#include "ih264e_statistics.h"
#include "ih264e_trace.h"
#include "ih264_list.h"
#include "ih264e_encode_header.h"
#include "ih264e_me.h"
#include "ime_defs.h"
#include "ime.h"
#include "ih264e_rate_control.h"
#include "ih264e_core_coding.h"
#include "ih264e_rc_mem_interface.h"
#include "ih264e_time_stamp.h"
#include "ih264e_debug.h"
#include "ih264e_process.h"
#include "ih264e_master.h"
#include "irc_rate_control_api.h"
#include "ime_statistics.h"

/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

/**
*******************************************************************************
*
* @brief
*  Used to get minimum level index for a given picture size
*
* @par Description:
*  Gets the minimum level index and then gets corresponding level.
*  Also used to ignore invalid levels like 2.3, 3.3 etc
*
* @param[in] level
*  Level of the stream
*
* @returns  Level index for a given level
*
* @remarks
*
*******************************************************************************
*/
WORD32 ih264e_get_min_level(WORD32 pic_size)
{
    WORD32 lvl_idx = MAX_LEVEL, i;

    for (i = 0; i < MAX_LEVEL; i++)
    {
        if (pic_size <= gai4_ih264_max_luma_pic_size[i])
        {
            lvl_idx = i;
            break;
        }
    }

    return gai4_ih264_levels[lvl_idx];
}

/**
*******************************************************************************
*
* @brief
*  Used to get level index for a given level
*
* @par Description:
*  Converts from level_idc (which is multiplied by 30) to an index that can be
*  used as a lookup. Also used to ignore invalid levels like 2.2 , 3.2 etc
*
* @param[in] level
*  Level of the stream
*
* @returns  Level index for a given level
*
* @remarks
*
*******************************************************************************
*/
WORD32 ih264e_get_lvl_idx(WORD32 level)
{
    WORD32 lvl_idx = 0;

    if (level < IH264_LEVEL_11)
    {
        lvl_idx = 0;
    }
    else if (level < IH264_LEVEL_12)
    {
        lvl_idx = 1;
    }
    else if (level < IH264_LEVEL_13)
    {
        lvl_idx = 2;
    }
    else if (level < IH264_LEVEL_20)
    {
        lvl_idx = 3;
    }
    else if (level < IH264_LEVEL_21)
    {
        lvl_idx = 4;
    }
    else if (level < IH264_LEVEL_22)
    {
        lvl_idx = 5;
    }
    else if (level < IH264_LEVEL_30)
    {
        lvl_idx = 6;
    }
    else if (level < IH264_LEVEL_31)
    {
        lvl_idx = 7;
    }
    else if (level < IH264_LEVEL_32)
    {
        lvl_idx = 8;
    }
    else if (level < IH264_LEVEL_40)
    {
        lvl_idx = 9;
    }
    else if (level < IH264_LEVEL_41)
    {
        lvl_idx = 10;
    }
    else if (level < IH264_LEVEL_42)
    {
        lvl_idx = 11;
    }
    else if (level < IH264_LEVEL_50)
    {
        lvl_idx = 12;
    }
    else if (level < IH264_LEVEL_51)
    {
        lvl_idx = 13;
    }
    else
    {
        lvl_idx = 14;
    }

    return (lvl_idx);
}

/**
*******************************************************************************
*
* @brief returns maximum number of pictures allowed in dpb for a given level
*
* @par Description:
*  For given width, height and level, number of pictures allowed in decoder
*  picture buffer is computed as per Annex A.3.1
*
* @param[in] level
*  level of the bit-stream
*
* @param[in] pic_size
*  width * height
*
* @returns  Number of buffers in DPB
*
* @remarks
*  From annexure A.3.1 of H264 specification,
*  max_dec_frame_buffering <= MaxDpbSize, where MaxDpbSize is equal to
*  Min( 1024 * MaxDPB / ( PicWidthInMbs * FrameHeightInMbs * 384 ), 16 ) and
*  MaxDPB is given in Table A-1 in units of 1024 bytes. However the MaxDPB size
*  presented in the look up table gas_ih264_lvl_tbl is in units of 512
*  bytes. Hence the expression is modified accordingly.
*
*******************************************************************************
*/
WORD32 ih264e_get_dpb_size(WORD32 level, WORD32 pic_size)
{
    /* dpb size */
    WORD32 max_dpb_size_bytes = 0;

    /* dec frame buffering */
    WORD32 max_dpb_size_frames = 0;

    /* temp var */
    WORD32 i;

    /* determine max luma samples */
    for (i = 0; i < 16; i++)
        if (level == (WORD32)gas_ih264_lvl_tbl[i].u4_level_idc)
            max_dpb_size_bytes = gas_ih264_lvl_tbl[i].u4_max_dpb_size;

    /* from Annexure A.3.1 h264 specification */
    max_dpb_size_frames =
                    MIN( 1024 * max_dpb_size_bytes / ( pic_size * 3 ), MAX_DPB_SIZE );

    return max_dpb_size_frames;
}

/**
*******************************************************************************
*
* @brief
*  Used to get reference picture buffer size for a given level and
*  and padding used
*
* @par Description:
*  Used to get reference picture buffer size for a given level and padding used
*  Each picture is padded on all four sides
*
* @param[in] pic_size
*  Number of luma samples (Width * Height)
*
* @param[in] level
*  Level
*
* @param[in] horz_pad
*  Total padding used in horizontal direction
*
* @param[in] vert_pad
*  Total padding used in vertical direction
*
* @returns  Total picture buffer size
*
* @remarks
*
*
*******************************************************************************
*/
WORD32 ih264e_get_total_pic_buf_size(WORD32 pic_size,
                                     WORD32 level,
                                     WORD32 horz_pad,
                                     WORD32 vert_pad,
                                     WORD32 num_ref_frames,
                                     WORD32 num_reorder_frames)
{
    WORD32 size;
    WORD32 num_luma_samples;
    WORD32 lvl_idx;
    WORD32 max_wd, min_ht;
    WORD32 num_samples;
    WORD32 max_num_bufs;
    WORD32 pad = MAX(horz_pad, vert_pad);
    UNUSED(pic_size);
    /*
     * If num_ref_frames and num_reorder_frmaes is specified
     * Use minimum value
     */
    max_num_bufs = (num_ref_frames + num_reorder_frames + MAX_CTXT_SETS);

    /* Get level index */
    lvl_idx = ih264e_get_lvl_idx(level);

    /* Maximum number of luma samples in a picture at given level */
    num_luma_samples = gai4_ih264_max_luma_pic_size[lvl_idx];

    /* Account for chroma */
    num_samples = num_luma_samples * 3 / 2;

    /* Maximum width of luma samples in a picture at given level */
    max_wd = gai4_ih264_max_wd_ht[lvl_idx];

    /* Minimum height of luma samples in a picture at given level */
    min_ht = gai4_ih264_min_wd_ht[lvl_idx];

    /* Allocation is required for
     * (Wd + horz_pad) * (Ht + vert_pad) * (2 * max_dpb_size + 1)
     *
     * Above expanded as
     * ((Wd * Ht) + (horz_pad * vert_pad) + Wd * vert_pad + Ht * horz_pad) * (2 * max_dpb_size + 1)
     * (Wd * Ht) * (2 * max_dpb_size + 1) + ((horz_pad * vert_pad) + Wd * vert_pad + Ht * horz_pad) * (2 * max_dpb_size + 1)
     * Now  max_dpb_size increases with smaller Wd and Ht, but Wd * ht * max_dpb_size will still be lesser or equal to max_wd * max_ht * dpb_size
     *
     * In the above equation (Wd * Ht) * (2 * max_dpb_size + 1) is accounted by using num_samples * (2 * max_dpb_size + 1) below
     *
     * For the padded area use MAX(horz_pad, vert_pad) as pad
     * ((pad * pad) + pad * (Wd + Ht)) * (2 * max_dpb_size + 1) has to accounted from the above for padding
     *
     * Since Width and Height can change worst Wd + Ht is when One of the dimensions is max and other is min
     * So use max_wd and min_ht
     */

    /* Number of bytes in reference pictures */
    size = num_samples * max_num_bufs;

    /* Account for padding area */
    size += ((pad * pad) + pad * (max_wd + min_ht)) * 3 / 2 * max_num_bufs;

    return size;
}

/**
*******************************************************************************
*
* @brief Returns MV bank buffer size for a given number of luma samples
*
* @par Description:
*  For given number of luma samples  one MV bank size is computed.
*  Each MV bank includes pu_map and enc_pu_t for all the min PUs(4x4) in a picture
*
* @param[in] num_luma_samples
*  Max number of luma pixels in the frame
*
* @returns  Total MV Bank size
*
* @remarks
*
*******************************************************************************
*/
WORD32 ih264e_get_pic_mv_bank_size(WORD32 num_luma_samples)
{
    /* mv bank buffer size */
    WORD32 mv_bank_size = 0;

    /* number of sub mb partitions possible */
    WORD32 num_pu = num_luma_samples / (MIN_PU_SIZE * MIN_PU_SIZE);

    /* number of mbs */
    WORD32 num_mb = num_luma_samples / (MB_SIZE * MB_SIZE);

    /* Size for storing enc_pu_t start index each MB */
    /* One extra entry is needed to compute number of PUs in the last MB */
    mv_bank_size += num_mb * sizeof(WORD32);

    /* Size for pu_map */
    mv_bank_size += num_pu;

    /* Size for storing enc_pu_t for each PU */
    mv_bank_size += num_pu * sizeof(enc_pu_t);

    return mv_bank_size;
}

/**
*******************************************************************************
*
* @brief
*  Function to initialize ps_pic_buf structs add pic buffers to
*  buffer manager in case of non-shared mode
*
* @par Description:
*  Function to initialize ps_pic_buf structs add pic buffers to
*  buffer manager in case of non-shared mode
*  To be called once per stream or for every reset
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @returns  error status
*
* @remarks
*
*******************************************************************************
*/
IH264E_ERROR_T ih264e_pic_buf_mgr_add_bufs(codec_t *ps_codec)
{
    /* error status */
    IH264E_ERROR_T ret = IH264E_SUCCESS;

    /* max ref buffer cnt */
    WORD32 max_num_bufs = ps_codec->i4_ref_buf_cnt;

    /* total size for pic buffers */
    WORD32 pic_buf_size_allocated = ps_codec->i4_total_pic_buf_size
                    - BUF_MGR_MAX_CNT * sizeof(pic_buf_t);

    /* temp var */
    UWORD8 *pu1_buf = (UWORD8 *) ps_codec->ps_pic_buf;
    pic_buf_t *ps_pic_buf = (pic_buf_t *) ps_codec->ps_pic_buf;
    WORD32 i;

    pu1_buf += BUF_MGR_MAX_CNT * sizeof(pic_buf_t);

    /* In case of non-shared mode, add picture buffers to buffer manager
     * In case of shared mode, buffers are added in the run-time
     */
    {
        WORD32 buf_ret;

        WORD32 luma_samples = (ps_codec->i4_rec_strd)
                        * (ps_codec->s_cfg.u4_ht + PAD_HT);

        WORD32 chroma_samples = luma_samples >> 1;

        /* Try and add as many buffers as possible for the memory that is allocated */
        /* If the number of buffers that can be added is less than max_num_bufs
         * return with an error */
        for (i = 0; i < max_num_bufs; i++)
        {
            pic_buf_size_allocated -= (luma_samples + chroma_samples);

            if (pic_buf_size_allocated < 0)
            {
                ps_codec->i4_error_code = IH264E_INSUFFICIENT_MEM_PICBUF;
                return IH264E_INSUFFICIENT_MEM_PICBUF;
            }

            ps_pic_buf->pu1_luma = pu1_buf + ps_codec->i4_rec_strd * PAD_TOP
                            + PAD_LEFT;
            pu1_buf += luma_samples;

            ps_pic_buf->pu1_chroma = pu1_buf
                            + ps_codec->i4_rec_strd * (PAD_TOP / 2)+ PAD_LEFT;
            pu1_buf += chroma_samples;

            buf_ret = ih264_buf_mgr_add((buf_mgr_t *) ps_codec->pv_ref_buf_mgr,
                                        ps_pic_buf, i);

            if (0 != buf_ret)
            {
                ps_codec->i4_error_code = IH264E_BUF_MGR_ERROR;
                return IH264E_BUF_MGR_ERROR;
            }
            pu1_buf += (HPEL_PLANES_CNT - 1) * (chroma_samples + luma_samples);
            ps_pic_buf++;
        }
    }

    return ret;
}

/**
*******************************************************************************
*
* @brief Function to add buffers to MV Bank buffer manager
*
* @par Description:
*  Function to add buffers to MV Bank buffer manager.  To be called once per
*  stream or for every reset
*
* @param[in] ps_codec
*  Pointer to codec context
*
* @returns  error status
*
* @remarks
*
*******************************************************************************
*/
IH264E_ERROR_T ih264e_mv_buf_mgr_add_bufs(codec_t *ps_codec)
{
    /* error status */
    IH264E_ERROR_T error_status = IH264E_SUCCESS;
    IH264_ERROR_T ret;

    /* max dpb size in frames */
    WORD32 max_dpb_size = 0;

    /* mv bank size for the entire dpb */
    WORD32 mv_bank_size_allocated = 0;

    /* mv bank size per pic */
    WORD32 pic_mv_bank_size = 0;

    /* mv buffer ptr */
    mv_buf_t *ps_mv_buf = NULL;

    /* num of luma samples */
    WORD32 num_luma_samples = ALIGN16(ps_codec->s_cfg.u4_wd)
                    * ALIGN16(ps_codec->s_cfg.u4_ht);

    /* number of mb's & frame partitions */
    WORD32 num_pu, num_mb;

    /* temp var */
    UWORD8 *pu1_buf = NULL;
    WORD32 i;

    /* Compute the number of MB Bank buffers needed */
    max_dpb_size = ps_codec->i4_ref_buf_cnt;

    /* allocate memory for mv buffer array */
    ps_codec->ps_mv_buf = ps_codec->pv_mv_bank_buf_base;
    pu1_buf = ps_codec->pv_mv_bank_buf_base;
    pu1_buf += BUF_MGR_MAX_CNT * sizeof(mv_buf_t);

    /********************************************************************/
    /* allocate memory for individual elements of mv buffer ptr         */
    /********************************************************************/
    mv_bank_size_allocated = ps_codec->i4_total_mv_bank_size
                    - (BUF_MGR_MAX_CNT * sizeof(mv_buf_t));

    /* compute MV bank size per picture */
    pic_mv_bank_size = ih264e_get_pic_mv_bank_size(num_luma_samples);

    num_pu = num_luma_samples / (MIN_PU_SIZE * MIN_PU_SIZE);
    num_mb = num_luma_samples / (MB_SIZE * MB_SIZE);
    i = 0;
    ps_mv_buf = ps_codec->pv_mv_bank_buf_base;

    while (i < max_dpb_size)
    {
        mv_bank_size_allocated -= pic_mv_bank_size;

        if (mv_bank_size_allocated < 0)
        {
            ps_codec->i4_error_code = IH264E_INSUFFICIENT_MEM_MVBANK;

            error_status = IH264E_INSUFFICIENT_MEM_MVBANK;

            return error_status;
        }

        ps_mv_buf->pu4_mb_pu_cnt = (UWORD32 *) pu1_buf;

        ps_mv_buf->pu1_pic_pu_map = (pu1_buf + num_mb * sizeof(WORD32));

        ps_mv_buf->ps_pic_pu = (enc_pu_t *) (pu1_buf + num_mb * sizeof(WORD32)
                        + num_pu);

        ret = ih264_buf_mgr_add((buf_mgr_t *) ps_codec->pv_mv_buf_mgr,
                                ps_mv_buf, i);

        if (IH264_SUCCESS != ret)
        {
            ps_codec->i4_error_code = IH264E_BUF_MGR_ERROR;
            error_status = IH264E_BUF_MGR_ERROR;
            return error_status;
        }

        pu1_buf += pic_mv_bank_size;
        ps_mv_buf++;
        i++;
    }

    return error_status;
}

/**
*******************************************************************************
*
* @brief Function to initialize quant params structure
*
* @par Description:
*  The forward quantization modules depends on qp/6, qp mod 6, forward scale
*  matrix, forward threshold matrix, weight list. The inverse quantization
*  modules depends on qp/6, qp mod 6, inverse scale matrix, weight list.
*  These params are initialized in this function.
*
* @param[in] ps_proc
*  pointer to process context
*
* @param[in] qp
*  quantization parameter
*
* @returns none
*
* @remarks
*
*******************************************************************************
*/
void ih264e_init_quant_params(process_ctxt_t *ps_proc, int qp)
{
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
    while (plane <= V)
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
        for (i = 0; i < 16; i++)
        {
            ps_qp_params->pu2_thres_mat[i] = pu2_thres_mat[i]
                            >> (8 - u4_qp_div6);
            ps_qp_params->pu2_weigh_mat[i] = 16;
        }

        /* qp dependent rounding constant */
        ps_qp_params->u4_dead_zone =
                        gu4_forward_quant_round_factor_4x4[u4_qp_div6];

        /* slice dependent rounding constant */
        if (ps_proc->i4_slice_type != ISLICE
                        && ps_proc->i4_slice_type != SISLICE)
        {
            ps_qp_params->u4_dead_zone >>= 1;
        }

        /* SATQD threshold for zero block prediction */
        if (ps_proc->ps_codec->s_cfg.u4_enable_satqd)
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
            ps_qp_params->pu2_sad_thrsh[8] = u4_satdq_t / MAX(MAX3(pu2_smat[0], pu2_smat[2], pu2_smat[8]), pu2_smat[10]);
        }

        /* inverse scale matrix */
        ps_qp_params->pu2_iscale_mat = pu2_iscale_mat + (u4_qp_mod6 * 16);

        plane += 1;
    }
    return ;
}

/**
*******************************************************************************
*
* @brief
*  Initialize AIR mb frame Map
*
* @par Description:
*  Initialize AIR mb frame map
*  MB frame map indicates which frame an Mb should be coded as intra according to AIR
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
IH264E_ERROR_T ih264e_init_air_map(codec_t *ps_codec)
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

    switch (air_mode)
    {
        case IVE_AIR_MODE_CYCLIC:

            for (curr_mb = 0; curr_mb < u4_mb_cnt; curr_mb++)
            {
                pu2_intr_rfrsh_map[curr_mb] = curr_mb % air_period;
            }
            break;

        case IVE_AIR_MODE_RANDOM:

            for (curr_mb = 0; curr_mb < u4_mb_cnt; curr_mb++)
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
IH264E_ERROR_T ih264e_codec_init(codec_t *ps_codec)
{
    /********************************************************************
     *                     INITIALIZE CODEC CONTEXT                     *
     ********************************************************************/
    /* encoder presets */
    if (ps_codec->s_cfg.u4_enc_speed_preset != IVE_CONFIG)
    {
        if (ps_codec->s_cfg.u4_enc_speed_preset == IVE_SLOWEST)
        {/* high quality */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;
            ps_codec->s_cfg.u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 1;
            ps_codec->luma_energy_compaction[1] =
                            ih264e_code_luma_intra_macroblock_4x4_rdopt_on;

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 1;

            /* deblocking off */
            ps_codec->s_cfg.u4_disable_deblock_level = DISABLE_DEBLK_LEVEL_0;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if (ps_codec->s_cfg.u4_enc_speed_preset == IVE_NORMAL)
        {/* normal */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;
            ps_codec->s_cfg.u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 1;

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 1;

            /* deblocking off */
            ps_codec->s_cfg.u4_disable_deblock_level = DISABLE_DEBLK_LEVEL_0;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if (ps_codec->s_cfg.u4_enc_speed_preset == IVE_FAST)
         {/* normal */
             /* enable diamond search */
             ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;
             ps_codec->s_cfg.u4_enable_fast_sad = 0;

             /* disable intra 4x4 */
             ps_codec->s_cfg.u4_enable_intra_4x4 = 0;

             /* sub pel off */
             ps_codec->s_cfg.u4_enable_hpel = 1;

             /* deblocking off */
             ps_codec->s_cfg.u4_disable_deblock_level = DISABLE_DEBLK_LEVEL_0;

             /* disabled intra inter gating in Inter slices */
             ps_codec->u4_inter_gate = 1;
         }
        else if (ps_codec->s_cfg.u4_enc_speed_preset == IVE_HIGH_SPEED)
        {/* fast */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;
            ps_codec->s_cfg.u4_enable_fast_sad = 0;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 0;

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 0;

            /* deblocking off */
            ps_codec->s_cfg.u4_disable_deblock_level = DISABLE_DEBLK_LEVEL_4;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 0;
        }
        else if (ps_codec->s_cfg.u4_enc_speed_preset == IVE_FASTEST)
        {/* fastest */
            /* enable diamond search */
            ps_codec->s_cfg.u4_me_speed_preset = DMND_SRCH;

            /* disable intra 4x4 */
            ps_codec->s_cfg.u4_enable_intra_4x4 = 0;

            /* sub pel off */
            ps_codec->s_cfg.u4_enable_hpel = 0;

            /* deblocking off */
            ps_codec->s_cfg.u4_disable_deblock_level = DISABLE_DEBLK_LEVEL_4;

            /* disabled intra inter gating in Inter slices */
            ps_codec->u4_inter_gate = 1;
        }
    }

    /*****************************************************************
     * Initialize AIR inside codec
     *****************************************************************/
    if (IVE_AIR_MODE_NONE != ps_codec->s_cfg.e_air_mode)
    {
        ih264e_init_air_map(ps_codec);

        ps_codec->i4_air_pic_cnt = -1;
    }

    /****************************************************/
    /*           INITIALIZE RATE CONTROL                */
    /****************************************************/
    {
        /* init qp */
        UWORD8 au1_init_qp[MAX_PIC_TYPE];

        /* min max qp */
        UWORD8 au1_min_max_qp[2 * MAX_PIC_TYPE];

        /* init i,p,b qp */
        au1_init_qp[0] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_i_qp];
        au1_init_qp[1] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_p_qp];
        au1_init_qp[2] = gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_b_qp];

        /* init min max qp */
        au1_min_max_qp[2 * I_PIC] =
                        gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_i_qp_min];
        au1_min_max_qp[2 * I_PIC + 1] =
                        gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_i_qp_max];

        au1_min_max_qp[2 * P_PIC] =
                        gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_p_qp_min];
        au1_min_max_qp[2 * P_PIC + 1] =
                        gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_p_qp_max];

        au1_min_max_qp[2 * B_PIC] =
                        gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_b_qp_min];
        au1_min_max_qp[2 * B_PIC + 1] =
                        gau1_h264_to_mpeg2_qmap[ps_codec->s_cfg.u4_b_qp_max];

        /* get rc mode */
        switch (ps_codec->s_cfg.e_rc_mode)
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
        ih264e_rc_init(ps_codec->s_rate_control.pps_rate_control_api,
                       ps_codec->s_rate_control.pps_frame_time,
                       ps_codec->s_rate_control.pps_time_stamp,
                       ps_codec->s_rate_control.pps_pd_frm_rate,
                       ps_codec->s_cfg.u4_max_framerate,
                       ps_codec->s_cfg.u4_src_frame_rate,
                       ps_codec->s_cfg.u4_tgt_frame_rate,
                       ps_codec->s_rate_control.e_rc_type,
                       ps_codec->s_cfg.u4_target_bitrate,
                       ps_codec->s_cfg.u4_max_bitrate,
                       ps_codec->s_cfg.u4_vbv_buffer_delay,
                       ps_codec->s_cfg.u4_i_frm_interval, au1_init_qp,
                       H264_ALLOC_INTER_FRM_INTV, au1_min_max_qp,
                       ps_codec->s_cfg.u4_max_level);
    }

    /* src stride */
    ps_codec->i4_src_strd = ps_codec->s_cfg.u4_strd;

    /* recon stride */
    ps_codec->i4_rec_strd = ALIGN16(ps_codec->s_cfg.u4_max_wd) + PAD_WD;

    /* max ref and reorder cnt */
    ps_codec->i4_ref_buf_cnt = ps_codec->s_cfg.u4_max_ref_cnt
                    + ps_codec->s_cfg.u4_max_reorder_cnt;
    ps_codec->i4_ref_buf_cnt += MAX_CTXT_SETS;

    DEBUG_HISTOGRAM_INIT();

    return IH264E_SUCCESS;
}

/**
*******************************************************************************
*
* @brief
*  Picture level initializations
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
IH264E_ERROR_T ih264e_pic_init(codec_t *ps_codec, inp_buf_t *ps_inp_buf)
{
    /* error status */
    IH264E_ERROR_T error_status = IH264E_SUCCESS;
    IH264_ERROR_T ret = IH264_SUCCESS;

    /* mv buff bank */
    mv_buf_t *ps_mv_buf = NULL;
    WORD32 cur_mv_bank_buf_id;

    /* recon buffer set */
    pic_buf_t *ps_cur_pic;
    WORD32 cur_pic_buf_id;
    UWORD8 *pu1_cur_pic_luma, *pu1_cur_pic_chroma;

    /* ref buffer set */
    pic_buf_t *ps_ref_pic;
    WORD32 ref_set_id;

    /* pic time stamp */
    UWORD32 u4_timestamp_high = ps_inp_buf->u4_timestamp_high;
    UWORD32 u4_timestamp_low = ps_inp_buf->u4_timestamp_low;

    /* indices to access curr/prev frame info */
    WORD32 ctxt_sel = ps_codec->i4_encode_api_call_cnt & 1;

    /* curr pic type */
    PIC_TYPE_T *pic_type = &ps_codec->pic_type;

    /* should src be skipped */
    WORD32 *skip_src = &ps_codec->s_rate_control.pre_encode_skip[ctxt_sel];

    /* Diamond search Iteration Max Cnt */
    UWORD32 u4_num_layers =
                    (ps_codec->s_cfg.u4_enc_speed_preset == IVE_FASTEST) ?
                                    (NUM_LAYERS >> 2) : NUM_LAYERS;

    /* enable fast sad */
    UWORD32 u4_enable_fast_sad = ps_codec->s_cfg.u4_enable_fast_sad;

    /********************************************************************/
    /*                     INITIALIZE CODEC CONTEXT                     */
    /********************************************************************/

    /* pre enc rc call */
    *skip_src = ih264e_set_rc_pic_params(ps_codec,
                                         ps_codec->i4_encode_api_call_cnt,
                                         (WORD32 *) pic_type);
    if (*skip_src == 1)
    {
        ps_codec->as_process[ctxt_sel * MAX_PROCESS_THREADS].s_inp_buf =
                        *ps_inp_buf;

        /* inform output bytes generated as zero */
        ps_codec->as_out_buf[ctxt_sel].s_bits_buf.u4_bytes = 0;

        return error_status;
    }

    /********************************************************************/
    /*                     Alternate reference frame                    */
    /********************************************************************/
    if (ps_codec->s_cfg.u4_enable_alt_ref)
    {
        if (PIC_IDR == *pic_type || PIC_I == *pic_type)
        {
            ps_codec->u4_is_curr_frm_ref = 1;
        }
        else
        {
            ps_codec->u4_is_curr_frm_ref = 1;
                if(ps_codec->i4_encode_api_call_cnt % (ps_codec->s_cfg.u4_enable_alt_ref + 1))
                    ps_codec->u4_is_curr_frm_ref = 0;
            }

        if ((ps_codec->u4_is_curr_frm_ref == 1) || (ps_codec->i4_frame_num < 0))
        {
            ps_codec->i4_frame_num++;
        }
    }
    else
    {
        ps_codec->u4_is_curr_frm_ref = 1;

        ps_codec->i4_frame_num++;
    }

    /* slice_type */
    ps_codec->i4_slice_type = PSLICE;

    if ((PIC_I == *pic_type) || (PIC_IDR == *pic_type))
    {
        ps_codec->i4_slice_type = ISLICE;
    }
    else if (PIC_P == *pic_type)
    {
        ps_codec->i4_slice_type = PSLICE;
    }

    /* is this an IDR pic */
    ps_codec->u4_is_idr = 0;

    if (PIC_IDR == *pic_type)
    {
        /* set idr flag */
        ps_codec->u4_is_idr = 1;

        /* reset frame num */
        ps_codec->i4_frame_num = 0;

        /* idr_pic_id */
        ps_codec->i4_idr_pic_id++;
    }

    /* set deblock disable flags based on disable deblock level */
    ps_codec->i4_disable_deblk_pic = 1;

    if (ps_codec->s_cfg.u4_disable_deblock_level == DISABLE_DEBLK_LEVEL_0)
    {
        /* enable deblocking */
        ps_codec->i4_disable_deblk_pic = 0;
    }
    else if (ps_codec->s_cfg.u4_disable_deblock_level == DISABLE_DEBLK_LEVEL_2)
    {
        /* enable deblocking after a period of frames */
        if (ps_codec->i4_disable_deblk_pic_cnt == DISABLE_DEBLOCK_INTERVAL
                        || ps_codec->i4_slice_type == ISLICE)
        {
            ps_codec->i4_disable_deblk_pic = 0;
        }
    }
    else if (ps_codec->s_cfg.u4_disable_deblock_level == DISABLE_DEBLK_LEVEL_3)
    {
        if (ps_codec->i4_slice_type == ISLICE)
        {
            ps_codec->i4_disable_deblk_pic = 0;
        }
    }

    if (ps_codec->i4_disable_deblk_pic)
    {
        ps_codec->i4_disable_deblk_pic_cnt++;
    }
    else
    {
        ps_codec->i4_disable_deblk_pic_cnt = 0;
    }

    /* In slice mode - lets not deblk mb edges that lie along slice boundaries */
    if (ps_codec->i4_disable_deblk_pic == 0)
    {
        if (ps_codec->s_cfg.e_slice_mode != IVE_SLICE_MODE_NONE)
        {
            ps_codec->i4_disable_deblk_pic = 2;
        }
    }

    /* error status */
    ps_codec->i4_error_code = IH264E_SUCCESS;

    /* populate header */
    if (ps_codec->i4_gen_header)
    {
        /* sps */
        sps_t *ps_sps = NULL;

        /* pps */
        pps_t *ps_pps = NULL;

        /*ps_codec->i4_pps_id ++;*/
        ps_codec->i4_pps_id %= MAX_PPS_CNT;

        /*ps_codec->i4_sps_id ++;*/
        ps_codec->i4_sps_id %= MAX_SPS_CNT;

        /* populate sps header */
        ps_sps = ps_codec->ps_sps_base + ps_codec->i4_sps_id;
        ih264e_populate_sps(ps_codec, ps_sps);

        /* populate pps header */
        ps_pps = ps_codec->ps_pps_base + ps_codec->i4_pps_id;
        ih264e_populate_pps(ps_codec, ps_pps);
    }

    /* Reference and MV bank Buffer Manager */
    {
        /* min pic cnt among the list of pics stored in ref list */
        WORD32 min_pic_cnt;

        /* max pic cnt among the list of pics stored in ref list */
        WORD32 max_pic_cnt;

        /* temp var */
        WORD32 i;

        ps_ref_pic = NULL;

        /* get reference picture when necessary */
        /* Only nearest picture encoded (max pic cnt) is used as reference */
        if ((*pic_type != PIC_IDR) && (*pic_type != PIC_I))
        {
            max_pic_cnt = ps_codec->as_ref_set[0].i4_pic_cnt;

            ps_ref_pic = ps_codec->as_ref_set[0].ps_pic_buf;

            /* loop through to get the max pic cnt among the list of pics stored in ref list */
            for (i = 1; i < ps_codec->i4_ref_buf_cnt; i++)
            {
                if (max_pic_cnt < ps_codec->as_ref_set[i].i4_pic_cnt)
                {
                    max_pic_cnt = ps_codec->as_ref_set[i].i4_pic_cnt;
                    ps_ref_pic = ps_codec->as_ref_set[i].ps_pic_buf;
                }
            }
        }

        /* get a location at which the curr pic info can be stored for future reference */
        ref_set_id = -1;

        for (i = 0; i < ps_codec->i4_ref_buf_cnt; i++)
        {
            if (-1 == ps_codec->as_ref_set[i].i4_pic_cnt)
            {
                ref_set_id = i;
                break;
            }
        }

        /* If all the entries in the ref_set array are filled, then remove the entry with least pic_cnt */
        if (ref_set_id == -1)
        {
            /* pic info */
            pic_buf_t *ps_cur_pic;

            /* mv info */
            mv_buf_t *ps_cur_mv_buf;

            ref_set_id = 0;
            min_pic_cnt = ps_codec->as_ref_set[0].i4_pic_cnt;

            /* loop through to get the min pic cnt among the list of pics stored in ref list */
            for (i = 1; i < ps_codec->i4_ref_buf_cnt; i++)
            {
                if (min_pic_cnt > ps_codec->as_ref_set[i].i4_pic_cnt)
                {
                    min_pic_cnt = ps_codec->as_ref_set[i].i4_pic_cnt;
                    ref_set_id = i;
                }
            }

            ps_cur_pic = ps_codec->as_ref_set[ref_set_id].ps_pic_buf;

            ps_cur_mv_buf = ps_codec->as_ref_set[ref_set_id].ps_mv_buf;

            /* release this frame from reference list */
            ih264_buf_mgr_release(ps_codec->pv_mv_buf_mgr,
                                  ps_cur_mv_buf->i4_buf_id, BUF_MGR_REF);

            ih264_buf_mgr_release(ps_codec->pv_ref_buf_mgr,
                                  ps_cur_pic->i4_buf_id, BUF_MGR_REF);
        }

        if (ps_codec->s_cfg.u4_enable_recon)
        {
            ret = ih264_buf_mgr_check_free((buf_mgr_t *)ps_codec->pv_ref_buf_mgr);

            if (ret != IH264_SUCCESS)
            {
                return IH264E_NO_FREE_RECONBUF;
            }
        }
    }

    {
        /*****************************************************************/
        /* Get free MV Bank to hold current picture's motion vector data */
        /* If there are no free buffers then return with an error code.  */
        /* If the buffer is to be freed by another thread, change the    */
        /* following to call thread yield and wait for buffer to be freed*/
        /*****************************************************************/
        ps_mv_buf = (mv_buf_t *) ih264_buf_mgr_get_next_free(
                        (buf_mgr_t *) ps_codec->pv_mv_buf_mgr,
                        &cur_mv_bank_buf_id);

        if (NULL == ps_mv_buf)
        {
            ps_codec->i4_error_code = IH264E_NO_FREE_MVBANK;
            return IH264E_NO_FREE_MVBANK;
        }

        /* mark the buffer as needed for reference if the curr pic is available for ref */
        if (ps_codec->u4_is_curr_frm_ref)
        {
            ih264_buf_mgr_set_status(ps_codec->pv_mv_buf_mgr,
                                     cur_mv_bank_buf_id, BUF_MGR_REF);
        }

        /* Set current ABS poc to ps_mv_buf, so that while freeing a reference buffer
         * corresponding mv buffer can be found by looping through ps_codec->ps_mv_buf array
         * and getting a buffer id to free
         */
        ps_mv_buf->i4_abs_poc = ps_codec->i4_abs_pic_order_cnt;

        ps_mv_buf->i4_buf_id = cur_mv_bank_buf_id;
    }

    {
        /*****************************************************************/
        /* Get free pic buf to hold current picture's recon data         */
        /* If there are no free buffers then return with an error code.  */
        /* If the buffer is to be freed by another thread, change the    */
        /* following to call thread yield and wait for buffer to be freed*/
        /*****************************************************************/
        ps_cur_pic = (pic_buf_t *) ih264_buf_mgr_get_next_free(
                        (buf_mgr_t *) ps_codec->pv_ref_buf_mgr,
                        &cur_pic_buf_id);

        if (NULL == ps_cur_pic)
        {
            ps_codec->i4_error_code = IH264E_NO_FREE_PICBUF;
            return IH264E_NO_FREE_PICBUF;
        }

        /* mark the buffer as needed for reference if the curr pic is available for ref */
        if (1 == ps_codec->u4_is_curr_frm_ref)
        {
            ih264_buf_mgr_set_status(ps_codec->pv_ref_buf_mgr, cur_pic_buf_id,
                                     BUF_MGR_REF);
        }

        /* Mark the current buffer as needed for IO if recon is enabled */
        if (1 == ps_codec->s_cfg.u4_enable_recon)
        {
            ih264_buf_mgr_set_status(ps_codec->pv_ref_buf_mgr, cur_pic_buf_id,
                                     BUF_MGR_IO);
        }

        /* Associate input timestamp with current buffer */
        ps_cur_pic->u4_timestamp_high = ps_inp_buf->u4_timestamp_high;
        ps_cur_pic->u4_timestamp_low = ps_inp_buf->u4_timestamp_low;

        ps_cur_pic->i4_abs_poc = ps_codec->i4_abs_pic_order_cnt;
        ps_cur_pic->i4_poc_lsb = ps_codec->i4_pic_order_cnt_lsb;

        ps_cur_pic->i4_buf_id = cur_pic_buf_id;

        pu1_cur_pic_luma = ps_cur_pic->pu1_luma;
        pu1_cur_pic_chroma = ps_cur_pic->pu1_chroma;
    }

    /* in case the current picture is used for reference then add it to the reference set */
    if (ps_codec->u4_is_curr_frm_ref
                    && ((*pic_type == PIC_IDR) || (*pic_type == PIC_I)
                                    || (*pic_type == PIC_P)))
    {
        ps_codec->as_ref_set[ref_set_id].i4_pic_cnt = ps_codec->i4_pic_cnt;

        /* TODO: Currently pic_cnt and poc are same - Once frame drops are introduced change appropriately */
        ps_codec->as_ref_set[ref_set_id].i4_poc = ps_codec->i4_pic_cnt;

        ps_codec->as_ref_set[ref_set_id].ps_mv_buf = ps_mv_buf;

        ps_codec->as_ref_set[ref_set_id].ps_pic_buf = ps_cur_pic;
    }

    /********************************************************************/
    /*                     INITIALIZE PROCESS CONTEXT                   */
    /********************************************************************/
    {
        /* temp var */
        WORD32 i, j = 0;

        /* curr proc ctxt */
        process_ctxt_t *ps_proc = NULL;

        j = ctxt_sel * MAX_PROCESS_THREADS;

        /* begin init */
        for (i = j; i < (j + MAX_PROCESS_THREADS); i++)
        {
            ps_proc = &ps_codec->as_process[i];

            /* luma src buffer */
            if (ps_codec->s_cfg.e_inp_color_fmt == IV_YUV_422ILE)
            {
                ps_proc->pu1_src_buf_luma_base = ps_codec->pu1_y_csc_buf_base;
            }
            else
            {
                ps_proc->pu1_src_buf_luma_base =
                                ps_inp_buf->s_raw_buf.apv_bufs[0];
            }

            /* chroma src buffer */
            if (ps_codec->s_cfg.e_inp_color_fmt == IV_YUV_422ILE
                            || ps_codec->s_cfg.e_inp_color_fmt == IV_YUV_420P)
            {
                ps_proc->pu1_src_buf_chroma_base =
                                ps_codec->pu1_uv_csc_buf_base;
            }
            else
            {
                ps_proc->pu1_src_buf_chroma_base =
                                ps_inp_buf->s_raw_buf.apv_bufs[1];
            }

            /* luma rec buffer */
            ps_proc->pu1_rec_buf_luma_base = pu1_cur_pic_luma;

            /* chroma rec buffer */
            ps_proc->pu1_rec_buf_chroma_base = pu1_cur_pic_chroma;

            /* src stride */
            ps_proc->i4_src_strd = ps_codec->i4_src_strd;

            /* rec stride */
            ps_proc->i4_rec_strd = ps_codec->i4_rec_strd;

            /* frame num */
            ps_proc->i4_frame_num = ps_codec->i4_frame_num;

            /* is idr */
            ps_proc->u4_is_idr = ps_codec->u4_is_idr;

            /* idr pic id */
            ps_proc->u4_idr_pic_id = ps_codec->i4_idr_pic_id;

            /* slice_type */
            ps_proc->i4_slice_type = ps_codec->i4_slice_type;

            /* Input width in mbs */
            ps_proc->i4_wd_mbs = ps_codec->s_cfg.i4_wd_mbs;

            /* Input height in mbs */
            ps_proc->i4_ht_mbs = ps_codec->s_cfg.i4_ht_mbs;

            /* Half x plane offset from pic buf */
            ps_proc->u4_half_x_offset = 0;

            /* Half y plane offset from half x plane */
            ps_proc->u4_half_y_offset = 0;

            /* Half x plane offset from half y plane */
            ps_proc->u4_half_xy_offset = 0;

            /* top row syntax elements */
            ps_proc->ps_top_row_mb_syntax_ele =
                            ps_proc->ps_top_row_mb_syntax_ele_base;

            ps_proc->pu1_top_mb_intra_modes =
                            ps_proc->pu1_top_mb_intra_modes_base;

            ps_proc->ps_top_row_pu = ps_proc->ps_top_row_pu_base;

            /* initialize quant params */
            ps_proc->u4_frame_qp = ps_codec->u4_frame_qp;
            ps_proc->u4_mb_qp = ps_codec->u4_frame_qp;
            ih264e_init_quant_params(ps_proc, ps_proc->u4_frame_qp);

            /* previous mb qp*/
            ps_proc->u4_mb_qp_prev = ps_proc->u4_frame_qp;

            /* Reset frame info */
            memset(&ps_proc->s_frame_info, 0, sizeof(frame_info_t));

            /* initialize proc, deblk and ME map */
            if (i == j)
            {
                /* row '-1' */
                memset(ps_proc->pu1_proc_map - ps_proc->i4_wd_mbs, 1, ps_proc->i4_wd_mbs);
                /* row 0 to ht in mbs */
                memset(ps_proc->pu1_proc_map, 0, ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs);

                /* row '-1' */
                memset(ps_proc->pu1_deblk_map - ps_proc->i4_wd_mbs, 1, ps_proc->i4_wd_mbs);
                /* row 0 to ht in mbs */
                memset(ps_proc->pu1_deblk_map, 0, ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs);

                /* row '-1' */
                memset(ps_proc->pu1_me_map - ps_proc->i4_wd_mbs, 1, ps_proc->i4_wd_mbs);
                /* row 0 to ht in mbs */
                memset(ps_proc->pu1_me_map, 0, ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs);

                /* at the start of air refresh period, reset intra coded map */
                if (IVE_AIR_MODE_NONE != ps_codec->s_cfg.e_air_mode)
                {
                    ps_codec->i4_air_pic_cnt = (ps_codec->i4_air_pic_cnt + 1)
                                    % ps_codec->s_cfg.u4_air_refresh_period;

                    if (!ps_codec->i4_air_pic_cnt)
                    {
                        memset(ps_proc->pu1_is_intra_coded, 0, ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs);
                    }
                }
            }

            /* deblock level */
            ps_proc->u4_disable_deblock_level = ps_codec->i4_disable_deblk_pic;

            /* slice index map */
            /* no slice */
            if (ps_codec->s_cfg.e_slice_mode == IVE_SLICE_MODE_NONE)
            {
                memset(ps_proc->pu1_slice_idx, 0, ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs);
            }
            /* generate slices for every 'n' rows, 'n' is given through slice param */
            else if (ps_codec->s_cfg.e_slice_mode == IVE_SLICE_MODE_BLOCKS)
            {
                /* slice idx map */
                UWORD8 *pu1_slice_idx = ps_proc->pu1_slice_idx;

                /* temp var */
                WORD32 i4_mb_y = 0, slice_idx = 0, cnt;

                while (i4_mb_y < ps_proc->i4_ht_mbs)
                {
                    if (i4_mb_y +(WORD32)ps_codec->s_cfg.u4_slice_param < ps_proc->i4_ht_mbs)
                    {
                        cnt = ps_codec->s_cfg.u4_slice_param * ps_proc->i4_wd_mbs;
                        i4_mb_y += ps_codec->s_cfg.u4_slice_param;
                    }
                    else
                    {
                        cnt = (ps_proc->i4_ht_mbs - i4_mb_y) * ps_proc->i4_wd_mbs;
                        i4_mb_y += (ps_proc->i4_ht_mbs - i4_mb_y);
                    }
                    memset(pu1_slice_idx, slice_idx, cnt);
                    slice_idx++;
                    pu1_slice_idx += cnt;
                }
            }

            /* Current MV Bank's buffer ID */
            ps_proc->i4_cur_mv_bank_buf_id = cur_mv_bank_buf_id;

            /* Pointer to current picture buffer structure */
            ps_proc->ps_cur_pic = ps_cur_pic;

            /* Pointer to current pictures mv buffers */
            ps_proc->ps_cur_mv_buf = ps_mv_buf;

            /* pointer to ref picture */
            ps_proc->ps_ref_pic = ps_ref_pic;

            if ((*pic_type != PIC_IDR) && (*pic_type != PIC_I))
            {
                /* ref pointer luma */
                ps_proc->pu1_ref_buf_luma_base = ps_ref_pic->pu1_luma;

                /* ref pointer chroma */
                ps_proc->pu1_ref_buf_chroma_base = ps_ref_pic->pu1_chroma;
            }

            /* Structure for current input buffer */
            ps_proc->s_inp_buf = *ps_inp_buf;

            /* Number of encode frame API calls made */
            ps_proc->i4_encode_api_call_cnt = ps_codec->i4_encode_api_call_cnt;

            /* Current Picture count */
            ps_proc->i4_pic_cnt = ps_codec->i4_pic_cnt;

            /* error status */
            ps_proc->i4_error_code = 0;

            /********************************************************************/
            /*                     INITIALIZE ENTROPY CONTEXT                   */
            /********************************************************************/
            {
                entropy_ctxt_t *ps_entropy = &ps_proc->s_entropy;

                /* start of frame */
                ps_entropy->i4_sof = 0;

                /* end of frame */
                ps_entropy->i4_eof = 0;

                /* generate header */
                ps_entropy->i4_gen_header = ps_codec->i4_gen_header;

                /* sps ref_set_id */
                ps_entropy->u4_sps_id = ps_codec->i4_sps_id;

                /* sps base */
                ps_entropy->ps_sps_base = ps_codec->ps_sps_base;

                /* sps id */
                ps_entropy->u4_pps_id = ps_codec->i4_pps_id;

                /* sps base */
                ps_entropy->ps_pps_base = ps_codec->ps_pps_base;

                /* slice map */
                ps_entropy->pu1_slice_idx = ps_proc->pu1_slice_idx;

                /* slice hdr base */
                ps_entropy->ps_slice_hdr_base = ps_proc->ps_slice_hdr_base;

                /* initialize entropy map */
                if (i == j)
                {
                    /* row '-1' */
                    memset(ps_entropy->pu1_entropy_map - ps_proc->i4_wd_mbs, 1, ps_proc->i4_wd_mbs);
                    /* row 0 to ht in mbs */
                    memset(ps_entropy->pu1_entropy_map, 0, ps_proc->i4_wd_mbs * ps_proc->i4_ht_mbs);
                }

                /* wd in mbs */
                ps_entropy->i4_wd_mbs = ps_proc->i4_wd_mbs;

                /* ht in mbs */
                ps_entropy->i4_ht_mbs = ps_proc->i4_ht_mbs;

                /* transform_8x8_mode_flag */
                ps_entropy->i1_transform_8x8_mode_flag = 0;

                /* entropy_coding_mode_flag */
                ps_entropy->u1_entropy_coding_mode_flag =
                                ps_codec->s_cfg.u4_entropy_coding_mode;

                /* error code */
                ps_entropy->i4_error_code = IH264E_SUCCESS;

                /* mb skip run */
                *(ps_proc->s_entropy.pi4_mb_skip_run) = 0;

                /* last frame to encode */
                ps_proc->s_entropy.u4_is_last = ps_inp_buf->u4_is_last;

                /* Current Picture count */
                ps_proc->s_entropy.i4_pic_cnt = ps_codec->i4_pic_cnt;

                /* time stamps */
                ps_entropy->u4_timestamp_low = u4_timestamp_low;
                ps_entropy->u4_timestamp_high = u4_timestamp_high;

                /* init frame statistics */
                ps_entropy->u4_header_bits[MB_TYPE_INTRA] = 0;
                ps_entropy->u4_header_bits[MB_TYPE_INTER] = 0;
                ps_entropy->u4_residue_bits[MB_TYPE_INTRA] = 0;
                ps_entropy->u4_residue_bits[MB_TYPE_INTER] = 0;
            }

            /********************************************************************/
            /*                     INITIALIZE DEBLOCK CONTEXT                   */
            /********************************************************************/
            {
                /* deblk ctxt */
                deblk_ctxt_t *ps_deblk = &ps_proc->s_deblk_ctxt;

                /* slice idx map */
                ps_deblk->pu1_slice_idx = ps_proc->pu1_slice_idx;
            }

            /********************************************************************/
            /*                     INITIALIZE ME CONTEXT                        */
            /********************************************************************/
            {
                /* me ctxt */
                me_ctxt_t *ps_me_ctxt = &ps_proc->s_me_ctxt;

                /* srch range x */
                ps_me_ctxt->ai2_srch_boundaries[0] =
                                ps_codec->s_cfg.u4_srch_rng_x;

                /* srch range y */
                ps_me_ctxt->ai2_srch_boundaries[1] =
                                ps_codec->s_cfg.u4_srch_rng_y;

                /* src stride */
                ps_me_ctxt->i4_src_strd = ps_codec->i4_src_strd;

                /* rec stride */
                ps_me_ctxt->i4_rec_strd = ps_codec->i4_rec_strd;

                /* Half x plane offset from pic buf */
                ps_me_ctxt->u4_half_x_offset = ps_proc->u4_half_x_offset;

                /* Half y plane offset from half x plane */
                ps_me_ctxt->u4_half_y_offset = ps_proc->u4_half_y_offset;

                /* Half x plane offset from half y plane */
                ps_me_ctxt->u4_half_xy_offset = ps_proc->u4_half_xy_offset;

                /* enable fast sad */
                ps_me_ctxt->u4_enable_fast_sad = u4_enable_fast_sad;

                /* half pel */
                ps_me_ctxt->u4_enable_hpel = ps_codec->s_cfg.u4_enable_hpel;

                /* Diamond search Iteration Max Cnt */
                ps_me_ctxt->u4_num_layers = u4_num_layers;

                /* me speed preset */
                ps_me_ctxt->u4_me_speed_preset =
                                ps_codec->s_cfg.u4_me_speed_preset;

                /* qp */
                ps_me_ctxt->u1_mb_qp = ps_codec->u4_frame_qp;

                if ((i == 0) && (0 == ps_codec->i4_pic_cnt))
                {
                    /* init mv bits tables */
                    ih264e_init_mv_bits(ps_me_ctxt);
                }
            }

            ps_proc->ps_ngbr_avbl = &(ps_proc->s_ngbr_avbl);

        }

        /* reset encoder header */
        ps_codec->i4_gen_header = 0;
    }

    /********************************************************************/
    /*                       ADD JOBS TO THE QUEUE                      */
    /********************************************************************/
    {
        /* job structures */
        job_t s_job;

        /* temp var */
        WORD32 i;

        /* job class */
        s_job.i4_cmd = CMD_PROCESS;

        /* number of mbs to be processed in the current job */
        s_job.i2_mb_cnt = ps_codec->s_cfg.i4_wd_mbs;

        /* job start index x */
        s_job.i2_mb_x = 0;

        /* proc base idx */
        s_job.i2_proc_base_idx = ctxt_sel ? (MAX_PROCESS_CTXT / 2) : 0;

        for (i = 0; i < (WORD32)ps_codec->s_cfg.i4_ht_mbs; i++)
        {
            /* job start index y */
            s_job.i2_mb_y = i;

            /* queue the job */
            ret = ih264_list_queue(ps_codec->pv_proc_jobq, &s_job, 1);
            if (ret != IH264_SUCCESS)
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

    return error_status;
}
