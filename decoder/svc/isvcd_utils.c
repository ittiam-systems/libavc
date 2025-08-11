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
 *  isvcd_utils.c
 *
 * @brief
 *  Contains routines that handle of start and end of pic processing
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <string.h>
#include "ih264_typedefs.h"
#include "ithread.h"
#include "ih264d_deblocking.h"
#include "ih264d_parse_slice.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_dpb_manager.h"
#include "ih264d_defs.h"
#include "isvcd_structs.h"
#include "ih264d_mem_request.h"
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_tables.h"
#include "ih264d_debug.h"
#include "ih264d_mb_utils.h"
#include "ih264d_error_handler.h"
#include "ih264d_dpb_manager.h"
#include "ih264d_utils.h"
#include "ih264d_defs.h"
#include "ih264d_tables.h"
#include "ih264d_inter_pred.h"
#include "ih264d_dpb_manager.h"
#include "iv.h"
#include "ivd.h"
#include "ih264d_format_conv.h"
#include "ih264_error.h"
#include "ih264_disp_mgr.h"
#include "ih264_buf_mgr.h"
#include "ih264d_utils.h"

WORD32 ih264d_init_dec_mb_grp(dec_struct_t *ps_dec);
/*!
**************************************************************************
* \if Function name : isvcd_free_dynamic_bufs \endif
*
* \brief
*    This function frees dynamic memory allocated by Decoder.
*
* \param ps_dec: Pointer to dec_struct_t.
*
* \return
*    Returns i4_status as returned by MemManager.
*
**************************************************************************
*/
WORD16 isvcd_free_dynamic_bufs(svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    WORD32 i;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    /* Free any avc dynamic buffers that are allocated */
    ih264d_free_dynamic_bufs(ps_dec);
    PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->pu1_crop_wnd_flag);
    PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->ps_inter_lyr_mb_prms_base);
    PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->ps_il_pred_mv_bank_buf_base);
    PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->pi2_il_residual_resample_luma_base);
    PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->pi2_il_residual_resample_chroma_base);
    PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->ps_svc_frm_mb_info);
    PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->pu2_frm_res_luma_csbp);
    PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->pu1_svc_base_mode_flag);

    memset(ps_dec->ps_pic_buf_base, 0, sizeof(struct pic_buffer_t) * (H264_MAX_REF_PICS * 2));
    for(i = 0; i < MAX_DISP_BUFS_NEW; i++)
    {
        ps_dec->apv_buf_id_pic_buf_map[i] = NULL;
    }
    return 0;
}

/*!
**************************************************************************
* \if Function name : isvcd_allocate_dynamic_bufs \endif
*
* \brief
*    This function allocates memory required by Decoder.
*
* \param ps_dec: Pointer to dec_struct_t.
*
* \return
*    Returns i4_status as returned by MemManager.
*
**************************************************************************
*/

WORD16 isvcd_allocate_dynamic_bufs(svc_dec_lyr_struct_t *ps_svc_lyr_dec)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    WORD16 i16_status = 0;
    UWORD8 uc_frmOrFld = (1 - ps_dec->ps_cur_sps->u1_frame_mbs_only_flag);
    dec_seq_params_t *ps_sps = ps_dec->ps_cur_sps;
    UWORD32 u4_total_mbs = ps_sps->u4_total_num_of_mbs << uc_frmOrFld;
    WORD32 size;
    void *pv_buf;
    void *pv_mem_ctxt = ps_dec->pv_mem_ctxt;
    size = u4_total_mbs;

    i16_status = ih264d_allocate_dynamic_bufs(ps_dec);

    if(i16_status != OK)
    {
        /* Free any dynamic buffers that are allocated */
        ih264d_free_dynamic_bufs(ps_dec);
        ps_dec->i4_error_code = IVD_MEM_ALLOC_FAILED;
        return IVD_MEM_ALLOC_FAILED;
    }
    if(u4_total_mbs == 0)
    {
        return IVD_MEM_ALLOC_FAILED;
    }

    /* Allocate frame level mb info */
    size = sizeof(dec_svc_mb_info_t) * u4_total_mbs;
    pv_buf = ps_dec->pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_svc_lyr_dec->ps_svc_frm_mb_info = pv_buf;
    memset(ps_svc_lyr_dec->ps_svc_frm_mb_info, 0, size);

    /* Allocate frame level residual luma csbp info */
    size = sizeof(UWORD16) * u4_total_mbs;
    pv_buf = ps_dec->pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_svc_lyr_dec->pu2_frm_res_luma_csbp = pv_buf;
    memset(ps_svc_lyr_dec->pu2_frm_res_luma_csbp, 0, size);
    ps_svc_lyr_dec->i4_frm_res_luma_csbp_stride = ps_dec->u2_frm_wd_in_mbs;

    /* Allocate frame level residual luma csbp info */
    size = sizeof(UWORD8) * u4_total_mbs;
    pv_buf = ps_dec->pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_svc_lyr_dec->pu1_svc_base_mode_flag = pv_buf;
    memset(ps_svc_lyr_dec->pu1_svc_base_mode_flag, 0, size);
    ps_svc_lyr_dec->i4_frm_svc_base_mode_cabac_stride = ps_dec->u2_frm_wd_in_mbs;
    ps_svc_lyr_dec->i4_frm_svc_base_mode_cabac_size = u4_total_mbs;

    /* Allocate frame level crop windows flags */
    size = sizeof(UWORD8) * u4_total_mbs;
    pv_buf = ps_dec->pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    ps_svc_lyr_dec->pu1_crop_wnd_flag = pv_buf;
    memset(ps_svc_lyr_dec->pu1_crop_wnd_flag, 0, size);

    /**********************************/
    /*Creation of Inter layer buffers */
    /**********************************/

    /* MB type buffer : one element per MB */
    size = (ps_dec->u2_frm_wd_in_mbs + 2) * (ps_dec->u2_frm_ht_in_mbs + 2) *
           sizeof(inter_lyr_mb_prms_t);
    pv_buf = ps_dec->pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, -1, size);
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_base = pv_buf;
    ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride = ps_dec->u2_frm_wd_in_mbs + 2;
    ps_svc_lyr_dec->ps_inter_lyr_mb_prms_frm_start =
        ps_svc_lyr_dec->ps_inter_lyr_mb_prms_base + 1 + ps_svc_lyr_dec->u2_inter_lyr_mb_prms_stride;

    ps_svc_lyr_dec->u4_inter_lyr_mb_prms_size = (ps_dec->u2_frm_wd_in_mbs + 2) *
                                                (ps_dec->u2_frm_ht_in_mbs + 2) *
                                                sizeof(inter_lyr_mb_prms_t);

    /* Luma Residual data at each layer : dafault 0*/
    size = ((ps_dec->u2_pic_wd + 4) * (ps_dec->u2_pic_ht + 4)) * sizeof(WORD16);
    pv_buf = ps_dec->pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svc_lyr_dec->pi2_il_residual_resample_luma_base = pv_buf;
    ps_svc_lyr_dec->u2_residual_resample_luma_stride = (ps_dec->u2_pic_wd + 4);
    ps_svc_lyr_dec->pi2_il_residual_resample_mb_luma_frm_start =
        ps_svc_lyr_dec->pi2_il_residual_resample_luma_base + 2 +
        (2 * ps_svc_lyr_dec->u2_residual_resample_luma_stride);
    ps_svc_lyr_dec->u4_residual_resample_luma_size =
        ((ps_dec->u2_pic_wd + 4) * (ps_dec->u2_pic_ht + 4)) * sizeof(WORD16);

    /* Chroma Residual data at each layer : dafault 0*/
    size = (((4 + ps_dec->u2_pic_wd) * ((4 + ps_dec->u2_pic_ht) >> 1)) * sizeof(WORD16));
    pv_buf = ps_dec->pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svc_lyr_dec->pi2_il_residual_resample_chroma_base = pv_buf;
    ps_svc_lyr_dec->u2_residual_resample_chroma_stride = (ps_dec->u2_pic_wd + 4);
    ps_svc_lyr_dec->pi2_il_residual_resample_mb_chroma_frm_start =
        ps_svc_lyr_dec->pi2_il_residual_resample_chroma_base + 2 +
        ps_svc_lyr_dec->u2_residual_resample_chroma_stride;
    ps_svc_lyr_dec->u4_residual_resample_chroma_size =
        (((4 + ps_dec->u2_pic_wd) * ((4 + ps_dec->u2_pic_ht) >> 1)) * sizeof(WORD16));

    /* mv bank buffer : 16 elements per MB: each at 4x4 block level */
    size = ((ps_dec->u2_pic_wd) * (ps_dec->u2_pic_ht >> 4)) * sizeof(mv_pred_t);
    pv_buf = ps_dec->pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svc_lyr_dec->ps_il_pred_mv_bank_buf_base = pv_buf;

    /*syntax for SVC related bin ctxt tables*/
    {
        bin_ctxt_model_t *const p_cabac_ctxt_table_t = ps_dec->p_cabac_ctxt_table_t;

        ps_svc_lyr_dec->ps_base_mode_flag = p_cabac_ctxt_table_t + CABAC_BASE_MODE_FLAG;
        ps_svc_lyr_dec->ps_motion_prediction_flag_l0 = p_cabac_ctxt_table_t + CABAC_MOT_PRED_FLAG0;
        ps_svc_lyr_dec->ps_motion_prediction_flag_l1 = p_cabac_ctxt_table_t + CABAC_MOT_PRED_FLAG1;
        ps_svc_lyr_dec->ps_residual_prediction_flag = p_cabac_ctxt_table_t + CABAC_RES_PRED_FLAG;
    }
    return (i16_status);
}

/*!
**************************************************************************
* \if Function name : isvcd_decode_pic_order_cnt \endif
*
* \brief
*    Calculates picture order count of picture.
*
* \return
*    Returns the pic order count of the picture to which current
*    Slice belongs.
*
**************************************************************************
*/
WORD32 isvcd_decode_pic_order_cnt(
    UWORD8 u1_is_idr_slice, UWORD32 u2_frame_num, pocstruct_t *ps_prev_poc, pocstruct_t *ps_cur_poc,
    dec_slice_params_t *ps_cur_slice, /*!< Pointer to current slice Params*/
    dec_pic_params_t *ps_pps, UWORD8 u1_nal_ref_idc, UWORD8 u1_bottom_field_flag,
    UWORD8 u1_field_pic_flag, WORD32 *pi4_poc, dec_struct_t *ps_dec)
{
    WORD64 i8_pic_msb;
    WORD32 i4_top_field_order_cnt = 0, i4_bottom_field_order_cnt = 0;
    dec_seq_params_t *ps_seq = ps_dec->ps_cur_sps;
    WORD32 i4_prev_frame_num_ofst;

    switch(ps_seq->u1_pic_order_cnt_type)
    {
        case 0:
            /* POC TYPE 0 */
            if(u1_is_idr_slice)
            {
                ps_prev_poc->i4_pic_order_cnt_msb = 0;
                ps_prev_poc->i4_pic_order_cnt_lsb = 0;
            }
            if(ps_prev_poc->u1_mmco_equalto5)
            {
                if(ps_prev_poc->u1_bot_field != 1)
                {
                    ps_prev_poc->i4_pic_order_cnt_msb = 0;
                    ps_prev_poc->i4_pic_order_cnt_lsb = ps_prev_poc->i4_top_field_order_count;
                }
                else
                {
                    ps_prev_poc->i4_pic_order_cnt_msb = 0;
                    ps_prev_poc->i4_pic_order_cnt_lsb = 0;
                }
            }

            if((ps_cur_poc->i4_pic_order_cnt_lsb < ps_prev_poc->i4_pic_order_cnt_lsb) &&
               ((ps_prev_poc->i4_pic_order_cnt_lsb - ps_cur_poc->i4_pic_order_cnt_lsb) >=
                (ps_seq->i4_max_pic_order_cntLsb >> 1)))
            {
                i8_pic_msb =
                    (WORD64) ps_prev_poc->i4_pic_order_cnt_msb + ps_seq->i4_max_pic_order_cntLsb;
            }
            else if((ps_cur_poc->i4_pic_order_cnt_lsb > ps_prev_poc->i4_pic_order_cnt_lsb) &&
                    ((ps_cur_poc->i4_pic_order_cnt_lsb - ps_prev_poc->i4_pic_order_cnt_lsb) >=
                     (ps_seq->i4_max_pic_order_cntLsb >> 1)))
            {
                i8_pic_msb =
                    (WORD64) ps_prev_poc->i4_pic_order_cnt_msb - ps_seq->i4_max_pic_order_cntLsb;
            }
            else
            {
                i8_pic_msb = ps_prev_poc->i4_pic_order_cnt_msb;
            }

            if(!u1_field_pic_flag || !u1_bottom_field_flag)
            {
                WORD64 i8_result = i8_pic_msb + ps_cur_poc->i4_pic_order_cnt_lsb;
                if(IS_OUT_OF_RANGE_S32(i8_result))
                {
                    return ERROR_INV_POC;
                }
                i4_top_field_order_cnt = (WORD32) i8_result;
            }

            if(!u1_field_pic_flag)
            {
                WORD64 i8_result =
                    (WORD64) i4_top_field_order_cnt + ps_cur_poc->i4_delta_pic_order_cnt_bottom;
                if(IS_OUT_OF_RANGE_S32(i8_result))
                {
                    return ERROR_INV_POC;
                }
                i4_bottom_field_order_cnt = (WORD32) i8_result;
            }
            else if(u1_bottom_field_flag)
            {
                WORD64 i8_result = i8_pic_msb + ps_cur_poc->i4_pic_order_cnt_lsb;
                if(IS_OUT_OF_RANGE_S32(i8_result))
                {
                    return ERROR_INV_POC;
                }
                i4_bottom_field_order_cnt = (WORD32) i8_result;
            }

            if(IS_OUT_OF_RANGE_S32(i8_pic_msb))
            {
                return ERROR_INV_POC;
            }
            ps_cur_poc->i4_pic_order_cnt_msb = (WORD32) i8_pic_msb;
            break;

        case 1:
        {
            /* POC TYPE 1 */
            UWORD8 i;
            WORD32 prev_frame_num;
            WORD32 frame_num_ofst;
            WORD32 abs_frm_num;
            WORD32 poc_cycle_cnt, frame_num_in_poc_cycle;
            WORD64 i8_expected_delta_poc_cycle;
            WORD32 expected_poc;
            WORD64 i8_result;

            prev_frame_num = (WORD32) ps_cur_slice->u2_frame_num;
            if(!u1_is_idr_slice)
            {
                if(ps_cur_slice->u1_mmco_equalto5)
                {
                    prev_frame_num = 0;
                    i4_prev_frame_num_ofst = 0;
                }
                else
                {
                    i4_prev_frame_num_ofst = ps_prev_poc->i4_prev_frame_num_ofst;
                }
            }
            else
                i4_prev_frame_num_ofst = 0;

            /* 1. Derivation for FrameNumOffset */
            if(u1_is_idr_slice)
            {
                frame_num_ofst = 0;
                ps_cur_poc->i4_delta_pic_order_cnt[0] = 0;
                ps_cur_poc->i4_delta_pic_order_cnt[1] = 0;
            }
            else if(prev_frame_num > ((WORD32) u2_frame_num))
            {
                WORD64 i8_result =
                    i4_prev_frame_num_ofst + (WORD64) ps_seq->u2_u4_max_pic_num_minus1 + 1;
                if(IS_OUT_OF_RANGE_S32(i8_result))
                {
                    return ERROR_INV_FRAME_NUM;
                }
                frame_num_ofst = (WORD32) i8_result;
            }
            else
                frame_num_ofst = i4_prev_frame_num_ofst;

            /* 2. Derivation for absFrameNum */
            if(0 != ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle)
            {
                WORD64 i8_result = frame_num_ofst + (WORD64) u2_frame_num;
                if(IS_OUT_OF_RANGE_S32(i8_result))
                {
                    return ERROR_INV_FRAME_NUM;
                }
                abs_frm_num = (WORD32) i8_result;
            }
            else
                abs_frm_num = 0;
            if((u1_nal_ref_idc == 0) && (abs_frm_num > 0)) abs_frm_num = abs_frm_num - 1;

            /* 4. expectedDeltaPerPicOrderCntCycle is derived as */
            i8_expected_delta_poc_cycle = 0;
            for(i = 0; i < ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle; i++)
            {
                i8_expected_delta_poc_cycle += ps_seq->i4_ofst_for_ref_frame[i];
            }

            /* 3. When absFrameNum > 0, picOrderCntCycleCnt and
            frame_num_in_poc_cycle are derived as : */
            /* 5. expectedPicOrderCnt is derived as : */
            if(abs_frm_num > 0)
            {
                poc_cycle_cnt =
                    DIV((abs_frm_num - 1), ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle);
                frame_num_in_poc_cycle =
                    MOD((abs_frm_num - 1), ps_seq->u1_num_ref_frames_in_pic_order_cnt_cycle);

                i8_result = poc_cycle_cnt * i8_expected_delta_poc_cycle;

                for(i = 0; i <= frame_num_in_poc_cycle; i++)
                {
                    i8_result = i8_result + ps_seq->i4_ofst_for_ref_frame[i];
                }

                if(IS_OUT_OF_RANGE_S32(i8_result)) return ERROR_INV_POC;

                expected_poc = (WORD32) i8_result;
            }
            else
                expected_poc = 0;

            if(u1_nal_ref_idc == 0)
            {
                i8_result = (WORD64) expected_poc + ps_seq->i4_ofst_for_non_ref_pic;

                if(IS_OUT_OF_RANGE_S32(i8_result)) return ERROR_INV_POC;

                expected_poc = (WORD32) i8_result;
            }

            /* 6. TopFieldOrderCnt or BottomFieldOrderCnt are derived as */
            if(!u1_field_pic_flag)
            {
                i8_result = (WORD64) expected_poc + ps_cur_poc->i4_delta_pic_order_cnt[0];

                if(IS_OUT_OF_RANGE_S32(i8_result)) return ERROR_INV_POC;
                i4_top_field_order_cnt = (WORD32) i8_result;

                i8_result = (WORD64) i4_top_field_order_cnt +
                            ps_seq->i4_ofst_for_top_to_bottom_field +
                            ps_cur_poc->i4_delta_pic_order_cnt[1];

                if(IS_OUT_OF_RANGE_S32(i8_result)) return ERROR_INV_POC;
                i4_bottom_field_order_cnt = (WORD32) i8_result;
            }
            else if(!u1_bottom_field_flag)
            {
                i8_result = (WORD64) expected_poc + ps_cur_poc->i4_delta_pic_order_cnt[0];

                if(IS_OUT_OF_RANGE_S32(i8_result)) return ERROR_INV_POC;
                i4_top_field_order_cnt = (WORD32) i8_result;
            }
            else
            {
                i8_result = (WORD64) expected_poc + ps_seq->i4_ofst_for_top_to_bottom_field +
                            ps_cur_poc->i4_delta_pic_order_cnt[0];

                if(IS_OUT_OF_RANGE_S32(i8_result)) return ERROR_INV_POC;
                i4_bottom_field_order_cnt = (WORD32) i8_result;
            }
            /* Copy the current POC info into Previous POC structure */
            ps_cur_poc->i4_prev_frame_num_ofst = frame_num_ofst;
        }

        break;
        case 2:
        {
            /* POC TYPE 2 */
            WORD32 prev_frame_num;
            WORD32 frame_num_ofst;
            WORD32 tmp_poc;

            prev_frame_num = (WORD32) ps_cur_slice->u2_frame_num;
            if(!u1_is_idr_slice)
            {
                if(ps_cur_slice->u1_mmco_equalto5)
                {
                    prev_frame_num = 0;
                    i4_prev_frame_num_ofst = 0;
                }
                else
                    i4_prev_frame_num_ofst = ps_prev_poc->i4_prev_frame_num_ofst;
            }
            else
                i4_prev_frame_num_ofst = 0;

            /* 1. Derivation for FrameNumOffset */
            if(u1_is_idr_slice)
            {
                frame_num_ofst = 0;
                ps_cur_poc->i4_delta_pic_order_cnt[0] = 0;
                ps_cur_poc->i4_delta_pic_order_cnt[1] = 0;
            }
            else if(prev_frame_num > ((WORD32) u2_frame_num))
            {
                WORD64 i8_result =
                    i4_prev_frame_num_ofst + (WORD64) ps_seq->u2_u4_max_pic_num_minus1 + 1;
                if(IS_OUT_OF_RANGE_S32(i8_result))
                {
                    return ERROR_INV_FRAME_NUM;
                }
                frame_num_ofst = (WORD32) i8_result;
            }
            else
                frame_num_ofst = i4_prev_frame_num_ofst;

            /* 2. Derivation for tempPicOrderCnt */
            if(u1_is_idr_slice)
                tmp_poc = 0;
            else if(u1_nal_ref_idc == 0)
            {
                WORD64 i8_result = ((frame_num_ofst + (WORD64) u2_frame_num) << 1) - 1;
                if(IS_OUT_OF_RANGE_S32(i8_result))
                {
                    return ERROR_INV_POC;
                }
                tmp_poc = (WORD32) i8_result;
            }
            else
            {
                WORD64 i8_result = (frame_num_ofst + (WORD64) u2_frame_num) << 1;
                if(IS_OUT_OF_RANGE_S32(i8_result))
                {
                    return ERROR_INV_POC;
                }
                tmp_poc = (WORD32) i8_result;
            }

            /* 6. TopFieldOrderCnt or BottomFieldOrderCnt are derived as */
            if(!u1_field_pic_flag)
            {
                i4_top_field_order_cnt = tmp_poc;
                i4_bottom_field_order_cnt = tmp_poc;
            }
            else if(!u1_bottom_field_flag)
                i4_top_field_order_cnt = tmp_poc;
            else
                i4_bottom_field_order_cnt = tmp_poc;

            /* Copy the current POC info into Previous POC structure */
            ps_prev_poc->i4_prev_frame_num_ofst = frame_num_ofst;
            ps_cur_poc->i4_prev_frame_num_ofst = frame_num_ofst;
        }
        break;
        default:
            return ERROR_INV_POC_TYPE_T;
            break;
    }

    if(!u1_field_pic_flag)  // or a complementary field pair
    {
        *pi4_poc = MIN(i4_top_field_order_cnt, i4_bottom_field_order_cnt);
        ps_pps->i4_top_field_order_cnt = i4_top_field_order_cnt;
        ps_pps->i4_bottom_field_order_cnt = i4_bottom_field_order_cnt;
    }
    else if(!u1_bottom_field_flag)
    {
        *pi4_poc = i4_top_field_order_cnt;
        ps_pps->i4_top_field_order_cnt = i4_top_field_order_cnt;
    }
    else
    {
        *pi4_poc = i4_bottom_field_order_cnt;
        ps_pps->i4_bottom_field_order_cnt = i4_bottom_field_order_cnt;
    }

    ps_pps->i4_avg_poc = *pi4_poc;

    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_decode_gaps_in_frame_num */
/*                                                                           */
/*  Description   : This function decodes gaps in frame number               */
/*                                                                           */
/*  Inputs        : ps_dec          Decoder parameters                       */
/*                  u2_frame_num   current frame number                     */
/*                                                                           */
/*  Globals       : None                                                     */
/*  Processing    : This functionality needs to be implemented               */
/*  Outputs       : None                                                     */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : Not implemented                                          */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 05 2002   NS              Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_decode_gaps_in_frame_num(dec_struct_t *ps_dec, UWORD16 u2_frame_num)
{
    UWORD32 u4_next_frm_num, u4_start_frm_num;
    UWORD32 u4_max_frm_num;
    pocstruct_t s_tmp_poc;
    WORD32 i4_poc;
    dec_slice_params_t *ps_cur_slice;

    dec_pic_params_t *ps_pic_params;
    WORD8 i1_gap_idx;
    WORD32 *i4_gaps_start_frm_num;
    dpb_manager_t *ps_dpb_mgr;
    WORD8 *pi1_gaps_per_seq;
    WORD32 ret;

    ps_cur_slice = ps_dec->ps_cur_slice;
    if(ps_cur_slice->u1_field_pic_flag)
    {
        if(ps_dec->u2_prev_ref_frame_num == u2_frame_num) return 0;
    }

    u4_next_frm_num = ps_dec->u2_prev_ref_frame_num + 1;
    u4_max_frm_num = ps_dec->ps_cur_sps->u2_u4_max_pic_num_minus1 + 1;

    if(u4_next_frm_num >= u4_max_frm_num)
    {
        u4_next_frm_num -= u4_max_frm_num;
    }

    if(u4_next_frm_num == u2_frame_num)
    {
        return (0);
    }

    if((ps_dec->u1_nal_unit_type == IDR_SLICE_NAL) && (u4_next_frm_num >= u2_frame_num))
    {
        return (0);
    }
    u4_start_frm_num = u4_next_frm_num;

    s_tmp_poc.i4_pic_order_cnt_lsb = 0;
    s_tmp_poc.i4_delta_pic_order_cnt_bottom = 0;
    s_tmp_poc.i4_pic_order_cnt_lsb = 0;
    s_tmp_poc.i4_delta_pic_order_cnt_bottom = 0;
    s_tmp_poc.i4_delta_pic_order_cnt[0] = 0;
    s_tmp_poc.i4_delta_pic_order_cnt[1] = 0;

    ps_cur_slice = ps_dec->ps_cur_slice;
    ps_pic_params = ps_dec->ps_cur_pps;

    ps_dpb_mgr = ps_dec->ps_dpb_mgr;

    /* Find a empty slot to store gap seqn info */
    i4_gaps_start_frm_num = ps_dpb_mgr->ai4_gaps_start_frm_num;
    for(i1_gap_idx = 0; i1_gap_idx < MAX_FRAMES; i1_gap_idx++)
    {
        if(INVALID_FRAME_NUM == i4_gaps_start_frm_num[i1_gap_idx]) break;
    }
    if(MAX_FRAMES == i1_gap_idx)
    {
        UWORD32 i4_error_code;
        i4_error_code = ERROR_DBP_MANAGER_T;
        return i4_error_code;
    }

    i4_poc = 0;
    i4_gaps_start_frm_num[i1_gap_idx] = u4_start_frm_num;
    ps_dpb_mgr->ai4_gaps_end_frm_num[i1_gap_idx] = u2_frame_num - 1;
    pi1_gaps_per_seq = ps_dpb_mgr->ai1_gaps_per_seq;
    pi1_gaps_per_seq[i1_gap_idx] = 0;
    while(u4_next_frm_num != u2_frame_num)
    {
        ih264d_delete_nonref_nondisplay_pics(ps_dpb_mgr);
        if(ps_pic_params->ps_sps->u1_pic_order_cnt_type)
        {
            /* allocate a picture buffer and insert it as ST node */
            ret =
                isvcd_decode_pic_order_cnt(0, u4_next_frm_num, &ps_dec->s_prev_pic_poc, &s_tmp_poc,
                                           ps_cur_slice, ps_pic_params, 1, 0, 0, &i4_poc, ps_dec);
            if(ret != OK) return ret;

            /* Display seq no calculations */
            if(i4_poc >= ps_dec->i4_max_poc) ps_dec->i4_max_poc = i4_poc;
            /* IDR Picture or POC wrap around */
            if(i4_poc == 0)
            {
                WORD64 i8_temp;
                i8_temp = (WORD64) ps_dec->i4_prev_max_display_seq + ps_dec->i4_max_poc +
                          ps_dec->u1_max_dec_frame_buffering + 1;
                /*If i4_prev_max_display_seq overflows integer range, reset it */
                ps_dec->i4_prev_max_display_seq =
                    IS_OUT_OF_RANGE_S32(i8_temp) ? 0 : (WORD32) i8_temp;
                ps_dec->i4_max_poc = 0;
            }

            ps_cur_slice->u1_mmco_equalto5 = 0;
            ps_cur_slice->u2_frame_num = u4_next_frm_num;
        }

        if(ps_dpb_mgr->i1_poc_buf_id_entries >= ps_dec->u1_max_dec_frame_buffering)
        {
            ret = ih264d_assign_display_seq(ps_dec);
            if(ret != OK) return ret;
        }

        {
            WORD64 i8_display_poc;
            i8_display_poc = (WORD64) ps_dec->i4_prev_max_display_seq + i4_poc;
            if(IS_OUT_OF_RANGE_S32(i8_display_poc))
            {
                ps_dec->i4_prev_max_display_seq = 0;
            }
        }
        ret = ih264d_insert_pic_in_display_list(ps_dec->ps_dpb_mgr, (WORD8) DO_NOT_DISP,
                                                (WORD32) (ps_dec->i4_prev_max_display_seq + i4_poc),
                                                u4_next_frm_num);
        if(ret != OK) return ret;

        pi1_gaps_per_seq[i1_gap_idx]++;
        ret = ih264d_do_mmco_for_gaps(ps_dpb_mgr, ps_dec->ps_cur_sps->u1_num_ref_frames);
        if(ret != OK) return ret;

        ih264d_delete_nonref_nondisplay_pics(ps_dpb_mgr);

        u4_next_frm_num++;
        if(u4_next_frm_num >= u4_max_frm_num)
        {
            u4_next_frm_num -= u4_max_frm_num;
        }
    }

    return OK;
}

/*!
****************************************************************************
*                                                                           
*  \if Function name : isvcd_init_dpb_ref_bufs \endif
*
*  \brief
*    Initializes the reference buffers
*
*  \return
*    None
*
*  \note
*    This function is called to initialize the reference buffers.
****************************************************************************
*/

void isvcd_init_dpb_ref_bufs(dec_struct_t *ps_dec)
{
    UWORD8 i;
    struct pic_buffer_t *ps_init_dpb;
    ps_init_dpb = ps_dec->ps_dpb_mgr->ps_init_dpb[0][0];
    for(i = 0; i < 2 * MAX_REF_BUFS; i++)
    {
        memset(ps_init_dpb, 0, sizeof(struct pic_buffer_t));
        ps_init_dpb->pu1_buf1 = NULL;
        ps_init_dpb->u1_long_term_frm_idx = MAX_REF_BUFS + 1;
        ps_dec->ps_dpb_mgr->ps_init_dpb[0][i] = ps_init_dpb;
        ps_dec->ps_dpb_mgr->ps_mod_dpb[0][i] = ps_init_dpb;
        ps_init_dpb++;
    }

    ps_init_dpb = ps_dec->ps_dpb_mgr->ps_init_dpb[1][0];
    for(i = 0; i < 2 * MAX_REF_BUFS; i++)
    {
        memset(ps_init_dpb, 0, sizeof(struct pic_buffer_t));
        ps_init_dpb->pu1_buf1 = NULL;
        ps_init_dpb->u1_long_term_frm_idx = MAX_REF_BUFS + 1;
        ps_dec->ps_dpb_mgr->ps_init_dpb[1][i] = ps_init_dpb;
        ps_dec->ps_dpb_mgr->ps_mod_dpb[1][i] = ps_init_dpb;
        ps_init_dpb++;
    }
}

/*!
**************************************************************************
* \if Function name : isvcd_init_pic \endif
*
* \brief
*    Initializes the picture.
*
* \return
*    0 on Success and Error code otherwise
*
* \note
*    This function is called when first slice of the
*    NON -IDR picture is encountered.
**************************************************************************
*/
WORD32 isvcd_init_pic(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_frame_num, WORD32 i4_poc,
                      dec_pic_params_t *ps_pps)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_seq_params_t *ps_seq = ps_dec->ps_cur_sps;
    prev_seq_params_t *ps_prev_seq_params = &ps_dec->s_prev_seq_params;
    WORD32 ret;

    ps_dec->ps_cur_slice->u2_frame_num = u2_frame_num;
    ps_dec->ps_cur_slice->i4_poc = i4_poc;
    ps_dec->ps_cur_pps = ps_pps;
    ps_dec->ps_cur_pps->pv_codec_handle = ps_dec;

    ps_dec->ps_dpb_mgr->i4_max_frm_num = ps_seq->u2_u4_max_pic_num_minus1 + 1;

    ps_dec->ps_dpb_mgr->u2_pic_ht = ps_dec->u2_pic_ht;
    ps_dec->ps_dpb_mgr->u2_pic_wd = ps_dec->u2_pic_wd;
    ps_dec->i4_pic_type = NA_SLICE;
    ps_dec->i4_frametype = IV_NA_FRAME;
    ps_dec->i4_content_type = IV_CONTENTTYPE_NA;

    /*--------------------------------------------------------------------*/
    /* Get the value of MaxMbAddress and frmheight in Mbs                 */
    /*--------------------------------------------------------------------*/
    ps_seq->u4_max_mb_addr =
        ((UWORD32)ps_seq->u2_frm_wd_in_mbs *
         ((UWORD32)ps_dec->u2_pic_ht >> (4 + ps_dec->ps_cur_slice->u1_field_pic_flag))) -
        1;
    ps_dec->u2_frm_ht_in_mbs = (ps_dec->u2_pic_ht >> (4 + ps_dec->ps_cur_slice->u1_field_pic_flag));

    /***************************************************************************/
    /* If change in Level or the required PicBuffers i4_size is more than the  */
    /* current one FREE the current PicBuffers and allocate affresh            */
    /***************************************************************************/
    if(!ps_dec->u1_init_dec_flag)
    {
        ps_dec->u1_max_dec_frame_buffering = ih264d_get_dpb_size(ps_seq);

        ps_dec->i4_display_delay = ps_dec->u1_max_dec_frame_buffering;
        if((1 == ps_seq->u1_vui_parameters_present_flag) &&
           (1 == ps_seq->s_vui.u1_bitstream_restriction_flag))
        {
            if(ps_seq->u1_frame_mbs_only_flag == 1)
                ps_dec->i4_display_delay = ps_seq->s_vui.u4_num_reorder_frames + 1;
            else
                ps_dec->i4_display_delay = ps_seq->s_vui.u4_num_reorder_frames * 2 + 2;
        }

        if(IVD_DECODE_FRAME_OUT == ps_dec->e_frm_out_mode) ps_dec->i4_display_delay = 0;

        if(ps_dec->u4_share_disp_buf == 0)
        {
            if(ps_seq->u1_frame_mbs_only_flag == 1)
                ps_dec->u1_pic_bufs = ps_dec->i4_display_delay + ps_seq->u1_num_ref_frames + 1;
            else
                ps_dec->u1_pic_bufs = ps_dec->i4_display_delay + ps_seq->u1_num_ref_frames * 2 + 2;
        }
        else
        {
            ps_dec->u1_pic_bufs = (WORD32) ps_dec->u4_num_disp_bufs;
        }

        /* Ensure at least two buffers are allocated */
        ps_dec->u1_pic_bufs = MAX(ps_dec->u1_pic_bufs, 2);

        if(ps_dec->u4_share_disp_buf == 0)
            ps_dec->u1_pic_bufs = MIN(ps_dec->u1_pic_bufs, (H264_MAX_REF_PICS * 2));

        ps_dec->u1_max_dec_frame_buffering =
            MIN(ps_dec->u1_max_dec_frame_buffering, ps_dec->u1_pic_bufs);

        /* Temporary hack to run Tractor Cav/Cab/MbAff Profiler streams  also for
         * CAFI1_SVA_C.264 in conformance*/
        if(ps_dec->u1_init_dec_flag)
        {
            ih264d_release_pics_in_dpb((void *) ps_dec, ps_dec->u1_pic_bufs);
            ih264d_release_display_bufs(ps_dec);
            ih264d_reset_ref_bufs(ps_dec->ps_dpb_mgr);
        }

        /*********************************************************************/
        /* Configuring decoder parameters based on level and then            */
        /* fresh pointer initialisation in decoder scratch and state buffers */
        /*********************************************************************/
        if(!ps_dec->u1_init_dec_flag || ((ps_seq->u1_level_idc < H264_LEVEL_3_0) ^
                                         (ps_prev_seq_params->u1_level_idc < H264_LEVEL_3_0)))
        {
            ret = ih264d_init_dec_mb_grp(ps_dec);
            if(ret != OK) return ret;
        }

        ret = isvcd_allocate_dynamic_bufs(ps_svc_lyr_dec);

        if(ret != OK)
        {
            /* Free any dynamic buffers that are allocated */
            isvcd_free_dynamic_bufs(ps_svc_lyr_dec);
            ps_dec->i4_error_code = IVD_MEM_ALLOC_FAILED;
            return IVD_MEM_ALLOC_FAILED;
        }

        ih264d_init_ref_bufs((dpb_manager_t *)ps_dec->ps_dpb_mgr);
        isvcd_init_dpb_ref_bufs(ps_dec);

        ret = ih264d_create_pic_buffers(ps_dec->u1_pic_bufs, ps_dec);
        if(ret != OK) return ret;

        ret = ih264d_create_mv_bank(ps_dec, ps_dec->u2_pic_wd, ps_dec->u2_pic_ht);
        if(ret != OK) return ret;

        /* In shared mode, set all of them as used by display */
        if(ps_dec->u4_share_disp_buf == 1)
        {
            WORD32 i;

            for(i = 0; i < ps_dec->u1_pic_bufs; i++)
            {
                ih264_buf_mgr_set_status((buf_mgr_t *) ps_dec->pv_pic_buf_mgr, i, BUF_MGR_IO);
            }
        }

        ps_dec->u1_init_dec_flag = 1;
        ps_prev_seq_params->u2_frm_wd_in_mbs = ps_seq->u2_frm_wd_in_mbs;
        ps_prev_seq_params->u1_level_idc = ps_seq->u1_level_idc;
        ps_prev_seq_params->u1_profile_idc = ps_seq->u1_profile_idc;
        ps_prev_seq_params->u2_frm_ht_in_mbs = ps_seq->u2_frm_ht_in_mbs;
        ps_prev_seq_params->u1_frame_mbs_only_flag = ps_seq->u1_frame_mbs_only_flag;
        ps_prev_seq_params->u1_direct_8x8_inference_flag = ps_seq->u1_direct_8x8_inference_flag;

        ps_dec->i4_cur_display_seq = 0;
        ps_dec->i4_prev_max_display_seq = 0;
        ps_dec->i4_max_poc = 0;

        {
            /* 0th entry of CtxtIncMbMap will be always be containing default values
            for CABAC context representing MB not available */
            ctxt_inc_mb_info_t *p_DefCtxt = ps_dec->p_ctxt_inc_mb_map - 1;
            UWORD8 *pu1_temp;
            WORD8 i;
            p_DefCtxt->u1_mb_type = CAB_SKIP;

            p_DefCtxt->u1_cbp = 0x0f;
            p_DefCtxt->u1_intra_chroma_pred_mode = 0;

            p_DefCtxt->u1_yuv_dc_csbp = 0x7;

            p_DefCtxt->u1_transform8x8_ctxt = 0;

            pu1_temp = (UWORD8 *) p_DefCtxt->i1_ref_idx;
            for(i = 0; i < 4; i++, pu1_temp++) (*pu1_temp) = 0;
            pu1_temp = (UWORD8 *) p_DefCtxt->u1_mv;
            for(i = 0; i < 16; i++, pu1_temp++) (*pu1_temp) = 0;
            ps_dec->ps_def_ctxt_mb_info = p_DefCtxt;
        }
    }
    /* reset DBP commands read u4_flag */
    ps_dec->ps_dpb_cmds->u1_dpb_commands_read = 0;

    return OK;
}
