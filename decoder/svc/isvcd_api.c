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
 *  isvcd_api.c
 *
 * @brief
 *  Contains all the API related functions
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - api_check_struct_sanity()
 *  - isvcd_set_processor()
 *  - isvcd_init_decoder()
 *  - isvcd_nal_parse_ctxt_free()
 *  - isvcd_residual_resample_ctxt_free()
 *  - isvcd_intra_resample_ctxt_free()
 *  - isvcd_mode_mv_resample_ctxt_free()
 *  - isvcd_free_static_bufs()
 *  - isvcd_nal_parse_ctxt_create()
 *  - isvcd_intra_resample_ctxt_create()
 *  - isvcd_residual_resample_ctxt_create()
 *  - isvcd_mode_mv_resample_ctxt_create()
 *  - isvcd_allocate_static_bufs()
 *  - isvcd_create()
 *  - isvcd_update_dqid()
 *  - isvcd_detect_res_change()
 *  - isvcd_parse_ref_pic_list_modify()
 *  - isvcd_parse_slice_hdr_refdq_id()
 *  - isvcd_get_ref_lyr_dqid()
 *  - isvcd_conceal_node_params()
 *  - isvcd_refine_dep_list()
 *  - isvcd_dec_non_vcl()
 *  - isvcd_seq_hdr_dec()
 *  - isvcd_pre_parse_refine_au()
 *  - isvcd_video_decode()
 *  - isvcd_set_display_frame()
 *  - isvcd_set_flush_mode()
 *  - isvcd_get_status()
 *  - isvcd_get_buf_info()
 *  - isvcd_set_params()
 *  - isvcd_set_target_layer()
 *  - isvcd_set_default_params()
 *  - isvcd_delete()
 *  - isvcd_reset()
 *  - isvcd_ctl()
 *  - isvcd_rel_display_frame()
 *  - isvcd_set_degrade()
 *  - isvcd_get_frame_dimensions()
 *  - isvcd_get_vui_params()
 *  - isvcd_get_sei_mdcv_params()
 *  - isvcd_get_sei_cll_params()
 *  - isvcd_get_sei_ave_params()
 *  - isvcd_get_sei_ccv_params()
 *  - isvcd_set_num_cores()
 *  - isvcd_fill_output_struct_from_context()
 *  - isvcd_api_function()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <string.h>
#include <limits.h>
#include <stddef.h>
#include <assert.h>
#include "ih264_defs.h"
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_tables.h"
#include "iv.h"
#include "ivd.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "ih264_debug.h"
#include "ih264d_inter_pred.h"
#include "isvcd_structs.h"
#include "ih264d_nal.h"
#include "ih264d_error_handler.h"
#include "ithread.h"
#include "ih264d_parse_slice.h"
#include "ih264d_function_selector.h"
#include "ih264_error.h"
#include "ih264_disp_mgr.h"
#include "ih264_buf_mgr.h"
#include "ih264d_deblocking.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_parse_cabac.h"
#include "ih264d_process_pslice.h"
#include "isvcd_process_epslice.h"
#include "ih264d_utils.h"
#include "ih264d_api_utils.h"
#include "ih264d_format_conv.h"
#include "ih264d_parse_headers.h"
#include "ih264d_thread_compute_bs.h"
#include "isvcd_utils.h"
#include "isvcd.h"
#include "isvcd_mode_mv_resamp.h"
#include "isvcd_parse_headers.h"
#include "isvcd_thread_compute_bs.h"
#include "isvcd_function_selector.h"
/*********************/
/* Codec Versioning  */
/*********************/
// Move this to where it is used
#define CODEC_NAME "H264VDEC"
#define CODEC_RELEASE_TYPE "production"
#define CODEC_RELEASE_VER "05.00"
#define CODEC_VENDOR "ITTIAM"
#define MAXVERSION_STRLEN 511
#ifdef ANDROID
#define VERSION(version_string, codec_name, codec_release_type, codec_release_ver, codec_vendor)  \
    snprintf(version_string, MAXVERSION_STRLEN, "@(#)Id:%s_%s Ver:%s Released by %s", codec_name, \
             codec_release_type, codec_release_ver, codec_vendor)
#else
#define VERSION(version_string, codec_name, codec_release_type, codec_release_ver, codec_vendor)  \
    snprintf(version_string, MAXVERSION_STRLEN,                                                   \
             "@(#)Id:%s_%s Ver:%s Released by %s Build: %s @ %s", codec_name, codec_release_type, \
             codec_release_ver, codec_vendor, __DATE__, __TIME__)
#endif

#define MIN_IN_BUFS 1
#define MIN_OUT_BUFS_420 3
#define MIN_OUT_BUFS_422ILE 1
#define MIN_OUT_BUFS_RGB565 1
#define MIN_OUT_BUFS_420SP 2

#define NUM_FRAMES_LIMIT_ENABLED 0

#if NUM_FRAMES_LIMIT_ENABLED
#define NUM_FRAMES_LIMIT 10000
#else
#define NUM_FRAMES_LIMIT 0x7FFFFFFF
#endif
WORD32 ih264d_get_version(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);
WORD32 ih264d_parse_sei(dec_struct_t *ps_dec, dec_bit_stream_t *ps_bitstrm);
WORD32 check_app_out_buf_size(dec_struct_t *ps_dec);
UWORD32 ih264d_get_extra_mem_external(UWORD32 width, UWORD32 height);
WORD32 isvcd_get_frame_dimensions(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);
WORD32 isvcd_get_vui_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

WORD32 isvcd_get_sei_mdcv_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);
WORD32 isvcd_get_sei_cll_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);
WORD32 isvcd_get_sei_ave_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);
WORD32 isvcd_get_sei_ccv_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);
WORD32 isvcd_set_num_cores(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

WORD32 ih264d_deblock_display(dec_struct_t *ps_dec);
WORD32 ih264d_get_display_frame(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op);

void ih264d_signal_decode_thread(dec_struct_t *ps_dec);

void ih264d_signal_bs_deblk_thread(dec_struct_t *ps_dec);
void ih264d_decode_picture_thread(dec_struct_t *ps_dec);
void isvcd_decode_picture_thread(svc_dec_lyr_struct_t *ps_svc_lyr_dec);

WORD32 isvcd_set_degrade(iv_obj_t *ps_codec_obj, void *pv_api_ip, void *pv_api_op);

void isvcd_fill_output_struct_from_context(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                           ivd_video_decode_op_t *ps_dec_op);

static IV_API_CALL_STATUS_T api_check_struct_sanity(iv_obj_t *ps_handle, void *pv_api_ip,
                                                    void *pv_api_op)
{
    IVD_API_COMMAND_TYPE_T e_cmd;
    UWORD32 *pu4_api_ip;
    UWORD32 *pu4_api_op;
    UWORD32 i;

    if(NULL == pv_api_op) return (IV_FAIL);

    if(NULL == pv_api_ip) return (IV_FAIL);

    pu4_api_ip = (UWORD32 *) pv_api_ip;
    pu4_api_op = (UWORD32 *) pv_api_op;
    e_cmd = *(pu4_api_ip + 1);

    /* error checks on handle */
    switch((WORD32) e_cmd)
    {
        case IVD_CMD_CREATE:
            break;

        case IVD_CMD_REL_DISPLAY_FRAME:
        case IVD_CMD_SET_DISPLAY_FRAME:
        case IVD_CMD_GET_DISPLAY_FRAME:
        case IVD_CMD_VIDEO_DECODE:
        case IVD_CMD_DELETE:
        case IVD_CMD_VIDEO_CTL:
            if(ps_handle == NULL)
            {
                *(pu4_api_op + 1) |= 1 << IVD_UNSUPPORTEDPARAM;
                *(pu4_api_op + 1) |= IVD_HANDLE_NULL;
                return IV_FAIL;
            }

            if(ps_handle->u4_size != sizeof(iv_obj_t))
            {
                *(pu4_api_op + 1) |= 1 << IVD_UNSUPPORTEDPARAM;
                *(pu4_api_op + 1) |= IVD_HANDLE_STRUCT_SIZE_INCORRECT;
                return IV_FAIL;
            }

            if(ps_handle->pv_fxns != isvcd_api_function)
            {
                *(pu4_api_op + 1) |= 1 << IVD_UNSUPPORTEDPARAM;
                *(pu4_api_op + 1) |= IVD_INVALID_HANDLE_NULL;
                return IV_FAIL;
            }

            if(ps_handle->pv_codec_handle == NULL)
            {
                *(pu4_api_op + 1) |= 1 << IVD_UNSUPPORTEDPARAM;
                *(pu4_api_op + 1) |= IVD_INVALID_HANDLE_NULL;
                return IV_FAIL;
            }
            break;
        default:
            *(pu4_api_op + 1) |= 1 << IVD_UNSUPPORTEDPARAM;
            *(pu4_api_op + 1) |= IVD_INVALID_API_CMD;
            return IV_FAIL;
    }

    switch((WORD32) e_cmd)
    {
        case IVD_CMD_CREATE:
        {
            isvcd_create_ip_t *ps_ip = (isvcd_create_ip_t *) pv_api_ip;
            isvcd_create_op_t *ps_op = (isvcd_create_op_t *) pv_api_op;

            ps_op->s_ivd_create_op_t.u4_error_code = 0;

            if((ps_ip->s_ivd_create_ip_t.u4_size > sizeof(isvcd_create_ip_t)) ||
               (ps_ip->s_ivd_create_ip_t.u4_size < sizeof(ivd_create_ip_t)))
            {
                ps_op->s_ivd_create_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_create_op_t.u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                H264_DEC_DEBUG_PRINT("\n");
                return (IV_FAIL);
            }

            if((ps_op->s_ivd_create_op_t.u4_size != sizeof(isvcd_create_op_t)) &&
               (ps_op->s_ivd_create_op_t.u4_size != sizeof(ivd_create_op_t)))
            {
                ps_op->s_ivd_create_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_create_op_t.u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                H264_DEC_DEBUG_PRINT("\n");
                return (IV_FAIL);
            }

            if((ps_ip->s_ivd_create_ip_t.e_output_format != IV_YUV_420P) &&
               (ps_ip->s_ivd_create_ip_t.e_output_format != IV_YUV_422ILE) &&
               (ps_ip->s_ivd_create_ip_t.e_output_format != IV_RGB_565) &&
               (ps_ip->s_ivd_create_ip_t.e_output_format != IV_YUV_420SP_UV) &&
               (ps_ip->s_ivd_create_ip_t.e_output_format != IV_YUV_420SP_VU))
            {
                ps_op->s_ivd_create_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_create_op_t.u4_error_code |= IVD_INIT_DEC_COL_FMT_NOT_SUPPORTED;
                H264_DEC_DEBUG_PRINT("\n");
                return (IV_FAIL);
            }
        }
        break;

        case IVD_CMD_GET_DISPLAY_FRAME:
        {
            isvcd_get_display_frame_ip_t *ps_ip = (isvcd_get_display_frame_ip_t *) pv_api_ip;
            isvcd_get_display_frame_op_t *ps_op = (isvcd_get_display_frame_op_t *) pv_api_op;

            ps_op->s_ivd_get_display_frame_op_t.u4_error_code = 0;

            if((ps_ip->s_ivd_get_display_frame_ip_t.u4_size !=
                sizeof(isvcd_get_display_frame_ip_t)) &&
               (ps_ip->s_ivd_get_display_frame_ip_t.u4_size != sizeof(ivd_get_display_frame_ip_t)))
            {
                ps_op->s_ivd_get_display_frame_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_get_display_frame_op_t.u4_error_code |=
                    IVD_IP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }

            if((ps_op->s_ivd_get_display_frame_op_t.u4_size !=
                sizeof(isvcd_get_display_frame_op_t)) &&
               (ps_op->s_ivd_get_display_frame_op_t.u4_size != sizeof(ivd_get_display_frame_op_t)))
            {
                ps_op->s_ivd_get_display_frame_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_get_display_frame_op_t.u4_error_code |=
                    IVD_OP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }
        }
        break;

        case IVD_CMD_REL_DISPLAY_FRAME:
        {
            isvcd_rel_display_frame_ip_t *ps_ip = (isvcd_rel_display_frame_ip_t *) pv_api_ip;
            isvcd_rel_display_frame_op_t *ps_op = (isvcd_rel_display_frame_op_t *) pv_api_op;

            ps_op->s_ivd_rel_display_frame_op_t.u4_error_code = 0;

            if((ps_ip->s_ivd_rel_display_frame_ip_t.u4_size !=
                sizeof(isvcd_rel_display_frame_ip_t)) &&
               (ps_ip->s_ivd_rel_display_frame_ip_t.u4_size != sizeof(ivd_rel_display_frame_ip_t)))
            {
                ps_op->s_ivd_rel_display_frame_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_rel_display_frame_op_t.u4_error_code |=
                    IVD_IP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }

            if((ps_op->s_ivd_rel_display_frame_op_t.u4_size !=
                sizeof(isvcd_rel_display_frame_op_t)) &&
               (ps_op->s_ivd_rel_display_frame_op_t.u4_size != sizeof(ivd_rel_display_frame_op_t)))
            {
                ps_op->s_ivd_rel_display_frame_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_rel_display_frame_op_t.u4_error_code |=
                    IVD_OP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }
        }
        break;

        case IVD_CMD_SET_DISPLAY_FRAME:
        {
            isvcd_set_display_frame_ip_t *ps_ip = (isvcd_set_display_frame_ip_t *) pv_api_ip;
            isvcd_set_display_frame_op_t *ps_op = (isvcd_set_display_frame_op_t *) pv_api_op;
            UWORD32 j;

            ps_op->s_ivd_set_display_frame_op_t.u4_error_code = 0;

            if((ps_ip->s_ivd_set_display_frame_ip_t.u4_size !=
                sizeof(isvcd_set_display_frame_ip_t)) &&
               (ps_ip->s_ivd_set_display_frame_ip_t.u4_size != sizeof(ivd_set_display_frame_ip_t)))
            {
                ps_op->s_ivd_set_display_frame_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_set_display_frame_op_t.u4_error_code |=
                    IVD_IP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }

            if((ps_op->s_ivd_set_display_frame_op_t.u4_size !=
                sizeof(isvcd_set_display_frame_op_t)) &&
               (ps_op->s_ivd_set_display_frame_op_t.u4_size != sizeof(ivd_set_display_frame_op_t)))
            {
                ps_op->s_ivd_set_display_frame_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_set_display_frame_op_t.u4_error_code |=
                    IVD_OP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }

            if(ps_ip->s_ivd_set_display_frame_ip_t.num_disp_bufs == 0)
            {
                ps_op->s_ivd_set_display_frame_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_set_display_frame_op_t.u4_error_code |= IVD_DISP_FRM_ZERO_OP_BUFS;
                return IV_FAIL;
            }

            for(j = 0; j < ps_ip->s_ivd_set_display_frame_ip_t.num_disp_bufs; j++)
            {
                if(ps_ip->s_ivd_set_display_frame_ip_t.s_disp_buffer[j].u4_num_bufs == 0)
                {
                    ps_op->s_ivd_set_display_frame_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                    ps_op->s_ivd_set_display_frame_op_t.u4_error_code |= IVD_DISP_FRM_ZERO_OP_BUFS;
                    return IV_FAIL;
                }

                for(i = 0; i < ps_ip->s_ivd_set_display_frame_ip_t.s_disp_buffer[j].u4_num_bufs;
                    i++)
                {
                    if(ps_ip->s_ivd_set_display_frame_ip_t.s_disp_buffer[j].pu1_bufs[i] == NULL)
                    {
                        ps_op->s_ivd_set_display_frame_op_t.u4_error_code |=
                            1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_set_display_frame_op_t.u4_error_code |=
                            IVD_DISP_FRM_OP_BUF_NULL;
                        return IV_FAIL;
                    }

                    if(ps_ip->s_ivd_set_display_frame_ip_t.s_disp_buffer[j]
                           .u4_min_out_buf_size[i] == 0)
                    {
                        ps_op->s_ivd_set_display_frame_op_t.u4_error_code |=
                            1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_set_display_frame_op_t.u4_error_code |=
                            IVD_DISP_FRM_ZERO_OP_BUF_SIZE;
                        return IV_FAIL;
                    }
                }
            }
        }
        break;

        case IVD_CMD_VIDEO_DECODE:
        {
            isvcd_video_decode_ip_t *ps_ip = (isvcd_video_decode_ip_t *) pv_api_ip;
            isvcd_video_decode_op_t *ps_op = (isvcd_video_decode_op_t *) pv_api_op;

            H264_DEC_DEBUG_PRINT("The input bytes is: %d",
                                 ps_ip->s_ivd_video_decode_ip_t.u4_num_Bytes);
            ps_op->s_ivd_video_decode_op_t.u4_error_code = 0;

            if(ps_ip->s_ivd_video_decode_ip_t.u4_size != sizeof(isvcd_video_decode_ip_t) &&
               ps_ip->s_ivd_video_decode_ip_t.u4_size !=
                   offsetof(ivd_video_decode_ip_t, s_out_buffer))
            {
                ps_op->s_ivd_video_decode_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_video_decode_op_t.u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }

            if(ps_op->s_ivd_video_decode_op_t.u4_size != sizeof(isvcd_video_decode_op_t) &&
               ps_op->s_ivd_video_decode_op_t.u4_size !=
                   offsetof(ivd_video_decode_op_t, u4_output_present))
            {
                ps_op->s_ivd_video_decode_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_video_decode_op_t.u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }
            {
                svc_dec_ctxt_t *ps_svcd_ctxt;
                svc_dec_lyr_struct_t *ps_svc_lyr_dec;
                dec_struct_t *ps_dec;
                ps_svcd_ctxt = (svc_dec_ctxt_t *) (ps_handle->pv_codec_handle);
                ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[0];
                ps_dec = &ps_svc_lyr_dec->s_dec;
                if(ps_dec->u1_enable_mb_info)
                {
                    if(!ps_ip->pu1_8x8_blk_qp_map && !ps_ip->pu1_8x8_blk_type_map)
                    {
                        ps_op->s_ivd_video_decode_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_video_decode_op_t.u4_error_code |=
                            IH264D_FRAME_INFO_OP_BUF_NULL;
                        return IV_FAIL;
                    }
                }
            }
        }
        break;

        case IVD_CMD_DELETE:
        {
            isvcd_delete_ip_t *ps_ip = (isvcd_delete_ip_t *) pv_api_ip;
            isvcd_delete_op_t *ps_op = (isvcd_delete_op_t *) pv_api_op;

            ps_op->s_ivd_delete_op_t.u4_error_code = 0;

            if(ps_ip->s_ivd_delete_ip_t.u4_size != sizeof(isvcd_delete_ip_t))
            {
                ps_op->s_ivd_delete_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_delete_op_t.u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }

            if(ps_op->s_ivd_delete_op_t.u4_size != sizeof(isvcd_delete_op_t))
            {
                ps_op->s_ivd_delete_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_op->s_ivd_delete_op_t.u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                return (IV_FAIL);
            }
        }
        break;

        case IVD_CMD_VIDEO_CTL:
        {
            UWORD32 *pu4_ptr_cmd;
            UWORD32 sub_command;

            pu4_ptr_cmd = (UWORD32 *) pv_api_ip;
            pu4_ptr_cmd += 2;
            sub_command = *pu4_ptr_cmd;

            switch(sub_command)
            {
                case IVD_CMD_CTL_SETPARAMS:
                {
                    isvcd_ctl_set_config_ip_t *ps_ip;
                    isvcd_ctl_set_config_op_t *ps_op;
                    ps_ip = (isvcd_ctl_set_config_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_set_config_op_t *) pv_api_op;

                    if(ps_ip->s_ivd_ctl_set_config_ip_t.u4_size !=
                       sizeof(isvcd_ctl_set_config_ip_t))
                    {
                        ps_op->s_ivd_ctl_set_config_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_set_config_op_t.u4_error_code |=
                            IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                }
                break;

                case IVD_CMD_CTL_SETDEFAULT:
                {
                    isvcd_ctl_set_config_op_t *ps_op;
                    ps_op = (isvcd_ctl_set_config_op_t *) pv_api_op;
                    if(ps_op->s_ivd_ctl_set_config_op_t.u4_size !=
                       sizeof(isvcd_ctl_set_config_op_t))
                    {
                        ps_op->s_ivd_ctl_set_config_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_set_config_op_t.u4_error_code |=
                            IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                }
                break;

                case IVD_CMD_CTL_GETPARAMS:
                {
                    isvcd_ctl_getstatus_ip_t *ps_ip;
                    isvcd_ctl_getstatus_op_t *ps_op;

                    ps_ip = (isvcd_ctl_getstatus_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_getstatus_op_t *) pv_api_op;
                    if(ps_ip->s_ivd_ctl_getstatus_ip_t.u4_size != sizeof(isvcd_ctl_getstatus_ip_t))
                    {
                        ps_op->s_ivd_ctl_getstatus_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_getstatus_op_t.u4_error_code |=
                            IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                    if(ps_op->s_ivd_ctl_getstatus_op_t.u4_size != sizeof(isvcd_ctl_getstatus_op_t))
                    {
                        ps_op->s_ivd_ctl_getstatus_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_getstatus_op_t.u4_error_code |=
                            IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                }
                break;

                case IVD_CMD_CTL_GETBUFINFO:
                {
                    isvcd_ctl_getbufinfo_ip_t *ps_ip;
                    isvcd_ctl_getbufinfo_op_t *ps_op;
                    ps_ip = (isvcd_ctl_getbufinfo_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_getbufinfo_op_t *) pv_api_op;

                    if(ps_ip->s_ivd_ctl_getbufinfo_ip_t.u4_size !=
                       sizeof(isvcd_ctl_getbufinfo_ip_t))
                    {
                        ps_op->s_ivd_ctl_getbufinfo_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_getbufinfo_op_t.u4_error_code |=
                            IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                    if(ps_op->s_ivd_ctl_getbufinfo_op_t.u4_size !=
                       sizeof(isvcd_ctl_getbufinfo_op_t))
                    {
                        ps_op->s_ivd_ctl_getbufinfo_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_getbufinfo_op_t.u4_error_code |=
                            IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                }
                break;

                case IVD_CMD_CTL_GETVERSION:
                {
                    isvcd_ctl_getversioninfo_ip_t *ps_ip;
                    isvcd_ctl_getversioninfo_op_t *ps_op;
                    ps_ip = (isvcd_ctl_getversioninfo_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_getversioninfo_op_t *) pv_api_op;
                    if(ps_ip->s_ivd_ctl_getversioninfo_ip_t.u4_size !=
                       sizeof(isvcd_ctl_getversioninfo_ip_t))
                    {
                        ps_op->s_ivd_ctl_getversioninfo_op_t.u4_error_code |=
                            1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_getversioninfo_op_t.u4_error_code |=
                            IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                    if(ps_op->s_ivd_ctl_getversioninfo_op_t.u4_size !=
                       sizeof(isvcd_ctl_getversioninfo_op_t))
                    {
                        ps_op->s_ivd_ctl_getversioninfo_op_t.u4_error_code |=
                            1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_getversioninfo_op_t.u4_error_code |=
                            IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                }
                break;

                case IVD_CMD_CTL_FLUSH:
                {
                    isvcd_ctl_flush_ip_t *ps_ip;
                    isvcd_ctl_flush_op_t *ps_op;
                    ps_ip = (isvcd_ctl_flush_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_flush_op_t *) pv_api_op;
                    if(ps_ip->s_ivd_ctl_flush_ip_t.u4_size != sizeof(isvcd_ctl_flush_ip_t))
                    {
                        ps_op->s_ivd_ctl_flush_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_flush_op_t.u4_error_code |=
                            IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                    if(ps_op->s_ivd_ctl_flush_op_t.u4_size != sizeof(isvcd_ctl_flush_op_t))
                    {
                        ps_op->s_ivd_ctl_flush_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_flush_op_t.u4_error_code |=
                            IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                }
                break;

                case IVD_CMD_CTL_RESET:
                {
                    isvcd_ctl_reset_ip_t *ps_ip;
                    isvcd_ctl_reset_op_t *ps_op;
                    ps_ip = (isvcd_ctl_reset_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_reset_op_t *) pv_api_op;
                    if(ps_ip->s_ivd_ctl_reset_ip_t.u4_size != sizeof(isvcd_ctl_reset_ip_t))
                    {
                        ps_op->s_ivd_ctl_reset_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_reset_op_t.u4_error_code |=
                            IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                    if(ps_op->s_ivd_ctl_reset_op_t.u4_size != sizeof(isvcd_ctl_reset_op_t))
                    {
                        ps_op->s_ivd_ctl_reset_op_t.u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->s_ivd_ctl_reset_op_t.u4_error_code |=
                            IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }
                }
                break;

                case IH264D_CMD_CTL_DEGRADE:
                {
                    isvcd_ctl_degrade_ip_t *ps_ip;
                    isvcd_ctl_degrade_op_t *ps_op;

                    ps_ip = (isvcd_ctl_degrade_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_degrade_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_degrade_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_degrade_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if((ps_ip->i4_degrade_pics < 0) || (ps_ip->i4_degrade_pics > 4) ||
                       (ps_ip->i4_nondegrade_interval < 0) || (ps_ip->i4_degrade_type < 0) ||
                       (ps_ip->i4_degrade_type > 15))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        return IV_FAIL;
                    }

                    break;
                }

                case IH264D_CMD_CTL_GET_BUFFER_DIMENSIONS:
                {
                    isvcd_ctl_get_frame_dimensions_ip_t *ps_ip;
                    isvcd_ctl_get_frame_dimensions_op_t *ps_op;

                    ps_ip = (isvcd_ctl_get_frame_dimensions_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_get_frame_dimensions_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_get_frame_dimensions_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_get_frame_dimensions_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    break;
                }
                case IH264D_CMD_CTL_GET_VUI_PARAMS:
                {
                    isvcd_ctl_get_vui_params_ip_t *ps_ip;
                    isvcd_ctl_get_vui_params_op_t *ps_op;

                    ps_ip = (isvcd_ctl_get_vui_params_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_get_vui_params_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_get_vui_params_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_get_vui_params_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    break;
                }
                case IH264D_CMD_CTL_GET_SEI_MDCV_PARAMS:
                {
                    isvcd_ctl_get_sei_mdcv_params_ip_t *ps_ip;
                    isvcd_ctl_get_sei_mdcv_params_op_t *ps_op;

                    ps_ip = (isvcd_ctl_get_sei_mdcv_params_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_get_sei_mdcv_params_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_get_sei_mdcv_params_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_get_sei_mdcv_params_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    break;
                }

                case IH264D_CMD_CTL_GET_SEI_CLL_PARAMS:
                {
                    isvcd_ctl_get_sei_cll_params_ip_t *ps_ip;
                    isvcd_ctl_get_sei_cll_params_op_t *ps_op;

                    ps_ip = (isvcd_ctl_get_sei_cll_params_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_get_sei_cll_params_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_get_sei_cll_params_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_get_sei_cll_params_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    break;
                }

                case IH264D_CMD_CTL_GET_SEI_AVE_PARAMS:
                {
                    isvcd_ctl_get_sei_ave_params_ip_t *ps_ip;
                    isvcd_ctl_get_sei_ave_params_op_t *ps_op;

                    ps_ip = (isvcd_ctl_get_sei_ave_params_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_get_sei_ave_params_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_get_sei_ave_params_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_get_sei_ave_params_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    break;
                }

                case IH264D_CMD_CTL_GET_SEI_CCV_PARAMS:
                {
                    isvcd_ctl_get_sei_ccv_params_ip_t *ps_ip;
                    isvcd_ctl_get_sei_ccv_params_op_t *ps_op;

                    ps_ip = (isvcd_ctl_get_sei_ccv_params_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_get_sei_ccv_params_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_get_sei_ccv_params_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_get_sei_ccv_params_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    break;
                }

                case IH264D_CMD_CTL_SET_NUM_CORES:
                {
                    isvcd_ctl_set_num_cores_ip_t *ps_ip;
                    isvcd_ctl_set_num_cores_op_t *ps_op;

                    ps_ip = (isvcd_ctl_set_num_cores_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_set_num_cores_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_set_num_cores_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_set_num_cores_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if((ps_ip->u4_num_cores != 1) && (ps_ip->u4_num_cores != 2) &&
                       (ps_ip->u4_num_cores != 3) && (ps_ip->u4_num_cores != 4))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        return IV_FAIL;
                    }
                    break;
                }
                case IH264D_CMD_CTL_SET_PROCESSOR:
                {
                    isvcd_ctl_set_processor_ip_t *ps_ip;
                    isvcd_ctl_set_processor_op_t *ps_op;

                    ps_ip = (isvcd_ctl_set_processor_ip_t *) pv_api_ip;
                    ps_op = (isvcd_ctl_set_processor_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_ctl_set_processor_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_ctl_set_processor_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    break;
                }

                case ISVCD_CMD_CTL_SET_TGT_LAYER:
                {
                    isvcd_set_target_layer_ip_t *ps_ip;
                    isvcd_set_target_layer_op_t *ps_op;

                    ps_ip = (isvcd_set_target_layer_ip_t *) pv_api_ip;
                    ps_op = (isvcd_set_target_layer_op_t *) pv_api_op;

                    if(ps_ip->u4_size != sizeof(isvcd_set_target_layer_ip_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_IP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    if(ps_ip->u1_tgt_dep_id > MAX_DEPENDENCY_ID)
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        return IV_FAIL;
                    }

                    if(ps_ip->u1_tgt_temp_id > MAX_TEMPORAL_ID)
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        return IV_FAIL;
                    }

                    if(ps_ip->u1_tgt_quality_id > MAX_QUALITY_ID)
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        return IV_FAIL;
                    }

                    if(ps_ip->u1_tgt_priority_id > MAX_PRIORITY_ID)
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        return IV_FAIL;
                    }

                    if(ps_op->u4_size != sizeof(isvcd_set_target_layer_op_t))
                    {
                        ps_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                        ps_op->u4_error_code |= IVD_OP_API_STRUCT_SIZE_INCORRECT;
                        return IV_FAIL;
                    }

                    break;
                }

                default:
                    *(pu4_api_op + 1) |= 1 << IVD_UNSUPPORTEDPARAM;
                    *(pu4_api_op + 1) |= IVD_UNSUPPORTED_API_CMD;
                    return IV_FAIL;
                    break;
            }
        }
        break;
    }

    return IV_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief
 *  Sets Processor type
 *
 * @par Description:
 *  Sets Processor type
 *
 * @param[in] ps_codec_obj
 *  Pointer to codec object at API level
 *
 * @param[in] pv_api_ip
 *  Pointer to input argument structure
 *
 * @param[out] pv_api_op
 *  Pointer to output argument structure
 *
 * @returns  Status
 *
 * @remarks
 *
 *
 *******************************************************************************
 */

WORD32 isvcd_set_processor(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_ctl_set_processor_ip_t *ps_ip;
    isvcd_ctl_set_processor_op_t *ps_op;
    UWORD8 u1_layer_id;
    svc_dec_lyr_struct_t *ps_codec;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ip = (isvcd_ctl_set_processor_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_set_processor_op_t *) pv_api_op;

    ps_svcd_ctxt->e_processor_arch = (IVD_ARCH_T) ps_ip->u4_arch;
    ps_svcd_ctxt->e_processor_soc = (IVD_SOC_T) ps_ip->u4_soc;

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_codec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_codec->s_dec.e_processor_arch = (IVD_ARCH_T) ps_ip->u4_arch;
        ps_codec->s_dec.e_processor_soc = (IVD_SOC_T) ps_ip->u4_soc;

        isvcd_init_function_ptr(ps_codec);
    }

    ps_op->u4_error_code = 0;
    return IV_SUCCESS;
}

/**************************************************************************
 * \if Function name : isvcd_init_decoder \endif
 *
 *
 * \brief
 *    Initializes the decoder
 *
 * \param apiVersion               : Version of the api being used.
 * \param errorHandlingMechanism   : Mechanism to be used for errror handling.
 * \param postFilteringType: Type of post filtering operation to be used.
 * \param uc_outputFormat: Format of the decoded picture [default 4:2:0].
 * \param uc_dispBufs: Number of Display Buffers.
 * \param p_NALBufAPI: Pointer to NAL Buffer API.
 * \param p_DispBufAPI: Pointer to Display Buffer API.
 * \param ih264d_dec_mem_manager  :Pointer to the function that will be called
 *by decoder for memory allocation and freeing.
 *
 * \return
 *    0 on Success and -1 on error
 *
 **************************************************************************
 */
void isvcd_init_decoder(svc_dec_lyr_struct_t *ps_dec_svc_lyr_params)
{
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) ps_dec_svc_lyr_params;
    dec_struct_t *ps_dec;
    dec_slice_params_t *ps_cur_slice;
    pocstruct_t *ps_prev_poc, *ps_cur_poc;
    size_t size;
    ps_dec = &ps_svc_lyr_dec->s_dec;

    size = sizeof(pred_info_t) * 2 * 32;
    memset(ps_dec->ps_pred, 0, size);

    size = sizeof(disp_mgr_t);
    memset(ps_dec->pv_disp_buf_mgr, 0, size);

    size = ih264_buf_mgr_size();
    memset(ps_dec->pv_pic_buf_mgr, 0, size);

    size = sizeof(dec_err_status_t);
    memset(ps_dec->ps_dec_err_status, 0, size);

    size = sizeof(sei);
    memset(ps_dec->ps_sei, 0, size);

    size = sizeof(sei);
    memset(ps_dec->ps_sei_parse, 0, size);

    size = sizeof(dpb_commands_t);
    memset(ps_dec->ps_dpb_cmds, 0, size);

    size = sizeof(dec_bit_stream_t);
    memset(ps_dec->ps_bitstrm, 0, size);

    size = sizeof(dec_slice_params_t);
    memset(ps_dec->ps_cur_slice, 0, size);

    size = MAX(sizeof(dec_seq_params_t), sizeof(dec_pic_params_t));
    memset(ps_dec->pv_scratch_sps_pps, 0, size);

    size = sizeof(dec_svc_seq_params_t);
    memset(ps_svc_lyr_dec->pv_scratch_subset_sps, 0, size);

    size = sizeof(ctxt_inc_mb_info_t);
    memset(ps_dec->ps_left_mb_ctxt_info, 0, size);

    size = (sizeof(neighbouradd_t) << 2);
    memset(ps_dec->ps_left_mvpred_addr, 0, size);

    size = ih264_buf_mgr_size();
    memset(ps_dec->pv_mv_buf_mgr, 0, size);

    /* Free any dynamic buffers that are allocated */
    isvcd_free_dynamic_bufs(ps_svc_lyr_dec);

    {
        UWORD8 i;
        struct pic_buffer_t *ps_init_dpb;
        ps_init_dpb = ps_dec->ps_dpb_mgr->ps_init_dpb[0][0];
        for(i = 0; i < 2 * MAX_REF_BUFS; i++)
        {
            ps_init_dpb->pu1_buf1 = NULL;
            ps_init_dpb->u1_long_term_frm_idx = MAX_REF_BUFS + 1;
            ps_dec->ps_dpb_mgr->ps_init_dpb[0][i] = ps_init_dpb;
            ps_dec->ps_dpb_mgr->ps_mod_dpb[0][i] = ps_init_dpb;
            ps_init_dpb++;
        }

        ps_init_dpb = ps_dec->ps_dpb_mgr->ps_init_dpb[1][0];
        for(i = 0; i < 2 * MAX_REF_BUFS; i++)
        {
            ps_init_dpb->pu1_buf1 = NULL;
            ps_init_dpb->u1_long_term_frm_idx = MAX_REF_BUFS + 1;
            ps_dec->ps_dpb_mgr->ps_init_dpb[1][i] = ps_init_dpb;
            ps_dec->ps_dpb_mgr->ps_mod_dpb[1][i] = ps_init_dpb;
            ps_init_dpb++;
        }
    }

    ps_cur_slice = ps_dec->ps_cur_slice;
    ps_dec->init_done = 0;

    ps_dec->u4_num_cores = 1;
    ps_dec->u2_pic_ht = ps_dec->u2_pic_wd = 0;

    ps_dec->u1_separate_parse = DEFAULT_SEPARATE_PARSE;
    ps_dec->u4_app_disable_deblk_frm = 0;
    ps_dec->i4_degrade_type = 0;
    ps_dec->i4_degrade_pics = 0;

    /* Initialization of function pointers ih264d_deblock_picture function*/
    ps_dec->p_DeblockPicture[0] = ih264d_deblock_picture_non_mbaff;
    ps_dec->p_DeblockPicture[1] = ih264d_deblock_picture_mbaff;

    ps_dec->s_cab_dec_env.pv_codec_handle = ps_dec;
    ps_dec->u4_num_fld_in_frm = 0;
    ps_dec->ps_dpb_mgr->pv_codec_handle = ps_dec;

    /* Initialize the sei validity u4_flag with zero indiacting sei is not valid*/
    ps_dec->ps_sei->u1_is_valid = 0;

    /* decParams Initializations */
    ps_dec->ps_cur_pps = NULL;
    ps_dec->ps_cur_sps = NULL;
    ps_dec->u1_init_dec_flag = 0;
    ps_dec->u1_first_slice_in_stream = 1;
    ps_dec->u1_last_pic_not_decoded = 0;
    ps_dec->u4_app_disp_width = 0;
    ps_dec->i4_header_decoded = 0;
    ps_dec->u4_total_frames_decoded = 0;

    ps_dec->i4_error_code = 0;
    ps_dec->i4_content_type = IV_CONTENTTYPE_NA;
    ps_dec->ps_cur_slice->u1_mbaff_frame_flag = 0;

    ps_dec->ps_dec_err_status->u1_err_flag = ACCEPT_ALL_PICS;
    ps_dec->ps_dec_err_status->u1_cur_pic_type = PIC_TYPE_UNKNOWN;
    ps_dec->ps_dec_err_status->u4_frm_sei_sync = SYNC_FRM_DEFAULT;
    ps_dec->ps_dec_err_status->u4_cur_frm = INIT_FRAME;
    ps_dec->ps_dec_err_status->u1_pic_aud_i = PIC_TYPE_UNKNOWN;

    ps_dec->u1_pr_sl_type = 0xFF;
    ps_dec->u2_mbx = 0xffff;
    ps_dec->u2_mby = 0;
    ps_dec->u2_total_mbs_coded = 0;

    /* POC initializations */
    ps_prev_poc = &ps_dec->s_prev_pic_poc;
    ps_cur_poc = &ps_dec->s_cur_pic_poc;
    ps_prev_poc->i4_pic_order_cnt_lsb = ps_cur_poc->i4_pic_order_cnt_lsb = 0;
    ps_prev_poc->i4_pic_order_cnt_msb = ps_cur_poc->i4_pic_order_cnt_msb = 0;
    ps_prev_poc->i4_delta_pic_order_cnt_bottom = ps_cur_poc->i4_delta_pic_order_cnt_bottom = 0;
    ps_prev_poc->i4_delta_pic_order_cnt[0] = ps_cur_poc->i4_delta_pic_order_cnt[0] = 0;
    ps_prev_poc->i4_delta_pic_order_cnt[1] = ps_cur_poc->i4_delta_pic_order_cnt[1] = 0;
    ps_prev_poc->u1_mmco_equalto5 = ps_cur_poc->u1_mmco_equalto5 = 0;
    ps_prev_poc->i4_top_field_order_count = ps_cur_poc->i4_top_field_order_count = 0;
    ps_prev_poc->i4_bottom_field_order_count = ps_cur_poc->i4_bottom_field_order_count = 0;
    ps_prev_poc->u1_bot_field = ps_cur_poc->u1_bot_field = 0;
    ps_prev_poc->u1_mmco_equalto5 = ps_cur_poc->u1_mmco_equalto5 = 0;
    ps_prev_poc->i4_prev_frame_num_ofst = ps_cur_poc->i4_prev_frame_num_ofst = 0;
    ps_cur_slice->u1_mmco_equalto5 = 0;
    ps_cur_slice->u2_frame_num = 0;

    ps_dec->i4_max_poc = 0;
    ps_dec->i4_prev_max_display_seq = 0;
    ps_dec->u1_recon_mb_grp = 4;
    ps_dec->i4_reorder_depth = -1;

    /* Field PIC initializations */
    ps_dec->u1_second_field = 0;
    ps_dec->s_prev_seq_params.u1_eoseq_pending = 0;

    /* Set the cropping parameters as zero */
    ps_dec->u2_crop_offset_y = 0;
    ps_dec->u2_crop_offset_uv = 0;

    /* The Initial Frame Rate Info is not Present */
    ps_dec->i4_vui_frame_rate = -1;
    ps_dec->i4_pic_type = NA_SLICE;
    ps_dec->i4_frametype = IV_NA_FRAME;
    ps_dec->i4_content_type = IV_CONTENTTYPE_NA;

    ps_dec->u1_res_changed = 0;

    ps_dec->u1_frame_decoded_flag = 0;

    /* Set the default frame seek mask mode */
    ps_dec->u4_skip_frm_mask = SKIP_NONE;

    /********************************************************/
    /* Initialize CAVLC residual decoding function pointers */
    /********************************************************/
    ps_dec->pf_cavlc_4x4res_block[0] = ih264d_cavlc_4x4res_block_totalcoeff_1;
    ps_dec->pf_cavlc_4x4res_block[1] = ih264d_cavlc_4x4res_block_totalcoeff_2to10;
    ps_dec->pf_cavlc_4x4res_block[2] = ih264d_cavlc_4x4res_block_totalcoeff_11to16;

    ps_dec->pf_cavlc_parse4x4coeff[0] = ih264d_cavlc_parse4x4coeff_n0to7;
    ps_dec->pf_cavlc_parse4x4coeff[1] = ih264d_cavlc_parse4x4coeff_n8;

    ps_dec->pf_cavlc_parse_8x8block[0] = ih264d_cavlc_parse_8x8block_none_available;
    ps_dec->pf_cavlc_parse_8x8block[1] = ih264d_cavlc_parse_8x8block_left_available;
    ps_dec->pf_cavlc_parse_8x8block[2] = ih264d_cavlc_parse_8x8block_top_available;
    ps_dec->pf_cavlc_parse_8x8block[3] = ih264d_cavlc_parse_8x8block_both_available;

    /***************************************************************************/
    /* Initialize Bs calculation function pointers for P and B, 16x16/non16x16 */
    /***************************************************************************/
    ps_dec->pf_fill_bs1[0][0] = ih264d_fill_bs1_16x16mb_pslice;
    ps_dec->pf_fill_bs1[0][1] = ih264d_fill_bs1_non16x16mb_pslice;

    ps_dec->pf_fill_bs1[1][0] = ih264d_fill_bs1_16x16mb_bslice;
    ps_dec->pf_fill_bs1[1][1] = ih264d_fill_bs1_non16x16mb_bslice;

    ps_dec->pf_fill_bs_xtra_left_edge[0] = ih264d_fill_bs_xtra_left_edge_cur_frm;
    ps_dec->pf_fill_bs_xtra_left_edge[1] = ih264d_fill_bs_xtra_left_edge_cur_fld;

    /* Initialize Reference Pic Buffers */
    ih264d_init_ref_bufs(ps_dec->ps_dpb_mgr);

    ps_dec->u2_prv_frame_num = 0;
    ps_dec->u1_top_bottom_decoded = 0;
    ps_dec->u1_dangling_field = 0;

    ps_dec->s_cab_dec_env.cabac_table = gau4_ih264d_cabac_table;

    ps_dec->pu1_left_mv_ctxt_inc = ps_dec->u1_left_mv_ctxt_inc_arr[0];
    ps_dec->pi1_left_ref_idx_ctxt_inc = &ps_dec->i1_left_ref_idx_ctx_inc_arr[0][0];
    ps_dec->pu1_left_yuv_dc_csbp = &ps_dec->u1_yuv_dc_csbp_topmb;

    /* ! */
    /* Initializing flush frame u4_flag */
    ps_dec->u1_flushfrm = 0;

    ps_dec->s_cab_dec_env.pv_codec_handle = (void *) ps_dec;
    ps_dec->ps_bitstrm->pv_codec_handle = (void *) ps_dec;
    ps_dec->ps_cur_slice->pv_codec_handle = (void *) ps_dec;
    ps_dec->ps_dpb_mgr->pv_codec_handle = (void *) ps_dec;

    memset(ps_dec->disp_bufs, 0, (MAX_DISP_BUFS_NEW) * sizeof(disp_buf_t));
    memset(ps_dec->u4_disp_buf_mapping, 0, (MAX_DISP_BUFS_NEW) * sizeof(UWORD32));
    memset(ps_dec->u4_disp_buf_to_be_freed, 0, (MAX_DISP_BUFS_NEW) * sizeof(UWORD32));
    memset(ps_dec->ps_cur_slice, 0, sizeof(dec_slice_params_t));

    ih264d_init_arch(ps_dec);
    isvcd_init_function_ptr(ps_svc_lyr_dec);
    ps_dec->e_frm_out_mode = IVD_DISPLAY_FRAME_OUT;
    ps_dec->init_done = 1;
    ps_svc_lyr_dec->u1_layer_identifier = BASE_LAYER;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_parse_ctxt_free                                */
/*                                                                           */
/*  Description   :this function is used to free the nal parse context       */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

void isvcd_nal_parse_ctxt_free(svc_dec_ctxt_t *ps_svcd_ctxt)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    void (*pf_aligned_free)(void *pv_mem_ctxt, void *pv_buf);
    void *pv_mem_ctxt;
    nal_parse_ctxt_t *ps_ctxt;
    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[0];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    pf_aligned_free = ps_dec->pf_aligned_free;

    pv_mem_ctxt = ps_dec->pv_mem_ctxt;
    ps_ctxt = (nal_parse_ctxt_t *) ps_svcd_ctxt->pv_nal_parse_ctxt;

    pf_aligned_free(pv_mem_ctxt, ps_ctxt->s_dqid_ctxt.ps_dqid_node);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->pv_nal_header_buf);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->pv_nal_unit);
    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->pv_vcl_nal_buff);
    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->pv_non_vcl_nal_buff);
    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->pv_nal_parse_ctxt);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_resample_ctxt_free                        */
/*                                                                           */
/*  Description   :this function is used to free the resd_resamp context     */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

void isvcd_residual_resample_ctxt_free(svc_dec_ctxt_t *ps_svcd_ctxt)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    void (*pf_aligned_free)(void *pv_mem_ctxt, void *pv_buf);
    void *pv_mem_ctxt;
    residual_sampling_ctxt_t *ps_ctxt;
    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[0];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    pf_aligned_free = ps_dec->pf_aligned_free;

    pv_mem_ctxt = ps_dec->pv_mem_ctxt;
    ps_ctxt = (residual_sampling_ctxt_t *) ps_svcd_ctxt->pv_residual_sample_ctxt;

    pf_aligned_free(pv_mem_ctxt, ps_ctxt->pi2_refarray_buffer);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->pu1_ref_x_ptr_incr);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_x_offset_length);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_y_offset_length);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_x_pos_phase);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_y_pos_phase);
    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->pv_residual_sample_ctxt);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_intra_resample_ctxt_free                           */
/*                                                                           */
/*  Description   :this function is used to free the intra_resamp context    */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/
void isvcd_intra_resample_ctxt_free(svc_dec_ctxt_t *ps_svcd_ctxt)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    void (*pf_aligned_free)(void *pv_mem_ctxt, void *pv_buf);
    void *pv_mem_ctxt;
    intra_sampling_ctxt_t *ps_ctxt;
    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[0];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    pf_aligned_free = ps_dec->pf_aligned_free;

    pv_mem_ctxt = ps_dec->pv_mem_ctxt;
    ps_ctxt = (intra_sampling_ctxt_t *) ps_svcd_ctxt->pv_intra_sample_ctxt;

    pf_aligned_free(pv_mem_ctxt, ps_ctxt->pu1_refarray_buffer);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->pu1_refarray_cb);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->pu1_refarray_cr);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->pi4_temp_interpolation_buffer);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_x_offset_length);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_y_offset_length);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_x_min_max);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_y_min_max);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_x_pos_phase);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_y_pos_phase);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.pi2_xd_index);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.pi2_yd_index);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.pi2_ya_index);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_seg_lookup_horz);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.ps_seg_lookup_vert);
    pf_aligned_free(pv_mem_ctxt, ps_ctxt->as_res_lyrs[0].s_luma_map_ctxt.pu1_refarray_x_idx);
    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->pv_intra_sample_ctxt);
    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->pv_ii_pred_ctxt);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_mode_mv_resample_ctxt_free                         */
/*                                                                           */
/*  Description   :this function is used to free the mv resamp context       */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/
void isvcd_mode_mv_resample_ctxt_free(svc_dec_ctxt_t *ps_svcd_ctxt)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    void (*pf_aligned_free)(void *pv_mem_ctxt, void *pv_buf);
    void *pv_mem_ctxt;
    mode_motion_ctxt_t *ps_mode_motion;

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[0];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    pf_aligned_free = ps_dec->pf_aligned_free;

    pv_mem_ctxt = ps_dec->pv_mem_ctxt;
    ps_mode_motion = (mode_motion_ctxt_t *) ps_svcd_ctxt->pv_mode_mv_sample_ctxt;

    pf_aligned_free(pv_mem_ctxt, ps_mode_motion->ps_motion_pred_struct);
    pf_aligned_free(pv_mem_ctxt, ps_mode_motion->as_res_lyr_mem[0].pi2_ref_loc_x);
    pf_aligned_free(pv_mem_ctxt, ps_mode_motion->as_res_lyr_mem[0].pi2_ref_loc_y);
    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->pv_ref_lyr_offset);
    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->pv_mode_mv_sample_ctxt);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_free_static_bufs                                   */
/*                                                                           */
/*  Description   :this function is used to free the static buffers          */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_free_static_bufs(iv_obj_t *dec_hdl)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;

    UWORD8 u1_layer_id;
    svc_dec_ctxt_t *ps_svcd_ctxt;

    void (*pf_aligned_free)(void *pv_mem_ctxt, void *pv_buf);
    void *pv_mem_ctxt;

    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    isvcd_intra_resample_ctxt_free(ps_svcd_ctxt);
    isvcd_residual_resample_ctxt_free(ps_svcd_ctxt);
    isvcd_mode_mv_resample_ctxt_free(ps_svcd_ctxt);
    isvcd_nal_parse_ctxt_free(ps_svcd_ctxt);

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_dec = &ps_svc_lyr_dec->s_dec;
        pf_aligned_free = ps_dec->pf_aligned_free;
        pv_mem_ctxt = ps_dec->pv_mem_ctxt;

#ifdef KEEP_THREADS_ACTIVE
        /* Wait for threads */
        ps_dec->i4_break_threads = 1;
        if(ps_dec->u4_dec_thread_created)
        {
            ithread_mutex_lock(ps_dec->apv_proc_start_mutex[0]);

            ps_dec->ai4_process_start[0] = PROC_START;

            ithread_cond_signal(ps_dec->apv_proc_start_condition[0]);

            ithread_mutex_unlock(ps_dec->apv_proc_start_mutex[0]);

            ithread_join(ps_dec->pv_dec_thread_handle, NULL);

            ps_dec->u4_dec_thread_created = 0;
        }

        if(ps_dec->u4_bs_deblk_thread_created)
        {
            ithread_mutex_lock(ps_dec->apv_proc_start_mutex[1]);

            ps_dec->ai4_process_start[1] = PROC_START;

            ithread_cond_signal(ps_dec->apv_proc_start_condition[1]);

            ithread_mutex_unlock(ps_dec->apv_proc_start_mutex[1]);

            ithread_join(ps_dec->pv_bs_deblk_thread_handle, NULL);

            ps_dec->u4_bs_deblk_thread_created = 0;
        }

        // destroy mutex and condition variable for both the threads
        // 1. ih264d_decode_picture_thread
        // 2. ih264d_recon_deblk_thread
        {
            UWORD32 i;
            for(i = 0; i < 2; i++)
            {
                ithread_cond_destroy(ps_dec->apv_proc_start_condition[i]);
                ithread_cond_destroy(ps_dec->apv_proc_done_condition[i]);

                ithread_mutex_destroy(ps_dec->apv_proc_start_mutex[i]);
                ithread_mutex_destroy(ps_dec->apv_proc_done_mutex[i]);
            }
        }
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->apv_proc_start_mutex[0]);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->apv_proc_start_condition[0]);
#endif
        if(0 == u1_layer_id)
        {
            UWORD8 u1_sps_ctr;
            PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_sps);
            for(u1_sps_ctr = 0; u1_sps_ctr < (2 * MAX_NUM_SEQ_PARAMS); u1_sps_ctr++)
            {
                if(NULL != ps_svcd_ctxt->ps_subset_sps[u1_sps_ctr].s_sps_svc_ext.ps_svc_vui_ext)
                {
                    PS_DEC_ALIGNED_FREE(
                        ps_dec,
                        ps_svcd_ctxt->ps_subset_sps[u1_sps_ctr].s_sps_svc_ext.ps_svc_vui_ext);
                }
            }
            PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->ps_subset_sps);
            PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_pps);
            PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_sei);
            PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_sei_parse);
        }

        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pv_dec_thread_handle);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pv_bs_deblk_thread_handle);

        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_dpb_mgr);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_pred);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pv_disp_buf_mgr);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_pic_buf_base);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_dec_err_status);

        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_dpb_cmds);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_bitstrm);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->ps_nal_svc_ext);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_cur_slice);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pv_scratch_sps_pps);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->pv_scratch_subset_sps);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pu1_bits_buf_static);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ppv_map_ref_idx_to_poc_base);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->p_cabac_ctxt_table_t);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_left_mb_ctxt_info);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pu1_ref_buff_base);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pi2_pred1);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pu1_temp_mc_buffer);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pu1_init_dpb_base);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pu4_mbaff_wt_mat);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pu4_wts_ofsts_mat);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_left_mvpred_addr);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->ps_col_mv_base);
        PS_DEC_ALIGNED_FREE(ps_dec, ps_svc_lyr_dec->pu1_ii_resamp_buffer_luma);

        if(NULL != ps_dec->pv_pic_buf_mgr)
        {
            if(u1_layer_id < ps_svcd_ctxt->u1_prev_num_res_layers)
            {
                if(((buf_mgr_t *) ps_dec->pv_pic_buf_mgr)->pv_mutex != NULL)
                    ih264_buf_mgr_free(ps_dec->pv_pic_buf_mgr);
            }
            PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pv_pic_buf_mgr);
        }
        if(NULL != ps_dec->pv_mv_buf_mgr)
        {
            if(u1_layer_id < ps_svcd_ctxt->u1_prev_num_res_layers)
            {
                if(((buf_mgr_t *) ps_dec->pv_mv_buf_mgr)->pv_mutex != NULL)
                    ih264_buf_mgr_free(ps_dec->pv_mv_buf_mgr);
            }
            PS_DEC_ALIGNED_FREE(ps_dec, ps_dec->pv_mv_buf_mgr);
        }
    }

    pf_aligned_free(pv_mem_ctxt, ps_svcd_ctxt->ps_svc_dec_lyr);
    pf_aligned_free(pv_mem_ctxt, dec_hdl->pv_codec_handle);

    if(dec_hdl)
    {
        pf_aligned_free(pv_mem_ctxt, dec_hdl);
    }

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_parse_init                                      */
/*                                                                           */
/*  Description   : Initiaization of allocation of memory                    */
/*                                                                           */
/*  Inputs        : pv_mem_rec - Allocated memory records                    */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : None                                                     */
/*                                                                           */
/*  Returns       : Module's handle                                          */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay           Draft                               */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_nal_parse_ctxt_create(svc_dec_ctxt_t *ps_svcd_ctxt, void *pv_api_ip, void *pv_api_op)
{
    isvcd_create_ip_t *ps_create_ip;
    void *pv_buf;
    void *(*pf_aligned_alloc)(void *pv_mem_ctxt, WORD32 alignment, WORD32 size);
    void *pv_mem_ctxt;
    WORD32 size;
    nal_parse_ctxt_t *ps_nal_parse_ctxt;
    UWORD8 *pu1_ptr;
    UNUSED(pv_api_op);

    ps_create_ip = (isvcd_create_ip_t *) pv_api_ip;

    pf_aligned_alloc = ps_create_ip->s_ivd_create_ip_t.pf_aligned_alloc;
    pv_mem_ctxt = ps_create_ip->s_ivd_create_ip_t.pv_mem_ctxt;

    /*-----------------------------------------------------------------------*/
    /* Handle                                                                */
    /*-----------------------------------------------------------------------*/
    size = sizeof(nal_parse_ctxt_t);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_nal_parse_ctxt = pv_buf;

    /* set the lowest dqid to -1 */
    ps_nal_parse_ctxt->i4_prev_dq_id = -1;

    /*-----------------------------------------------------------------------*/
    /* DQID list buffer and initialization of vcl node buffer context        */
    /*-----------------------------------------------------------------------*/
    {
        WORD32 i4_lyr_idx;
        WORD32 i4_max_num_lyrs;
        vcl_node_t *ps_vcl_node;
        dqid_node_t *ps_dqid_node;
        dqid_ctxt_t *ps_dqid_ctxt;

        size = sizeof(vcl_node_t);
        size += sizeof(dqid_node_t);
        size *= MAX_NUM_RES_LYRS;

        ps_dqid_ctxt = &ps_nal_parse_ctxt->s_dqid_ctxt;

        ps_dqid_ctxt->i4_max_num_lyrs = MAX_NUM_RES_LYRS;

        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);

        ps_dqid_ctxt->ps_dqid_node = pv_buf;
        ps_dqid_node = ps_dqid_ctxt->ps_dqid_node;

        i4_max_num_lyrs = ps_dqid_ctxt->i4_max_num_lyrs;

        pu1_ptr = pv_buf;
        pu1_ptr += sizeof(dqid_node_t) * i4_max_num_lyrs;
        ps_vcl_node = (vcl_node_t *) pu1_ptr;

        for(i4_lyr_idx = 0; i4_lyr_idx < i4_max_num_lyrs; i4_lyr_idx++)
        {
            ps_dqid_node->ps_vcl_node = ps_vcl_node;

            /* Loop updates */
            ps_vcl_node += 1;
            ps_dqid_node += 1;
        } /* Loop over all the layers */
    }

    /*-----------------------------------------------------------------------*/
    /* Common memory                                                         */
    /*-----------------------------------------------------------------------*/
    size = UP_ALIGN_8(HEADER_BUFFER_LEN_BEFORE_EP);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_nal_parse_ctxt->pv_nal_header_buf = (void *) pv_buf;

    /*-----------------------------------------------------------------------*/
    /* Layer params memory                                                   */
    /*-----------------------------------------------------------------------*/
    size = sizeof(nal_unit_t);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_nal_parse_ctxt->pv_nal_unit = pv_buf;

    size = MAX_VCL_NAL_BUFF_SIZE * sizeof(UWORD8);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->pv_vcl_nal_buff = pv_buf;

    size = MAX_NON_VCL_NAL_BUFF_SIZE * sizeof(UWORD8);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->pv_non_vcl_nal_buff = pv_buf;

    /*-----------------------------------------------------------------------*/
    /* Registering the seq and pic prms buffer pointers                      */
    /*-----------------------------------------------------------------------*/
    if(NULL == ps_svcd_ctxt->ps_sps || NULL == ps_svcd_ctxt->ps_pps)
    {
        return IV_FAIL;
    }

    ps_svcd_ctxt->pv_nal_parse_ctxt = ps_nal_parse_ctxt;
    ps_nal_parse_ctxt->pv_seq_prms = ps_svcd_ctxt->ps_sps;
    ps_nal_parse_ctxt->pv_pic_prms = ps_svcd_ctxt->ps_pps;

    /* register VCL and NON VCL buffer pointers */
    if(NULL == ps_svcd_ctxt->pv_vcl_nal_buff || NULL == ps_svcd_ctxt->pv_non_vcl_nal_buff)
    {
        return IV_FAIL;
    }

    ps_nal_parse_ctxt->pv_vcl_nal_buf = (UWORD8 *) ps_svcd_ctxt->pv_vcl_nal_buff;
    ps_nal_parse_ctxt->pv_non_vcl_nal_buf = (UWORD8 *) ps_svcd_ctxt->pv_non_vcl_nal_buff;
    isvcd_nal_parse_reset_ctxt(ANNEX_B, PARTIAL_INPUT_MODE, ps_nal_parse_ctxt);

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_intra_resample_ctxt_create                         */
/*                                                                           */
/*  Description   :this function is used to create intra_resamp context      */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_intra_resample_ctxt_create(svc_dec_ctxt_t *ps_svcd_ctxt, void *pv_api_ip,
                                        void *pv_api_op)
{
    isvcd_create_ip_t *ps_create_ip;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    void *pv_buf;
    UWORD8 u1_layer_id;
    void *(*pf_aligned_alloc)(void *pv_mem_ctxt, WORD32 alignment, WORD32 size);
    void *pv_mem_ctxt;
    WORD32 size;
    intra_inter_pred_ctxt_t *ps_ii_pred_ctxt;

    intra_sampling_ctxt_t *ps_ctxt;
    UNUSED(pv_api_op);
    ps_create_ip = (isvcd_create_ip_t *) pv_api_ip;

    pf_aligned_alloc = ps_create_ip->s_ivd_create_ip_t.pf_aligned_alloc;
    pv_mem_ctxt = ps_create_ip->s_ivd_create_ip_t.pv_mem_ctxt;

    {
        intra_samp_lyr_ctxt *ps_lyr_ctxt;

        /* allocate context structure */
        size = ((sizeof(intra_sampling_ctxt_t) + 127) >> 7) << 7;
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_ctxt = pv_buf;

        /* luma reference array buffer  */
        size = REF_ARRAY_WIDTH * REF_ARRAY_HEIGHT * sizeof(UWORD8);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_ctxt->pu1_refarray_buffer = pv_buf;

        /* cb reference array buffer */
        size = REF_ARRAY_WIDTH * REF_ARRAY_HEIGHT * sizeof(UWORD8);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_ctxt->pu1_refarray_cb = pv_buf;

        /* cr reference array buffer */
        size = ((DYADIC_REF_W_C + 2) * (DYADIC_REF_H_C + 2) * sizeof(UWORD8));
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_ctxt->pu1_refarray_cr = pv_buf;

        /* Temp Intermediate Buffer */
        size = INTERMEDIATE_BUFF_WIDTH * INTERMEDIATE_BUFF_HEIGHT * sizeof(WORD32);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_ctxt->pi4_temp_interpolation_buffer = pv_buf;

        /****************** projected locations buffers ******************/
        {
            intra_samp_map_ctxt_t *ps_luma_map;
            intra_samp_map_ctxt_t *ps_chroma_map;
            WORD32 i4_lyr_id;
            ref_mb_map_t *ps_off_len_map;
            ref_pixel_map_t *ps_pos_phase_map;
            ref_min_max_map_t *ps_min_max;
            WORD16 *pi2_mem;
            UWORD8 *pu1_mem;
            seg_lookup_desc_t *ps_seg_lookup;

            /****************** Horz offset length ******************/

            size = (H264_MAX_FRAME_WIDTH >> 4) * MAX_NUM_RES_LYRS * 2 * sizeof(ref_mb_map_t);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            ps_off_len_map = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->ps_x_offset_length = ps_off_len_map;
                ps_off_len_map += (H264_MAX_FRAME_WIDTH >> 4);
                ps_chroma_map->ps_x_offset_length = ps_off_len_map;
                ps_off_len_map += (H264_MAX_FRAME_WIDTH >> 4);

            } /* end of loop over resolution layers */

            /****************** Vert offset length ******************/
            size = (H264_MAX_FRAME_HEIGHT >> 4) * MAX_NUM_RES_LYRS * 2 * sizeof(ref_mb_map_t);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            ps_off_len_map = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->ps_y_offset_length = ps_off_len_map;
                ps_off_len_map += (H264_MAX_FRAME_HEIGHT >> 4);
                ps_chroma_map->ps_y_offset_length = ps_off_len_map;
                ps_off_len_map += (H264_MAX_FRAME_HEIGHT >> 4);

            } /* end of loop over resolution layers */

            /****************** Horz Min Max Pos ******************/

            size = (H264_MAX_FRAME_WIDTH >> 4) * MAX_NUM_RES_LYRS * 2 * sizeof(ref_mb_map_t);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            ps_min_max = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->ps_x_min_max = ps_min_max;
                ps_min_max += (H264_MAX_FRAME_WIDTH >> 4);
                ps_chroma_map->ps_x_min_max = ps_min_max;
                ps_min_max += (H264_MAX_FRAME_WIDTH >> 4);
            } /* end of loop over resolution layers */

            /****************** Vert Min Max Pos ******************/
            size = (H264_MAX_FRAME_HEIGHT >> 4) * MAX_NUM_RES_LYRS * 2 * sizeof(ref_mb_map_t);

            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            ps_min_max = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->ps_y_min_max = ps_min_max;
                ps_min_max += (H264_MAX_FRAME_HEIGHT >> 4);
                ps_chroma_map->ps_y_min_max = ps_min_max;
                ps_min_max += (H264_MAX_FRAME_HEIGHT >> 4);

            } /* end of loop over resolution layers */

            /****************** Horz position phase ******************/
            size = (H264_MAX_FRAME_WIDTH) *MAX_NUM_RES_LYRS * 2 * sizeof(ref_pixel_map_t);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            ps_pos_phase_map = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->ps_x_pos_phase = ps_pos_phase_map;
                ps_pos_phase_map += (H264_MAX_FRAME_WIDTH);
                ps_chroma_map->ps_x_pos_phase = ps_pos_phase_map;
                ps_pos_phase_map += (H264_MAX_FRAME_WIDTH);

            } /* end of loop over resolution layers */

            /****************** Vert position phase ******************/

            size = (H264_MAX_FRAME_HEIGHT) *MAX_NUM_RES_LYRS * 2 * sizeof(ref_pixel_map_t);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            ps_pos_phase_map = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->ps_y_pos_phase = ps_pos_phase_map;
                ps_pos_phase_map += (H264_MAX_FRAME_HEIGHT);
                ps_chroma_map->ps_y_pos_phase = ps_pos_phase_map;
                ps_pos_phase_map += (H264_MAX_FRAME_HEIGHT);

            } /* end of loop over resolution layers */

            /**************** XD Index ******************************/
            size = (MB_WIDTH) *MAX_NUM_RES_LYRS * 2 * sizeof(WORD16);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            pi2_mem = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->pi2_xd_index = pi2_mem;
                pi2_mem += MB_WIDTH;
                ps_chroma_map->pi2_xd_index = pi2_mem;
                pi2_mem += MB_WIDTH;

            } /* end of loop over resolution layers */

            /**************** YD Index ******************************/
            size = (MB_HEIGHT) *MAX_NUM_RES_LYRS * 2 * sizeof(WORD16);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            pi2_mem = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->pi2_yd_index = pi2_mem;
                pi2_mem += MB_HEIGHT;
                ps_chroma_map->pi2_yd_index = pi2_mem;
                pi2_mem += MB_HEIGHT;

            } /* end of loop over resolution layers */

            /**************** YA Index ******************************/
            size = MB_HEIGHT * MAX_NUM_RES_LYRS * 2 * sizeof(WORD16);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            pi2_mem = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->pi2_ya_index = pi2_mem;
                pi2_mem += MB_HEIGHT;
                ps_chroma_map->pi2_ya_index = pi2_mem;
                pi2_mem += MB_HEIGHT;

            } /* end of loop over resolution layers */

            /**************** Horizontal segment lookup **************************/
            /* (MB_WIDTH x seg_lookup_desc_t) x (num layers - 1)   (for luma   )*/
            /* (BLOCK_WIDTH x seg_lookup_desc_t) x (num layers - 1) (for chroma )*/
            size = (MB_WIDTH * sizeof(seg_lookup_desc_t)) * MAX_NUM_RES_LYRS;

            size += (BLOCK_WIDTH * sizeof(seg_lookup_desc_t)) * MAX_NUM_RES_LYRS;

            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            ps_seg_lookup = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->ps_seg_lookup_horz = ps_seg_lookup;
                ps_seg_lookup += MB_WIDTH;
                ps_chroma_map->ps_seg_lookup_horz = ps_seg_lookup;
                ps_seg_lookup += BLOCK_WIDTH;

            } /* end of loop over resolution layers */

            /**************** Vertical segment lookup ****************************/
            /* (MB_HEIGHT x seg_lookup_desc_t) x (num layers - 1)    (for luma  )*/
            /* (BLOCK_HEIGHT x seg_lookup_desc_t) x (num layers - 1) (for chroma)*/
            size = (MB_HEIGHT * sizeof(seg_lookup_desc_t)) * MAX_NUM_RES_LYRS;

            size += (BLOCK_HEIGHT * sizeof(seg_lookup_desc_t)) * MAX_NUM_RES_LYRS;

            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            ps_seg_lookup = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->ps_seg_lookup_vert = ps_seg_lookup;
                ps_seg_lookup += MB_HEIGHT;
                ps_chroma_map->ps_seg_lookup_vert = ps_seg_lookup;
                ps_seg_lookup += BLOCK_HEIGHT;

            } /* end of loop over resolution layers */

            /**************** X and Y Reference Array Index lookup ***************/
            /* (MAX_REF_IDX_ARRAY) x (num layers - 1)     (for luma  x-index)     */
            /* (MAX_REF_IDX_ARRAY) x (num layers - 1)     (for luma  y-index)     */
            /* (MAX_REF_IDX_ARRAY) x (num layers - 1)     (for chroma x-index)    */
            /* (MAX_REF_IDX_ARRAY) x (num layers - 1)     (for chroma y-index)    */
            /*********************************************************************/
            size = (MAX_REF_IDX_ARRAY * MAX_NUM_RES_LYRS * 4);
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);
            pu1_mem = pv_buf;

            /* loop over num layers -1 */
            for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
            {
                /* derive the layer map ctxt */
                ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
                ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
                ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

                /* initialise the pointers */
                ps_luma_map->pu1_refarray_x_idx = pu1_mem;
                pu1_mem += MAX_REF_IDX_ARRAY;

                ps_luma_map->pu1_refarray_y_idx = pu1_mem;
                pu1_mem += MAX_REF_IDX_ARRAY;

                ps_chroma_map->pu1_refarray_x_idx = pu1_mem;
                pu1_mem += MAX_REF_IDX_ARRAY;

                ps_chroma_map->pu1_refarray_y_idx = pu1_mem;
                pu1_mem += MAX_REF_IDX_ARRAY;

            } /* end of loop over resolution layers */
        }

        size = ((sizeof(intra_inter_pred_ctxt_t) + 127) >> 7) << 7;
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_ii_pred_ctxt = pv_buf;
    }

    ps_svcd_ctxt->pv_intra_sample_ctxt = ps_ctxt;
    ps_svcd_ctxt->pv_ii_pred_ctxt = ps_ii_pred_ctxt;

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_svc_lyr_dec->pv_intra_sample_ctxt = ps_svcd_ctxt->pv_intra_sample_ctxt;
        ps_svc_lyr_dec->pv_ii_pred_ctxt = ps_svcd_ctxt->pv_ii_pred_ctxt;
    }

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_resample_ctxt_create                      */
/*                                                                           */
/*  Description   :this function is used to create resd_resamp context       */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   ittiam                creation                       */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_residual_resample_ctxt_create(svc_dec_ctxt_t *ps_svcd_ctxt, void *pv_api_ip,
                                           void *pv_api_op)
{
    isvcd_create_ip_t *ps_create_ip;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    void *pv_buf;
    UWORD8 u1_layer_id;
    void *(*pf_aligned_alloc)(void *pv_mem_ctxt, WORD32 alignment, WORD32 size);
    void *pv_mem_ctxt;
    WORD32 size;

    residual_sampling_ctxt_t *ps_ctxt;
    res_lyr_ctxt *ps_lyr_ctxt;
    UNUSED(pv_api_op);
    ps_create_ip = (isvcd_create_ip_t *) pv_api_ip;

    pf_aligned_alloc = ps_create_ip->s_ivd_create_ip_t.pf_aligned_alloc;
    pv_mem_ctxt = ps_create_ip->s_ivd_create_ip_t.pv_mem_ctxt;

    /* allocate context structure */
    size = ((sizeof(residual_sampling_ctxt_t) + 127) >> 7) << 7;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_ctxt = pv_buf;

    /* reference array buffer  */
    size = REF_ARRAY_WIDTH_RES_SAMP * REF_ARRAY_HEIGHT_RES_SAMP * sizeof(WORD16);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_ctxt->pi2_refarray_buffer = pv_buf;

    /* reference array pointer increment buffer */
    {
        WORD32 i4_size;

        i4_size = REF_ARRAY_WIDTH_RES_SAMP * REF_ARRAY_HEIGHT_RES_SAMP * sizeof(UWORD8);
        size = REF_ARRAY_WIDTH_RES_SAMP * REF_ARRAY_HEIGHT_RES_SAMP * 2 * sizeof(UWORD8);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_ctxt->pu1_ref_x_ptr_incr = pv_buf;
        ps_ctxt->pu1_ref_y_ptr_incr = ps_ctxt->pu1_ref_x_ptr_incr + i4_size;
    }

    /****************** projected locations buffers ******************/
    {
        residual_samp_map_ctxt_t *ps_luma_map;
        residual_samp_map_ctxt_t *ps_chroma_map;
        WORD32 i4_lyr_id;
        ref_mb_map_t *ps_off_len_map;
        ref_pixel_map_t *ps_pos_phase_map;

        /****************** Horz offset length ******************/
        size = (H264_MAX_FRAME_WIDTH >> 4) * MAX_NUM_RES_LYRS * 2 * sizeof(ref_mb_map_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_off_len_map = pv_buf;

        /* loop over num layers -1 */
        for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
        {
            /* derive the layer map ctxt */
            ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
            ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
            ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

            /* initialise the pointers */
            ps_luma_map->ps_x_offset_length = ps_off_len_map;
            ps_off_len_map += (H264_MAX_FRAME_WIDTH >> 4);
            ps_chroma_map->ps_x_offset_length = ps_off_len_map;
            ps_off_len_map += (H264_MAX_FRAME_WIDTH >> 4);

        } /* end of loop over resolution layers */

        /****************** Vert offset length ******************/
        size = (H264_MAX_FRAME_HEIGHT >> 4) * MAX_NUM_RES_LYRS * 2 * sizeof(ref_mb_map_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_off_len_map = pv_buf;

        /* loop over num layers -1 */
        for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
        {
            /* derive the layer map ctxt */
            ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
            ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
            ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

            /* initialise the pointers */
            ps_luma_map->ps_y_offset_length = ps_off_len_map;
            ps_off_len_map += (H264_MAX_FRAME_HEIGHT >> 4);
            ps_chroma_map->ps_y_offset_length = ps_off_len_map;
            ps_off_len_map += (H264_MAX_FRAME_HEIGHT >> 4);

        } /* end of loop over resolution layers */

        /****************** Horz position phase ******************/
        size = H264_MAX_FRAME_WIDTH * MAX_NUM_RES_LYRS * 2 * sizeof(ref_pixel_map_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_pos_phase_map = pv_buf;

        /* loop over num layers -1 */
        for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
        {
            /* derive the layer map ctxt */
            ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
            ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
            ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

            /* initialise the pointers */
            ps_luma_map->ps_x_pos_phase = ps_pos_phase_map;
            ps_pos_phase_map += H264_MAX_FRAME_WIDTH;
            ps_chroma_map->ps_x_pos_phase = ps_pos_phase_map;
            ps_pos_phase_map += H264_MAX_FRAME_WIDTH;

        } /* end of loop over resolution layers */

        /****************** Vert position phase ******************/

        size = H264_MAX_FRAME_HEIGHT * MAX_NUM_RES_LYRS * 2 * sizeof(ref_pixel_map_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_pos_phase_map = pv_buf;

        /* loop over num layers -1 */
        for(i4_lyr_id = 0; i4_lyr_id < MAX_NUM_RES_LYRS; i4_lyr_id++)
        {
            /* derive the layer map ctxt */
            ps_lyr_ctxt = &ps_ctxt->as_res_lyrs[i4_lyr_id];
            ps_luma_map = &ps_lyr_ctxt->s_luma_map_ctxt;
            ps_chroma_map = &ps_lyr_ctxt->s_chroma_map_ctxt;

            /* initialise the pointers */
            ps_luma_map->ps_y_pos_phase = ps_pos_phase_map;
            ps_pos_phase_map += H264_MAX_FRAME_HEIGHT;
            ps_chroma_map->ps_y_pos_phase = ps_pos_phase_map;
            ps_pos_phase_map += H264_MAX_FRAME_HEIGHT;

        } /* end of loop over resolution layers */
    }

    ps_svcd_ctxt->pv_residual_sample_ctxt = ps_ctxt;

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_svc_lyr_dec->pv_residual_sample_ctxt = ps_svcd_ctxt->pv_residual_sample_ctxt;
    }
    return IV_SUCCESS;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_mode_mv_resample_ctxt_create                       */
/*                                                                           */
/*  Description   :this function is used to create mv_resamp context         */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_mode_mv_resample_ctxt_create(svc_dec_ctxt_t *ps_svcd_ctxt, void *pv_api_ip,
                                          void *pv_api_op)
{
    isvcd_create_ip_t *ps_create_ip;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    void *pv_buf;
    WORD16 *pi2_mem;
    UWORD8 u1_layer_id;
    void *(*pf_aligned_alloc)(void *pv_mem_ctxt, WORD32 alignment, WORD32 size);
    void *pv_mem_ctxt;
    WORD32 size, i4_res_id;
    ref_lyr_scaled_offset_t *ps_ref_pic_offsets;
    mode_motion_ctxt_t *ps_mode_motion;
    mode_motion_lyr_ctxt *ps_lyr_mem;
    UNUSED(pv_api_op);
    ps_create_ip = (isvcd_create_ip_t *) pv_api_ip;

    pf_aligned_alloc = ps_create_ip->s_ivd_create_ip_t.pf_aligned_alloc;
    pv_mem_ctxt = ps_create_ip->s_ivd_create_ip_t.pv_mem_ctxt;

    size = ((sizeof(mode_motion_ctxt_t) + 127) >> 7) << 7;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_mode_motion = pv_buf;

    /* motion pred structure */
    size = 2 * NUM_MB_PARTS * NUM_SUB_MB_PARTS * sizeof(mv_pred_t);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_mode_motion->ps_motion_pred_struct = (mv_pred_t *) pv_buf;

    /* projected locations X */
    size = H264_MAX_FRAME_WIDTH * MAX_NUM_RES_LYRS * sizeof(WORD16);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    pi2_mem = (WORD16 *) pv_buf;

    /* loop over NUM resolution layers */
    for(i4_res_id = 0; i4_res_id < MAX_NUM_RES_LYRS; i4_res_id++)
    {
        ps_lyr_mem = &ps_mode_motion->as_res_lyr_mem[i4_res_id];

        /* initialise the pointers */
        ps_lyr_mem->pi2_ref_loc_x = pi2_mem;

        /* increment the buffer pointer */
        pi2_mem += H264_MAX_FRAME_WIDTH;

    } /* end of loop over num resolution layers */

    /* projected locations Y */
    size = H264_MAX_FRAME_HEIGHT * MAX_NUM_RES_LYRS * sizeof(WORD16);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    pi2_mem = (WORD16 *) pv_buf;
    /* loop over NUM resolution layers */
    for(i4_res_id = 0; i4_res_id < MAX_NUM_RES_LYRS; i4_res_id++)
    {
        ps_lyr_mem = &ps_mode_motion->as_res_lyr_mem[i4_res_id];

        /* initialise the pointers */
        ps_lyr_mem->pi2_ref_loc_y = pi2_mem;
        /* increment the buffer pointer */
        pi2_mem += H264_MAX_FRAME_HEIGHT;

    } /* end of loop over num resolution layers */

    size = sizeof(ref_lyr_scaled_offset_t) * MAX_NUM_RES_LYRS * MAX_NUM_PIC_BUFS;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->pv_ref_lyr_offset = pv_buf;

    /* loop over NUM resolution layers */
    ps_ref_pic_offsets = (ref_lyr_scaled_offset_t *) ps_svcd_ctxt->pv_ref_lyr_offset;

    for(i4_res_id = 0; i4_res_id < MAX_NUM_RES_LYRS; i4_res_id++)
    {
        ps_lyr_mem = &ps_mode_motion->as_res_lyr_mem[i4_res_id];

        /* store the current resolution layer pic offset start pointer */
        ps_lyr_mem->ps_ref_pic_lyr_offsets = ps_ref_pic_offsets + (i4_res_id * MAX_NUM_PIC_BUFS);

    } /* end of loop over num resolution layers */

    ps_svcd_ctxt->pv_mode_mv_sample_ctxt = ps_mode_motion;

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_svc_lyr_dec->pv_mode_mv_sample_ctxt = ps_svcd_ctxt->pv_mode_mv_sample_ctxt;
        ps_svc_lyr_dec->pv_ref_lyr_offset = ps_svcd_ctxt->pv_ref_lyr_offset;
    }
    return IV_SUCCESS;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_allocate_static_bufs                               */
/*                                                                           */
/*  Description   : allocates static buffers                                 */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_allocate_static_bufs(iv_obj_t **dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_create_ip_t *ps_create_ip;
    isvcd_create_op_t *ps_create_op;
    void *pv_buf;
    UWORD8 *pu1_buf;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    void *(*pf_aligned_alloc)(void *pv_mem_ctxt, WORD32 alignment, WORD32 size);
    void (*pf_aligned_free)(void *pv_mem_ctxt, void *pv_buf);
    void *pv_mem_ctxt;
    WORD32 size;
    UWORD8 u1_layer_id, u1_sps_ctr;
    UWORD8 u1_chroma_format;
    WORD32 ret;

    ps_create_ip = (isvcd_create_ip_t *) pv_api_ip;
    ps_create_op = (isvcd_create_op_t *) pv_api_op;

    ps_create_op->s_ivd_create_op_t.u4_error_code = 0;
    pf_aligned_alloc = ps_create_ip->s_ivd_create_ip_t.pf_aligned_alloc;
    pf_aligned_free = ps_create_ip->s_ivd_create_ip_t.pf_aligned_free;
    pv_mem_ctxt = ps_create_ip->s_ivd_create_ip_t.pv_mem_ctxt;
    u1_chroma_format = (UWORD8) (ps_create_ip->s_ivd_create_ip_t.e_output_format);

    if((u1_chroma_format != IV_YUV_420P) && (u1_chroma_format != IV_YUV_420SP_UV) &&
       (u1_chroma_format != IV_YUV_420SP_VU))
    {
        ps_create_op->s_ivd_create_op_t.pv_handle = NULL;

        return IV_FAIL;
    }

    /* Initialize return handle to NULL */
    ps_create_op->s_ivd_create_op_t.pv_handle = NULL;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, sizeof(iv_obj_t));
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, sizeof(iv_obj_t));
    *dec_hdl = (iv_obj_t *) pv_buf;
    ps_create_op->s_ivd_create_op_t.pv_handle = *dec_hdl;

    (*dec_hdl)->pv_codec_handle = NULL;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, sizeof(svc_dec_ctxt_t));
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    (*dec_hdl)->pv_codec_handle = (svc_dec_ctxt_t *) pv_buf;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) pv_buf;

    memset(ps_svcd_ctxt, 0, sizeof(svc_dec_ctxt_t));

    ps_svcd_ctxt->u1_prev_num_res_layers = UINT8_MAX;
    ps_svcd_ctxt->u1_pre_parse_in_flush = 1;
    /* set default to maximum values supported */
    ps_svcd_ctxt->u1_tgt_dep_id = MAX_DEPENDENCY_ID;
    ps_svcd_ctxt->u1_tgt_quality_id = MAX_QUALITY_ID;
    ps_svcd_ctxt->u1_tgt_temp_id = MAX_TEMPORAL_ID;
    ps_svcd_ctxt->u1_tgt_priority_id = MAX_PRIORITY_ID;

    /* two sets of MAX_NUM_SEQ_PARAMS are created one for sps-base layer;  one for
     * subset_sps- enhancement*/
    size = ((sizeof(dec_seq_params_t)) * MAX_NUM_SEQ_PARAMS * 2);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->ps_sps = pv_buf;

    /* two sets of MAX_NUM_SEQ_PARAMS are created one for sps-base layer;  one for
     * subset_sps- enhancement*/
    size = ((sizeof(dec_svc_seq_params_t)) * MAX_NUM_SEQ_PARAMS * 2);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->ps_subset_sps = pv_buf;

    for(u1_sps_ctr = 0; u1_sps_ctr < (2 * MAX_NUM_SEQ_PARAMS); u1_sps_ctr++)
    {
        ps_svcd_ctxt->ps_subset_sps[u1_sps_ctr].ps_seq = &ps_svcd_ctxt->ps_sps[u1_sps_ctr];
    }

    size = sizeof(sei);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->ps_sei = (sei *) pv_buf;

    size = sizeof(sei);
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->ps_sei_parse = (sei *) pv_buf;

    size = (sizeof(dec_pic_params_t)) * MAX_NUM_PIC_PARAMS;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->ps_pps = pv_buf;

    size = (sizeof(svc_dec_lyr_struct_t)) * MAX_NUM_RES_LYRS;
    pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
    RETURN_IF((NULL == pv_buf), IV_FAIL);
    memset(pv_buf, 0, size);
    ps_svcd_ctxt->ps_svc_dec_lyr = pv_buf;
    ps_svcd_ctxt->u1_target_layer_id = 0;
    ps_svcd_ctxt->u1_cur_layer_id = 0;
    ps_svcd_ctxt->i4_eos_flag = 0;

    ret = isvcd_mode_mv_resample_ctxt_create(ps_svcd_ctxt, pv_api_ip, pv_api_op);
    if(ret != IV_SUCCESS)
    {
        return ret;
    }
    ret = isvcd_intra_resample_ctxt_create(ps_svcd_ctxt, pv_api_ip, pv_api_op);
    if(ret != IV_SUCCESS)
    {
        return ret;
    }
    ret = isvcd_residual_resample_ctxt_create(ps_svcd_ctxt, pv_api_ip, pv_api_op);
    if(ret != IV_SUCCESS)
    {
        return ret;
    }
    ret = isvcd_nal_parse_ctxt_create(ps_svcd_ctxt, pv_api_ip, pv_api_op);
    if(ret != IV_SUCCESS)
    {
        return ret;
    }
    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_dec = &ps_svc_lyr_dec->s_dec;

        ps_svc_lyr_dec->ps_svcd_ctxt = ps_svcd_ctxt;
        ps_svc_lyr_dec->u1_layer_id = u1_layer_id;
        ps_svc_lyr_dec->u1_dyadic_flag = 1;
        ps_svc_lyr_dec->u1_restricted_res_change_flag = 1;
        ps_svc_lyr_dec->u1_base_res_flag = 1;
        ps_svc_lyr_dec->u1_ref_layer_id = 0;
        ps_svc_lyr_dec->ps_dec_svc_ref_layer =
            &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svc_lyr_dec->u1_ref_layer_id];
        ps_svc_lyr_dec->u4_pps_id_for_layer = UINT32_MAX;

#ifndef LOGO_EN
        ps_dec->u4_share_disp_buf = ps_create_ip->s_ivd_create_ip_t.u4_share_disp_buf;
#else
        ps_dec->u4_share_disp_buf = 0;
#endif

        ps_dec->u1_chroma_format = (UWORD8) (ps_create_ip->s_ivd_create_ip_t.e_output_format);

        if((ps_dec->u1_chroma_format != IV_YUV_420P) &&
           (ps_dec->u1_chroma_format != IV_YUV_420SP_UV) &&
           (ps_dec->u1_chroma_format != IV_YUV_420SP_VU))
        {
            ps_dec->u4_share_disp_buf = 0;
        }

        ps_dec->u1_enable_mb_info = ps_create_ip->u4_enable_frame_info;
        ps_dec->pf_aligned_alloc = pf_aligned_alloc;
        ps_dec->pf_aligned_free = pf_aligned_free;
        ps_dec->pv_mem_ctxt = pv_mem_ctxt;

        ps_dec->ps_sps = ps_svcd_ctxt->ps_sps;
        ps_svc_lyr_dec->ps_subset_sps = ps_svcd_ctxt->ps_subset_sps;
        ps_dec->ps_pps = ps_svcd_ctxt->ps_pps;
        ps_dec->ps_sei = ps_svcd_ctxt->ps_sei;
        ps_dec->ps_sei_parse = ps_svcd_ctxt->ps_sei_parse;

        size = ithread_get_handle_size();
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pv_dec_thread_handle = pv_buf;

        size = ithread_get_handle_size();
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pv_bs_deblk_thread_handle = pv_buf;

#ifdef KEEP_THREADS_ACTIVE
        {
            UWORD32 i;
            /* Request memory to hold mutex (start/done) for both threads */
            size = ithread_get_mutex_lock_size() << 2;
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 8, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);

            // init mutex variable for both the threads
            // 1. ih264d_decode_picture_thread
            // 2. ih264d_recon_deblk_thread
            for(i = 0; i < 2; i++)
            {
                WORD32 ret;
                WORD32 mutex_size = ithread_get_mutex_lock_size();

                ps_dec->apv_proc_start_mutex[i] = (UWORD8 *) pv_buf + (2 * i * mutex_size);
                ps_dec->apv_proc_done_mutex[i] = (UWORD8 *) pv_buf + ((2 * i + 1) * mutex_size);

                ret = ithread_mutex_init(ps_dec->apv_proc_start_mutex[0]);
                RETURN_IF((ret != IV_SUCCESS), ret);

                ret = ithread_mutex_init(ps_dec->apv_proc_done_mutex[i]);
                RETURN_IF((ret != IV_SUCCESS), ret);
            }

            size = ithread_get_cond_struct_size() << 2;
            pv_buf = pf_aligned_alloc(pv_mem_ctxt, 8, size);
            RETURN_IF((NULL == pv_buf), IV_FAIL);
            memset(pv_buf, 0, size);

            // init condition variable for both the threads
            for(i = 0; i < 2; i++)
            {
                WORD32 ret;
                WORD32 cond_size = ithread_get_cond_struct_size();
                ps_dec->apv_proc_start_condition[i] = (UWORD8 *) pv_buf + (2 * i * cond_size);
                ps_dec->apv_proc_done_condition[i] = (UWORD8 *) pv_buf + ((2 * i + 1) * cond_size);

                ret = ithread_cond_init(ps_dec->apv_proc_start_condition[i]);
                RETURN_IF((ret != IV_SUCCESS), ret);

                ret = ithread_cond_init(ps_dec->apv_proc_done_condition[i]);
                RETURN_IF((ret != IV_SUCCESS), ret);
            }
        }
#endif
        size = sizeof(dpb_manager_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_dpb_mgr = pv_buf;

        size = sizeof(pred_info_t) * 2 * 32;
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_pred = pv_buf;

        size = sizeof(disp_mgr_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pv_disp_buf_mgr = pv_buf;

        size = ih264_buf_mgr_size();
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pv_pic_buf_mgr = pv_buf;

        size = sizeof(struct pic_buffer_t) * (H264_MAX_REF_PICS * 2);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_pic_buf_base = pv_buf;

        size = sizeof(dec_err_status_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_dec_err_status = (dec_err_status_t *) pv_buf;

        size = sizeof(dpb_commands_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_dpb_cmds = (dpb_commands_t *) pv_buf;

        size = sizeof(dec_bit_stream_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_bitstrm = (dec_bit_stream_t *) pv_buf;

        size = sizeof(dec_nal_unit_svc_ext_params_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_svc_lyr_dec->ps_nal_svc_ext = (dec_nal_unit_svc_ext_params_t *) pv_buf;

        size = sizeof(dec_slice_params_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_cur_slice = (dec_slice_params_t *) pv_buf;

        size = MAX(sizeof(dec_seq_params_t), sizeof(dec_pic_params_t));
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pv_scratch_sps_pps = pv_buf;

        size = sizeof(dec_svc_seq_params_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_svc_lyr_dec->pv_scratch_subset_sps = pv_buf;

        ps_dec->u4_static_bits_buf_size = 256000;
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, ps_dec->u4_static_bits_buf_size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, ps_dec->u4_static_bits_buf_size);
        ps_dec->pu1_bits_buf_static = pv_buf;

        size = ((TOTAL_LIST_ENTRIES + PAD_MAP_IDX_POC) * sizeof(void *));
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        ps_dec->ppv_map_ref_idx_to_poc_base = pv_buf;
        memset(ps_dec->ppv_map_ref_idx_to_poc_base, 0, size);

        ps_dec->ppv_map_ref_idx_to_poc = ps_dec->ppv_map_ref_idx_to_poc_base + OFFSET_MAP_IDX_POC;

        size = (sizeof(bin_ctxt_model_t) * NUM_CABAC_CTXTS_SVC);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->p_cabac_ctxt_table_t = pv_buf;

        size = sizeof(ctxt_inc_mb_info_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_left_mb_ctxt_info = pv_buf;

        size = MAX_REF_BUF_SIZE * 2;
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pu1_ref_buff_base = pv_buf;
        ps_dec->pu1_ref_buff = ps_dec->pu1_ref_buff_base + MAX_REF_BUF_SIZE;

        size = ((sizeof(WORD16)) * PRED_BUFFER_WIDTH * PRED_BUFFER_HEIGHT * 2);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pi2_pred1 = pv_buf;

        size = sizeof(UWORD8) * (MB_LUM_SIZE);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pu1_temp_mc_buffer = pv_buf;

        size = 8 * MAX_REF_BUFS * sizeof(struct pic_buffer_t);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);

        ps_dec->pu1_init_dpb_base = pv_buf;
        pu1_buf = pv_buf;
        ps_dec->ps_dpb_mgr->ps_init_dpb[0][0] = (struct pic_buffer_t *) pu1_buf;

        pu1_buf += size / 2;
        ps_dec->ps_dpb_mgr->ps_init_dpb[1][0] = (struct pic_buffer_t *) pu1_buf;

        size = (sizeof(UWORD32) * 2 * 3 * ((MAX_FRAMES << 1) * (MAX_FRAMES << 1)) * 2);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pu4_mbaff_wt_mat = pv_buf;

        size = sizeof(UWORD32) * 2 * 3 * ((MAX_FRAMES << 1) * (MAX_FRAMES << 1));
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pu4_wts_ofsts_mat = pv_buf;

        size = (sizeof(neighbouradd_t) << 2);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->ps_left_mvpred_addr = pv_buf;

        size = ih264_buf_mgr_size();
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        memset(pv_buf, 0, size);
        ps_dec->pv_mv_buf_mgr = pv_buf;

        size = sizeof(col_mv_buf_t) * (H264_MAX_REF_PICS * 2);
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        ps_dec->ps_col_mv_base = pv_buf;
        memset(ps_dec->ps_col_mv_base, 0, size);

        size = ((MB_SIZE * MB_SIZE * 3) >> 1) + MB_SIZE;
        pv_buf = pf_aligned_alloc(pv_mem_ctxt, 128, size);
        RETURN_IF((NULL == pv_buf), IV_FAIL);
        ps_svc_lyr_dec->pu1_ii_resamp_buffer_luma = pv_buf;
        ps_svc_lyr_dec->pu1_ii_resamp_buffer_chroma =
            ps_svc_lyr_dec->pu1_ii_resamp_buffer_luma + (MB_SIZE * MB_SIZE);
        memset(ps_svc_lyr_dec->pu1_ii_resamp_buffer_luma, 0, size);

        isvcd_init_decoder(ps_svc_lyr_dec);
    }
    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_create                                             */
/*                                                                           */
/*  Description   : creates decoder                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore                                              */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_create(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_create_ip_t *ps_create_ip;
    isvcd_create_op_t *ps_create_op;
    WORD32 ret;

    ps_create_ip = (isvcd_create_ip_t *) pv_api_ip;
    ps_create_op = (isvcd_create_op_t *) pv_api_op;

    ps_create_op->s_ivd_create_op_t.u4_error_code = 0;
    dec_hdl = NULL;
    ret = isvcd_allocate_static_bufs(&dec_hdl, pv_api_ip, pv_api_op);

    /* If allocation of some buffer fails, then free buffers allocated till then */
    if(IV_FAIL == ret)
    {
        if(dec_hdl)
        {
            if(dec_hdl->pv_codec_handle)
            {
                isvcd_free_static_bufs(dec_hdl);
            }
            else
            {
                void (*pf_aligned_free)(void *pv_mem_ctxt, void *pv_buf);
                void *pv_mem_ctxt;

                pf_aligned_free = ps_create_ip->s_ivd_create_ip_t.pf_aligned_free;
                pv_mem_ctxt = ps_create_ip->s_ivd_create_ip_t.pv_mem_ctxt;
                pf_aligned_free(pv_mem_ctxt, dec_hdl);
            }
        }
        ps_create_op->s_ivd_create_op_t.u4_error_code = IVD_MEM_ALLOC_FAILED;
        ps_create_op->s_ivd_create_op_t.u4_error_code |= 1 << IVD_FATALERROR;
        return IV_FAIL;
    }

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_update_dqid                                        */
/*                                                                           */
/*  Description   : Updates the DQID list based on reference layer DQID      */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. Reference layer DQID                                  */
/*                  2. current layer's vcl node structure                    */
/*                  3. pointer to store the bottom layer VCL node            */
/*  Globals       : None                                                     */
/*  Processing    : 1. Searches for a layer with reference layer DQID        */
/*                  2. Updates the bottom and top nodes of current layer and */
/*                     reference layer vcl nodes respectively                */
/*                                                                           */
/*  Outputs       : Updates top and bottom node field of vcl nodes of current*/
/*                  layer and reference layer respectively                   */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*****************************************************************************/
WORD32 isvcd_update_dqid(WORD32 i4_ref_lyr_dqid, vcl_node_t *ps_cur_lyr_node,
                         vcl_node_t **pps_bot_lyr_node)
{
    vcl_node_t *ps_vcl_node;

    /* sanity checks */
    if((NULL == ps_cur_lyr_node) || (NULL == pps_bot_lyr_node))
    {
        return NOT_OK;
    }

    ps_vcl_node = ps_cur_lyr_node->ps_bot_node;
    while(NULL != ps_vcl_node)
    {
        WORD32 i4_dqid;

        i4_dqid = (ps_vcl_node->i4_dependency_id << 4) + ps_vcl_node->i4_quality_id;

        /* if reference layer DQ ID matches */
        /* or reference layer is a layer below reference dq id layer */
        if((i4_dqid == i4_ref_lyr_dqid) ||
           (ps_vcl_node->i4_quality_id < (i4_ref_lyr_dqid & 0x0F)) ||
           (ps_vcl_node->i4_dependency_id < (i4_ref_lyr_dqid >> 4)))
        {
            break;
        }
        ps_vcl_node = ps_vcl_node->ps_bot_node;
    }

    /* Update the top and bottom node of ref layer and current layer nodes */

    if(NULL != ps_vcl_node)
    {
        ps_cur_lyr_node->ps_bot_node = ps_vcl_node;
        ps_vcl_node->ps_top_node = ps_cur_lyr_node;
    }

    /* Update pointer to bottom VCL node */
    *pps_bot_lyr_node = ps_vcl_node;
    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_detect_res_change                                  */
/*                                                                           */
/*  Description   : This function detects the resolution change              */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. Pointer to Current SPS                                */
/*                  2. Pointer to prevoius SPS                               */
/*  Globals       : None                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : SVCD_TRUE if different resolution else SVCD_FALSE        */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijayakumar      Draft                               */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_detect_res_change(dec_seq_params_t *ps_curr_sps, dec_seq_params_t *ps_prev_sps,
                               dec_svc_seq_params_t *ps_curr_subset_sps,
                               dec_svc_seq_params_t *ps_prev_subset_sps)
{
    UWORD16 u2_scaled_ref_width_sps;
    UWORD16 u2_scaled_ref_ht_sps;
    UNUSED(ps_prev_subset_sps);

    if(NULL == ps_prev_sps)
    {
        /* indicates bottom most layer in Access unit */
        return (SVCD_FALSE);
    }
    /* Check for the ESS idc */
    if(2 == ps_curr_subset_sps->s_sps_svc_ext.u1_extended_spatial_scalability_idc)
    {
        return (SVCD_TRUE);
    }

    /* Calculate the scaled reference width and height */
    u2_scaled_ref_width_sps = (ps_curr_sps->u2_frm_wd_in_mbs << 4);
    u2_scaled_ref_width_sps -=
        (ps_curr_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_left_offset +
         ps_curr_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_right_offset);

    u2_scaled_ref_ht_sps = (ps_curr_sps->u2_frm_ht_in_mbs << 4);
    u2_scaled_ref_ht_sps -=
        (ps_curr_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_top_offset +
         ps_curr_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_bottom_offset);

    /* Check for frame width being different */
    if(u2_scaled_ref_width_sps != (ps_prev_sps->u2_frm_wd_in_mbs << 4))
    {
        return (SVCD_TRUE);
    }

    /* Check for frame height being different */
    if(u2_scaled_ref_ht_sps != (ps_prev_sps->u2_frm_ht_in_mbs << 4))
    {
        return (SVCD_TRUE);
    }

    /* check for crop offset not MB aligned */
    if((0 != (ps_curr_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_left_offset & 15)) ||
       (0 != (ps_curr_subset_sps->s_sps_svc_ext.i4_seq_scaled_ref_layer_top_offset & 15)))
    {
        return (SVCD_TRUE);
    }

    /* check for chroma Phase Y being different */
    if(ps_curr_subset_sps->s_sps_svc_ext.u1_chroma_phase_x_plus1_flag !=
       ps_curr_subset_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_x_plus1_flag)
    {
        return (SVCD_TRUE);
    }

    /* check for chroma Phase Y being different */
    if(ps_curr_subset_sps->s_sps_svc_ext.u1_chroma_phase_y_plus1 !=
       ps_curr_subset_sps->s_sps_svc_ext.u1_seq_ref_layer_chroma_phase_y_plus1)
    {
        return (SVCD_TRUE);
    }

    /* If none of the above are true then there is no resolution change */
    return (SVCD_FALSE);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_ref_pic_list_modify                          */
/*                                                                           */
/*  Description   : Parses the reference picture modification related        */
/*                  syntax elements                                          */
/*                                                                           */
/*  Inputs        : 1. stream context structure                              */
/*                  2. slice prms structure                                  */
/*  Globals       : None                                                     */
/*  Processing    : Parses the syntax elements                               */
/*                                                                           */
/*  Outputs       : Updated stream buffer context                            */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_parse_ref_pic_list_modify(dec_bit_stream_t *ps_bitstrm,
                                       dec_slice_params_t *ps_slice_prms,
                                       dec_seq_params_t *ps_curr_sps)
{
    WORD32 i4_mod_flag;
    UWORD16 ui_nextUev;
    WORD32 i4_num_sets_ctr = 0;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;

    if(I_SLICE != ps_slice_prms->u1_slice_type)
    {
        /* ref_pic_list_modification_flag_l0 */
        i4_mod_flag = ih264d_get_bit_h264(ps_bitstrm);

        if(0 != i4_mod_flag)
        {
            WORD32 i4_mod_pic_num_idc;

            i4_num_sets_ctr = 0;
            do
            {
                /* modification_of_pic_nums_idc */
                i4_mod_pic_num_idc = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

                if((i4_mod_pic_num_idc > 3) || (i4_mod_pic_num_idc < 0))
                {
                    return ERROR_INV_SLICE_HDR_T;
                }
                if(3 != i4_mod_pic_num_idc)
                {
                    /* i4_mod_pic_num_idc = 0,1 ==> abs_diff_pic_num_minus1 */
                    /* i4_mod_pic_num_idc = 2 ==> long_term_pic_num */

                    ui_nextUev = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                    if(ui_nextUev > (ps_curr_sps->u2_u4_max_pic_num_minus1 + 1))
                        return ERROR_INV_SLICE_HDR_T;
                }

                i4_num_sets_ctr++;

                /* if the number of commands recieved exceeds max limit */
                if((H264_MAX_REF_PICS) == i4_num_sets_ctr) break;

            } while(3 != i4_mod_pic_num_idc);
        }

        /*********** if (I_SLICE != u1_slice_type) ***************************/
    }

    if(B_SLICE != ps_slice_prms->u1_slice_type)
    {
        return (OK);
    }

    /* ref_pic_list_modification_flag_l1 */
    i4_mod_flag = ih264d_get_bit_h264(ps_bitstrm);

    if(0 != i4_mod_flag)
    {
        WORD32 i4_mod_pic_num_idc;

        do
        {
            /* modification_of_pic_nums_idc */
            i4_mod_pic_num_idc = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

            if((i4_mod_pic_num_idc > 3) || (i4_mod_pic_num_idc < 0))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            if(3 != i4_mod_pic_num_idc)
            {
                /* i4_mod_pic_num_idc = 0,1 ==> abs_diff_pic_num_minus1 */
                /* i4_mod_pic_num_idc = 2 ==> long_term_pic_num */

                ui_nextUev = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                if(ui_nextUev > (ps_curr_sps->u2_u4_max_pic_num_minus1 + 1))
                    return ERROR_INV_SLICE_HDR_T;
            }

            i4_num_sets_ctr++;

            /* if the number of commands recieved exceeds max limit */
            if((H264_MAX_REF_PICS) == i4_num_sets_ctr) break;

        } while(3 != i4_mod_pic_num_idc);
    }

    return (OK);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_slice_hdr_refdq_id                           */
/*                                                                           */
/*  Description   : this function decodes the slice until Inter layer deblk  */
/*                  parameters are parsed                                    */
/*                                                                           */
/*  Inputs        : 1. pointer to VCL node of given layer                    */
/*                  2. pointer toi slice params structure                    */
/*                  3. pointer to layer params strcuture                     */
/*                  4. pointer to stream context structure                   */
/*                  5. pointer to store the reference DQID                   */
/*                  6. pointer to array of SPS                               */
/*                  7. pointer to array of PPS                               */
/*                  8. no inter layer pred flag of current slice             */
/*  Globals       : none                                                     */
/*  Processing    : prases syntax elements                                   */
/*                                                                           */
/*  Outputs       : reference layer DQ ID                                    */
/*  Returns       : Error code                                               */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_parse_slice_hdr_refdq_id(vcl_node_t *ps_vcl_node, dec_slice_params_t *ps_slice_prms,
                                      dec_bit_stream_t *ps_bitstrm, WORD32 *pi4_ref_dq_id,
                                      dec_seq_params_t *ps_sps, dec_pic_params_t *ps_pps,
                                      WORD32 i4_no_int_lyr_pred, svc_dec_ctxt_t *ps_svcd_ctxt)
{
    UWORD8 u1_pps_id;
    WORD32 i_temp;
    UWORD32 u4_temp;
    WORD32 i4_nal_unit_type;
    WORD32 i4_nal_ref_idc, i4_quality_id;
    WORD32 i4_use_ref_base, i4_idr_pic_flag;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    dec_svc_seq_params_t *ps_subset_sps = NULL;
    WORD32 ret = OK;

    i4_nal_unit_type = ps_vcl_node->i4_nal_unit_type;
    i4_nal_ref_idc = ps_vcl_node->i4_nal_ref_idc;
    i4_quality_id = ps_vcl_node->i4_quality_id;
    i4_use_ref_base = ps_vcl_node->i4_use_ref_base;
    i4_idr_pic_flag = ps_vcl_node->i4_idr_pic_flag;

    /*-----------------------------------------------------------------------*/
    /*--------------------- first mb in slice -------------------------------*/
    /*-----------------------------------------------------------------------*/
    ps_slice_prms->u2_first_mb_in_slice = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    /*-----------------------------------------------------------------------*/
    /*---------------------------- slice type -------------------------------*/
    /*-----------------------------------------------------------------------*/
    ps_slice_prms->u1_slice_type = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(ps_slice_prms->u1_slice_type > 4)
    {
        ps_slice_prms->u1_slice_type -= 5;
    }

    /*-----------------------------------------------------------------------*/
    /*----------------------------- PPS id ----------------------------------*/
    /*-----------------------------------------------------------------------*/
    u1_pps_id = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    /* set correspoding sps and pps id also */
    ps_pps += u1_pps_id;
    if(FALSE == ps_pps->u1_is_valid)
    {
        return ERROR_INV_SLICE_HDR_T;
    }
    ps_sps = ps_pps->ps_sps;
    ps_subset_sps = &ps_svcd_ctxt->ps_subset_sps[ps_sps->u1_seq_parameter_set_id];
    if(CODED_SLICE_EXTENSION_NAL == i4_nal_unit_type)
    {
        ps_sps += MAX_NUM_SEQ_PARAMS;
        ps_subset_sps =
            &ps_svcd_ctxt->ps_subset_sps[MAX_NUM_SEQ_PARAMS + ps_sps->u1_seq_parameter_set_id];
    }
    /*-----------------------------------------------------------------------*/
    /*--------------------------- frm num -----------------------------------*/
    /*-----------------------------------------------------------------------*/
    if(!ps_sps) return ERROR_INV_SLICE_HDR_T;
    if(FALSE == ps_sps->u1_is_valid) return ERROR_INV_SLICE_HDR_T;

    ps_slice_prms->u2_frame_num = ih264d_get_bits_h264(ps_bitstrm, ps_sps->u1_bits_in_frm_num);

    /*-----------------------------------------------------------------------*/
    /*------------------ field pic flag and bottom field flag ---------------*/
    /*-----------------------------------------------------------------------*/
    if(SVCD_TRUE != ps_sps->u1_frame_mbs_only_flag)
    {
        return ERROR_INV_SLICE_HDR_T;
    }
    /*-----------------------------------------------------------------------*/
    /*--------------------------- IDR pic id --------------------------------*/
    /*-----------------------------------------------------------------------*/
    if(SVCD_TRUE == i4_idr_pic_flag)
    {
        UWORD32 u4_idr_pic_id = 0;
        u4_idr_pic_id = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_idr_pic_id > 65535) return ERROR_INV_SLICE_HDR_T;
    }

    /*-----------------------------------------------------------------------*/
    /*----------------- poc lsb and delts_poc_bottom ------------------------*/
    /*-----------------------------------------------------------------------*/
    if(0 == ps_sps->u1_pic_order_cnt_type)
    {
        i_temp = ih264d_get_bits_h264(ps_bitstrm, ps_sps->u1_log2_max_pic_order_cnt_lsb_minus);

        if(i_temp < 0 || i_temp >= ps_sps->i4_max_pic_order_cntLsb) return ERROR_INV_SLICE_HDR_T;
        if(SVCD_TRUE == ps_pps->u1_pic_order_present_flag)
        {
            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }
    }

    /*-----------------------------------------------------------------------*/
    /*---------------- delta_poc_count[0] and [1] ---------------------------*/
    /*-----------------------------------------------------------------------*/
    if((1 == ps_sps->u1_pic_order_cnt_type) && (!ps_sps->u1_delta_pic_order_always_zero_flag))
    {
        i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_pps->u1_pic_order_present_flag)
        {
            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }
    }

    /*-----------------------------------------------------------------------*/
    /*---------------------- redundant pic cnt ------------------------------*/
    /*-----------------------------------------------------------------------*/
    if(ps_pps->u1_redundant_pic_cnt_present_flag)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(u4_temp > MAX_REDUNDANT_PIC_CNT) return ERROR_INV_SLICE_HDR_T;
    }
    /*-----------------------------------------------------------------------*/
    /*-----------------Direct_spatial_mv_pred_flag --------------------------*/
    /*-----------------num ref active override flag -------------------------*/
    /*-----------------num_ref_idx_active_l0&1 ------------------------------*/
    /*-----------------------------------------------------------------------*/
    if(0 == i4_quality_id)
    {
        if(B_SLICE == ps_slice_prms->u1_slice_type)
        {
            ps_slice_prms->u1_direct_spatial_mv_pred_flag = ih264d_get_bit_h264(ps_bitstrm);
        }

        if((P_SLICE == ps_slice_prms->u1_slice_type) || (B_SLICE == ps_slice_prms->u1_slice_type))
        {
            WORD8 i1_over_ride_flag;
            i1_over_ride_flag = ih264d_get_bit_h264(ps_bitstrm);

            ps_slice_prms->u1_num_ref_idx_active_override_flag = i1_over_ride_flag;

            if(SVCD_TRUE == i1_over_ride_flag)
            {
                UWORD8 u8_ref_idx_l0;
                u8_ref_idx_l0 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                if(u8_ref_idx_l0 > H264_MAX_REF_PICS)
                {
                    return ERROR_NUM_REF;
                }
                if(B_SLICE == ps_slice_prms->u1_slice_type)
                {
                    UWORD8 u8_ref_idx_l1;
                    u8_ref_idx_l1 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                    if(u8_ref_idx_l1 > H264_MAX_REF_PICS)
                    {
                        return ERROR_NUM_REF;
                    }
                }
            }
        }

        /*-----------------------------------------------------------------------*/
        /*---------------------- ref pic list modification ----------------------*/
        /*-----------------------------------------------------------------------*/
        {
            ret = isvcd_parse_ref_pic_list_modify(ps_bitstrm, ps_slice_prms, ps_sps);
            if(OK != ret) return ret;
        }

        if(((1 == ps_pps->u1_wted_pred_flag) && (P_SLICE == ps_slice_prms->u1_slice_type)) ||
           ((B_SLICE == ps_slice_prms->u1_slice_type) && (1 == ps_pps->u1_wted_bipred_idc)))
        {
            if((ps_slice_prms->u1_num_ref_idx_lx_active[0] >= H264_MAX_REF_IDX) ||
               (ps_slice_prms->u1_num_ref_idx_lx_active[1] >= H264_MAX_REF_IDX))
            {
                return ERROR_NUM_REF;
            }
            /*-------------------------------------------------------------------*/
            /*------------------------- Pred weight table -----------------------*/
            /*-------------------------------------------------------------------*/
            if(CODED_SLICE_EXTENSION_NAL == i4_nal_unit_type)
            {
                WORD32 i4_base_pred_wt_tbl_flag = 1;

                /* base_pred_weight_table_flag */
                if(0 == i4_no_int_lyr_pred)
                {
                    i4_base_pred_wt_tbl_flag = ih264d_get_bit_h264(ps_bitstrm);
                }

                if((1 == i4_no_int_lyr_pred) || (0 == i4_base_pred_wt_tbl_flag))
                {
                    ret = ih264d_parse_pred_weight_table(ps_slice_prms, ps_bitstrm);
                    if(ret != OK) return ret;
                }
            }
            else
            {
                ret = ih264d_parse_pred_weight_table(ps_slice_prms, ps_bitstrm);
                if(ret != OK) return ret;
            }
        }

        /*-----------------------------------------------------------------------*/
        /*------------------------- ref pic marking -----------------------------*/
        /*-----------------------------------------------------------------------*/
        if(0 != i4_nal_ref_idc)
        {
            dec_struct_t *ps_dec;
            svc_dec_lyr_struct_t *ps_svc_lyr_dec;
            UWORD8 u1_store_ref_base_pic;
            ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr;
            ps_dec = &ps_svc_lyr_dec->s_dec;
            {
                dec_seq_params_t *ps_sps_tmp = ps_pps->ps_sps;

                ps_dec->u1_nal_unit_type = i4_nal_unit_type;
                ps_svc_lyr_dec->ps_nal_svc_ext->u1_idr_flag = i4_idr_pic_flag;
                ps_dec->ps_cur_sps = ps_sps;
                ps_dec->ps_cur_pps = ps_pps;
                ps_pps->ps_sps = ps_sps;

                if(ps_svc_lyr_dec->ps_nal_svc_ext->u1_idr_flag)
                    ps_dec->u1_nal_unit_type = IDR_SLICE_NAL;

                i_temp = ih264d_read_mmco_commands(ps_dec);
                ps_pps->ps_sps = ps_sps_tmp;
                ps_dec->u1_nal_unit_type = i4_nal_unit_type;
                if(i_temp < 0)
                {
                    return ERROR_DBP_MANAGER_T;
                }
                ps_dec->u4_bitoffset = i_temp;
            }

            if(0 == ps_subset_sps->s_sps_svc_ext.u1_slice_header_restriction_flag)
            {
                /* store_ref_base_pic_flag */
                u1_store_ref_base_pic = ih264d_get_bit_h264(ps_bitstrm);
                if(0 != u1_store_ref_base_pic)
                {
                    return ERROR_INV_SLICE_HDR_T;
                }

                if(((1 == i4_use_ref_base) || (1 == u1_store_ref_base_pic)) &&
                   (SVCD_FALSE == i4_idr_pic_flag))
                {
                    i_temp = isvcd_dec_ref_base_pic_marking(
                        &ps_svc_lyr_dec->s_svc_slice_params.s_ref_base_pic_marking_svc_ext,
                        ps_bitstrm);
                    if(i_temp != OK)
                    {
                        return i_temp;
                    }
                }
                /******* End of if (SVC_VCL_NAL == i4_nal_unit_type) *********/
            }
            /******** End of if(0 != i4_nal_ref_idc) *************************/
        }
        /************* End of if(0 == i4_quality_id) *************************/
    }

    /*-----------------------------------------------------------------------*/
    /*--------------------------- cabac int idc -----------------------------*/
    /*-----------------------------------------------------------------------*/
    if((ps_pps->u1_entropy_coding_mode == CABAC) && (I_SLICE != ps_slice_prms->u1_slice_type))
    {
        ps_slice_prms->u1_cabac_init_idc = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if(ps_slice_prms->u1_cabac_init_idc > MAX_CABAC_INIT_IDC)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
    }

    /*-----------------------------------------------------------------------*/
    /*--------------------------- slice qp delta ----------------------------*/
    /*-----------------------------------------------------------------------*/
    {
        WORD8 i1_slice_qp_delta;
        i1_slice_qp_delta = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        i1_slice_qp_delta += ps_pps->u1_pic_init_qp;
        if((i1_slice_qp_delta < MIN_H264_QP) || (i1_slice_qp_delta > MAX_H264_QP))
        {
            return ERROR_INV_RANGE_QP_T;
        }
        ps_slice_prms->u1_slice_qp = (UWORD8) i1_slice_qp_delta;
    }

    /*-----------------------------------------------------------------------*/
    /*--------------------------- disable dblk filter idc -------------------*/
    /*-----------------------------------------------------------------------*/
    /* Set to default value */

    ps_slice_prms->u1_disable_dblk_filter_idc = 0;
    if(SVCD_TRUE == ps_pps->u1_deblocking_filter_parameters_present_flag)
    {
        ps_slice_prms->u1_disable_dblk_filter_idc = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_slice_prms->u1_disable_dblk_filter_idc > SLICE_BOUNDARY_DBLK_DISABLED)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        /*-------------------------------------------------------------------*/
        /*--------------------------- slice_alpha_c0_offset_div2 ------------*/
        /*--------------------------- slice_beta_offset_div2 ----------------*/
        /*-------------------------------------------------------------------*/
        if(1 != ps_slice_prms->u1_disable_dblk_filter_idc)
        {
            /* slice_alpha_c0_offset_div2 */
            ps_slice_prms->i1_slice_alpha_c0_offset = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            if((MIN_DBLK_FIL_OFF > ps_slice_prms->i1_slice_alpha_c0_offset) ||
               (ps_slice_prms->i1_slice_alpha_c0_offset > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice_prms->i1_slice_beta_offset = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
            if((MIN_DBLK_FIL_OFF > ps_slice_prms->i1_slice_beta_offset) ||
               (ps_slice_prms->i1_slice_beta_offset > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
        }
    }

    *pi4_ref_dq_id = -1;

    if((0 == i4_no_int_lyr_pred) && (0 == i4_quality_id))
    {
        WORD32 i4_inter_lyr_dblk_idc;
        WORD32 i4_inter_lyr_alpha_c0_offset;
        WORD32 i4_inter_lyr_beta_offset;

        /*-------------------------------------------------------------------*/
        /*--------------------------- ref_layer_dq_id -----------------------*/
        /*-------------------------------------------------------------------*/
        *pi4_ref_dq_id = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(*pi4_ref_dq_id > MAX_REF_DEP_ID)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        /* ------------------------------------------- */
        /* ---- Inter layer de-blocking parameters ---- */
        /* ------------------------------------------- */
        i4_inter_lyr_dblk_idc = 0;
        i4_inter_lyr_alpha_c0_offset = 0;
        i4_inter_lyr_beta_offset = 0;
        if(SVCD_TRUE ==
           ps_subset_sps->s_sps_svc_ext.u1_inter_layer_deblocking_filter_control_present_flag)
        {
            i4_inter_lyr_dblk_idc = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

            if(i4_inter_lyr_dblk_idc > 6)
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            if(1 != i4_inter_lyr_dblk_idc)
            {
                /* Alpha Offset */
                i4_inter_lyr_alpha_c0_offset = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                if(i4_inter_lyr_alpha_c0_offset > 6 || i4_inter_lyr_alpha_c0_offset < -6)
                {
                    return ERROR_INV_SLICE_HDR_T;
                }
                i4_inter_lyr_alpha_c0_offset <<= 1;

                /* Beta Offset */
                i4_inter_lyr_beta_offset = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
                if(i4_inter_lyr_beta_offset > 6 || i4_inter_lyr_beta_offset < -6)
                {
                    return ERROR_INV_SLICE_HDR_T;
                }
                i4_inter_lyr_beta_offset <<= 1;
            }
        }
        ps_vcl_node->i4_inter_lyr_dblk_idc = i4_inter_lyr_dblk_idc;
        ps_vcl_node->i4_inter_lyr_beta_offset = i4_inter_lyr_beta_offset;
        ps_vcl_node->i4_inter_lyr_alpha_c0_offset = i4_inter_lyr_alpha_c0_offset;
    }

    return (0);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_ref_lyr_dqid                                   */
/*                                                                           */
/*  Description   : Parses the slice header till ref lyr dqid.               */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. vcl node                                              */
/*                  2. sequence prms set base buffer pointer                 */
/*                  3. picture prms set base buffer pointer                  */
/*                  4. Target layer flag                                     */
/*                  5. DPB command context structure                         */
/*                  6. Reference layer DQID                                  */
/*  Globals       : None                                                     */
/*  Processing    : 1. Parses the slice header till ref lyr DQID             */
/*                  2. If current layer is target layer then it calculates   */
/*                     poc and gaps in frame number                          */
/*                                                                           */
/*  Outputs       : Updates 1) ref dqid variable 2) dpb command context      */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           Draft                                */
/*****************************************************************************/
WORD32 isvcd_get_ref_lyr_dqid(vcl_node_t *ps_vcl_node, dec_seq_params_t *ps_sps,
                              dec_pic_params_t *ps_pps, WORD32 *pi4_ref_lyr_dqid,
                              WORD32 i4_prev_au_dqid, WORD32 *pi4_err_code,
                              svc_dec_ctxt_t *ps_svcd_ctxt)
{
    WORD32 i4_status;
    WORD32 ai4_ref_dq_id[2] = {0};
    WORD32 i4_num_slc_dec;

    /* local structures */
    dec_slice_params_t s_slice_prms = {0};

    /* vcl buffer */
    vcl_buf_hdr_t *ps_vcl_buf;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr;
    UNUSED(i4_prev_au_dqid);
    ps_dec = &ps_svc_lyr_dec->s_dec;
    /* Sanity checks */
    if((NULL == ps_vcl_node) || (NULL == ps_sps) || (NULL == ps_pps) || (NULL == pi4_ref_lyr_dqid))
    {
        return NOT_OK;
    }

    i4_num_slc_dec = 0;
    ps_vcl_buf = ps_vcl_node->ps_first_vcl_nal;
    i4_status = NOT_OK;

    while(NULL != ps_vcl_buf)
    {
        WORD32 i4_error;

        /* Fill the stream context structure */
        ps_dec->ps_bitstrm->u4_ofst = 0;
        ps_dec->ps_bitstrm->pu4_buffer =
            (UWORD32 *) ((UWORD8 *) ps_vcl_buf + ps_vcl_buf->i4_buf_offset +
                         ps_vcl_buf->i4_slice_offset);
        ps_dec->ps_bitstrm->u4_max_ofst = ps_vcl_buf->u4_max_bits;

        /* call the function which decodes the slice header */
        i4_error = isvcd_parse_slice_hdr_refdq_id(ps_vcl_node, &s_slice_prms, ps_dec->ps_bitstrm,
                                                  &ai4_ref_dq_id[i4_num_slc_dec], ps_sps, ps_pps,
                                                  ps_vcl_buf->i4_no_int_lyr_pred, ps_svcd_ctxt);

        /* store the first error encountered */
        if(0 == *pi4_err_code)
        {
            *pi4_err_code = i4_error;
        }
        if(i4_error != 0)
        {
            /* check on the Error returned */
            return NOT_OK;
        }

        /* set the return status */
        i4_status = OK;
        break;

        /* go to the next slice header */
        ps_vcl_buf = ps_vcl_buf->ps_next;
    }

    /* set the appropriate reference dqid of the first slice */
    *pi4_ref_lyr_dqid = ai4_ref_dq_id[0];

    return (i4_status);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_conceal_node_params                                */
/*                                                                           */
/*  Description   : This function detects the resolution change              */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. Pointer to Current SPS                                */
/*                  2. Pointer to prevoius SPS                               */
/*  Globals       : None                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : SVCD_TRUE if different resolution else SVCD_FALSE        */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijayakumar     Draft                                */
/*                                                                           */
/*****************************************************************************/
void isvcd_conceal_node_params(vcl_nal_t *ps_vcl_nal, prev_au_prms_t *ps_prev_au_prms)
{
    vcl_node_t *ps_node;
    WORD32 i4_conceal_lyrs;
    WORD32 i4_no_gaps_flag;

    /* get the bottom node */
    ps_node = ps_vcl_nal->ps_bot_node;
    i4_conceal_lyrs = SVCD_FALSE;
    i4_no_gaps_flag = SVCD_FALSE;

    /* loop over all nodes present in the current AU */
    while(NULL != ps_node)
    {
        WORD32 i4_dep_id = 0;
        WORD32 i4_qua_id = 0;
        UWORD16 u2_frm_num_dep = 0;
        WORD32 i4_idr_pic_flag = 0;
        WORD32 i4_idr_pic_num = 0;
        WORD32 i4_nal_ref_idc = 0;
        WORD32 i4_poc_syntax = 0;
        WORD32 i4_qua_zero_lyr_sts = 0;

        i4_dep_id = ps_node->i4_dependency_id;
        i4_qua_id = ps_node->i4_quality_id;

        /* reset the quality 0 layer updated status */
        if(0 == i4_qua_id)
        {
            i4_qua_zero_lyr_sts = SVCD_FALSE;
        }

        /* process the quality id 0 layers */
        if((0 == i4_qua_id) && (NULL != ps_node->ps_first_vcl_nal))
        {
            /* if current and previous are reference pictures */
            if((0 != ps_prev_au_prms[i4_dep_id].i4_nal_ref_id) && (0 != ps_node->i4_nal_ref_idc))
            {
                if(ps_prev_au_prms[i4_dep_id].u2_frm_num == ps_node->u2_frm_num)
                {
                    /* frame number is concealed */
                    ps_node->u2_frm_num++;
                    i4_conceal_lyrs = SVCD_TRUE;
                }
                else if((SVCD_TRUE == i4_conceal_lyrs) || (SVCD_TRUE == i4_no_gaps_flag))
                {
                    /* if the current au frm_num is less than prev */
                    /* or the difference is greater than 1         */
                    if((ps_prev_au_prms[i4_dep_id].u2_frm_num > ps_node->u2_frm_num) ||
                       ((ps_node->u2_frm_num - ps_prev_au_prms[i4_dep_id].u2_frm_num) > 1))
                    {
                        /* frame number is concealed */
                        ps_node->u2_frm_num = ps_prev_au_prms[i4_dep_id].u2_frm_num + 1;
                    }
                }

                /* set the no gaps flag */
                if(1 == (ps_node->u2_frm_num - ps_prev_au_prms[i4_dep_id].u2_frm_num))
                {
                    i4_no_gaps_flag = SVCD_TRUE;
                }
            }

            /* store the final frame number */
            u2_frm_num_dep = ps_node->u2_frm_num;
            i4_idr_pic_flag = ps_node->i4_idr_pic_flag;
            i4_idr_pic_num = ps_node->i4_idr_pic_num;
            i4_nal_ref_idc = ps_node->i4_nal_ref_idc;
            i4_poc_syntax = ps_node->i4_poc_syntax;
            i4_qua_zero_lyr_sts = SVCD_TRUE;
        }
        else
        {
            if(SVCD_TRUE == i4_qua_zero_lyr_sts)
            {
                /* for higher quality layers store the same value */
                /* present in the quality id 0 layer              */
                ps_node->u2_frm_num = u2_frm_num_dep;
                ps_node->i4_idr_pic_flag = i4_idr_pic_flag;
                ps_node->i4_idr_pic_num = i4_idr_pic_num;
                ps_node->i4_nal_ref_idc = i4_nal_ref_idc;
                ps_node->i4_poc_syntax = i4_poc_syntax;
            }
        }

        /* get the upper node pointer */
        ps_node = ps_node->ps_top_node;
    }
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_refine_dep_list                                    */
/*                                                                           */
/*  Description   : Refines the DQID list based on reference layer DQID      */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. vcl nal structure (input and output)                  */
/*                  2. sps prms base buffer pointer                          */
/*                  3. pps prms base buffer pointer                          */
/*                  4. pointer to array in which the dep id should be stored */
/*                  5. pointer to array having prev AU ref dq id             */
/*                  6. pointer to init params structure                      */
/*                  7. pointer to store the Error code                       */
/*  Globals       : None                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : Updates the vcl nal structure                            */
/*                  Also determines frm_num and poc                          */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           Draft                                */
/*****************************************************************************/
WORD32 isvcd_refine_dep_list(void *pv_out_vcl_ctxt, dec_seq_params_t *ps_sps,
                             dec_svc_seq_params_t *ps_subset_sps, dec_pic_params_t *ps_pps,
                             WORD32 *pi4_dep_id_map, prev_au_prms_t *ps_prev_au_prms,
                             prev_au_sps_pps_t *ps_pps_sps_prev, WORD32 *pi4_err_code,
                             svc_dec_ctxt_t *ps_svcd_ctxt)
{
    vcl_nal_t *ps_vcl_nal;
    vcl_node_t *ps_vcl_node;
    WORD32 i4_idr_pic_flag;
    WORD32 i4_nal_ref_idc;
    WORD32 i4_idr_pic_num;
    WORD32 i4_num_res_lyrs_bup;
    WORD32 i4_restore_prms_flag;
    vcl_node_t *ps_node_bup;
    WORD32 ai4_dep_id[MAX_NUM_RES_LYRS] = {0};

    /* used for checking the init prms */
    dec_seq_params_t *ps_sps_tgt_minus1_lyr = NULL;
    dec_seq_params_t *ps_sps_tgt_minus2_lyr = NULL;
    UNUSED(pi4_err_code);
    /* sanity checks */
    if((NULL == pv_out_vcl_ctxt) || (NULL == ps_sps) || (NULL == ps_pps))
    {
        return NOT_OK;
    }

    ps_vcl_nal = (vcl_nal_t *) pv_out_vcl_ctxt;

    /*  no node is present */
    if(NULL == ps_vcl_nal->ps_bot_node)
    {
        return (NOT_OK);
    }

    /* set the single layer flag if top node and bottom node are same */
    if((ps_vcl_nal->ps_bot_node == ps_vcl_nal->ps_top_node) &&
       (0 == ps_vcl_nal->ps_bot_node->i4_dependency_id))
    {
    }
    else
    {
        /* call the function which corrects the frame number of each node */
        /* based on previous access unit frame number                     */
        isvcd_conceal_node_params(ps_vcl_nal, ps_prev_au_prms);
    }
    /* get the top most node */
    ps_vcl_node = ps_vcl_nal->ps_top_node;

    /* get the IDR picture flag for top most layer in current AU */
    /* if not valid then set the value present in the first valid node */
    {
        vcl_node_t *ps_node;
        WORD32 i4_node_present_flag;

        ps_node = ps_vcl_node;
        i4_node_present_flag = SVCD_FALSE;

        /* store default values */
        i4_idr_pic_flag = SVCD_FALSE;
        i4_nal_ref_idc = 0;
        i4_idr_pic_num = 0;

        /* loop until valid node */
        while(NULL != ps_node)
        {
            if(NULL != ps_node->ps_first_vcl_nal)
            {
                i4_idr_pic_flag = ps_node->i4_idr_pic_flag;
                i4_nal_ref_idc = ps_node->i4_nal_ref_idc;
                i4_idr_pic_num = ps_node->i4_idr_pic_num;
                i4_node_present_flag = SVCD_TRUE;
                break;
            }
            else if(SVCD_TRUE == ps_node->i4_idr_pic_flag)
            {
                i4_idr_pic_flag = ps_node->i4_idr_pic_flag;
                i4_nal_ref_idc = ps_node->i4_nal_ref_idc;
                i4_idr_pic_num = ps_node->i4_idr_pic_num;
                i4_node_present_flag = SVCD_TRUE;
                break;
            }
            /* point to next node */
            ps_node = ps_node->ps_bot_node;
        }

        /* alteast one node should be present */
        if(SVCD_FALSE == i4_node_present_flag)
        {
            return (NOT_OK);
        }
    }

    /* initially the access unit is considered to have a single resolution */
    ai4_dep_id[0] = 0;
    ps_vcl_nal->i4_num_res_lyrs = 1;
    i4_restore_prms_flag = SVCD_FALSE;

    /*-----------------------------------------------------------------------*/
    /* loop until all the nodes are processed                                */
    /*-----------------------------------------------------------------------*/
    while(NULL != ps_vcl_node)
    {
        WORD32 i4_ref_lyr_dqid, i4_status;
        vcl_node_t *ps_bot_vcl_node;
        WORD32 i4_res_chnge_flag = SVCD_FALSE;
        WORD32 i4_dep_id, i4_qua_id;
        WORD32 i4_prev_sps_pps_valid;
        WORD32 i4_prev_au_prms_valid;

        /* set the reference layer DQID to -1 */
        i4_ref_lyr_dqid = -1;

        /* get the current layer dependency and quality id */
        i4_dep_id = ps_vcl_node->i4_dependency_id;
        i4_qua_id = ps_vcl_node->i4_quality_id;

        /* get the valid status of prev access unit params */
        i4_prev_au_prms_valid = ps_prev_au_prms[i4_dep_id].i4_updated_sts;
        i4_prev_sps_pps_valid = ps_pps_sps_prev[(i4_dep_id << 4) + i4_qua_id].i4_updated_sts;

        /* missing layer handling */
        if(NULL == ps_vcl_node->ps_first_vcl_nal)
        {
            /* store the params appropriately */
            ps_vcl_node->i4_idr_pic_flag = i4_idr_pic_flag;
            ps_vcl_node->i4_nal_ref_idc = i4_nal_ref_idc;
            ps_vcl_node->i4_idr_pic_num = i4_idr_pic_num;
            ps_vcl_node->i4_num_slices = 0;
            ps_vcl_node->i4_use_ref_base = 0;
            ps_vcl_node->i4_temporal_id = 0;

            if((0 != i4_dep_id) || (0 != i4_qua_id))
            {
                ps_vcl_node->i4_nal_unit_type = CODED_SLICE_EXTENSION_NAL;
                ps_vcl_node->u1_acc_no_int_pred = 0;
            }
            else if(SVCD_TRUE == i4_idr_pic_flag)
            {
                ps_vcl_node->i4_nal_unit_type = IDR_SLICE_NAL;
                ps_vcl_node->u1_acc_no_int_pred = 1;
            }
            else
            {
                ps_vcl_node->i4_nal_unit_type = SLICE_NAL;
                ps_vcl_node->u1_acc_no_int_pred = 1;
            }

            if(SVCD_FALSE == i4_idr_pic_flag)
            {
                /* pick the other params form previous access unit */
                if(SVCD_TRUE == i4_prev_sps_pps_valid)
                {
                    ps_vcl_node->u1_pps_id =
                        ps_pps_sps_prev[(i4_dep_id << 4) + i4_qua_id].u1_pps_id;

                    ps_vcl_node->u1_sps_id =
                        ps_pps_sps_prev[(i4_dep_id << 4) + i4_qua_id].u1_sps_id;
                }

                if(SVCD_TRUE == i4_prev_au_prms_valid)
                {
                    if(0 == ps_vcl_node->i4_nal_ref_idc)
                    {
                        ps_vcl_node->u2_frm_num = ps_prev_au_prms[i4_dep_id].u2_frm_num;
                    }
                    else
                    {
                        ps_vcl_node->u2_frm_num = ps_prev_au_prms[i4_dep_id].u2_frm_num + 1;
                    }
                }
            }
        }

        /* SPS id cannot change unless its an IDR pic */
        if(SVCD_FALSE == ps_vcl_node->i4_idr_pic_flag)
        {
            if(SVCD_TRUE == i4_prev_sps_pps_valid)
            {
                /* store the SPS id of the current layer */
                ps_vcl_node->u1_sps_id = ps_pps_sps_prev[(i4_dep_id << 4) + i4_qua_id].u1_sps_id;
            }
        }

        /* store the PPS id and SPS id of the current layer */
        ps_pps_sps_prev[(i4_dep_id << 4) + i4_qua_id].u1_pps_id = ps_vcl_node->u1_pps_id;
        ps_pps_sps_prev[(i4_dep_id << 4) + i4_qua_id].u1_sps_id = ps_vcl_node->u1_sps_id;
        ps_pps_sps_prev[(i4_dep_id << 4) + i4_qua_id].i4_updated_sts = SVCD_TRUE;

        /* handling of no_inter_layer_pred_flag 1 cases */
        if((1 == ps_vcl_node->u1_acc_no_int_pred) && (NULL != ps_vcl_node->ps_bot_node))
        {
            if(SVCD_TRUE == i4_idr_pic_flag)
            {
                /* take a back up of the parameters till the current node. */
                /* these parameters will be restored at the end of loop */

                if(SVCD_FALSE == i4_restore_prms_flag)
                {
                    /* get the number of resolution detected so far */
                    i4_num_res_lyrs_bup = ps_vcl_nal->i4_num_res_lyrs;

                    ps_node_bup = ps_vcl_node;

                    /* set the restore params flag */
                    i4_restore_prms_flag = SVCD_TRUE;
                }
            }
            else
            {
                ps_vcl_node->i4_ref_dq_id = -1;
                ps_vcl_node->i4_res_change_flag = i4_res_chnge_flag;

                /* store the reference DQID for current dependency */
                ps_prev_au_prms[i4_dep_id].i4_ref_dq_id = -1;
                ps_prev_au_prms[i4_dep_id].u2_frm_num = ps_vcl_node->u2_frm_num;
                ps_prev_au_prms[i4_dep_id].i4_nal_ref_id = ps_vcl_node->i4_nal_ref_idc;

                /* the bottom node is set to NULL */
                ps_vcl_node->ps_bot_node = NULL;
                break;
            }
        }

        /* derive the reference layer DQID for quality id equal to 0 */
        if(0 == i4_qua_id)
        {
            dec_seq_params_t *ps_curr_sps;
            dec_svc_seq_params_t *ps_curr_subset_sps;

            /* derive current SPS */
            ps_curr_sps = ps_sps + ps_vcl_node->u1_sps_id;
            ps_curr_subset_sps = ps_subset_sps + ps_vcl_node->u1_sps_id;

            {
                WORD32 i4_max_frm_num;

                /* get the maximum value of frame number */
                i4_max_frm_num = (1 << (ps_curr_sps->u1_bits_in_frm_num + 1));
                ps_vcl_node->u2_frm_num = ps_vcl_node->u2_frm_num % i4_max_frm_num;
                if(SVCD_TRUE == ps_vcl_node->i4_idr_pic_flag)
                {
                    /* if idr then frm num should be 0 */
                    ps_vcl_node->u2_frm_num = 0;
                }
            }

            /* store default params to inter layer deblocking params */
            ps_vcl_node->i4_inter_lyr_dblk_idc = 0;
            ps_vcl_node->i4_inter_lyr_beta_offset = 0;
            ps_vcl_node->i4_inter_lyr_alpha_c0_offset = 0;
            /* No SEI support for scalability info*/
            i4_status = NOT_OK;

            /* if no inter layer pred flag is present set the   */
            /* status to fail since the slices will not contain */
            /* reference layer Dqid                             */
            if(1 == ps_vcl_node->u1_acc_no_int_pred)
            {
                i4_status = NOT_OK;
            }
            else
            {
                WORD32 *pi4_ref_dq_id;
                WORD32 i4_ref_dq_id_temp;

                /* check if the SEI message has given the ref_dq_id */
                if(NOT_OK == i4_status)
                {
                    pi4_ref_dq_id = &i4_ref_lyr_dqid;
                }
                else
                {
                    pi4_ref_dq_id = &i4_ref_dq_id_temp;
                }

                i4_status = isvcd_get_ref_lyr_dqid(ps_vcl_node, ps_sps, ps_pps, pi4_ref_dq_id,
                                                   ps_prev_au_prms[i4_dep_id].i4_ref_dq_id,
                                                   &ps_svcd_ctxt->i4_error_code, ps_svcd_ctxt);
            }

            /* no slice in the layer has been successfully decoded */
            if(NOT_OK == i4_status)
            {
                /* check for IDR picture */
                if(SVCD_TRUE == i4_idr_pic_flag)
                {
                    /* set the next lower layer as the reference layer */
                    if(NULL != ps_vcl_node->ps_bot_node)
                    {
                        i4_ref_lyr_dqid = ps_vcl_node->ps_bot_node->i4_dependency_id << 4;

                        i4_ref_lyr_dqid += ps_vcl_node->ps_bot_node->i4_quality_id;
                    }
                    else
                    {
                        i4_ref_lyr_dqid = -1;
                    }
                }
                else
                {
                    /* take the reference dq id from previous access unit */
                    i4_ref_lyr_dqid = ps_prev_au_prms[i4_dep_id].i4_ref_dq_id;
                }
            }

            /* Update the DQID list based on ref DQID.     */
            /* This routine also updates the ref_dq_id     */
            /* in case the actual layer is completely lost */
            i4_status = isvcd_update_dqid(i4_ref_lyr_dqid, ps_vcl_node, &ps_bot_vcl_node);

            if(!(OK == i4_status))
            {
                return i4_status;
            }

            /* store the reference DQID for current depedency and */
            /* quality id 0 layer                                 */
            ps_prev_au_prms[i4_dep_id].i4_ref_dq_id = i4_ref_lyr_dqid;
            ps_prev_au_prms[i4_dep_id].i4_nal_ref_id = ps_vcl_node->i4_nal_ref_idc;
            ps_prev_au_prms[i4_dep_id].u2_frm_num = ps_vcl_node->u2_frm_num;
            ps_prev_au_prms[i4_dep_id].i4_updated_sts = SVCD_TRUE;

            /* ------- Detect Resolution Change ---------------- */
            {
                dec_seq_params_t *ps_lower_sps = NULL;
                dec_svc_seq_params_t *ps_lower_subset_sps = NULL;

                if(NULL != ps_bot_vcl_node)
                {
                    if((NULL != ps_bot_vcl_node->ps_first_vcl_nal) ||
                       (SVCD_TRUE == i4_idr_pic_flag))
                    {
                        /* get the SPS of layer */
                        ps_lower_sps = ps_sps + ps_bot_vcl_node->u1_sps_id;
                        ps_lower_subset_sps = ps_subset_sps + ps_bot_vcl_node->u1_sps_id;
                    }
                    else
                    {
                        /* if the bottom layer is completely missed */
                        WORD32 i4_bot_dep_id, i4_bot_qua_id;
                        UWORD8 u1_sps_id = 0;

                        /* sps id is picked from previous access unit */
                        i4_bot_dep_id = ps_bot_vcl_node->i4_dependency_id;
                        i4_bot_qua_id = ps_bot_vcl_node->i4_quality_id;

                        if(SVCD_TRUE ==
                           ps_pps_sps_prev[(i4_bot_dep_id << 4) + i4_bot_qua_id].i4_updated_sts)
                        {
                            u1_sps_id =
                                ps_pps_sps_prev[(i4_bot_dep_id << 4) + i4_bot_qua_id].u1_sps_id;
                        }
                        else
                        {
                            /* should not enter here */
                            return NOT_OK;
                        }

                        /* get the SPS of lower layer */
                        ps_lower_sps = ps_sps + u1_sps_id;
                        ps_lower_subset_sps = ps_subset_sps + u1_sps_id;
                    }
                }

                /* call the function which detects resolution change */
                i4_res_chnge_flag = isvcd_detect_res_change(
                    ps_curr_sps, ps_lower_sps, ps_curr_subset_sps, ps_lower_subset_sps);

                /* if a resolution exists below current resolution */
                if(SVCD_TRUE == i4_res_chnge_flag)
                {
                    /* if current picture id IDR */
                    if(SVCD_TRUE == i4_idr_pic_flag)
                    {
                        /* store the depedency id of bottom most layer in current resolution */
                        ai4_dep_id[ps_vcl_nal->i4_num_res_lyrs - 1] = i4_dep_id;
                    }

                    /* increment the num resolution layer counter */
                    ps_vcl_nal->i4_num_res_lyrs++;

                    /* store the SPS of target -1 and -2 resolution layers */
                    if(2 == ps_vcl_nal->i4_num_res_lyrs)
                    {
                        ps_sps_tgt_minus1_lyr = ps_curr_sps;
                    }
                    else if(3 == ps_vcl_nal->i4_num_res_lyrs)
                    {
                        ps_sps_tgt_minus2_lyr = ps_curr_sps;
                    }
                    else if(ps_vcl_nal->i4_num_res_lyrs > MAX_NUM_RES_LYRS)
                    {
                        return NOT_OK;
                    }
                }
            }

            /* -------- end of resolution change detection -------- */
        }
        else
        {
            i4_ref_lyr_dqid = (i4_dep_id << 4);
            i4_ref_lyr_dqid += (i4_qua_id - 1);

            /* Update the DQID list based on ref DQID.     */
            /* This routine also updates the ref_dq_id     */
            /* in case the actual layer is completely lost */
            i4_status = isvcd_update_dqid(i4_ref_lyr_dqid, ps_vcl_node, &ps_bot_vcl_node);

            if(!(OK == i4_status))
            {
                return i4_status;
            }
            if(SVCD_TRUE == ps_vcl_node->i4_idr_pic_flag)
            {
                /* if idr then frm num should be 0 */
                ps_vcl_node->u2_frm_num = 0;
            }
        }

        /* Update resolution change flag inside VCL    */
        /* node structure. This parameter is later used*/
        /* in detecting the top most layer in the      */
        /* resolution currently being decoded          */
        ps_vcl_node->i4_res_change_flag = i4_res_chnge_flag;
        ps_vcl_node->i4_ref_dq_id = i4_ref_lyr_dqid;

        /* go to the next node */
        ps_vcl_node = ps_bot_vcl_node;
    }

    /* update the Dependency array for each resolution */
    if(SVCD_TRUE == i4_idr_pic_flag)
    {
        WORD32 i4_idx;

        ai4_dep_id[ps_vcl_nal->i4_num_res_lyrs - 1] = 0;

        /* loop over number of resolutions detected */
        for(i4_idx = 0; i4_idx < ps_vcl_nal->i4_num_res_lyrs; i4_idx++)
        {
            pi4_dep_id_map[i4_idx] = ai4_dep_id[ps_vcl_nal->i4_num_res_lyrs - 1 - i4_idx];
        }
    }

    if(SVCD_TRUE == i4_restore_prms_flag)
    {
        /* restore the number of resolutions */
        ps_vcl_nal->i4_num_res_lyrs = i4_num_res_lyrs_bup;

        ps_vcl_node = ps_node_bup;

        /* set the bottom node to NULL */
        ps_vcl_node->ps_bot_node = NULL;

        ps_vcl_node->i4_ref_dq_id = -1;
        ps_vcl_node->i4_res_change_flag = SVCD_FALSE;

        /* store the reference DQID for current dependency */
        ps_prev_au_prms[ps_vcl_node->i4_dependency_id].i4_ref_dq_id = -1;

        ps_prev_au_prms[ps_vcl_node->i4_dependency_id].u2_frm_num = ps_vcl_node->u2_frm_num;

        ps_prev_au_prms[ps_vcl_node->i4_dependency_id].i4_nal_ref_id = ps_vcl_node->i4_nal_ref_idc;
    }

    /* Finally update the bottom most node in the current access unit */
    ps_vcl_node = ps_vcl_nal->ps_top_node;

    while(NULL != ps_vcl_node->ps_bot_node)
    {
        ps_vcl_node = ps_vcl_node->ps_bot_node;
    }

    ps_vcl_nal->ps_bot_node = ps_vcl_node;

    /* check on validity of Target Layer -1 and -2 dimensions */
    if((NULL != ps_sps_tgt_minus1_lyr) && (0 == ps_sps_tgt_minus1_lyr->u1_is_valid))
    {
        if((H264_MAX_FRAME_WIDTH < (WORD32) (ps_sps_tgt_minus1_lyr->u2_frm_wd_in_mbs << 4)) ||
           (H264_MAX_FRAME_HEIGHT < (WORD32) (ps_sps_tgt_minus1_lyr->u2_frm_ht_in_mbs << 4)))
        {
            return NOT_OK;
        }
    }

    if((NULL != ps_sps_tgt_minus2_lyr) && (0 == ps_sps_tgt_minus2_lyr->u1_is_valid))
    {
        if((H264_MAX_FRAME_WIDTH < (WORD32) (ps_sps_tgt_minus2_lyr->u2_frm_wd_in_mbs << 4)) ||
           (H264_MAX_FRAME_HEIGHT < (WORD32) (ps_sps_tgt_minus2_lyr->u2_frm_ht_in_mbs << 4)))
        {
            return NOT_OK;
        }
    }

    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_dec_non_vcl                                        */
/*                                                                           */
/*  Description   : this function decodes the NON VCL NAL units              */
/*                                                                           */
/*                                                                           */
/*  Inputs        : pv_out_non_vcl : pointer to the structure containing     */
/*                                  NON VCL NAL units                        */
/*                  ps_seq_params : pointer to array of SPS structures       */
/*                  ps_pic_params : pointer to array of PPS structures       */
/*                  ps_sei_ctxt : pointer to array of SEI structures         */
/*  Globals       : none                                                     */
/*  Processing    : it decodes the units unitl all the units are             */
/*                  decoded                                                  */
/*  Outputs       : decoded parameters in appropriate structures             */
/*  Returns       : Success or Faliure                                       */
/*                                                                           */
/*  Issues        :                                                          */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_dec_non_vcl(void *pv_out_non_vcl, void *pv_seq_params, void *pv_pic_params,
                         svc_dec_ctxt_t *ps_svcd_ctxt)
{
    /* local varibles */
    non_vcl_nal_t *ps_non_vcl;
    WORD32 i4_unit_indx;
    non_vcl_buf_hdr_t *ps_non_vcl_buf;
    WORD32 i_status = OK;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    dec_bit_stream_t *ps_bitstrm;

    if((NULL == pv_out_non_vcl) || (NULL == pv_seq_params) || (NULL == pv_pic_params))
    {
        return NOT_OK;
    }
    UNUSED(pv_seq_params);
    UNUSED(pv_pic_params);

    /* currently SEI decoding is not supported */
    /* derive the local variables */
    ps_non_vcl = (non_vcl_nal_t *) pv_out_non_vcl;
    ps_non_vcl_buf = ps_non_vcl->ps_first_non_vcl_nal;
    if(NULL == ps_non_vcl_buf) return (NOT_OK);

    /* loop until all NON VCL NAL are decoded */
    for(i4_unit_indx = 0; i4_unit_indx < ps_non_vcl->i4_num_non_vcl_nals; i4_unit_indx++)
    {
        UWORD32 u4_nal_unit_type;
        ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr;
        ps_dec = &ps_svc_lyr_dec->s_dec;
        if(NULL == ps_non_vcl_buf) return (NOT_OK);
        /* get the current NAL unit type */
        u4_nal_unit_type = (UWORD32) ps_non_vcl_buf->i4_nal_unit_type;
        if(u4_nal_unit_type > MAX_SVC_NAL_UNIT_TYPE) return (NOT_OK);
        ps_dec->u1_nal_unit_type = u4_nal_unit_type;

        ps_dec->ps_bitstrm->pu4_buffer =
            (UWORD32 *) ((UWORD8 *) ps_non_vcl_buf + ps_non_vcl_buf->i4_buf_offset);
        ps_dec->ps_bitstrm->u4_ofst = 0;
        ps_dec->ps_bitstrm->u4_max_ofst = isvcd_nal_rbsp_to_sodb(
            (UWORD8 *) ps_dec->ps_bitstrm->pu4_buffer, ps_non_vcl_buf->i4_buf_size, 0);
        if(ps_dec->ps_bitstrm->u4_max_ofst <= 0) return (NOT_OK);

        ps_bitstrm = ps_dec->ps_bitstrm;

        /* call the processing module based on nal unit type */
        switch(u4_nal_unit_type)
        {
            case SEQ_PARAM_NAL:

                i_status = isvcd_parse_sps(ps_svc_lyr_dec, ps_bitstrm);

                if(!i_status)
                {
                    ps_dec->i4_header_decoded |= 0x1;
                    ps_svcd_ctxt->u4_num_sps_ctr++;
                }

                if(i_status) return i_status;

                break;
            case SUBSET_SPS_NAL:

                i_status = isvcd_parse_subset_sps(ps_svc_lyr_dec, ps_bitstrm);

                if(!i_status)
                {
                    ps_svcd_ctxt->u4_num_sps_ctr++;
                    ps_dec->i4_header_decoded |= 0x1;
                }
                if(i_status) return i_status;

                break;

            case PIC_PARAM_NAL:

                i_status = isvcd_parse_pps(ps_svc_lyr_dec, ps_bitstrm);
                if(i_status == ERROR_INV_SPS_PPS_T) return i_status;
                if(!i_status)
                {
                    ps_dec->i4_header_decoded |= 0x2;
                    ps_svcd_ctxt->u4_num_pps_ctr++;
                }
                break;
            case SEI_NAL:
            {
                i_status = ih264d_parse_sei_message(ps_dec, ps_bitstrm);
                ih264d_parse_sei(ps_dec, ps_bitstrm);
            }
            break;
            default:
                /* no other NON VCL UNIT is supported */
                break;
        }

        /* get the next non vcl bufffer */
        ps_non_vcl_buf = ps_non_vcl_buf->ps_next;

    } /* end of loop over all NAL units */

    return (OK);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_seq_hdr_dec                                        */
/*                                                                           */
/*  Description   : This function decodes sequence header, which includes    */
/*                  non VCL NAL before the first VCL unit                    */
/*  Inputs        : Decoder context, inbufs, place holder for number of bytes*/
/*                  consumed and number of packets consumed                  */
/*  Globals       : None                                                     */
/*  Processing    : 1. Parse non VCL units before first VCL unit             */
/*                  2. Decode parsed non VCL units                           */
/*  Outputs       : Decoded header                                           */
/*  Returns       : OK or NOT_OK                                             */
/*                                                                           */
/*  Issues        : no known issues                                          */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Creation                             */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_seq_hdr_dec(svc_dec_ctxt_t *ps_svcd_ctxt, ivd_video_decode_ip_t *ps_in_bufs,
                         UWORD32 *pu4_bytes_consumed)
{
    WORD32 i4_status;

    /* Decode all non VCL NAL till first VCL NAL is encountered */
    ps_svcd_ctxt->s_non_vcl_nal.i4_num_non_vcl_nals = 0;
    i4_status = isvcd_nal_parse_non_vcl_nal(
        ps_svcd_ctxt->pv_nal_parse_ctxt, ps_in_bufs->pv_stream_buffer, &ps_svcd_ctxt->s_non_vcl_nal,
        pu4_bytes_consumed, &ps_in_bufs->u4_num_Bytes);

    /* Note: The bitstream extraction module expects updated  */
    /* pointer whenever a new call to this module has been    */
    /* made. Hence the buffer pointer has to be incremented   */
    /* by bytes consumed                                      */
    ps_in_bufs->u4_num_Bytes -= *pu4_bytes_consumed;

    /* ------------------------------------------------------ */
    /* Decoding of non VCL data. As current implementation it */
    /* decodes the followings:                                */
    /*          1. Sequence parameter set                     */
    /*          2. Picture parameter set                      */
    /*          3. SEI message                                */
    /* ------------------------------------------------------ */
    isvcd_dec_non_vcl(&ps_svcd_ctxt->s_non_vcl_nal, ps_svcd_ctxt->ps_sps, ps_svcd_ctxt->ps_pps,
                      ps_svcd_ctxt);

    return (i4_status);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pre_parse_refine_au                               */
/*                                                                           */
/*  Description   : This function process a decode call                      */
/*  Inputs        : ps_dec_ctxt : decoder context structure                  */
/*                  ps_in_bufs : input buffer descriptor                     */
/*                  pu4_bytes_consumed : pointer to store the bytes consumed */
/*                  pi4_packets_consumed : pointer to store the packets      */
/*                                   consumed                                */
/*  Globals       : None                                                     */
/*  Processing    : It calls the NAL parse module to parse the input stream  */
/*                  if a picture boundary is detected it calls the           */
/*                  Dependency list refiniment and Picture Decode routines   */
/*  Outputs       : Decoded picture                                          */
/*  Returns       : OK or NOT_OK                                             */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Creation                             */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_pre_parse_refine_au(svc_dec_ctxt_t *ps_svcd_ctxt, ivd_video_decode_ip_t *ps_in_bufs,
                                 UWORD32 *pu4_bytes_consumed)
{
    WORD32 i4_status, i4_non_vcl_status;
    UWORD32 u4_bytes_consumed = 0;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr;
    ps_dec = &ps_svc_lyr_dec->s_dec;

    /* Sequence header decode:                                 */
    /* If sequence header is not decoded then decode  the seq  */
    /* uence header                                            */

    if(SVCD_FALSE == ps_dec->i4_header_decoded)
    {
        i4_status = isvcd_seq_hdr_dec(ps_svcd_ctxt, ps_in_bufs, &u4_bytes_consumed);

        if((VCL_NAL_FOUND_TRUE == i4_status) && (ps_svcd_ctxt->u4_num_sps_ctr != 0) &&
           (ps_svcd_ctxt->u4_num_pps_ctr != 0))
        {
            /* set the header decoded flag */
            ps_dec->i4_header_decoded = 3;
        }
    }
    *pu4_bytes_consumed = u4_bytes_consumed;
    if(1 == ps_dec->i4_decode_header)
    {
        return OK;
    }
    /* Bit-stream Parsing. It performs following tasks:        */
    /*          1. NAL hader decoder                           */
    /*          2. Emulation prevention and byte swap          */
    /*             (During this process data to moved to output*/
    /*              buffer)                                    */
    /*          3. Dependency list creation based on NAL header*/
    /*          4. Detection of picture boundary               */
    /* NOTE1:                                                  */
    /*       Output buffers for VCL and non VCL data are       */
    /*       different. VCL data can be retrieved through      */
    /*       dependency list. Whereas non VCL data is stored in*/
    /*       one single buffer, which is accessed through NON  */
    /*       VCL structure                                     */
    /* NOTE2:Partial input case for nal parsing requires a     */
    /*       flush API to be called when end of bitstream      */
    /*       occurs                                            */

    if(SVCD_FALSE == ps_svcd_ctxt->i4_eos_flag)
    {
        if(ps_dec->i4_header_decoded == 3)
        {
            i4_status = isvcd_nal_parse_vcl_nal_partial(
                ps_svcd_ctxt->pv_nal_parse_ctxt, ps_in_bufs->pv_stream_buffer,
                &ps_svcd_ctxt->s_non_vcl_nal, &ps_svcd_ctxt->s_vcl_nal, &u4_bytes_consumed,
                &ps_in_bufs->u4_num_Bytes);
        }
        else
        {
            return NOT_OK;
        }
    }
    else
    {
        void *pv_nal_parse_ctxt;
        pv_nal_parse_ctxt = ps_svcd_ctxt->pv_nal_parse_ctxt;

        i4_status = isvcd_nal_parse_partial_signal_eos(pv_nal_parse_ctxt, &ps_svcd_ctxt->s_vcl_nal,
                                                       &ps_svcd_ctxt->s_non_vcl_nal);

        u4_bytes_consumed = 0;
    }

    *pu4_bytes_consumed += u4_bytes_consumed;

    /* Picture Boundary detected: Go ahead and do the decoding  */
    /* Picture boundary not detected: Otherwsie retrun from this*/
    /* function and update the bytes consumed variable. This    */
    /* should be repeated till we get a picture boundary        */

    if(PIC_BOUNDARY_FALSE == i4_status)
    {
        return (NOT_OK);
    }

    else if(FLUSH_DECODED_PICTURE == i4_status)
    {
        /* No more data is expected to come. Pictures decoded   */
        /* so far needs to be sent for display                  */
        return (FLUSH);
    }

    if(PIC_BOUNDARY_TRUE != i4_status)
    {
        return (NOT_OK);
    }

    /* check if the application has set any of the skip modes       */
    /* add the support for P and B skip modes                       */
    /* if(ps_dec_ctxt->s_dyn_prms.u1_frame_skip_mode)               */

    /* Parse slice header to decode reference layer dQId and refine */
    /* the dependency list                                          */
    /* NOTE: Yes, this processing could be moved into NAL parsing   */
    /*       routine to avoid unneccessary emulation prevention and */
    /*       byte swapping over discardable data. This Optimization */
    /*       has been deferred for some time. In future if we found */
    /*       that there are many such streams which doesn't set     */
    /*       'discard_flag' correctly in NAL header, we will take a */
    /*       hit to optimize it.                                    */

    /* At present this routine also performs the following          */
    /* 1. Refine DQID list based on reference layer DQID            */
    /* 2. Calculates the POC for the target layer                   */

    {
        i4_status = isvcd_refine_dep_list(
            &ps_svcd_ctxt->s_vcl_nal, ps_svcd_ctxt->ps_sps, ps_svcd_ctxt->ps_subset_sps,
            ps_svcd_ctxt->ps_pps, &ps_svcd_ctxt->ai4_dq_id_map[0], &ps_svcd_ctxt->as_au_prms_dep[0],
            &ps_svcd_ctxt->as_pps_sps_prev_au[0], &ps_svcd_ctxt->i4_error_code, ps_svcd_ctxt);
    }

    if(0 != ps_svcd_ctxt->s_non_vcl_nal.i4_num_non_vcl_nals)
    {
        /* Decoding of non VCL data. In current implementation it  */
        /* decodes the followings:                                 */
        /*          1. Sequence parameter set                      */
        /*          2. Picture parameter set                       */
        /*          3. SEI message                                 */
        i4_non_vcl_status = isvcd_dec_non_vcl(&ps_svcd_ctxt->s_non_vcl_nal, ps_svcd_ctxt->ps_sps,
                                              ps_svcd_ctxt->ps_pps, ps_svcd_ctxt);

        if(OK != i4_non_vcl_status) return i4_non_vcl_status;
    }
    if(OK != i4_status) return (i4_status);
    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_video_decode                                      */
/*                                                                           */
/*  Description   :  handle video decode API command                         */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_video_decode(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    dec_struct_t *ps_dec;
    dec_struct_t *ps_dec_zero_lyr;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_zero_dec;

    svc_dec_ctxt_t *ps_svcd_ctxt;
    WORD32 i4_err_status = 0;

    UWORD32 bytes_consumed = 0;
    WORD32 ret = 0, api_ret_value = IV_SUCCESS;
    isvcd_video_decode_ip_t *ps_h264d_dec_ip;
    isvcd_video_decode_op_t *ps_h264d_dec_op;
    ivd_video_decode_ip_t *ps_dec_ip;
    ivd_video_decode_op_t *ps_dec_op;
    UWORD8 u1_res_id;

    ithread_set_name((void *) "Parse_thread");

    ps_svcd_ctxt = (svc_dec_ctxt_t *) (dec_hdl->pv_codec_handle);
    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[0];
    ps_dec = &ps_svc_lyr_dec->s_dec;

    ps_h264d_dec_ip = (isvcd_video_decode_ip_t *) pv_api_ip;
    ps_h264d_dec_op = (isvcd_video_decode_op_t *) pv_api_op;
    ps_dec_ip = &ps_h264d_dec_ip->s_ivd_video_decode_ip_t;
    ps_dec_op = &ps_h264d_dec_op->s_ivd_video_decode_op_t;

    {
        UWORD32 u4_size;
        u4_size = ps_dec_op->u4_size;
        memset(ps_h264d_dec_op, 0, sizeof(isvcd_video_decode_op_t));
        ps_dec_op->u4_size = u4_size;
    }

    ps_dec->pv_dec_out = ps_dec_op;
    if(ps_dec->init_done != 1)
    {
        return IV_FAIL;
    }

    /*Data memory barries instruction,so that bitstream write by the application
     * is complete*/
    DATA_SYNC();

    if(0 == ps_dec->u1_flushfrm)
    {
        if(ps_dec_ip->pv_stream_buffer == NULL)
        {
            ps_dec_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
            ps_dec_op->u4_error_code |= IVD_DEC_FRM_BS_BUF_NULL;
            return IV_FAIL;
        }
        if(ps_dec_ip->u4_num_Bytes <= 16)
        {
            ps_dec_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
            ps_dec_op->u4_error_code |= IVD_DEC_NUMBYTES_INV;
            return IV_FAIL;
        }
    }
#ifdef KEEP_THREADS_ACTIVE
    {
        UWORD32 i;
        ps_dec->i4_break_threads = 0;
        for(i = 0; i < 2; i++)
        {
            ret = ithread_mutex_lock(ps_dec->apv_proc_start_mutex[i]);
            RETURN_IF((ret != IV_SUCCESS), ret);

            ps_dec->ai4_process_start[i] = PROC_INIT;

            ret = ithread_mutex_unlock(ps_dec->apv_proc_start_mutex[i]);
            RETURN_IF((ret != IV_SUCCESS), ret);
        }
    }
#else
    ps_dec->u4_dec_thread_created = 0;
    ps_dec->u4_bs_deblk_thread_created = 0;
#endif
    ps_dec_op->u4_num_bytes_consumed = 0;
    ps_dec_op->i4_reorder_depth = -1;
    ps_dec_op->i4_display_index = DEFAULT_POC;

    ps_dec->ps_out_buffer = NULL;
    if(ps_dec_ip->u4_size >= offsetof(ivd_video_decode_ip_t, s_out_buffer))
        ps_dec->ps_out_buffer = &ps_dec_ip->s_out_buffer;

    if(0 == ps_dec->u4_share_disp_buf && ps_dec->i4_decode_header == 0)
    {
        UWORD32 i;
        if((ps_dec->ps_out_buffer->u4_num_bufs == 0) ||
           (ps_dec->ps_out_buffer->u4_num_bufs > IVD_VIDDEC_MAX_IO_BUFFERS))
        {
            ps_dec_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
            ps_dec_op->u4_error_code |= IVD_DISP_FRM_ZERO_OP_BUFS;
            return IV_FAIL;
        }

        for(i = 0; i < ps_dec->ps_out_buffer->u4_num_bufs; i++)
        {
            if(ps_dec->ps_out_buffer->pu1_bufs[i] == NULL)
            {
                ps_dec_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_dec_op->u4_error_code |= IVD_DISP_FRM_OP_BUF_NULL;
                return IV_FAIL;
            }

            if(ps_dec->ps_out_buffer->u4_min_out_buf_size[i] == 0)
            {
                ps_dec_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
                ps_dec_op->u4_error_code |= IVD_DISP_FRM_ZERO_OP_BUF_SIZE;
                return IV_FAIL;
            }
        }
    }

    if(ps_dec->u4_total_frames_decoded >= NUM_FRAMES_LIMIT)
    {
        ps_dec_op->u4_error_code = ERROR_FRAME_LIMIT_OVER;
        return IV_FAIL;
    }

    ps_dec_op->u4_error_code = 0;
    ps_dec_op->e_pic_type = IV_NA_FRAME;
    ps_dec_op->u4_output_present = 0;
    ps_dec_op->u4_frame_decoded_flag = 0;

    /* In case the decoder is not in flush mode(in shared mode),
     then decoder has to pick up a buffer to write current frame.
     Check if a frame is available in such cases */
    if(ps_dec->u1_init_dec_flag == 1 && ps_dec->u4_share_disp_buf == 1 && ps_dec->u1_flushfrm == 0)
    {
        UWORD32 i;
        WORD32 disp_avail = 0, free_id;

        /* Check if at least one buffer is available with the codec */
        /* If not then return to application with error */
        for(i = 0; i < ps_dec->u1_pic_bufs; i++)
        {
            if(0 == ps_dec->u4_disp_buf_mapping[i] || 1 == ps_dec->u4_disp_buf_to_be_freed[i])
            {
                disp_avail = 1;
                break;
            }
        }

        if(0 == disp_avail)
        {
            /* If something is queued for display wait for that buffer to be returned
             */

            ps_dec_op->u4_error_code = IVD_DEC_REF_BUF_NULL;
            ps_dec_op->u4_error_code |= (1 << IVD_UNSUPPORTEDPARAM);
            return (IV_FAIL);
        }

        while(1)
        {
            pic_buffer_t *ps_pic_buf;
            ps_pic_buf = (pic_buffer_t *) ih264_buf_mgr_get_next_free(
                (buf_mgr_t *) ps_dec->pv_pic_buf_mgr, &free_id);

            if(ps_pic_buf == NULL)
            {
                UWORD32 display_queued = 0;

                /* check if any buffer was given for display which is not returned yet */
                for(i = 0; i < (MAX_DISP_BUFS_NEW); i++)
                {
                    if(0 != ps_dec->u4_disp_buf_mapping[i])
                    {
                        display_queued = 1;
                        break;
                    }
                }
                /* If some buffer is queued for display, then codec has to singal an
                 error and wait for that buffer to be returned. If nothing is queued for
                 display then codec has ownership of all display buffers and it can
                 reuse any of the existing buffers and continue decoding */

                if(1 == display_queued)
                {
                    /* If something is queued for display wait for that buffer to be
                     * returned */
                    ps_dec_op->u4_error_code = IVD_DEC_REF_BUF_NULL;
                    ps_dec_op->u4_error_code |= (1 << IVD_UNSUPPORTEDPARAM);
                    return (IV_FAIL);
                }
            }
            else
            {
                /* If the buffer is with display, then mark it as in use and then look
                 * for a buffer again */
                if(1 == ps_dec->u4_disp_buf_mapping[free_id])
                {
                    ih264_buf_mgr_set_status((buf_mgr_t *) ps_dec->pv_pic_buf_mgr, free_id,
                                             BUF_MGR_IO);
                }
                else
                {
                    /**
                     *  Found a free buffer for present call. Release it now.
                     *  Will be again obtained later.
                     */
                    ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_pic_buf_mgr, free_id,
                                          BUF_MGR_IO);
                    break;
                }
            }
        }
    }

    if(ps_dec->u1_enable_mb_info && (ps_dec->i4_header_decoded & DECODED_SPS_MASK))
    {
        UWORD32 blk_qp_map_size = ps_h264d_dec_ip->u4_8x8_blk_qp_map_size;
        UWORD32 blk_type_map_size = ps_h264d_dec_ip->u4_8x8_blk_type_map_size;
        UWORD32 blk_8x8_map_size = ps_dec->u4_total_mbs << 2;
        if((ps_h264d_dec_ip->pu1_8x8_blk_qp_map && blk_qp_map_size < blk_8x8_map_size) ||
           (ps_h264d_dec_ip->pu1_8x8_blk_type_map && blk_type_map_size < blk_8x8_map_size))
        {
            ps_dec_op->u4_error_code |= 1 << IVD_UNSUPPORTEDPARAM;
            ps_dec_op->u4_error_code |= IH264D_INSUFFICIENT_METADATA_BUFFER;
            return IV_FAIL;
        }
    }

    if(ps_dec->u1_flushfrm && (1 == ps_svcd_ctxt->u1_pre_parse_in_flush))
    {
        if(ps_dec->u1_init_dec_flag == 0)
        {
            ps_dec->u1_flushfrm = 0;
            return (IV_FAIL);
        }

        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->s_vcl_nal.i4_num_res_lyrs - 1];
        ps_dec = &ps_svc_lyr_dec->s_dec;
        ps_dec->u4_fmt_conv_cur_row = 0;
        ps_dec->u4_output_present = 0;
        ps_dec->s_disp_op.u4_error_code = 1;

        ps_dec->ps_out_buffer = NULL;
        if(ps_dec_ip->u4_size >= offsetof(ivd_video_decode_ip_t, s_out_buffer))
        {
            ps_dec->ps_out_buffer = &ps_dec_ip->s_out_buffer;
        }
        ih264d_get_next_display_field(ps_dec, ps_dec->ps_out_buffer, &(ps_dec->s_disp_op));
        if(0 == ps_dec->s_disp_op.u4_error_code)
        {
            /* check output buffer size given by the application */
            if(check_app_out_buf_size(ps_dec) != IV_SUCCESS)
            {
                ps_dec_op->u4_error_code = IVD_DISP_FRM_ZERO_OP_BUF_SIZE;
                return (IV_FAIL);
            }

            ps_dec->u4_fmt_conv_cur_row = 0;
            ps_dec->u4_fmt_conv_num_rows = ps_dec->s_disp_frame_info.u4_y_ht;
            ih264d_format_convert(ps_dec, &(ps_dec->s_disp_op), ps_dec->u4_fmt_conv_cur_row,
                                  ps_dec->u4_fmt_conv_num_rows);
            ps_dec->u4_fmt_conv_cur_row += ps_dec->u4_fmt_conv_num_rows;
            ps_dec->u4_output_present = 1;
            if(ps_dec->u1_enable_mb_info)
            {
                UWORD32 disp_buf_id = ps_dec->s_disp_op.u4_disp_buf_id;
                if(ps_h264d_dec_ip->pu1_8x8_blk_qp_map)
                {
                    ps_h264d_dec_op->pu1_8x8_blk_qp_map = ps_h264d_dec_ip->pu1_8x8_blk_qp_map;
                    ps_h264d_dec_op->u4_8x8_blk_qp_map_size = ps_dec->u4_total_mbs << 2;
                    ih264_memcpy(ps_h264d_dec_op->pu1_8x8_blk_qp_map,
                                 ps_dec->as_buf_id_info_map[disp_buf_id].pu1_qp_map,
                                 ps_dec->u4_total_mbs << 2);
                }
                if(ps_h264d_dec_ip->pu1_8x8_blk_type_map)
                {
                    ps_h264d_dec_op->pu1_8x8_blk_type_map = ps_h264d_dec_ip->pu1_8x8_blk_type_map;
                    ps_h264d_dec_op->u4_8x8_blk_type_map_size = ps_dec->u4_total_mbs << 2;
                    ih264_memcpy(ps_h264d_dec_op->pu1_8x8_blk_type_map,
                                 ps_dec->as_buf_id_info_map[disp_buf_id].pu1_mb_type_map,
                                 ps_dec->u4_total_mbs << 2);
                }
            }
        }
        ih264d_export_sei_params(&ps_dec_op->s_sei_decode_op, ps_dec);

        ih264d_release_display_field(ps_dec, &(ps_dec->s_disp_op));

        ps_dec_op->u4_pic_wd = (UWORD32) ps_dec->u2_disp_width;
        ps_dec_op->u4_pic_ht = (UWORD32) ps_dec->u2_disp_height;
        ps_dec_op->i4_reorder_depth = ps_dec->i4_reorder_depth;
        ps_dec_op->i4_display_index = ps_dec->i4_display_index;
        ps_dec_op->u4_new_seq = 0;

        ps_dec_op->u4_output_present = ps_dec->u4_output_present;
        ps_dec_op->u4_progressive_frame_flag = ps_dec->s_disp_op.u4_progressive_frame_flag;
        ps_dec_op->e_output_format = ps_dec->s_disp_op.e_output_format;
        ps_dec_op->s_disp_frm_buf = ps_dec->s_disp_op.s_disp_frm_buf;
        ps_dec_op->e4_fld_type = ps_dec->s_disp_op.e4_fld_type;
        ps_dec_op->u4_ts = ps_dec->s_disp_op.u4_ts;
        ps_dec_op->u4_disp_buf_id = ps_dec->s_disp_op.u4_disp_buf_id;

        /*In the case of flush ,since no frame is decoded set pic type as invalid*/
        ps_dec_op->u4_is_ref_flag = UINT32_MAX;
        ps_dec_op->e_pic_type = IV_NA_FRAME;
        ps_dec_op->u4_frame_decoded_flag = 0;

        if(0 == ps_dec->s_disp_op.u4_error_code)
        {
            return (IV_SUCCESS);
        }
        else
            return (IV_FAIL);
    }

    if(ps_dec->u1_res_changed == 1)
    {
        /*if resolution has changed and all buffers have been flushed, reset
         * decoder*/
        if(((buf_mgr_t *) ps_dec->pv_pic_buf_mgr)->pv_mutex != NULL)
            ih264_buf_mgr_free(ps_dec->pv_pic_buf_mgr);
        if(((buf_mgr_t *) ps_dec->pv_mv_buf_mgr)->pv_mutex != NULL)
            ih264_buf_mgr_free(ps_dec->pv_mv_buf_mgr);

        isvcd_init_decoder(ps_svc_lyr_dec);
    }

    DEBUG_THREADS_PRINTF(" Starting process call\n");

    {
        vcl_node_t *ps_cur_node;
        UWORD8 u1_num_res_lyrs;
        vcl_buf_hdr_t *ps_vcl_buf;
        UWORD8 flush_decode = 1;
        ps_svcd_ctxt->u1_pre_parse_in_flush = 0;

        ret = isvcd_pre_parse_refine_au(ps_svcd_ctxt, ps_dec_ip, &ps_dec_op->u4_num_bytes_consumed);
        ps_svcd_ctxt->u1_pre_parse_in_flush = (ret == FLUSH);

        if(ret != OK)
        {
            UWORD32 error = ih264d_map_error((UWORD32) ret);
            if(ret != NOT_OK)
            {
                ps_dec_op->u4_error_code = error | ret;
            }
            if((ps_dec_op->u4_error_code >> IVD_FATALERROR) & 1)
            {
                ps_svcd_ctxt->u1_exit_till_next_IDR = 1;
            }
            api_ret_value = IV_FAIL;
            if((ret == IVD_RES_CHANGED) || (ret == IVD_MEM_ALLOC_FAILED) ||
               (ret == ERROR_UNAVAIL_PICBUF_T) || (ret == ERROR_UNAVAIL_MVBUF_T) ||
               (ret == ERROR_INV_SPS_PPS_T) || (ret == ERROR_FEATURE_UNAVAIL) ||
               (ret == IVD_STREAM_WIDTH_HEIGHT_NOT_SUPPORTED) ||
               (ret == IVD_DISP_FRM_ZERO_OP_BUF_SIZE))
            {
                ps_dec->u4_slice_start_code_found = 0;
            }
            if((ret == ERROR_INCOMPLETE_FRAME) || (ret == ERROR_DANGLING_FIELD_IN_PIC))
            {
                api_ret_value = IV_FAIL;
            }

            if(ret == ERROR_IN_LAST_SLICE_OF_PIC)
            {
                api_ret_value = IV_FAIL;
            }
        }

        if(NOT_OK == ret)
        {
            if(ps_dec->u4_pic_buf_got == 0)
            {
                ps_dec->i4_error_code = ERROR_START_CODE_NOT_FOUND;
                ps_dec_op->u4_error_code |= 1 << IVD_INSUFFICIENTDATA;

                isvcd_fill_output_struct_from_context(ps_svc_lyr_dec, ps_dec_op);

                ps_dec_op->u4_error_code = ps_dec->i4_error_code;
                ps_dec_op->u4_frame_decoded_flag = 0;
                return (IV_FAIL);
            }
            return (IV_SUCCESS);
        }

        u1_num_res_lyrs = ps_svcd_ctxt->s_vcl_nal.i4_num_res_lyrs;

        /* error concelment: exit till next IDR if any of Non Target layers are
         * corrupted */
        {
            ps_cur_node = ps_svcd_ctxt->s_vcl_nal.ps_bot_node;

            if(NULL != ps_cur_node)
            {
                if(!ps_cur_node->i4_idr_pic_flag)
                {
                    if(u1_num_res_lyrs != ps_svcd_ctxt->u1_prev_num_res_layers)
                    {
                        ps_svcd_ctxt->u1_exit_till_next_IDR = 1;
                        ps_dec_op->u4_error_code = ERROR_UNKNOWN_NAL;
                        return IV_FAIL;
                    }
                }
                else
                {
                    if(u1_num_res_lyrs != ps_svcd_ctxt->u1_prev_num_res_layers)
                    {
                        ps_svcd_ctxt->u1_prev_num_res_layers = u1_num_res_lyrs;
                    }
                }
            }
        }
        if(ps_svcd_ctxt->u1_prev_num_res_layers != u1_num_res_lyrs && (u1_num_res_lyrs != 0))
        {
            ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr + u1_num_res_lyrs - 1;
            ps_dec = &ps_svc_lyr_dec->s_dec;

            if(ps_dec->u1_init_dec_flag == 1)
            {
                ih264d_release_pics_in_dpb((void *) ps_dec, ps_dec->u1_pic_bufs);
                ih264d_release_display_bufs(ps_dec);
                ih264_disp_mgr_init((disp_mgr_t *) ps_dec->pv_disp_buf_mgr);

                ih264_buf_mgr_reset(ps_dec->pv_pic_buf_mgr);
                ih264_buf_mgr_reset(ps_dec->pv_mv_buf_mgr);
                ih264d_init_ref_bufs(ps_dec->ps_dpb_mgr);
            }

            // ps_svcd_ctxt->u1_prev_num_res_layers = u1_num_res_lyrs;
        }
        ps_svcd_ctxt->u1_parse_nal_unit_error = 0;

        if((1 == ps_svcd_ctxt->u1_exit_till_next_IDR) &&
           (ps_svcd_ctxt->s_vcl_nal.ps_bot_node != NULL))
        {
            if(1 == ps_svcd_ctxt->s_vcl_nal.ps_bot_node->i4_idr_pic_flag)
            {
                ps_svcd_ctxt->u1_exit_till_next_IDR = 0;

                for(u1_res_id = 0; u1_res_id < u1_num_res_lyrs; u1_res_id++)
                {
                    ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr + u1_res_id;
                    ps_dec = &ps_svc_lyr_dec->s_dec;
                    ih264_buf_mgr_reset(ps_dec->pv_pic_buf_mgr);
                    ih264_buf_mgr_reset(ps_dec->pv_mv_buf_mgr);
                }
            }
            else
            {
                ps_dec_op->u4_error_code = ERROR_UNKNOWN_NAL;
                return IV_FAIL;
            }
        }

        if((0 == ps_dec->i4_decode_header) && (OK == ret))
        {
            flush_decode = 0;
            ps_cur_node = ps_svcd_ctxt->s_vcl_nal.ps_bot_node;
            ps_svc_lyr_zero_dec = ps_svcd_ctxt->ps_svc_dec_lyr;
            ps_dec_zero_lyr = &ps_svc_lyr_zero_dec->s_dec;
            /* master loop */

            for(u1_res_id = 0; u1_res_id < u1_num_res_lyrs; u1_res_id++)
            {
                UWORD8 u1_layer_nal_data_present = 0;
                ps_svcd_ctxt->u1_cur_layer_id = u1_res_id;
                ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr + u1_res_id;
                ps_svc_lyr_dec->u1_res_init_done = 0;
                ps_svc_lyr_dec->u1_first_mb_addr_check = 1;
                ps_dec = &ps_svc_lyr_dec->s_dec;

                ps_dec->i4_decode_header = ps_dec_zero_lyr->i4_decode_header;
                ps_dec->i4_header_decoded = ps_dec_zero_lyr->i4_header_decoded;
                ps_dec->u1_pic_decode_done = 0;
                ps_dec->u4_fmt_conv_cur_row = 0;

                ps_dec->u4_output_present = 0;
                ps_dec->s_disp_op.u4_error_code = 1;
                ps_dec->u4_fmt_conv_num_rows = FMT_CONV_NUM_ROWS;
                ps_dec->u4_ts = ps_dec_ip->u4_ts;
                ps_dec->i4_frametype = IV_NA_FRAME;
                ps_dec->i4_content_type = IV_CONTENTTYPE_NA;

                ps_dec->u4_slice_start_code_found = 0;
                ps_dec->u2_cur_mb_addr = 0;
                ps_dec->u2_total_mbs_coded = 0;
                ps_dec->u2_cur_slice_num = 0;
                ps_dec->cur_dec_mb_num = 0;
                ps_dec->cur_recon_mb_num = 0;
                ps_dec->u4_first_slice_in_pic = 1;
                ps_dec->u1_slice_header_done = 0;
                ps_dec->u1_dangling_field = 0;

                ps_dec->u4_dec_thread_created = 0;
                ps_dec->u4_bs_deblk_thread_created = 0;
                ps_dec->u4_cur_bs_mb_num = 0;
                ps_dec->u4_cur_deblk_mb_num = 0;
                ps_dec->u4_start_recon_deblk = 0;
                ps_dec->u4_sps_cnt_in_process = 0;
                ps_dec->u4_pic_buf_got = 0;
                ps_dec->pv_dec_out = ps_dec_op;

                if(ps_dec_ip->u4_size >= offsetof(ivd_video_decode_ip_t, s_out_buffer))
                    ps_dec->ps_out_buffer = &ps_dec_ip->s_out_buffer;

                ps_dec->u1_nal_unit_type = ps_cur_node->i4_nal_unit_type;
                ps_dec->u1_separate_parse = 0;
                if(u1_res_id == (u1_num_res_lyrs - 1))
                {
                    ps_svc_lyr_dec->u1_layer_identifier = TARGET_LAYER;
                    if(ps_dec->u4_num_cores >= 2)
                    {
                        ps_dec->u4_num_cores = 2;
                        ps_dec->u1_separate_parse = 1;
                    }
                }
                else if(u1_res_id == 0)
                {
                    ps_svc_lyr_dec->u1_layer_identifier = BASE_LAYER;
                    ps_dec->u1_separate_parse = 0;
                    ps_dec->u4_num_cores = 1;
                }
                else if(u1_res_id != 0)
                {
                    ps_svc_lyr_dec->u1_layer_identifier = MEDIAL_ENHANCEMENT_LAYER;
                    ps_dec->u1_separate_parse = 0;
                    ps_dec->u4_num_cores = 1;
                }
                else
                {
                    return IV_FAIL;
                }

                ps_svc_lyr_dec->u1_base_res_flag = (0 == u1_res_id);
                ps_svc_lyr_dec->ps_nal_svc_ext->u1_idr_flag = ps_cur_node->i4_idr_pic_flag;
                ps_svc_lyr_dec->ps_nal_svc_ext->u1_dependency_id = ps_cur_node->i4_dependency_id;
                ps_svc_lyr_dec->ps_nal_svc_ext->u1_priority_id = ps_cur_node->i4_priority_id;
                ps_svc_lyr_dec->ps_nal_svc_ext->u1_no_inter_layer_pred_flag =
                    ps_cur_node->u1_acc_no_int_pred;

                ps_svc_lyr_dec->ps_nal_svc_ext->u1_quality_id = ps_cur_node->i4_quality_id;
                ps_svc_lyr_dec->ps_nal_svc_ext->u1_temporal_id = ps_cur_node->i4_temporal_id;

                ps_svc_lyr_dec->ps_nal_svc_ext->u1_use_ref_base_pic_flag =
                    ps_cur_node->i4_use_ref_base;
                ps_svc_lyr_dec->ps_nal_svc_ext->u1_discardable_flag = 0;
                ps_svc_lyr_dec->ps_nal_svc_ext->u1_svc_ext_flag = (u1_res_id > 1);
                ps_svc_lyr_dec->u4_pps_id_for_layer = UINT32_MAX;
                ps_vcl_buf = ps_cur_node->ps_first_vcl_nal;
                ps_svc_lyr_dec->u1_error_in_cur_frame = 0;

                /* Only for Non target Layers*/
                if(NULL != ps_cur_node->ps_top_node)
                {
                    ps_svc_lyr_dec->u1_inter_lyr_disable_dblk_filter_idc =
                        ps_cur_node->ps_top_node->i4_inter_lyr_dblk_idc;
                    ps_svc_lyr_dec->i1_inter_lyr_slice_alpha_c0_offset =
                        ps_cur_node->ps_top_node->i4_inter_lyr_alpha_c0_offset;
                    ps_svc_lyr_dec->i1_inter_lyr_slice_beta_offset =
                        ps_cur_node->ps_top_node->i4_inter_lyr_beta_offset;
                }

                while(NULL != ps_vcl_buf)
                {
                    u1_layer_nal_data_present = 1;
                    ps_dec->ps_bitstrm->u4_ofst = 0;
                    ps_dec->ps_bitstrm->pu4_buffer =
                        (UWORD32 *) ((UWORD8 *) ps_vcl_buf + ps_vcl_buf->i4_buf_offset +
                                     ps_vcl_buf->i4_slice_offset);

                    ps_dec->ps_bitstrm->u4_max_ofst = ps_vcl_buf->u4_max_bits;

                    ps_dec_op->u4_frame_decoded_flag = 0;
                    ret = isvcd_parse_nal_unit(ps_svc_lyr_dec, ps_cur_node->i4_nal_ref_idc);
                    if(ret != OK)
                    {
                        ps_svcd_ctxt->u1_parse_nal_unit_error = 1;
                        break;
                    }

                    /* go to the next slice */
                    ps_vcl_buf = ps_vcl_buf->ps_next;
                }
                /* error concelment: exit till next IDR if a Layer data is missing */
                if(0 == u1_layer_nal_data_present)
                {
                    ps_svcd_ctxt->u1_exit_till_next_IDR = 1;
                    ps_dec_op->u4_error_code = ERROR_UNKNOWN_NAL;
                    return IV_FAIL;
                }
                /* error concelment: exit till next IDR if any of Non Target layers are
                 * corrupted */
                if((ret != OK) && (u1_res_id != (u1_num_res_lyrs - 1)))
                {
                    ps_svcd_ctxt->u1_exit_till_next_IDR = 1;
                    ps_dec_op->u4_error_code = ERROR_UNKNOWN_NAL;
                    return IV_FAIL;
                }

                if((ret != OK) && (u1_res_id == (u1_num_res_lyrs - 1)))
                {
                    ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr + u1_num_res_lyrs - 1;
                    ps_dec = &ps_svc_lyr_dec->s_dec;

                    if((0 == ps_svcd_ctxt->u4_num_sps_ctr) || (0 == ps_svcd_ctxt->u4_num_pps_ctr) ||
                       (NULL == ps_dec->ps_cur_pps) || (ps_svc_lyr_dec->u1_res_init_done == 0))
                    {
                        ps_svcd_ctxt->u1_exit_till_next_IDR = 1;
                        ps_dec_op->u4_error_code = ERROR_UNKNOWN_NAL;
                        ih264d_signal_decode_thread(ps_dec);
                        return IV_FAIL;
                    }
                }
                ps_cur_node = ps_cur_node->ps_top_node;

                if((ps_dec->u4_pic_buf_got == 1) && (ret != IVD_MEM_ALLOC_FAILED) &&
                   ps_dec->u2_total_mbs_coded < ps_dec->u2_frm_ht_in_mbs * ps_dec->u2_frm_wd_in_mbs)
                {
                    // last slice - missing/corruption
                    WORD32 num_mb_skipped;
                    WORD32 prev_slice_err;
                    pocstruct_t temp_poc;
                    WORD32 ret1;
                    WORD32 ht_in_mbs;
                    ht_in_mbs = ps_dec->u2_pic_ht >> (4 + ps_dec->ps_cur_slice->u1_field_pic_flag);
                    num_mb_skipped =
                        (ht_in_mbs * ps_dec->u2_frm_wd_in_mbs) - ps_dec->u2_total_mbs_coded;

                    if(ps_dec->u4_first_slice_in_pic && (ps_dec->u4_pic_buf_got == 0))
                        prev_slice_err = 1;
                    else
                        prev_slice_err = 2;

                    if(ps_dec->u2_total_mbs_coded == 0)
                    {
                        prev_slice_err = 1;
                    }
                    ret1 = isvcd_mark_err_slice_skip(
                        ps_svc_lyr_dec, num_mb_skipped, ps_dec->u1_nal_unit_type == IDR_SLICE_NAL,
                        ps_dec->ps_cur_slice->u2_frame_num, &temp_poc, prev_slice_err);

                    if((ret1 == ERROR_UNAVAIL_PICBUF_T) || (ret1 == ERROR_UNAVAIL_MVBUF_T) ||
                       (ret1 == ERROR_INV_SPS_PPS_T) || (ret1 == ERROR_CORRUPTED_SLICE) ||
                       (ret == NOT_OK))
                    {
                        ret = ret1;
                    }
                }

                if((ret == IVD_RES_CHANGED) || (ret == IVD_MEM_ALLOC_FAILED) ||
                   (ret == ERROR_UNAVAIL_PICBUF_T) || (ret == ERROR_UNAVAIL_MVBUF_T) ||
                   (ret == ERROR_INV_SPS_PPS_T) || (ret == ERROR_CORRUPTED_SLICE) ||
                   (ret == IVD_DISP_FRM_ZERO_OP_BUF_SIZE) || (ret == NOT_OK))
                {
                    ps_svcd_ctxt->u1_exit_till_next_IDR = 1;
                    /* signal the decode thread */
                    ih264d_signal_decode_thread(ps_dec);
                    /* dont consume bitstream for change in resolution case */
                    if(ret == IVD_RES_CHANGED)
                    {
                        ps_dec_op->u4_num_bytes_consumed -= bytes_consumed;
                    }
                    return IV_FAIL;
                }

                /* Multi thread - for target Layer decoding*/
                if((ps_dec->u1_separate_parse) &&
                   (ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER) &&
                   (0 == ps_svc_lyr_dec->u1_error_in_cur_frame))
                {
                    /* If Format conversion is not complete,
                     complete it here */
                    if(ps_dec->u4_num_cores == 2)
                    {
                        /*do deblocking of all mbs*/
                        if((ps_dec->u4_nmb_deblk == 0) && (ps_dec->u4_start_recon_deblk == 1) &&
                           (ps_dec->ps_cur_sps->u1_mb_aff_flag == 0))
                        {
                            UWORD8 u1_end_of_row = 0;
                            UWORD32 u4_max_addr;
                            tfr_ctxt_t s_tfr_ctxt = {0};
                            tfr_ctxt_t *ps_tfr_cxt = &s_tfr_ctxt;
                            pad_mgr_t *ps_pad_mgr = &ps_dec->s_pad_mgr;
                            UWORD32 u4_slice_end = 0;

                            /*BS is done for all mbs while parsing*/
                            u4_max_addr = (ps_dec->u2_frm_wd_in_mbs * ps_dec->u2_frm_ht_in_mbs) - 1;
                            /* BS is moved post recon gen in SVC ext*/

                            ih264d_init_deblk_tfr_ctxt(ps_dec, ps_pad_mgr, ps_tfr_cxt,
                                                       ps_dec->u2_frm_wd_in_mbs, 0);

                            {
                                while(u4_slice_end != 1)
                                {
                                    dec_mb_info_t *p_cur_mb;
                                    WORD32 i, bs_mb_grp;
                                    bs_mb_grp = ps_dec->cur_dec_mb_num - ps_dec->u4_cur_bs_mb_num;

                                    for(i = 0; i < bs_mb_grp; i++)
                                    {
                                        p_cur_mb =
                                            &ps_dec->ps_frm_mb_info[ps_dec->u4_cur_bs_mb_num];

                                        DEBUG_THREADS_PRINTF("ps_dec->u4_cur_bs_mb_num = %d\n",
                                                             ps_dec->u4_cur_bs_mb_num);
                                        isvcd_compute_bs_non_mbaff_thread(ps_svc_lyr_dec, p_cur_mb,
                                                                          ps_dec->u4_cur_bs_mb_num);

                                        ps_dec->u4_cur_bs_mb_num++;
                                        ps_dec->u4_bs_cur_slice_num_mbs++;
                                    }
                                    if(ps_dec->u4_cur_bs_mb_num > u4_max_addr)
                                    {
                                        u4_slice_end = 1;
                                        u1_end_of_row = 1;
                                    }
                                    /*deblock MB group*/
                                    {
                                        UWORD32 u4_num_mbs;

                                        if(ps_dec->u4_cur_bs_mb_num > ps_dec->u4_cur_deblk_mb_num)
                                        {
                                            if(u1_end_of_row)
                                            {
                                                u4_num_mbs = ps_dec->u4_cur_bs_mb_num -
                                                             ps_dec->u4_cur_deblk_mb_num;
                                            }
                                            else
                                            {
                                                u4_num_mbs = ps_dec->u4_cur_bs_mb_num -
                                                             ps_dec->u4_cur_deblk_mb_num - 1;
                                            }
                                        }
                                        else
                                            u4_num_mbs = 0;

                                        ih264d_check_mb_map_deblk(ps_dec, u4_num_mbs, ps_tfr_cxt,
                                                                  0);
                                    }
                                }
                            }
                        }
                    }

                    /*signal the decode thread*/
                    ih264d_signal_decode_thread(ps_dec);
                }
                else if((ps_dec->u1_separate_parse) &&
                        (ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER))
                {
                    /*signal the decode thread*/
                    ih264d_signal_decode_thread(ps_dec);
                }

                DATA_SYNC();

                if((ps_dec_op->u4_error_code & 0xff) != ERROR_DYNAMIC_RESOLUTION_NOT_SUPPORTED)
                {
                    ps_dec_op->u4_pic_wd = (UWORD32) ps_dec->u2_disp_width;
                    ps_dec_op->u4_pic_ht = (UWORD32) ps_dec->u2_disp_height;
                    ps_dec_op->i4_reorder_depth = ps_dec->i4_reorder_depth;
                }

                // Report if header (sps and pps) has not been decoded yet
                if(ps_dec->i4_decode_header == 1 && ps_dec->i4_header_decoded != 3)
                {
                    ps_dec_op->u4_error_code |= (1 << IVD_INSUFFICIENTDATA);
                    api_ret_value = IV_FAIL;
                }

                if((ps_dec->u4_pic_buf_got == 1) && (ERROR_DANGLING_FIELD_IN_PIC != i4_err_status))
                {
                    /* For field pictures, set bottom and top picture decoded u4_flag correctly */

                    if(ps_dec->ps_cur_slice->u1_field_pic_flag)
                    {
                        if(1 == ps_dec->ps_cur_slice->u1_bottom_field_flag)
                        {
                            ps_dec->u1_top_bottom_decoded |= BOT_FIELD_ONLY;
                        }
                        else
                        {
                            ps_dec->u1_top_bottom_decoded |= TOP_FIELD_ONLY;
                        }
                    }
                    else
                    {
                        ps_dec->u1_top_bottom_decoded = TOP_FIELD_ONLY | BOT_FIELD_ONLY;
                    }

                    /* if new frame in not found (if we are still getting slices from
                     * previous frame) ih264d_deblock_display is not called. Such frames
                     * will not be added to reference /display
                     */
                    if((ps_dec->ps_dec_err_status->u1_err_flag & REJECT_CUR_PIC) == 0)
                    {
                        /* Calling Function to deblock Picture and Display */
                        ret = ih264d_deblock_display(ps_dec);
                    }

                    /*set to complete ,as we dont support partial frame decode*/
                    if(ps_dec->i4_header_decoded == 3)
                    {
                        ps_dec->u2_total_mbs_coded = ps_dec->ps_cur_sps->u2_max_mb_addr + 1;
                    }

                    /*Update the i4_frametype at the end of picture*/
                    if(ps_dec->ps_cur_slice->u1_nal_unit_type == IDR_SLICE_NAL)
                    {
                        ps_dec->i4_frametype = IV_IDR_FRAME;
                    }
                    else if(ps_dec->i4_pic_type == B_SLICE)
                    {
                        ps_dec->i4_frametype = IV_B_FRAME;
                    }
                    else if(ps_dec->i4_pic_type == P_SLICE)
                    {
                        ps_dec->i4_frametype = IV_P_FRAME;
                    }
                    else if(ps_dec->i4_pic_type == I_SLICE)
                    {
                        ps_dec->i4_frametype = IV_I_FRAME;
                    }
                    else
                    {
                        H264_DEC_DEBUG_PRINT("Shouldn't come here\n");
                    }

                    // Update the content type
                    ps_dec->i4_content_type = ps_dec->ps_cur_slice->u1_field_pic_flag;

                    ps_dec->u4_total_frames_decoded = ps_dec->u4_total_frames_decoded + 2;
                    ps_dec->u4_total_frames_decoded =
                        ps_dec->u4_total_frames_decoded - ps_dec->ps_cur_slice->u1_field_pic_flag;
                }

                /* In case the decoder is configured to run in low delay mode,
                 * then get display buffer and then format convert.
                 * Note in this mode, format conversion does not run paralelly in a
                 * thread and adds to the codec cycles
                 */
                if((IVD_DECODE_FRAME_OUT == ps_dec->e_frm_out_mode) && ps_dec->u1_init_dec_flag)
                {
                    ih264d_get_next_display_field(ps_dec, ps_dec->ps_out_buffer,
                                                  &(ps_dec->s_disp_op));

                    if(0 == ps_dec->s_disp_op.u4_error_code)
                    {
                        ps_dec->u4_fmt_conv_cur_row = 0;
                        ps_dec->u4_output_present = 1;
                    }
                    else
                    {
                        ps_dec->u4_output_present = 0;
                    }
                }

                isvcd_fill_output_struct_from_context(ps_svc_lyr_dec, ps_dec_op);

                /* If Format conversion is not complete,
                 complete it here */
                /* For Non -target Layers , Buffers are retrived but not displayed*/

                if((ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER) &&
                   ps_dec->u4_output_present &&
                   (ps_dec->u4_fmt_conv_cur_row < ps_dec->s_disp_frame_info.u4_y_ht))
                {
                    ps_dec->u4_fmt_conv_num_rows =
                        ps_dec->s_disp_frame_info.u4_y_ht - ps_dec->u4_fmt_conv_cur_row;
                    ih264d_format_convert(ps_dec, &(ps_dec->s_disp_op), ps_dec->u4_fmt_conv_cur_row,
                                          ps_dec->u4_fmt_conv_num_rows);
                    ps_dec->u4_fmt_conv_cur_row += ps_dec->u4_fmt_conv_num_rows;
                }

                ih264d_release_display_field(ps_dec, &(ps_dec->s_disp_op));

                if(ps_dec->i4_decode_header == 1 && (ps_dec->i4_header_decoded & 1) == 1)
                {
                    ps_dec_op->u4_progressive_frame_flag = 1;
                    if((NULL != ps_dec->ps_cur_sps) && (1 == (ps_dec->ps_cur_sps->u1_is_valid)))
                    {
                        if((0 == ps_dec->ps_sps->u1_frame_mbs_only_flag) &&
                           (0 == ps_dec->ps_sps->u1_mb_aff_flag))
                            ps_dec_op->u4_progressive_frame_flag = 0;
                    }
                }

                if((TOP_FIELD_ONLY | BOT_FIELD_ONLY) == ps_dec->u1_top_bottom_decoded)
                {
                    ps_dec->u1_top_bottom_decoded = 0;
                }
                /*--------------------------------------------------------------------*/
                /* Do End of Pic processing.                                          */
                /* Should be called only if frame was decoded in previous process call*/
                /*--------------------------------------------------------------------*/
                if(ps_dec->u4_pic_buf_got == 1)
                {
                    if(1 == ps_dec->u1_last_pic_not_decoded)
                    {
                        ret = ih264d_end_of_pic_dispbuf_mgr(ps_dec);

                        if(ret != OK) return ret;

                        ret = ih264d_end_of_pic(ps_dec);
                        if(ret != OK) return ret;
                    }
                    else
                    {
                        ret = ih264d_end_of_pic(ps_dec);
                        if(ret != OK) return ret;
                    }
                }

                if(ps_dec->u1_enable_mb_info && ps_dec->u4_output_present)
                {
                    UWORD32 disp_buf_id = ps_dec->s_disp_op.u4_disp_buf_id;
                    if(ps_h264d_dec_ip->pu1_8x8_blk_qp_map)
                    {
                        ps_h264d_dec_op->pu1_8x8_blk_qp_map = ps_h264d_dec_ip->pu1_8x8_blk_qp_map;
                        ps_h264d_dec_op->u4_8x8_blk_qp_map_size = ps_dec->u4_total_mbs << 2;
                        ih264_memcpy(ps_h264d_dec_op->pu1_8x8_blk_qp_map,
                                     ps_dec->as_buf_id_info_map[disp_buf_id].pu1_qp_map,
                                     ps_dec->u4_total_mbs << 2);
                    }
                    if(ps_h264d_dec_ip->pu1_8x8_blk_type_map)
                    {
                        ps_h264d_dec_op->pu1_8x8_blk_type_map =
                            ps_h264d_dec_ip->pu1_8x8_blk_type_map;
                        ps_h264d_dec_op->u4_8x8_blk_type_map_size = ps_dec->u4_total_mbs << 2;
                        ih264_memcpy(ps_h264d_dec_op->pu1_8x8_blk_type_map,
                                     ps_dec->as_buf_id_info_map[disp_buf_id].pu1_mb_type_map,
                                     ps_dec->u4_total_mbs << 2);
                    }
                }
                /*Data memory barrier instruction,so that yuv write by the library is
                 * complete*/
                DATA_SYNC();

                H264_DEC_DEBUG_PRINT("The num bytes consumed: %d\n",
                                     ps_dec_op->u4_num_bytes_consumed);
            }
        }
        /* highest layer for flush validation */

        if((ps_dec->u1_flushfrm) && (1 == flush_decode))
        {
            u1_res_id = u1_num_res_lyrs - 1;
            ps_svc_lyr_dec = ps_svcd_ctxt->ps_svc_dec_lyr + u1_res_id;
            ps_dec = &ps_svc_lyr_dec->s_dec;

            ih264d_get_next_display_field(ps_dec, ps_dec->ps_out_buffer, &(ps_dec->s_disp_op));
            if(0 == ps_dec->s_disp_op.u4_error_code)
            {
                /* check output buffer size given by the application */
                if(check_app_out_buf_size(ps_dec) != IV_SUCCESS)
                {
                    ps_dec_op->u4_error_code = IVD_DISP_FRM_ZERO_OP_BUF_SIZE;
                    return (IV_FAIL);
                }

                ps_dec->u4_fmt_conv_cur_row = 0;
                ps_dec->u4_fmt_conv_num_rows = ps_dec->s_disp_frame_info.u4_y_ht;
                ih264d_format_convert(ps_dec, &(ps_dec->s_disp_op), ps_dec->u4_fmt_conv_cur_row,
                                      ps_dec->u4_fmt_conv_num_rows);
                ps_dec->u4_fmt_conv_cur_row += ps_dec->u4_fmt_conv_num_rows;
                ps_dec->u4_output_present = 1;
            }
            else
            {
                ps_dec->u4_output_present = 0;
            }
            ih264d_export_sei_params(&ps_dec_op->s_sei_decode_op, ps_dec);

            ih264d_release_display_field(ps_dec, &(ps_dec->s_disp_op));

            ps_dec_op->u4_pic_wd = (UWORD32) ps_dec->u2_disp_width;
            ps_dec_op->u4_pic_ht = (UWORD32) ps_dec->u2_disp_height;
            ps_dec_op->i4_reorder_depth = ps_dec->i4_reorder_depth;
            ps_dec_op->i4_display_index = ps_dec->i4_display_index;

            ps_dec_op->u4_new_seq = 0;
            ps_dec_op->u4_output_present = (ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER)
                                               ? ps_dec->u4_output_present
                                               : 0;
            ps_dec_op->u4_progressive_frame_flag = ps_dec->s_disp_op.u4_progressive_frame_flag;
            ps_dec_op->e_output_format = ps_dec->s_disp_op.e_output_format;
            ps_dec_op->s_disp_frm_buf = ps_dec->s_disp_op.s_disp_frm_buf;
            ps_dec_op->e4_fld_type = ps_dec->s_disp_op.e4_fld_type;
            ps_dec_op->u4_ts = ps_dec->s_disp_op.u4_ts;
            ps_dec_op->u4_disp_buf_id = ps_dec->s_disp_op.u4_disp_buf_id;

            /*In the case of flush ,since no frame is decoded set pic type as invalid*/
            ps_dec_op->u4_is_ref_flag = UINT32_MAX;
            ps_dec_op->e_pic_type = IV_NA_FRAME;
            ps_dec_op->u4_frame_decoded_flag = 0;

            if(0 == ps_dec->s_disp_op.u4_error_code)
            {
                return (IV_SUCCESS);
            }
            else
                return (IV_FAIL);
        }
    }

    if((ps_dec_op->u4_error_code & 0xff) != ERROR_DYNAMIC_RESOLUTION_NOT_SUPPORTED)
    {
        ps_dec_op->u4_pic_wd = (UWORD32) ps_dec->u2_disp_width;
        ps_dec_op->u4_pic_ht = (UWORD32) ps_dec->u2_disp_height;
        ps_dec_op->i4_reorder_depth = ps_dec->i4_reorder_depth;
    }
    return api_ret_value;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_set_display_frame                                 */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_set_display_frame(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    UWORD32 u4_disp_buf_size[3] = {0};
    UWORD32 u4_num_disp_bufs;
    ivd_set_display_frame_ip_t *dec_disp_ip;
    ivd_set_display_frame_op_t *dec_disp_op;
    UWORD32 i;
    dec_struct_t *ps_dec;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;

    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;

    dec_disp_ip = (ivd_set_display_frame_ip_t *) pv_api_ip;
    dec_disp_op = (ivd_set_display_frame_op_t *) pv_api_op;
    dec_disp_op->u4_error_code = 0;

    ps_dec->u4_num_disp_bufs = 0;
    if(ps_dec->u4_share_disp_buf)
    {
        UWORD32 u4_num_bufs = dec_disp_ip->num_disp_bufs;

        u4_num_bufs = MIN(u4_num_bufs, MAX_DISP_BUFS_NEW);
        ps_dec->u4_num_disp_bufs = u4_num_bufs;

        /* Get the number and sizes of the first buffer. Compare this with the
         * rest to make sure all the buffers are of the same size.
         */
        u4_num_disp_bufs = dec_disp_ip->s_disp_buffer[0].u4_num_bufs;

        u4_disp_buf_size[0] = dec_disp_ip->s_disp_buffer[0].u4_min_out_buf_size[0];
        u4_disp_buf_size[1] = dec_disp_ip->s_disp_buffer[0].u4_min_out_buf_size[1];
        u4_disp_buf_size[2] = dec_disp_ip->s_disp_buffer[0].u4_min_out_buf_size[2];

        for(i = 0; i < u4_num_bufs; i++)
        {
            if(dec_disp_ip->s_disp_buffer[i].u4_num_bufs != u4_num_disp_bufs)
            {
                return IV_FAIL;
            }

            if((dec_disp_ip->s_disp_buffer[i].u4_min_out_buf_size[0] != u4_disp_buf_size[0]) ||
               (dec_disp_ip->s_disp_buffer[i].u4_min_out_buf_size[1] != u4_disp_buf_size[1]) ||
               (dec_disp_ip->s_disp_buffer[i].u4_min_out_buf_size[2] != u4_disp_buf_size[2]))
            {
                return IV_FAIL;
            }

            ps_dec->disp_bufs[i].u4_num_bufs = dec_disp_ip->s_disp_buffer[i].u4_num_bufs;

            ps_dec->disp_bufs[i].buf[0] = dec_disp_ip->s_disp_buffer[i].pu1_bufs[0];
            ps_dec->disp_bufs[i].buf[1] = dec_disp_ip->s_disp_buffer[i].pu1_bufs[1];
            ps_dec->disp_bufs[i].buf[2] = dec_disp_ip->s_disp_buffer[i].pu1_bufs[2];

            ps_dec->disp_bufs[i].u4_bufsize[0] =
                dec_disp_ip->s_disp_buffer[i].u4_min_out_buf_size[0];
            ps_dec->disp_bufs[i].u4_bufsize[1] =
                dec_disp_ip->s_disp_buffer[i].u4_min_out_buf_size[1];
            ps_dec->disp_bufs[i].u4_bufsize[2] =
                dec_disp_ip->s_disp_buffer[i].u4_min_out_buf_size[2];
        }
    }
    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : ih264d_set_flush_mode_svt_ext                            */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Globals       : <Does it use any global variables?>                      */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_set_flush_mode(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    UWORD8 u1_layer_id;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    ivd_ctl_flush_op_t *ps_ctl_op = (ivd_ctl_flush_op_t *) pv_api_op;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ctl_op->u4_error_code = 0;

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[0];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    if(0 == ps_dec->i4_decode_header)
    {
        ps_svcd_ctxt->i4_eos_flag = 1;
    }

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_dec = &ps_svc_lyr_dec->s_dec;
        UNUSED(pv_api_ip);

        /* Signal flush frame control call */
        ps_dec->u1_flushfrm = 1;

        if(ps_dec->u1_init_dec_flag == 1)
        {
            ih264d_release_pics_in_dpb((void *) ps_dec, ps_dec->u1_pic_bufs);
            ih264d_release_display_bufs(ps_dec);
        }

        ps_ctl_op->u4_error_code = 0;

        /* Ignore dangling fields during flush */
        ps_dec->u1_top_bottom_decoded = 0;
    }

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_status                                         */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Globals       : <Does it use any global variables?>                      */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_get_status(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    UWORD32 i;
    dec_struct_t *ps_dec;
    UWORD32 pic_wd, pic_ht;
    ivd_ctl_getstatus_op_t *ps_ctl_op = (ivd_ctl_getstatus_op_t *) pv_api_op;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    UNUSED(pv_api_ip);
    ps_ctl_op->u4_error_code = 0;

    if((NULL != ps_dec->ps_cur_sps) && (1 == (ps_dec->ps_cur_sps->u1_is_valid)))
    {
        ps_ctl_op->u4_pic_ht = ps_dec->u2_disp_height;
        ps_ctl_op->u4_pic_wd = ps_dec->u2_disp_width;

        if(0 == ps_dec->u4_share_disp_buf)
        {
            pic_wd = ps_dec->u2_disp_width;
            pic_ht = ps_dec->u2_disp_height;
        }
        else
        {
            pic_wd = ps_dec->u2_frm_wd_y;
            pic_ht = ps_dec->u2_frm_ht_y;
        }
    }
    else
    {
        pic_wd = 0;
        pic_ht = 0;
        ps_ctl_op->u4_pic_ht = pic_wd;
        ps_ctl_op->u4_pic_wd = pic_ht;

        if(1 == ps_dec->u4_share_disp_buf)
        {
            pic_wd += (PAD_LEN_Y_H << 1);
            pic_ht += (PAD_LEN_Y_V << 2);
        }
    }

    if(ps_dec->u4_app_disp_width > pic_wd) pic_wd = ps_dec->u4_app_disp_width;
    if(0 == ps_dec->u4_share_disp_buf)
        ps_ctl_op->u4_num_disp_bufs = 1;
    else
    {
        if((NULL != ps_dec->ps_cur_sps) && (1 == (ps_dec->ps_cur_sps->u1_is_valid)))
        {
            if((ps_dec->ps_cur_sps->u1_vui_parameters_present_flag == 1) &&
               (1 == ps_dec->ps_cur_sps->s_vui.u1_bitstream_restriction_flag))
            {
                ps_ctl_op->u4_num_disp_bufs = ps_dec->ps_cur_sps->s_vui.u4_num_reorder_frames + 1;
            }
            else
            {
                /*if VUI is not present assume maximum possible refrence frames for the
                 * level, as max reorder frames*/
                ps_ctl_op->u4_num_disp_bufs = ih264d_get_dpb_size(ps_dec->ps_cur_sps);
            }

            ps_ctl_op->u4_num_disp_bufs += ps_dec->ps_cur_sps->u1_num_ref_frames + 1;
        }
        else
        {
            ps_ctl_op->u4_num_disp_bufs = 32;
        }
        ps_ctl_op->u4_num_disp_bufs = MAX(ps_ctl_op->u4_num_disp_bufs, 6);
        ps_ctl_op->u4_num_disp_bufs = MIN(ps_ctl_op->u4_num_disp_bufs, 32);
    }

    ps_ctl_op->u4_error_code = ps_dec->i4_error_code;
    ps_ctl_op->u4_frame_rate = 0;
    ps_ctl_op->u4_bit_rate = 0;
    ps_ctl_op->e_content_type = ps_dec->i4_content_type;
    ps_ctl_op->e_output_chroma_format = ps_dec->u1_chroma_format;
    ps_ctl_op->u4_min_num_in_bufs = MIN_IN_BUFS;

    if(ps_dec->u1_chroma_format == IV_YUV_420P)
    {
        ps_ctl_op->u4_min_num_out_bufs = MIN_OUT_BUFS_420;
    }
    else if(ps_dec->u1_chroma_format == IV_YUV_422ILE)
    {
        ps_ctl_op->u4_min_num_out_bufs = MIN_OUT_BUFS_422ILE;
    }
    else if(ps_dec->u1_chroma_format == IV_RGB_565)
    {
        ps_ctl_op->u4_min_num_out_bufs = MIN_OUT_BUFS_RGB565;
    }
    else if((ps_dec->u1_chroma_format == IV_YUV_420SP_UV) ||
            (ps_dec->u1_chroma_format == IV_YUV_420SP_VU))
    {
        ps_ctl_op->u4_min_num_out_bufs = MIN_OUT_BUFS_420SP;
    }
    else
    {
        // Invalid chroma format; Error code may be updated, verify in testing if needed
        ps_ctl_op->u4_error_code = ERROR_FEATURE_UNAVAIL;
        return IV_FAIL;
    }

    for(i = 0; i < ps_ctl_op->u4_min_num_in_bufs; i++)
    {
        ps_ctl_op->u4_min_in_buf_size[i] = MAX(256000, pic_wd * pic_ht * 3 / 2);
    }

    if(ps_dec->u1_chroma_format == IV_YUV_420P)
    {
        ps_ctl_op->u4_min_out_buf_size[0] = (pic_wd * pic_ht);
        ps_ctl_op->u4_min_out_buf_size[1] = (pic_wd * pic_ht) >> 2;
        ps_ctl_op->u4_min_out_buf_size[2] = (pic_wd * pic_ht) >> 2;
    }
    else if(ps_dec->u1_chroma_format == IV_YUV_422ILE)
    {
        ps_ctl_op->u4_min_out_buf_size[0] = (pic_wd * pic_ht) * 2;
        ps_ctl_op->u4_min_out_buf_size[1] = ps_ctl_op->u4_min_out_buf_size[2] = 0;
    }
    else if(ps_dec->u1_chroma_format == IV_RGB_565)
    {
        ps_ctl_op->u4_min_out_buf_size[0] = (pic_wd * pic_ht) * 2;
        ps_ctl_op->u4_min_out_buf_size[1] = ps_ctl_op->u4_min_out_buf_size[2] = 0;
    }
    else if((ps_dec->u1_chroma_format == IV_YUV_420SP_UV) ||
            (ps_dec->u1_chroma_format == IV_YUV_420SP_VU))
    {
        ps_ctl_op->u4_min_out_buf_size[0] = (pic_wd * pic_ht);
        ps_ctl_op->u4_min_out_buf_size[1] = (pic_wd * pic_ht) >> 1;
        ps_ctl_op->u4_min_out_buf_size[2] = 0;
    }

    ps_dec->u4_num_disp_bufs_requested = ps_ctl_op->u4_num_disp_bufs;
    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name :    isvcd_get_buf_info                                    */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Globals       : <Does it use any global variables?>                      */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_buf_info(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    UWORD8 i = 0;  // Default for 420P format
    UWORD16 pic_wd, pic_ht;
    ivd_ctl_getbufinfo_op_t *ps_ctl_op = (ivd_ctl_getbufinfo_op_t *) pv_api_op;
    UWORD32 au4_min_out_buf_size[IVD_VIDDEC_MAX_IO_BUFFERS] = {0};
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    UNUSED(pv_api_ip);

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    ps_ctl_op->u4_error_code = 0;

    ps_ctl_op->u4_min_num_in_bufs = MIN_IN_BUFS;
    ps_ctl_op->u4_num_disp_bufs = 1;
    pic_wd = 0;
    pic_ht = 0;

    if(ps_dec->i4_header_decoded == 3)
    {
        if(0 == ps_dec->u4_share_disp_buf)
        {
            pic_wd = ps_dec->u2_disp_width;
            pic_ht = ps_dec->u2_disp_height;
        }
        else
        {
            pic_wd = ps_dec->u2_frm_wd_y;
            pic_ht = ps_dec->u2_frm_ht_y;
        }
    }

    for(i = 0; i < ps_ctl_op->u4_min_num_in_bufs; i++)
    {
        ps_ctl_op->u4_min_in_buf_size[i] = MAX(256000, pic_wd * pic_ht * 3 / 2);
    }
    if((WORD32) ps_dec->u4_app_disp_width > pic_wd) pic_wd = ps_dec->u4_app_disp_width;

    if(0 == ps_dec->u4_share_disp_buf)
        ps_ctl_op->u4_num_disp_bufs = 1;
    else
    {
        if((NULL != ps_dec->ps_cur_sps) && (1 == (ps_dec->ps_cur_sps->u1_is_valid)))
        {
            if((ps_dec->ps_cur_sps->u1_vui_parameters_present_flag == 1) &&
               (1 == ps_dec->ps_cur_sps->s_vui.u1_bitstream_restriction_flag))
            {
                ps_ctl_op->u4_num_disp_bufs = ps_dec->ps_cur_sps->s_vui.u4_num_reorder_frames + 1;
            }
            else
            {
                /*if VUI is not present assume maximum possible refrence frames for the
                 * level, as max reorder frames*/
                ps_ctl_op->u4_num_disp_bufs = ih264d_get_dpb_size(ps_dec->ps_cur_sps);
            }

            ps_ctl_op->u4_num_disp_bufs += ps_dec->ps_cur_sps->u1_num_ref_frames + 1;
        }
        else
        {
            ps_ctl_op->u4_num_disp_bufs = 32;
        }

        ps_ctl_op->u4_num_disp_bufs = MAX(ps_ctl_op->u4_num_disp_bufs, 6);
        ps_ctl_op->u4_num_disp_bufs = MIN(ps_ctl_op->u4_num_disp_bufs, 32);
    }

    ps_ctl_op->u4_min_num_out_bufs =
        ih264d_get_outbuf_size(pic_wd, pic_ht, ps_dec->u1_chroma_format, &au4_min_out_buf_size[0]);

    for(i = 0; i < ps_ctl_op->u4_min_num_out_bufs; i++)
    {
        ps_ctl_op->u4_min_out_buf_size[i] = au4_min_out_buf_size[i];
    }

    ps_dec->u4_num_disp_bufs_requested = ps_ctl_op->u4_num_disp_bufs;

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_set_params                                         */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_set_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    WORD32 ret = IV_SUCCESS;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    WORD32 u1_layer_id;

    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;
    ps_svcd_ctxt->i4_eos_flag = 0;
    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        isvcd_ctl_set_config_ip_t *ps_h264d_ctl_ip = (isvcd_ctl_set_config_ip_t *) pv_api_ip;
        isvcd_ctl_set_config_op_t *ps_h264d_ctl_op = (isvcd_ctl_set_config_op_t *) pv_api_op;
        ivd_ctl_set_config_ip_t *ps_ctl_ip = &ps_h264d_ctl_ip->s_ivd_ctl_set_config_ip_t;
        ivd_ctl_set_config_op_t *ps_ctl_op = &ps_h264d_ctl_op->s_ivd_ctl_set_config_op_t;

        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_dec = &ps_svc_lyr_dec->s_dec;

        ps_dec->u1_flushfrm = 0;
        ps_dec->u4_skip_frm_mask = 0;
        ps_ctl_op->u4_error_code = 0;

        if(ps_ctl_ip->e_frm_skip_mode != IVD_SKIP_NONE)
        {
            ps_ctl_op->u4_error_code = (1 << IVD_UNSUPPORTEDPARAM);
            ret = IV_FAIL;
        }

        if(ps_ctl_ip->u4_disp_wd >= ps_dec->u2_disp_width)
        {
            ps_dec->u4_app_disp_width = ps_ctl_ip->u4_disp_wd;
        }
        else if(0 == ps_dec->i4_header_decoded)
        {
            ps_dec->u4_app_disp_width = ps_ctl_ip->u4_disp_wd;
        }
        else if(ps_ctl_ip->u4_disp_wd == 0)
        {
            ps_dec->u4_app_disp_width = 0;
        }
        else
        {
            /*
             * Set the display width to zero. This will ensure that the wrong value we
             * had stored (0xFFFFFFFF) does not propogate.
             */
            ps_dec->u4_app_disp_width = 0;
            ps_ctl_op->u4_error_code |= (1 << IVD_UNSUPPORTEDPARAM);
            ps_ctl_op->u4_error_code |= ERROR_DISP_WIDTH_INVALID;
            ret = IV_FAIL;
        }

        if(ps_ctl_ip->e_vid_dec_mode == IVD_DECODE_FRAME)
            ps_dec->i4_decode_header = 0;
        else if(ps_ctl_ip->e_vid_dec_mode == IVD_DECODE_HEADER)
            ps_dec->i4_decode_header = 1;
        else
        {
            ps_ctl_op->u4_error_code = (1 << IVD_UNSUPPORTEDPARAM);
            ps_dec->i4_decode_header = 1;
            ret = IV_FAIL;
        }
        ps_dec->e_frm_out_mode = IVD_DISPLAY_FRAME_OUT;

        if((ps_ctl_ip->e_frm_out_mode != IVD_DECODE_FRAME_OUT) &&
           (ps_ctl_ip->e_frm_out_mode != IVD_DISPLAY_FRAME_OUT))
        {
            ps_ctl_op->u4_error_code = (1 << IVD_UNSUPPORTEDPARAM);
            ret = IV_FAIL;
        }
        ps_dec->e_frm_out_mode = ps_ctl_ip->e_frm_out_mode;
    }
    return ret;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_set_target_layer                                   */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         05 04 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_set_target_layer(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    WORD32 ret = IV_SUCCESS;

    isvcd_set_target_layer_ip_t *ps_ip;
    isvcd_set_target_layer_op_t *ps_op;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ip = (isvcd_set_target_layer_ip_t *) pv_api_ip;
    ps_op = (isvcd_set_target_layer_op_t *) pv_api_op;

    ps_svcd_ctxt->u1_tgt_dep_id = ps_ip->u1_tgt_dep_id;
    ps_svcd_ctxt->u1_tgt_quality_id = ps_ip->u1_tgt_quality_id;
    ps_svcd_ctxt->u1_tgt_temp_id = ps_ip->u1_tgt_temp_id;
    ps_svcd_ctxt->u1_tgt_priority_id = ps_ip->u1_tgt_priority_id;

    ret = isvcd_nal_parse_set_target_attr(ps_ip->u1_tgt_quality_id, ps_ip->u1_tgt_dep_id,
                                          ps_ip->u1_tgt_temp_id, ps_ip->u1_tgt_priority_id,
                                          ps_svcd_ctxt->pv_nal_parse_ctxt);
    ps_op->u4_error_code = 0;

    return ret;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_set_default_params                                 */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Copied from set_params               */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_set_default_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    WORD32 ret = IV_SUCCESS;
    UWORD8 u1_layer_id;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ivd_ctl_set_config_op_t *ps_ctl_op = (ivd_ctl_set_config_op_t *) pv_api_op;

    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;
    UNUSED(pv_api_ip);

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_dec = &ps_svc_lyr_dec->s_dec;

        ps_dec->u4_app_disp_width = 0;
        ps_dec->u4_skip_frm_mask = 0;
        ps_dec->i4_decode_header = 1;
    }
    ps_ctl_op->u4_error_code = 0;

    return ret;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_delete                                            */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Globals       : <Does it use any global variables?>                      */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_delete(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    isvcd_delete_ip_t *ps_ip = (isvcd_delete_ip_t *) pv_api_ip;
    isvcd_delete_op_t *ps_op = (isvcd_delete_op_t *) pv_api_op;

    UWORD8 u1_layer_id;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;
    UNUSED(ps_ip);

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        isvcd_free_dynamic_bufs(ps_svc_lyr_dec);
    }
    isvcd_free_static_bufs(dec_hdl);
    ps_op->s_ivd_delete_op_t.u4_error_code = 0;

    return IV_SUCCESS;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_reset                                             */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Globals       : <Does it use any global variables?>                      */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_reset(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    ivd_ctl_reset_op_t *ps_ctl_op = (ivd_ctl_reset_op_t *) pv_api_op;
    UWORD8 u1_layer_id;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;
    UNUSED(pv_api_ip);
    ps_ctl_op->u4_error_code = 0;

    ps_svcd_ctxt->i4_eos_flag = 0;
    ps_svcd_ctxt->u4_num_sps_ctr = 0;
    ps_svcd_ctxt->u4_num_pps_ctr = 0;
    ps_svcd_ctxt->u1_pre_parse_in_flush = 1;
    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_dec = &ps_svc_lyr_dec->s_dec;
        if(ps_dec != NULL)
        {
            if(((buf_mgr_t *) ps_dec->pv_pic_buf_mgr)->pv_mutex != NULL)
                ih264_buf_mgr_free(ps_dec->pv_pic_buf_mgr);
            if(((buf_mgr_t *) ps_dec->pv_mv_buf_mgr)->pv_mutex != NULL)
                ih264_buf_mgr_free(ps_dec->pv_mv_buf_mgr);

            isvcd_init_decoder(ps_svc_lyr_dec);
            ps_dec->u1_flushfrm = 0;
        }
        else
        {
            H264_DEC_DEBUG_PRINT("\nReset called without Initializing the decoder\n");
            ps_ctl_op->u4_error_code = ERROR_INIT_NOT_DONE;
        }
    }
    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name :  isvcd_ctl                                               */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_ctl(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    ivd_ctl_set_config_ip_t *ps_ctl_ip;
    ivd_ctl_set_config_op_t *ps_ctl_op;
    WORD32 ret = IV_SUCCESS;
    UWORD32 subcommand;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;
    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    if(ps_dec->init_done != 1)
    {
        return IV_FAIL;
    }
    ps_ctl_ip = (ivd_ctl_set_config_ip_t *) pv_api_ip;
    ps_ctl_op = (ivd_ctl_set_config_op_t *) pv_api_op;
    ps_ctl_op->u4_error_code = 0;
    subcommand = ps_ctl_ip->e_sub_cmd;

    switch(subcommand)
    {
        case IVD_CMD_CTL_GETPARAMS:
            ret = isvcd_get_status(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IVD_CMD_CTL_SETPARAMS:
            ret = isvcd_set_params(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IVD_CMD_CTL_RESET:
            ret = isvcd_reset(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IVD_CMD_CTL_SETDEFAULT:
            ret = isvcd_set_default_params(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IVD_CMD_CTL_FLUSH:
            ret = isvcd_set_flush_mode(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IVD_CMD_CTL_GETBUFINFO:
            ret = isvcd_get_buf_info(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IVD_CMD_CTL_GETVERSION:
            ret = ih264d_get_version(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IH264D_CMD_CTL_DEGRADE:
            ret = isvcd_set_degrade(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;

        case IH264D_CMD_CTL_SET_NUM_CORES:
            ret = isvcd_set_num_cores(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IH264D_CMD_CTL_GET_BUFFER_DIMENSIONS:
            ret = isvcd_get_frame_dimensions(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IH264D_CMD_CTL_GET_VUI_PARAMS:
            ret = isvcd_get_vui_params(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IH264D_CMD_CTL_GET_SEI_MDCV_PARAMS:
            ret = isvcd_get_sei_mdcv_params(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IH264D_CMD_CTL_GET_SEI_CLL_PARAMS:
            ret = isvcd_get_sei_cll_params(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IH264D_CMD_CTL_GET_SEI_AVE_PARAMS:
            ret = isvcd_get_sei_ave_params(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IH264D_CMD_CTL_GET_SEI_CCV_PARAMS:
            ret = isvcd_get_sei_ccv_params(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IH264D_CMD_CTL_SET_PROCESSOR:
            ret = isvcd_set_processor(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case ISVCD_CMD_CTL_SET_TGT_LAYER:
            ret = isvcd_set_target_layer(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        default:
            H264_DEC_DEBUG_PRINT("\ndo nothing\n");
            break;
    }

    return ret;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name :   isvcd_rel_display_frame                                */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_rel_display_frame(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    ivd_rel_display_frame_ip_t *ps_rel_ip;
    ivd_rel_display_frame_op_t *ps_rel_op;
    UWORD32 buf_released = 0;

    UWORD32 u4_ts = 0;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    UWORD8 u1_layer_id;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_rel_ip = (ivd_rel_display_frame_ip_t *) pv_api_ip;
    ps_rel_op = (ivd_rel_display_frame_op_t *) pv_api_op;
    ps_rel_op->u4_error_code = 0;
    u4_ts = ps_rel_ip->u4_disp_buf_id;

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_dec = &ps_svc_lyr_dec->s_dec;
        if(0 == ps_dec->u4_share_disp_buf)
        {
            ps_dec->u4_disp_buf_mapping[u4_ts] = 0;
            ps_dec->u4_disp_buf_to_be_freed[u4_ts] = 0;
            return IV_SUCCESS;
        }

        if(ps_dec->pv_pic_buf_mgr != NULL)
        {
            if(1 == ps_dec->u4_disp_buf_mapping[u4_ts])
            {
                ih264_buf_mgr_release((buf_mgr_t *) ps_dec->pv_pic_buf_mgr,
                                      ps_rel_ip->u4_disp_buf_id, BUF_MGR_IO);
                ps_dec->u4_disp_buf_mapping[u4_ts] = 0;
                buf_released = 1;
            }
        }

        if((1 == ps_dec->u4_share_disp_buf) && (0 == buf_released))
            ps_dec->u4_disp_buf_to_be_freed[u4_ts] = 1;
    }
    return IV_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief
 *  Sets degrade params
 *
 * @par Description:
 *  Sets degrade params.
 *  Refer to ih264d_ctl_degrade_ip_t definition for details
 *
 * @param[in] ps_codec_obj
 *  Pointer to codec object at API level
 *
 * @param[in] pv_api_ip
 *  Pointer to input argument structure
 *
 * @param[out] pv_api_op
 *  Pointer to output argument structure
 *
 * @returns  Status
 *
 * @remarks
 *
 *
 *******************************************************************************
 */

WORD32 isvcd_set_degrade(iv_obj_t *ps_codec_obj, void *pv_api_ip, void *pv_api_op)
{
    isvcd_ctl_degrade_ip_t *ps_ip;
    isvcd_ctl_degrade_op_t *ps_op;
    dec_struct_t *ps_codec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    UWORD8 u1_layer_id;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) ps_codec_obj->pv_codec_handle;

    ps_ip = (isvcd_ctl_degrade_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_degrade_op_t *) pv_api_op;

    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_codec = &ps_svc_lyr_dec->s_dec;
        ps_codec->i4_degrade_type = ps_ip->i4_degrade_type;
        ps_codec->i4_nondegrade_interval = ps_ip->i4_nondegrade_interval;
        ps_codec->i4_degrade_pics = ps_ip->i4_degrade_pics;

        ps_codec->i4_degrade_pic_cnt = 0;
    }
    ps_op->u4_error_code = 0;

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_frame_dimensions                               */
/*                                                                           */
/*  Description   : gets the frame wd and ht and the buffer sizes            */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_frame_dimensions(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_ctl_get_frame_dimensions_ip_t *ps_ip;
    isvcd_ctl_get_frame_dimensions_op_t *ps_op;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    UWORD32 disp_wd, disp_ht, buffer_wd, buffer_ht, x_offset, y_offset;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;

    ps_ip = (isvcd_ctl_get_frame_dimensions_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_get_frame_dimensions_op_t *) pv_api_op;
    UNUSED(ps_ip);
    if((NULL != ps_dec->ps_cur_sps) && (1 == (ps_dec->ps_cur_sps->u1_is_valid)))
    {
        disp_wd = ps_dec->u2_disp_width;
        disp_ht = ps_dec->u2_disp_height;

        if(0 == ps_dec->u4_share_disp_buf)
        {
            buffer_wd = disp_wd;
            buffer_ht = disp_ht;
        }
        else
        {
            buffer_wd = ps_dec->u2_frm_wd_y;
            buffer_ht = ps_dec->u2_frm_ht_y;
        }
    }
    else
    {
        disp_wd = 0;
        disp_ht = 0;

        if(0 == ps_dec->u4_share_disp_buf)
        {
            buffer_wd = disp_wd;
            buffer_ht = disp_ht;
        }
        else
        {
            buffer_wd = ALIGN16(disp_wd) + (PAD_LEN_Y_H << 1);
            buffer_ht = ALIGN16(disp_ht) + (PAD_LEN_Y_V << 2);
        }
    }
    if(ps_dec->u4_app_disp_width > buffer_wd) buffer_wd = ps_dec->u4_app_disp_width;

    if(0 == ps_dec->u4_share_disp_buf)
    {
        x_offset = 0;
        y_offset = 0;
    }
    else
    {
        y_offset = (PAD_LEN_Y_V << 1);
        x_offset = PAD_LEN_Y_H;

        if((NULL != ps_dec->ps_sps) && (1 == (ps_dec->ps_sps->u1_is_valid)) &&
           (0 != ps_dec->u2_crop_offset_y))
        {
            y_offset += ps_dec->u2_crop_offset_y / ps_dec->u2_frm_wd_y;
            x_offset += ps_dec->u2_crop_offset_y % ps_dec->u2_frm_wd_y;
        }
    }

    ps_op->u4_disp_wd[0] = disp_wd;
    ps_op->u4_disp_ht[0] = disp_ht;
    ps_op->u4_buffer_wd[0] = buffer_wd;
    ps_op->u4_buffer_ht[0] = buffer_ht;
    ps_op->u4_x_offset[0] = x_offset;
    ps_op->u4_y_offset[0] = y_offset;

    ps_op->u4_disp_wd[1] = ps_op->u4_disp_wd[2] = ((ps_op->u4_disp_wd[0] + 1) >> 1);
    ps_op->u4_disp_ht[1] = ps_op->u4_disp_ht[2] = ((ps_op->u4_disp_ht[0] + 1) >> 1);
    ps_op->u4_buffer_wd[1] = ps_op->u4_buffer_wd[2] = (ps_op->u4_buffer_wd[0] >> 1);
    ps_op->u4_buffer_ht[1] = ps_op->u4_buffer_ht[2] = (ps_op->u4_buffer_ht[0] >> 1);
    ps_op->u4_x_offset[1] = ps_op->u4_x_offset[2] = (ps_op->u4_x_offset[0] >> 1);
    ps_op->u4_y_offset[1] = ps_op->u4_y_offset[2] = (ps_op->u4_y_offset[0] >> 1);

    if((ps_dec->u1_chroma_format == IV_YUV_420SP_UV) ||
       (ps_dec->u1_chroma_format == IV_YUV_420SP_VU))
    {
        ps_op->u4_disp_wd[2] = 0;
        ps_op->u4_disp_ht[2] = 0;
        ps_op->u4_buffer_wd[2] = 0;
        ps_op->u4_buffer_ht[2] = 0;
        ps_op->u4_x_offset[2] = 0;
        ps_op->u4_y_offset[2] = 0;

        ps_op->u4_disp_wd[1] <<= 1;
        ps_op->u4_buffer_wd[1] <<= 1;
        ps_op->u4_x_offset[1] <<= 1;
    }

    return IV_SUCCESS;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_vui_params                                     */
/*                                                                           */
/*  Description   : gets the VUI params                                      */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : success or failure                                       */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_get_vui_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_ctl_get_vui_params_ip_t *ps_ip;
    isvcd_ctl_get_vui_params_op_t *ps_op;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    dec_seq_params_t *ps_sps;
    vui_t *ps_vui;
    UWORD32 u4_size;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ip = (isvcd_ctl_get_vui_params_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_get_vui_params_op_t *) pv_api_op;
    UNUSED(ps_ip);

    u4_size = ps_op->u4_size;
    memset(ps_op, 0, sizeof(isvcd_ctl_get_vui_params_op_t));
    ps_op->u4_size = u4_size;

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    if(NULL == ps_dec->ps_cur_sps)
    {
        ps_op->u4_error_code = ERROR_VUI_PARAMS_NOT_FOUND;
        return IV_FAIL;
    }
    ps_sps = ps_dec->ps_cur_sps;

    if((0 == ps_sps->u1_is_valid) || (0 == ps_sps->u1_vui_parameters_present_flag))
    {
        ps_op->u4_error_code = ERROR_VUI_PARAMS_NOT_FOUND;
        return IV_FAIL;
    }

    ps_vui = &ps_sps->s_vui;

    ps_op->u1_aspect_ratio_idc = ps_vui->u1_aspect_ratio_idc;
    ps_op->u2_sar_width = ps_vui->u2_sar_width;
    ps_op->u2_sar_height = ps_vui->u2_sar_height;
    ps_op->u1_overscan_appropriate_flag = ps_vui->u1_overscan_appropriate_flag;
    ps_op->u1_video_format = ps_vui->u1_video_format;
    ps_op->u1_video_full_range_flag = ps_vui->u1_video_full_range_flag;
    ps_op->u1_colour_primaries = ps_vui->u1_colour_primaries;
    ps_op->u1_tfr_chars = ps_vui->u1_tfr_chars;
    ps_op->u1_matrix_coeffs = ps_vui->u1_matrix_coeffs;
    ps_op->u1_cr_top_field = ps_vui->u1_cr_top_field;
    ps_op->u1_cr_bottom_field = ps_vui->u1_cr_bottom_field;
    ps_op->u4_num_units_in_tick = ps_vui->u4_num_units_in_tick;
    ps_op->u4_time_scale = ps_vui->u4_time_scale;
    ps_op->u1_fixed_frame_rate_flag = ps_vui->u1_fixed_frame_rate_flag;
    ps_op->u1_nal_hrd_params_present = ps_vui->u1_nal_hrd_params_present;
    ps_op->u1_vcl_hrd_params_present = ps_vui->u1_vcl_hrd_params_present;
    ps_op->u1_low_delay_hrd_flag = ps_vui->u1_low_delay_hrd_flag;
    ps_op->u1_pic_struct_present_flag = ps_vui->u1_pic_struct_present_flag;
    ps_op->u1_bitstream_restriction_flag = ps_vui->u1_bitstream_restriction_flag;
    ps_op->u1_mv_over_pic_boundaries_flag = ps_vui->u1_mv_over_pic_boundaries_flag;
    ps_op->u4_max_bytes_per_pic_denom = ps_vui->u4_max_bytes_per_pic_denom;
    ps_op->u4_max_bits_per_mb_denom = ps_vui->u4_max_bits_per_mb_denom;
    ps_op->u4_log2_max_mv_length_horz = ps_vui->u4_log2_max_mv_length_horz;
    ps_op->u4_log2_max_mv_length_vert = ps_vui->u4_log2_max_mv_length_vert;
    ps_op->u4_num_reorder_frames = ps_vui->u4_num_reorder_frames;
    ps_op->u4_max_dec_frame_buffering = ps_vui->u4_max_dec_frame_buffering;

    return IV_SUCCESS;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_sei_mdcv_params                               */
/*                                                                           */
/*  Description   : This function populates SEI mdcv message in              */
/*                     output structure                                      */
/*  Inputs        : iv_obj_t decoder handle                                  */
/*                : pv_api_ip pointer to input structure                     */
/*                : pv_api_op pointer to output structure                    */
/*  Outputs       :                                                          */
/*  Returns       : returns 0; 1 with error code when MDCV is not present    */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_sei_mdcv_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_ctl_get_sei_mdcv_params_ip_t *ps_ip;
    isvcd_ctl_get_sei_mdcv_params_op_t *ps_op;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    sei_mdcv_params_t *ps_sei_mdcv;
    WORD32 i4_count;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ip = (isvcd_ctl_get_sei_mdcv_params_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_get_sei_mdcv_params_op_t *) pv_api_op;
    UNUSED(ps_ip);

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    if(0 == ps_dec->s_sei_export.u1_sei_mdcv_params_present_flag)
    {
        ps_op->u4_error_code = ERROR_SEI_MDCV_PARAMS_NOT_FOUND;
        return IV_FAIL;
    }
    ps_sei_mdcv = &ps_dec->s_sei_export.s_sei_mdcv_params;

    for(i4_count = 0; i4_count < NUM_SEI_MDCV_PRIMARIES; i4_count++)
    {
        ps_op->au2_display_primaries_x[i4_count] = ps_sei_mdcv->au2_display_primaries_x[i4_count];
        ps_op->au2_display_primaries_y[i4_count] = ps_sei_mdcv->au2_display_primaries_y[i4_count];
    }

    ps_op->u2_white_point_x = ps_sei_mdcv->u2_white_point_x;
    ps_op->u2_white_point_y = ps_sei_mdcv->u2_white_point_y;
    ps_op->u4_max_display_mastering_luminance = ps_sei_mdcv->u4_max_display_mastering_luminance;
    ps_op->u4_min_display_mastering_luminance = ps_sei_mdcv->u4_min_display_mastering_luminance;

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_sei_cll_params                                */
/*                                                                           */
/*  Description   : This function populates SEI cll message in               */
/*                     output structure                                      */
/*  Inputs        : iv_obj_t decoder handle                                  */
/*                : pv_api_ip pointer to input structure                     */
/*                : pv_api_op pointer to output structure                    */
/*  Outputs       :                                                          */
/*  Returns       : returns 0; 1 with error code when CLL is not present     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_sei_cll_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_ctl_get_sei_cll_params_ip_t *ps_ip;
    isvcd_ctl_get_sei_cll_params_op_t *ps_op;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    sei_cll_params_t *ps_sei_cll;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ip = (isvcd_ctl_get_sei_cll_params_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_get_sei_cll_params_op_t *) pv_api_op;
    UNUSED(ps_ip);

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    if(0 == ps_dec->s_sei_export.u1_sei_cll_params_present_flag)
    {
        ps_op->u4_error_code = ERROR_SEI_CLL_PARAMS_NOT_FOUND;
        return IV_FAIL;
    }
    ps_sei_cll = &ps_dec->s_sei_export.s_sei_cll_params;

    ps_op->u2_max_content_light_level = ps_sei_cll->u2_max_content_light_level;
    ps_op->u2_max_pic_average_light_level = ps_sei_cll->u2_max_pic_average_light_level;

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_sei_ave_params                                */
/*                                                                           */
/*  Description   : This function populates SEI ave message in               */
/*                     output structure                                      */
/*  Inputs        : iv_obj_t decoder handle                                  */
/*                : pv_api_ip pointer to input structure                     */
/*                : pv_api_op pointer to output structure                    */
/*  Outputs       :                                                          */
/*  Returns       : returns 0; 1 with error code when AVE is not present     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_sei_ave_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_ctl_get_sei_ave_params_ip_t *ps_ip;
    isvcd_ctl_get_sei_ave_params_op_t *ps_op;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    sei_ave_params_t *ps_sei_ave;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ip = (isvcd_ctl_get_sei_ave_params_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_get_sei_ave_params_op_t *) pv_api_op;
    UNUSED(ps_ip);

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    if(0 == ps_dec->s_sei_export.u1_sei_ave_params_present_flag)
    {
        ps_op->u4_error_code = ERROR_SEI_AVE_PARAMS_NOT_FOUND;
        return IV_FAIL;
    }
    ps_sei_ave = &ps_dec->s_sei_export.s_sei_ave_params;

    ps_op->u4_ambient_illuminance = ps_sei_ave->u4_ambient_illuminance;
    ps_op->u2_ambient_light_x = ps_sei_ave->u2_ambient_light_x;
    ps_op->u2_ambient_light_y = ps_sei_ave->u2_ambient_light_y;

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_sei_ccv_params                                */
/*                                                                           */
/*  Description   : This function populates SEI mdcv message in              */
/*                     output structure                                      */
/*  Inputs        : iv_obj_t decoder handle                                  */
/*                : pv_api_ip pointer to input structure                     */
/*                : pv_api_op pointer to output structure                    */
/*  Outputs       :                                                          */
/*  Returns       : returns 0; 1 with error code when CCV is not present    */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_sei_ccv_params(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    isvcd_ctl_get_sei_ccv_params_ip_t *ps_ip;
    isvcd_ctl_get_sei_ccv_params_op_t *ps_op;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    sei_ccv_params_t *ps_sei_ccv;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    WORD32 i4_count;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ip = (isvcd_ctl_get_sei_ccv_params_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_get_sei_ccv_params_op_t *) pv_api_op;
    UNUSED(ps_ip);

    ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[ps_svcd_ctxt->u1_target_layer_id];
    ps_dec = &ps_svc_lyr_dec->s_dec;
    if(0 == ps_dec->s_sei_export.u1_sei_ccv_params_present_flag)
    {
        ps_op->u4_error_code = ERROR_SEI_CCV_PARAMS_NOT_FOUND;
        return IV_FAIL;
    }
    ps_sei_ccv = &ps_dec->s_sei_export.s_sei_ccv_params;
    ps_op->u1_ccv_cancel_flag = ps_sei_ccv->u1_ccv_cancel_flag;

    if(0 == ps_op->u1_ccv_cancel_flag)
    {
        ps_op->u1_ccv_persistence_flag = ps_sei_ccv->u1_ccv_persistence_flag;
        ps_op->u1_ccv_primaries_present_flag = ps_sei_ccv->u1_ccv_primaries_present_flag;
        ps_op->u1_ccv_min_luminance_value_present_flag =
            ps_sei_ccv->u1_ccv_min_luminance_value_present_flag;
        ps_op->u1_ccv_max_luminance_value_present_flag =
            ps_sei_ccv->u1_ccv_max_luminance_value_present_flag;
        ps_op->u1_ccv_avg_luminance_value_present_flag =
            ps_sei_ccv->u1_ccv_avg_luminance_value_present_flag;
        ps_op->u1_ccv_reserved_zero_2bits = ps_sei_ccv->u1_ccv_reserved_zero_2bits;

        if(1 == ps_sei_ccv->u1_ccv_primaries_present_flag)
        {
            for(i4_count = 0; i4_count < NUM_SEI_CCV_PRIMARIES; i4_count++)
            {
                ps_op->ai4_ccv_primaries_x[i4_count] = ps_sei_ccv->ai4_ccv_primaries_x[i4_count];
                ps_op->ai4_ccv_primaries_y[i4_count] = ps_sei_ccv->ai4_ccv_primaries_y[i4_count];
            }
        }

        if(1 == ps_sei_ccv->u1_ccv_min_luminance_value_present_flag)
        {
            ps_op->u4_ccv_min_luminance_value = ps_sei_ccv->u4_ccv_min_luminance_value;
        }
        if(1 == ps_sei_ccv->u1_ccv_max_luminance_value_present_flag)
        {
            ps_op->u4_ccv_max_luminance_value = ps_sei_ccv->u4_ccv_max_luminance_value;
        }
        if(1 == ps_sei_ccv->u1_ccv_avg_luminance_value_present_flag)
        {
            ps_op->u4_ccv_avg_luminance_value = ps_sei_ccv->u4_ccv_avg_luminance_value;
        }
    }

    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_set_num_cores                                      */
/*                                                                           */
/*  Description   : This function sets the no of cores which decoder         */
/*                     can use for decoding                                  */
/*  Inputs        : iv_obj_t decoder handle                                  */
/*                : pv_api_ip pointer to input structure                     */
/*                : pv_api_op pointer to output structure                    */
/*  Outputs       :                                                          */
/*  Returns       : returns 0; 1 with error code                             */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_set_num_cores(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    UWORD8 u1_layer_id;
    isvcd_ctl_set_num_cores_ip_t *ps_ip;
    isvcd_ctl_set_num_cores_op_t *ps_op;
    dec_struct_t *ps_dec;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec;
    svc_dec_ctxt_t *ps_svcd_ctxt;
    ps_svcd_ctxt = (svc_dec_ctxt_t *) dec_hdl->pv_codec_handle;

    ps_ip = (isvcd_ctl_set_num_cores_ip_t *) pv_api_ip;
    ps_op = (isvcd_ctl_set_num_cores_op_t *) pv_api_op;
    ps_op->u4_error_code = 0;
    for(u1_layer_id = 0; u1_layer_id < MAX_NUM_RES_LYRS; u1_layer_id++)
    {
        ps_svc_lyr_dec = &ps_svcd_ctxt->ps_svc_dec_lyr[u1_layer_id];
        ps_dec = &ps_svc_lyr_dec->s_dec;
        ps_dec->u4_num_cores = ps_ip->u4_num_cores;

        if(ps_dec->u4_num_cores == 1)
        {
            ps_dec->u1_separate_parse = 0;
        }
        else
        {
            ps_dec->u1_separate_parse = 1;
        }

        /*using only upto three threads currently*/
        if(ps_dec->u4_num_cores > 3) ps_dec->u4_num_cores = 3;
    }
    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_fill_output_struct_from_context                    */
/*                                                                           */
/*  Description   : This function fills the output structure from the        */
/*                     svc layer context                                     */
/*  Inputs        : iv_obj_t decoder handle                                  */
/*                : ps_svc_lyr_dec pointer to svc layer context              */
/*                : ps_dec_op pointer to output structure                    */
/*  Outputs       :                                                          */
/*  Returns       :                                                          */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
void isvcd_fill_output_struct_from_context(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                           ivd_video_decode_op_t *ps_dec_op)
{
    dec_struct_t *ps_dec;
    ps_dec = &ps_svc_lyr_dec->s_dec;
    if((ps_dec_op->u4_error_code & 0xff) != ERROR_DYNAMIC_RESOLUTION_NOT_SUPPORTED)
    {
        ps_dec_op->u4_pic_wd = (UWORD32) ps_dec->u2_disp_width;
        ps_dec_op->u4_pic_ht = (UWORD32) ps_dec->u2_disp_height;
    }
    ps_dec_op->i4_reorder_depth = ps_dec->i4_reorder_depth;
    ps_dec_op->i4_display_index = ps_dec->i4_display_index;
    ps_dec_op->e_pic_type = ps_dec->i4_frametype;

    ps_dec_op->u4_new_seq = 0;
    ps_dec_op->u4_output_present =
        (ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER) ? ps_dec->u4_output_present : 0;
    ps_dec_op->u4_progressive_frame_flag = ps_dec->s_disp_op.u4_progressive_frame_flag;

    ps_dec_op->u4_is_ref_flag = 1;
    if(ps_dec_op->u4_frame_decoded_flag)
    {
        if(ps_dec->ps_cur_slice->u1_nal_ref_idc == 0) ps_dec_op->u4_is_ref_flag = 0;
    }

    ps_dec_op->e_output_format = ps_dec->s_disp_op.e_output_format;
    ps_dec_op->s_disp_frm_buf = ps_dec->s_disp_op.s_disp_frm_buf;
    ps_dec_op->e4_fld_type = ps_dec->s_disp_op.e4_fld_type;
    ps_dec_op->u4_ts = ps_dec->s_disp_op.u4_ts;
    ps_dec_op->u4_disp_buf_id = ps_dec->s_disp_op.u4_disp_buf_id;

    ih264d_export_sei_params(&ps_dec_op->s_sei_decode_op, ps_dec);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_api_function                                       */
/*                                                                           */
/*  Description   :                                                          */
/*                                                                           */
/*  Inputs        :iv_obj_t decoder handle                                   */
/*                :pv_api_ip pointer to input structure                      */
/*                :pv_api_op pointer to output structure                     */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
IV_API_CALL_STATUS_T isvcd_api_function(iv_obj_t *dec_hdl, void *pv_api_ip, void *pv_api_op)
{
    UWORD32 command;
    UWORD32 *pu2_ptr_cmd;
    UWORD32 u4_api_ret;
    IV_API_CALL_STATUS_T e_status;
    e_status = api_check_struct_sanity(dec_hdl, pv_api_ip, pv_api_op);

    if(e_status != IV_SUCCESS)
    {
        UWORD32 *ptr_err;

        ptr_err = (UWORD32 *) pv_api_op;
        UNUSED(ptr_err);
        H264_DEC_DEBUG_PRINT("error code = %d\n", *(ptr_err + 1));
        return IV_FAIL;
    }

    pu2_ptr_cmd = (UWORD32 *) pv_api_ip;
    pu2_ptr_cmd++;

    command = *pu2_ptr_cmd;
    switch(command)
    {
        case IVD_CMD_CREATE:
            u4_api_ret = isvcd_create(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        case IVD_CMD_DELETE:
            u4_api_ret = isvcd_delete(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;

        case IVD_CMD_VIDEO_DECODE:
            u4_api_ret = isvcd_video_decode(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;

        case IVD_CMD_GET_DISPLAY_FRAME:
            u4_api_ret = ih264d_get_display_frame(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);

            break;

        case IVD_CMD_SET_DISPLAY_FRAME:
            u4_api_ret = isvcd_set_display_frame(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);

            break;

        case IVD_CMD_REL_DISPLAY_FRAME:
            u4_api_ret = isvcd_rel_display_frame(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;

        case IVD_CMD_VIDEO_CTL:
            u4_api_ret = isvcd_ctl(dec_hdl, (void *) pv_api_ip, (void *) pv_api_op);
            break;
        default:
            u4_api_ret = IV_FAIL;
            break;
    }

    return u4_api_ret;
}
