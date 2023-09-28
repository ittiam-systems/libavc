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
/*****************************************************************************/
/*                                                                           */
/*  File Name         : app.h                                                */
/*                                                                           */
/*  Description       : This file contains all the necessary structure and   */
/*                      enumeration definitions needed for the Application   */
/*                                                                           */
/*  List of Functions :                                                      */
/*                                                                           */
/*  Issues / Problems : None                                                 */
/*                                                                           */
/*  Revision History  :                                                      */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         26 08 2010   Ittiam          Draft                                */
/*                                                                           */
/*****************************************************************************/

#ifndef _SVCE_APP_H_
#define _SVCE_APP_H_

#include <stdbool.h>
#include <sys/time.h>

#include "iv2.h"
#include "ive2.h"

/*****************************************************************************/
/* Function Macros                                                           */
/*****************************************************************************/
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define MIN(a, b) ((a) < (b)) ? (a) : (b)

#define ALIGN16(x) ((((x) + 15) >> 4) << 4)

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

#define DEFAULT_NUM_INPUT_BUFS 32
#define DEFAULT_MAX_INPUT_BUFS 32

#define DEFAULT_NUM_OUTPUT_BUFS 32
#define DEFAULT_MAX_OUTPUT_BUFS 32

#define DEFAULT_NUM_RECON_BUFS 32
#define DEFAULT_MAX_RECON_BUFS DEFAULT_NUM_RECON_BUFS

#define DEFAULT_NUM_NALU_INFO_BUFS 32
#define DEFAULT_MAX_NALU_INFO_BUFS 32

#define LEN_STATUS_BUFFER (10 * 1024)
#define MAX_VBV_BUFF_SIZE (120 * 16384)
#define MAX_NUM_IO_BUFS 3

#define DEFAULT_MAX_REF_FRM 2
#define DEFAULT_MAX_REORDER_FRM 0
#define DEFAULT_QP_MIN 4
#define DEFAULT_QP_MAX 51
#define DEFAULT_MAX_BITRATE 240000000
#define DEFAULT_NUM_BFRAMES 0
#define DEFAULT_MAX_SRCH_RANGE_X 256
#define DEFAULT_MAX_SRCH_RANGE_Y 256
#define DEFAULT_MAX_FRAMERATE 120000
#define DEFAULT_NUM_CORES 1
#define DEFAULT_NUM_CORES_PRE_ENC 0
#define DEFAULT_FPS 30
#define DEFAULT_ENC_SPEED 100
#define DEFAULT_MEM_REC_CNT 0
#define DEFAULT_RECON_ENABLE 0
#define DEFAULT_NALU_INFO_EXPORT_ENABLE 0
#define DEFAULT_CHKSUM_ENABLE 0
#define DEFAULT_START_FRM 0
#define DEFAULT_NUM_FRMS 100
#define DEFAULT_INP_COLOR_FMT IV_YUV_420P
#define DEFAULT_RECON_COLOR_FMT IV_YUV_420P
#define DEFAULT_LOOPBACK 0
#define DEFAULT_SRC_FRAME_RATE 50
#define DEFAULT_TGT_FRAME_RATE 50
#define DEFAULT_MAX_WD 1920
#define DEFAULT_MAX_HT 1920
#define DEFAULT_MAX_LEVEL 51
#define DEFAULT_STRIDE 0
#define DEFAULT_WD 1920
#define DEFAULT_HT 1080
#define DEFAULT_PSNR_ENABLE 0
#define DEFAULT_ME_SPEED 100
#define DEFAULT_ENABLE_FAST_SAD 0
#define DEFAULT_ENABLE_ALT_REF 0
#define DEFAULT_RC 0
#define DEFAULT_BITRATE 6000000
#define DEFAULT_I_QP 25
#define DEFAULT_I_QP_MAX DEFAULT_QP_MAX
#define DEFAULT_I_QP_MIN 0
#define DEFAULT_P_QP 28
#define DEFAULT_P_QP_MAX DEFAULT_QP_MAX
#define DEFAULT_P_QP_MIN 0
#define DEFAULT_B_QP 28
#define DEFAULT_B_QP_MAX DEFAULT_QP_MAX
#define DEFAULT_B_QP_MIN 0
#define DEFAULT_AIR 0
#define DEFAULT_AIR_REFRESH_PERIOD 30
#define DEFAULT_SRCH_RNG_X 64
#define DEFAULT_SRCH_RNG_Y 48
#define DEFAULT_I_INTERVAL 50
#define DEFAULT_IDR_INTERVAL 100
#define DEFAULT_B_FRAMES 0
#define DEFAULT_DISABLE_DEBLK_LEVEL 4
#define DEFAULT_HPEL 1
#define DEFAULT_QPEL 1
#define DEFAULT_I4 1
#define DEFAULT_EPROFILE IV_PROFILE_BASE
#define DEFAULT_SLICE_MODE 0
#define DEFAULT_SLICE_PARAM 256
#define DEFAULT_ENTROPY_CODING_MODE 1
#define DEFAULT_NUM_TEMPORAL_LAYERS 1
#define DEFAULT_NUM_SPATIAL_LAYERS 1
#define DEFAULT_SPATIAL_RES_RATIO 2.0

#define DEFAULT_MAX_DISPLAY_MASTERING_LUMINANCE 50000
#define DEFAULT_MIN_DISPLAY_MASTERING_LUMINANCE 1

#define STRLENGTH 500

/* specifies the number of colour primary components of the mastering
   display */
#define NUM_SEI_MDCV_PRIMARIES 3

/* specifies the number of colour primary components of the nominal
   content colour volume */
#define NUM_SEI_CCV_PRIMARIES 3

/*****************************************************************************/
/*  profile Macros                                                           */
/*****************************************************************************/
#ifdef PROFILE_ENABLE
#ifdef WINDOWS
typedef LARGE_INTEGER TIMER;
#else
typedef struct timeval TIMER;
#endif
#else
typedef int32_t TIMER;
#endif

#ifdef PROFILE_ENABLE
#ifdef WINDOWS
#define GETTIME(timer) QueryPerformanceCounter(timer);
#else
#define GETTIME(timer) gettimeofday(timer, NULL);
#endif

#ifdef WINDOWS
#define ELAPSEDTIME(s_start_timer, s_end_timer, s_elapsed_time, frequency)                     \
    {                                                                                          \
        TIMER s_temp_time;                                                                     \
        s_temp_time.LowPart = s_end_timer.LowPart - s_start_timer.LowPart;                     \
        s_elapsed_time =                                                                       \
            (UWORD32) (((DOUBLE) s_temp_time.LowPart / (DOUBLE) frequency.LowPart) * 1000000); \
    }
#else
#define ELAPSEDTIME(s_start_timer, s_end_timer, s_elapsed_time, frequency)     \
    s_elapsed_time = ((s_end_timer.tv_sec - s_start_timer.tv_sec) * 1000000) + \
                     (s_end_timer.tv_usec - s_start_timer.tv_usec);
#endif
#else
#define GETTIME(timer)
#define ELAPSEDTIME(s_start_timer, s_end_timer, s_elapsed_time, frequency)
#endif

/*****************************************************************************/
/*  Structure definitions                                                    */
/*****************************************************************************/
typedef struct
{
    UWORD8 *pu1_buf;
    UWORD32 u4_buf_size;
    UWORD32 u4_timestamp_low;
    UWORD32 u4_timestamp_high;
    UWORD32 u4_is_free;
    void *pv_mb_info;
    void *pv_pic_info;
} input_buf_t;

typedef struct
{
    UWORD8 *pu1_buf;
    UWORD32 u4_buf_size;
    UWORD32 u4_timestamp_low;
    UWORD32 u4_timestamp_high;
    UWORD32 u4_is_free;
} output_buf_t;

typedef struct
{
    UWORD8 *pu1_buf;
    UWORD32 u4_buf_size;
    UWORD32 u4_timestamp_low;
    UWORD32 u4_timestamp_high;
    UWORD32 u4_is_free;
} recon_buf_t;

typedef struct nalu_info_buf_t
{
    UWORD8 *pu1_buf;

    UWORD32 u4_buf_size;

    bool b_is_free;
} nalu_info_buf_t;

typedef struct
{
    iv_obj_t *ps_enc;
    iv_mem_rec_t *ps_mem_rec;
    UWORD32 u4_num_mem_rec;
    UWORD32 u4_recon_enable;
    UWORD32 u4_chksum_enable;
    UWORD32 u4_nalu_info_export_enable;
    UWORD32 u4_mb_info_type;
    UWORD32 u4_pic_info_type;
    UWORD32 u4_mb_info_size;
    UWORD32 u4_pic_info_size;
    UWORD32 u4_start_frm;
    UWORD32 u4_max_num_frms;
    UWORD32 u4_total_bytes;
    UWORD32 u4_pics_cnt;
    IV_COLOR_FORMAT_T e_inp_color_fmt;
    IV_COLOR_FORMAT_T e_recon_color_fmt;
    IV_ARCH_T e_arch;
    IV_SOC_T e_soc;

    WORD32 header_generated;
    void *pv_codec_obj;

    UWORD32 u4_num_cores;
    UWORD32 u4_pre_enc_me;
    UWORD32 u4_pre_enc_ipe;
    CHAR ac_ip_fname[STRLENGTH];
    CHAR ac_op_fname[STRLENGTH];
    CHAR ac_recon_fname[STRLENGTH];
    CHAR ac_nalu_info_csv_fname[STRLENGTH];
    CHAR ac_chksum_fname[STRLENGTH];
    CHAR ac_mb_info_fname[STRLENGTH];
    CHAR ac_pic_info_fname[STRLENGTH];

    FILE *fp_ip;
    FILE *fp_op;
    FILE *fp_recon;
    FILE *fp_nalu_info;
    FILE *fp_chksum;
    FILE *fp_psnr_ip;
    FILE *fp_mb_info;
    FILE *fp_pic_info;
    FILE *fp_dump_op;

    UWORD32 u4_loopback;
    UWORD32 u4_max_frame_rate;
    UWORD32 u4_src_frame_rate;
    UWORD32 u4_tgt_frame_rate;
    UWORD32 u4_max_wd;
    UWORD32 u4_max_ht;
    UWORD32 u4_max_level;

    UWORD32 u4_strd;

    UWORD32 u4_wd;
    UWORD32 u4_ht;

    UWORD32 u4_enc_wd;
    UWORD32 u4_enc_ht;

    UWORD32 u4_psnr_enable;

    UWORD32 u4_enc_speed;
    UWORD32 u4_me_speed;
    UWORD32 u4_enable_fast_sad;
    UWORD32 u4_enable_alt_ref;
    UWORD32 u4_rc;
    UWORD32 *pu4_max_bitrate;
    UWORD32 *pu4_bitrate;
    UWORD32 *pu4_i_qp, *pu4_i_qp_max, *pu4_i_qp_min;
    UWORD32 *pu4_p_qp, *pu4_p_qp_max, *pu4_p_qp_min;
    UWORD32 *pu4_b_qp, *pu4_b_qp_max, *pu4_b_qp_min;
    UWORD32 u4_air;
    UWORD32 u4_air_refresh_period;
    UWORD32 u4_srch_rng_x;
    UWORD32 u4_srch_rng_y;
    UWORD32 u4_i_interval;
    UWORD32 u4_idr_interval;
    UWORD32 u4_b_frames;
    UWORD32 u4_num_bframes;
    UWORD32 u4_disable_deblock_level;
    UWORD32 u4_hpel;
    UWORD32 u4_qpel;
    UWORD32 u4_enable_intra_4x4;
    IV_PROFILE_T e_profile;

    UWORD32 u4_slice_mode;
    UWORD32 u4_slice_param;
    UWORD32 u4_entropy_coding_mode;

    void *pv_input_thread_handle;
    void *pv_output_thread_handle;
    void *pv_recon_thread_handle;

    isvce_ctl_getbufinfo_op_t s_get_buf_info_op;
    input_buf_t as_input_buf[DEFAULT_MAX_INPUT_BUFS];
    output_buf_t as_output_buf[DEFAULT_MAX_OUTPUT_BUFS];
    recon_buf_t as_recon_buf[DEFAULT_MAX_RECON_BUFS];
    nalu_info_buf_t as_nalu_info_bufs[DEFAULT_MAX_NALU_INFO_BUFS];

    DOUBLE adbl_psnr[3];
    UWORD32 u4_psnr_cnt;
    UWORD8 *pu1_psnr_buf;
    UWORD8 u4_psnr_buf_size;

    UWORD32 *pu4_vbv_buffer_delay;

    TIMER enc_start_time;
    TIMER enc_last_time;
    WORD32 avg_time;

    UWORD32 u4_sei_mdcv_params_present_flag;
    UWORD32 au4_display_primaries_x[NUM_SEI_MDCV_PRIMARIES];
    UWORD32 au4_display_primaries_y[NUM_SEI_MDCV_PRIMARIES];
    UWORD32 u4_white_point_x;
    UWORD32 u4_white_point_y;
    UWORD32 u4_max_display_mastering_luminance;
    UWORD32 u4_min_display_mastering_luminance;

    UWORD32 u4_sei_cll_params_present_flag;
    UWORD32 u4_max_content_light_level;
    UWORD32 u4_max_pic_average_light_level;

    UWORD32 u4_sei_ave_params_present_flag;
    UWORD32 u4_ambient_illuminance;
    UWORD32 u4_ambient_light_x;
    UWORD32 u4_ambient_light_y;

    UWORD32 u4_sei_ccv_params_present_flag;
    UWORD32 u4_ccv_cancel_flag;
    UWORD32 u4_ccv_persistence_flag;
    UWORD32 u4_ccv_primaries_present_flag;
    UWORD32 u4_ccv_min_luminance_value_present_flag;
    UWORD32 u4_ccv_max_luminance_value_present_flag;
    UWORD32 u4_ccv_avg_luminance_value_present_flag;
    UWORD32 u4_ccv_reserved_zero_2bits;
    WORD32 ai4_ccv_primaries_x[NUM_SEI_CCV_PRIMARIES];
    WORD32 ai4_ccv_primaries_y[NUM_SEI_CCV_PRIMARIES];
    UWORD32 u4_ccv_min_luminance_value;
    UWORD32 u4_ccv_max_luminance_value;
    UWORD32 u4_ccv_avg_luminance_value;
    UWORD32 u4_use_default_vui;

    isvce_ctl_set_sei_mdcv_params_ip_t s_sei_mdcv_params;
    isvce_ctl_set_sei_cll_params_ip_t s_sei_cll_params;
    isvce_ctl_set_sei_ave_params_ip_t s_sei_ave_params;

    UWORD8 u1_num_temporal_layers;
    UWORD8 u1_num_spatial_layers;
    DOUBLE d_spatial_res_ratio;

} app_ctxt_t;

/*****************************************************************************/
/*  Function Declarations                                                    */
/*****************************************************************************/
void codec_exit(CHAR *pc_err_message);
void allocate_input(app_ctxt_t *ps_app_ctxt);
void allocate_output(app_ctxt_t *ps_app_ctxt);
void allocate_recon(app_ctxt_t *ps_app_ctxt);

IV_STATUS_T read_input(FILE *fp, iv_raw_buf_t *ps_raw_buf);
IV_STATUS_T write_recon(FILE *fp, iv_raw_buf_t *ps_raw_buf);
IV_STATUS_T write_output(FILE *fp, UWORD8 *pu1_buf, WORD32 num_bytes);

IV_STATUS_T read_mb_info(app_ctxt_t *ps_app_ctxt, void *pv_mb_info);
IV_STATUS_T read_pic_info(app_ctxt_t *ps_app_ctxt, void *pv_pic_info);

void *isvca_aligned_malloc(WORD32 alignment, WORD32 size);
void isvca_aligned_free(void *pv_buf);

void free_input(app_ctxt_t *ps_app_ctxt);
void free_recon(app_ctxt_t *ps_app_ctxt);
void free_output(app_ctxt_t *ps_app_ctxt);

void init_raw_buf_descr(app_ctxt_t *ps_app_ctxt, iv_raw_buf_t *ps_raw_buf, UWORD8 *pu1_buf,
                        IV_COLOR_FORMAT_T e_color_fmt);

#ifndef MD5_DISABLE
void calc_md5_cksum(UWORD8 *pu1_inbuf, UWORD32 u4_stride, UWORD32 u4_width, UWORD32 u4_height,
                    UWORD8 *pu1_cksum_p);
#else
#define calc_md5_cksum(a, b, c, d, e)
#endif

#endif
