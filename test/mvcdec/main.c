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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>

#include "ih264_typedefs.h"
#include "iv.h"
#include "ivd.h"
#include "imvcd.h"
#include "ih264_debug.h"
#include "ih264d.h"
#include "ithread.h"

/* Constants */
#define DEFAULT_NON_DEGRADE_INTERVAL 4

#define MAX_DISP_BUFFERS 64

#define STRLENGTH 1000

#define PEAK_WINDOW_SIZE 8

#define DEFAULT_COLOR_FORMAT IV_YUV_420P

#define DEFAULT_NUM_CORES 1

#define MAX_ARG_SHORTNAME_LENGTH 4

#define MAX_ARG_NAME_LENGTH 128

#define MAX_ARG_DESC_LENGTH 512

#define MAX_NUM_VIEWS 6

#undef PROFILE_ENABLE

/* Macro functions */
#ifdef PROFILE_ENABLE
#define GETTIME(timer) gettimeofday(timer, NULL);

#define ELAPSEDTIME(s_start_timer, s_end_timer, s_elapsed_time, frequency)     \
    s_elapsed_time = ((s_end_timer.tv_sec - s_start_timer.tv_sec) * 1000000) + \
                     (s_end_timer.tv_usec - s_start_timer.tv_usec);
#else
#define GETTIME(timer)
#define ELAPSEDTIME(s_start_timer, s_end_timer, s_elapsed_time, frequency)
#endif

/* Typedefs */
typedef enum ARGUMENT_T
{
    INVALID,
    HELP,
    INPUT_FILE,
    OUTPUT,
    CHKSUM,
    SAVE_OUTPUT,
    SAVE_CHKSUM,
    NUM_FRAMES,
    NUM_CORES,
    DISABLE_DEBLOCK_LEVEL,
    LOOPBACK,
    CONFIG,
    DEGRADE_TYPE,
    DEGRADE_PICS,
    ARCH
} ARGUMENT_T;

typedef enum COMPONENT_TYPES_T
{
    Y = 0,
    UV = 1,
    U = 1,
    V = 2,
    NUM_SP_COMPONENTS = 2,
    NUM_COMPONENTS = 3
} COMPONENT_TYPES_T;

#ifdef PROFILE_ENABLE
typedef struct timeval TIMER;
#else
typedef WORD32 TIMER;
#endif

typedef struct mvc_app_files_t
{
    UWORD8 au1_ip_fname[STRLENGTH];

    UWORD8 au1_op_fname[STRLENGTH];

    UWORD8 au1_op_chksum_fname[STRLENGTH];

    UWORD32 au4_disp_frm_id_queue[MAX_DISP_BUFFERS];
} mvc_app_files_t;

typedef struct mvc_dec_ctx_t
{
    mvc_app_files_t s_mvc_app_files;

    imvcd_get_buf_info_op_t s_disp_buf_props;

    iv_yuv_buf_t as_view_disp_bufs[MAX_NUM_VIEWS];

    iv_obj_t *ps_codec_obj;

    IV_COLOR_FORMAT_T e_output_chroma_format;

    IVD_ARCH_T e_arch;

    IVD_SOC_T e_soc;

    UWORD32 u4_max_frm_ts;

    UWORD32 u4_disable_dblk_level;

    WORD32 i4_degrade_type;

    WORD32 i4_degrade_pics;

    UWORD32 u4_num_cores;

    UWORD32 u4_disp_delay;

    UWORD32 u4_fps;

    UWORD32 u4_pic_wd;

    UWORD32 u4_pic_ht;

    UWORD32 u4_loopback;

    UWORD32 u4_file_save_flag;

    UWORD32 u4_chksum_save_flag;

    UWORD8 u1_quit;

} mvc_dec_ctx_t;

typedef struct argument_t
{
    UWORD8 au1_argument_shortname[MAX_ARG_SHORTNAME_LENGTH];

    UWORD8 au1_argument_name[MAX_ARG_NAME_LENGTH];

    ARGUMENT_T e_argument;

    UWORD8 au1_description[MAX_ARG_DESC_LENGTH];

} argument_t;

/* Function declarations */
#ifndef MD5_DISABLE
void calc_md5_cksum(UWORD8 *pu1_inbuf, UWORD32 u4_stride, UWORD32 u4_width, UWORD32 u4_height,
                    UWORD8 *pu1_cksum_p);
#else
#define calc_md5_cksum(a, b, c, d, e)
#endif

static inline void *mvcd_aligned_malloc(void *pv_ctxt, WORD32 alignment, WORD32 i4_size)
{
    void *buf = NULL;
    (void) pv_ctxt;
    if(0 != posix_memalign(&buf, alignment, i4_size))
    {
        return NULL;
    }
    return buf;
}

static inline void mvcd_aligned_free(void *pv_ctxt, void *pv_buf)
{
    (void) pv_ctxt;
    free(pv_buf);
    return;
}

static const argument_t argument_mapping[] = {
    {"-h", "--help", HELP, "Print this help\n"},
    {"-c", "--config", CONFIG, "config file (Default: test.cfg)\n"},
    {"-i", "--input", INPUT_FILE, "Input file\n"},
    {"-o", "--output", OUTPUT, "Output file\n"},
    {"--", "--chksum", CHKSUM, "Output MD5 Checksum file\n"},
    {"-s", "--save_output", SAVE_OUTPUT, "Save Output file\n"},
    {"--", "--save_chksum", SAVE_CHKSUM, "Save Check sum file\n"},
    {"-n", "--num_frames", NUM_FRAMES, "Number of frames to be decoded\n"},
    {"--", "--num_cores", NUM_CORES, "Number of cores to be used\n"},
    {"--", "--disable_deblock_level", DISABLE_DEBLOCK_LEVEL,
     "Disable deblocking level : 0 to 4 - 0 Enable deblocking 4 Disable "
     "deblocking completely\n"},
    {"--", "--u4_loopback", LOOPBACK, "Enable playback in a loop\n"},
    {"--", "--degrade_type", DEGRADE_TYPE,
     "Degrade type : 0: No degrade 0th bit set : Disable SAO 1st bit set : "
     "Disable deblocking 2nd bit set : Faster inter prediction filters 3rd bit "
     "set : Fastest inter prediction filters\n"},
    {"--", "--degrade_pics", DEGRADE_PICS,
     "Degrade pics : 0 : No degrade  1 : Only on non-reference frames  2 : Do "
     "not degrade every 4th or key frames  3 : All non-key frames  4 : All "
     "frames\n"},
    {"--", "--arch", ARCH,
     "Set Architecture. Supported values - ARM_A9Q, ARMV8_GENERIC, "
     "X86_GENERIC, X86_SSE4 \n"},
};

#if ANDROID_NDK
/*****************************************************************************/
/*                                                                           */
/*  Function Name : raise                                                    */
/*                                                                           */
/*  Description   : Needed as a workaround when the application is built in  */
/*                  Android NDK. This is an exception to be called for divide*/
/*                  by zero error                                            */
/*                                                                           */
/*  Inputs        : a                                                        */
/*  Globals       :                                                          */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       :                                                          */
/*  Returns       :                                                          */
/*                                                                           */
/*  Issues        :                                                          */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes                              */
/*         07 09 2012   100189          Initial Version                      */
/*                                                                           */
/*****************************************************************************/
static int raise(int a)
{
    printf("Divide by zero\n");
    return 0;
}
#endif

static void mvcd_print_usage(void)
{
    WORD32 i = 0;
    WORD32 num_entries = sizeof(argument_mapping) / sizeof(argument_t);

    printf("\nUsage:\n");

    while(i < num_entries)
    {
        printf("%-32s\t %s", argument_mapping[i].au1_argument_name,
               argument_mapping[i].au1_description);
        i++;
    }
}

static ARGUMENT_T mvcd_get_argument(char *name)
{
    WORD32 i = 0;
    WORD32 num_entries = sizeof(argument_mapping) / sizeof(argument_t);

    while(i < num_entries)
    {
        if((0 == strcmp((char *) argument_mapping[i].au1_argument_name, name)) ||
           ((0 == strcmp((char *) argument_mapping[i].au1_argument_shortname, name)) &&
            (0 != strcmp((char *) argument_mapping[i].au1_argument_shortname, "--"))))
        {
            return argument_mapping[i].e_argument;
        }
        i++;
    }

    return INVALID;
}

static void mvcd_exit(UWORD8 *pu1_err_message)
{
    printf("%s\n", pu1_err_message);
    exit(-1);
}

static void mvcd_parse_argument(mvc_dec_ctx_t *ps_app_ctx, char *argument, char *value)
{
    ARGUMENT_T arg;

    arg = mvcd_get_argument((char *) argument);

    switch(arg)
    {
        case HELP:
        {
            mvcd_print_usage();

            exit(-1);
        }
        case INPUT_FILE:
        {
            sscanf(value, "%s", ps_app_ctx->s_mvc_app_files.au1_ip_fname);

            break;
        }
        case OUTPUT:
        {
            sscanf(value, "%s", ps_app_ctx->s_mvc_app_files.au1_op_fname);

            break;
        }
        case CHKSUM:
        {
            sscanf(value, "%s", ps_app_ctx->s_mvc_app_files.au1_op_chksum_fname);

            break;
        }
        case SAVE_OUTPUT:
        {
            sscanf(value, "%u", &ps_app_ctx->u4_file_save_flag);

            break;
        }
        case SAVE_CHKSUM:
        {
            sscanf(value, "%u", &ps_app_ctx->u4_chksum_save_flag);

            break;
        }
        case NUM_FRAMES:
        {
            sscanf(value, "%u", &ps_app_ctx->u4_max_frm_ts);

            break;
        }
        case NUM_CORES:
        {
            sscanf(value, "%u", &ps_app_ctx->u4_num_cores);

            break;
        }
        case DEGRADE_PICS:
        {
            sscanf(value, "%d", &ps_app_ctx->i4_degrade_pics);

            break;
        }
        case DEGRADE_TYPE:
        {
            sscanf(value, "%d", &ps_app_ctx->i4_degrade_type);

            break;
        }
        case LOOPBACK:
        {
            sscanf(value, "%u", &ps_app_ctx->u4_loopback);

            break;
        }
        case ARCH:
        {
            if((strcmp(value, "ARM_A9Q")) == 0)
            {
                ps_app_ctx->e_arch = ARCH_ARM_A9Q;
            }
            else if((strcmp(value, "X86_GENERIC")) == 0)
            {
                ps_app_ctx->e_arch = ARCH_X86_GENERIC;
            }
            else if((strcmp(value, "X86_SSSE3")) == 0)
            {
                ps_app_ctx->e_arch = ARCH_X86_SSSE3;
            }
            else if((strcmp(value, "X86_SSE42")) == 0)
            {
                ps_app_ctx->e_arch = ARCH_X86_SSE42;
            }
            else if((strcmp(value, "ARMV8_GENERIC")) == 0)
            {
                ps_app_ctx->e_arch = ARCH_ARMV8_GENERIC;
            }
            else
            {
                printf("\nInvalid Arch. Setting it to ARCH_ARMV8_GENERIC\n");
                ps_app_ctx->e_arch = ARCH_ARMV8_GENERIC;
            }

            break;
        }
        case DISABLE_DEBLOCK_LEVEL:
        {
            sscanf(value, "%u", &ps_app_ctx->u4_disable_dblk_level);

            break;
        }
        default:
        {
            printf("Ignoring argument :  %s\n", argument);

            break;
        }
    }
}

static void mvcd_read_cfg_file(mvc_dec_ctx_t *ps_app_ctx, FILE *ps_cfg_file)
{
    char line[STRLENGTH];
    char description[STRLENGTH];
    char value[STRLENGTH];
    char argument[STRLENGTH];
    void *ret;

    while(0 == feof(ps_cfg_file))
    {
        line[0] = '\0';
        ret = fgets(line, STRLENGTH, ps_cfg_file);
        if(NULL == ret) break;
        argument[0] = '\0';
        /* Reading Input File Name */
        sscanf(line, "%s %s %s", argument, value, description);
        if(argument[0] == '\0') continue;

        mvcd_parse_argument(ps_app_ctx, argument, value);
    }
}

static void mvcd_get_view_file_name(const UWORD8 *pu1_default_name, UWORD8 *pu1_view_file_name,
                             UWORD16 u2_view_id)
{
    CHAR *apc_sub_str[2];
    CHAR ac_string[STRLENGTH];

    const CHAR ac_delimiters[] = ".";

    strcpy(ac_string, (char *) pu1_default_name);

    apc_sub_str[0] = strtok(ac_string, ac_delimiters);
    apc_sub_str[1] = strtok(NULL, ac_delimiters);

    ASSERT(NULL == strtok(NULL, ac_delimiters));
    ASSERT((strlen(apc_sub_str[0]) + strlen(apc_sub_str[1]) + 3) < STRLENGTH);

    sprintf((char *) pu1_view_file_name, "%s_%d.%s", apc_sub_str[0], u2_view_id, apc_sub_str[1]);
}

static IV_API_CALL_STATUS_T mvcd_in_buf_alloc(mvc_dec_ctx_t *ps_app_ctxt, UWORD8 **ppu1_bs_buf)
{
    ppu1_bs_buf[0] =
        (UWORD8 *) malloc(ps_app_ctxt->s_disp_buf_props.s_ivd_op.u4_min_in_buf_size[0]);

    if(ppu1_bs_buf[0] == NULL)
    {
        return IV_FAIL;
    }

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T mvcd_in_buf_free(UWORD8 *pu1_bs_buf)
{
    free(pu1_bs_buf);

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T mvcd_out_buf_alloc(mvc_dec_ctx_t *ps_app_ctxt, ivd_out_bufdesc_t *ps_out_buf)
{
    UWORD32 i;

    UWORD16 u2_num_views = ps_app_ctxt->s_disp_buf_props.s_mvc_buf_info.u2_num_views;

    if(ps_app_ctxt->s_disp_buf_props.s_ivd_op.u4_min_num_out_bufs < (NUM_COMPONENTS * u2_num_views))
    {
        return IV_FAIL;
    }

    ps_out_buf->u4_num_bufs = ps_app_ctxt->s_disp_buf_props.s_ivd_op.u4_min_num_out_bufs;

    for(i = 0; i < ps_out_buf->u4_num_bufs; i++)
    {
        ps_out_buf->u4_min_out_buf_size[i] =
            ps_app_ctxt->s_disp_buf_props.s_ivd_op.u4_min_out_buf_size[i];
        ps_out_buf->pu1_bufs[i] =
            (UWORD8 *) malloc(ps_app_ctxt->s_disp_buf_props.s_ivd_op.u4_min_out_buf_size[i]);

        if(ps_out_buf->pu1_bufs[i] == NULL)
        {
            return IV_FAIL;
        }
    }

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T mvcd_out_buf_free(ivd_out_bufdesc_t *ps_out_buf)
{
    UWORD32 i;

    for(i = 0; i < ps_out_buf->u4_num_bufs; i++)
    {
        free(ps_out_buf->pu1_bufs[i]);
    }

    return IV_SUCCESS;
}

static void mvcd_output_write_stall(UWORD8 *au1_fname, UWORD32 u4_cur_frm_idx)
{
    FILE *fp_fast_file = NULL;

    const UWORD8 threshold = 64;
    char past_fname[1000];

    if(u4_cur_frm_idx >= threshold)
    {
        sprintf(past_fname, (char *) au1_fname, u4_cur_frm_idx - threshold);
        do
        {
            fp_fast_file = fopen(past_fname, "rb");

            if(fp_fast_file != NULL)
            {
                fclose(fp_fast_file);
                /* Wait until the resource is released by a third party app*/
                ithread_sleep(5000);
            }
            else
            {
                break;
            }
        } while(1);
    }
}

static IV_API_CALL_STATUS_T mvcd_create_decoder(mvc_dec_ctx_t *ps_app_ctxt)
{
    imvcd_create_ip_t s_create_ip;
    imvcd_create_op_t s_create_op;

    IV_API_CALL_STATUS_T ret;

    s_create_ip.s_ivd_ip.e_cmd = IVD_CMD_CREATE;
    s_create_ip.s_ivd_ip.e_output_format = ps_app_ctxt->e_output_chroma_format;
    s_create_ip.s_ivd_ip.pf_aligned_alloc = mvcd_aligned_malloc;
    s_create_ip.s_ivd_ip.pf_aligned_free = mvcd_aligned_free;
    s_create_ip.s_ivd_ip.u4_share_disp_buf = 0;
    s_create_ip.s_ivd_ip.pv_mem_ctxt = NULL;

    s_create_ip.s_ivd_ip.u4_size = sizeof(s_create_ip.s_ivd_ip);
    s_create_op.s_ivd_op.u4_size = sizeof(s_create_op.s_ivd_op);

    ret = imvcd_api_function(NULL, &s_create_ip, &s_create_op);

    if(ret != IV_SUCCESS)
    {
        return ret;
    }

    ps_app_ctxt->ps_codec_obj = (iv_obj_t *) s_create_op.s_ivd_op.pv_handle;

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T mvcd_set_decode_mode(mvc_dec_ctx_t *ps_app_ctxt,
                                          IVD_VIDEO_DECODE_MODE_T e_decode_mode)
{
    imvcd_set_config_ip_t s_ctl_ip;
    imvcd_set_config_op_t s_ctl_op;

    s_ctl_ip.s_ivd_ip.u4_size = sizeof(s_ctl_ip.s_ivd_ip);
    s_ctl_op.s_ivd_op.u4_size = sizeof(s_ctl_op.s_ivd_op);
    s_ctl_ip.s_ivd_ip.e_cmd = IVD_CMD_VIDEO_CTL;
    s_ctl_ip.s_ivd_ip.e_sub_cmd = (WORD32) IVD_CMD_CTL_SETPARAMS;
    s_ctl_ip.s_ivd_ip.e_frm_out_mode = IVD_DISPLAY_FRAME_OUT;
    s_ctl_ip.s_ivd_ip.e_frm_skip_mode = IVD_SKIP_NONE;
    s_ctl_ip.s_ivd_ip.e_vid_dec_mode = e_decode_mode;

    return imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_ctl_ip, &s_ctl_op);
}

static IV_API_CALL_STATUS_T mvcd_set_num_cores(mvc_dec_ctx_t *ps_app_ctxt)
{
    imvcd_set_num_cores_ip_t s_ctl_ip;
    imvcd_set_num_cores_op_t s_ctl_op;

    s_ctl_ip.u4_size = sizeof(s_ctl_ip);
    s_ctl_op.u4_size = sizeof(s_ctl_op);
    s_ctl_ip.e_cmd = IVD_CMD_VIDEO_CTL;
    s_ctl_ip.e_sub_cmd = (WORD32) IMVCD_CTL_SET_NUM_CORES;
    s_ctl_ip.u4_num_cores = ps_app_ctxt->u4_num_cores;

    return imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_ctl_ip, &s_ctl_op);
}

static IV_API_CALL_STATUS_T mvcd_set_arch(mvc_dec_ctx_t *ps_app_ctxt)
{
    imvcd_set_arch_ip_t s_ctl_ip;
    imvcd_set_arch_op_t s_ctl_op;

    s_ctl_ip.u4_size = sizeof(s_ctl_ip);
    s_ctl_op.u4_size = sizeof(s_ctl_op);
    s_ctl_ip.e_cmd = IVD_CMD_VIDEO_CTL;
    s_ctl_ip.e_sub_cmd = (WORD32) IMVCD_CTL_SET_PROCESSOR;
    s_ctl_ip.e_arch = ps_app_ctxt->e_arch;
    s_ctl_ip.e_soc = ps_app_ctxt->e_soc;

    return imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_ctl_ip, &s_ctl_op);
}

static IV_API_CALL_STATUS_T mvcd_dump_output(mvc_dec_ctx_t *ps_app_ctxt,
                                             iv_yuv_buf_t *ps_view_disp_bufs, FILE **pps_op_file,
                                             FILE **pps_op_chksum_file, UWORD16 u2_num_views)
{
    UWORD32 i, j;

    UWORD32 u4_file_save = ps_app_ctxt->u4_file_save_flag;
    UWORD32 u4_chksum_save = ps_app_ctxt->u4_chksum_save_flag;

    if(!u4_file_save && !u4_chksum_save)
    {
        return IV_SUCCESS;
    }

    if(NULL == ps_view_disp_bufs->pv_y_buf)
    {
        return IV_FAIL;
    }

    for(i = 0; i < u2_num_views; i++)
    {
        iv_yuv_buf_t *ps_view_buf = &ps_view_disp_bufs[i];

        if(ps_app_ctxt->e_output_chroma_format == IV_YUV_420P)
        {
            if(u4_file_save)
            {
                UWORD8 *pu1_buf = (UWORD8 *) ps_view_buf->pv_y_buf;
                UWORD16 u2_width = ps_view_buf->u4_y_wd;
                UWORD16 u2_height = ps_view_buf->u4_y_ht;

                for(j = 0; j < u2_height; j++)
                {
                    fwrite(pu1_buf, 1, u2_width, pps_op_file[i]);

                    pu1_buf += ps_view_buf->u4_y_strd;
                }

                pu1_buf = (UWORD8 *) ps_view_buf->pv_u_buf;
                u2_width = ps_view_buf->u4_u_wd;
                u2_height = ps_view_buf->u4_u_ht;

                for(j = 0; j < u2_height; j++)
                {
                    fwrite(pu1_buf, 1, u2_width, pps_op_file[i]);

                    pu1_buf += ps_view_buf->u4_u_strd;
                }

                pu1_buf = (UWORD8 *) ps_view_buf->pv_v_buf;
                u2_width = ps_view_buf->u4_v_wd;
                u2_height = ps_view_buf->u4_v_ht;

                for(j = 0; j < u2_height; j++)
                {
                    fwrite(pu1_buf, 1, u2_width, pps_op_file[i]);

                    pu1_buf += ps_view_buf->u4_v_strd;
                }
            }

#ifndef MD5_DISABLE
            if(u4_chksum_save)
            {
                UWORD8 au1_cksum[16];

                UWORD8 *pu1_buf = (UWORD8 *) ps_view_buf->pv_y_buf;
                UWORD16 u2_width = ps_view_buf->u4_y_wd;
                UWORD16 u2_height = ps_view_buf->u4_y_ht;

                calc_md5_cksum(pu1_buf, ps_view_buf->u4_y_strd, u2_width, u2_height, au1_cksum);

                fwrite(au1_cksum, sizeof(UWORD8), 16, pps_op_chksum_file[i]);

                pu1_buf = (UWORD8 *) ps_view_buf->pv_u_buf;
                u2_width = ps_view_buf->u4_u_wd;
                u2_height = ps_view_buf->u4_u_ht;

                calc_md5_cksum(pu1_buf, ps_view_buf->u4_u_strd, u2_width, u2_height, au1_cksum);

                fwrite(au1_cksum, sizeof(UWORD8), 16, pps_op_chksum_file[i]);

                pu1_buf = (UWORD8 *) ps_view_buf->pv_v_buf;
                u2_width = ps_view_buf->u4_v_wd;
                u2_height = ps_view_buf->u4_v_ht;

                calc_md5_cksum(pu1_buf, ps_view_buf->u4_v_strd, u2_width, u2_height, au1_cksum);

                fwrite(au1_cksum, sizeof(UWORD8), 16, pps_op_chksum_file[i]);
            }
#endif
        }
        else
        {
            return IV_FAIL;
        }

        fflush(pps_op_file[i]);
        fflush(pps_op_chksum_file[i]);
    }

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T mvcd_flush(mvc_dec_ctx_t *ps_app_ctxt, ivd_out_bufdesc_t *ps_out_buf,
                                FILE **pps_op_file, FILE **pps_op_chksum_file, UWORD8 *pu1_bs_buf,
                                UWORD32 *pu4_op_frm_ts, UWORD32 u4_ip_frm_ts,
                                UWORD32 u4_bytes_remaining, UWORD16 u2_num_views)
{
    IV_API_CALL_STATUS_T ret = IV_SUCCESS;

    do
    {
        imvcd_flush_dec_ip_t s_ctl_ip;
        imvcd_flush_dec_op_t s_ctl_op;

        if(*(pu4_op_frm_ts) >= u4_ip_frm_ts)
        {
            break;
        }

        s_ctl_ip.s_ivd_ip.u4_size = sizeof(s_ctl_ip.s_ivd_ip);
        s_ctl_op.s_ivd_op.u4_size = sizeof(s_ctl_op.s_ivd_op);
        s_ctl_ip.s_ivd_ip.e_cmd = IVD_CMD_VIDEO_CTL;
        s_ctl_ip.s_ivd_ip.e_sub_cmd = (WORD32) IVD_CMD_CTL_FLUSH;

        ret = imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_ctl_ip, &s_ctl_op);

        if(IV_SUCCESS == ret)
        {
            imvcd_video_decode_ip_t s_video_decode_ip;
            imvcd_video_decode_op_t s_video_decode_op;

            s_video_decode_ip.s_ivd_ip.e_cmd = IVD_CMD_VIDEO_DECODE;
            s_video_decode_ip.s_ivd_ip.u4_ts = u4_ip_frm_ts;
            s_video_decode_ip.s_ivd_ip.pv_stream_buffer = pu1_bs_buf;
            s_video_decode_ip.s_ivd_ip.u4_num_Bytes = u4_bytes_remaining;
            s_video_decode_ip.s_ivd_ip.s_out_buffer = ps_out_buf[0];
            s_video_decode_op.ps_view_disp_bufs = ps_app_ctxt->as_view_disp_bufs;

            s_video_decode_ip.s_ivd_ip.u4_size = sizeof(s_video_decode_ip.s_ivd_ip);
            s_video_decode_op.s_ivd_op.u4_size = sizeof(s_video_decode_op.s_ivd_op);

            ret = imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_video_decode_ip,
                                     &s_video_decode_op);

            if(1 == s_video_decode_op.s_ivd_op.u4_output_present)
            {
                ret = mvcd_dump_output(ps_app_ctxt, s_video_decode_op.ps_view_disp_bufs,
                                       pps_op_file, pps_op_chksum_file, u2_num_views);

                (*pu4_op_frm_ts)++;
            }
            else
            {
                break;
            }
        }
    } while(IV_SUCCESS == ret);

    return ret;
}

static IV_API_CALL_STATUS_T mvcd_decode_header(mvc_dec_ctx_t *ps_app_ctxt, ivd_out_bufdesc_t *ps_out_buf,
                                        FILE *ps_ip_file, IVD_ERROR_CODES_T *pe_error_code,
                                        UWORD32 *pu4_num_bytes_dec, UWORD32 u4_ip_frm_ts)
{
    IV_API_CALL_STATUS_T ret;

    UWORD32 u4_ip_buf_len;
    UWORD8 *pu1_bs_buf;

    imvcd_video_decode_ip_t s_video_decode_ip = {0};
    imvcd_video_decode_op_t s_video_decode_op = {0};

    UWORD32 u4_file_pos = 0;
    UWORD32 u4_num_bytes_dec = 0;

    ret = mvcd_set_decode_mode(ps_app_ctxt, IVD_DECODE_HEADER);

    if(ret != IV_SUCCESS)
    {
        pe_error_code[0] = IVD_INIT_DEC_FAILED;

        return ret;
    }

    /* Allocate input buffer for header */
    u4_ip_buf_len = 256 * 1024;
    pu1_bs_buf = (UWORD8 *) malloc(u4_ip_buf_len);

    if(pu1_bs_buf == NULL)
    {
        pe_error_code[0] = IVD_MEM_ALLOC_FAILED;

        return IV_FAIL;
    }

    do
    {
        WORD32 u4_numbytes;
        UWORD32 u4_bytes_remaining;

        fseek(ps_ip_file, u4_file_pos, SEEK_SET);
        u4_numbytes = u4_ip_buf_len;

        u4_bytes_remaining = fread(pu1_bs_buf, sizeof(UWORD8), u4_numbytes, ps_ip_file);

        if(0 == u4_bytes_remaining)
        {
            pe_error_code[0] = IVD_UNEXPECTED_END_OF_STREAM;

            return IV_FAIL;
        }

        s_video_decode_ip.s_ivd_ip.e_cmd = IVD_CMD_VIDEO_DECODE;
        s_video_decode_ip.s_ivd_ip.u4_ts = u4_ip_frm_ts;
        s_video_decode_ip.s_ivd_ip.pv_stream_buffer = pu1_bs_buf;
        s_video_decode_ip.s_ivd_ip.u4_num_Bytes = u4_bytes_remaining;
        s_video_decode_ip.s_ivd_ip.s_out_buffer = ps_out_buf[0];
        s_video_decode_op.ps_view_disp_bufs = ps_app_ctxt->as_view_disp_bufs;

        s_video_decode_ip.s_ivd_ip.u4_size = sizeof(s_video_decode_ip.s_ivd_ip);
        s_video_decode_op.s_ivd_op.u4_size = sizeof(s_video_decode_op.s_ivd_op);

        ret = imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_video_decode_ip, &s_video_decode_op);

        if(ret != IV_SUCCESS)
        {
            pe_error_code[0] = s_video_decode_op.s_ivd_op.u4_error_code;

            return ret;
        }

        u4_num_bytes_dec += s_video_decode_op.s_ivd_op.u4_num_bytes_consumed;

#ifndef PROFILE_ENABLE
        printf("%d\n", s_video_decode_op.s_ivd_op.u4_num_bytes_consumed);
#endif
    } while(ret != IV_SUCCESS);

    ps_app_ctxt->u4_pic_wd = s_video_decode_op.s_ivd_op.u4_pic_wd;
    ps_app_ctxt->u4_pic_ht = s_video_decode_op.s_ivd_op.u4_pic_ht;

    pu4_num_bytes_dec[0] = u4_num_bytes_dec;
    pe_error_code[0] = IVD_ERROR_NONE;

    free(pu1_bs_buf);

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T mvcd_get_buf_info(mvc_dec_ctx_t *ps_app_ctxt)
{
    imvcd_get_buf_info_ip_t s_ctl_ip;
    imvcd_get_buf_info_op_t s_ctl_op;

    IV_API_CALL_STATUS_T e_retval;

    s_ctl_ip.s_ivd_ip.u4_size = sizeof(s_ctl_ip.s_ivd_ip);
    s_ctl_op.s_ivd_op.u4_size = sizeof(s_ctl_op.s_ivd_op);
    s_ctl_ip.s_ivd_ip.e_cmd = IVD_CMD_VIDEO_CTL;
    s_ctl_ip.s_ivd_ip.e_sub_cmd = (WORD32) IVD_CMD_CTL_GETBUFINFO;

    e_retval = imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_ctl_ip, &s_ctl_op);

    ps_app_ctxt->s_disp_buf_props = s_ctl_op;

    return e_retval;
}

static IV_API_CALL_STATUS_T mvcd_set_degrade_type(mvc_dec_ctx_t *ps_app_ctxt)
{
    imvcd_set_degrade_mode_ip_t s_ctl_ip;
    imvcd_set_degrade_mode_op_t s_ctl_op;

    s_ctl_ip.u4_size = sizeof(s_ctl_ip);
    s_ctl_op.u4_size = sizeof(s_ctl_op);
    s_ctl_ip.e_cmd = IVD_CMD_VIDEO_CTL;
    s_ctl_ip.e_sub_cmd = (WORD32) IMVCD_CTL_DEGRADE;
    s_ctl_ip.i4_degrade_type = ps_app_ctxt->i4_degrade_type;
    s_ctl_ip.i4_degrade_pics = ps_app_ctxt->i4_degrade_pics;
    s_ctl_ip.i4_nondegrade_interval = DEFAULT_NON_DEGRADE_INTERVAL;

    return imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_ctl_ip, &s_ctl_op);
}

static IV_API_CALL_STATUS_T mvcd_decode_frame(mvc_dec_ctx_t *ps_app_ctxt, ivd_out_bufdesc_t *ps_out_buf,
                                       FILE *ps_ip_file, FILE **pps_op_file,
                                       FILE **pps_op_chksum_file, UWORD8 *pu1_bs_buf,
                                       IVD_ERROR_CODES_T *pe_error_code,
                                       UWORD32 *pu4_num_bytes_consumed, UWORD32 *pu4_file_pos,
                                       UWORD32 *pu4_ip_frm_ts, UWORD32 *pu4_op_frm_ts,
                                       UWORD32 u4_ip_buf_len)
{
    imvcd_video_decode_ip_t s_video_decode_ip;
    imvcd_video_decode_op_t s_video_decode_op;
#ifdef PROFILE_ENABLE
    TIMER s_start_timer;
    TIMER s_end_timer;
#ifdef WINDOWS_TIMER
    TIMER frequency;
#endif

    UWORD32 au4_peak_window[PEAK_WINDOW_SIZE];
    UWORD32 s_elapsed_time;
#endif
    UWORD32 u4_max_op_frm_ts;
    UWORD32 u4_bytes_remaining;
    UWORD32 u4_numbytes;

    IV_API_CALL_STATUS_T ret;

    UWORD32 u4_num_bytes_dec = 0;
    UWORD16 u2_num_views = ps_app_ctxt->s_disp_buf_props.s_mvc_buf_info.u2_num_views;
#ifdef PROFILE_ENABLE
    UWORD32 u4_frm_cnt = 0;
    UWORD32 u4_tot_cycles = 0;
    UWORD32 u4_tot_fmt_cycles = 0;
    UWORD32 u4_peak_window_idx = 0;
    UWORD32 u4_peak_avg_max = 0;
#endif

#ifdef PROFILE_ENABLE
    memset(au4_peak_window, 0, sizeof(WORD32) * PEAK_WINDOW_SIZE);

#ifdef WINDOWS_TIMER
    QueryPerformanceFrequency(&frequency);
#endif
#endif

    /*************************************************************************/
    /* Set the decoder in frame decode mode. It was set in header decode     */
    /* mode earlier                                                          */
    /*************************************************************************/
    ret = mvcd_set_decode_mode(ps_app_ctxt, IVD_DECODE_FRAME);

    if(IV_SUCCESS != ret)
    {
        pe_error_code[0] = IVD_INIT_DEC_FAILED;

        return ret;
    }

    /*************************************************************************/
    /* If required disable deblocking and sao at given level                 */
    /*************************************************************************/
    ret = mvcd_set_degrade_type(ps_app_ctxt);

    if(IV_SUCCESS != ret)
    {
        pe_error_code[0] = IVD_INIT_DEC_FAILED;

        return ret;
    }

    u4_max_op_frm_ts = ps_app_ctxt->u4_max_frm_ts + ps_app_ctxt->u4_disp_delay;

    if(u4_max_op_frm_ts < ps_app_ctxt->u4_disp_delay)
    {
        /* clip as overflow has occured*/
        u4_max_op_frm_ts = UINT32_MAX;
    }

    if(IV_SUCCESS != ret)
    {
        pe_error_code[0] = IVD_INIT_DEC_FAILED;

        return ret;
    }

    while(pu4_op_frm_ts[0] < u4_max_op_frm_ts)
    {
        fseek(ps_ip_file, pu4_file_pos[0], SEEK_SET);
        u4_numbytes = u4_ip_buf_len;

        u4_bytes_remaining = fread(pu1_bs_buf, sizeof(UWORD8), u4_numbytes, ps_ip_file);

        if(0 == u4_bytes_remaining)
        {
            if(ps_app_ctxt->u4_loopback)
            {
                pu4_file_pos[0] = 0;

                fseek(ps_ip_file, pu4_file_pos[0], SEEK_SET);
                u4_numbytes = u4_ip_buf_len;

                u4_bytes_remaining = fread(pu1_bs_buf, sizeof(UWORD8), u4_numbytes, ps_ip_file);
            }
            else
            {
                break;
            }
        }

        s_video_decode_ip.s_ivd_ip.e_cmd = IVD_CMD_VIDEO_DECODE;
        s_video_decode_ip.s_ivd_ip.u4_ts = pu4_ip_frm_ts[0];
        s_video_decode_ip.s_ivd_ip.pv_stream_buffer = pu1_bs_buf;
        s_video_decode_ip.s_ivd_ip.u4_num_Bytes = u4_bytes_remaining;
        s_video_decode_ip.s_ivd_ip.s_out_buffer = ps_out_buf[0];
        s_video_decode_op.ps_view_disp_bufs = ps_app_ctxt->as_view_disp_bufs;

        s_video_decode_ip.s_ivd_ip.u4_size = sizeof(s_video_decode_ip.s_ivd_ip);
        s_video_decode_op.s_ivd_op.u4_size = sizeof(s_video_decode_op.s_ivd_op);

        /****************/
        /* Video Decode */
        /****************/
        GETTIME(&s_start_timer);

        ret = imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_video_decode_ip, &s_video_decode_op);

        GETTIME(&s_end_timer);
        ELAPSEDTIME(s_start_timer, s_end_timer, s_elapsed_time, frequency);

#ifdef PROFILE_ENABLE
        {
            UWORD32 peak_avg, id;

            u4_tot_cycles += s_elapsed_time;
            au4_peak_window[u4_peak_window_idx++] = s_elapsed_time;

            if(u4_peak_window_idx == PEAK_WINDOW_SIZE)
            {
                u4_peak_window_idx = 0;
            }

            peak_avg = 0;
            for(id = 0; id < PEAK_WINDOW_SIZE; id++)
            {
                peak_avg += au4_peak_window[id];
            }

            peak_avg /= PEAK_WINDOW_SIZE;
            if(peak_avg > u4_peak_avg_max) u4_peak_avg_max = peak_avg;
            u4_frm_cnt++;

            printf(
                "FrameNum: %4d TimeTaken(microsec): %6d AvgTime: %6d "
                "PeakAvgTimeMax: "
                "%6d Output: %2d NumBytes: %6d \n",
                u4_frm_cnt, s_elapsed_time, u4_tot_cycles / u4_frm_cnt, u4_peak_avg_max,
                s_video_decode_op.s_ivd_op.u4_output_present,
                s_video_decode_op.s_ivd_op.u4_num_bytes_consumed);
        }
#else
        printf("%d\n", s_video_decode_op.s_ivd_op.u4_num_bytes_consumed);
#endif

        if(ret != IV_SUCCESS)
        {
            pe_error_code[0] = s_video_decode_op.s_ivd_op.u4_error_code;

            return ret;
        }

        u4_num_bytes_dec = s_video_decode_op.s_ivd_op.u4_num_bytes_consumed;

        pu4_file_pos[0] += u4_num_bytes_dec;
        pu4_num_bytes_consumed[0] += u4_num_bytes_dec;
        pu4_ip_frm_ts[0]++;

        if(s_video_decode_op.s_ivd_op.u4_output_present)
        {
            char cur_fname[1000];

            char *extn = NULL;

            /* The objective is to dump the decoded frames into separate files
             * instead of dumping all the frames in one common file. Also, the
             * number of dumped frames at any given instance of time cannot exceed
             * 'frame_memory'
             */
            if(ps_app_ctxt->u4_file_save_flag)
            {
                /* Locate the position of extension yuv */
                extn = strstr((char *) ps_app_ctxt->s_mvc_app_files.au1_op_fname, "%d");

                if(extn != NULL)
                {
                    mvcd_output_write_stall(ps_app_ctxt->s_mvc_app_files.au1_op_fname,
                                            pu4_op_frm_ts[0]);

                    sprintf(cur_fname, (char *) ps_app_ctxt->s_mvc_app_files.au1_op_fname,
                            pu4_op_frm_ts[0]);

                    pps_op_file[0] = fopen(cur_fname, "wb");

                    if(NULL == pps_op_file[0])
                    {
                        pe_error_code[0] = IVD_MEM_ALLOC_FAILED;

                        return IV_FAIL;
                    }
                }
            }

            if(u2_num_views > 1)
            {
                if(s_video_decode_op.ps_view_disp_bufs[0].u4_y_wd !=
                   s_video_decode_op.ps_view_disp_bufs[1].u4_y_wd)
                {
                    pe_error_code[0] = IVD_MEM_ALLOC_FAILED;

                    return IV_FAIL;
                }

                if(s_video_decode_op.ps_view_disp_bufs[0].u4_y_ht !=
                   s_video_decode_op.ps_view_disp_bufs[1].u4_y_ht)
                {
                    pe_error_code[0] = IVD_MEM_ALLOC_FAILED;

                    return IV_FAIL;
                }
            }

            mvcd_dump_output(ps_app_ctxt, s_video_decode_op.ps_view_disp_bufs, pps_op_file,
                             pps_op_chksum_file, u2_num_views);

            pu4_op_frm_ts[0]++;
        }
        else if((s_video_decode_op.s_ivd_op.u4_error_code >> IVD_FATALERROR) & 1)
        {
            pe_error_code[0] = s_video_decode_op.s_ivd_op.u4_error_code;

            return IV_FAIL;
        }
    }

#ifdef PROFILE_ENABLE
    printf("Summary\n");
    printf("Input filename                  : %s\n", ps_app_ctxt->s_mvc_app_files.au1_ip_fname);
    printf("Output Width                    : %-4d\n", ps_app_ctxt->u4_pic_wd);
    printf("Output Height                   : %-4d\n", ps_app_ctxt->u4_pic_ht);

    if(u4_frm_cnt)
    {
        double avg = u4_tot_cycles / u4_frm_cnt;
        double bytes_avg = pu4_num_bytes_consumed[0] / u4_frm_cnt;
        double bitrate = (bytes_avg * 8 * ps_app_ctxt->u4_fps) / 1000000;

        printf("Bitrate @ %2d u4_fps(mbps)          : %-6.2f\n", ps_app_ctxt->u4_fps, bitrate);
        printf("Average decode time(micro sec)  : %-6d\n", (WORD32) avg);
        printf("Avg Peak decode time(%2d frames) : %-6d\n", PEAK_WINDOW_SIZE,
               (WORD32) u4_peak_avg_max);

        avg = (u4_tot_cycles + u4_tot_fmt_cycles) * 1.0 / u4_frm_cnt;
        printf("FPS achieved (with format conv) : %-3.2f\n", 1000000 / avg);
    }
#endif

    return IV_SUCCESS;
}

static IV_API_CALL_STATUS_T mvcd_delete_decoder(mvc_dec_ctx_t *ps_app_ctxt)
{
    imvcd_delete_ip_t s_delete_ip;
    imvcd_delete_op_t s_delete_op;

    s_delete_ip.s_ivd_ip.e_cmd = IVD_CMD_DELETE;

    s_delete_ip.s_ivd_ip.u4_size = sizeof(s_delete_ip.s_ivd_ip);
    s_delete_op.s_ivd_op.u4_size = sizeof(s_delete_op.s_ivd_op);

    return imvcd_api_function(ps_app_ctxt->ps_codec_obj, &s_delete_ip, &s_delete_op);
}

int main(WORD32 argc, char *argv[])
{
    mvc_dec_ctx_t s_app_ctxt;
    ivd_out_bufdesc_t s_out_buf;

    IV_API_CALL_STATUS_T ret;
    IVD_ERROR_CODES_T e_error_code;

    UWORD8 au1_cfg_fname[STRLENGTH];
    UWORD8 au1_error_str[STRLENGTH];
    UWORD32 i;
    UWORD32 u4_ip_buf_len;
    UWORD32 u4_total_bytes_consumed;
    UWORD8 au1_view_file_names[STRLENGTH];
    UWORD16 u2_num_views;
    UWORD32 u4_num_header_bytes;

    FILE *ps_cfg_file = NULL;
    FILE *ps_ip_file = NULL;
    FILE *aps_op_file[MAX_NUM_VIEWS] = {NULL};
    FILE *aps_op_chksum_file[MAX_NUM_VIEWS] = {NULL};

    UWORD8 *pu1_bs_buf = NULL;
    UWORD32 u4_file_pos = 0;
    UWORD32 u4_ip_frm_ts = 0, u4_op_frm_ts = 0;
    UWORD32 u4_bytes_remaining = 0;

    /* Usage */
    if(argc < 2)
    {
        mvcd_print_usage();

        exit(-1);
    }
    else if(argc == 2)
    {
        strcpy((char *) au1_cfg_fname, argv[1]);
    }

    /***********************************************************************/
    /*                  Initialize Application parameters                  */
    /***********************************************************************/
    strcpy((char *) s_app_ctxt.s_mvc_app_files.au1_ip_fname, "\0");
    s_app_ctxt.u4_disp_delay = 0;
    s_app_ctxt.u4_max_frm_ts = 100;
    s_app_ctxt.u4_loopback = 0;
    u4_file_pos = 0;
    u4_total_bytes_consumed = 0;
    u4_ip_frm_ts = 0;
    u4_op_frm_ts = 0;

    s_app_ctxt.u4_num_cores = DEFAULT_NUM_CORES;
    s_app_ctxt.i4_degrade_type = 0;
    s_app_ctxt.i4_degrade_pics = 0;
    s_app_ctxt.e_arch = ARCH_X86_SSE42;
    s_app_ctxt.e_soc = SOC_GENERIC;
    s_app_ctxt.u1_quit = 0;
    s_app_ctxt.u4_disable_dblk_level = 0;
    s_app_ctxt.u4_file_save_flag = 0;
    s_app_ctxt.u4_chksum_save_flag = 0;
    s_app_ctxt.e_output_chroma_format = DEFAULT_COLOR_FORMAT;

    /*************************************************************************/
    /* Parse arguments                                                       */
    /*************************************************************************/

    /* Read command line arguments */
    if(argc > 2)
    {
        for(i = 1; i < (UWORD32) argc; i += 2)
        {
            if(CONFIG == mvcd_get_argument(argv[i]))
            {
                strcpy((char *) au1_cfg_fname, argv[i + 1]);

                if((ps_cfg_file = fopen((char *) au1_cfg_fname, "r")) == NULL)
                {
                    sprintf((char *) au1_error_str, "Could not open Configuration file %s",
                            au1_cfg_fname);

                    mvcd_exit(au1_error_str);
                }
                mvcd_read_cfg_file(&s_app_ctxt, ps_cfg_file);
                fclose(ps_cfg_file);
            }
            else
            {
                mvcd_parse_argument(&s_app_ctxt, argv[i], argv[i + 1]);
            }
        }
    }
    else
    {
        if((ps_cfg_file = fopen((char *) au1_cfg_fname, "r")) == NULL)
        {
            sprintf((char *) au1_error_str, "Could not open Configuration file %s", au1_cfg_fname);
            mvcd_exit(au1_error_str);
        }

        mvcd_read_cfg_file(&s_app_ctxt, ps_cfg_file);

        fclose(ps_cfg_file);
    }

    if(strcmp((char *) s_app_ctxt.s_mvc_app_files.au1_ip_fname, "\0") == 0)
    {
        printf("\nNo input file given for decoding\n");

        exit(-1);
    }

    /***********************************************************************/
    /*          create the file object for input file                      */
    /***********************************************************************/
    ps_ip_file = fopen((char *) s_app_ctxt.s_mvc_app_files.au1_ip_fname, "rb");

    if(NULL == ps_ip_file)
    {
        sprintf((char *) au1_error_str, "Could not open input file %s",
                s_app_ctxt.s_mvc_app_files.au1_ip_fname);

        mvcd_exit(au1_error_str);
    }

    /***********************************************************************/
    /*          create the file object for output file                     */
    /***********************************************************************/
    /* If the filename does not contain %d, then output will be dumped to
       a single file and it is opened here */
    if((1 == s_app_ctxt.u4_file_save_flag) &&
       (strstr((char *) s_app_ctxt.s_mvc_app_files.au1_op_fname, "%d") == NULL))
    {
        mvcd_get_view_file_name(s_app_ctxt.s_mvc_app_files.au1_op_fname, au1_view_file_names, 0);

        aps_op_file[0] = fopen((char *) au1_view_file_names, "wb");

        if(NULL == aps_op_file[0])
        {
            const CHAR au1_explanatory_string[] = "Could not open output file ";
            UWORD8 au1_error_str[STRLENGTH + sizeof(au1_explanatory_string) + 1];

            sprintf((char *) au1_error_str, "%s%s", au1_explanatory_string, au1_view_file_names);
            mvcd_exit(au1_error_str);
        }
    }

    /***********************************************************************/
    /*          create the file object for check sum file                  */
    /***********************************************************************/
    if((1 == s_app_ctxt.u4_chksum_save_flag) &&
       (strstr((char *) s_app_ctxt.s_mvc_app_files.au1_op_chksum_fname, "%d") == NULL))
    {
        mvcd_get_view_file_name(s_app_ctxt.s_mvc_app_files.au1_op_chksum_fname, au1_view_file_names,
                                0);
        aps_op_chksum_file[0] = fopen((char *) au1_view_file_names, "wb");

        if(NULL == aps_op_chksum_file[0])
        {
            const CHAR au1_explanatory_string[] = "Could not open check sum file ";
            UWORD8 au1_error_str[STRLENGTH + sizeof(au1_explanatory_string) + 1];

            sprintf((char *) au1_error_str, "%s%s", au1_explanatory_string, au1_view_file_names);
            mvcd_exit(au1_error_str);
        }
    }

    /***************************/
    /* Create decoder instance */
    /***************************/
    ret = mvcd_create_decoder(&s_app_ctxt);

    if(ret != IV_SUCCESS)
    {
        sprintf((char *) au1_error_str, "Error in Create %8x\n", ret);

        mvcd_exit(au1_error_str);
    }

    /**************/
    /* set params */
    /**************/
    ret = mvcd_set_decode_mode(&s_app_ctxt, IVD_DECODE_HEADER);

    if(ret != IV_SUCCESS)
    {
        sprintf((char *) au1_error_str, "\nError in Dec mode");
        mvcd_exit(au1_error_str);
    }

    /*************************************************************************/
    /* set num of cores                                                      */
    /*************************************************************************/
    ret = mvcd_set_num_cores(&s_app_ctxt);

    if(ret != IV_SUCCESS)
    {
        sprintf((char *) au1_error_str, "\nError in setting number of cores");
        mvcd_exit(au1_error_str);
    }

    /*************************************************************************/
    /* set processsor                                                        */
    /*************************************************************************/
    ret = mvcd_set_arch(&s_app_ctxt);

    if(ret != IV_SUCCESS)
    {
        sprintf((char *) au1_error_str, "\nError in setting Processor type");
        mvcd_exit(au1_error_str);
    }

    /*****************/
    /* Header Decode */
    /*****************/
    ret = mvcd_decode_header(&s_app_ctxt, &s_out_buf, ps_ip_file, &e_error_code,
                             &u4_num_header_bytes, u4_ip_frm_ts);

    if(ret != IV_SUCCESS)
    {
        sprintf((char *) au1_error_str, "Error in header decode 0x%x\n", e_error_code);
        mvcd_exit(au1_error_str);
    }

    u4_file_pos += u4_num_header_bytes;
    u4_total_bytes_consumed += u4_num_header_bytes;

    /****************/
    /* Get buf info */
    /****************/
    ret = mvcd_get_buf_info(&s_app_ctxt);

    if(ret != IV_SUCCESS)
    {
        sprintf((char *) au1_error_str, "Error in Get Buf Info %x",
                s_app_ctxt.s_disp_buf_props.s_ivd_op.u4_error_code);
        mvcd_exit(au1_error_str);
    }

    u2_num_views = s_app_ctxt.s_disp_buf_props.s_mvc_buf_info.u2_num_views;
    ASSERT(u2_num_views <= MAX_NUM_VIEWS);

    ret = mvcd_out_buf_alloc(&s_app_ctxt, &s_out_buf);

    if(IV_SUCCESS != ret)
    {
        sprintf((char *) au1_error_str, "Error in Out buf alloc\n");
        mvcd_exit(au1_error_str);
    }

    /***************************************************************************/
    /*   create the file object for output file for other views(if present) */
    /***************************************************************************/
    for(i = 1; i < u2_num_views; i++)
    {
        if((1 == s_app_ctxt.u4_file_save_flag) &&
           (strstr((char *) s_app_ctxt.s_mvc_app_files.au1_op_fname, "%d") == NULL))
        {
            mvcd_get_view_file_name(s_app_ctxt.s_mvc_app_files.au1_op_fname, au1_view_file_names,
                                    i);

            aps_op_file[i] = fopen((char *) au1_view_file_names, "wb");

            if(NULL == aps_op_file[i])
            {
                const UWORD8 au1_explanatory_string[] = "Could not open output file ";
                UWORD8 au1_error_str[STRLENGTH + sizeof(au1_explanatory_string) + 1];

                sprintf((char *) au1_error_str, "%s%s", au1_explanatory_string,
                        au1_view_file_names);
                mvcd_exit(au1_error_str);
            }
        }

        if((1 == s_app_ctxt.u4_chksum_save_flag) &&
           (strstr((char *) s_app_ctxt.s_mvc_app_files.au1_op_chksum_fname, "%d") == NULL))
        {
            mvcd_get_view_file_name(s_app_ctxt.s_mvc_app_files.au1_op_chksum_fname,
                                    au1_view_file_names, i);
            aps_op_chksum_file[i] = fopen((char *) au1_view_file_names, "wb");

            if(NULL == aps_op_chksum_file[i])
            {
                const UWORD8 au1_explanatory_string[] = "Could not open check sum file ";
                UWORD8 au1_error_str[STRLENGTH + sizeof(au1_explanatory_string) + 1];

                sprintf((char *) au1_error_str, "%s%s", au1_explanatory_string,
                        au1_view_file_names);
                mvcd_exit(au1_error_str);
            }
        }
    }

    u4_ip_buf_len = s_app_ctxt.s_disp_buf_props.s_ivd_op.u4_min_in_buf_size[0];
    ret = mvcd_in_buf_alloc(&s_app_ctxt, &pu1_bs_buf);

    if(IV_SUCCESS != ret)
    {
        sprintf((char *) au1_error_str, "Error in In buf alloc\n");
        mvcd_exit(au1_error_str);
    }

    ret = mvcd_decode_frame(&s_app_ctxt, &s_out_buf, ps_ip_file, aps_op_file, aps_op_chksum_file,
                            pu1_bs_buf, &e_error_code, &u4_total_bytes_consumed, &u4_file_pos,
                            &u4_ip_frm_ts, &u4_op_frm_ts, u4_ip_buf_len);

    if(ret != IV_SUCCESS)
    {
        sprintf((char *) au1_error_str, "Error in Frame decode 0x%x\n", e_error_code);
        mvcd_exit(au1_error_str);
    }

    /***********************************************************************/
    /*      To get the last decoded frames, call process with NULL input    */
    /***********************************************************************/
    ret = mvcd_flush(&s_app_ctxt, &s_out_buf, aps_op_file, aps_op_chksum_file, pu1_bs_buf,
                     &u4_op_frm_ts, u4_ip_frm_ts, u4_bytes_remaining, u2_num_views);

    if(ret != IV_SUCCESS)
    {
        sprintf((char *) au1_error_str, "\nError in Setting the decoder in flush mode");
        mvcd_exit(au1_error_str);
    }

    s_app_ctxt.u1_quit = 1;

    ret = mvcd_delete_decoder(&s_app_ctxt);

    if(IV_SUCCESS != ret)
    {
        sprintf((char *) au1_error_str, "Error in Codec delete");
        mvcd_exit(au1_error_str);
    }

    /***********************************************************************/
    /*              Close all the files and free all the memory            */
    /***********************************************************************/
    fclose(ps_ip_file);

    for(i = 0; i < u2_num_views; i++)
    {
        if(aps_op_file[i])
        {
            fclose(aps_op_file[i]);
        }

        if(aps_op_chksum_file[i])
        {
            fclose(aps_op_chksum_file[i]);
        }
    }

    mvcd_in_buf_free(pu1_bs_buf);

    mvcd_out_buf_free(&s_out_buf);

    return (0);
}
