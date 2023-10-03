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
*  recon.c
*
* @brief
*  Contains functions necessary for managing recon buffers
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/* User include files */
#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "ih264e.h"
#include "app.h"

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/


/*****************************************************************************/
/*  Macros                                                                   */
/*****************************************************************************/


/*****************************************************************************/
/*  Function Definitions                                                     */
/*****************************************************************************/

/**
**************************************************************************
* @brief Write recon buffers to a file
**************************************************************************
*/
IV_STATUS_T write_recon(FILE *fp, iv_raw_buf_t *ps_raw_buf)
{
    WORD32 num_comp = 2;
    WORD32 comp_idx;

    if(IV_YUV_420P == ps_raw_buf->e_color_fmt)
        num_comp = 3;

    for(comp_idx = 0; comp_idx < num_comp; comp_idx++)
    {
        WORD32 wd = ps_raw_buf->au4_wd[comp_idx];
        WORD32 ht = ps_raw_buf->au4_ht[comp_idx];
        UWORD8 *pu1_buf = ps_raw_buf->apv_bufs[comp_idx];
        WORD32 i;

        for(i = 0; i < ht; i++)
        {
            WORD32 bytes = fwrite(pu1_buf, sizeof(UWORD8), wd, fp);
            if(bytes != wd)
            {
                return (IV_FAIL);
            }
            pu1_buf += wd;
        }
    }
    fflush(fp);
    return IV_SUCCESS;
}

/**
**************************************************************************
* @brief Allocate space for recon buffers
**************************************************************************
*/
void allocate_recon(app_ctxt_t *ps_app_ctxt)
{
    WORD32 num_bufs = DEFAULT_NUM_RECON_BUFS;
    /* Size of buffer for YUV420 */
    WORD32 luma_size = ps_app_ctxt->u4_max_wd * ps_app_ctxt->u4_max_ht;
    WORD32 chroma_size = (luma_size) / 4;
    WORD32 pic_size = luma_size + chroma_size * 2;
    WORD32 i;

    for(i = 0; i < num_bufs; i++)
    {
        UWORD8 *pu1_buf = (UWORD8 *)ih264a_aligned_malloc(16, pic_size);

        if(NULL == pu1_buf)
        {
            CHAR ac_error[STRLENGTH];
            sprintf(ac_error, "Allocation failed for recon buffer of size %d\n",
                    pic_size);
            codec_exit(ac_error);
        }
        ps_app_ctxt->as_recon_buf[i].pu1_buf = pu1_buf;
        ps_app_ctxt->as_recon_buf[i].u4_buf_size = pic_size;
        ps_app_ctxt->as_recon_buf[i].u4_is_free = 1;
    }
    if(ps_app_ctxt->u4_psnr_enable)
    {
        UWORD8 *pu1_buf = (UWORD8 *)ih264a_aligned_malloc(16, pic_size);

        if(NULL == pu1_buf)
        {
            CHAR ac_error[STRLENGTH];
            sprintf(ac_error, "Allocation failed for recon buffer of size %d\n",
                    pic_size);
            codec_exit(ac_error);
        }
        ps_app_ctxt->pu1_psnr_buf = pu1_buf;
        ps_app_ctxt->u4_psnr_buf_size = pic_size;
    }
}

/**
**************************************************************************
* @brief free recon buffers
**************************************************************************
*/
void free_recon(app_ctxt_t *ps_app_ctxt)
{
    WORD32 num_bufs = DEFAULT_NUM_RECON_BUFS;
    WORD32 i;

    for(i = 0; i < num_bufs; i++)
    {
        ih264a_aligned_free(ps_app_ctxt->as_recon_buf[i].pu1_buf);
    }
    if(ps_app_ctxt->u4_psnr_enable)
    {
        ih264a_aligned_free(ps_app_ctxt->pu1_psnr_buf);
    }
}

/**
**************************************************************************
* @brief initialize raw buffer descriptor
* All the pointers and dimensions are initialized here to support change in
* resolution from the application
**************************************************************************
*/
void init_raw_buf_descr(app_ctxt_t *ps_app_ctxt,
                        iv_raw_buf_t *ps_raw_buf,
                        UWORD8 *pu1_buf,
                        IV_COLOR_FORMAT_T e_color_fmt)
{
    WORD32 luma_size = ps_app_ctxt->u4_max_wd * ps_app_ctxt->u4_max_ht;
    WORD32 chroma_size = (luma_size) / 4;

    ps_raw_buf->u4_size = sizeof(iv_raw_buf_t);
    ps_raw_buf->e_color_fmt = e_color_fmt;

    ps_raw_buf->apv_bufs[0] = pu1_buf;
    pu1_buf += luma_size;
    ps_raw_buf->apv_bufs[1] = pu1_buf;
    pu1_buf += chroma_size;
    ps_raw_buf->apv_bufs[2] = NULL;
    if(IV_YUV_420P == e_color_fmt)
    {
        ps_raw_buf->apv_bufs[2] = pu1_buf;
    }

    ps_raw_buf->au4_wd[0] = ps_app_ctxt->u4_wd;
    ps_raw_buf->au4_ht[0] = ps_app_ctxt->u4_ht;
    ps_raw_buf->au4_strd[0] = ps_app_ctxt->u4_wd;
    /* Initialize for 420SP */
    {
        ps_raw_buf->au4_wd[1] = ps_app_ctxt->u4_wd;
        ps_raw_buf->au4_wd[2] = 0;

        ps_raw_buf->au4_ht[1] = ps_app_ctxt->u4_ht / 2;
        ps_raw_buf->au4_ht[2] = 0;

        ps_raw_buf->au4_strd[1] = ps_app_ctxt->u4_wd;
        ps_raw_buf->au4_strd[2] = 0;
    }
    if(IV_YUV_420P == e_color_fmt)
    {
        ps_raw_buf->au4_wd[1] = ps_app_ctxt->u4_wd / 2;
        ps_raw_buf->au4_wd[2] = ps_app_ctxt->u4_wd / 2;

        ps_raw_buf->au4_ht[1] = ps_app_ctxt->u4_ht / 2;
        ps_raw_buf->au4_ht[2] = ps_app_ctxt->u4_ht / 2;

        ps_raw_buf->au4_strd[1] = ps_app_ctxt->u4_wd / 2;
        ps_raw_buf->au4_strd[2] = ps_app_ctxt->u4_wd / 2;
    }
    /* If stride is not initialized, then use width as stride */
    if(0 == ps_raw_buf->au4_strd[0])
    {
        ps_raw_buf->au4_strd[0] = ps_raw_buf->au4_wd[0];
        ps_raw_buf->au4_strd[1] = ps_raw_buf->au4_wd[1];
        ps_raw_buf->au4_strd[2] = ps_raw_buf->au4_wd[2];
    }
}

