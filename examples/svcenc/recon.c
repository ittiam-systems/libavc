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
#include "isvce.h"
#include "app.h"

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*  Macros                                                                   */
/*****************************************************************************/

/*****************************************************************************/
/*  Function Declarations                                                    */
/*****************************************************************************/

IV_STATUS_T write_recon(FILE *fp, iv_raw_buf_t *ps_raw_buf)
{
    WORD32 bytes;
    WORD32 wd, ht;
    UWORD8 *pu1_buf;
    WORD32 i;
    WORD32 comp;
    WORD32 num_comp;

    num_comp = 2;
    if(IV_YUV_420P == ps_raw_buf->e_color_fmt) num_comp = 3;

    for(comp = 0; comp < num_comp; comp++)
    {
        wd = ps_raw_buf->au4_wd[comp];
        ht = ps_raw_buf->au4_ht[comp];
        pu1_buf = ps_raw_buf->apv_bufs[comp];
        for(i = 0; i < ht; i++)
        {
            bytes = (WORD32) fwrite(pu1_buf, sizeof(UWORD8), wd, fp);
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

void allocate_recon(app_ctxt_t *ps_app_ctxt)
{
    WORD32 num_bufs;
    WORD32 pic_size;
    WORD32 luma_size;
    WORD32 chroma_size;
    WORD32 i;
    UWORD8 *pu1_buf;

    num_bufs = DEFAULT_NUM_RECON_BUFS;

    assert(ps_app_ctxt->s_get_buf_info_op.u4_rec_comp_cnt == 3);

    /* Size of buffer for YUV420/420SP */
    luma_size = ps_app_ctxt->s_get_buf_info_op.au4_min_rec_buf_size[0];
    chroma_size = ps_app_ctxt->s_get_buf_info_op.au4_min_rec_buf_size[1] +
                  ps_app_ctxt->s_get_buf_info_op.au4_min_rec_buf_size[2];
    pic_size = luma_size + chroma_size;

    for(i = 0; i < num_bufs; i++)
    {
        pu1_buf = (UWORD8 *) isvca_aligned_malloc(16, pic_size);

        if(NULL == pu1_buf)
        {
            CHAR ac_error[STRLENGTH];
            snprintf(ac_error, sizeof(ac_error) - 1,
                     "Allocation failed for recon buffer of size %d\n", pic_size);
            codec_exit(ac_error);
        }

        ps_app_ctxt->as_recon_buf[i].pu1_buf = pu1_buf;
        ps_app_ctxt->as_recon_buf[i].u4_buf_size = pic_size;
        ps_app_ctxt->as_recon_buf[i].u4_is_free = 1;
    }

    if(ps_app_ctxt->u4_psnr_enable)
    {
        pu1_buf = (UWORD8 *) isvca_aligned_malloc(16, pic_size);

        if(NULL == pu1_buf)
        {
            CHAR ac_error[STRLENGTH];
            snprintf(ac_error, sizeof(ac_error) - 1,
                     "Allocation failed for recon buffer of size %d\n", pic_size);
            codec_exit(ac_error);
        }

        ps_app_ctxt->pu1_psnr_buf = pu1_buf;
        ps_app_ctxt->u4_psnr_buf_size = pic_size;
    }
}

void free_recon(app_ctxt_t *ps_app_ctxt)
{
    WORD32 num_bufs;
    WORD32 i;

    num_bufs = DEFAULT_NUM_RECON_BUFS;

    for(i = 0; i < num_bufs; i++)
    {
        isvca_aligned_free(ps_app_ctxt->as_recon_buf[i].pu1_buf);
    }

    if(ps_app_ctxt->u4_psnr_enable)
    {
        isvca_aligned_free(ps_app_ctxt->pu1_psnr_buf);
    }
}

void init_raw_buf_descr(app_ctxt_t *ps_app_ctxt, iv_raw_buf_t *ps_raw_buf, UWORD8 *pu1_buf,
                        IV_COLOR_FORMAT_T e_color_fmt)
{
    WORD32 luma_size;
    WORD32 au4_chroma_sizes[2];

    assert(IV_YUV_420P == e_color_fmt);

    luma_size = ps_app_ctxt->s_get_buf_info_op.au4_min_rec_buf_size[0];
    au4_chroma_sizes[0] = ps_app_ctxt->s_get_buf_info_op.au4_min_rec_buf_size[1];
    au4_chroma_sizes[1] = ps_app_ctxt->s_get_buf_info_op.au4_min_rec_buf_size[2];

    ps_raw_buf->apv_bufs[0] = pu1_buf;
    pu1_buf += luma_size;

    ps_raw_buf->apv_bufs[1] = pu1_buf;
    pu1_buf += au4_chroma_sizes[0];

    ps_raw_buf->apv_bufs[2] = NULL;
    if(IV_YUV_420P == e_color_fmt)
    {
        ps_raw_buf->apv_bufs[2] = pu1_buf;
    }

    ps_raw_buf->e_color_fmt = e_color_fmt;
    ps_raw_buf->au4_wd[0] = ps_app_ctxt->u4_enc_wd;
    ps_raw_buf->au4_ht[0] = ps_app_ctxt->u4_enc_ht;
    ps_raw_buf->au4_strd[0] = ps_app_ctxt->u4_enc_wd;

    /* Initialize for 420SP */
    {
        ps_raw_buf->au4_wd[1] = ps_app_ctxt->u4_enc_wd;
        ps_raw_buf->au4_wd[2] = 0;

        ps_raw_buf->au4_ht[1] = ps_app_ctxt->u4_enc_ht / 2;
        ps_raw_buf->au4_ht[2] = 0;

        ps_raw_buf->au4_strd[1] = ps_app_ctxt->u4_enc_wd;
        ps_raw_buf->au4_strd[2] = 0;
    }

    if(IV_YUV_420P == e_color_fmt)
    {
        ps_raw_buf->au4_wd[1] = ps_app_ctxt->u4_enc_wd / 2;
        ps_raw_buf->au4_wd[2] = ps_app_ctxt->u4_enc_wd / 2;

        ps_raw_buf->au4_ht[1] = ps_app_ctxt->u4_enc_ht / 2;
        ps_raw_buf->au4_ht[2] = ps_app_ctxt->u4_enc_ht / 2;

        ps_raw_buf->au4_strd[1] = ps_app_ctxt->u4_enc_wd / 2;
        ps_raw_buf->au4_strd[2] = ps_app_ctxt->u4_enc_wd / 2;
    }
    /* If stride is not initialized, then use width as stride */
    if(0 == ps_raw_buf->au4_strd[0])
    {
        ps_raw_buf->au4_strd[0] = ps_raw_buf->au4_wd[0];
        ps_raw_buf->au4_strd[1] = ps_raw_buf->au4_wd[1];
        ps_raw_buf->au4_strd[2] = ps_raw_buf->au4_wd[2];
    }

    ps_raw_buf->u4_size = sizeof(iv_raw_buf_t);
}
