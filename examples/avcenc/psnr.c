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
*  psnr.c
*
* @brief
*  Contains functions necessary for computing psnr
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* User include files */
#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "ih264e.h"
#include "app.h"
#include "psnr.h"


/*****************************************************************************/
/*  Function Definition                                                      */
/*****************************************************************************/

/**
**************************************************************************
* @brief Initialize PSNR for the Y, U, V component
**************************************************************************
*/
void init_psnr(app_ctxt_t *ps_app_ctxt)
{
    ps_app_ctxt->adbl_psnr[0] = 0;
    ps_app_ctxt->adbl_psnr[1] = 0;
    ps_app_ctxt->adbl_psnr[2] = 0;
    ps_app_ctxt->u4_psnr_cnt = 0;
}

/**
**************************************************************************
* @brief Computes PSNR for the Y, U, V component
**************************************************************************
*/
void compute_psnr(app_ctxt_t *ps_app_ctxt,
                  iv_raw_buf_t *ps_buf1,
                  iv_raw_buf_t *ps_buf2)
{
    WORD32 i, j;
    WORD32 comp;
    DOUBLE df_psnr[3];
    WORD32 wd, ht, strd1, strd2;
    UWORD8 *pu1_buf1, *pu1_buf2;
    WORD32 incr1, incr2;

    printf("\nPicNum %4d\t ", ps_app_ctxt->u4_psnr_cnt);

    for(comp = 0; comp < 3; comp++)
    {
        df_psnr[comp] = 0;
        pu1_buf1 = (UWORD8 *)ps_buf1->apv_bufs[comp];
        pu1_buf2 = (UWORD8 *)ps_buf2->apv_bufs[comp];
        wd = ps_buf1->au4_wd[comp];
        ht = ps_buf1->au4_ht[comp];
        strd1 = ps_buf1->au4_strd[comp] - ps_buf1->au4_wd[comp];
        strd2 = ps_buf2->au4_strd[comp] - ps_buf2->au4_wd[comp];
        incr1 = 1;
        incr2 = 1;

        if((IV_YUV_420SP_UV == ps_buf1->e_color_fmt)
                        || (IV_YUV_420SP_VU == ps_buf1->e_color_fmt))
        {
            switch(comp)
            {
                case 0:
                    pu1_buf1 = ps_buf1->apv_bufs[0];
                    break;
                case 1:
                    if(IV_YUV_420SP_UV == ps_buf1->e_color_fmt)
                        pu1_buf1 = (UWORD8 *)ps_buf1->apv_bufs[1];
                    else
                        pu1_buf1 = (UWORD8 *)ps_buf1->apv_bufs[1] + 1;
                    incr1 = 2;
                    wd = ps_buf1->au4_wd[0] >> 1;
                    ht = ps_buf1->au4_ht[0] >> 1;
                    break;
                case 2:
                    if(IV_YUV_420SP_UV == ps_buf1->e_color_fmt)
                        pu1_buf1 = (UWORD8 *)ps_buf1->apv_bufs[1] + 1;
                    else
                        pu1_buf1 = ps_buf1->apv_bufs[1];
                    incr1 = 2;
                    wd = ps_buf1->au4_wd[0] >> 1;
                    ht = ps_buf1->au4_ht[0] >> 1;
                    strd1 = ps_buf1->au4_strd[1] - ps_buf1->au4_wd[1];
                    break;
            }
        }
        if((IV_YUV_420SP_UV == ps_buf2->e_color_fmt)
                        || (IV_YUV_420SP_VU == ps_buf2->e_color_fmt))
        {
            switch(comp)
            {
                case 0:
                    pu1_buf2 = ps_buf2->apv_bufs[0];
                    break;
                case 1:
                    if(IV_YUV_420SP_UV == ps_buf2->e_color_fmt)
                        pu1_buf2 = ps_buf2->apv_bufs[1];
                    else
                        pu1_buf2 = (UWORD8 *)ps_buf2->apv_bufs[1] + 1;
                    incr2 = 2;
                    wd = ps_buf1->au4_wd[0] >> 1;
                    ht = ps_buf1->au4_ht[0] >> 1;

                    break;
                case 2:
                    if(IV_YUV_420SP_UV == ps_buf2->e_color_fmt)
                        pu1_buf2 = (UWORD8 *)ps_buf2->apv_bufs[1] + 1;
                    else
                        pu1_buf2 = ps_buf2->apv_bufs[1];
                    incr2 = 2;
                    wd = ps_buf1->au4_wd[0] >> 1;
                    ht = ps_buf1->au4_ht[0] >> 1;
                    strd2 = ps_buf2->au4_strd[1] - ps_buf2->au4_wd[1];

                    break;
            }
        }

        for(i = 0; i < ht; i++)
        {
            for(j = 0; j < wd; j++)
            {
                WORD32 diff;
                diff = (*pu1_buf1 - *pu1_buf2);
                pu1_buf1 += incr1;
                pu1_buf2 += incr2;
                df_psnr[comp] += diff * diff;
            }
            pu1_buf1 += strd1;
            pu1_buf2 += strd2;
        }
        df_psnr[comp] /= (wd * ht);
        if(df_psnr[comp])
            df_psnr[comp] = 20 * log10(255 / sqrt(df_psnr[comp]));
        else
            df_psnr[comp] = 100;

        ps_app_ctxt->adbl_psnr[comp] += df_psnr[comp];
        switch(comp)
        {
            case 0:
                printf("Y :");
                break;
            case 1:
                printf("U :");
                break;
            case 2:
                printf("V :");
                break;
            default:
                break;
        }
        printf("%2.2f\t", df_psnr[comp]);
    }

    ps_app_ctxt->u4_psnr_cnt++;
}

/**
**************************************************************************
* @brief Prints PSNR for the Y, U, V component
**************************************************************************
*/
void print_average_psnr(app_ctxt_t *ps_app_ctxt)
{
    printf("\n");

    printf("Avg PSNR Y                      : %-2.2f\n",
           (ps_app_ctxt->adbl_psnr[0] / ps_app_ctxt->u4_psnr_cnt));
    printf("Avg PSNR U                      : %-2.2f\n",
           (ps_app_ctxt->adbl_psnr[1] / ps_app_ctxt->u4_psnr_cnt));
    printf("Avg PSNR V                      : %-2.2f\n",
           (ps_app_ctxt->adbl_psnr[2] / ps_app_ctxt->u4_psnr_cnt));
}

