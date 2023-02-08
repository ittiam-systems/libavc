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
 *  isvcd_pred_residual_recon.c
 *
 * @brief
 *  Contains definition of functions for svc inverse quantization inverse
 *    transformation and resd comp
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_pred_residual_recon_chroma_8x8()
 *  - isvcd_residual_chroma_cb_cr_8x8()
 *  - isvcd_pred_residual_recon_chroma_4x4()
 *  - isvcd_pred_residual_recon_16x16()
 *  - isvcd_pred_residual_recon_4x4()
 *  - isvcd_pred_residual_recon_8x8()
 *  - isvcd_residual_luma_4x4()
 *  - isvcd_residual_luma_8x8()
 *  - isvcd_residual_luma_16x16()
 *
 * @remarks
 *   None
 *
 *******************************************************************************
 */

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* User include files */
#include "ih264_typedefs.h"
#include "ih264_defs.h"
#include "ih264_trans_macros.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_trans_data.h"
#include "ih264_size_defs.h"
#include "ih264_structs.h"
#include "isvcd_pred_residual_recon.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_chroma_8x8                      */
/*                                                                           */
/*  Description   : this function computes the recon from                    */
/*                  the residual and pred buffer                             */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

void isvcd_pred_residual_recon_chroma_8x8(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                          WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    UWORD8 *pu1_pred_ptr = pu1_pred;
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    UWORD8 *pu1_out_ptr = pu1_out;
    WORD16 i, j;
    WORD16 i_macro;

    for(i = 0; i < 8; i++)
    {
        pu1_pred_ptr = pu1_pred;
        pi2_rsd_ptr = pi2_rsd;
        pu1_out = pu1_out_ptr;

        for(j = 0; j < 8; j++)
        {
            i_macro = *pu1_pred_ptr + *pi2_rsd_ptr;
            *pu1_out = CLIP_U8(i_macro);
            pu1_pred_ptr += pred_strd;
            pi2_rsd_ptr += rsd_strd;
            pu1_out += out_strd;
        }

        pu1_out_ptr += 2;  // Interleaved store for output
        pu1_pred += 2;     // Interleaved load for pred buffer
        pi2_rsd += 2;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_chroma_cb_cr_8x8                           */
/*                                                                           */
/*  Description   : this function computes the nnz from the resd             */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_residual_chroma_cb_cr_8x8(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    WORD16 *pi2_rsd_ptr_Cb = pi2_rsd;
    WORD16 *pi2_rsd_ptr_Cr = pi2_rsd + 1;
    WORD16 i, j;
    WORD32 i4_nnz = 0, ai4_nnz_Cb[2][2] = {0}, ai4_nnz_Cr[2][2] = {0};

    for(i = 0; i < 8; i++)
    {
        pi2_rsd_ptr_Cb = pi2_rsd;
        pi2_rsd_ptr_Cr = pi2_rsd + 1;

        for(j = 0; j < 8; j++)
        {
            ai4_nnz_Cb[j >> 2][i >> 2] |= !!(*pi2_rsd_ptr_Cb);
            ai4_nnz_Cr[j >> 2][i >> 2] |= !!(*pi2_rsd_ptr_Cr);
            pi2_rsd_ptr_Cb += rsd_strd;
            pi2_rsd_ptr_Cr += rsd_strd;
        }
        pi2_rsd += 2;
    }
    i4_nnz = ai4_nnz_Cr[0][0] | (ai4_nnz_Cr[1][0] << 2);
    i4_nnz |= (ai4_nnz_Cr[0][1] << 1) | (ai4_nnz_Cr[1][1] << 3);
    i4_nnz <<= 4;
    i4_nnz |= ai4_nnz_Cb[0][0] | (ai4_nnz_Cb[1][0] << 2);
    i4_nnz |= (ai4_nnz_Cb[0][1] << 1) | (ai4_nnz_Cb[1][1] << 3);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_chroma_4x4                      */
/*                                                                           */
/*  Description   : this function computes the recon from                    */
/*                  the residual and pred buffer                             */
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
void isvcd_pred_residual_recon_chroma_4x4(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                          WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    UWORD8 *pu1_pred_ptr = pu1_pred;
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    UWORD8 *pu1_out_ptr = pu1_out;
    WORD16 i, j;
    WORD16 i_macro;

    for(i = 0; i < 4; i++)
    {
        pu1_pred_ptr = pu1_pred;
        pi2_rsd_ptr = pi2_rsd;
        pu1_out = pu1_out_ptr;

        for(j = 0; j < 4; j++)
        {
            i_macro = *pu1_pred_ptr + *pi2_rsd_ptr;
            *pu1_out = CLIP_U8(i_macro);
            pu1_pred_ptr += pred_strd;
            pi2_rsd_ptr += rsd_strd;
            pu1_out += out_strd;
        }

        pu1_out_ptr += 2;  // Interleaved store for output
        pu1_pred += 2;     // Interleaved load for pred buffer
        pi2_rsd += 2;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_16x16                           */
/*                                                                           */
/*  Description   : this function computes the recon from                    */
/*                  the residual and pred buffer                             */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_pred_residual_recon_16x16(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                       WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    WORD32 i4_nnz = 0, i4_nnz_blk[4][4] = {0};
    UWORD8 *pu1_pred_ptr = pu1_pred;
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    UWORD8 *pu1_out_ptr = pu1_out;
    WORD16 i, j;
    WORD16 i_macro;

    for(i = 0; i < 16; i++)
    {
        pu1_pred_ptr = pu1_pred;
        pi2_rsd_ptr = pi2_rsd;
        pu1_out = pu1_out_ptr;

        for(j = 0; j < 16; j++)
        {
            i_macro = *pi2_rsd_ptr;
            i4_nnz_blk[j >> 2][i >> 2] |= !!i_macro;
            i_macro += *pu1_pred_ptr;
            *pu1_out = CLIP_U8(i_macro);
            pu1_pred_ptr += pred_strd;
            pi2_rsd_ptr += rsd_strd;
            pu1_out += out_strd;
        }

        pu1_out_ptr++;
        pi2_rsd++;
        pu1_pred++;
    }

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            i4_nnz |= (i4_nnz_blk[j][i]) << (i + (j << 2));
        }
    }

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_4x4                             */
/*                                                                           */
/*  Description   : this function computes the recon from                    */
/*                  the residual and pred buffer                             */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_pred_residual_recon_4x4(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                     WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    WORD32 i4_nnz_blk[4][4] = {0};
    UWORD8 *pu1_pred_ptr = pu1_pred;
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    UWORD8 *pu1_out_ptr = pu1_out;
    WORD16 i, j;
    WORD16 i_macro;

    for(i = 0; i < 4; i++)
    {
        pu1_pred_ptr = pu1_pred;
        pi2_rsd_ptr = pi2_rsd;
        pu1_out = pu1_out_ptr;

        for(j = 0; j < 4; j++)
        {
            i_macro = *pi2_rsd_ptr;
            i4_nnz_blk[j >> 2][i >> 2] |= !!i_macro;
            i_macro += *pu1_pred_ptr;
            *pu1_out = CLIP_U8(i_macro);
            pu1_pred_ptr += pred_strd;
            pi2_rsd_ptr += rsd_strd;
            pu1_out += out_strd;
        }

        pu1_out_ptr++;
        pi2_rsd++;
        pu1_pred++;
    }

    return i4_nnz_blk[0][0];
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_8x8                             */
/*                                                                           */
/*  Description   : this function computes the recon from                    */
/*                  the residual and pred buffer                             */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_pred_residual_recon_8x8(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                     WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    WORD32 i4_nnz = 0, i4_nnz_blk[4][4] = {0};
    UWORD8 *pu1_pred_ptr = pu1_pred;
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    UWORD8 *pu1_out_ptr = pu1_out;
    WORD16 i, j;
    WORD16 i_macro;

    for(i = 0; i < 8; i++)
    {
        pu1_pred_ptr = pu1_pred;
        pi2_rsd_ptr = pi2_rsd;
        pu1_out = pu1_out_ptr;

        for(j = 0; j < 8; j++)
        {
            i_macro = *pi2_rsd_ptr;
            i4_nnz_blk[j >> 2][i >> 2] |= !!i_macro;
            i_macro += *pu1_pred_ptr;
            *pu1_out = CLIP_U8(i_macro);
            pu1_pred_ptr += pred_strd;
            pi2_rsd_ptr += rsd_strd;
            pu1_out += out_strd;
        }

        pu1_out_ptr++;
        pi2_rsd++;
        pu1_pred++;
    }

    i4_nnz = i4_nnz_blk[0][0] | (i4_nnz_blk[1][0] << 4);
    i4_nnz |= (i4_nnz_blk[0][1] << 1) | (i4_nnz_blk[1][1] << 5);

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_4x4                                   */
/*                                                                           */
/*  Description   : this function computes the nnz from resd                 */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_residual_luma_4x4(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    WORD32 i4_nnz_blk[4][4] = {0};
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    WORD16 i, j;
    WORD16 i_macro;

    for(i = 0; i < 4; i++)
    {
        pi2_rsd_ptr = pi2_rsd;

        for(j = 0; j < 4; j++)
        {
            i_macro = *pi2_rsd_ptr;
            i4_nnz_blk[j >> 2][i >> 2] |= !!i_macro;
            pi2_rsd_ptr += rsd_strd;
        }
        pi2_rsd++;
    }
    return i4_nnz_blk[0][0];
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_8x8                                   */
/*                                                                           */
/*  Description   : this function computes the nnz from resd                 */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_residual_luma_8x8(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    WORD32 i4_nnz = 0, i4_nnz_blk[4][4] = {0};
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    WORD16 i, j;
    WORD16 i_macro;

    for(i = 0; i < 8; i++)
    {
        pi2_rsd_ptr = pi2_rsd;

        for(j = 0; j < 8; j++)
        {
            i_macro = *pi2_rsd_ptr;
            i4_nnz_blk[j >> 2][i >> 2] |= !!i_macro;
            pi2_rsd_ptr += rsd_strd;
        }
        pi2_rsd++;
    }

    i4_nnz = i4_nnz_blk[0][0] | (i4_nnz_blk[1][0] << 4);
    i4_nnz |= (i4_nnz_blk[0][1] << 1) | (i4_nnz_blk[1][1] << 5);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_16x16                                 */
/*                                                                           */
/*  Description   : this function computes the nnz from resd                 */
/*                                                                           */
/*  Inputs        :                                                          */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : nnz                                                      */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 11 2021   Kishore               creation                       */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_residual_luma_16x16(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    WORD32 i4_nnz = 0, i4_nnz_blk[4][4] = {0};
    WORD16 *pi2_rsd_ptr = pi2_rsd;
    WORD16 i, j;
    WORD16 i_macro;

    for(i = 0; i < 16; i++)
    {
        pi2_rsd_ptr = pi2_rsd;

        for(j = 0; j < 16; j++)
        {
            i_macro = *pi2_rsd_ptr;
            i4_nnz_blk[j >> 2][i >> 2] |= !!i_macro;
            pi2_rsd_ptr += rsd_strd;
        }
        pi2_rsd++;
    }

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            i4_nnz |= (i4_nnz_blk[j][i]) << (i + (j << 2));
        }
    }
    return i4_nnz;
}
