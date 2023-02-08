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
 *  isvcd_pred_residual_recon_sse42.c
 *
 * @brief
 *  Contains function definitions for pred_residual and recon transform
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_pred_residual_recon_4x4_sse42()
 *  - isvcd_pred_residual_recon_8x8_sse42()
 *  - isvcd_pred_residual_recon_16x16_sse42()
 *  - isvcd_pred_residual_recon_chroma_4x4_sse42()
 *  - isvcd_pred_residual_recon_chroma_8x8_sse42()
 *  - isvcd_residual_luma_4x4_sse42()
 *  - isvcd_residual_luma_8x8_sse42()
 *  - isvcd_residual_luma_16x16_sse42()
 *  - isvcd_residual_chroma_cb_cr_8x8_sse42()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
/* User include files */
#include <immintrin.h>
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
/*  Function Name : isvcd_pred_residual_recon_4x4_sse42                       */
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

WORD32 isvcd_pred_residual_recon_4x4_sse42(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                           WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    __m128i pred_16x8b_0, pred_8x16b_0, rsd_8x16b_0, out_8x16b_0, out_16x8b_0;
    __m128i pred_16x8b_1, pred_8x16b_1, rsd_8x16b_1, out_8x16b_1, out_16x8b_1;
    __m128i pred_16x8b_2, pred_8x16b_2, rsd_8x16b_2, out_8x16b_2, out_16x8b_2;
    __m128i pred_16x8b_3, pred_8x16b_3, rsd_8x16b_3, out_8x16b_3, out_16x8b_3;
    __m128i rsd_8x16b_01, rsd_8x16b_23;

    __m128i zero_8x16b = _mm_setzero_si128();
    WORD32 i4_nnz, row_01, row_23;

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + (pred_strd << 1)));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + (pred_strd << 1) + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + (rsd_strd << 1)));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + (rsd_strd << 1) + rsd_strd));

    rsd_8x16b_01 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);

    row_01 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01, zero_8x16b));  // return 1 if all zeros, else 0
    row_23 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23, zero_8x16b));

    out_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    out_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    out_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    out_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);

    out_16x8b_0 = _mm_packus_epi16(out_8x16b_0, zero_8x16b);
    out_16x8b_1 = _mm_packus_epi16(out_8x16b_1, zero_8x16b);
    out_16x8b_2 = _mm_packus_epi16(out_8x16b_2, zero_8x16b);
    out_16x8b_3 = _mm_packus_epi16(out_8x16b_3, zero_8x16b);

    *((WORD32 *) (pu1_out)) = _mm_cvtsi128_si32(out_16x8b_0);
    *((WORD32 *) (pu1_out + out_strd)) = _mm_cvtsi128_si32(out_16x8b_1);
    *((WORD32 *) (pu1_out + (out_strd << 1))) = _mm_cvtsi128_si32(out_16x8b_2);
    *((WORD32 *) (pu1_out + (out_strd * 3))) = _mm_cvtsi128_si32(out_16x8b_3);
    i4_nnz = !(row_01 && row_23);

    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_8x8_sse42                       */
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

WORD32 isvcd_pred_residual_recon_8x8_sse42(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                           WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    __m128i pred_16x8b_0, pred_8x16b_0, rsd_8x16b_0, out_8x16b_0, out_16x8b_0;
    __m128i pred_16x8b_1, pred_8x16b_1, rsd_8x16b_1, out_8x16b_1, out_16x8b_1;
    __m128i pred_16x8b_2, pred_8x16b_2, rsd_8x16b_2, out_8x16b_2, out_16x8b_2;
    __m128i pred_16x8b_3, pred_8x16b_3, rsd_8x16b_3, out_8x16b_3, out_16x8b_3;
    __m128i pred_16x8b_4, pred_8x16b_4, rsd_8x16b_4, out_8x16b_4, out_16x8b_4;
    __m128i pred_16x8b_5, pred_8x16b_5, rsd_8x16b_5, out_8x16b_5, out_16x8b_5;
    __m128i pred_16x8b_6, pred_8x16b_6, rsd_8x16b_6, out_8x16b_6, out_16x8b_6;
    __m128i pred_16x8b_7, pred_8x16b_7, rsd_8x16b_7, out_8x16b_7, out_16x8b_7;
    __m128i rsd_8x16b_01_b0, rsd_8x16b_23_b0, rsd_8x16b_45_b2, rsd_8x16b_67_b2;
    __m128i rsd_8x16b_01_b1, rsd_8x16b_23_b1, rsd_8x16b_45_b3, rsd_8x16b_67_b3;

    WORD32 row_01_b0, row_23_b0, row_45_b2, row_67_b2;
    WORD32 row_01_b1, row_23_b1, row_45_b3, row_67_b3;
    WORD32 i4_nnz, i4_nnz_b0, i4_nnz_b1, i4_nnz_b2, i4_nnz_b3;

    __m128i zero_8x16b = _mm_setzero_si128();

    WORD32 pred_strd2 = (pred_strd << 1);
    WORD32 pred_strd4 = (pred_strd << 2);
    WORD32 rsd_strd2 = (rsd_strd << 1);
    WORD32 rsd_strd4 = (rsd_strd << 2);
    WORD32 out_strd2 = (out_strd << 1);
    WORD32 out_strd4 = (out_strd << 2);

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2 + pred_strd));
    pred_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4));
    pred_16x8b_5 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd));
    pred_16x8b_6 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2));
    pred_16x8b_7 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2 + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);
    pred_8x16b_4 = _mm_cvtepu8_epi16(pred_16x8b_4);
    pred_8x16b_5 = _mm_cvtepu8_epi16(pred_16x8b_5);
    pred_8x16b_6 = _mm_cvtepu8_epi16(pred_16x8b_6);
    pred_8x16b_7 = _mm_cvtepu8_epi16(pred_16x8b_7);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    out_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    out_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    out_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    out_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);
    out_8x16b_4 = _mm_add_epi16(pred_8x16b_4, rsd_8x16b_4);
    out_8x16b_5 = _mm_add_epi16(pred_8x16b_5, rsd_8x16b_5);
    out_8x16b_6 = _mm_add_epi16(pred_8x16b_6, rsd_8x16b_6);
    out_8x16b_7 = _mm_add_epi16(pred_8x16b_7, rsd_8x16b_7);

    out_16x8b_0 = _mm_packus_epi16(out_8x16b_0, zero_8x16b);
    out_16x8b_1 = _mm_packus_epi16(out_8x16b_1, zero_8x16b);
    out_16x8b_2 = _mm_packus_epi16(out_8x16b_2, zero_8x16b);
    out_16x8b_3 = _mm_packus_epi16(out_8x16b_3, zero_8x16b);
    out_16x8b_4 = _mm_packus_epi16(out_8x16b_4, zero_8x16b);
    out_16x8b_5 = _mm_packus_epi16(out_8x16b_5, zero_8x16b);
    out_16x8b_6 = _mm_packus_epi16(out_8x16b_6, zero_8x16b);
    out_16x8b_7 = _mm_packus_epi16(out_8x16b_7, zero_8x16b);

    _mm_storel_epi64((__m128i *) (pu1_out), out_16x8b_0);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd), out_16x8b_1);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2), out_16x8b_2);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2 + out_strd), out_16x8b_3);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4), out_16x8b_4);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd), out_16x8b_5);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2), out_16x8b_6);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2 + out_strd), out_16x8b_7);

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0));
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 1;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 4;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 5;

    i4_nnz = (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_16x16_sse42                     */
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

WORD32 isvcd_pred_residual_recon_16x16_sse42(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                             WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    __m128i pred_16x8b_0, pred_8x16b_0, rsd_8x16b_0, out_8x16b_0, out_16x8b_0;
    __m128i pred_16x8b_1, pred_8x16b_1, rsd_8x16b_1, out_8x16b_1, out_16x8b_1;
    __m128i pred_16x8b_2, pred_8x16b_2, rsd_8x16b_2, out_8x16b_2, out_16x8b_2;
    __m128i pred_16x8b_3, pred_8x16b_3, rsd_8x16b_3, out_8x16b_3, out_16x8b_3;
    __m128i pred_16x8b_4, pred_8x16b_4, rsd_8x16b_4, out_8x16b_4, out_16x8b_4;
    __m128i pred_16x8b_5, pred_8x16b_5, rsd_8x16b_5, out_8x16b_5, out_16x8b_5;
    __m128i pred_16x8b_6, pred_8x16b_6, rsd_8x16b_6, out_8x16b_6, out_16x8b_6;
    __m128i pred_16x8b_7, pred_8x16b_7, rsd_8x16b_7, out_8x16b_7, out_16x8b_7;
    __m128i rsd_8x16b_01_b0, rsd_8x16b_23_b0, rsd_8x16b_45_b2, rsd_8x16b_67_b2;
    __m128i rsd_8x16b_01_b1, rsd_8x16b_23_b1, rsd_8x16b_45_b3, rsd_8x16b_67_b3;

    WORD32 row_01_b0, row_23_b0, row_45_b2, row_67_b2;
    WORD32 row_01_b1, row_23_b1, row_45_b3, row_67_b3;
    WORD32 i4_nnz, i4_nnz_b0, i4_nnz_b1, i4_nnz_b2, i4_nnz_b3;

    __m128i zero_8x16b = _mm_setzero_si128();

    WORD32 pred_strd2 = (pred_strd << 1);
    WORD32 pred_strd4 = (pred_strd << 2);
    WORD32 rsd_strd2 = (rsd_strd << 1);
    WORD32 rsd_strd4 = (rsd_strd << 2);
    WORD32 out_strd2 = (out_strd << 1);
    WORD32 out_strd4 = (out_strd << 2);

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2 + pred_strd));
    pred_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4));
    pred_16x8b_5 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd));
    pred_16x8b_6 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2));
    pred_16x8b_7 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2 + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);
    pred_8x16b_4 = _mm_cvtepu8_epi16(pred_16x8b_4);
    pred_8x16b_5 = _mm_cvtepu8_epi16(pred_16x8b_5);
    pred_8x16b_6 = _mm_cvtepu8_epi16(pred_16x8b_6);
    pred_8x16b_7 = _mm_cvtepu8_epi16(pred_16x8b_7);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    out_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    out_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    out_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    out_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);
    out_8x16b_4 = _mm_add_epi16(pred_8x16b_4, rsd_8x16b_4);
    out_8x16b_5 = _mm_add_epi16(pred_8x16b_5, rsd_8x16b_5);
    out_8x16b_6 = _mm_add_epi16(pred_8x16b_6, rsd_8x16b_6);
    out_8x16b_7 = _mm_add_epi16(pred_8x16b_7, rsd_8x16b_7);

    out_16x8b_0 = _mm_packus_epi16(out_8x16b_0, zero_8x16b);
    out_16x8b_1 = _mm_packus_epi16(out_8x16b_1, zero_8x16b);
    out_16x8b_2 = _mm_packus_epi16(out_8x16b_2, zero_8x16b);
    out_16x8b_3 = _mm_packus_epi16(out_8x16b_3, zero_8x16b);
    out_16x8b_4 = _mm_packus_epi16(out_8x16b_4, zero_8x16b);
    out_16x8b_5 = _mm_packus_epi16(out_8x16b_5, zero_8x16b);
    out_16x8b_6 = _mm_packus_epi16(out_8x16b_6, zero_8x16b);
    out_16x8b_7 = _mm_packus_epi16(out_8x16b_7, zero_8x16b);

    _mm_storel_epi64((__m128i *) (pu1_out), out_16x8b_0);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd), out_16x8b_1);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2), out_16x8b_2);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2 + out_strd), out_16x8b_3);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4), out_16x8b_4);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd), out_16x8b_5);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2), out_16x8b_6);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2 + out_strd), out_16x8b_7);

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0));
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 1;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 4;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 5;

    i4_nnz = (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);

    pu1_pred += 8;
    pi2_rsd += 8;
    pu1_out += 8;

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2 + pred_strd));
    pred_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4));
    pred_16x8b_5 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd));
    pred_16x8b_6 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2));
    pred_16x8b_7 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2 + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);
    pred_8x16b_4 = _mm_cvtepu8_epi16(pred_16x8b_4);
    pred_8x16b_5 = _mm_cvtepu8_epi16(pred_16x8b_5);
    pred_8x16b_6 = _mm_cvtepu8_epi16(pred_16x8b_6);
    pred_8x16b_7 = _mm_cvtepu8_epi16(pred_16x8b_7);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    out_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    out_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    out_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    out_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);
    out_8x16b_4 = _mm_add_epi16(pred_8x16b_4, rsd_8x16b_4);
    out_8x16b_5 = _mm_add_epi16(pred_8x16b_5, rsd_8x16b_5);
    out_8x16b_6 = _mm_add_epi16(pred_8x16b_6, rsd_8x16b_6);
    out_8x16b_7 = _mm_add_epi16(pred_8x16b_7, rsd_8x16b_7);

    out_16x8b_0 = _mm_packus_epi16(out_8x16b_0, zero_8x16b);
    out_16x8b_1 = _mm_packus_epi16(out_8x16b_1, zero_8x16b);
    out_16x8b_2 = _mm_packus_epi16(out_8x16b_2, zero_8x16b);
    out_16x8b_3 = _mm_packus_epi16(out_8x16b_3, zero_8x16b);
    out_16x8b_4 = _mm_packus_epi16(out_8x16b_4, zero_8x16b);
    out_16x8b_5 = _mm_packus_epi16(out_8x16b_5, zero_8x16b);
    out_16x8b_6 = _mm_packus_epi16(out_8x16b_6, zero_8x16b);
    out_16x8b_7 = _mm_packus_epi16(out_8x16b_7, zero_8x16b);

    _mm_storel_epi64((__m128i *) (pu1_out), out_16x8b_0);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd), out_16x8b_1);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2), out_16x8b_2);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2 + out_strd), out_16x8b_3);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4), out_16x8b_4);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd), out_16x8b_5);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2), out_16x8b_6);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2 + out_strd), out_16x8b_7);

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0)) << 2;
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 3;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 6;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 7;

    i4_nnz |= (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);

    pu1_pred -= 8;
    pi2_rsd -= 8;
    pu1_out -= 8;

    pu1_pred += (pred_strd << 3);
    pi2_rsd += (rsd_strd << 3);
    pu1_out += (out_strd << 3);

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2 + pred_strd));
    pred_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4));
    pred_16x8b_5 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd));
    pred_16x8b_6 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2));
    pred_16x8b_7 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2 + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);
    pred_8x16b_4 = _mm_cvtepu8_epi16(pred_16x8b_4);
    pred_8x16b_5 = _mm_cvtepu8_epi16(pred_16x8b_5);
    pred_8x16b_6 = _mm_cvtepu8_epi16(pred_16x8b_6);
    pred_8x16b_7 = _mm_cvtepu8_epi16(pred_16x8b_7);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    out_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    out_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    out_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    out_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);
    out_8x16b_4 = _mm_add_epi16(pred_8x16b_4, rsd_8x16b_4);
    out_8x16b_5 = _mm_add_epi16(pred_8x16b_5, rsd_8x16b_5);
    out_8x16b_6 = _mm_add_epi16(pred_8x16b_6, rsd_8x16b_6);
    out_8x16b_7 = _mm_add_epi16(pred_8x16b_7, rsd_8x16b_7);

    out_16x8b_0 = _mm_packus_epi16(out_8x16b_0, zero_8x16b);
    out_16x8b_1 = _mm_packus_epi16(out_8x16b_1, zero_8x16b);
    out_16x8b_2 = _mm_packus_epi16(out_8x16b_2, zero_8x16b);
    out_16x8b_3 = _mm_packus_epi16(out_8x16b_3, zero_8x16b);
    out_16x8b_4 = _mm_packus_epi16(out_8x16b_4, zero_8x16b);
    out_16x8b_5 = _mm_packus_epi16(out_8x16b_5, zero_8x16b);
    out_16x8b_6 = _mm_packus_epi16(out_8x16b_6, zero_8x16b);
    out_16x8b_7 = _mm_packus_epi16(out_8x16b_7, zero_8x16b);

    _mm_storel_epi64((__m128i *) (pu1_out), out_16x8b_0);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd), out_16x8b_1);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2), out_16x8b_2);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2 + out_strd), out_16x8b_3);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4), out_16x8b_4);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd), out_16x8b_5);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2), out_16x8b_6);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2 + out_strd), out_16x8b_7);

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0)) << 8;
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 9;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 12;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 13;

    i4_nnz |= (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);

    pu1_pred += 8;
    pi2_rsd += 8;
    pu1_out += 8;

    pred_16x8b_0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_16x8b_1 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd));
    pred_16x8b_2 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2));
    pred_16x8b_3 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd2 + pred_strd));
    pred_16x8b_4 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4));
    pred_16x8b_5 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd));
    pred_16x8b_6 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2));
    pred_16x8b_7 = _mm_loadu_si128((__m128i *) (pu1_pred + pred_strd4 + pred_strd2 + pred_strd));

    pred_8x16b_0 = _mm_cvtepu8_epi16(pred_16x8b_0);
    pred_8x16b_1 = _mm_cvtepu8_epi16(pred_16x8b_1);
    pred_8x16b_2 = _mm_cvtepu8_epi16(pred_16x8b_2);
    pred_8x16b_3 = _mm_cvtepu8_epi16(pred_16x8b_3);
    pred_8x16b_4 = _mm_cvtepu8_epi16(pred_16x8b_4);
    pred_8x16b_5 = _mm_cvtepu8_epi16(pred_16x8b_5);
    pred_8x16b_6 = _mm_cvtepu8_epi16(pred_16x8b_6);
    pred_8x16b_7 = _mm_cvtepu8_epi16(pred_16x8b_7);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    out_8x16b_0 = _mm_add_epi16(pred_8x16b_0, rsd_8x16b_0);
    out_8x16b_1 = _mm_add_epi16(pred_8x16b_1, rsd_8x16b_1);
    out_8x16b_2 = _mm_add_epi16(pred_8x16b_2, rsd_8x16b_2);
    out_8x16b_3 = _mm_add_epi16(pred_8x16b_3, rsd_8x16b_3);
    out_8x16b_4 = _mm_add_epi16(pred_8x16b_4, rsd_8x16b_4);
    out_8x16b_5 = _mm_add_epi16(pred_8x16b_5, rsd_8x16b_5);
    out_8x16b_6 = _mm_add_epi16(pred_8x16b_6, rsd_8x16b_6);
    out_8x16b_7 = _mm_add_epi16(pred_8x16b_7, rsd_8x16b_7);

    out_16x8b_0 = _mm_packus_epi16(out_8x16b_0, zero_8x16b);
    out_16x8b_1 = _mm_packus_epi16(out_8x16b_1, zero_8x16b);
    out_16x8b_2 = _mm_packus_epi16(out_8x16b_2, zero_8x16b);
    out_16x8b_3 = _mm_packus_epi16(out_8x16b_3, zero_8x16b);
    out_16x8b_4 = _mm_packus_epi16(out_8x16b_4, zero_8x16b);
    out_16x8b_5 = _mm_packus_epi16(out_8x16b_5, zero_8x16b);
    out_16x8b_6 = _mm_packus_epi16(out_8x16b_6, zero_8x16b);
    out_16x8b_7 = _mm_packus_epi16(out_8x16b_7, zero_8x16b);

    _mm_storel_epi64((__m128i *) (pu1_out), out_16x8b_0);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd), out_16x8b_1);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2), out_16x8b_2);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd2 + out_strd), out_16x8b_3);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4), out_16x8b_4);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd), out_16x8b_5);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2), out_16x8b_6);
    _mm_storel_epi64((__m128i *) (pu1_out + out_strd4 + out_strd2 + out_strd), out_16x8b_7);

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0)) << 10;
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 11;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 14;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 15;

    i4_nnz |= (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_chroma_4x4_sse42                */
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

void isvcd_pred_residual_recon_chroma_4x4_sse42(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                                WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    __m128i src_r0, src_r1, src_r2, src_r3;
    __m128i pred_r0, pred_r1, pred_r2, pred_r3;
    __m128i pred0, pred1, pred2, pred3;
    __m128i rsd_r0, rsd_r1, rsd_r2, rsd_r3;
    __m128i zero_16x8b;  // all bits reset to zero
    __m128i chroma_mask_even;
    __m128i chroma_mask_odd;

    zero_16x8b = _mm_setzero_si128();

    rsd_r0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_r1 = _mm_loadu_si128((__m128i *) (pi2_rsd + (1 * rsd_strd)));
    rsd_r2 = _mm_loadu_si128((__m128i *) (pi2_rsd + (2 * rsd_strd)));
    rsd_r3 = _mm_loadu_si128((__m128i *) (pi2_rsd + (3 * rsd_strd)));

    pred_r0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred_r1 = _mm_loadu_si128((__m128i *) (pu1_pred + (1 * pred_strd)));
    pred_r2 = _mm_loadu_si128((__m128i *) (pu1_pred + (2 * pred_strd)));
    pred_r3 = _mm_loadu_si128((__m128i *) (pu1_pred + (3 * pred_strd)));

    src_r0 = _mm_loadu_si128((__m128i *) (pu1_out));
    src_r1 = _mm_loadu_si128((__m128i *) (pu1_out + (1 * out_strd)));
    src_r2 = _mm_loadu_si128((__m128i *) (pu1_out + (2 * out_strd)));
    src_r3 = _mm_loadu_si128((__m128i *) (pu1_out + (3 * out_strd)));

    pred0 = _mm_cvtepu8_epi16(pred_r0);
    pred1 = _mm_cvtepu8_epi16(pred_r1);
    pred2 = _mm_cvtepu8_epi16(pred_r2);
    pred3 = _mm_cvtepu8_epi16(pred_r3);

    pred0 = _mm_add_epi16(pred0, rsd_r0);
    pred1 = _mm_add_epi16(pred1, rsd_r1);
    pred2 = _mm_add_epi16(pred2, rsd_r2);
    pred3 = _mm_add_epi16(pred3, rsd_r3);

    pred0 = _mm_packus_epi16(pred0, zero_16x8b);
    pred1 = _mm_packus_epi16(pred1, zero_16x8b);
    pred2 = _mm_packus_epi16(pred2, zero_16x8b);
    pred3 = _mm_packus_epi16(pred3, zero_16x8b);

    chroma_mask_even = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
                                    0x00, 0xff, 0x00, 0xff, 0x00, 0xff);
    chroma_mask_odd = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x00, 0xff,
                                   0x00, 0xff, 0x00, 0xff, 0x00);

    src_r0 = _mm_and_si128(src_r0, chroma_mask_odd);  // 0 src1 0 src2 0 ...
    src_r1 = _mm_and_si128(src_r1, chroma_mask_odd);
    src_r2 = _mm_and_si128(src_r2, chroma_mask_odd);
    src_r3 = _mm_and_si128(src_r3, chroma_mask_odd);

    pred0 = _mm_and_si128(pred0, chroma_mask_even);  // val 0 val 0 ..
    pred1 = _mm_and_si128(pred1, chroma_mask_even);
    pred2 = _mm_and_si128(pred2, chroma_mask_even);
    pred3 = _mm_and_si128(pred3, chroma_mask_even);

    src_r0 = _mm_add_epi8(src_r0, pred0);  // macro  src1 macro src2 macro ...
    src_r1 = _mm_add_epi8(src_r1, pred1);
    src_r2 = _mm_add_epi8(src_r2, pred2);
    src_r3 = _mm_add_epi8(src_r3, pred3);

    _mm_storel_epi64((__m128i *) (&pu1_out[0]), src_r0);
    _mm_storel_epi64((__m128i *) (&pu1_out[out_strd]), src_r1);
    _mm_storel_epi64((__m128i *) (&pu1_out[2 * out_strd]), src_r2);
    _mm_storel_epi64((__m128i *) (&pu1_out[3 * out_strd]), src_r3);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pred_residual_recon_chroma_8x8_sse42                */
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

void isvcd_pred_residual_recon_chroma_8x8_sse42(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                                WORD32 pred_strd, WORD32 rsd_strd, WORD32 out_strd)
{
    __m128i src_r0, src_r1, src_r2, src_r3, src_r4, src_r5, src_r6, src_r7;
    __m128i pred0, pred1, pred2, pred3, pred4, pred5, pred6, pred7;
    __m128i rsd_r0, rsd_r1, rsd_r2, rsd_r3, rsd_r4, rsd_r5, rsd_r6, rsd_r7;
    __m128i zero_16x8b;  // all bits reset to zero
    __m128i chroma_mask_even;
    __m128i chroma_mask_odd;

    zero_16x8b = _mm_setzero_si128();

    rsd_r0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_r1 = _mm_loadu_si128((__m128i *) (pi2_rsd + (1 * rsd_strd)));
    rsd_r2 = _mm_loadu_si128((__m128i *) (pi2_rsd + (2 * rsd_strd)));
    rsd_r3 = _mm_loadu_si128((__m128i *) (pi2_rsd + (3 * rsd_strd)));
    rsd_r4 = _mm_loadu_si128((__m128i *) (pi2_rsd + (4 * rsd_strd)));
    rsd_r5 = _mm_loadu_si128((__m128i *) (pi2_rsd + (5 * rsd_strd)));
    rsd_r6 = _mm_loadu_si128((__m128i *) (pi2_rsd + (6 * rsd_strd)));
    rsd_r7 = _mm_loadu_si128((__m128i *) (pi2_rsd + (7 * rsd_strd)));

    pred0 = _mm_loadu_si128((__m128i *) (pu1_pred));
    pred1 = _mm_loadu_si128((__m128i *) (pu1_pred + (1 * pred_strd)));
    pred2 = _mm_loadu_si128((__m128i *) (pu1_pred + (2 * pred_strd)));
    pred3 = _mm_loadu_si128((__m128i *) (pu1_pred + (3 * pred_strd)));
    pred4 = _mm_loadu_si128((__m128i *) (pu1_pred + (4 * pred_strd)));
    pred5 = _mm_loadu_si128((__m128i *) (pu1_pred + (5 * pred_strd)));
    pred6 = _mm_loadu_si128((__m128i *) (pu1_pred + (6 * pred_strd)));
    pred7 = _mm_loadu_si128((__m128i *) (pu1_pred + (7 * pred_strd)));

    src_r0 = _mm_loadu_si128((__m128i *) (pu1_out));
    src_r1 = _mm_loadu_si128((__m128i *) (pu1_out + (1 * out_strd)));
    src_r2 = _mm_loadu_si128((__m128i *) (pu1_out + (2 * out_strd)));
    src_r3 = _mm_loadu_si128((__m128i *) (pu1_out + (3 * out_strd)));
    src_r4 = _mm_loadu_si128((__m128i *) (pu1_out + (4 * out_strd)));
    src_r5 = _mm_loadu_si128((__m128i *) (pu1_out + (5 * out_strd)));
    src_r6 = _mm_loadu_si128((__m128i *) (pu1_out + (6 * out_strd)));
    src_r7 = _mm_loadu_si128((__m128i *) (pu1_out + (7 * out_strd)));

    pred0 = _mm_cvtepu8_epi16(pred0);
    pred1 = _mm_cvtepu8_epi16(pred1);
    pred2 = _mm_cvtepu8_epi16(pred2);
    pred3 = _mm_cvtepu8_epi16(pred3);
    pred4 = _mm_cvtepu8_epi16(pred4);
    pred5 = _mm_cvtepu8_epi16(pred5);
    pred6 = _mm_cvtepu8_epi16(pred6);
    pred7 = _mm_cvtepu8_epi16(pred7);

    pred0 = _mm_add_epi16(pred0, rsd_r0);
    pred1 = _mm_add_epi16(pred1, rsd_r1);
    pred2 = _mm_add_epi16(pred2, rsd_r2);
    pred3 = _mm_add_epi16(pred3, rsd_r3);
    pred4 = _mm_add_epi16(pred4, rsd_r4);
    pred5 = _mm_add_epi16(pred5, rsd_r5);
    pred6 = _mm_add_epi16(pred6, rsd_r6);
    pred7 = _mm_add_epi16(pred7, rsd_r7);

    pred0 = _mm_packus_epi16(pred0, zero_16x8b);
    pred1 = _mm_packus_epi16(pred1, zero_16x8b);
    pred2 = _mm_packus_epi16(pred2, zero_16x8b);
    pred3 = _mm_packus_epi16(pred3, zero_16x8b);
    pred4 = _mm_packus_epi16(pred4, zero_16x8b);
    pred5 = _mm_packus_epi16(pred5, zero_16x8b);
    pred6 = _mm_packus_epi16(pred6, zero_16x8b);
    pred7 = _mm_packus_epi16(pred7, zero_16x8b);

    chroma_mask_even = _mm_set_epi8(0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff,
                                    0x00, 0xff, 0x00, 0xff, 0x00, 0xff);
    chroma_mask_odd = _mm_set_epi8(0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff,
                                   0x00, 0xff, 0x00, 0xff, 0x00);

    src_r0 = _mm_and_si128(src_r0, chroma_mask_odd);  // 0 src1 0 src2 0 ...
    src_r1 = _mm_and_si128(src_r1, chroma_mask_odd);
    src_r2 = _mm_and_si128(src_r2, chroma_mask_odd);
    src_r3 = _mm_and_si128(src_r3, chroma_mask_odd);
    src_r4 = _mm_and_si128(src_r4, chroma_mask_odd);
    src_r5 = _mm_and_si128(src_r5, chroma_mask_odd);
    src_r6 = _mm_and_si128(src_r6, chroma_mask_odd);
    src_r7 = _mm_and_si128(src_r7, chroma_mask_odd);

    pred0 = _mm_and_si128(pred0, chroma_mask_even);  // val 0 val 0 ..
    pred1 = _mm_and_si128(pred1, chroma_mask_even);
    pred2 = _mm_and_si128(pred2, chroma_mask_even);
    pred3 = _mm_and_si128(pred3, chroma_mask_even);
    pred4 = _mm_and_si128(pred4, chroma_mask_even);
    pred5 = _mm_and_si128(pred5, chroma_mask_even);
    pred6 = _mm_and_si128(pred6, chroma_mask_even);
    pred7 = _mm_and_si128(pred7, chroma_mask_even);

    src_r0 = _mm_add_epi8(src_r0, pred0);  // macro  src1 macro src2 macro ...
    src_r1 = _mm_add_epi8(src_r1, pred1);
    src_r2 = _mm_add_epi8(src_r2, pred2);
    src_r3 = _mm_add_epi8(src_r3, pred3);
    src_r4 = _mm_add_epi8(src_r4, pred4);
    src_r5 = _mm_add_epi8(src_r5, pred5);
    src_r6 = _mm_add_epi8(src_r6, pred6);
    src_r7 = _mm_add_epi8(src_r7, pred7);

    _mm_storel_epi64((__m128i *) (&pu1_out[0]), src_r0);
    _mm_storel_epi64((__m128i *) (&pu1_out[out_strd]), src_r1);
    _mm_storel_epi64((__m128i *) (&pu1_out[2 * out_strd]), src_r2);
    _mm_storel_epi64((__m128i *) (&pu1_out[3 * out_strd]), src_r3);
    _mm_storel_epi64((__m128i *) (&pu1_out[4 * out_strd]), src_r4);
    _mm_storel_epi64((__m128i *) (&pu1_out[5 * out_strd]), src_r5);
    _mm_storel_epi64((__m128i *) (&pu1_out[6 * out_strd]), src_r6);
    _mm_storel_epi64((__m128i *) (&pu1_out[7 * out_strd]), src_r7);

    /* load and repeat for the last 4 elements interleaved in the row */

    rsd_r0 = _mm_loadu_si128((__m128i *) (pi2_rsd + 8));
    rsd_r1 = _mm_loadu_si128((__m128i *) (pi2_rsd + (1 * rsd_strd) + 8));
    rsd_r2 = _mm_loadu_si128((__m128i *) (pi2_rsd + (2 * rsd_strd) + 8));
    rsd_r3 = _mm_loadu_si128((__m128i *) (pi2_rsd + (3 * rsd_strd) + 8));
    rsd_r4 = _mm_loadu_si128((__m128i *) (pi2_rsd + (4 * rsd_strd) + 8));
    rsd_r5 = _mm_loadu_si128((__m128i *) (pi2_rsd + (5 * rsd_strd) + 8));
    rsd_r6 = _mm_loadu_si128((__m128i *) (pi2_rsd + (6 * rsd_strd) + 8));
    rsd_r7 = _mm_loadu_si128((__m128i *) (pi2_rsd + (7 * rsd_strd) + 8));

    pred0 = _mm_loadu_si128((__m128i *) (pu1_pred + 8));
    pred1 = _mm_loadu_si128((__m128i *) (pu1_pred + (1 * pred_strd) + 8));
    pred2 = _mm_loadu_si128((__m128i *) (pu1_pred + (2 * pred_strd) + 8));
    pred3 = _mm_loadu_si128((__m128i *) (pu1_pred + (3 * pred_strd) + 8));
    pred4 = _mm_loadu_si128((__m128i *) (pu1_pred + (4 * pred_strd) + 8));
    pred5 = _mm_loadu_si128((__m128i *) (pu1_pred + (5 * pred_strd) + 8));
    pred6 = _mm_loadu_si128((__m128i *) (pu1_pred + (6 * pred_strd) + 8));
    pred7 = _mm_loadu_si128((__m128i *) (pu1_pred + (7 * pred_strd) + 8));

    src_r0 = _mm_loadu_si128((__m128i *) (pu1_out + 8));
    src_r1 = _mm_loadu_si128((__m128i *) (pu1_out + (1 * out_strd) + 8));
    src_r2 = _mm_loadu_si128((__m128i *) (pu1_out + (2 * out_strd) + 8));
    src_r3 = _mm_loadu_si128((__m128i *) (pu1_out + (3 * out_strd) + 8));
    src_r4 = _mm_loadu_si128((__m128i *) (pu1_out + (4 * out_strd) + 8));
    src_r5 = _mm_loadu_si128((__m128i *) (pu1_out + (5 * out_strd) + 8));
    src_r6 = _mm_loadu_si128((__m128i *) (pu1_out + (6 * out_strd) + 8));
    src_r7 = _mm_loadu_si128((__m128i *) (pu1_out + (7 * out_strd) + 8));

    pred0 = _mm_cvtepu8_epi16(pred0);
    pred1 = _mm_cvtepu8_epi16(pred1);
    pred2 = _mm_cvtepu8_epi16(pred2);
    pred3 = _mm_cvtepu8_epi16(pred3);
    pred4 = _mm_cvtepu8_epi16(pred4);
    pred5 = _mm_cvtepu8_epi16(pred5);
    pred6 = _mm_cvtepu8_epi16(pred6);
    pred7 = _mm_cvtepu8_epi16(pred7);

    pred0 = _mm_add_epi16(pred0, rsd_r0);
    pred1 = _mm_add_epi16(pred1, rsd_r1);
    pred2 = _mm_add_epi16(pred2, rsd_r2);
    pred3 = _mm_add_epi16(pred3, rsd_r3);
    pred4 = _mm_add_epi16(pred4, rsd_r4);
    pred5 = _mm_add_epi16(pred5, rsd_r5);
    pred6 = _mm_add_epi16(pred6, rsd_r6);
    pred7 = _mm_add_epi16(pred7, rsd_r7);

    pred0 = _mm_packus_epi16(pred0, zero_16x8b);
    pred1 = _mm_packus_epi16(pred1, zero_16x8b);
    pred2 = _mm_packus_epi16(pred2, zero_16x8b);
    pred3 = _mm_packus_epi16(pred3, zero_16x8b);
    pred4 = _mm_packus_epi16(pred4, zero_16x8b);
    pred5 = _mm_packus_epi16(pred5, zero_16x8b);
    pred6 = _mm_packus_epi16(pred6, zero_16x8b);
    pred7 = _mm_packus_epi16(pred7, zero_16x8b);

    src_r0 = _mm_and_si128(src_r0, chroma_mask_odd);  // 0 src1 0 src2 0 ...
    src_r1 = _mm_and_si128(src_r1, chroma_mask_odd);
    src_r2 = _mm_and_si128(src_r2, chroma_mask_odd);
    src_r3 = _mm_and_si128(src_r3, chroma_mask_odd);
    src_r4 = _mm_and_si128(src_r4, chroma_mask_odd);
    src_r5 = _mm_and_si128(src_r5, chroma_mask_odd);
    src_r6 = _mm_and_si128(src_r6, chroma_mask_odd);
    src_r7 = _mm_and_si128(src_r7, chroma_mask_odd);

    pred0 = _mm_and_si128(pred0, chroma_mask_even);  // val 0 val 0 ..
    pred1 = _mm_and_si128(pred1, chroma_mask_even);
    pred2 = _mm_and_si128(pred2, chroma_mask_even);
    pred3 = _mm_and_si128(pred3, chroma_mask_even);
    pred4 = _mm_and_si128(pred4, chroma_mask_even);
    pred5 = _mm_and_si128(pred5, chroma_mask_even);
    pred6 = _mm_and_si128(pred6, chroma_mask_even);
    pred7 = _mm_and_si128(pred7, chroma_mask_even);

    src_r0 = _mm_add_epi8(src_r0, pred0);  // macro  src1 macro src2 macro ...
    src_r1 = _mm_add_epi8(src_r1, pred1);
    src_r2 = _mm_add_epi8(src_r2, pred2);
    src_r3 = _mm_add_epi8(src_r3, pred3);
    src_r4 = _mm_add_epi8(src_r4, pred4);
    src_r5 = _mm_add_epi8(src_r5, pred5);
    src_r6 = _mm_add_epi8(src_r6, pred6);
    src_r7 = _mm_add_epi8(src_r7, pred7);

    _mm_storel_epi64((__m128i *) (&pu1_out[0] + 8), src_r0);
    _mm_storel_epi64((__m128i *) (&pu1_out[out_strd] + 8), src_r1);
    _mm_storel_epi64((__m128i *) (&pu1_out[(2 * out_strd)] + 8), src_r2);
    _mm_storel_epi64((__m128i *) (&pu1_out[(3 * out_strd)] + 8), src_r3);
    _mm_storel_epi64((__m128i *) (&pu1_out[(4 * out_strd)] + 8), src_r4);
    _mm_storel_epi64((__m128i *) (&pu1_out[(5 * out_strd)] + 8), src_r5);
    _mm_storel_epi64((__m128i *) (&pu1_out[(6 * out_strd)] + 8), src_r6);
    _mm_storel_epi64((__m128i *) (&pu1_out[(7 * out_strd)] + 8), src_r7);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_4x4_sse42                             */
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

WORD32 isvcd_residual_luma_4x4_sse42(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    __m128i rsd_8x16b_0;
    __m128i rsd_8x16b_1;
    __m128i rsd_8x16b_2;
    __m128i rsd_8x16b_3;
    __m128i rsd_8x16b_01, rsd_8x16b_23;

    __m128i zero_8x16b = _mm_setzero_si128();
    WORD32 i4_nnz, row_01, row_23;

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + (rsd_strd << 1)));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + (rsd_strd << 1) + rsd_strd));

    rsd_8x16b_01 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);

    row_01 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01, zero_8x16b));  // return 1 if all zeros, else 0
    row_23 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23, zero_8x16b));

    i4_nnz = !(row_01 && row_23);
    return i4_nnz;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_8x8_sse42                             */
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

WORD32 isvcd_residual_luma_8x8_sse42(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    __m128i rsd_8x16b_0;
    __m128i rsd_8x16b_1;
    __m128i rsd_8x16b_2;
    __m128i rsd_8x16b_3;
    __m128i rsd_8x16b_4;
    __m128i rsd_8x16b_5;
    __m128i rsd_8x16b_6;
    __m128i rsd_8x16b_7;
    __m128i rsd_8x16b_01_b0, rsd_8x16b_23_b0, rsd_8x16b_45_b2, rsd_8x16b_67_b2;
    __m128i rsd_8x16b_01_b1, rsd_8x16b_23_b1, rsd_8x16b_45_b3, rsd_8x16b_67_b3;

    WORD32 row_01_b0, row_23_b0, row_45_b2, row_67_b2;
    WORD32 row_01_b1, row_23_b1, row_45_b3, row_67_b3;
    WORD32 i4_nnz, i4_nnz_b0, i4_nnz_b1, i4_nnz_b2, i4_nnz_b3;

    __m128i zero_8x16b = _mm_setzero_si128();

    WORD32 rsd_strd2 = (rsd_strd << 1);
    WORD32 rsd_strd4 = (rsd_strd << 2);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0));
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 1;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 4;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 5;

    i4_nnz = (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_luma_16x16_sse42                           */
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

WORD32 isvcd_residual_luma_16x16_sse42(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    __m128i rsd_8x16b_0;
    __m128i rsd_8x16b_1;
    __m128i rsd_8x16b_2;
    __m128i rsd_8x16b_3;
    __m128i rsd_8x16b_4;
    __m128i rsd_8x16b_5;
    __m128i rsd_8x16b_6;
    __m128i rsd_8x16b_7;
    __m128i rsd_8x16b_01_b0, rsd_8x16b_23_b0, rsd_8x16b_45_b2, rsd_8x16b_67_b2;
    __m128i rsd_8x16b_01_b1, rsd_8x16b_23_b1, rsd_8x16b_45_b3, rsd_8x16b_67_b3;

    WORD32 row_01_b0, row_23_b0, row_45_b2, row_67_b2;
    WORD32 row_01_b1, row_23_b1, row_45_b3, row_67_b3;
    WORD32 i4_nnz, i4_nnz_b0, i4_nnz_b1, i4_nnz_b2, i4_nnz_b3;

    __m128i zero_8x16b = _mm_setzero_si128();

    WORD32 rsd_strd2 = (rsd_strd << 1);
    WORD32 rsd_strd4 = (rsd_strd << 2);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0));
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 1;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 4;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 5;

    i4_nnz = (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);

    pi2_rsd += 8;

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0)) << 2;
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 3;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 6;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 7;

    i4_nnz |= (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);

    pi2_rsd -= 8;
    pi2_rsd += (rsd_strd << 3);

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0)) << 8;
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 9;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 12;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 13;

    i4_nnz |= (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);

    pi2_rsd += 8;

    rsd_8x16b_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_2 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_3 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_4 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_5 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_6 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_7 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_01_b0 = _mm_unpacklo_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b0 = _mm_unpacklo_epi64(rsd_8x16b_2, rsd_8x16b_3);
    rsd_8x16b_01_b1 = _mm_unpackhi_epi64(rsd_8x16b_0, rsd_8x16b_1);
    rsd_8x16b_23_b1 = _mm_unpackhi_epi64(rsd_8x16b_2, rsd_8x16b_3);

    rsd_8x16b_45_b2 = _mm_unpacklo_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b2 = _mm_unpacklo_epi64(rsd_8x16b_6, rsd_8x16b_7);
    rsd_8x16b_45_b3 = _mm_unpackhi_epi64(rsd_8x16b_4, rsd_8x16b_5);
    rsd_8x16b_67_b3 = _mm_unpackhi_epi64(rsd_8x16b_6, rsd_8x16b_7);

    row_01_b0 = _mm_test_all_ones(
        _mm_cmpeq_epi16(rsd_8x16b_01_b0, zero_8x16b));  // return 1 if all zeros, else 0
    row_23_b0 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b0, zero_8x16b));
    row_01_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_01_b1, zero_8x16b));
    row_23_b1 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_23_b1, zero_8x16b));
    row_45_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b2, zero_8x16b));
    row_67_b2 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b2, zero_8x16b));
    row_45_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_45_b3, zero_8x16b));
    row_67_b3 = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_67_b3, zero_8x16b));

    i4_nnz_b0 = (!(row_01_b0 && row_23_b0)) << 10;
    i4_nnz_b1 = (!(row_01_b1 && row_23_b1)) << 11;
    i4_nnz_b2 = (!(row_45_b2 && row_67_b2)) << 14;
    i4_nnz_b3 = (!(row_45_b3 && row_67_b3)) << 15;

    i4_nnz |= (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);
    return i4_nnz;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_residual_chroma_cb_cr_8x8_sse42                     */
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

WORD32 isvcd_residual_chroma_cb_cr_8x8_sse42(WORD16 *pi2_rsd, WORD32 rsd_strd)
{
    __m128i rsd_8x16b_r0_0, rsd_8x16b_r0_1, mix_8x16b_r01_0_l, mix_8x16b_r01_1_l,
        rsd_8x16b_r01_b0_cb, rsd_8x16b_r01_b1_cb;
    __m128i rsd_8x16b_r1_0, rsd_8x16b_r1_1, mix_8x16b_r23_0_l, mix_8x16b_r23_1_l,
        rsd_8x16b_r01_b0_cr, rsd_8x16b_r01_b1_cr;
    __m128i rsd_8x16b_r2_0, rsd_8x16b_r2_1, mix_8x16b_r45_0_l, mix_8x16b_r45_1_l,
        rsd_8x16b_r23_b0_cb, rsd_8x16b_r23_b1_cb;
    __m128i rsd_8x16b_r3_0, rsd_8x16b_r3_1, mix_8x16b_r67_0_l, mix_8x16b_r67_1_l,
        rsd_8x16b_r23_b0_cr, rsd_8x16b_r23_b1_cr;
    __m128i rsd_8x16b_r4_0, rsd_8x16b_r4_1, mix_8x16b_r01_0_h, mix_8x16b_r01_1_h,
        rsd_8x16b_r45_b2_cb, rsd_8x16b_r45_b3_cb;
    __m128i rsd_8x16b_r5_0, rsd_8x16b_r5_1, mix_8x16b_r23_0_h, mix_8x16b_r23_1_h,
        rsd_8x16b_r45_b2_cr, rsd_8x16b_r45_b3_cr;
    __m128i rsd_8x16b_r6_0, rsd_8x16b_r6_1, mix_8x16b_r45_0_h, mix_8x16b_r45_1_h,
        rsd_8x16b_r67_b2_cb, rsd_8x16b_r67_b3_cb;
    __m128i rsd_8x16b_r7_0, rsd_8x16b_r7_1, mix_8x16b_r67_0_h, mix_8x16b_r67_1_h,
        rsd_8x16b_r67_b2_cr, rsd_8x16b_r67_b3_cr;

    WORD32 r01_b0_cb, r01_b0_cr;
    WORD32 r23_b0_cb, r23_b0_cr;
    WORD32 r01_b1_cb, r01_b1_cr;
    WORD32 r23_b1_cb, r23_b1_cr;
    WORD32 r45_b2_cb, r45_b2_cr;
    WORD32 r67_b2_cb, r67_b2_cr;
    WORD32 r45_b3_cb, r45_b3_cr;
    WORD32 r67_b3_cb, r67_b3_cr;

    WORD32 i4_nnz, i4_nnz_b0, i4_nnz_b1, i4_nnz_b2, i4_nnz_b3;

    __m128i zero_8x16b = _mm_setzero_si128();

    WORD32 rsd_strd2 = (rsd_strd << 1);
    WORD32 rsd_strd4 = (rsd_strd << 2);

    rsd_8x16b_r0_0 = _mm_loadu_si128((__m128i *) (pi2_rsd));
    rsd_8x16b_r1_0 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd));
    rsd_8x16b_r2_0 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2));
    rsd_8x16b_r3_0 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd));
    rsd_8x16b_r4_0 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4));
    rsd_8x16b_r5_0 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd));
    rsd_8x16b_r6_0 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2));
    rsd_8x16b_r7_0 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd));

    rsd_8x16b_r0_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + 8));
    rsd_8x16b_r1_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd + 8));
    rsd_8x16b_r2_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + 8));
    rsd_8x16b_r3_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd2 + rsd_strd + 8));
    rsd_8x16b_r4_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + 8));
    rsd_8x16b_r5_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd + 8));
    rsd_8x16b_r6_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + 8));
    rsd_8x16b_r7_1 = _mm_loadu_si128((__m128i *) (pi2_rsd + rsd_strd4 + rsd_strd2 + rsd_strd + 8));

    mix_8x16b_r01_0_l =
        _mm_unpacklo_epi16(rsd_8x16b_r0_0, rsd_8x16b_r1_0);  // a0, b0 a1 b1 a2 b2 a3 b3
    mix_8x16b_r23_0_l = _mm_unpacklo_epi16(rsd_8x16b_r2_0, rsd_8x16b_r3_0);
    mix_8x16b_r45_0_l = _mm_unpacklo_epi16(rsd_8x16b_r4_0, rsd_8x16b_r5_0);
    mix_8x16b_r67_0_l = _mm_unpacklo_epi16(rsd_8x16b_r6_0, rsd_8x16b_r7_0);
    mix_8x16b_r01_0_h =
        _mm_unpackhi_epi16(rsd_8x16b_r0_0, rsd_8x16b_r1_0);  // a4 b4 a5 b5 a6 b6 a7 b7
    mix_8x16b_r23_0_h = _mm_unpackhi_epi16(rsd_8x16b_r2_0, rsd_8x16b_r3_0);
    mix_8x16b_r45_0_h = _mm_unpackhi_epi16(rsd_8x16b_r4_0, rsd_8x16b_r5_0);
    mix_8x16b_r67_0_h = _mm_unpackhi_epi16(rsd_8x16b_r6_0, rsd_8x16b_r7_0);

    mix_8x16b_r01_1_l =
        _mm_unpacklo_epi16(rsd_8x16b_r0_1, rsd_8x16b_r1_1);  // a8, b8 a9 b9 a10 b10 a11 b11
    mix_8x16b_r23_1_l = _mm_unpacklo_epi16(rsd_8x16b_r2_1, rsd_8x16b_r3_1);
    mix_8x16b_r45_1_l = _mm_unpacklo_epi16(rsd_8x16b_r4_1, rsd_8x16b_r5_1);
    mix_8x16b_r67_1_l = _mm_unpacklo_epi16(rsd_8x16b_r6_1, rsd_8x16b_r7_1);
    mix_8x16b_r01_1_h =
        _mm_unpackhi_epi16(rsd_8x16b_r0_1, rsd_8x16b_r1_1);  // a12 b12 a13 b13 a14 b14 a15 b15
    mix_8x16b_r23_1_h = _mm_unpackhi_epi16(rsd_8x16b_r2_1, rsd_8x16b_r3_1);
    mix_8x16b_r45_1_h = _mm_unpackhi_epi16(rsd_8x16b_r4_1, rsd_8x16b_r5_1);
    mix_8x16b_r67_1_h = _mm_unpackhi_epi16(rsd_8x16b_r6_1, rsd_8x16b_r7_1);

    mix_8x16b_r01_0_l = _mm_shuffle_epi32(mix_8x16b_r01_0_l, 0b11011000);  // a0b0 a2b2 a1b1 a3b3
    mix_8x16b_r23_0_l = _mm_shuffle_epi32(mix_8x16b_r23_0_l, 0b11011000);  // c0d0
    mix_8x16b_r45_0_l = _mm_shuffle_epi32(mix_8x16b_r45_0_l, 0b11011000);  // e0f0
    mix_8x16b_r67_0_l = _mm_shuffle_epi32(mix_8x16b_r67_0_l, 0b11011000);  // g0h0
    mix_8x16b_r01_0_h = _mm_shuffle_epi32(mix_8x16b_r01_0_h, 0b11011000);  // a4b4 a6b6 a5b5 a7b7
    mix_8x16b_r23_0_h = _mm_shuffle_epi32(mix_8x16b_r23_0_h, 0b11011000);  // c4d4
    mix_8x16b_r45_0_h = _mm_shuffle_epi32(mix_8x16b_r45_0_h, 0b11011000);  // e4f4
    mix_8x16b_r67_0_h = _mm_shuffle_epi32(mix_8x16b_r67_0_h, 0b11011000);  // g4h4

    mix_8x16b_r01_1_l = _mm_shuffle_epi32(mix_8x16b_r01_1_l, 0b11011000);
    mix_8x16b_r23_1_l = _mm_shuffle_epi32(mix_8x16b_r23_1_l, 0b11011000);
    mix_8x16b_r45_1_l = _mm_shuffle_epi32(mix_8x16b_r45_1_l, 0b11011000);
    mix_8x16b_r67_1_l = _mm_shuffle_epi32(mix_8x16b_r67_1_l, 0b11011000);
    mix_8x16b_r01_1_h = _mm_shuffle_epi32(mix_8x16b_r01_1_h, 0b11011000);
    mix_8x16b_r23_1_h = _mm_shuffle_epi32(mix_8x16b_r23_1_h, 0b11011000);
    mix_8x16b_r45_1_h = _mm_shuffle_epi32(mix_8x16b_r45_1_h, 0b11011000);
    mix_8x16b_r67_1_h = _mm_shuffle_epi32(mix_8x16b_r67_1_h, 0b11011000);

    rsd_8x16b_r01_b0_cb =
        _mm_unpacklo_epi64(mix_8x16b_r01_0_l, mix_8x16b_r01_0_h);  // a0b0 a2b2 a4b4 a6b6
    rsd_8x16b_r01_b0_cr =
        _mm_unpackhi_epi64(mix_8x16b_r01_0_l, mix_8x16b_r01_0_h);  // a1b1 a3b3 a5b5 a7b7
    rsd_8x16b_r23_b0_cb = _mm_unpacklo_epi64(mix_8x16b_r23_0_l, mix_8x16b_r23_0_h);  //
    rsd_8x16b_r23_b0_cr = _mm_unpackhi_epi64(mix_8x16b_r23_0_l, mix_8x16b_r23_0_h);
    rsd_8x16b_r45_b2_cb = _mm_unpacklo_epi64(mix_8x16b_r45_0_l, mix_8x16b_r45_0_h);
    rsd_8x16b_r45_b2_cr = _mm_unpackhi_epi64(mix_8x16b_r45_0_l, mix_8x16b_r45_0_h);
    rsd_8x16b_r67_b2_cb = _mm_unpacklo_epi64(mix_8x16b_r67_0_l, mix_8x16b_r67_0_h);
    rsd_8x16b_r67_b2_cr = _mm_unpackhi_epi64(mix_8x16b_r67_0_l, mix_8x16b_r67_0_h);

    rsd_8x16b_r01_b1_cb =
        _mm_unpacklo_epi64(mix_8x16b_r01_1_l, mix_8x16b_r01_1_h);  // a8b8 a10b10 a12b12 a14b14
    rsd_8x16b_r01_b1_cr =
        _mm_unpackhi_epi64(mix_8x16b_r01_1_l, mix_8x16b_r01_1_h);  // a9b9 a11b11 a13b13 a15b15
    rsd_8x16b_r23_b1_cb = _mm_unpacklo_epi64(mix_8x16b_r23_1_l, mix_8x16b_r23_1_h);
    rsd_8x16b_r23_b1_cr = _mm_unpackhi_epi64(mix_8x16b_r23_1_l, mix_8x16b_r23_1_h);
    rsd_8x16b_r45_b3_cb = _mm_unpacklo_epi64(mix_8x16b_r45_1_l, mix_8x16b_r45_1_h);
    rsd_8x16b_r45_b3_cr = _mm_unpackhi_epi64(mix_8x16b_r45_1_l, mix_8x16b_r45_1_h);
    rsd_8x16b_r67_b3_cb = _mm_unpacklo_epi64(mix_8x16b_r67_1_l, mix_8x16b_r67_1_h);
    rsd_8x16b_r67_b3_cr = _mm_unpackhi_epi64(mix_8x16b_r67_1_l, mix_8x16b_r67_1_h);

    r01_b0_cb = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r01_b0_cb, zero_8x16b));
    r23_b0_cb = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r23_b0_cb, zero_8x16b));
    r01_b1_cb = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r01_b1_cb, zero_8x16b));
    r23_b1_cb = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r23_b1_cb, zero_8x16b));
    r45_b2_cb = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r45_b2_cb, zero_8x16b));
    r67_b2_cb = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r67_b2_cb, zero_8x16b));
    r45_b3_cb = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r45_b3_cb, zero_8x16b));
    r67_b3_cb = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r67_b3_cb, zero_8x16b));

    r01_b0_cr = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r01_b0_cr, zero_8x16b));
    r23_b0_cr = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r23_b0_cr, zero_8x16b));
    r01_b1_cr = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r01_b1_cr, zero_8x16b));
    r23_b1_cr = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r23_b1_cr, zero_8x16b));
    r45_b2_cr = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r45_b2_cr, zero_8x16b));
    r67_b2_cr = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r67_b2_cr, zero_8x16b));
    r45_b3_cr = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r45_b3_cr, zero_8x16b));
    r67_b3_cr = _mm_test_all_ones(_mm_cmpeq_epi16(rsd_8x16b_r67_b3_cr, zero_8x16b));

    i4_nnz_b0 = (!(r01_b0_cr && r23_b0_cr));
    i4_nnz_b1 = (!(r01_b1_cr && r23_b1_cr)) << 1;
    i4_nnz_b2 = (!(r45_b2_cr && r67_b2_cr)) << 2;
    i4_nnz_b3 = (!(r45_b3_cr && r67_b3_cr)) << 3;

    i4_nnz = (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);
    i4_nnz = i4_nnz << 4;

    i4_nnz_b0 = (!(r01_b0_cb && r23_b0_cb));
    i4_nnz_b1 = (!(r01_b1_cb && r23_b1_cb)) << 1;
    i4_nnz_b2 = (!(r45_b2_cb && r67_b2_cb)) << 2;
    i4_nnz_b3 = (!(r45_b3_cb && r67_b3_cb)) << 3;

    i4_nnz |= (i4_nnz_b0 | i4_nnz_b1 | i4_nnz_b2 | i4_nnz_b3);
    return i4_nnz;
}
