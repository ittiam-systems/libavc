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
 *  isvc_inter_pred_filters.h
 *
 * @brief
 *  Declarations of functions used for inter prediction
 *
 * @author
 *  Ittiam
 *
 * @par List of Functions:
 *  -ih264_inter_pred_luma_copy
 *  -ih264_interleave_copy
 *  -ih264_inter_pred_luma_horz
 *  -ih264_inter_pred_luma_vert
 *  -ih264_inter_pred_luma_horz_hpel_vert_hpel
 *  -ih264_inter_pred_luma_vert_qpel
 *  -ih264_inter_pred_luma_horz_qpel
 *  -ih264_inter_pred_luma_horz_qpel_vert_qpel
 *  -ih264_inter_pred_luma_horz_qpel_vert_hpel
 *  -ih264_inter_pred_luma_horz_hpel_vert_qpel
 *  -ih264_inter_pred_luma_bilinear
 *  -ih264_inter_pred_chroma
 *  -ih264_inter_pred_luma_copy_a9q
 *  -ih264_interleave_copy_a9
 *  -ih264_inter_pred_luma_horz_a9q
 *  -ih264_inter_pred_luma_vert_a9q
 *  -ih264_inter_pred_luma_bilinear_a9q
 *  -ih264_inter_pred_luma_horz_hpel_vert_hpel_a9q
 *  -ih264_inter_pred_luma_horz_qpel_a9q
 *  -ih264_inter_pred_luma_vert_qpel_a9q
 *  -ih264_inter_pred_luma_horz_qpel_vert_qpel_a9q
 *  -ih264_inter_pred_luma_horz_qpel_vert_hpel_a9q
 *  -ih264_inter_pred_luma_horz_hpel_vert_qpel_a9q
 *  -ih264_inter_pred_chroma_a9q
 *  -ih264_inter_pred_luma_copy_av8
 *  -ih264_interleave_copy_av8
 *  -ih264_inter_pred_luma_horz_av8
 *  -ih264_inter_pred_luma_vert_av8
 *  -ih264_inter_pred_luma_bilinear_av8
 *  -ih264_inter_pred_luma_horz_hpel_vert_hpel_av8
 *  -ih264_inter_pred_luma_horz_qpel_av8
 *  -ih264_inter_pred_luma_vert_qpel_av8
 *  -ih264_inter_pred_luma_horz_qpel_vert_qpel_av8
 *  -ih264_inter_pred_luma_horz_qpel_vert_hpel_av8
 *  -ih264_inter_pred_luma_horz_hpel_vert_qpel_av8
 *  -ih264_inter_pred_chroma_av8
 *  -ih264_inter_pred_chroma_dx_zero_av8
 *  -ih264_inter_pred_chroma_dy_zero_av8
 *  -ih264_inter_pred_luma_copy_ssse3
 *  -ih264_inter_pred_luma_copy_ssse3
 *  -ih264_inter_pred_luma_horz_ssse3
 *  -ih264_inter_pred_luma_vert_ssse3
 *  -ih264_inter_pred_luma_bilinear_ssse3
 *  -ih264_inter_pred_luma_horz_hpel_vert_hpel_ssse3
 *  -ih264_inter_pred_luma_horz_qpel_ssse3
 *  -ih264_inter_pred_luma_vert_qpel_ssse3
 *  -ih264_inter_pred_luma_horz_qpel_vert_qpel_ssse3
 *  -ih264_inter_pred_luma_horz_qpel_vert_hpel_ssse3
 *  -ih264_inter_pred_luma_horz_hpel_vert_qpel_ssse3
 *  -ih264_inter_pred_chroma_ssse3
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVC_INTER_PRED_FILTERS_H_
#define _ISVC_INTER_PRED_FILTERS_H_

/*****************************************************************************/
/* Constant Data variables                                                   */
/*****************************************************************************/

extern const WORD32 ih264_g_six_tap[3]; /* coefficients for 6 tap filtering*/

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/

typedef void FT_INTER_PRED_LUMA(UWORD8 *pu1_src, UWORD8 *pu1_dst, WORD32 src_strd, WORD32 dst_strd,
                                WORD32 ht, WORD32 wd, UWORD8 *pu1_tmp, WORD32 dydx);

typedef void FT_INTERLEAVE_COPY(UWORD8 *pu1_src, UWORD8 *pu1_dst, WORD32 src_strd, WORD32 dst_strd,
                                WORD32 ht, WORD32 wd);

typedef void FT_INTER_PRED_LUMA_BILINEAR(UWORD8 *pu1_src1, UWORD8 *pu1_src2, UWORD8 *pu1_dst,
                                         WORD32 src_strd1, WORD32 src_strd2, WORD32 dst_strd,
                                         WORD32 height, WORD32 width);

typedef void FT_INTER_PRED_CHROMA(UWORD8 *pu1_src, UWORD8 *pu1_dst, WORD32 src_strd,
                                  WORD32 dst_strd, WORD32 dx, WORD32 dy, WORD32 ht, WORD32 wd);

/* No NEON Declarations */

FT_INTER_PRED_LUMA ih264_inter_pred_luma_copy;

FT_INTERLEAVE_COPY ih264_interleave_copy;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_vert;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_hpel_vert_hpel;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_vert_qpel;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_vert_qpel;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_vert_hpel;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_hpel_vert_qpel;

FT_INTER_PRED_LUMA_BILINEAR ih264_inter_pred_luma_bilinear;

FT_INTER_PRED_CHROMA ih264_inter_pred_chroma;

/* A9 NEON Declarations */
FT_INTER_PRED_LUMA ih264_inter_pred_luma_copy_a9q;

FT_INTERLEAVE_COPY ih264_interleave_copy_a9;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_a9q;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_vert_a9q;

FT_INTER_PRED_LUMA_BILINEAR ih264_inter_pred_luma_bilinear_a9q;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_hpel_vert_hpel_a9q;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_a9q;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_vert_qpel_a9q;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_vert_qpel_a9q;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_vert_hpel_a9q;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_hpel_vert_qpel_a9q;

FT_INTER_PRED_CHROMA ih264_inter_pred_chroma_a9q;

/* AV8 NEON Declarations */
FT_INTER_PRED_LUMA ih264_inter_pred_luma_copy_av8;

FT_INTERLEAVE_COPY ih264_interleave_copy_av8;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_av8;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_vert_av8;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_hpel_vert_hpel_av8;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_av8;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_vert_qpel_av8;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_vert_qpel_av8;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_vert_hpel_av8;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_hpel_vert_qpel_av8;

FT_INTER_PRED_CHROMA ih264_inter_pred_chroma_av8;

FT_INTER_PRED_CHROMA ih264_inter_pred_chroma_dx_zero_av8;

FT_INTER_PRED_CHROMA ih264_inter_pred_chroma_dy_zero_av8;

/* SSSE3 Intrinsic Declarations */
FT_INTER_PRED_LUMA ih264_inter_pred_luma_copy_ssse3;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_ssse3;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_vert_ssse3;

FT_INTER_PRED_LUMA_BILINEAR ih264_inter_pred_luma_bilinear_ssse3;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_hpel_vert_hpel_ssse3;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_ssse3;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_vert_qpel_ssse3;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_vert_qpel_ssse3;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_qpel_vert_hpel_ssse3;

FT_INTER_PRED_LUMA ih264_inter_pred_luma_horz_hpel_vert_qpel_ssse3;

FT_INTER_PRED_CHROMA ih264_inter_pred_chroma_ssse3;

/** Nothing past this point */

#endif
