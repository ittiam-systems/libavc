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
 *  isvcd_pred_residual_recon.h
 *
 * @brief
 *  Contains declarations for forward and inverse transform paths for SVC
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_PRED_RESDUAL_RECON_
#define _ISVCD_PRED_RESDUAL_RECON_

/*Function prototype declarations*/

typedef WORD32 ih264_residual_ft(WORD16 *pi2_rsd, WORD32 rsd_stride);

typedef WORD32 ih264_residual_chroma_ft(WORD16 *pi2_rsd, WORD32 rsd_stride);

typedef WORD32 ih264_pred_residual_recon_ft(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                            WORD32 pred_strd, WORD32 rsd_stride, WORD32 out_strd);

typedef void ih264_pred_residual_recon_chroma_ft(UWORD8 *pu1_pred, WORD16 *pi2_rsd, UWORD8 *pu1_out,
                                                 WORD32 pred_strd, WORD32 rsd_stride,
                                                 WORD32 out_strd);

ih264_residual_ft isvcd_residual_luma_4x4;
ih264_residual_ft isvcd_residual_luma_8x8;
ih264_residual_ft isvcd_residual_luma_16x16;
ih264_residual_chroma_ft isvcd_residual_chroma_cb_cr_8x8;

ih264_residual_ft isvcd_residual_luma_4x4_sse42;
ih264_residual_ft isvcd_residual_luma_8x8_sse42;
ih264_residual_ft isvcd_residual_luma_16x16_sse42;
ih264_residual_chroma_ft isvcd_residual_chroma_cb_cr_8x8_sse42;

ih264_residual_ft isvcd_residual_luma_4x4_neonintr;
ih264_residual_ft isvcd_residual_luma_8x8_neonintr;
ih264_residual_ft isvcd_residual_luma_16x16_neonintr;
ih264_residual_chroma_ft isvcd_residual_chroma_cb_cr_8x8_neonintr;

ih264_pred_residual_recon_ft isvcd_pred_residual_recon_4x4;
ih264_pred_residual_recon_ft isvcd_pred_residual_recon_8x8;
ih264_pred_residual_recon_ft isvcd_pred_residual_recon_16x16;
ih264_pred_residual_recon_chroma_ft isvcd_pred_residual_recon_chroma_4x4;
ih264_pred_residual_recon_chroma_ft isvcd_pred_residual_recon_chroma_8x8;

ih264_pred_residual_recon_ft isvcd_pred_residual_recon_4x4_neonintr;
ih264_pred_residual_recon_ft isvcd_pred_residual_recon_8x8_neonintr;
ih264_pred_residual_recon_ft isvcd_pred_residual_recon_16x16_neonintr;
ih264_pred_residual_recon_chroma_ft isvcd_pred_residual_recon_chroma_4x4_neonintr;
ih264_pred_residual_recon_chroma_ft isvcd_pred_residual_recon_chroma_8x8_neonintr;

ih264_pred_residual_recon_ft isvcd_pred_residual_recon_4x4_sse42;
ih264_pred_residual_recon_ft isvcd_pred_residual_recon_8x8_sse42;
ih264_pred_residual_recon_ft isvcd_pred_residual_recon_16x16_sse42;
ih264_pred_residual_recon_chroma_ft isvcd_pred_residual_recon_chroma_4x4_sse42;
ih264_pred_residual_recon_chroma_ft isvcd_pred_residual_recon_chroma_8x8_sse42;

#endif /* _ISVCD_PRED_RESDUAL_RECON_ */