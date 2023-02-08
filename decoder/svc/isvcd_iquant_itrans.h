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
 *  isvcd_iquant_itrans.h
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

#ifndef ISVCD_ITRANS_IQUANT_
#define ISVCD_ITRANS_IQUANT_

/*Function prototype declarations*/

typedef void ih264_iquant_itrans_ft(WORD16 *pi2_src, WORD16 *pi2_out, WORD32 out_strd,
                                    const UWORD16 *pu2_iscale_mat, const UWORD16 *pu2_weigh_mat,
                                    UWORD32 qp_div, WORD16 *pi2_tmp, WORD32 iq_start_idx,
                                    WORD16 *pi2_dc_ld_addr);

typedef void ih264_iquant_itrans_chroma_ft(WORD16 *pi2_src, WORD16 *pi2_out, WORD32 out_strd,
                                           const UWORD16 *pu2_iscal_mat,
                                           const UWORD16 *pu2_weigh_mat, UWORD32 u4_qp_div_6,
                                           WORD16 *pi2_tmp, WORD16 *pi2_dc_src);

ih264_iquant_itrans_ft isvcd_iquant_itrans_4x4;
ih264_iquant_itrans_ft isvcd_iquant_itrans_8x8;
ih264_iquant_itrans_ft isvcd_iquant_itrans_4x4_dc;
ih264_iquant_itrans_ft isvcd_iquant_itrans_8x8_dc;
ih264_iquant_itrans_chroma_ft isvcd_iquant_itrans_chroma_4x4;
ih264_iquant_itrans_chroma_ft isvcd_iquant_itrans_chroma_4x4_dc;

ih264_iquant_itrans_ft isvcd_iquant_itrans_4x4_dc_neonintr;
ih264_iquant_itrans_ft isvcd_iquant_itrans_8x8_dc_neonintr;
ih264_iquant_itrans_ft isvcd_iquant_itrans_4x4_neonintr;
ih264_iquant_itrans_ft isvcd_iquant_itrans_8x8_neonintr;
ih264_iquant_itrans_chroma_ft isvcd_iquant_itrans_chroma_4x4_dc_neonintr;
ih264_iquant_itrans_chroma_ft isvcd_iquant_itrans_chroma_4x4_neonintr;

ih264_iquant_itrans_ft isvcd_iquant_itrans_4x4_dc_sse42;
ih264_iquant_itrans_ft isvcd_iquant_itrans_8x8_dc_sse42;
ih264_iquant_itrans_ft isvcd_iquant_itrans_4x4_sse42;
ih264_iquant_itrans_ft isvcd_iquant_itrans_8x8_sse42;
ih264_iquant_itrans_chroma_ft isvcd_iquant_itrans_chroma_4x4_dc_sse42;
ih264_iquant_itrans_chroma_ft isvcd_iquant_itrans_chroma_4x4_sse42;

#endif /* ISVCD_ITRANS_IQUANT_ */