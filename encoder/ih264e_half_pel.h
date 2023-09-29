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
*  ih264e_half_pel.h
*
* @brief
*  Contains declarations of subpel functions used by the encoder
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_HALF_PEL_H_
#define _IH264E_HALF_PEL_H_

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/
/*
 * Dimensions of subpel plane buffers
 */
#define HP_PL_WD  MB_SIZE + 1
#define HP_PL_HT  MB_SIZE + 1

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

typedef void ih264e_sixtapfilter_horz_ft(UWORD8 *pu1_src,
                                         UWORD8 *pu1_dst,
                                         WORD32 src_strd,
                                         WORD32 dst_strd);

typedef void ih264e_sixtap_filter_2dvh_vert_ft(UWORD8 *pu1_src,
                                               UWORD8 *pu1_dst1,
                                               UWORD8 *pu1_dst2,
                                               WORD32 src_strd,
                                               WORD32 dst_strd,
                                               WORD32 *pi4_pred,
                                               WORD32 i4_pred_strd);

/* C Declarations */
ih264e_sixtapfilter_horz_ft ih264e_sixtapfilter_horz;
ih264e_sixtap_filter_2dvh_vert_ft ih264e_sixtap_filter_2dvh_vert;

/* A9 Declarations */
ih264e_sixtapfilter_horz_ft ih264e_sixtapfilter_horz_a9q;
ih264e_sixtap_filter_2dvh_vert_ft ih264e_sixtap_filter_2dvh_vert_a9q;

/* AV8 Declarations */
ih264e_sixtapfilter_horz_ft ih264e_sixtapfilter_horz_av8;
ih264e_sixtap_filter_2dvh_vert_ft ih264e_sixtap_filter_2dvh_vert_av8;

/* SSSE3 Declarations */
ih264e_sixtapfilter_horz_ft ih264e_sixtapfilter_horz_ssse3;
ih264e_sixtap_filter_2dvh_vert_ft ih264e_sixtap_filter_2dvh_vert_ssse3;

#endif /* _IH264E_HALF_PEL_H_ */
