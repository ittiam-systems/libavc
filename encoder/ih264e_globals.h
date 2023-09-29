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
*  ih264e_globals.h
*
* @brief
*  Contains declarations of global variables used in the encoder
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_GLOBALS_H_
#define _IH264E_GLOBALS_H_


/*****************************************************************************/
/* Global Declarations                                                       */
/*****************************************************************************/

extern const UWORD8 gu1_qp_lambdaIP[52];
extern const UWORD8 gu1_qp_lambdaB[52];
extern const UWORD8 gu1_qp0[52];
extern const UWORD8 u1_uev_codelength[32];
extern const UWORD8 gu1_coeff_cost[6];
extern const UWORD8 gu1_luma_scan_order[16];
extern const UWORD8 gu1_chroma_scan_order[15];
extern const UWORD8 gu1_luma_scan_order_dc[16];
extern const UWORD8 gu1_chroma_scan_order_dc[4];
extern const WORD8 gi1_mv_pred_condition[8];
extern const UWORD8 gau1_h264_to_mpeg2_qmap[H264_QP_ELEM];
extern const UWORD8 gau1_mpeg2_to_h264_qmap[MPEG2_QP_ELEM];


#endif /* _IH264E_GLOBALS_H_ */
