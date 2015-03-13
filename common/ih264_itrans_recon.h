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
*  ih264_itrans_recon.h
*
* @brief
*  Contains function declarations for inverse transform  and reconstruction of
*  the quantized macro blocks
*
* @author
*  Ittiam
*
* @par List of Functions:
*  - ih264_itrans_recon_ft
*  - ih264_itrans_recon_4x4
*  - ih264_itrans_recon_8x8
*  - ih264_itrans_recon_4x4_a9
*
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef IH264_ITRANS_RECON_H_
#define IH264_ITRANS_RECON_H_

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/

typedef void ih264_itrans_recon_ft(WORD16 *pi2_src,
                                   UWORD8 *pu1_pred,
                                   UWORD8 *pu1_recon,
                                   WORD32 src_strd,
                                   WORD32 pred_strd,
                                   WORD32 dst_strd,
                                   UWORD32 q_lev,
                                   WORD32 *pi4_tmp);

/*C declarations*/

ih264_itrans_recon_ft ih264_itrans_recon_4x4;

ih264_itrans_recon_ft ih264_itrans_recon_8x8;

/*A9 declarations */

ih264_itrans_recon_ft ih264_itrans_recon_4x4_a9;

#endif /* IH264_ITRANS_RECON_H_ */
