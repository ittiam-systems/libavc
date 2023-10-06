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
*  ih264_trans_data.h
*
* @brief
*  Contains declaration of global variables for H264 transform, qunat and
*  inverse quant
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264_TRANS_DATA_H_
#define _IH264_TRANS_DATA_H_

/*****************************************************************************/
/* Extern global declarations                                                */
/*****************************************************************************/

extern const UWORD16 gu2_quant_scale_matrix_4x4[96];

extern const UWORD32 gu4_forward_quant_round_factor_4x4[9];

extern const UWORD16 gu2_forward_quant_threshold_4x4[96];

extern const UWORD16 gu2_quant_scale_matrix_8x8 [384];

extern const UWORD8 gu1_qpc_fqpi[52];

#endif /* _IH264_TRANS_DATA_H_ */
