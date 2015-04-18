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
******************************************************************************
* @file ih264_cabac_tables.h
*
* @brief
*  This file contains enumerations, macros and extern declarations of H264
*  cabac tables
*
* @author
*  Ittiam
*
* @remarks
*  none
******************************************************************************
*/

#ifndef IH264_CABAC_TABLES_H_
#define IH264_CABAC_TABLES_H_

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

/**
******************************************************************************
 *  @brief  maximum range of cabac_init_idc (0-2)
******************************************************************************
 */
#define IH264_NUM_CABAC_INIT_IDC_PLUS_ONE   4

/**
******************************************************************************
 *  @brief  max range of qps in H264 (0-51)
******************************************************************************
 */
#define IH264_MAX_QP      52

/**
******************************************************************************
 *  @brief  max range of cabac contexts in H264 (0-459)
******************************************************************************
 */
#define IH264_NUM_CABAC_CTXTS 460

/*****************************************************************************/
/* Extern global declarations                                                */
/*****************************************************************************/

/**
 ******************************************************************************
 * @breif  Table for rangeTabLPS depending on pStateIdx and qCodIRangeIdx
 * input   : pStateIdx(0-63) and qCodIRangeIdx(0-3) [(Range >> 6) & 0x3]
 * output  : RLps
 *
 * @remarks See Table 9-35 of H264 spec for rangeTabLPS
 *******************************************************************************
 */
extern const UWORD8 gau1_ih264_cabac_rlps[64][4];


/**
 ******************************************************************************
 * @breif  probability+MPS state transition tables based on cur State and bin
 * input  : curpState[bits7-2]  | curMPS[bit1] | decodedBin[bit0]
 * output : nextpState[bits6-1] | nextMPS[bit0]
 * @remarks Modified form of Table-9-36 State Transition table in H264 spec
 ******************************************************************************
 */
extern const UWORD8 gau1_ih264_next_state[128*2];


/**
 ******************************************************************************
 * @brief  Init context tables for all combinations of qp and cabac_init_idc
 * @remarks Packing format MPS in lsb and pState in bits[1-6]
 ******************************************************************************
 */
extern const UWORD8 gau1_ih264_cab_ctxts[IH264_NUM_CABAC_INIT_IDC_PLUS_ONE][IH264_MAX_QP][IH264_NUM_CABAC_CTXTS];


#endif /* IH264_CABAC_TABLES_H_ */
