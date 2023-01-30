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
 *  isvce_cabac_structs.h
 *
 * @brief
 *  This file contains cabac related macros, enums, tables and function
 *declarations.
 *
 * @author
 *  Doney Alex
 *
 * @remarks
 *  none
 *
 *******************************************************************************
 */

#ifndef _ISVCE_CABAC_H_
#define _ISVCE_CABAC_H_

#include "ih264e_cabac.h"
#include "isvce_cabac_structs.h"
#include "isvce_defs.h"
#include "isvce_structs.h"

/*****************************************************************************/
/* Function Declarations                                                 */
/*****************************************************************************/

/**
 *******************************************************************************
 *
 * @brief
 * Initialize default context values and pointers.
 *
 * @param[in] ps_ent_ctxt
 *  Pointer to entropy context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
void isvce_init_cabac_table(isvce_entropy_ctxt_t *ps_ent_ctxt);

/**
 *******************************************************************************
 *
 * @brief
 * Initialize cabac context: Intitalize all contest with init values given in
 *the spec. Called at the beginning of entropy coding of each slice for CABAC
 *encoding.
 *
 * @param[in] ps_ent_ctxt
 *  Pointer to entropy context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
extern void isvce_init_cabac_ctxt(isvce_entropy_ctxt_t *ps_ent_ctxt, slice_header_t *ps_slice_hdr);

/**
 *******************************************************************************
 *
 * @brief
 *  k-th order Exp-Golomb (UEGk) binarization process: Implements concatenated
 *   unary/ k-th order Exp-Golomb  (UEGk) binarization process,
 *   where k = 0 as defined in 9.3.2.3 of  ITU_T_H264-201402
 *
 * @param[in] i2_sufs
 *  Suffix bit string
 *
 * @param[in] pi1_bins_len
 *  Pointer to length of the string
 *
 * @returns Binarized value
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
UWORD32 isvce_cabac_UEGk0_binarization(WORD16 i2_sufs, WORD8 *pi1_bins_len);

/**
 *******************************************************************************
 *
 * @brief
 *  Get cabac context for the MB :calculates the pointers to Top and   left
 *          cabac neighbor context depending upon neighbor  availability.
 *
 * @param[in] ps_ent_ctxt
 *  Pointer to entropy context structure
 *
 * @param[in] u4_mb_type
 *  Type of MB
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
void isvce_get_cabac_context(isvce_entropy_ctxt_t *ps_ent_ctxt, WORD32 u4_mb_type);

/**
 *******************************************************************************
 * @brief
 *  flushing at termination: Explained in flowchart 9-12(ITU_T_H264-201402).
 *
 *  @param[in]   ps_cabac_ctxt
 *  pointer to cabac context (handle)
 *
 * @returns  none
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
void isvce_cabac_flush(isvce_cabac_ctxt_t *ps_cabac_ctxt);

/**
 ******************************************************************************
 *
 *  @brief Puts new byte (and outstanding bytes) into bitstream after cabac
 *         renormalization
 *
 *  @par   Description
 *  1. Extract the leading byte of low(L)
 *  2. If leading byte=0xff increment outstanding bytes and return
 *     (as the actual bits depend on carry propogation later)
 *  3. If leading byte is not 0xff check for any carry propogation
 *  4. Insert the carry (propogated in previous byte) along with outstanding
 *     bytes (if any) and leading byte
 *
 *
 *  @param[inout]   ps_cabac_ctxt
 *  pointer to cabac context (handle)
 *
 *  @return
 *
 ******************************************************************************
 */
void isvce_cabac_put_byte(isvce_cabac_ctxt_t *ps_cabac_ctxt);

/**
 ******************************************************************************
 *
 *  @brief Codes a bin based on probablilty and mps packed context model
 *
 *  @par   Description
 *  1. Apart from encoding bin, context model is updated as per state transition
 *  2. Range and Low renormalization is done based on bin and original state
 *  3. After renorm bistream is updated (if required)
 *
 *  @param[inout]   ps_cabac
 *  pointer to cabac context (handle)
 *
 *  @param[in]   bin
 *  bin(boolean) to be encoded
 *
 *  @param[in]  pu1_bin_ctxts
 *  index of cabac context model containing pState[bits 5-0] | MPS[bit6]
 *
 *  @return
 *
 ******************************************************************************
 */
void isvce_cabac_encode_bin(isvce_cabac_ctxt_t *ps_cabac, WORD32 bin,
                            bin_ctxt_model *pu1_bin_ctxts);

/**
 *******************************************************************************
 *
 * @brief
 *  Encoding process for a binary decision :implements encoding process of a
 decision
 *  as defined in 9.3.4.2 . This function encodes multiple bins, of a symbol.
 Implements
 *  flowchart Figure 9-7( ITU_T_H264-201402)
 *
 * @param[in] u4_bins
 * array of bin values
 *
 * @param[in] i1_bins_len
 *  Length of bins, maximum 32
 *
 * @param[in] u4_ctx_inc
 *  CtxInc, byte0- bin0, byte1-bin1 ..
 *
 * @param[in] i1_valid_len
 *  valid length of bins, after that CtxInc is constant
 *
 * @param[in] pu1_bin_ctxt_type
 *  Pointer to binary contexts

 * @param[in] ps_cabac
 *  Pointer to cabac_context_structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
void isvce_encode_decision_bins(UWORD32 u4_bins, WORD8 i1_bins_len, UWORD32 u4_ctx_inc,
                                WORD8 i1_valid_len, bin_ctxt_model *pu1_bin_ctxt_type,
                                isvce_cabac_ctxt_t *ps_cabac);

/**
 *******************************************************************************
 * @brief
 *  Encoding process for a binary decision before termination:Encoding process
 *  of a termination(9.3.4.5 :ITU_T_H264-201402) . Explained in flowchart 9-11.
 *
 * @param[in] ps_cabac
 *  Pointer to cabac structure
 *
 * @param[in] term_bin
 *  Symbol value, end of slice or not, term_bin is binary
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
void isvce_cabac_encode_terminate(isvce_cabac_ctxt_t *ps_cabac, WORD32 term_bin);

/**
 *******************************************************************************
 * @brief
 * Bypass encoding process for binary decisions:  Explained (9.3.4.4
 *:ITU_T_H264-201402) , flowchart 9-10.
 *
 *  @param[in]  ps_cabac : pointer to cabac context (handle)
 *
 *  @param[in]   bin :  bypass bin(0/1) to be encoded
 *
 *  @returns
 *
 *  @remarks
 *  None
 *
 *******************************************************************************
 */

void isvce_cabac_encode_bypass_bin(isvce_cabac_ctxt_t *ps_cabac, WORD32 bin);

/**
 ******************************************************************************
 *
 *  @brief Encodes a series of bypass bins (FLC bypass bins)
 *
 *  @par   Description
 *  This function is more optimal than calling isvce_cabac_encode_bypass_bin()
 *  in a loop as cabac low, renorm and generating the stream (8bins at a time)
 *  can be done in one operation
 *
 *  @param[inout]ps_cabac
 *   pointer to cabac context (handle)
 *
 *  @param[in]   u4_bins
 *   syntax element to be coded (as FLC bins)
 *
 *  @param[in]   num_bins
 *   This is the FLC length for u4_sym
 *
 *  @return
 *
 ******************************************************************************
 */

void isvce_cabac_encode_bypass_bins(isvce_cabac_ctxt_t *ps_cabac, UWORD32 u4_bins, WORD32 num_bins);

/**
 *******************************************************************************
 *
 * @brief
 *  This function generates CABAC coded bit stream for an Intra Slice.
 *
 * @description
 *  The mb syntax layer for intra slices constitutes luma mb mode, luma sub
 *modes (if present), mb qp delta, coded block pattern, chroma mb mode and
 *  luma/chroma residue. These syntax elements are written as directed by table
 *  7.3.5 of h264 specification.
 *
 * @param[in] ps_ent_ctxt
 *  pointer to entropy context
 *
 * @returns error code
 *
 * @remarks none
 *
 *******************************************************************************
 */
IH264E_ERROR_T isvce_write_islice_mb_cabac(isvce_entropy_ctxt_t *ps_ent_ctxt);

/**
 *******************************************************************************
 *
 * @brief
 *  This function generates CABAC coded bit stream for Inter slices
 *
 * @description
 *  The mb syntax layer for inter slices constitutes luma mb mode, luma sub
 *modes (if present), mb qp delta, coded block pattern, chroma mb mode and
 *  luma/chroma residue. These syntax elements are written as directed by table
 *  7.3.5 of h264 specification
 *
 * @param[in] ps_ent_ctxt
 *  pointer to entropy context
 *
 * @returns error code
 *
 * @remarks none
 *
 *******************************************************************************
 */
IH264E_ERROR_T isvce_write_pslice_mb_cabac(isvce_entropy_ctxt_t *ps_ent_ctxt);

/**
 *******************************************************************************
 *
 * @brief
 *  This function generates CABAC coded bit stream for B slices
 *
 * @description
 *  The mb syntax layer for inter slices constitutes luma mb mode,
 *  mb qp delta, coded block pattern, chroma mb mode and
 *  luma/chroma residue. These syntax elements are written as directed by table
 *  7.3.5 of h264 specification
 *
 * @param[in] ps_ent_ctxt
 *  pointer to entropy context
 *
 * @returns error code
 *
 * @remarks none
 *
 *******************************************************************************
 */
IH264E_ERROR_T isvce_write_bslice_mb_cabac(isvce_entropy_ctxt_t *ps_ent_ctxt);

#if ENABLE_RE_ENC_AS_SKIP
IH264E_ERROR_T isvce_reencode_as_skip_frame_cabac(isvce_entropy_ctxt_t *ps_ent_ctxt);
#endif

#endif
