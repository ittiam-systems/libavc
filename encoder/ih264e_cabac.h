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
*  ih264e_cabac.h
*
* @brief
*  This file contains declarations necessary for cabac encoding
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_CABAC_H_
#define _IH264E_CABAC_H_

/*****************************************************************************/
/* Macros                                                                    */
/*****************************************************************************/

/**
*******************************************************************************
*  @brief Bit precision of cabac engine;
*******************************************************************************
*/
#define CABAC_BITS  9

/**
*******************************************************************************
*  @macro Reverse bits in an unsigned integer
*******************************************************************************
*/
#define REV(u4_input, u4_output)                 \
{                                                \
    UWORD32 u4_temp = (u4_input);                \
    WORD8 i;                                     \
    u4_output = 0;                               \
    for (i = 0; i < 32; i++)                     \
    {                                            \
        u4_output = (u4_output << 1) +           \
                        ((u4_temp >> i) & 0x01); \
    }                                            \
}

/**
*******************************************************************************
*! Bit manipulation macros
*******************************************************************************
*/
#define SETBIT(a, i)   ((a) |= (1 << (i)))
#define CLEARBIT(a, i) ((a) &= ~(1 << (i)))

/**
*******************************************************************************
*! Cabac module expect atlesat MIN_STREAM_SIZE_MB bytes left in stream buffer
*! for encoding an MB
*******************************************************************************
*/
#define MIN_STREAM_SIZE_MB   1024


/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

void ih264e_init_cabac_table(entropy_ctxt_t *ps_ent_ctxt);

void ih264e_init_cabac_ctxt(entropy_ctxt_t *ps_ent_ctxt);

UWORD32 ih264e_cabac_UEGk0_binarization(WORD16 i2_sufs, WORD8 *pi1_bins_len);

void ih264e_get_cabac_context(entropy_ctxt_t *ps_ent_ctxt, WORD32 u4_mb_type);

IH264E_ERROR_T ih264e_cabac_flush(cabac_ctxt_t *ps_cabac_ctxt);

IH264E_ERROR_T ih264e_cabac_put_byte(cabac_ctxt_t *ps_cabac_ctxt);

void ih264e_cabac_encode_bin(cabac_ctxt_t *ps_cabac, WORD32 bin,
                             bin_ctxt_model *pu1_bin_ctxts);

void ih264e_encode_decision_bins(UWORD32 u4_bins, WORD8 i1_bins_len,
                                 UWORD32 u4_ctx_inc, WORD8 i1_valid_len,
                                 bin_ctxt_model *pu1_bin_ctxt_type,
                                 cabac_ctxt_t *ps_cabac);

void ih264e_cabac_encode_terminate(cabac_ctxt_t *ps_cabac, WORD32 term_bin);

void ih264e_cabac_encode_bypass_bin(cabac_ctxt_t *ps_cabac, WORD32 bin);

void ih264e_cabac_encode_bypass_bins(cabac_ctxt_t *ps_cabac, UWORD32 u4_bins,
                                     WORD32 num_bins);

IH264E_ERROR_T ih264e_write_islice_mb_cabac(entropy_ctxt_t *ps_ent_ctxt);

IH264E_ERROR_T ih264e_write_pslice_mb_cabac(entropy_ctxt_t *ps_ent_ctxt);

IH264E_ERROR_T ih264e_write_bslice_mb_cabac(entropy_ctxt_t *ps_ent_ctxt);

#endif /* _IH264E_CABAC_H_ */
