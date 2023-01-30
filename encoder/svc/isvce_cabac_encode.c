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
*  isvce_cabac.c
*
* @brief
*  Contains all functions to encode in CABAC entropy mode
*
*
* @author
* Doney Alex
*
* @par List of Functions:
*
*
* @remarks
*  None
*
*******************************************************************************
*/

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <string.h>

/* User include files */
#include "ih264e_config.h"
#include "ih264_typedefs.h"
#include "iv2.h"
#include "ive2.h"
#include "ih264_debug.h"
#include "isvc_defs.h"
#include "isvce_defs.h"
#include "isvc_macros.h"
#include "ih264e_error.h"
#include "ih264e_bitstream.h"
#include "ime_distortion_metrics.h"
#include "ime_defs.h"
#include "ime_structs.h"
#include "ih264_error.h"
#include "isvc_structs.h"
#include "isvc_trans_quant_itrans_iquant.h"
#include "isvc_inter_pred_filters.h"
#include "isvc_mem_fns.h"
#include "ih264_padding.h"
#include "ih264_platform_macros.h"
#include "ih264_intra_pred_filters.h"
#include "ih264_deblk_edge_filters.h"
#include "isvc_cabac_tables.h"
#include "irc_cntrl_param.h"
#include "irc_frame_info_collector.h"
#include "isvce_rate_control.h"
#include "isvce_cabac_structs.h"
#include "isvce_structs.h"
#include "isvce_cabac.h"
#include "isvce_encode_header.h"
#include "ih264_cavlc_tables.h"
#include "isvce_cavlc.h"
#include "ih264e_statistics.h"
#include "ih264e_trace.h"
#include "isvce_cabac_utils.h"
#include "isvce_utils.h"

/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

/**
 *******************************************************************************
 *
 * @brief
 *  Encodes mb_skip_flag  using CABAC entropy coding mode.
 *
 * @param[in] u1_mb_skip_flag
 *  mb_skip_flag
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 * @param[in] u4_ctxidx_offset
 *  ctxIdxOffset for mb_skip_flag context
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_enc_mb_skip(UWORD8 u1_mb_skip_flag, isvce_cabac_ctxt_t *ps_cabac_ctxt,
                                    UWORD32 u4_ctxidx_offset)
{
    UWORD8 u4_ctx_inc;
    WORD8 a, b;
    a = ((ps_cabac_ctxt->ps_left_ctxt_mb_info->u1_mb_type & CAB_SKIP_MASK) ? 0 : 1);
    b = ((ps_cabac_ctxt->ps_top_ctxt_mb_info->u1_mb_type & CAB_SKIP_MASK) ? 0 : 1);

    u4_ctx_inc = a + b;
    /* Encode the bin */
    isvce_cabac_encode_bin(ps_cabac_ctxt, (UWORD32) u1_mb_skip_flag,
                           ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctxidx_offset + u4_ctx_inc);
}

/* ! < Table 9-36 – Binarization for macroblock types in I slices  in
 * ITU_T_H264-201402 Bits 0-7 : binarised value Bits 8-15: length of binary
 * sequence
 */
static const UWORD32 u4_mb_type_intra[26] = {0x0100, 0x0620, 0x0621, 0x0622, 0x0623, 0x0748, 0x0749,
                                             0x074a, 0x074b, 0x074c, 0x074d, 0x074e, 0x074f, 0x0628,
                                             0x0629, 0x062a, 0x062b, 0x0758, 0x0759, 0x075a, 0x075b,
                                             0x075c, 0x075d, 0x075e, 0x075f, 0x0203};

/* CtxInc for mb types */
static const UWORD32 u4_mb_ctxinc[2][26] = {
    /* Intra CtxInc's */
    {0x00,     0x03467,  0x03467,  0x03467,  0x03467,  0x034567, 0x034567, 0x034567, 0x034567,
     0x034567, 0x034567, 0x034567, 0x034567, 0x03467,  0x03467,  0x03467,  0x03467,  0x034567,
     0x034567, 0x034567, 0x034567, 0x034567, 0x034567, 0x034567, 0x034567, 0x00},
    /* Inter CtxInc's */
    {0x00,      0x001233,  0x001233,  0x001233,  0x001233,  0x0012233, 0x0012233,
     0x0012233, 0x0012233, 0x0012233, 0x0012233, 0x0012233, 0x0012233, 0x001233,
     0x001233,  0x001233,  0x001233,  0x0012233, 0x0012233, 0x0012233, 0x0012233,
     0x0012233, 0x0012233, 0x0012233, 0x0012233, 0x00}};

/**
 *******************************************************************************
 *
 * @brief
 *  Encodes mb_type for an intra MB.
 *
 * @param[in] u4_slice_type
 *  slice type
 *
 * @param[in] u4_intra_mb_type
 *  MB type (Table 7-11)
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 ** @param[in] u4_ctxidx_offset
 *  ctxIdxOffset for mb_type context
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

static void isvce_cabac_enc_intra_mb_type(UWORD32 u4_slice_type, UWORD32 u4_intra_mb_type,
                                          isvce_cabac_ctxt_t *ps_cabac_ctxt,
                                          UWORD32 u4_ctx_idx_offset)
{
    encoding_envirnoment_t *ps_cab_enc_env = &(ps_cabac_ctxt->s_cab_enc_env);
    bin_ctxt_model *pu1_mb_bin_ctxt, *pu1_bin_ctxt;
    UWORD8 u1_bin;
    isvce_mb_info_ctxt_t *ps_left_ctxt = ps_cabac_ctxt->ps_left_ctxt_mb_info;
    isvce_mb_info_ctxt_t *ps_top_ctxt = ps_cabac_ctxt->ps_top_ctxt_mb_info;
    UWORD32 u4_bins;
    UWORD32 u4_ctx_inc;
    WORD8 i1_bins_len;
    UWORD32 u4_code_int_range;
    UWORD32 u4_code_int_low;
    UWORD16 u2_quant_code_int_range;
    UWORD16 u4_code_int_range_lps;
    WORD8 i;
    UWORD8 u1_ctx_inc;
    UWORD32 u4_table_val;

    pu1_mb_bin_ctxt = ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctx_idx_offset;

    u4_bins = u4_mb_type_intra[u4_intra_mb_type];
    i1_bins_len = (WORD8) ((u4_bins >> 8) & 0x0f);
    u4_ctx_inc = u4_mb_ctxinc[(u4_slice_type != ISLICE)][u4_intra_mb_type];
    u1_ctx_inc = 0;
    if(u4_slice_type == ISLICE)
    {
        if(ps_left_ctxt != ps_cabac_ctxt->ps_def_ctxt_mb_info)
            u1_ctx_inc += ((ps_left_ctxt->u1_mb_type != CAB_I4x4) ? 1 : 0);
        if(ps_top_ctxt != ps_cabac_ctxt->ps_def_ctxt_mb_info)
            u1_ctx_inc += ((ps_top_ctxt->u1_mb_type != CAB_I4x4) ? 1 : 0);

        u4_ctx_inc = (u4_ctx_inc | (u1_ctx_inc << ((i1_bins_len - 1) << 2)));
    }
    else
    {
        pu1_mb_bin_ctxt += 3;
        if(u4_slice_type == BSLICE) pu1_mb_bin_ctxt += 2;
    }

    u4_code_int_range = ps_cab_enc_env->u4_code_int_range;
    u4_code_int_low = ps_cab_enc_env->u4_code_int_low;

    for(i = (i1_bins_len - 1); i >= 0; i--)
    {
        WORD32 shift;

        u1_ctx_inc = ((u4_ctx_inc >> (i << 2)) & 0x0f);
        u1_bin = ((u4_bins >> i) & 0x01);
        /* Encode the bin */
        pu1_bin_ctxt = pu1_mb_bin_ctxt + u1_ctx_inc;
        if(i != (i1_bins_len - 2))
        {
            WORD8 i1_mps = !!((*pu1_bin_ctxt) & (0x40));
            WORD8 i1_state = (*pu1_bin_ctxt) & 0x3F;

            u2_quant_code_int_range = ((u4_code_int_range >> 6) & 0x03);
            u4_table_val = gau4_isvc_cabac_table[i1_state][u2_quant_code_int_range];
            u4_code_int_range_lps = u4_table_val & 0xFF;

            u4_code_int_range -= u4_code_int_range_lps;
            if(u1_bin != i1_mps)
            {
                u4_code_int_low += u4_code_int_range;
                u4_code_int_range = u4_code_int_range_lps;
                if(i1_state == 0)
                {
                    /* MPS(CtxIdx) = 1 - MPS(CtxIdx) */
                    i1_mps = 1 - i1_mps;
                }

                i1_state = (u4_table_val >> 15) & 0x3F;
            }
            else
            {
                i1_state = (u4_table_val >> 8) & 0x3F;
            }

            (*pu1_bin_ctxt) = (i1_mps << 6) | i1_state;
        }
        else
        {
            u4_code_int_range -= 2;
        }

        /* Renormalize */
        /*****************************************************************/
        /* Renormalization; calculate bits generated based on range(R)   */
        /* Note : 6 <= R < 512; R is 2 only for terminating encode       */
        /*****************************************************************/
        GETRANGE(shift, u4_code_int_range);
        shift = 9 - shift;
        u4_code_int_low <<= shift;
        u4_code_int_range <<= shift;

        /* bits to be inserted in the bitstream */
        ps_cab_enc_env->u4_bits_gen += shift;
        ps_cab_enc_env->u4_code_int_range = u4_code_int_range;
        ps_cab_enc_env->u4_code_int_low = u4_code_int_low;

        /* generate stream when a byte is ready */
        if(ps_cab_enc_env->u4_bits_gen > CABAC_BITS)
        {
            isvce_cabac_put_byte(ps_cabac_ctxt);
            u4_code_int_range = ps_cab_enc_env->u4_code_int_range;
            u4_code_int_low = ps_cab_enc_env->u4_code_int_low;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief
 *  Encodes prev_intra4x4_pred_mode_flag and
 *  rem_intra4x4_pred_mode using CABAC entropy coding mode
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 *  @param[in] pu1_intra_4x4_modes
 *  Pointer to array containing prev_intra4x4_pred_mode_flag and
 *  rem_intra4x4_pred_mode
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_enc_4x4mb_modes(isvce_cabac_ctxt_t *ps_cabac_ctxt,
                                        UWORD8 *pu1_intra_4x4_modes)
{
    WORD32 i;
    WORD8 byte;
    for(i = 0; i < 16; i += 2)
    {
        /* sub blk idx 1 */
        byte = pu1_intra_4x4_modes[i >> 1];
        if(byte & 0x1)
        {
            isvce_cabac_encode_bin(
                ps_cabac_ctxt, 1,
                ps_cabac_ctxt->au1_cabac_ctxt_table + PREV_INTRA4X4_PRED_MODE_FLAG);
        }
        else
        {
            /* Binarization is FL and Cmax=7 */
            isvce_encode_decision_bins(
                byte & 0xF, 4, 0x05554, 4,
                ps_cabac_ctxt->au1_cabac_ctxt_table + REM_INTRA4X4_PRED_MODE - 5, ps_cabac_ctxt);
        }
        /* sub blk idx 2 */
        byte >>= 4;
        if(byte & 0x1)
        {
            isvce_cabac_encode_bin(
                ps_cabac_ctxt, 1,
                ps_cabac_ctxt->au1_cabac_ctxt_table + PREV_INTRA4X4_PRED_MODE_FLAG);
        }
        else
        {
            isvce_encode_decision_bins(
                byte & 0xF, 4, 0x05554, 4,
                ps_cabac_ctxt->au1_cabac_ctxt_table + REM_INTRA4X4_PRED_MODE - 5, ps_cabac_ctxt);
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief
 *  Encodes chroma  intrapred mode for the MB.
 *
 * @param[in] u1_chroma_pred_mode
 *  Chroma intr prediction mode
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_enc_chroma_predmode(UWORD8 u1_chroma_pred_mode,
                                            isvce_cabac_ctxt_t *ps_cabac_ctxt)
{
    WORD8 i1_temp;
    isvce_mb_info_ctxt_t *ps_curr_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info;
    isvce_mb_info_ctxt_t *ps_left_ctxt = ps_cabac_ctxt->ps_left_ctxt_mb_info;
    isvce_mb_info_ctxt_t *ps_top_ctxt = ps_cabac_ctxt->ps_top_ctxt_mb_info;
    UWORD32 u4_bins = 0;
    WORD8 i1_bins_len = 1;
    UWORD32 u4_ctx_inc = 0;
    UWORD8 a, b;
    a = ((ps_left_ctxt->u1_intrapred_chroma_mode != 0) ? 1 : 0);
    b = ((ps_top_ctxt->u1_intrapred_chroma_mode != 0) ? 1 : 0);

    /* Binarization is TU and Cmax=3 */
    ps_curr_ctxt->u1_intrapred_chroma_mode = u1_chroma_pred_mode;

    u4_ctx_inc = a + b;
    u4_ctx_inc = (u4_ctx_inc | 0x330);
    if(u1_chroma_pred_mode)
    {
        u4_bins = 1;
        i1_temp = u1_chroma_pred_mode;
        i1_temp--;
        /* Put a stream of 1's of length Chromaps_pred_mode_ctxt value */
        while(i1_temp)
        {
            u4_bins = (u4_bins | (1 << i1_bins_len));
            i1_bins_len++;
            i1_temp--;
        }
        /* If Chromaps_pred_mode_ctxt < Cmax i.e 3. Terminate put a zero */
        if(u1_chroma_pred_mode < 3)
        {
            i1_bins_len++;
        }
    }

    isvce_encode_decision_bins(u4_bins, i1_bins_len, u4_ctx_inc, 3,
                               ps_cabac_ctxt->au1_cabac_ctxt_table + INTRA_CHROMA_PRED_MODE,
                               ps_cabac_ctxt);
}

/**
 *******************************************************************************
 *
 * @brief
 *  Encodes CBP for the MB.
 *
 * @param[in] u1_cbp
 *  CBP for the MB
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_enc_cbp(UWORD32 u4_cbp, isvce_cabac_ctxt_t *ps_cabac_ctxt)
{
    isvce_mb_info_ctxt_t *ps_left_ctxt = ps_cabac_ctxt->ps_left_ctxt_mb_info;
    isvce_mb_info_ctxt_t *ps_top_ctxt = ps_cabac_ctxt->ps_top_ctxt_mb_info;
    WORD8 i2_cbp_chroma, i, j;
    UWORD8 u1_ctxt_inc, u1_bin;
    UWORD8 a, b;
    UWORD32 u4_ctx_inc;
    UWORD32 u4_bins;
    WORD8 i1_bins_len;

    /* CBP Luma, FL, Cmax = 15, L = 4 */
    u4_ctx_inc = 0;
    u4_bins = 0;
    i1_bins_len = 5;
    for(i = 0; i < 4; i++)
    {
        /* calulate ctxtInc, depending on neighbour availability */
        /* u1_ctxt_inc = CondTerm(A) + 2 * CondTerm(B);
         A: Left block and B: Top block */

        /* Check for Top availability */
        if(i >> 1)
        {
            j = i - 2;
            /* Top is available always and it's current MB */
            b = (((u4_cbp >> j) & 0x01) != 0 ? 0 : 1);
        }
        else
        {
            /* for blocks whose top reference is in another MB */
            {
                j = i + 2;
                b = ((ps_top_ctxt->u1_cbp >> j) & 0x01) ? 0 : 1;
            }
        }

        /* Check for Left availability */
        if(i & 0x01)
        {
            /* Left is available always and it's current MB */
            j = i - 1;
            a = (((u4_cbp >> j) & 0x01) != 0 ? 0 : 1);
        }
        else
        {
            {
                j = i + 1;
                a = ((ps_left_ctxt->u1_cbp >> j) & 0x01) ? 0 : 1;
            }
        }
        u1_ctxt_inc = a + 2 * b;
        u1_bin = ((u4_cbp >> i) & 0x01);
        u4_ctx_inc = (u4_ctx_inc | (u1_ctxt_inc << (i << 2)));
        u4_bins = (u4_bins | (u1_bin << i));
    }

    /* CBP Chroma, TU, Cmax = 2 */
    i2_cbp_chroma = u4_cbp >> 4;
    /* calulate ctxtInc, depending on neighbour availability */
    a = (ps_left_ctxt->u1_cbp > 15) ? 1 : 0;
    b = (ps_top_ctxt->u1_cbp > 15) ? 1 : 0;

    u1_ctxt_inc = a + 2 * b;
    if(i2_cbp_chroma)
    {
        u4_ctx_inc = u4_ctx_inc | ((4 + u1_ctxt_inc) << 16);
        u4_bins = (u4_bins | 0x10);
        /* calulate ctxtInc, depending on neighbour availability */
        a = (ps_left_ctxt->u1_cbp > 31) ? 1 : 0;
        b = (ps_top_ctxt->u1_cbp > 31) ? 1 : 0;
        u1_ctxt_inc = a + 2 * b;
        u4_ctx_inc = u4_ctx_inc | ((8 + u1_ctxt_inc) << 20);
        u4_bins = (u4_bins | (((i2_cbp_chroma >> 1) & 0x01) << i1_bins_len));
        i1_bins_len++;
    }
    else
    {
        u4_ctx_inc = (u4_ctx_inc | ((4 + u1_ctxt_inc) << 16));
    }
    isvce_encode_decision_bins(u4_bins, i1_bins_len, u4_ctx_inc, 8,
                               ps_cabac_ctxt->au1_cabac_ctxt_table + CBP_LUMA, ps_cabac_ctxt);
}

/**
 *******************************************************************************
 *
 * @brief
 *  Encodes mb_qp_delta for the MB.
 *
 * @param[in] i1_mb_qp_delta
 *  mb_qp_delta
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_enc_mb_qp_delta(WORD8 i1_mb_qp_delta, isvce_cabac_ctxt_t *ps_cabac_ctxt)
{
    UWORD8 u1_code_num;
    UWORD8 u1_ctxt_inc;

    UWORD32 u4_bins;
    WORD8 i1_bins_len;

    /* Range of ps_mb_qp_delta_ctxt= -26 to +25 inclusive */
    ASSERT((i1_mb_qp_delta < 26) && (i1_mb_qp_delta > -27));

    /* if ps_mb_qp_delta_ctxt=0, then codeNum=0 */
    u1_code_num = 0;
    if(i1_mb_qp_delta > 0)
    {
        u1_code_num = (i1_mb_qp_delta << 1) - 1;
    }
    else if(i1_mb_qp_delta < 0)
    {
        u1_code_num = (ABS(i1_mb_qp_delta)) << 1;
    }

    u4_bins = 0;
    i1_bins_len = 1;

    u1_ctxt_inc = !!ps_cabac_ctxt->i1_prevps_mb_qp_delta_ctxt;

    if(u1_code_num == 0)
    {
        isvce_encode_decision_bins(u4_bins, i1_bins_len, u1_ctxt_inc, 3,
                                   ps_cabac_ctxt->au1_cabac_ctxt_table + MB_QP_DELTA,
                                   ps_cabac_ctxt);
    }
    else
    {
        u4_bins = 1;
        u1_code_num--;

        if(u1_code_num == 0)
        {
            i1_bins_len++;

            isvce_encode_decision_bins(u4_bins, i1_bins_len, u1_ctxt_inc | 0x20, 3,
                                       ps_cabac_ctxt->au1_cabac_ctxt_table + MB_QP_DELTA,
                                       ps_cabac_ctxt);
        }
        else
        {
            u4_bins = (u4_bins | (1 << i1_bins_len));
            i1_bins_len++;
            u1_code_num--;

            /* BinIdx from b2 onwards */
            if(u1_code_num < 30)
            {
                /* maximum i1_bins_len = 31 */
                while(u1_code_num)
                {
                    u4_bins = (u4_bins | (1 << i1_bins_len));
                    i1_bins_len++;
                    u1_code_num--;
                };

                i1_bins_len++;

                isvce_encode_decision_bins(u4_bins, i1_bins_len, u1_ctxt_inc | 0x320, 2,
                                           ps_cabac_ctxt->au1_cabac_ctxt_table + MB_QP_DELTA,
                                           ps_cabac_ctxt);
            }
            else
            {
                /* maximum i1_bins_len = 53 */
                u4_bins = 0xffffffff;
                i1_bins_len = 32;
                u1_code_num -= 30;

                isvce_encode_decision_bins(u4_bins, i1_bins_len, u1_ctxt_inc | 0x320, 2,
                                           ps_cabac_ctxt->au1_cabac_ctxt_table + MB_QP_DELTA,
                                           ps_cabac_ctxt);

                u4_bins = 0;
                i1_bins_len = 0;

                while(u1_code_num)
                {
                    u4_bins = (u4_bins | (1 << i1_bins_len));
                    i1_bins_len++;
                    u1_code_num--;
                };

                i1_bins_len++;

                isvce_encode_decision_bins(u4_bins, i1_bins_len, 0x333, 1,
                                           ps_cabac_ctxt->au1_cabac_ctxt_table + MB_QP_DELTA,
                                           ps_cabac_ctxt);
            }
        }
    }
}

/**
 *******************************************************************************
 * @brief
 *  Encodes 4residual_block_cabac as defined in 7.3.5.3.3.
 *
 * @param[in] pi2_res_block
 *  pointer to the array of residues
 *
 * @param[in]  u1_nnz
 *  Number of non zero coeffs in the block
 *
 * @param[in] u1_max_num_coeffs
 *  Max number of coeffs that can be there in the block
 *
 * @param[in] u2_sig_coeff_map
 *  Significant coeff map
 *
 * @param[in] u4_ctx_cat_offset
 *  ctxIdxOffset for  absolute value contexts
 *
 * @param[in]  pu1_ctxt_sig_coeff
 *  Pointer to residual state variables
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_write_coeff4x4(WORD16 *pi2_res_block, UWORD8 u1_nnz,
                                       UWORD8 u1_max_num_coeffs, UWORD16 u2_sig_coeff_map,
                                       UWORD32 u4_ctx_cat_offset,
                                       bin_ctxt_model *pu1_ctxt_sig_coeff,
                                       isvce_cabac_ctxt_t *ps_cabac_ctxt)
{
    WORD8 i;
    WORD16 *pi16_coeffs;
    UWORD32 u4_sig_coeff, u4_bins;
    UWORD32 u4_ctx_inc;
    UWORD8 u1_last_sig_coef_index = (31 - CLZ(u2_sig_coeff_map));

    /* Always put Coded Block Flag as 1 */

    pi16_coeffs = pi2_res_block;
    {
        bin_ctxt_model *pu1_bin_ctxt;
        UWORD8 u1_bin, uc_last;

        i = 0;
        pu1_bin_ctxt = pu1_ctxt_sig_coeff;
        u4_sig_coeff = 0;
        u1_bin = 1;
        if((u1_last_sig_coef_index))
        {
            u1_bin = !!(u2_sig_coeff_map & 01);
        }
        uc_last = 1;

        do
        {
            /* Encode Decision */
            isvce_cabac_encode_bin(ps_cabac_ctxt, u1_bin, pu1_bin_ctxt);

            if(u1_bin & uc_last)
            {
                u4_sig_coeff = (u4_sig_coeff | (1 << i));
                pu1_bin_ctxt = pu1_ctxt_sig_coeff + i + LAST_SIGNIFICANT_COEFF_FLAG_FRAME -
                               SIGNIFICANT_COEFF_FLAG_FRAME;
                u1_bin = (i == u1_last_sig_coef_index);
                uc_last = 0;
            }
            else
            {
                i = i + 1;
                pu1_bin_ctxt = pu1_ctxt_sig_coeff + i;
                u1_bin = (i == u1_last_sig_coef_index);
                uc_last = 1;
                if((i != u1_last_sig_coef_index))
                {
                    u1_bin = !!((u2_sig_coeff_map >> i) & 01);
                }
            }
        } while(!((i > u1_last_sig_coef_index) || (i > (u1_max_num_coeffs - 1))));
    }

    /* Encode coeff_abs_level_minus1 and coeff_sign_flag */
    {
        UWORD8 u1_sign;
        UWORD16 u2_abs_level;
        UWORD8 u1_abs_level_equal1 = 1, u1_abs_level_gt1 = 0;
        UWORD8 u1_ctx_inc;
        UWORD8 u1_coff;
        WORD16 i2_sufs;
        WORD8 i1_bins_len;
        i = u1_last_sig_coef_index;
        pi16_coeffs = pi2_res_block + u1_nnz - 1;
        do
        {
            {
                u4_sig_coeff = u4_sig_coeff & ((1 << i) - 1);
                u4_bins = 0;
                u4_ctx_inc = 0;
                i1_bins_len = 1;
                /* Encode the AbsLevelMinus1 */
                u2_abs_level = ABS(*(pi16_coeffs)) - 1;
                /* CtxInc for bin0 */
                u4_ctx_inc = MIN(u1_abs_level_equal1, 4);
                /* CtxInc for remaining */
                u1_ctx_inc = 5 + MIN(u1_abs_level_gt1, 4);
                u4_ctx_inc = u4_ctx_inc + (u1_ctx_inc << 4);
                if(u2_abs_level)
                {
                    u1_abs_level_gt1++;
                    u1_abs_level_equal1 = 0;
                }
                if(!u1_abs_level_gt1) u1_abs_level_equal1++;

                u1_coff = 14;
                if(u2_abs_level >= u1_coff)
                {
                    /* Prefix TU i.e string of 14 1's */
                    u4_bins = 0x3fff;
                    i1_bins_len = 14;
                    isvce_encode_decision_bins(
                        u4_bins, i1_bins_len, u4_ctx_inc, 1,
                        ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctx_cat_offset, ps_cabac_ctxt);

                    /* Suffix, uses EncodeBypass */
                    i2_sufs = u2_abs_level - u1_coff;

                    u4_bins = isvce_cabac_UEGk0_binarization(i2_sufs, &i1_bins_len);

                    isvce_cabac_encode_bypass_bins(ps_cabac_ctxt, u4_bins, i1_bins_len);
                }
                else
                {
                    /* Prefix only */
                    u4_bins = (1 << u2_abs_level) - 1;
                    i1_bins_len = u2_abs_level + 1;
                    /* Encode Terminating bit */
                    isvce_encode_decision_bins(
                        u4_bins, i1_bins_len, u4_ctx_inc, 1,
                        ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctx_cat_offset, ps_cabac_ctxt);
                }
            }
            /* encode coeff_sign_flag[i] */
            u1_sign = ((*pi16_coeffs) < 0) ? 1 : 0;
            isvce_cabac_encode_bypass_bin(ps_cabac_ctxt, u1_sign);
            i = CLZ(u4_sig_coeff);
            i = 31 - i;
            pi16_coeffs--;
        } while(u4_sig_coeff);
    }
}

/**
 *******************************************************************************
 * @brief
 * Write DC coeffs for intra predicted luma block
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
static void isvce_cabac_encode_residue_luma_dc(isvce_entropy_ctxt_t *ps_ent_ctxt)
{
    /* CABAC context */
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;
    tu_sblk_coeff_data_t *ps_mb_coeff_data;

    /* packed residue */
    void *pv_mb_coeff_data = ps_ent_ctxt->pv_mb_coeff_data;
    UWORD16 u2_sig_coeff_map;
    WORD16 *pi2_res_block;
    UWORD8 u1_nnz;
    UWORD8 u1_cbf;
    isvce_mb_info_ctxt_t *ps_top_ctxt = ps_cabac_ctxt->ps_top_ctxt_mb_info;
    isvce_mb_info_ctxt_t *p_CurCtxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info;

    PARSE_COEFF_DATA_BLOCK_4x4(pv_mb_coeff_data, ps_mb_coeff_data, u1_nnz, u2_sig_coeff_map,
                               pi2_res_block);

    u1_cbf = !!(u1_nnz);

    {
        UWORD32 u4_ctx_inc;
        UWORD8 u1_a, u1_b;

        u1_a = ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] & 0x1;
        u1_b = ps_top_ctxt->u1_yuv_dc_csbp & 0x1;
        u4_ctx_inc = u1_a + (u1_b << 1);

        isvce_cabac_encode_bin(
            ps_cabac_ctxt, u1_cbf,
            ps_cabac_ctxt->au1_cabac_ctxt_table + CBF + (LUMA_DC_CTXCAT << 2) + u4_ctx_inc);
    }

    /* Write coded_block_flag */
    if(u1_cbf)
    {
        isvce_cabac_write_coeff4x4(pi2_res_block, u1_nnz, 15, u2_sig_coeff_map,
                                   COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_0_OFFSET,
                                   ps_cabac_ctxt->au1_cabac_ctxt_table +
                                       SIGNIFICANT_COEFF_FLAG_FRAME + SIG_COEFF_CTXT_CAT_0_OFFSET,
                                   ps_cabac_ctxt);

        ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] |= 0x1;
        p_CurCtxt->u1_yuv_dc_csbp |= 0x1;
    }
    else
    {
        ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
        p_CurCtxt->u1_yuv_dc_csbp &= 0x6;
    }

    ps_ent_ctxt->pv_mb_coeff_data = pv_mb_coeff_data;
}

/**
 *******************************************************************************
 * @brief
 * Write chroma residues to the bitstream
 *
 * @param[in] ps_ent_ctxt
 *  Pointer to entropy context structure
 *
 * @param[in] u1_chroma_cbp
 * coded block pattern, chroma
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_write_chroma_residue(isvce_entropy_ctxt_t *ps_ent_ctxt,
                                             UWORD8 u1_chroma_cbp)
{
    /* CABAC context */
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;
    tu_sblk_coeff_data_t *ps_mb_coeff_data;
    /* packed residue */
    void *pv_mb_coeff_data = ps_ent_ctxt->pv_mb_coeff_data;
    UWORD16 u2_sig_coeff_map;
    UWORD8 u1_nnz;
    isvce_mb_info_ctxt_t *ps_top_ctxt_mb_info, *ps_curr_ctxt;

    ps_top_ctxt_mb_info = ps_cabac_ctxt->ps_top_ctxt_mb_info;
    ps_curr_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info;

    /********************/
    /* Write Chroma DC */
    /********************/
    {
        WORD16 *pi2_res_block;
        UWORD8 u1_left_dc_csbp, u1_top_dc_csbp, u1_uv, u1_cbf;

        u1_left_dc_csbp = (ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0]) >> 1;
        u1_top_dc_csbp = (ps_top_ctxt_mb_info->u1_yuv_dc_csbp) >> 1;

        for(u1_uv = 0; u1_uv < 2; u1_uv++)
        {
            PARSE_COEFF_DATA_BLOCK_4x4(pv_mb_coeff_data, ps_mb_coeff_data, u1_nnz, u2_sig_coeff_map,
                                       pi2_res_block);
            u1_cbf = !!(u1_nnz);
            {
                UWORD8 u1_a, u1_b;
                UWORD32 u4_ctx_inc;
                u1_a = (u1_left_dc_csbp >> u1_uv) & 0x01;
                u1_b = (u1_top_dc_csbp >> u1_uv) & 0x01;
                u4_ctx_inc = (u1_a + (u1_b << 1));

                isvce_cabac_encode_bin(ps_cabac_ctxt, u1_cbf,
                                       ps_cabac_ctxt->au1_cabac_ctxt_table + CBF +
                                           (CHROMA_DC_CTXCAT << 2) + u4_ctx_inc);
            }

            if(u1_cbf)
            {
                isvce_cabac_write_coeff4x4(pi2_res_block, u1_nnz, 3, u2_sig_coeff_map,
                                           COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_3_OFFSET,
                                           ps_cabac_ctxt->au1_cabac_ctxt_table +
                                               SIGNIFICANT_COEFF_FLAG_FRAME +
                                               SIG_COEFF_CTXT_CAT_3_OFFSET,
                                           ps_cabac_ctxt);

                SETBIT(u1_top_dc_csbp, u1_uv);
                SETBIT(u1_left_dc_csbp, u1_uv);
            }
            else
            {
                CLEARBIT(u1_top_dc_csbp, u1_uv);
                CLEARBIT(u1_left_dc_csbp, u1_uv);
            }
        }
        /*************************************************************/
        /*      Update the DC csbp                                   */
        /*************************************************************/
        ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x1;
        ps_curr_ctxt->u1_yuv_dc_csbp &= 0x1;
        ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] |= (u1_left_dc_csbp << 1);
        ps_curr_ctxt->u1_yuv_dc_csbp |= (u1_top_dc_csbp << 1);
    }
    /*******************/
    /* Write Chroma AC */
    /*******************/
    {
        if(u1_chroma_cbp == 2)
        {
            UWORD8 u1_uv_blkno, u1_left_ac_csbp, u1_top_ac_csbp;
            WORD16 *pi2_res_block;
            u1_left_ac_csbp = ps_cabac_ctxt->pu1_left_uv_ac_csbp[0];
            u1_top_ac_csbp = ps_top_ctxt_mb_info->u1_yuv_ac_csbp >> 4;

            for(u1_uv_blkno = 0; u1_uv_blkno < 8; u1_uv_blkno++)
            {
                UWORD8 u1_cbf;
                UWORD8 u1_b2b0, u1_b2b1;
                PARSE_COEFF_DATA_BLOCK_4x4(pv_mb_coeff_data, ps_mb_coeff_data, u1_nnz,
                                           u2_sig_coeff_map, pi2_res_block);

                u1_cbf = !!(u1_nnz);
                u1_b2b0 = ((u1_uv_blkno & 0x4) >> 1) | (u1_uv_blkno & 0x1);
                u1_b2b1 = ((u1_uv_blkno & 0x4) >> 1) | ((u1_uv_blkno & 0x2) >> 1);

                {
                    UWORD8 u1_a, u1_b;
                    UWORD32 u4_ctx_inc;
                    /* write coded_block_flag */
                    u1_a = (u1_left_ac_csbp >> u1_b2b1) & 0x1;
                    u1_b = (u1_top_ac_csbp >> u1_b2b0) & 0x1;
                    u4_ctx_inc = u1_a + (u1_b << 1);

                    isvce_cabac_encode_bin(ps_cabac_ctxt, u1_cbf,
                                           ps_cabac_ctxt->au1_cabac_ctxt_table + CBF +
                                               (CHROMA_AC_CTXCAT << 2) + u4_ctx_inc);
                }
                if(u1_cbf)
                {
                    isvce_cabac_write_coeff4x4(
                        pi2_res_block, u1_nnz, 14, u2_sig_coeff_map,
                        COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_4_OFFSET,
                        ps_cabac_ctxt->au1_cabac_ctxt_table + +SIGNIFICANT_COEFF_FLAG_FRAME +
                            SIG_COEFF_CTXT_CAT_4_OFFSET,
                        ps_cabac_ctxt);

                    SETBIT(u1_left_ac_csbp, u1_b2b1);
                    SETBIT(u1_top_ac_csbp, u1_b2b0);
                }
                else
                {
                    CLEARBIT(u1_left_ac_csbp, u1_b2b1);
                    CLEARBIT(u1_top_ac_csbp, u1_b2b0);
                }
            }
            /*************************************************************/
            /*      Update the AC csbp                                   */
            /*************************************************************/
            ps_cabac_ctxt->pu1_left_uv_ac_csbp[0] = u1_left_ac_csbp;
            ps_curr_ctxt->u1_yuv_ac_csbp &= 0x0f;
            ps_curr_ctxt->u1_yuv_ac_csbp |= (u1_top_ac_csbp << 4);
        }
        else
        {
            ps_cabac_ctxt->pu1_left_uv_ac_csbp[0] = 0;
            ps_curr_ctxt->u1_yuv_ac_csbp &= 0xf;
        }
    }
    ps_ent_ctxt->pv_mb_coeff_data = pv_mb_coeff_data;
}

/**
 *******************************************************************************
 * @brief
 * Encodes Residues for the MB as defined in 7.3.5.3
 *
 * @param[in] ps_ent_ctxt
 *  Pointer to entropy context structure
 *
 * @param[in] u1_cbp
 * coded block pattern
 *
 * @param[in] u1_ctx_cat
 * Context category, LUMA_AC_CTXCAT or LUMA_4x4_CTXCAT
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_encode_residue(isvce_entropy_ctxt_t *ps_ent_ctxt, UWORD32 u4_cbp,
                                       UWORD8 u1_ctx_cat)
{
    /* CABAC context */
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;

    tu_sblk_coeff_data_t *ps_mb_coeff_data;
    /* packed residue */
    void *pv_mb_coeff_data = ps_ent_ctxt->pv_mb_coeff_data;
    UWORD16 u2_sig_coeff_map;
    UWORD8 u1_nnz;
    isvce_mb_info_ctxt_t *ps_curr_ctxt;
    isvce_mb_info_ctxt_t *ps_top_ctxt;
    UWORD8 u1_left_ac_csbp;
    UWORD8 u1_top_ac_csbp;
    UWORD32 u4_ctx_idx_offset_sig_coef, u4_ctx_idx_offset_abs_lvl;
    ps_curr_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info;
    ps_top_ctxt = ps_cabac_ctxt->ps_top_ctxt_mb_info;
    u1_left_ac_csbp = ps_cabac_ctxt->pu1_left_y_ac_csbp[0];
    u1_top_ac_csbp = ps_top_ctxt->u1_yuv_ac_csbp;

    if(u4_cbp & 0xf)
    {
        /*  Write luma residue  */
        UWORD8 u1_offset;
        WORD16 *pi2_res_block;
        UWORD8 u1_subblk_num;
        if(u1_ctx_cat == LUMA_AC_CTXCAT)
        {
            u1_offset = 1;
            u4_ctx_idx_offset_sig_coef = SIG_COEFF_CTXT_CAT_1_OFFSET;
            u4_ctx_idx_offset_abs_lvl = COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_1_OFFSET;
        }
        else
        {
            u1_offset = 0;
            u4_ctx_idx_offset_sig_coef = SIG_COEFF_CTXT_CAT_2_OFFSET;
            u4_ctx_idx_offset_abs_lvl = COEFF_ABS_LEVEL_MINUS1 + COEFF_ABS_LEVEL_CAT_2_OFFSET;
        }

        for(u1_subblk_num = 0; u1_subblk_num < 16; u1_subblk_num++)
        {
            UWORD8 u1_b0, u1_b1, u1_b2, u1_b3, u1_b2b0, u1_b3b1, u1_b3b2;
            u1_b0 = (u1_subblk_num & 0x1);
            u1_b1 = (u1_subblk_num & 0x2) >> 1;
            u1_b2 = (u1_subblk_num & 0x4) >> 2;
            u1_b3 = (u1_subblk_num & 0x8) >> 3;
            u1_b2b0 = (u1_b2 << 1) | (u1_b0);
            u1_b3b1 = (u1_b3 << 1) | (u1_b1);
            u1_b3b2 = (u1_b3 << 1) | (u1_b2);

            if(!((u4_cbp >> u1_b3b2) & 0x1))
            {
                /* ---------------------------------------------------------- */
                /* The current block is not coded so skip all the sub block */
                /* and set the pointer of scan level, csbp accrodingly      */
                /* ---------------------------------------------------------- */
                CLEARBIT(u1_top_ac_csbp, u1_b2b0);
                CLEARBIT(u1_top_ac_csbp, (u1_b2b0 + 1));
                CLEARBIT(u1_left_ac_csbp, u1_b3b1);
                CLEARBIT(u1_left_ac_csbp, (u1_b3b1 + 1));

                u1_subblk_num += 3;
            }
            else
            {
                UWORD8 u1_csbf;

                PARSE_COEFF_DATA_BLOCK_4x4(pv_mb_coeff_data, ps_mb_coeff_data, u1_nnz,
                                           u2_sig_coeff_map, pi2_res_block);

                u1_csbf = !!(u1_nnz);
                {
                    UWORD8 u1_a, u1_b;
                    UWORD32 u4_ctx_inc;
                    u1_b = (u1_top_ac_csbp >> u1_b2b0) & 0x01;
                    u1_a = (u1_left_ac_csbp >> u1_b3b1) & 0x01;
                    u4_ctx_inc = u1_a + (u1_b << 1);

                    /* Encode the bin */
                    isvce_cabac_encode_bin(
                        ps_cabac_ctxt, u1_csbf,
                        ps_cabac_ctxt->au1_cabac_ctxt_table + CBF + (u1_ctx_cat << 2) + u4_ctx_inc);
                }
                /**************************/
                /* Write coded_block_flag */
                /**************************/
                if(u1_csbf)
                {
                    isvce_cabac_write_coeff4x4(pi2_res_block, u1_nnz, (UWORD8) (15 - u1_offset),
                                               u2_sig_coeff_map, u4_ctx_idx_offset_abs_lvl,
                                               ps_cabac_ctxt->au1_cabac_ctxt_table +
                                                   SIGNIFICANT_COEFF_FLAG_FRAME +
                                                   u4_ctx_idx_offset_sig_coef,
                                               ps_cabac_ctxt);

                    SETBIT(u1_top_ac_csbp, u1_b2b0);
                    SETBIT(u1_left_ac_csbp, u1_b3b1);
                }
                else
                {
                    CLEARBIT(u1_top_ac_csbp, u1_b2b0);
                    CLEARBIT(u1_left_ac_csbp, u1_b3b1);
                }
            }
        }
        /**************************************************************************/
        /*                   Update the AC csbp                                   */
        /**************************************************************************/
        ps_cabac_ctxt->pu1_left_y_ac_csbp[0] = u1_left_ac_csbp & 0xf;
        u1_top_ac_csbp &= 0x0f;
        ps_curr_ctxt->u1_yuv_ac_csbp &= 0xf0;
        ps_curr_ctxt->u1_yuv_ac_csbp |= u1_top_ac_csbp;
    }
    else
    {
        ps_cabac_ctxt->pu1_left_y_ac_csbp[0] = 0;
        ps_curr_ctxt->u1_yuv_ac_csbp &= 0xf0;
    }

    /*     Write chroma residue */

    ps_ent_ctxt->pv_mb_coeff_data = pv_mb_coeff_data;
    {
        UWORD8 u1_cbp_chroma;
        u1_cbp_chroma = u4_cbp >> 4;
        if(u1_cbp_chroma)
        {
            isvce_cabac_write_chroma_residue(ps_ent_ctxt, u1_cbp_chroma);
        }
        else
        {
            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x1;
            ps_curr_ctxt->u1_yuv_dc_csbp &= 0x1;
            ps_cabac_ctxt->pu1_left_uv_ac_csbp[0] = 0;
            ps_curr_ctxt->u1_yuv_ac_csbp &= 0xf;
        }
    }
}

/**
 *******************************************************************************
 * @brief
 * Encodes a Motion vector (9.3.3.1.1.7 )
 *
 * @param[in] u1_mvd
 *  Motion vector to be encoded
 *
 * @param[in] u4_ctx_idx_offset
 * *  ctxIdxOffset for MV_X or MV_Ycontext
 *
 * @param[in]  ui2_abs_mvd
 * sum of absolute value of corresponding neighboring motion vectors
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_enc_ctx_mvd(WORD16 u1_mvd, UWORD32 u4_ctx_idx_offset, UWORD16 ui2_abs_mvd,
                                    isvce_cabac_ctxt_t *ps_cabac_ctxt)
{
    UWORD8 u1_bin, u1_ctxt_inc;
    WORD8 k = 3, u1_coff = 9;
    WORD16 i2_abs_mvd, i2_sufs;
    UWORD32 u4_ctx_inc;
    UWORD32 u4_bins;
    WORD8 i1_bins_len;

    /* if mvd < u1_coff
     only Prefix
     else
     Prefix + Suffix

     encode sign bit

     Prefix TU encoding Cmax =u1_coff and Suffix 3rd order Exp-Golomb
     */

    if(ui2_abs_mvd < 3)
        u4_ctx_inc = 0;
    else if(ui2_abs_mvd > 32)
        u4_ctx_inc = 2;
    else
        u4_ctx_inc = 1;

    u4_bins = 0;
    i1_bins_len = 1;

    if(u1_mvd == 0)
    {
        isvce_cabac_encode_bin(
            ps_cabac_ctxt, 0, ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctx_idx_offset + u4_ctx_inc);
    }
    else
    {
        i2_abs_mvd = ABS(u1_mvd);
        if(i2_abs_mvd >= u1_coff)
        {
            /* Prefix TU i.e string of 9 1's */
            u4_bins = 0x1ff;
            i1_bins_len = 9;
            u4_ctx_inc = (u4_ctx_inc | 0x065430);

            isvce_encode_decision_bins(u4_bins, i1_bins_len, u4_ctx_inc, 4,
                                       ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctx_idx_offset,
                                       ps_cabac_ctxt);

            /* Suffix, uses EncodeBypass */
            u4_bins = 0;
            i1_bins_len = 0;
            i2_sufs = i2_abs_mvd - u1_coff;
            while(1)
            {
                if(i2_sufs >= (1 << k))
                {
                    u4_bins = (u4_bins | (1 << (31 - i1_bins_len)));
                    i1_bins_len++;
                    i2_sufs = i2_sufs - (1 << k);
                    k++;
                }
                else
                {
                    i1_bins_len++;
                    while(k--)
                    {
                        u1_bin = ((i2_sufs >> k) & 0x01);
                        u4_bins = (u4_bins | (u1_bin << (31 - i1_bins_len)));
                        i1_bins_len++;
                    }
                    break;
                }
            }
            u4_bins >>= (32 - i1_bins_len);
            isvce_cabac_encode_bypass_bins(ps_cabac_ctxt, u4_bins, i1_bins_len);
        }
        else
        {
            /* Prefix only */
            /* b0 */
            u4_bins = 1;
            i2_abs_mvd--;
            u1_ctxt_inc = 3;
            while(i2_abs_mvd)
            {
                i2_abs_mvd--;
                u4_bins = (u4_bins | (1 << i1_bins_len));
                if(u1_ctxt_inc <= 6)
                {
                    u4_ctx_inc = (u4_ctx_inc | (u1_ctxt_inc << (i1_bins_len << 2)));
                    u1_ctxt_inc++;
                }
                i1_bins_len++;
            }
            /* Encode Terminating bit */
            if(i1_bins_len <= 4) u4_ctx_inc = (u4_ctx_inc | (u1_ctxt_inc << (i1_bins_len << 2)));
            i1_bins_len++;
            isvce_encode_decision_bins(u4_bins, i1_bins_len, u4_ctx_inc, 4,
                                       ps_cabac_ctxt->au1_cabac_ctxt_table + u4_ctx_idx_offset,
                                       ps_cabac_ctxt);
        }
        /* sign bit, uses EncodeBypass */
        if(u1_mvd > 0)
            isvce_cabac_encode_bypass_bin(ps_cabac_ctxt, 0);
        else
            isvce_cabac_encode_bypass_bin(ps_cabac_ctxt, 1);
    }
}

/**
 *******************************************************************************
 * @brief
 * Encodes all motion vectors for a P16x16 MB
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 * @param[in] pi2_mv_ptr
 * Pointer to array of motion vectors
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_enc_mvds_p16x16(isvce_cabac_ctxt_t *ps_cabac_ctxt, WORD16 *pi2_mv_ptr)
{
    UWORD8 u1_abs_mvd_x, u1_abs_mvd_y;
    UWORD8 *pu1_top_mv_ctxt, *pu1_lft_mv_ctxt;
    WORD16 u2_mv;
    u1_abs_mvd_x = 0;
    u1_abs_mvd_y = 0;
    pu1_top_mv_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_mv[0];
    pu1_lft_mv_ctxt = ps_cabac_ctxt->pu1_left_mv_ctxt_inc[0];
    {
        UWORD16 u2_abs_mvd_x_a, u2_abs_mvd_x_b, u2_abs_mvd_y_a, u2_abs_mvd_y_b;
        u2_abs_mvd_x_b = (UWORD16) pu1_top_mv_ctxt[0];
        u2_abs_mvd_y_b = (UWORD16) pu1_top_mv_ctxt[1];
        u2_abs_mvd_x_a = (UWORD16) pu1_lft_mv_ctxt[0];
        u2_abs_mvd_y_a = (UWORD16) pu1_lft_mv_ctxt[1];
        u2_mv = *(pi2_mv_ptr++);

        isvce_cabac_enc_ctx_mvd(u2_mv, MVD_X, (UWORD16) (u2_abs_mvd_x_a + u2_abs_mvd_x_b),
                                ps_cabac_ctxt);

        u1_abs_mvd_x = CLIP3(0, 127, ABS(u2_mv));
        u2_mv = *(pi2_mv_ptr++);

        isvce_cabac_enc_ctx_mvd(u2_mv, MVD_Y, (UWORD16) (u2_abs_mvd_y_a + u2_abs_mvd_y_b),
                                ps_cabac_ctxt);

        u1_abs_mvd_y = CLIP3(0, 127, ABS(u2_mv));
    }
    /***************************************************************/
    /* Store abs_mvd_values cabac contexts                         */
    /***************************************************************/
    pu1_top_mv_ctxt[0] = pu1_lft_mv_ctxt[0] = u1_abs_mvd_x;
    pu1_top_mv_ctxt[1] = pu1_lft_mv_ctxt[1] = u1_abs_mvd_y;
}

/**
 *******************************************************************************
 * @brief
 * Encodes all motion vectors for a B MB (Assues that mbype is B_L0_16x16,
 *B_L1_16x16 or B_Bi_16x16
 *
 * @param[in] ps_cabac_ctxt
 *  Pointer to cabac context structure
 *
 * @param[in] pi2_mv_ptr
 * Pointer to array of motion vectors
 *
 * @returns
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */
static void isvce_cabac_enc_mvds_b16x16(isvce_cabac_ctxt_t *ps_cabac_ctxt, WORD16 *pi2_mv_ptr,
                                        WORD32 i4_mb_part_pred_mode)
{
    /* Encode the differential component of the motion vectors */

    {
        UWORD8 u1_abs_mvd_x, u1_abs_mvd_y;
        UWORD8 *pu1_top_mv_ctxt, *pu1_lft_mv_ctxt;
        WORD16 u2_mv;
        u1_abs_mvd_x = 0;
        u1_abs_mvd_y = 0;
        pu1_top_mv_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_mv[0];
        pu1_lft_mv_ctxt = ps_cabac_ctxt->pu1_left_mv_ctxt_inc[0];
        if(i4_mb_part_pred_mode != L1)
        {
            UWORD16 u2_abs_mvd_x_a, u2_abs_mvd_x_b, u2_abs_mvd_y_a, u2_abs_mvd_y_b;
            u2_abs_mvd_x_b = (UWORD16) pu1_top_mv_ctxt[0];
            u2_abs_mvd_y_b = (UWORD16) pu1_top_mv_ctxt[1];
            u2_abs_mvd_x_a = (UWORD16) pu1_lft_mv_ctxt[0];
            u2_abs_mvd_y_a = (UWORD16) pu1_lft_mv_ctxt[1];
            u2_mv = pi2_mv_ptr[0];

            isvce_cabac_enc_ctx_mvd(u2_mv, MVD_X, (UWORD16) (u2_abs_mvd_x_a + u2_abs_mvd_x_b),
                                    ps_cabac_ctxt);

            u1_abs_mvd_x = CLIP3(0, 127, ABS(u2_mv));
            u2_mv = pi2_mv_ptr[1];

            isvce_cabac_enc_ctx_mvd(u2_mv, MVD_Y, (UWORD16) (u2_abs_mvd_y_a + u2_abs_mvd_y_b),
                                    ps_cabac_ctxt);

            u1_abs_mvd_y = CLIP3(0, 127, ABS(u2_mv));
        }

        /***************************************************************/
        /* Store abs_mvd_values cabac contexts                         */
        /***************************************************************/
        pu1_top_mv_ctxt[0] = pu1_lft_mv_ctxt[0] = u1_abs_mvd_x;
        pu1_top_mv_ctxt[1] = pu1_lft_mv_ctxt[1] = u1_abs_mvd_y;

        u1_abs_mvd_x = 0;
        u1_abs_mvd_y = 0;
        if(i4_mb_part_pred_mode != L0)
        {
            UWORD16 u2_abs_mvd_x_a, u2_abs_mvd_x_b, u2_abs_mvd_y_a, u2_abs_mvd_y_b;
            u2_abs_mvd_x_b = (UWORD16) pu1_top_mv_ctxt[2];
            u2_abs_mvd_y_b = (UWORD16) pu1_top_mv_ctxt[3];
            u2_abs_mvd_x_a = (UWORD16) pu1_lft_mv_ctxt[2];
            u2_abs_mvd_y_a = (UWORD16) pu1_lft_mv_ctxt[3];
            u2_mv = pi2_mv_ptr[2];

            isvce_cabac_enc_ctx_mvd(u2_mv, MVD_X, (UWORD16) (u2_abs_mvd_x_a + u2_abs_mvd_x_b),
                                    ps_cabac_ctxt);

            u1_abs_mvd_x = CLIP3(0, 127, ABS(u2_mv));
            u2_mv = pi2_mv_ptr[3];

            isvce_cabac_enc_ctx_mvd(u2_mv, MVD_Y, (UWORD16) (u2_abs_mvd_y_a + u2_abs_mvd_y_b),
                                    ps_cabac_ctxt);

            u1_abs_mvd_y = CLIP3(0, 127, ABS(u2_mv));
        }
        /***************************************************************/
        /* Store abs_mvd_values cabac contexts                         */
        /***************************************************************/
        pu1_top_mv_ctxt[2] = pu1_lft_mv_ctxt[2] = u1_abs_mvd_x;
        pu1_top_mv_ctxt[3] = pu1_lft_mv_ctxt[3] = u1_abs_mvd_y;
    }
}

static FORCEINLINE void isvce_mb_ctxt_update(isvce_cabac_ctxt_t *ps_cabac_ctxt,
                                             isvce_mb_info_ctxt_t *ps_curr_ctxt,
                                             WORD8 i1_mb_qp_delta, UWORD8 u1_cbp,
                                             UWORD8 u1_base_mode_flag, MBTYPES_T e_mb_type)
{
    UWORD8 u1_is_intra_mb = (e_mb_type == I16x16) || (e_mb_type == I8x8) || (e_mb_type == I4x4);
    UWORD8 u1_is_skip_mb = (e_mb_type == PSKIP) || (e_mb_type == BSKIP);
    UWORD8 u1_is_direct_mb = (e_mb_type == BDIRECT);

    ps_curr_ctxt->u1_cbp = u1_cbp;
    ps_curr_ctxt->u1_base_mode_flag = u1_base_mode_flag;

    if(u1_is_intra_mb || u1_is_skip_mb || u1_is_direct_mb || u1_base_mode_flag)
    {
        memset(ps_curr_ctxt->u1_mv, 0, 16);
        memset(ps_cabac_ctxt->pu1_left_mv_ctxt_inc, 0, 16);
    }

    if((0 == u1_cbp) && (e_mb_type != I16x16))
    {
        ps_curr_ctxt->u1_yuv_ac_csbp = 0;
        ps_curr_ctxt->u1_yuv_dc_csbp = 0;

        ps_cabac_ctxt->pu1_left_uv_ac_csbp[0] = 0;
        ps_cabac_ctxt->pu1_left_y_ac_csbp[0] = 0;
        ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] = 0;
    }

    if(u1_is_skip_mb)
    {
        ps_cabac_ctxt->i1_prevps_mb_qp_delta_ctxt = 0;
    }
    else if((I16x16 != e_mb_type) && (0 == u1_cbp))
    {
        ps_cabac_ctxt->i1_prevps_mb_qp_delta_ctxt = 0;
    }
    else if(0 == i1_mb_qp_delta)
    {
        ps_cabac_ctxt->i1_prevps_mb_qp_delta_ctxt = 0;
    }
    else
    {
        ps_cabac_ctxt->i1_prevps_mb_qp_delta_ctxt = 1;
    }

    if(!u1_is_intra_mb || u1_base_mode_flag)
    {
        ps_curr_ctxt->u1_intrapred_chroma_mode = 0;
    }
}

/**
 *******************************************************************************
 *
 * @brief
 *  This function generates CABAC coded bit stream for an Intra Slice.
 *
 * @description
 *  The mb syntax layer for intra slices constitutes luma mb mode, mb qp delta,
 *coded block pattern, chroma mb mode and luma/chroma residue. These syntax
 *elements are written as directed by table 7.3.5 of h264 specification.
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
IH264E_ERROR_T isvce_write_islice_mb_cabac(isvce_entropy_ctxt_t *ps_ent_ctxt)
{
    isvce_mb_info_ctxt_t *ps_curr_ctxt;

    WORD32 mb_tpm, mb_type, chroma_intra_mode, luma_intra_mode;
    UWORD8 u1_cbp, u1_cbp_l, u1_cbp_c;
    WORD8 mb_qp_delta;
    WORD32 bitstream_start_offset, bitstream_end_offset;
    UWORD8 u1_base_mode_flag;

    bitstrm_t *ps_bitstream = ps_ent_ctxt->ps_bitstrm;
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;
    svc_slice_header_t *ps_svc_slice_header =
        ps_ent_ctxt->ps_svc_slice_hdr_base +
        (ps_ent_ctxt->i4_cur_slice_idx % SVC_MAX_SLICE_HDR_CNT);
    isvce_mb_hdr_common_t *ps_mb_hdr = (isvce_mb_hdr_common_t *) ps_ent_ctxt->pv_mb_header_data;

    UWORD8 *pu1_byte = ps_ent_ctxt->pv_mb_header_data;

    if((ps_bitstream->u4_strm_buf_offset + MIN_STREAM_SIZE_MB) >= ps_bitstream->u4_max_strm_size)
    {
        /* return without corrupting the buffer beyond its size */
        return (IH264E_BITSTREAM_BUFFER_OVERFLOW);
    }

    mb_tpm = ps_mb_hdr->u1_mb_type_mode;
    u1_base_mode_flag = ps_mb_hdr->u1_base_mode_flag;
    u1_cbp = ps_mb_hdr->u1_cbp;
    u1_cbp_c = (u1_cbp >> 4);
    u1_cbp_l = (u1_cbp & 0xF);
    mb_type = mb_tpm & 0xF;

    isvce_get_cabac_context(ps_ent_ctxt, mb_type);
    ps_curr_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info;

    bitstream_start_offset = isvce_get_num_bits(ps_bitstream);

    if(mb_type == I16x16)
    {
        luma_intra_mode = ((mb_tpm >> 4) & 3) + 1 + (u1_cbp_c << 2) + (u1_cbp_l == 15) * 12;
    }
    else
    {
        luma_intra_mode = 0;
    }

    chroma_intra_mode = (mb_tpm >> 6);

    if(ps_ent_ctxt->u1_spatial_layer_id && ps_svc_slice_header->i1_adaptive_base_mode_flag)
    {
        isvce_cabac_enc_base_mode_flag(ps_cabac_ctxt, u1_base_mode_flag);
    }

    if(!u1_base_mode_flag)
    {
        isvce_cabac_enc_intra_mb_type(ISLICE, luma_intra_mode, ps_cabac_ctxt, MB_TYPE_I_SLICE);

        if(mb_type == I4x4)
        {
            isvce_mb_hdr_i4x4_t *ps_mb_hdr_i4x4 =
                (isvce_mb_hdr_i4x4_t *) ps_ent_ctxt->pv_mb_header_data;

            isvce_cabac_enc_4x4mb_modes(ps_cabac_ctxt, ps_mb_hdr_i4x4->au1_sub_blk_modes);
        }

        isvce_cabac_enc_chroma_predmode(chroma_intra_mode, ps_cabac_ctxt);
    }

    if(u1_base_mode_flag || (mb_type != I16x16))
    {
        isvce_cabac_enc_cbp(u1_cbp, ps_cabac_ctxt);
    }

    if((u1_cbp > 0) || (mb_type == I16x16))
    {
        mb_qp_delta =
            ((WORD16) ps_mb_hdr->u1_mb_qp) - ((WORD16) ps_ent_ctxt->ps_mb_qp_ctxt->u1_cur_mb_qp);

        isvce_cabac_enc_mb_qp_delta(mb_qp_delta, ps_cabac_ctxt);
        ps_ent_ctxt->ps_mb_qp_ctxt->u1_cur_mb_qp = ps_mb_hdr->u1_mb_qp;

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_header_bits[0] += bitstream_end_offset - bitstream_start_offset;
        bitstream_start_offset = bitstream_end_offset;

        if(mb_type == I16x16)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I16x16;

            isvce_cabac_encode_residue_luma_dc(ps_ent_ctxt);

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_AC_CTXCAT);

            pu1_byte += sizeof(isvce_mb_hdr_i16x16_t);
        }
        else if(mb_type == I4x4)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I4x4;

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_i4x4_t);
        }
        else if(mb_type == BASE_MODE)
        {
            ps_curr_ctxt->u1_mb_type = CAB_P | CAB_NON_BD16x16;

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_base_mode_t);
        }

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_residue_bits[0] += bitstream_end_offset - bitstream_start_offset;
    }
    else
    {
        mb_qp_delta = 0;

        if(mb_type == I16x16)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I16x16;

            pu1_byte += sizeof(isvce_mb_hdr_i16x16_t);
        }
        else if(mb_type == I4x4)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I4x4;

            pu1_byte += sizeof(isvce_mb_hdr_i4x4_t);
        }
        else if(mb_type == BASE_MODE)
        {
            ps_curr_ctxt->u1_mb_type = CAB_P | CAB_NON_BD16x16;

            pu1_byte += sizeof(isvce_mb_hdr_base_mode_t);
        }

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_header_bits[0] += bitstream_end_offset - bitstream_start_offset;
    }

    isvce_mb_ctxt_update(ps_cabac_ctxt, ps_curr_ctxt, mb_qp_delta, u1_cbp, u1_base_mode_flag,
                         mb_type);

    ps_ent_ctxt->pv_mb_header_data = pu1_byte;

    return IH264E_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief
 *  This function generates CABAC coded bit stream for Inter slices
 *
 * @description
 *  The mb syntax layer for inter slices constitutes luma mb mode, mb qp delta,
 *coded block pattern, chroma mb mode and luma/chroma residue. These syntax
 *elements are written as directed by table 7.3.5 of h264 specification
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
IH264E_ERROR_T isvce_write_pslice_mb_cabac(isvce_entropy_ctxt_t *ps_ent_ctxt)
{
    isvce_mb_info_ctxt_t *ps_curr_ctxt;

    WORD32 mb_tpm, mb_type, chroma_intra_mode, luma_intra_mode;
    UWORD8 u1_cbp, u1_cbp_l, u1_cbp_c;
    WORD8 mb_qp_delta;
    WORD32 bitstream_start_offset, bitstream_end_offset;
    UWORD8 u1_base_mode_flag;
    UWORD8 u1_is_intra_mb;

    bitstrm_t *ps_bitstream = ps_ent_ctxt->ps_bitstrm;
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;
    svc_slice_header_t *ps_svc_slice_header =
        ps_ent_ctxt->ps_svc_slice_hdr_base +
        (ps_ent_ctxt->i4_cur_slice_idx % SVC_MAX_SLICE_HDR_CNT);
    isvce_mb_hdr_common_t *ps_mb_hdr = (isvce_mb_hdr_common_t *) ps_ent_ctxt->pv_mb_header_data;

    UWORD8 *pu1_byte = ps_ent_ctxt->pv_mb_header_data;

    if((ps_bitstream->u4_strm_buf_offset + MIN_STREAM_SIZE_MB) >= ps_bitstream->u4_max_strm_size)
    {
        /* return without corrupting the buffer beyond its size */
        return IH264E_BITSTREAM_BUFFER_OVERFLOW;
    }

    /* mb header info */
    mb_tpm = ps_mb_hdr->u1_mb_type_mode;
    u1_base_mode_flag = ps_mb_hdr->u1_base_mode_flag;
    u1_cbp = ps_mb_hdr->u1_cbp;
    u1_cbp_c = (u1_cbp >> 4);
    u1_cbp_l = (u1_cbp & 0xF);

    /* mb type */
    mb_type = mb_tpm & 0xF;
    u1_is_intra_mb = (mb_type == I16x16) || (mb_type == I8x8) || (mb_type == I4x4);

    /* CABAC contexts for the MB */
    isvce_get_cabac_context(ps_ent_ctxt, mb_type);
    ps_curr_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info;

    /* Starting bitstream offset for header in bits */
    bitstream_start_offset = isvce_get_num_bits(ps_bitstream);

    /* Encode mb_skip_flag */
    isvce_cabac_enc_mb_skip(mb_type == PSKIP, ps_cabac_ctxt, MB_SKIP_FLAG_P_SLICE);

    if(mb_type == PSKIP)
    {
        ps_curr_ctxt->u1_mb_type = CAB_P_SKIP;

        ps_ent_ctxt->pi4_mb_skip_run[0]++;

        isvce_mb_ctxt_update(ps_cabac_ctxt, ps_curr_ctxt, 0, 0, 0, PSKIP);

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_header_bits[!u1_is_intra_mb] +=
            bitstream_end_offset - bitstream_start_offset;

        pu1_byte += sizeof(isvce_mb_hdr_pskip_t);
        ps_ent_ctxt->pv_mb_header_data = pu1_byte;

        return IH264E_SUCCESS;
    }

    if(ps_ent_ctxt->u1_spatial_layer_id && ps_svc_slice_header->i1_adaptive_base_mode_flag)
    {
        isvce_cabac_enc_base_mode_flag(ps_cabac_ctxt, u1_base_mode_flag);
    }

    if(!u1_base_mode_flag)
    {
        if(u1_is_intra_mb)
        {
            if(mb_type == I16x16)
            {
                luma_intra_mode = ((mb_tpm >> 4) & 3) + 1 + (u1_cbp_c << 2) + (u1_cbp_l == 15) * 12;
            }
            else
            {
                luma_intra_mode = 0;
            }

            isvce_cabac_encode_bin(ps_cabac_ctxt, 1,
                                   ps_cabac_ctxt->au1_cabac_ctxt_table + MB_TYPE_P_SLICE);

            isvce_cabac_enc_intra_mb_type(PSLICE, (UWORD8) luma_intra_mode, ps_cabac_ctxt,
                                          MB_TYPE_P_SLICE);

            if(mb_type == I4x4)
            {
                isvce_mb_hdr_i4x4_t *ps_mb_hdr_i4x4 =
                    (isvce_mb_hdr_i4x4_t *) ps_ent_ctxt->pv_mb_header_data;

                isvce_cabac_enc_4x4mb_modes(ps_cabac_ctxt, ps_mb_hdr_i4x4->au1_sub_blk_modes);
            }

            chroma_intra_mode = (mb_tpm >> 6);

            isvce_cabac_enc_chroma_predmode(chroma_intra_mode, ps_cabac_ctxt);
        }
        else
        {
            UWORD32 u4_ctx_inc_p;

            isvce_mb_hdr_p16x16_t *ps_mb_hdr_p16x16 =
                (isvce_mb_hdr_p16x16_t *) ps_ent_ctxt->pv_mb_header_data;

            WORD16 *pi2_mv_ptr = (WORD16 *) ps_mb_hdr_p16x16->ai2_mvd;

            /* Encoding mb_type as P16x16 */
            u4_ctx_inc_p = (0x010 + ((2) << 8));

            isvce_encode_decision_bins(0, 3, u4_ctx_inc_p, 3,
                                       &(ps_cabac_ctxt->au1_cabac_ctxt_table[MB_TYPE_P_SLICE]),
                                       ps_cabac_ctxt);

            if(ps_ent_ctxt->u1_spatial_layer_id &&
               ps_svc_slice_header->i1_adaptive_motion_prediction_flag)
            {
                isvce_cabac_enc_motion_prediction_flag(ps_cabac_ctxt, ps_mb_hdr_p16x16->u1_mvp_idx,
                                                       1);
            }

            isvce_cabac_enc_mvds_p16x16(ps_cabac_ctxt, pi2_mv_ptr);
        }
    }

    if(ps_ent_ctxt->u1_spatial_layer_id && (u1_base_mode_flag || !u1_is_intra_mb) &&
       ps_svc_slice_header->i1_adaptive_residual_prediction_flag)
    {
        isvce_cabac_enc_residual_prediction_flag(ps_cabac_ctxt, u1_base_mode_flag,
                                                 ps_mb_hdr->u1_residual_prediction_flag);
    }

    if(u1_base_mode_flag || (mb_type != I16x16))
    {
        isvce_cabac_enc_cbp(u1_cbp, ps_cabac_ctxt);
    }

    if((u1_cbp > 0) || (mb_type == I16x16))
    {
        mb_qp_delta =
            ((WORD16) ps_mb_hdr->u1_mb_qp) - ((WORD16) ps_ent_ctxt->ps_mb_qp_ctxt->u1_cur_mb_qp);

        isvce_cabac_enc_mb_qp_delta(mb_qp_delta, ps_cabac_ctxt);
        ps_ent_ctxt->ps_mb_qp_ctxt->u1_cur_mb_qp = ps_mb_hdr->u1_mb_qp;

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_header_bits[!u1_is_intra_mb] +=
            bitstream_end_offset - bitstream_start_offset;

        bitstream_start_offset = bitstream_end_offset;

        if(mb_type == I16x16)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I16x16;

            isvce_cabac_encode_residue_luma_dc(ps_ent_ctxt);

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_AC_CTXCAT);

            pu1_byte += sizeof(isvce_mb_hdr_i16x16_t);
        }
        else if(mb_type == I4x4)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I4x4;

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_i4x4_t);
        }
        else if(mb_type == P16x16)
        {
            ps_curr_ctxt->u1_mb_type = (CAB_P | CAB_NON_BD16x16);

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_p16x16_t);
        }
        else if(mb_type == BASE_MODE)
        {
            ps_curr_ctxt->u1_mb_type = (CAB_P | CAB_NON_BD16x16);

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_base_mode_t);
        }

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_residue_bits[!u1_is_intra_mb] +=
            bitstream_end_offset - bitstream_start_offset;
    }
    else
    {
        mb_qp_delta = 0;

        if(mb_type == I16x16)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I16x16;

            pu1_byte += sizeof(isvce_mb_hdr_i16x16_t);
        }
        else if(mb_type == I4x4)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I4x4;

            pu1_byte += sizeof(isvce_mb_hdr_i4x4_t);
        }
        else if(mb_type == P16x16)
        {
            ps_curr_ctxt->u1_mb_type = (CAB_P | CAB_NON_BD16x16);

            pu1_byte += sizeof(isvce_mb_hdr_p16x16_t);
        }
        else if(mb_type == BASE_MODE)
        {
            ps_curr_ctxt->u1_mb_type = (CAB_P | CAB_NON_BD16x16);

            pu1_byte += sizeof(isvce_mb_hdr_base_mode_t);
        }

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_header_bits[!u1_is_intra_mb] +=
            bitstream_end_offset - bitstream_start_offset;
    }

    isvce_mb_ctxt_update(ps_cabac_ctxt, ps_curr_ctxt, mb_qp_delta, u1_cbp, u1_base_mode_flag,
                         mb_type);

    ps_ent_ctxt->pv_mb_header_data = pu1_byte;

    return IH264E_SUCCESS;
}

/* ! < Table 9-37 – Binarization for macroblock types in B slices  in
 * ITU_T_H264-201402 Bits 0-7 : binarised value Bits 8-15: length of binary
 * sequence */

static const UWORD32 u4_b_mb_type[27] = {
    0x0100, 0x0301, 0x0305, 0x0603, 0x0623, 0x0613, 0x0633, 0x060b, 0x062b, 0x061b, 0x063b, 0x061f,
    0x0707, 0x0747, 0x0727, 0x0767, 0x0717, 0x0757, 0x0737, 0x0777, 0x070f, 0x074f, 0x063f};
/* CtxInc for mb types in B slices */
static const UWORD32 ui_b_mb_type_ctx_inc[27] = {
    0x00,       0x0530,     0x0530,     0x0555430,  0x0555430,  0x0555430,  0x0555430,  0x0555430,
    0x0555430,  0x0555430,  0x0555430,  0x0555430,  0x05555430, 0x05555430, 0x05555430, 0x05555430,
    0x05555430, 0x05555430, 0x05555430, 0x05555430, 0x05555430, 0x05555430, 0x0555430};

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
IH264E_ERROR_T isvce_write_bslice_mb_cabac(isvce_entropy_ctxt_t *ps_ent_ctxt)
{
    isvce_mb_info_ctxt_t *ps_curr_ctxt;

    WORD32 mb_tpm, mb_type, chroma_intra_mode, luma_intra_mode;
    UWORD8 u1_cbp, u1_cbp_l, u1_cbp_c;
    WORD8 mb_qp_delta;
    WORD32 bitstream_start_offset, bitstream_end_offset;
    UWORD8 u1_base_mode_flag;
    UWORD8 u1_is_intra_mb;

    bitstrm_t *ps_bitstream = ps_ent_ctxt->ps_bitstrm;
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;
    svc_slice_header_t *ps_svc_slice_header =
        ps_ent_ctxt->ps_svc_slice_hdr_base +
        (ps_ent_ctxt->i4_cur_slice_idx % SVC_MAX_SLICE_HDR_CNT);
    isvce_mb_hdr_common_t *ps_mb_hdr = (isvce_mb_hdr_common_t *) ps_ent_ctxt->pv_mb_header_data;

    UWORD8 *pu1_byte = ps_ent_ctxt->pv_mb_header_data;

    if((ps_bitstream->u4_strm_buf_offset + MIN_STREAM_SIZE_MB) >= ps_bitstream->u4_max_strm_size)
    {
        /* return without corrupting the buffer beyond its size */
        return (IH264E_BITSTREAM_BUFFER_OVERFLOW);
    }

    /* mb header info */
    mb_tpm = ps_mb_hdr->u1_mb_type_mode;
    u1_base_mode_flag = ps_mb_hdr->u1_base_mode_flag;
    u1_cbp = ps_mb_hdr->u1_cbp;
    u1_cbp_c = (u1_cbp >> 4);
    u1_cbp_l = (u1_cbp & 0xF);

    /* mb type */
    mb_type = mb_tpm & 0xF;
    u1_is_intra_mb = (mb_type == I16x16) || (mb_type == I8x8) || (mb_type == I4x4);

    /* CABAC contexts for the MB */
    isvce_get_cabac_context(ps_ent_ctxt, mb_type);
    ps_curr_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info;

    /* Starting bitstream offset for header in bits */
    bitstream_start_offset = isvce_get_num_bits(ps_bitstream);

    /* Encode mb_skip_flag */
    isvce_cabac_enc_mb_skip(mb_type == BSKIP, ps_cabac_ctxt, MB_SKIP_FLAG_B_SLICE);

    if(mb_type == BSKIP)
    {
        ps_curr_ctxt->u1_mb_type = CAB_B_SKIP;

        ps_ent_ctxt->pi4_mb_skip_run[0]++;

        isvce_mb_ctxt_update(ps_cabac_ctxt, ps_curr_ctxt, 0, 0, 0, BSKIP);

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_header_bits[!u1_is_intra_mb] +=
            bitstream_end_offset - bitstream_start_offset;

        pu1_byte += sizeof(isvce_mb_hdr_bskip_t);
        ps_ent_ctxt->pv_mb_header_data = pu1_byte;

        return IH264E_SUCCESS;
    }

    if(ps_ent_ctxt->u1_spatial_layer_id && ps_svc_slice_header->i1_adaptive_base_mode_flag)
    {
        isvce_cabac_enc_base_mode_flag(ps_cabac_ctxt, u1_base_mode_flag);
    }

    if(!u1_base_mode_flag)
    {
        if(u1_is_intra_mb)
        {
            if(mb_type == I16x16)
            {
                luma_intra_mode = ((mb_tpm >> 4) & 3) + 1 + (u1_cbp_c << 2) + (u1_cbp_l == 15) * 12;
            }
            else
            {
                luma_intra_mode = 0;
            }

            {
                isvce_mb_info_ctxt_t *ps_left_ctxt = ps_cabac_ctxt->ps_left_ctxt_mb_info;
                isvce_mb_info_ctxt_t *ps_top_ctxt = ps_cabac_ctxt->ps_top_ctxt_mb_info;

                UWORD32 u4_ctx_inc = 0;

                if(ps_left_ctxt != ps_cabac_ctxt->ps_def_ctxt_mb_info)
                {
                    u4_ctx_inc +=
                        ((ps_left_ctxt->u1_mb_type & CAB_BD16x16_MASK) != CAB_BD16x16) ? 1 : 0;
                }

                if(ps_top_ctxt != ps_cabac_ctxt->ps_def_ctxt_mb_info)
                {
                    u4_ctx_inc +=
                        ((ps_top_ctxt->u1_mb_type & CAB_BD16x16_MASK) != CAB_BD16x16) ? 1 : 0;
                }

                /* Intra Prefix Only "111101" */
                u4_ctx_inc = (u4_ctx_inc | 0x05555430);
                isvce_encode_decision_bins(0x2f, 6, u4_ctx_inc, 3,
                                           ps_cabac_ctxt->au1_cabac_ctxt_table + MB_TYPE_B_SLICE,
                                           ps_cabac_ctxt);

                isvce_cabac_enc_intra_mb_type(BSLICE, (UWORD8) luma_intra_mode, ps_cabac_ctxt,
                                              MB_TYPE_B_SLICE);
            }

            if(mb_type == I4x4)
            {
                isvce_mb_hdr_i4x4_t *ps_mb_hdr_i4x4 =
                    (isvce_mb_hdr_i4x4_t *) ps_ent_ctxt->pv_mb_header_data;

                isvce_cabac_enc_4x4mb_modes(ps_cabac_ctxt, ps_mb_hdr_i4x4->au1_sub_blk_modes);
            }

            chroma_intra_mode = (mb_tpm >> 6);

            isvce_cabac_enc_chroma_predmode(chroma_intra_mode, ps_cabac_ctxt);
        }
        else if(mb_type == BDIRECT)
        {
            /* Encoding mb_type as B_Direct_16x16 */
            {
                isvce_mb_info_ctxt_t *ps_left_ctxt = ps_cabac_ctxt->ps_left_ctxt_mb_info;
                isvce_mb_info_ctxt_t *ps_top_ctxt = ps_cabac_ctxt->ps_top_ctxt_mb_info;

                UWORD32 u4_ctx_inc = 0;

                if(ps_left_ctxt != ps_cabac_ctxt->ps_def_ctxt_mb_info)
                {
                    u4_ctx_inc +=
                        ((ps_left_ctxt->u1_mb_type & CAB_BD16x16_MASK) != CAB_BD16x16) ? 1 : 0;
                }

                if(ps_top_ctxt != ps_cabac_ctxt->ps_def_ctxt_mb_info)
                {
                    u4_ctx_inc +=
                        ((ps_top_ctxt->u1_mb_type & CAB_BD16x16_MASK) != CAB_BD16x16) ? 1 : 0;
                }

                /* Encode the bin */
                isvce_cabac_encode_bin(
                    ps_cabac_ctxt, 0,
                    ps_cabac_ctxt->au1_cabac_ctxt_table + MB_TYPE_B_SLICE + u4_ctx_inc);
            }
        }
        else
        {
            WORD32 i;

            isvce_mb_hdr_b16x16_t *ps_mb_hdr_b16x16 =
                (isvce_mb_hdr_b16x16_t *) ps_ent_ctxt->pv_mb_header_data;

            WORD16 *pi2_mv_ptr = (WORD16 *) ps_mb_hdr_b16x16->ai2_mvd;
            WORD32 i4_mb_part_pred_mode = (mb_tpm >> 4);
            UWORD32 u4_mb_type = mb_type - B16x16 + B_L0_16x16 + i4_mb_part_pred_mode;

            /* Encoding mb_type as B16x16 */
            {
                isvce_mb_info_ctxt_t *ps_left_ctxt = ps_cabac_ctxt->ps_left_ctxt_mb_info;
                isvce_mb_info_ctxt_t *ps_top_ctxt = ps_cabac_ctxt->ps_top_ctxt_mb_info;
                UWORD32 u4_ctx_inc = 0;

                UWORD32 u4_mb_type_bins = u4_b_mb_type[u4_mb_type];
                UWORD32 u4_bin_len = (u4_mb_type_bins >> 8) & 0x0F;
                u4_mb_type_bins = u4_mb_type_bins & 0xFF;

                if(ps_left_ctxt != ps_cabac_ctxt->ps_def_ctxt_mb_info)
                    u4_ctx_inc +=
                        ((ps_left_ctxt->u1_mb_type & CAB_BD16x16_MASK) != CAB_BD16x16) ? 1 : 0;
                if(ps_top_ctxt != ps_cabac_ctxt->ps_def_ctxt_mb_info)
                    u4_ctx_inc +=
                        ((ps_top_ctxt->u1_mb_type & CAB_BD16x16_MASK) != CAB_BD16x16) ? 1 : 0;

                u4_ctx_inc = u4_ctx_inc | ui_b_mb_type_ctx_inc[u4_mb_type];

                isvce_encode_decision_bins(u4_mb_type_bins, u4_bin_len, u4_ctx_inc, u4_bin_len,
                                           &(ps_cabac_ctxt->au1_cabac_ctxt_table[MB_TYPE_B_SLICE]),
                                           ps_cabac_ctxt);
            }

            for(i = 0; i < NUM_PRED_DIRS; i++)
            {
                PRED_MODE_T e_pred_mode = (PRED_MODE_T) i;
                PRED_MODE_T e_cmpl_pred_mode = (e_pred_mode == L0) ? L1 : L0;

                if(((PRED_MODE_T) i4_mb_part_pred_mode) != e_pred_mode)
                {
                    if(ps_svc_slice_header->i1_adaptive_motion_prediction_flag &&
                       ps_ent_ctxt->u1_spatial_layer_id)
                    {
                        isvce_cabac_enc_motion_prediction_flag(
                            ps_cabac_ctxt, ps_mb_hdr_b16x16->au1_mvp_idx[e_cmpl_pred_mode],
                            e_cmpl_pred_mode == L0);
                    }
                }
            }

            isvce_cabac_enc_mvds_b16x16(ps_cabac_ctxt, pi2_mv_ptr, i4_mb_part_pred_mode);
        }
    }

    if(ps_svc_slice_header->i1_adaptive_residual_prediction_flag &&
       ps_ent_ctxt->u1_spatial_layer_id && (u1_base_mode_flag || !u1_is_intra_mb))
    {
        isvce_cabac_enc_residual_prediction_flag(ps_cabac_ctxt, u1_base_mode_flag,
                                                 ps_mb_hdr->u1_residual_prediction_flag);
    }

    if(u1_base_mode_flag || (mb_type != I16x16))
    {
        isvce_cabac_enc_cbp(u1_cbp, ps_cabac_ctxt);
    }

    if((u1_cbp > 0) || (mb_type == I16x16))
    {
        mb_qp_delta =
            ((WORD16) ps_mb_hdr->u1_mb_qp) - ((WORD16) ps_ent_ctxt->ps_mb_qp_ctxt->u1_cur_mb_qp);

        isvce_cabac_enc_mb_qp_delta(mb_qp_delta, ps_cabac_ctxt);
        ps_ent_ctxt->ps_mb_qp_ctxt->u1_cur_mb_qp = ps_mb_hdr->u1_mb_qp;

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_header_bits[!u1_is_intra_mb] +=
            bitstream_end_offset - bitstream_start_offset;
        bitstream_start_offset = bitstream_end_offset;

        if(mb_type == I16x16)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I16x16;

            isvce_cabac_encode_residue_luma_dc(ps_ent_ctxt);

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_AC_CTXCAT);

            pu1_byte += sizeof(isvce_mb_hdr_i16x16_t);
        }
        else if(mb_type == I4x4)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I4x4;

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_i4x4_t);
        }
        else if(mb_type == B16x16)
        {
            ps_curr_ctxt->u1_mb_type = CAB_NON_BD16x16;

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_b16x16_t);
        }
        else if(mb_type == BDIRECT)
        {
            ps_curr_ctxt->u1_mb_type = CAB_BD16x16;

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_b16x16_t);
        }
        else if(mb_type == BASE_MODE)
        {
            ps_curr_ctxt->u1_mb_type = CAB_NON_BD16x16;

            isvce_cabac_encode_residue(ps_ent_ctxt, u1_cbp, LUMA_4X4_CTXCAT);

            ps_cabac_ctxt->pu1_left_yuv_dc_csbp[0] &= 0x6;
            ps_cabac_ctxt->ps_curr_ctxt_mb_info->u1_yuv_dc_csbp &= 0x6;

            pu1_byte += sizeof(isvce_mb_hdr_base_mode_t);
        }

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_residue_bits[!u1_is_intra_mb] +=
            bitstream_end_offset - bitstream_start_offset;
    }
    else
    {
        mb_qp_delta = 0;

        if(mb_type == I16x16)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I16x16;

            pu1_byte += sizeof(isvce_mb_hdr_i16x16_t);
        }
        else if(mb_type == I4x4)
        {
            ps_curr_ctxt->u1_mb_type = CAB_I4x4;

            pu1_byte += sizeof(isvce_mb_hdr_i4x4_t);
        }
        else if(mb_type == B16x16)
        {
            ps_curr_ctxt->u1_mb_type = CAB_NON_BD16x16;

            pu1_byte += sizeof(isvce_mb_hdr_b16x16_t);
        }
        else if(mb_type == BDIRECT)
        {
            ps_curr_ctxt->u1_mb_type = CAB_BD16x16;

            pu1_byte += sizeof(isvce_mb_hdr_b16x16_t);
        }
        else if(mb_type == BDIRECT)
        {
            ps_curr_ctxt->u1_mb_type = CAB_NON_BD16x16;

            pu1_byte += sizeof(isvce_mb_hdr_base_mode_t);
        }

        bitstream_end_offset = isvce_get_num_bits(ps_bitstream);
        ps_ent_ctxt->u4_header_bits[!u1_is_intra_mb] +=
            bitstream_end_offset - bitstream_start_offset;
    }

    isvce_mb_ctxt_update(ps_cabac_ctxt, ps_curr_ctxt, mb_qp_delta, u1_cbp, u1_base_mode_flag,
                         mb_type);

    ps_ent_ctxt->pv_mb_header_data = pu1_byte;

    return IH264E_SUCCESS;
}

#if ENABLE_RE_ENC_AS_SKIP
IH264E_ERROR_T isvce_reencode_as_skip_frame_cabac(isvce_entropy_ctxt_t *ps_ent_ctxt)
{
    isvce_cabac_ctxt_t *ps_cabac_ctxt = ps_ent_ctxt->ps_cabac;
    bitstrm_t *ps_bitstrm = ps_ent_ctxt->ps_bitstrm;
    bitstrm_t *ps_bitstrm_after_slice_hdr = ps_ent_ctxt->ps_bitstrm_after_slice_hdr;

    isvce_mb_info_ctxt_t *ps_curr_ctxt;

    slice_header_t *ps_slice_header =
        (ps_ent_ctxt->u1_spatial_layer_id == 0)
            ? &ps_ent_ctxt->ps_slice_hdr_base[ps_ent_ctxt->i4_cur_slice_idx % SVC_MAX_SLICE_HDR_CNT]
            : &ps_ent_ctxt
                   ->ps_svc_slice_hdr_base[ps_ent_ctxt->i4_cur_slice_idx % SVC_MAX_SLICE_HDR_CNT]
                   .s_slice_header;

    /* total mb cnt */
    UWORD32 i4_wd_mbs = ps_ent_ctxt->i4_wd_mbs;
    UWORD32 i4_ht_mbs = ps_ent_ctxt->i4_ht_mbs;
    UWORD8 i, j;

    isvce_init_cabac_ctxt(ps_ent_ctxt, ps_slice_header);

    ps_bitstrm->i4_bits_left_in_cw = ps_bitstrm_after_slice_hdr->i4_bits_left_in_cw;
    ps_bitstrm->u4_cur_word = ps_bitstrm_after_slice_hdr->u4_cur_word;
    ps_bitstrm->u4_strm_buf_offset = ps_bitstrm_after_slice_hdr->u4_strm_buf_offset;
    ps_bitstrm->i4_zero_bytes_run = ps_bitstrm_after_slice_hdr->i4_zero_bytes_run;

    for(i = 0; i < i4_ht_mbs; i++)
    {
        for(j = 0; j < i4_wd_mbs; j++)
        {
            MBTYPES_T mb_type = PSKIP;

            ps_ent_ctxt->i4_mb_x = j;
            ps_ent_ctxt->i4_mb_y = i;

            isvce_get_cabac_context(ps_ent_ctxt, mb_type);
            ps_curr_ctxt = ps_cabac_ctxt->ps_curr_ctxt_mb_info;

            isvce_cabac_enc_mb_skip(mb_type == PSKIP, ps_cabac_ctxt, MB_SKIP_FLAG_P_SLICE);

            if(mb_type == PSKIP)
            {
                ps_curr_ctxt->u1_mb_type = CAB_P_SKIP;
                isvce_mb_ctxt_update(ps_cabac_ctxt, ps_curr_ctxt, 0, 0, 0, PSKIP);
            }

            if(j == i4_wd_mbs - 1 && i == i4_ht_mbs - 1)
            {
                isvce_cabac_encode_terminate(ps_cabac_ctxt, 1);
            }
            else
            {
                isvce_cabac_encode_terminate(ps_cabac_ctxt, 0);
            }
        }
    }
    return IH264E_SUCCESS;
}
#endif
