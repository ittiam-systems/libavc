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
* @file
*  ih264e_encode_header.h
*
* @brief
*  This file contains structures and interface prototypes for h264 bitstream
*  header encoding
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_ENCODE_HEADER_H_
#define _IH264E_ENCODE_HEADER_H_

/*****************************************************************************/
/* Function Macros                                                           */
/*****************************************************************************/

/**
******************************************************************************
 *  @brief   Macro to put a code with specified number of bits into the
 *           bitstream
******************************************************************************
 */
#define PUT_BITS(ps_bitstrm, code_val, code_len, ret_val, syntax_string)     \
        {                                                                    \
            ENTROPY_TRACE(syntax_string, code_val);                          \
            ret_val = ih264e_put_bits((ps_bitstrm), (code_val), (code_len)); \
            if(ret_val != IH264E_SUCCESS)                                    \
            {                                                                \
                return ret_val;                                              \
            }                                                                \
        }

/**
******************************************************************************
 *  @brief   Macro to put a code with specified number of bits into the
 *           bitstream using 0th order exponential Golomb encoding for
 *           signed numbers
******************************************************************************
 */
#define PUT_BITS_UEV(ps_bitstrm, code_val, ret_val, syntax_string)           \
        {                                                                    \
            ENTROPY_TRACE(syntax_string, code_val);                          \
            ret_val = ih264e_put_uev((ps_bitstrm), (code_val));              \
            if(ret_val != IH264E_SUCCESS)                                    \
            {                                                                \
                return ret_val;                                              \
            }                                                                \
        }
/**
******************************************************************************
 *  @brief   Macro to put a code with specified number of bits into the
 *           bitstream using 0th order exponential Golomb encoding for
 *           signed numbers
******************************************************************************
 */
#define PUT_BITS_SEV(ps_bitstrm, code_val, ret_val, syntax_string)           \
        {                                                                    \
            ENTROPY_TRACE(syntax_string, code_val);                          \
            ret_val = ih264e_put_sev((ps_bitstrm), (code_val));              \
            if(ret_val != IH264E_SUCCESS)                                    \
            {                                                                \
                return ret_val;                                              \
            }                                                                \
        }

/**
******************************************************************************
 *  @brief   Macro to set active entropy threads to zero and return
 *           in case of errors
******************************************************************************
 */
#define RETURN_ENTROPY_IF_ERROR(ps_codec, ps_entropy, ctxt_sel)              \
        if(ps_entropy->i4_error_code != IH264E_SUCCESS)                      \
        {                                                                    \
            DATA_SYNC();                                                     \
            return ps_entropy->i4_error_code;                                \
        }

/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

WORD32 ih264e_generate_sps(bitstrm_t *ps_bitstrm, sps_t *ps_sps, vui_t *ps_vui);

WORD32 ih264e_generate_pps(bitstrm_t *ps_bitstrm, pps_t *ps_pps, sps_t *ps_sps);

IH264E_ERROR_T ih264e_generate_sei(bitstrm_t *ps_bitstrm,
                                   sei_params_t *ps_sei,
                                   UWORD32 u4_insert_per_idr);

WORD32 ih264e_generate_slice_header(bitstrm_t *ps_bitstrm,
                                    slice_header_t *ps_slice_hdr,
                                    pps_t *ps_pps,
                                    sps_t *ps_sps);

IH264E_ERROR_T ih264e_populate_sps(codec_t *ps_codec, sps_t *ps_sps);

IH264E_ERROR_T ih264e_populate_pps(codec_t *ps_codec, pps_t *ps_pps);

WORD32 ih264e_populate_slice_header(process_ctxt_t *ps_proc,
                                    slice_header_t *ps_slice_hdr,
                                    pps_t *ps_pps,
                                    sps_t *ps_sps);

IH264E_ERROR_T ih264e_add_filler_nal_unit(bitstrm_t *ps_bitstrm,
                                          WORD32 insert_fill_bytes);


#endif /* _IH264E_ENCODE_HEADER_H_ */
