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
******************************************************************************
* @file
*  isvce_encode_header.h
*
* @brief
*  This file contains structures and interface prototypes for h264 bitstream
*  header encoding
*
* @author
*  ittiam
*
* @remarks
*  None
*
*******************************************************************************
*/

#ifndef _ISVCE_ENCODE_HEADER_H_
#define _ISVCE_ENCODE_HEADER_H_

#include "ih264_typedefs.h"

/* Dependencies of ih264e_bitstream.h */
#include "ih264e_error.h"

#include "ih264e_bitstream.h"
#include "ih264e_trace.h"
#include "isvce_structs.h"

/**
******************************************************************************
*  @brief   Macro to put a code with specified number of bits into the
*           bitstream
******************************************************************************
*/
#define PUT_BITS(ps_bitstrm, code_val, code_len, ret_val, syntax_string) \
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
#define PUT_BITS_UEV(ps_bitstrm, code_val, ret_val, syntax_string) \
    {                                                              \
        ENTROPY_TRACE(syntax_string, code_val);                    \
        ret_val = ih264e_put_uev((ps_bitstrm), (code_val));        \
        if(ret_val != IH264E_SUCCESS)                              \
        {                                                          \
            return ret_val;                                        \
        }                                                          \
    }
/**
******************************************************************************
*  @brief   Macro to put a code with specified number of bits into the
*           bitstream using 0th order exponential Golomb encoding for
*           signed numbers
******************************************************************************
*/
#define PUT_BITS_SEV(ps_bitstrm, code_val, ret_val, syntax_string) \
    {                                                              \
        ENTROPY_TRACE(syntax_string, code_val);                    \
        ret_val = ih264e_put_sev((ps_bitstrm), (code_val));        \
        if(ret_val != IH264E_SUCCESS)                              \
        {                                                          \
            return ret_val;                                        \
        }                                                          \
    }

/**
******************************************************************************
*  @brief   Macro to set active entropy threads to zero and return
*           in case of errors
******************************************************************************
*/
#define RETURN_ENTROPY_IF_ERROR(ps_codec, ps_entropy, ctxt_sel) \
    if(ps_entropy->i4_error_code != IH264E_SUCCESS)             \
    {                                                           \
        DATA_SYNC();                                            \
        ps_codec->au4_entropy_thread_active[ctxt_sel] = 0;      \
        return ps_entropy->i4_error_code;                       \
    }

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/
extern WORD32 ih264e_generate_nal_unit_header(bitstrm_t *ps_bitstrm, WORD32 nal_unit_type,
                                              WORD32 nal_ref_idc);

extern WORD32 ih264e_generate_vui(bitstrm_t *ps_bitstrm, vui_t *ps_vui);

extern IH264E_ERROR_T ih264e_generate_sei(bitstrm_t *ps_bitstrm, sei_params_t *ps_sei,
                                          UWORD32 u4_insert_per_idr);

extern IH264E_ERROR_T ih264e_add_filler_nal_unit(bitstrm_t *ps_bitstrm, WORD32 insert_fill_bytes);

/**
******************************************************************************
*
* @brief Generates SPS (Sequence Parameter Set)
*
* @par   Description
*  This function generates Sequence Parameter Set header as per the spec
*
* @param[in]   ps_bitstrm
*  pointer to bitstream context (handle)
*
* @param[in]   ps_sps
*  pointer to structure containing SPS data
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_generate_sps(bitstrm_t *ps_bitstrm, sps_t *ps_sps, NAL_UNIT_TYPE_T nal_type);

/**
******************************************************************************
*
* @brief Generates PPS (Picture Parameter Set)
*
* @par   Description
*  Generate Picture Parameter Set as per Section 7.3.2.2
*
* @param[in]   ps_bitstrm
*  pointer to bitstream context (handle)
*
* @param[in]   ps_pps
*  pointer to structure containing PPS data
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_generate_pps(bitstrm_t *ps_bitstrm, pps_t *ps_pps, sps_t *ps_sps);

/**
******************************************************************************
*
* @brief Generates Slice Header
*
* @par   Description
*  Generate Slice Header as per Section 7.3.5.1
*
* @param[inout]   ps_bitstrm
*  pointer to bitstream context for generating slice header
*
* @param[in]   ps_slice_hdr
*  pointer to slice header params
*
* @param[in]   ps_pps
*  pointer to pps params referred by slice
*
* @param[in]   ps_sps
*  pointer to sps params referred by slice
*
* @param[out]   ps_dup_bit_strm_ent_offset
*  Bitstream struct to store bitstream state
*
* @param[out]   pu4_first_slice_start_offset
*  first slice offset is returned
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_generate_slice_header(bitstrm_t *ps_bitstrm, slice_header_t *ps_slice_hdr,
                                   pps_t *ps_pps, sps_t *ps_sps, UWORD8 u1_idr_flag);
/**
******************************************************************************
*
* @brief Populates sps structure
*
* @par   Description
*  Populates sps structure for its use in header generation
*
* @param[in]   ps_codec
*  pointer to encoder context
*
* @param[out]  ps_sps
*  pointer to sps params that needs to be populated
*
* @return      success or failure error code
*
******************************************************************************
*/
IH264E_ERROR_T isvce_populate_sps(isvce_codec_t *ps_codec, sps_t *ps_sps, UWORD8 u1_sps_id,
                                  UWORD8 u1_profile_idc, isvce_inp_buf_t *ps_inp_buf,
                                  UWORD8 u1_spatial_layer_id);

/**
******************************************************************************
*
* @brief Populates pps structure
*
* @par   Description
*  Populates pps structure for its use in header generation
*
* @param[in]   ps_codec
*  pointer to encoder context
*
* @param[out]  ps_pps
*  pointer to pps params that needs to be populated
*
* @return      success or failure error code
*
******************************************************************************
*/
IH264E_ERROR_T isvce_populate_pps(isvce_codec_t *ps_codec, pps_t *ps_pps, UWORD8 u1_sps_id,
                                  UWORD8 u1_pps_id, UWORD8 u1_spatial_layer_id);

/**
******************************************************************************
*
* @brief Populates slice header structure
*
* @par   Description
*  Populates slice header structure for its use in header generation
*
* @param[in]  ps_proc
*  pointer to proc context
*
* @param[out]  ps_slice_hdr
*  pointer to slice header structure that needs to be populated
*
* @param[in]  ps_pps
*  pointer to pps params structure referred by the slice
*
* @param[in]   ps_sps
*  pointer to sps params referred by the pps
*
* @return      success or failure error code
*
******************************************************************************
*/
WORD32 isvce_populate_slice_header(isvce_process_ctxt_t *ps_proc, slice_header_t *ps_slice_hdr,
                                   pps_t *ps_pps, sps_t *ps_sps, UWORD8 u1_is_idr);

extern WORD32 isvce_populate_svc_nalu_extension(isvce_process_ctxt_t *ps_proc,
                                                svc_nalu_ext_t *ps_svc_nalu_ext,
                                                NAL_UNIT_TYPE_T nalu_type, UWORD8 u1_idr_flag);

extern WORD32 isvce_generate_svc_nalu_extension(bitstrm_t *ps_bitstrm,
                                                svc_nalu_ext_t *ps_svc_nalu_ext, UWORD8 u1_nalu_id);

extern WORD32 isvce_populate_svc_slice(isvce_process_ctxt_t *ps_proc,
                                       svc_slice_header_t *ps_svc_slice_hdr, pps_t *ps_pps,
                                       subset_sps_t *ps_subset_sps,
                                       svc_nalu_ext_t *ps_svc_nalu_ext);

extern WORD32 isvce_populate_subset_sps(isvce_codec_t *ps_codec, subset_sps_t *ps_subset_sps,
                                        UWORD8 u1_sps_id, isvce_inp_buf_t *ps_inp_buf,
                                        UWORD8 u1_spatial_layer_id);

extern WORD32 isvce_generate_prefix_nal(bitstrm_t *ps_bitstrm, svc_nalu_ext_t *ps_svc_nalu_ext,
                                        slice_header_t *ps_slice_header,
                                        UWORD8 u1_max_num_ref_frames, UWORD8 u1_num_spatial_layers);

extern WORD32 isvce_generate_slice_header_svc(bitstrm_t *ps_bitstrm, pps_t *ps_pps,
                                              svc_nalu_ext_t *ps_svc_nalu_ext,
                                              svc_slice_header_t *ps_svc_slice_hdr,
                                              subset_sps_t *ps_subset_sps);

extern WORD32 isvce_generate_subset_sps(bitstrm_t *ps_bitstrm, subset_sps_t *ps_subset_sps);

#endif
