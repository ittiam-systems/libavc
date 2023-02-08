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
/*!
 **************************************************************************
 * \file isvcd_nal_parse.h
 *
 * \brief
 *    Contains routines that resample for SVC resampling
 *
 * Detailed_description
 *
 * \date
 *
 *
 * \author : Kishore
 **************************************************************************
 */

#ifndef _SVCD_BITSTREAM_EXTRACT_H_
#define _SVCD_BITSTREAM_EXTRACT_H_

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

/*****************************************************************************/
/* Function Macros                                                           */
/*****************************************************************************/

/*****************************************************************************/
/* Typedefs                                                                  */
/*****************************************************************************/

/*****************************************************************************/
/* Enums                                                                     */
/*****************************************************************************/

typedef enum
{
    FULL_INPUT_MODE = 0,
    PARTIAL_INPUT_MODE
} NAL_PARSE_INPUT_MODE_T;

typedef enum
{
    VCL_NAL_FOUND_FALSE = 0,
    VCL_NAL_FOUND_TRUE
} EXTRACT_NON_VCL_NAL_RETURN_STATUS_T;

typedef enum
{
    PIC_BOUNDARY_FALSE,
    PIC_BOUNDARY_TRUE,
    FLUSH_DECODED_PICTURE
} EXTRACT_VCL_NAL_RETURN_STATUS_T;

/*****************************************************************************/
/* Structure                                                                 */
/*****************************************************************************/

/*****************************************************************************/
/* Extern Variable Declarations                                              */
/*****************************************************************************/

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/

void isvcd_nal_parse_reset_ctxt(WORD32 i4_input_bitstream_mode, WORD32 i4_input_mode,
                                void *pv_nal_parse_ctxt);

WORD32 isvcd_nal_parse_set_target_attr(WORD32 i4_target_quality_id, WORD32 i4_target_dependency_id,
                                       WORD32 i4_target_temporal_id, WORD32 i4_target_priority_id,
                                       void *pv_nal_parse_ctxt);

WORD32 isvcd_nal_parse_vcl_nal_partial(void *pv_nal_parse_ctxt, UWORD8 *pu1_stream_buffer,
                                       void *pv_out_non_vcl_nal, void *pv_out_vcl_nal,
                                       UWORD32 *pu4_bytes_consumed, UWORD32 *pu4_num_bytes);

WORD32 isvcd_nal_parse_non_vcl_nal(void *pv_nal_parse_ctxt, UWORD8 *pu1_stream_buffer,
                                   void *pv_out_non_vcl_nal, UWORD32 *pu4_bytes_consumed,
                                   UWORD32 *pu4_num_bytes);

WORD32 isvcd_nal_parse_partial_signal_eos(void *pv_nal_parse_ctxt, void *pv_out_vcl_nal,
                                          void *pv_out_non_vcl_nal);

#endif /* _SVCD_BITSTREAM_EXTRACT_H_ */
