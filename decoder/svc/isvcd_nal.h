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
 * \file isvcd_nal.h
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

#ifndef _SVCD_NAL_H_
#define _SVCD_NAL_H_

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

#define START_CODE_NOT_FOUND -1
#define END_OF_STREAM_BUFFER -2
#define END_OF_STREAM -1

#define SC_NOT_FOUND (-1)
#define SC_FOUND 1

#define NUM_OF_ZERO_BYTES_BEFORE_START_CODE (2)
#define START_CODE_BYTE (0x01)

#define EMULATION_PREVENTION_BYTE (0x03)

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
    NON_VCL_NAL,
    VCL_NAL
} DERIVED_NAL_UNIT_TYPE_T;

typedef enum
{
    NAL_START = 0,
    FIND_NAL_END,
    NAL_END
} NAL_BOUND_DETECT_STATE_T;

typedef enum
{
    STUFFED_BYTE = 0,
    /* Should be used for reset purposes */
    NOT_STUFFED_BYTE
} EMULATION_STATE_T;

typedef enum
{
    NAL_INSUFFICIENT_DATA = (WORD32) 0x80000000,
    NAL_CORRUPT_DATA = (WORD32) 0x80000001
} NAL_PARSE_ERR_CODES_T;

/*****************************************************************************/
/* Structure                                                                 */
/*****************************************************************************/
typedef struct
{
    WORD32 i4_nal_ref_idc;           /*!< NAL ref idc  - decoded prm from the bitstream */

    WORD32 i4_nal_unit_type;         /*!< NAL unit type - decoded prm from the
             bitstream                                                                */

    WORD32 i4_priority_id;           /*!< Priority id of NAL -  decoded prm from the
               bitstream. If not present then set to 0                                  */

    WORD32 i4_dependency_id;         /*!< dependency id of NAL -  decoded prm from the
             bitstream. If not present then set to 0                                  */

    WORD32 i4_quality_id;            /*!< Quality id of NAL - decoded prm from the
                bitstream. If not present then set to 0                                  */

    WORD32 i4_temporal_id;           /*!< Temporal id of NAL - decoded prm from the
               bitstream. If not present then set to 0                                  */

    WORD32 i4_no_int_lyr_pred;       /*!< No inter layer predictiion flag of NAL -
           decoded prm from the bitstream. if not present then it is set to 0       */

    WORD32 i4_use_ref_base_pic_flag; /*!<Use ref_base_pic flag  - decoded from
     the bitstream. If not present then it is set to 0                        */

    WORD32 i4_discard_flag;          /*!<Discard_flag - decoded from the bitstream.
              If not present then it is set to 0                                       */

    WORD32 i4_derived_nal_type;      /*! Derived NAL type - Place holder for enum
          which indicates whether the current NAL is VCL NAL or NON_VCL_NAL        */

    WORD32 i4_idr_pic_flag;          /*! Derived prm - Indicates whether current NAL
             is idr VCL NAL or not. SVCD_TRUE means IDR PIC NAL otherwsie not         */

    WORD32 i4_nal_header_len;        /*! lenght of NAL header in bytes              */

    WORD32 i4_dqid;                  /*! Derived parameter. It has same meaning as DQID in the
                      SVC standard                                                             */

    UWORD32 u4_first_mb_addr;        /*!< It shall hold the value of first MB address
                                     of the VCL NAL unit (excluding the prefix NAL
                                     unit) */

    UWORD8 u1_sps_id;                /*!< It shall have the value of SPS id for the VCL NAL
                                     unit (excluding the prefix NAL unit)*/

    UWORD8 u1_pps_id;                /*!< It shall have the value of PPS id for the VCL NAL
                                     unit (excluding the prefix NAL unit)*/

    UWORD16 u2_frm_num;              /*!< It shall have the value of frame number for the VCL
                                     NAL unit (excluding the prefix NAL unit)*/

    UWORD32 i4_idr_pic_num;          /*!< It shall have the value of IDR picture number when
                                    "i4_idr_pic_flag" is SVCD_TRUE for VCL NAL unit.
                                    (excluding prefix NAL unit) */

    WORD32 i4_poc_lsb;               /*!< It shall have the value of "picture order cnt lsb"
                                     when picture order count type is 0 for VCL NAL unit
                                     (excluding the prefix NAL unit). When not present in the
                                     bitstream, it shall be set to 0 */

    WORD32 i4_delta_poc_bot;         /*!< It shall have the value of "delta picture order
                               cnt bottom" when picture order count type is 0 for VCL NAL
                               unit (excluding the prefix NAL unit). When not present in
                               the bitstream, it shall be set to 0 */

    WORD32 ai4_delta_poc[2];         /*!< It shall have the value of "delta picture order
                               cnt bottom" when picture order count type is 1 for VCL NAL
                               unit (excluding the prefix NAL unit). When not present in
                               the bitstream, it shall be set to 0 */

} nal_prms_t;

/*****************************************************************************/
/* Extern Variable Declarations                                              */
/*****************************************************************************/

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/

void isvcd_nal_buf_reset(void *ps_nal_buf);

UWORD32 isvcd_nal_rbsp_to_sodb(UWORD8 *pu1_buf, WORD32 i4_nal_len_in_bytes, UWORD8 u1_ecd_mode);

WORD32 isvcd_get_first_start_code(UWORD8 *pu1_stream_buffer, UWORD32 *pu4_bytes_consumed,
                                  UWORD32 *pu4_num_bytes);

WORD32 isvcd_nal_find_start_code(UWORD8 *pu1_buf_start, WORD32 i4_cur_pos, WORD32 i4_max_num_bytes,
                                 WORD32 *pi4_zero_cnt, UWORD32 *pu4_bytes_consumed);

WORD32 isvcd_get_annex_b_nal_unit(UWORD8 *pu1_buf_start, WORD32 i4_cur_pos, WORD32 i4_max_num_bytes,
                                  WORD32 *pi4_state, WORD32 *pi4_zero_byte_cnt,
                                  UWORD32 *pu4_bytes_consumed, void *pv_nal_unit,
                                  WORD32 *pi4_more_data_flag);

void isvcd_reset_emulation_ctxt(void *pv_emulation_ctxt);

UWORD32 isvcd_nal_byte_swap_emulation(UWORD32 *pu4_out_stream, UWORD32 *pu4_out_len,
                                      UWORD8 *pu1_in_stream, UWORD32 u4_in_len, WORD32 i4_0s_bfr_sc,
                                      void *pv_emulation_ctxt);

void isvcd_dec_nal_hdr(void *pv_buf_ptr, WORD32 i4_buf_size, void *pv_nal_header_buf,
                       void *pv_nal_prms, void *pv_prefix_nal_buf, void *pv_prefix_nal_prms,
                       UWORD32 *pu4_err_code);

void isvcd_set_default_nal_prms(void *pv_nal_prms);

WORD32 isvcd_parse_part_slice_hdr(UWORD8 *pu1_input_buf, WORD32 i4_input_buf_size,
                                  UWORD8 *pu1_temp_buf, void *ps_sps, void *ps_pps,
                                  nal_prms_t *ps_nal_prms, UWORD32 *pu4_err_code,
                                  WORD32 *pi4_sps_pps_status);

WORD32 isvcd_discard_nal(void *pv_nal_prms, void *pv_app_attr, void *pv_int_attr,
                         WORD32 i4_update_flag);

#endif /* _SVCD_NAL_H_ */
