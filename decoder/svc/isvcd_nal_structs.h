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
 * \file isvcd_nal_structs.h
 *
 * \brief
 *    Contains the definitions of structures used in NAL processing
 *
 * Detailed_description
 *
 * \date
 *
 *
 * \author : Kishore
 **************************************************************************
 */

#ifndef _SVCD_NAL_STRUCTS_H_
#define _SVCD_NAL_STRUCTS_H_

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

#define MAX_NUM_SLICE_GRPS_IN_LYR 8

#define MAX_NUM_PACKETS_IN_NAL                    \
    100 /*! Maximum number of packets that can be \
  present in a NAL unit */

/*****************************************************************************/
/* Function Macros                                                           */
/*****************************************************************************/
#define SLICE_PARSE_ERR_HDLR(i4_error, u4_err_code, pu4_err_code) \
    if(0 != i4_error)                                             \
    {                                                             \
        *pu4_err_code = u4_err_code;                              \
        return (NOT_OK);                                          \
    }

#define SLICE_PARSE_ERR_HDLR1(i4_error, u4_err_code, pu4_err_code, pi4_sps_pps, i4_sps_pps_err) \
    if(0 != i4_error)                                                                           \
    {                                                                                           \
        *pu4_err_code = u4_err_code;                                                            \
        *pi4_sps_pps = i4_sps_pps_err;                                                          \
        return (NOT_OK);                                                                        \
    }

/*****************************************************************************/
/* Typedefs                                                                  */
/*****************************************************************************/

/*****************************************************************************/
/* Enums                                                                     */
/*****************************************************************************/

typedef enum
{
    NAL_DISCARD_RESET_STATE,
    NAL_DISCARD_ACTIVE_STATE
} NAL_DISCARD_STATE_MACHINE_T;

/*****************************************************************************/
/* Structure                                                                 */
/*****************************************************************************/

typedef struct
{
    WORD32 i4_dependency_id; /*!< Target dependency id */
    WORD32 i4_quality_id;    /*!< Target quality id */
    WORD32 i4_temporal_id;   /*!< Target temporal id */
    WORD32 i4_priority_id;   /*!< Target priority id */
} target_lyr_attr_t;

typedef struct
{
    WORD32 i4_num_bufs;  /*!< Number of buffers that comprises this NAL unit.
      In case of Annex B based input, this value will always be 1. Otherwise
      (RFC - input), this value indicates number of packets in the NAL unit
      in the current process call                                              */

    UWORD8 *pu1_bufs;    /*!< Nal unit buffer pointer                           */
    WORD32 i4_buf_sizes; /*!< Nal unit buffer size                            */

} nal_unit_t;

typedef struct
{
    WORD32 i4_state;          /*!< State of emulation prevention in the state machine */
    WORD32 i4_zeroes_cnt;     /*!< Number of consecutive zeroes counter */
    UWORD32 u4_word;          /*!< Place holder for WORD - output */
    UWORD32 u4_bytes_in_word; /*!< Number of bytes in the WORD */
} emulation_prevent_ctxt_t;

typedef struct vcl_buf_hdr_t
{
    struct vcl_buf_hdr_t *ps_next; /*!< Pointer to next VCL NAL buffer header.
                          This will be poiting to next slice in the layer. This
                          field shall be set to NULL for the last slice
                          (VCL NAL) in a layer */

    UWORD32 u4_max_bits;           /*!< Total number of SODB bits present in the VCL
                                   NAL (slice)*/

    WORD32 i4_buf_offset;          /*!< This is the offset from the start of the
                                   VCL NAL header to start of SODB data of VCL NAL.
                                   This shall be multiple of 4*/

    WORD32 i4_slice_offset;        /*!< It is the offset from start of VCL NAL data
                                   to start of slice data. A value of 0 shall
                                   indicate that prefix NAL unit is
                                   not present
                                   Note: If prefix NAL unit is present then it will
                                   be present between the offsets "i4_buf_offset" to
                                   "i4_slice_offset"*/

    UWORD32 u4_prefix_nal_bits;    /*!< Total number of SODB bits present in the VCL
                            NAL (slice). This shall have valid value when
                           "i4_slice_offset" has non zero value
                            */

    WORD32 i4_no_int_lyr_pred;     /*!< The value of no inter layer prediction
                            which is decoded from the NAL header
                            */

    WORD32 i4_first_mb_addr;       /*!< The 'first_mb_address' syntax
                                   element value decoded form the slice
                                   header present in the VCL NAL
                                   */
} vcl_buf_hdr_t;

typedef struct non_vcl_buf_hdr_t
{
    struct non_vcl_buf_hdr_t *ps_next; /*!< This shall point to start of next NON
                          VCL buffer header that is extracted from the bitstream.
                          It shall be set to NULL for the last non VCL NAL
                           */

    WORD32 i4_nal_unit_type;           /*!< NAL unit type that is decoded from the NAL
                                       header
                                      */

    WORD32 i4_buf_offset;              /*!< This is the offset from the start of the
                                            VCL NAL header to start of SODB data of VCL NAL.
                                            This shall be multiple of 4
                                         */

    WORD32 i4_buf_size;                /*!< Size of the NON VCL SODB data in bytes
                                        */

} non_vcl_buf_hdr_t;
typedef struct vcl_node_t
{
    struct vcl_node_t *ps_top_node; /*!< Pointer to top node present in the DQID
                                    list. This node is actually is a layer using
                                    the current layer as a reference layer. Value
                                    of NULL shall indicate that no more layers
                                    with DQID higher than current layer is present
                                    in the current access unit */

    struct vcl_node_t *ps_bot_node; /*!< Pointer to bottom node present in the
                                      DQID list. This node is actually the
                                      reference layer of the current layer. Value
                                      of NULL shall indicate that no more layers
                                      with DQID lower than current layer is
                                      present in the current access unit */

    /*------ info part -------*/

    WORD32 i4_quality_id;    /*!< Quality id of the layer */

    WORD32 i4_dependency_id; /*!< Dependency id of the layer */

    WORD32 i4_temporal_id;   /*!< Temporal id of the layer */

    WORD32 i4_priority_id;   /*!< Priority id of the layer */

    WORD32 i4_idr_pic_flag;  /*!< Flag indicating whether current layer is
                                  Idr picture or not. SVCD_TRUE shall indicate
                                  the idr picture
                              */

    WORD32 i4_nal_unit_type; /*!< NAL unit type of all slices in the current
                                  picture
                              */

    WORD32
    i4_nal_ref_idc;                      /*!< NAL ref idc of all slices in the current picture */

    WORD32 i4_use_ref_base;              /*!< Use ref base flag of NAL header. */

    UWORD8 u1_sps_id;                    /*!< It shall have the value of SPS id used by this layer.
                                             It's range is [0,63]
                                          */

    UWORD8 u1_pps_id;                    /*!< It shall have the value of PPS id used by this layer.
                                             It's range is [0,255]
                                          */
    UWORD8 u1_acc_no_int_pred;           /*! The value of accumulated no inter layer
                                            prediction flag. This value will be "logical and
                                            " of no inter layer prediction flag of all the
                                            slices in the corresponding DQID
                                          */

    UWORD16 u2_frm_num;                  /*!< It is the value of frame number of the layer.
                                          */

    UWORD32 i4_idr_pic_num;              /*!< It shall have the value of IDR picture number when
                                            "i4_idr_pic_flag" is SVCD_TRUE
                                         */

    WORD32 i4_poc_syntax;                /*!< It shall have either "picture order cnt lsb" or
                                         "delta picture order cnt [0]" that is decoded from the
                                         slice header. When picture order coutn type is 0 then
                                         this field holds "picture order cnt lsb" and when
                                         picture order cnt type is 1 then this field holds
                                         "delta picture order cnt [0]". This field will not
                                         have a valid value when picture order cnt type is 2
                                         */

    WORD32 i4_res_change_flag;           /*!< 'SpatialResolutionChangeFlag' for the
                                              layer as specified in section G.7.4.3.4
                                              Value 'SVCD_TRUE' indicates that parameter
                                              'SpatialResolutionChangeFlag' is set to 1.
                                              Otherwise it shall be set to SVCD_FALSE.
                                              @sa SVCD_BOOL_T
                                          */
    WORD32 i4_ref_dq_id;                 /*!< reference layer DQid for current layer */

    WORD32 i4_num_slices;                /*!< Number of slices in the current layer */

    WORD32 i4_inter_lyr_dblk_idc;        /*!< Deblock filter idc for inter layer
                                        deblocking. This shall be valid only
                                        for layer with quality id = 0.
                                        */
    WORD32 i4_inter_lyr_alpha_c0_offset; /*!< Alpha C0 offset for inter layer
                                 deblocking. This shall be valid only for
                                 layers with quality id = 0.
                                 */
    WORD32 i4_inter_lyr_beta_offset;     /*!< Beta offset for inter layer
                                     deblocking. This shall be valid only for
                                     layers with quality id = 0.
                                     */

    vcl_buf_hdr_t *ps_first_vcl_nal;     /*!< This shall point to start of the
                                   VCL NAL header of the first slice (VCL NAL)
                                   in a layer.
                                   */

    vcl_buf_hdr_t *aps_start_addr_slices[MAX_NUM_SLICE_GRPS_IN_LYR]; /*!< array
                                to hold the start address of first slice
                                of each slice group. the address will be
                                linked to each other within a slice
                                group based on first MB address
                                each entry will be pointing to the slice
                                which will be decoded next in the slice
                                group
                                */

} vcl_node_t;
/*****************************************************************************/
/* Extern Variable Declarations                                              */
/*****************************************************************************/

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/

#endif /* _SVCD_NAL_STRUCTS_H_ */
