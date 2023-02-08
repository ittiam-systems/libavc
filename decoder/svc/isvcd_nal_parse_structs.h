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
 * \file isvcd_nal_parse_structs.h
 *
 * \brief
 *    Contains the definitions of structures used in the
 *      bitstream extraction module
 * Detailed_description
 *
 * \date
 *
 *
 * \author : Kishore
 **************************************************************************
 */

#ifndef _SVCD_NAL_PARSE_STRUCTS_H_
#define _SVCD_NAL_PARSE_STRUCTS_H_

/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/*****************************************************************************/
/* Constant Macros                                                           */
/*****************************************************************************/

#define HEADER_BUFFER_LEN_BEFORE_EP 32

#define MAX_NAL_HEADER_SIZE 4

#define UP_ALIGN_8(x) (((((UWORD64) x) + 7) >> 3) << 3)
#define ALIGN_4(x) (((x) + 3) & (~3))

/*--------------------------------------------------------------------------*/
/* The start address of any VCL or NON VCL internal buffers (input to        */
/* emulation prevention should be aligned to 4 byte boundary                */
/*--------------------------------------------------------------------------*/

#define BUFFER_ALIGN_4 4

#define FIRST_PASS 0
#define SECOND_PASS 1
/*****************************************************************************/
/* Function Macros                                                           */
/*****************************************************************************/

static __inline UWORD32 GET_NAL_BUF_INC(WORD32 i4_derived_nal_type)
{
    UWORD32 u4_buf_inc;

    if(VCL_NAL == i4_derived_nal_type)
    {
        u4_buf_inc = sizeof(vcl_buf_hdr_t);
    }
    else
    {
        u4_buf_inc = sizeof(non_vcl_buf_hdr_t);
    }

    u4_buf_inc = UP_ALIGN_8(u4_buf_inc);
    return (u4_buf_inc);
}

static __inline void UPDATE_NAL_BUF_PTR(UWORD8 **ppu1_buf, WORD32 i4_derived_nal_type,
                                        UWORD32 *pu4_bytes_left)
{
    UWORD8 *pu1_buf_ptr;
    UWORD64 u4_inc;

    /* Align the start of the structure */

    pu1_buf_ptr = *ppu1_buf;

    /* Account for the vcl or non-vcl header */
    u4_inc = GET_NAL_BUF_INC(i4_derived_nal_type);
    u4_inc = UP_ALIGN_8(u4_inc);
    pu1_buf_ptr += u4_inc;

    /* Update the pointers */
    *pu4_bytes_left -= u4_inc;
    *ppu1_buf = pu1_buf_ptr;
}

/*****************************************************************************/
/* Typedefs                                                                  */
/*****************************************************************************/

/*****************************************************************************/
/* Enums                                                                     */
/*****************************************************************************/

typedef enum
{
    NAL_PARSE_HANDLE = 0,
    NAL_PARSE_DQID_LIST_MEM,
    NAL_PARSE_CMN_MEM,
    NAL_PARSE_NAL_UNIT_MEM,
    NAL_PARSE_NUM_MEM_TABS
} BITSTREAM_EXTRACT_MEMTABS_T;

typedef enum
{
    PIC_BOUND_DQID = 0,      /* Second slice has lower DQID than the first slice */
    PIC_BOUND_SLICE_PRMS = 1 /* Second slice has different slice prms as */
                             /* as compared to first slice               */
} PIC_BOUNDARY_TYPES_T;

/*****************************************************************************/
/* Structure                                                                 */
/*****************************************************************************/

typedef struct
{
    vcl_node_t *ps_vcl_node; /*!< The pointer to VCL NAL node buffer.
                              */
    UWORD8 u1_valid_flag;    /*!< This flag shall indicate that the occupancy of
                             vcl node buffer. SVCD_TRUE shall indicate that the vcl
                             node buffer is occupied. @sa SVCD_BOOL_T
                             */

    UWORD8 u1_dqid;          /*!< The value of DQID assigned for this structure.
                             The range is [0,127] and is computed as
                             (Dependency id * 16 + Quality id )
                             */
    WORD32 i4_poc_lsb;       /*!< It shall have the value of "picture order cnt lsb"
                             when picture order count type is 0 for the layer. When not
                             present in the bitstream, it shall be set to 0*/

    WORD32
    i4_delta_poc_bot; /*!< It shall have the value of "delta picture order cnt
                 bottom" when picture order count type is 0 for VCL NAL unit.
                 When not present in the bitstream, it shall be set to 0*/

    WORD32
    ai4_delta_poc[2]; /*!< It shall have the value of "delta picture order cnt
                   bottom" when picture order count type is 1 for VCL NAL
                   unit.
                   When not present in the bitstream, itshall be set to 0 */
} dqid_node_t;

typedef struct
{
    WORD32 i4_max_num_lyrs;    /*!< Maximum number of layers that will be
                       present in a access unit. This will determine the
                       length of the VCL NAL node buffer. This parameter
                       is configurable during instance creation time
                       */

    dqid_node_t *ps_dqid_node; /*!< Pointer to start of VCL NAL node buffer.
                                */

} dqid_ctxt_t;

typedef struct
{
    WORD32 i4_valid_flag; /*!< It shall indicate the validity of contents of
                          this buffer structure. SVCD_TRUE shall indicate
                          that the contents of this structure is valid.
                          @sa SVCD_BOOL_T
                          */

    UWORD8 *pu1_buf;      /*!< It shall point to start of SODB data of the NAL.
                          It should be 8 byte aligned.
                          */

    UWORD32 u4_max_bits;  /*!< The length of SODB data of NAL in bits. This
                          should be set properly by taking care of whether NAL
                          is coded in CAVLC or CABAC or NAL is a NON VCL NAL
                          */

    WORD32 i4_buf_size;   /*!< The size of SODB data of NAL in bytes
                           */
} nal_buf_t;

typedef struct
{
    /*----------------------------------------------------*/
    /*---------- Mode of operation -----------------------*/
    /*----------------------------------------------------*/

    WORD32 i4_input_bitstream_mode; /*!< RFC or Annex B   */

    /*----------------------------------------------------*/
    /*---------- NAL boundary detection ------------------*/
    /*----------------------------------------------------*/

    WORD32 i4_find_nal_state;   /*!< state of NAL boundary
                       detection logic */

    WORD32 i4_zero_byte_cnt;    /*< Number of zero bytes consumed */

    WORD32 i4_dec_frst_sc_flag; /*!< A flag to decode
                    the start code only. A value of SVCD_TRUE
                    shall indicate that start code shall be
                    decoded.@sa SVCD_BOOL_T */

    /*----------------------------------------------------*/
    /*--------- Emulation prevention info ----------------*/
    /*----------------------------------------------------*/

    emulation_prevent_ctxt_t s_emulation_ctxt;

    /*----------------------------------------------------*/
    /*--------- Picture boundary detetction info ---------*/
    /*----------------------------------------------------*/

    WORD32 i4_is_frst_vcl_nal_in_au; /*!< Indicates whether
                     current NAL is first NAL in the current
                     Access unit. This is needed for detecting
                     picture boundary in partial input mode of
                     operation */

    UWORD32 u4_time_stamp_lsb;       /*!< Holds the LSB of time stamp of the
                         first NAL unit in the  access unit.
                         Used for RFC based bitstreams */

    WORD32 i4_time_stamp_msb;        /*!< Holds the MSB of time stamp of the
                         first NAL unit in the  access unit.
                         Used for RFC based bitstreams */

    WORD32 i4_prev_dq_id;            /*!< Holds the value of DQID of
                               last NAL unit parsed. this is used for
                               detetecting the picture boundary.
                               in Annex B mode of input bitstream */
    WORD32 i4_idr_pic_err_flag;      /*!< place to hold the
                          IDR status of current AU
                          */

    /*----------------------------------------------------*/
    /*-------- DQID node context -------------------------*/
    /*----------------------------------------------------*/

    dqid_ctxt_t s_dqid_ctxt;

    /*----------------------------------------------------*/
    /*-------- VCL and NON VCL buf info ------------------*/
    /*----------------------------------------------------*/

    void *pv_non_vcl_nal_buf;               /*!< Start address of NON VCL
                                  NAL buffer */

    void *pv_vcl_nal_buf;                   /*!< Start address of VCL NAL
                                    buffer */

    UWORD32 u4_bytes_left_vcl;              /*!< number of bytes left in the
                                            VCL buffer
                                            */
    UWORD32 u4_bytes_left_non_vcl;          /*!< number of bytes left in the
                                        NON VCL buffer
                                        */

    UWORD8 *pu1_non_vcl_nal_buf;            /*!< Current position of
                               non VCL NAL buffer pointer */

    UWORD8 *pu1_vcl_nal_buf;                /*!< Current position of VCL NAL
                                   buffer pointer */

    WORD32 i4_num_non_vcl_nals;             /*!< Number of non vcl nals */

    nal_buf_t s_prefix_nal_buf;             /*!< NAL buffer structure
                            of prefix NAL unit */

    nal_buf_t s_nal_buf;                    /*!< NAL buffer structure of .
                                     active NAL unit ( which is not a prefix
                                     NAL unit) */

    nal_buf_t *ps_nal_buf;                  /*!< It shall point to active
                                  NAL buffer structure. It shall point to
                                  either "s_prefix_nal_buf" or "s_nal_buf"*/

    vcl_buf_hdr_t *ps_prev_vcl_buf;         /*!< It shall point
                          to vcl buffer header of the previous
                          slice of a layer */
    non_vcl_buf_hdr_t *ps_prev_non_vcl_buf; /*!< It shall
                      point to buffer header of the previous
                      non vcl nal present in the bitstream */

    /*----------------------------------------------------*/
    /*-------- NAL structure and NAL buffer --------------*/
    /*----------------------------------------------------*/

    void *pv_nal_unit;
    void *pv_nal_header_buf;
    nal_prms_t s_nal_prms;
    nal_prms_t s_prefix_nal_prms;

    /*----------------------------------------------------*/
    /*-------------- Target layer info -------------------*/
    /*----------------------------------------------------*/

    target_lyr_attr_t s_app_attr; /*!< This structure shall have
                                  the target layer attributes set
                                  by the application */

    target_lyr_attr_t s_int_attr; /*!< This structure shall have
                                  the target layer attributes set
                                  by the module. At any state, the
                                  module tries to attain the values
                                  of application attributes at
                                  the IDR pictures */

    WORD32 i4_tgt_lyr_update;     /*!< It is a flag which
                                  indicates whether the internal target layer
                                  attributes has to be updated or not. A value
                                  of SVCD_TRUE shall indicate that the target
                                  layer attributes shall be updated.
                                  @sa SVCD_BOOL_T */

    /*----------------------------------------------------*/
    /*-------- other NAL info ----------------------------*/
    /*----------------------------------------------------*/

    WORD32 i4_discard_nal_flag;
    WORD32 i4_nal_type;

    /*----------------------------------------------------*/
    /*---------- Seq, pic prms buffers -------------------*/
    /*----------------------------------------------------*/

    void *pv_seq_prms;
    void *pv_pic_prms;

    /*----------------------------------------------------*/
    /*----------        Others          ------------------*/
    /*----------------------------------------------------*/

    WORD32 i4_eos_flag;    /*!< Flush mode related parameter.
                   This is used by the module during the
                   flush call. SVCD_TRUE shall indicate that
                   current  end of bitstream has occurred.
                   @sa SVCD_BOOL_T; */

    UWORD32 u4_bytes_left; /*!< This field has number of bytes not
                    consumed by the NAL parse module in the
                    previous call because of insufficient bitstream
                    to decode the slice and NAL header. */
    UWORD8 u1_pic_boundary_aud_flag;
} nal_parse_ctxt_t;

/*****************************************************************************/
/* Extern Variable Declarations                                              */
/*****************************************************************************/

/*****************************************************************************/
/* Extern Function Declarations                                              */
/*****************************************************************************/

#endif /* _SVCD_NAL_PARSE_STRUCTS_H_ */
