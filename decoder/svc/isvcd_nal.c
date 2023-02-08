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
 * \file isvcd_nal.c
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
/*****************************************************************************/
/*                                                                           */
/*  File Name         : isvcd_nal.c                                          */
/*                                                                           */
/*  Description       : Contains fucntions which help in NAL extraction from */
/*                      the bitstream                                        */
/*                                                                           */
/*  List of Functions : isvcd_nal_find_start_code,                           */
/*                      isvcd_get_annex_b_nal_unit,                          */
/*                      isvcd_get_rfc_nal_unit,                              */
/*                      isvcd_nal_rbsp_to_sodb,                              */
/*                      isvcd_reset_emulation_ctxt,                          */
/*                      isvcd_nal_byte_swap_emulation,                       */
/*                      isvcd_set_default_nal_header_prms,                   */
/*                      isvcd_dec_nal_hdr,                                   */
/*                      isvcd_parse_part_slice_hdr,                          */
/*                      isvcd_get_int_tgt_lyr_attr,                          */
/*                      isvcd_discard_nal                                    */
/*                                                                           */
/*  Issues / Problems : None                                                 */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          14 09 2021   Kishore         Draft                               */
/*                                                                           */
/*****************************************************************************/
/*****************************************************************************/
/* File Includes                                                             */
/*****************************************************************************/

/* System include files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stddef.h>
#include <assert.h>

/* standard interface include files */
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_tables.h"
#include "iv.h"
#include "ivd.h"
#include "ih264d_defs.h"
#include "ih264_debug.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_inter_pred.h"
#include "isvcd_structs.h"
#include "ih264d_nal.h"
#include "ih264d_error_handler.h"
#include "ih264d_defs.h"

/*****************************************************************************/
/*Extern Variable Declarations                                               */
/*****************************************************************************/

/*****************************************************************************/
/* Global Variable Definitions                                               */
/*****************************************************************************/

/*****************************************************************************/
/* Static Global Variable Definitions                                        */
/*****************************************************************************/

/*****************************************************************************/
/* Static function Definitions                                               */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_reset_nal_buf                                      */
/*                                                                           */
/*  Description   : Performs the reset of NAL buffer structure               */
/*  Inputs        : 1. Pointer to NAL buffer structure                       */
/*  Globals       : None                                                     */
/*  Processing    : Updates different fields of the structure                */
/*  Outputs       : None                                                     */
/*  Returns       :                                                          */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/
void isvcd_nal_buf_reset(void *pv_nal_buf)
{
    nal_buf_t *ps_nal_buf = pv_nal_buf;

    ps_nal_buf->i4_valid_flag = SVCD_FALSE;
    ps_nal_buf->i4_buf_size = 0;
    ps_nal_buf->u4_max_bits = 0;
    ps_nal_buf->pu1_buf = NULL;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name :svcd_nal_find_start_code                                  */
/*                                                                           */
/*  Description   : Finds the position of the start code in the stream       */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. Pointer to buffer start                               */
/*                  2. start position                                        */
/*                  3. Maximum number of bytes in the buffer                 */
/*                  4. pointer to zero byte count                            */
/*                  5. pointer to bytes consumed variable                    */
/*  Globals       :                                                          */
/*  Processing    : Searches for the start code in the bitstream and updates */
/*                  consumed variable                                        */
/*                                                                           */
/*  Outputs       : Bytes consumed variable                                  */
/*  Returns       : If start code is found then it returns SC_FOUND otherwise*/
/*                  it returns SC_NOT_FOUND                                  */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_nal_find_start_code(UWORD8 *pu1_buf_start, WORD32 i4_cur_pos, WORD32 i4_max_num_bytes,
                                 WORD32 *pi4_zero_cnt, UWORD32 *pu4_bytes_consumed)
{
    UWORD8 *pu1_buf = pu1_buf_start + i4_cur_pos;
    WORD32 i4_i;

    for(i4_i = 0; i4_i < (i4_max_num_bytes - i4_cur_pos); i4_i++)
    {
        /*-------------------------------------------------------------------*/
        /* If zero increment the zero byte counter                           */
        /*-------------------------------------------------------------------*/
        if(0 == *pu1_buf)
        {
            (*pi4_zero_cnt)++;
        }

        /*-------------------------------------------------------------------*/
        /* If start code found then increment the byte consumed and return   */
        /*-------------------------------------------------------------------*/
        else if(0x01 == *pu1_buf && *pi4_zero_cnt >= NUM_OF_ZERO_BYTES_BEFORE_START_CODE)
        {
            (*pu4_bytes_consumed)++;
            return (SC_FOUND);
        }
        /*-------------------------------------------------------------------*/
        /* If non zero byte and value is not equal to 1 a then reset zero    */
        /* byte counter                                                      */
        /*-------------------------------------------------------------------*/
        else
        {
            *pi4_zero_cnt = 0;
        }

        (*pu4_bytes_consumed)++;
        pu1_buf++;
    }

    return (SC_NOT_FOUND);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_first_start_code                               */
/*                                                                           */
/*  Description   : Searches for the first start code in the bitstream       */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. input buffer structure                                */
/*                  2. Bytes consumed variable                               */
/*  Globals       : None                                                     */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : Updates bytes consumed variable                          */
/*  Returns       : Start code is found or not                               */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_get_first_start_code(UWORD8 *pu1_stream_buffer, UWORD32 *pu4_bytes_consumed,
                                  UWORD32 *pu4_num_bytes)
{
    WORD32 i4_zero_cnt = 0, i4_status;
    UWORD32 u4_bytes_consumed_temp = 0;

    i4_status = isvcd_nal_find_start_code(pu1_stream_buffer, 0, *pu4_num_bytes, &i4_zero_cnt,
                                          &u4_bytes_consumed_temp);

    /*-----------------------------------------------------------------------*/
    /* If start code is not found then return and start searching for it     */
    /* again in the next process call. This process is repeated till we      */
    /* get a start code                                                      */
    /*-----------------------------------------------------------------------*/
    if(SC_NOT_FOUND == i4_status)
    {
        *pu4_bytes_consumed += u4_bytes_consumed_temp;
        return (i4_status);
    }
    else
    {
        /*-------------------------------------------------------------------*/
        /* If start code found then proceed with bitstream extraction        */
        /*-------------------------------------------------------------------*/
        *pu4_bytes_consumed += u4_bytes_consumed_temp;
        return (i4_status);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_annex_b_nal_unit                               */
/*                                                                           */
/*  Description   : This function gets one NAL unit from the Annex B based   */
/*                  input bitstream                                          */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. Input buffer pointer                                  */
/*                  2. Current position in the input buffer                  */
/*                  3. Input buffer size                                     */
/*                  4. Pointer to state of NAL boundary detection variable   */
/*                  5. Pointer to bytes consumed variable                    */
/*                  6. pointer to nal structure                              */
/*  Globals       :                                                          */
/*  Processing    : This fucntion searches for start code from the current   */
/*                  position and once gets one start code it searches for    */
/*                  another start code to get a NAL unit.                    */
/*                                                                           */
/*  Outputs       : Updates the state of NAL boundary detection logic        */
/*                  Updates the bytes consumed variable from 0 to bytes      */
/*                  consumed in this call                                    */
/*  Returns       : start of nal flag                                        */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_get_annex_b_nal_unit(UWORD8 *pu1_buf_start, WORD32 i4_cur_pos, WORD32 i4_max_num_bytes,
                                  WORD32 *pi4_state, WORD32 *pi4_zero_byte_cnt,
                                  UWORD32 *pu4_bytes_consumed, void *pv_nal_unit,
                                  WORD32 *pi4_more_data_flag)
{
    nal_unit_t *ps_nal_unit = (nal_unit_t *) pv_nal_unit;
    WORD32 i4_status, i4_nal_start_flag = SVCD_FALSE;

    /*-----------------------------------------------------------------------*/
    /* Initialization                                                        */
    /*-----------------------------------------------------------------------*/
    *pu4_bytes_consumed = 0;
    *pi4_more_data_flag = SVCD_TRUE;

    /*------------------------ check ----------------------------------------*/
    /* Assumptions is that this fucntion should not be called with this state*/
    /* hence it is responsibility of the caller to reset the state after the */
    /* NAL_END.                                                              */
    /*-----------------------------------------------------------------------*/
    if(NAL_END == *pi4_state)
    {
        return i4_nal_start_flag;
    }

    /*-----------------------------------------------------------------------*/
    /* ps_nal_unit->apu1_bufs[0] is expected to point to start of buffer of  */
    /* current NAL unit of the current process call. If a NAL unit is frag-  */
    /* -mented across multiple process call then this buffer should point to */
    /* start address of buffers. But when start of NAL is present in the     */
    /* buffer of current process call then ps_nal_unit->apu1_bufs[0] is      */
    /* expected to point to start adress of NAL unit (should be pointing to) */
    /* NAL header)                                                           */
    /*-----------------------------------------------------------------------*/
    ps_nal_unit->pu1_bufs = pu1_buf_start + i4_cur_pos;

    if(NAL_START == *pi4_state)
    {
        if(0 != *pi4_zero_byte_cnt)
        {
            return i4_nal_start_flag;
        }
        i4_nal_start_flag = SVCD_TRUE;
        ps_nal_unit->i4_num_bufs = 1;
        ps_nal_unit->i4_buf_sizes = 0;
        *pi4_state = FIND_NAL_END;
    }

    i4_status = isvcd_nal_find_start_code(pu1_buf_start, i4_cur_pos, i4_max_num_bytes,
                                          pi4_zero_byte_cnt, pu4_bytes_consumed);

    if(SC_NOT_FOUND == i4_status)
    {
        /*-------------------------------------------------------------------*/
        /* If start code is not found then there are 2 possibilities         */
        /* 1. We are in the middle of decoding the start code. This means    */
        /*    that we might have decoded the one or 2 zeroes of the start    */
        /*    code. In such cases, we should not consume these bytes. Though */
        /*    doing so we might encounter spurious cases where 0's are not   */
        /*    actually corresponds to start code but these will not harm us  */
        /* 2. Not of above case. Straightforward one                         */
        /*-------------------------------------------------------------------*/
        ps_nal_unit->i4_buf_sizes = *pu4_bytes_consumed;
        *pi4_more_data_flag = SVCD_FALSE;

        return (i4_nal_start_flag);
    }
    else
    {
        /*-------------------------------------------------------------------*/
        /* If NAL END is found then increment the bytes consumed appropriatly*/
        /* reset the zero byte counter                                       */
        /*-------------------------------------------------------------------*/
        *pi4_state = NAL_END;
        ps_nal_unit->i4_buf_sizes = *pu4_bytes_consumed - 1;
        *pi4_zero_byte_cnt = 0;
        return (i4_nal_start_flag);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_rbsp_to_sodb                                   */
/*                                                                           */
/*  Description   : Converts the RBSP data to SODB data                      */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. Input buffer containing the NAL unit                  */
/*                  2. Length of NAL unit (in bytes)                         */
/*  Globals       : None                                                     */
/*  Processing    : Finds the RBSP stop bit, if present then finds the length*/
/*                  of SODB data                                             */
/*                                                                           */
/*  Outputs       :                                                          */
/*  Returns       : Number of bits in the SODB data                          */
/*                                                                           */
/*  Issues        :                                                          */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/

UWORD32 isvcd_nal_rbsp_to_sodb(UWORD8 *pu1_buf, WORD32 i4_nal_len_in_bytes, UWORD8 u1_ecd_mode)
{
    UWORD32 u4_last_word_pos;
    UWORD32 u4_word, u4_max_bit_offset;
    UWORD8 i4_num_bits;
    WORD32 i4_i;
    WORD64 i8_nal_len;
    UWORD32 *pu4_buf;

    if(0 >= i4_nal_len_in_bytes)
    {
        return (0);
    }

    /* Get offset in bits */
    i8_nal_len = (WORD64) i4_nal_len_in_bytes << 3;
    u4_max_bit_offset = (UWORD32) i8_nal_len;

    /* If NAL is coded in CABAC then SODB */
    /* length has to account for CABAC    */
    /* ZERO WORDS also                    */
    if(1 == u1_ecd_mode)
    {
        return (u4_max_bit_offset);
    }

    /* Calculate the position of last word */
    u4_last_word_pos = i4_nal_len_in_bytes >> 2;

    /* Load the last word                 */
    i4_i = i4_nal_len_in_bytes & 0x03;
    if(0 != i4_i)
    {
        pu4_buf = (UWORD32 *) pu1_buf;
        pu4_buf += u4_last_word_pos;
        u4_word = *pu4_buf;
        i4_num_bits = i4_i << 3;
        u4_word >>= (32 - i4_num_bits);
    }
    else
    {
        pu4_buf = (UWORD32 *) pu1_buf;
        pu4_buf += (u4_last_word_pos - 1);
        u4_word = *pu4_buf;
        i4_num_bits = 32;
    }

    /* Search for RBSP stop bit          */
    do
    {
        for(i4_i = 0; (i4_i < i4_num_bits) && !CHECKBIT(u4_word, i4_i); i4_i++)
            ;

        u4_max_bit_offset -= i4_i;

        /* RBSP stop bit is found then   */
        /* come out of the loop          */
        if(0 != CHECKBIT(u4_word, i4_i))
        {
            /* Remove RBSP stop bit */
            u4_max_bit_offset -= 1;
            break;
        }

        pu4_buf -= 1;
        u4_word = *pu4_buf;
        i4_num_bits = 32;
    } while(u4_max_bit_offset > 0);

    return (u4_max_bit_offset);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_reset_emulation_ctxt                               */
/*                                                                           */
/*  Description   : Resets the emulation prevention context structure        */
/*                                                                           */
/*  Inputs        : pv_emulation_ctxt - pointer to emulation prevention      */
/*                      context structure                                    */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : None                                                     */
/*                                                                           */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay      Draft                                    */
/*                                                                           */
/*****************************************************************************/

void isvcd_reset_emulation_ctxt(void *pv_emulation_ctxt)
{
    emulation_prevent_ctxt_t *ps_emulation_ctxt = (emulation_prevent_ctxt_t *) pv_emulation_ctxt;

    /*! Reset the emulation prevention context */
    ps_emulation_ctxt->i4_state = NOT_STUFFED_BYTE;
    ps_emulation_ctxt->i4_zeroes_cnt = 0;
    ps_emulation_ctxt->u4_bytes_in_word = 0;
    ps_emulation_ctxt->u4_word = 0;
}

/****************************************************************************/
/*                                                                          */
/* Function Name  : isvcd_nal_byte_swap_emulation                           */
/*                                                                          */
/* Description    : This function is does byte swap or emulation or both    */
/*                  in the stream.                                          */
/*                                                                          */
/* Inputs         :  pu4_out_stream : Pointer to bitstream out buffer       */
/*                   pu4_out_len    : Pointer to variable for out len       */
/*                   pu1_in_stream  : Pointer to bitstream in buffer        */
/*                   u4_in_len      : Input bitstream buffer length         */
/*                   u4_prev_0s     : In case of fragemented NAL 0s in last */
/*                                    fragmented unit                       */
/*                   u4_0s_bfr_sc   : Number of zeros before start code     */
/*                   u4_bytes       : Number of bytes in last fragmented    */
/*                                    word                                  */
/*                   u4_word        : Last fragmented word                  */
/*                                                                          */
/* Globals        :  None                                                   */
/*                                                                          */
/* Processing     :  It has three mode of operations                        */
/*                   1. Byte Swap and Emulation for H.264 WMV9 AP DEC       */
/*                      supports both fragmented and non fragmented packets */
/*                      set u4_prev_0s = last valid zeros for this operation*/
/*                   2. Byte Swap only for MPEG2 and MPEG4 WMV9 MP DEC      */
/*                      supports both fragmented and non fragmented packets */
/*                      set u4_prev_0s = 0 and  u4_0s_bfr_sc = u4_in_len    */
/*                   3. Annex B stream                                      */
/*                      only non fragmented                                 */
/*                      set u4_prev_0s = 0 for this operation               */
/* Outputs        :  pu4_out_len output length of the bit stream            */
/*                                                                          */
/* Returns        :  Number of zeros in case of framented start code        */
/*                                                                          */
/* Known Issues   :                                                         */
/*                                                                          */
/* Revision History                                                         */
/*                                                                          */
/*      DD MM YY            Author        Changes                           */
/*      06 09 2021          Vijay                                           */
/****************************************************************************/
UWORD32 isvcd_nal_byte_swap_emulation(UWORD32 *pu4_out_stream, UWORD32 *pu4_out_len,
                                      UWORD8 *pu1_in_stream, UWORD32 u4_in_len, WORD32 i4_0s_bfr_sc,
                                      void *pv_emulation_ctxt)
{
    UWORD32 u4_i, u4_num_bytes, u4_offset;
    UWORD8 u1_cur_byte;
    emulation_prevent_ctxt_t *ps_emulation_ctxt = (emulation_prevent_ctxt_t *) pv_emulation_ctxt;

    u4_offset = ps_emulation_ctxt->u4_bytes_in_word;
    u4_num_bytes = ps_emulation_ctxt->u4_bytes_in_word;

    for(u4_i = 0; u4_i < u4_in_len; u4_i++)
    {
        UWORD8 u1_cur_byte_emu, u1_cur_byte_sc;
        UWORD64 u8_sft_word;

        u1_cur_byte = *pu1_in_stream++;
        u1_cur_byte_emu = (EMULATION_PREVENTION_BYTE == u1_cur_byte);
        u1_cur_byte_sc = (START_CODE_BYTE == u1_cur_byte);

        if((ps_emulation_ctxt->i4_zeroes_cnt >= i4_0s_bfr_sc) & (u1_cur_byte_emu | u1_cur_byte_sc) &
           (NOT_STUFFED_BYTE == ps_emulation_ctxt->i4_state))
        {
            if(u1_cur_byte_sc)
            {
                break;
            }
            ps_emulation_ctxt->i4_zeroes_cnt = 0;
            ps_emulation_ctxt->i4_state = STUFFED_BYTE;
            continue;
        }

        u8_sft_word = (UWORD64) ps_emulation_ctxt->u4_word << 8;
        ps_emulation_ctxt->u4_word = (UWORD32) (u8_sft_word | u1_cur_byte);
        ps_emulation_ctxt->u4_bytes_in_word++;
        u4_num_bytes++;
        ps_emulation_ctxt->i4_zeroes_cnt++;
        if(u1_cur_byte != 0x00)
        {
            ps_emulation_ctxt->i4_zeroes_cnt = 0;
        }

        if((u4_num_bytes & 0x03) == 0x00)
        {
            *pu4_out_stream = ps_emulation_ctxt->u4_word;
            ps_emulation_ctxt->u4_bytes_in_word = 0;
            pu4_out_stream++;
        }

        ps_emulation_ctxt->i4_state = NOT_STUFFED_BYTE;
    }

    if(ps_emulation_ctxt->u4_bytes_in_word)
    {
        UWORD64 temp_out_stream = (UWORD64) ps_emulation_ctxt->u4_word
                                  << ((4 - ps_emulation_ctxt->u4_bytes_in_word) << 3);
        *pu4_out_stream = (UWORD32) temp_out_stream;
    }

    *pu4_out_len = (u4_num_bytes - u4_offset);
    return ((u4_num_bytes & 0xFFFFFFFC));
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_set_default_nal_header_prms                        */
/*                                                                           */
/*  Description   : Sets the members of NAL header structures to default     */
/*                  values                                                   */
/*                                                                           */
/*  Inputs        : pv_nal_prms - pointer nal header prms structure          */
/*                  i4_temp_id - default value of temporal id                */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : None                                                     */
/*                                                                           */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay      Draft                                    */
/*                                                                           */
/*****************************************************************************/
void isvcd_set_default_nal_prms(void *pv_nal_prms)
{
    nal_prms_t *ps_nal_prms;
    ps_nal_prms = (nal_prms_t *) pv_nal_prms;

    /* Set default values */
    ps_nal_prms->i4_dependency_id = 0;
    ps_nal_prms->i4_derived_nal_type = 0xFF;
    ps_nal_prms->i4_idr_pic_flag = SVCD_FALSE;
    ps_nal_prms->i4_nal_header_len = 0;
    ps_nal_prms->i4_nal_ref_idc = 0xFF;
    ps_nal_prms->i4_nal_unit_type = 0xFF;
    ps_nal_prms->i4_no_int_lyr_pred = 1;
    ps_nal_prms->i4_priority_id = 0;
    ps_nal_prms->i4_quality_id = 0;
    ps_nal_prms->i4_discard_flag = 0;
    ps_nal_prms->i4_dqid = 0;
    ps_nal_prms->i4_use_ref_base_pic_flag = 0;
    ps_nal_prms->i4_temporal_id = 0;
    ps_nal_prms->i4_idr_pic_num = 0;
    ps_nal_prms->u2_frm_num = 0;
    ps_nal_prms->i4_poc_lsb = 0;
    ps_nal_prms->i4_delta_poc_bot = 0;
    ps_nal_prms->ai4_delta_poc[0] = 0;
    ps_nal_prms->ai4_delta_poc[1] = 0;
    ps_nal_prms->u1_pps_id = 0;
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_dec_nal_hdr                                        */
/*                                                                           */
/*  Description   : None                                                     */
/*                                                                           */
/*  Inputs        : pv_buf_ptr - Pointer to buffer constaining start of NAL  */
/*                  pv_nal_header_buf - Temporray working buffer             */
/*                  pv_nal_prms - Pointer to nal header prms                 */
/*                      structure                                            */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : None                                                     */
/*                                                                           */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay      Draft                                    */
/*                                                                           */
/*****************************************************************************/
void isvcd_dec_nal_hdr(void *pv_buf_ptr, WORD32 i4_buf_size, void *pv_nal_header_buf,
                       void *pv_nal_prms, void *pv_prefix_nal_buf, void *pv_prefix_nal_prms,
                       UWORD32 *pu4_err_code)
{
    nal_prms_t *ps_nal_prms;
    nal_prms_t *ps_prefix_nal_prms;
    nal_buf_t *ps_prefix_nal_buf;
    dec_bit_stream_t s_stream_ctxt = {0};
    WORD32 i4_forbidden_zero_bit;

    /* byte swapping */
    UWORD8 *pu1_buf = (UWORD8 *) pv_nal_header_buf;
    UWORD8 *pu1_src = (UWORD8 *) pv_buf_ptr;

    ps_nal_prms = (nal_prms_t *) pv_nal_prms;
    ps_prefix_nal_prms = (nal_prms_t *) pv_prefix_nal_prms;
    ps_prefix_nal_buf = (nal_buf_t *) pv_prefix_nal_buf;

    /* The NAL header syntax elements are read through bitstream fucntions.  */
    /* Hence bitstream context structure initializaton is needed before      */
    /* parsing from the bitstream                                            */
    /* Also bitstream fucntions assume the buffer is byteswapped. Hence the  */
    /* byte swapping is also done for 4 bytes                                */
    s_stream_ctxt.u4_ofst = 0;
    s_stream_ctxt.pu4_buffer = pv_nal_header_buf;
    s_stream_ctxt.u4_max_ofst = (i4_buf_size << 3);

    *pu4_err_code = 0;

    /* Check the size of bitstream buffer */
    if(s_stream_ctxt.u4_max_ofst < 8)
    {
        *pu4_err_code = (UWORD32) NAL_INSUFFICIENT_DATA;
        return;
    }

    if(s_stream_ctxt.u4_max_ofst >= 32)
    {
        *pu1_buf++ = *(pu1_src + 3);
        *pu1_buf++ = *(pu1_src + 2);
        *pu1_buf++ = *(pu1_src + 1);
        *pu1_buf++ = *pu1_src;
    }
    else
    {
        *pu1_buf++ = *pu1_src;
    }

    /*-----------------------------------------------------------------------*/
    /*! Parse the NAL header and update the NAL header structure members     */
    /*-----------------------------------------------------------------------*/
    /* Read forbidden 0 bit */
    i4_forbidden_zero_bit = ih264d_get_bit_h264(&s_stream_ctxt);

    if(0 != i4_forbidden_zero_bit)
    {
        *pu4_err_code = (UWORD32) NAL_CORRUPT_DATA;
        return;
    }

    /*---------------- Read NAL ref idc -----------------------------*/
    ps_nal_prms->i4_nal_ref_idc = ih264d_get_bits_h264(&s_stream_ctxt, 2);

    /*----------------- Read NAL type -------------------------------*/
    ps_nal_prms->i4_nal_unit_type = ih264d_get_bits_h264(&s_stream_ctxt, 5);
    if(ps_nal_prms->i4_nal_unit_type > CODED_SLICE_EXTENSION_NAL)
    {
        *pu4_err_code = (UWORD32) NAL_CORRUPT_DATA;
        return;
    }
    if(ACCESS_UNIT_DELIMITER_RBSP == ps_nal_prms->i4_nal_unit_type)
    {
        ps_nal_prms->i4_derived_nal_type = NON_VCL_NAL;
        return;
    }

    /* set idr pic flag */
    if(IDR_SLICE_NAL == ps_nal_prms->i4_nal_unit_type)
    {
        ps_nal_prms->i4_idr_pic_flag = SVCD_TRUE;
    }
    else
    {
        ps_nal_prms->i4_idr_pic_flag = SVCD_FALSE;
    }

    /*----------------- Read SVC extension NAL header ---------------*/
    if(CODED_SLICE_EXTENSION_NAL == ps_nal_prms->i4_nal_unit_type ||
       PREFIX_UNIT_NAL == ps_nal_prms->i4_nal_unit_type)
    {
        WORD32 i4_svc_extension_flag, i4_idr_flag;

        /* check the size of the buffer */
        if(s_stream_ctxt.u4_max_ofst < 32)
        {
            *pu4_err_code = (UWORD32) NAL_INSUFFICIENT_DATA;
            return;
        }

        i4_svc_extension_flag = ih264d_get_bit_h264(&s_stream_ctxt);
        UNUSED(i4_svc_extension_flag);

        i4_idr_flag = ih264d_get_bit_h264(&s_stream_ctxt);

        /* Set idr pic flag based on idr flag */
        if(1 == i4_idr_flag)
        {
            ps_nal_prms->i4_idr_pic_flag = SVCD_TRUE;
        }
        else
        {
            ps_nal_prms->i4_idr_pic_flag = SVCD_FALSE;
        }

        /* parse priorit id */
        ps_nal_prms->i4_priority_id = ih264d_get_bits_h264(&s_stream_ctxt, 6);

        /* parse the no inter layer prediction flag */
        ps_nal_prms->i4_no_int_lyr_pred = ih264d_get_bit_h264(&s_stream_ctxt);

        /* parse dependency id */
        ps_nal_prms->i4_dependency_id = ih264d_get_bits_h264(&s_stream_ctxt, 3);

        /* parse quality id */
        ps_nal_prms->i4_quality_id = ih264d_get_bits_h264(&s_stream_ctxt, 4);

        if((ps_nal_prms->i4_quality_id > 0) || (ps_nal_prms->i4_dependency_id > 2))
        {
            *pu4_err_code = (UWORD32) NAL_CORRUPT_DATA;
            return;
        }
        /* parse temporal id */
        ps_nal_prms->i4_temporal_id = ih264d_get_bits_h264(&s_stream_ctxt, 3);

        /* parse use ref base pic flag */
        ps_nal_prms->i4_use_ref_base_pic_flag = ih264d_get_bit_h264(&s_stream_ctxt);

        if(0 != ps_nal_prms->i4_use_ref_base_pic_flag)
        {
            *pu4_err_code = (UWORD32) NAL_CORRUPT_DATA;
            return;
        }
        /* parse discrad flag */
        ps_nal_prms->i4_discard_flag = ih264d_get_bit_h264(&s_stream_ctxt);

        /* parse the reserved bits */
        ih264d_get_bits_h264(&s_stream_ctxt, 3);
    }

    /* update NAL hedaer length in bytes */
    ps_nal_prms->i4_nal_header_len = s_stream_ctxt.u4_ofst >> 3;

    /*************************************************************************/
    /* PREFIX NAL UNIT ASSOCIATION WITH ASSOCIATED NAL UNIT                  */
    /*************************************************************************/

    /* if current NAL is not a AVC NAL unit then */
    /* discard the prefix NAL unit if present    */
    if(CODED_SLICE_EXTENSION_NAL == ps_nal_prms->i4_nal_unit_type)
    {
        isvcd_nal_buf_reset(ps_prefix_nal_buf);
    }

    if(SVCD_TRUE == ps_prefix_nal_buf->i4_valid_flag)
    {
        /* Copy the required parameters from the prefix NAL unit */
        ps_nal_prms->i4_dependency_id = ps_prefix_nal_prms->i4_dependency_id;
        ps_nal_prms->i4_quality_id = ps_prefix_nal_prms->i4_quality_id;
        ps_nal_prms->i4_priority_id = ps_prefix_nal_prms->i4_priority_id;
        ps_nal_prms->i4_temporal_id = ps_prefix_nal_prms->i4_temporal_id;
        ps_nal_prms->i4_no_int_lyr_pred = ps_prefix_nal_prms->i4_no_int_lyr_pred;
        ps_nal_prms->i4_use_ref_base_pic_flag = ps_prefix_nal_prms->i4_use_ref_base_pic_flag;
        ps_nal_prms->i4_discard_flag = ps_prefix_nal_prms->i4_discard_flag;
    }

    /*-----------------------------------------------------------------------*/
    /* Set the derived NAL unit type and also update the DQID for VCL NAL    */
    /*  units                                                                */
    /*-----------------------------------------------------------------------*/
    if(CODED_SLICE_EXTENSION_NAL == ps_nal_prms->i4_nal_unit_type ||
       SLICE_NAL == ps_nal_prms->i4_nal_unit_type ||
       IDR_SLICE_NAL == ps_nal_prms->i4_nal_unit_type ||
       PREFIX_UNIT_NAL == ps_nal_prms->i4_nal_unit_type)
    {
        ps_nal_prms->i4_derived_nal_type = VCL_NAL;

        /* calculate the DQID and modified DQID */
        ps_nal_prms->i4_dqid = (ps_nal_prms->i4_dependency_id << 4) + ps_nal_prms->i4_quality_id;
    }
    else
    {
        ps_nal_prms->i4_derived_nal_type = NON_VCL_NAL;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_part_slice_hdr                                */
/*                                                                           */
/*  Description   : This routine parses the slice till POC parameters        */
/*                                                                           */
/*  Inputs        : 1. Pointer to input bitstream                            */
/*                  2. Temporary input buffer                                */
/*                  3. PPS start buffer                                      */
/*                  4. SPS start buffer                                      */
/*                  5. Pointer to NAL paramter structure                     */
/*                  6. Place holder for error code                           */
/*  Globals       : None                                                     */
/*  Processing    : Parses the slice header                                  */
/*                                                                           */
/*  Outputs       : Updated NAL prms structure                               */
/*                  Updated error code                                       */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : Does not support interlaced content                      */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_parse_part_slice_hdr(UWORD8 *pu1_input_buf, WORD32 i4_input_buf_size,
                                  UWORD8 *pu1_temp_buf, void *pv_sps, void *pv_pps,
                                  nal_prms_t *ps_nal_prms, UWORD32 *pu4_err_code,
                                  WORD32 *pi4_sps_pps_status)
{
    UWORD32 u4_slice_type;
    dec_seq_params_t *ps_sps = (dec_seq_params_t *) pv_sps;
    dec_pic_params_t *ps_pps = (dec_pic_params_t *) pv_pps;
    dec_bit_stream_t s_stream_ctxt = {0};
    dec_bit_stream_t *ps_stream_ctxt;
    UWORD32 *pu4_bitstrm_buf;
    UWORD32 *pu4_bitstrm_ofst;

    *pi4_sps_pps_status = NAL_CORRUPT_DATA;
    /* Perform the emulation prevention and byte swap */
    {
        emulation_prevent_ctxt_t s_emulation_ctxt = {0};
        WORD32 i4_size, i4_temp;

        isvcd_reset_emulation_ctxt((void *) &s_emulation_ctxt);
        i4_size = MIN(i4_input_buf_size, HEADER_BUFFER_LEN_BEFORE_EP);

        isvcd_nal_byte_swap_emulation((UWORD32 *) pu1_temp_buf, (UWORD32 *) &i4_temp, pu1_input_buf,
                                      (UWORD32) i4_size, NUM_OF_ZERO_BYTES_BEFORE_START_CODE,
                                      &s_emulation_ctxt);

        /* Initialize the stream context structure */
        s_stream_ctxt.pu4_buffer = (UWORD32 *) pu1_temp_buf;
        s_stream_ctxt.u4_ofst = 0;
        s_stream_ctxt.u4_max_ofst = (i4_size << 3);
    }

    ps_stream_ctxt = &s_stream_ctxt;

    /* Parse the first mb address in slice */
    pu4_bitstrm_buf = ps_stream_ctxt->pu4_buffer;
    pu4_bitstrm_ofst = &ps_stream_ctxt->u4_ofst;
    ps_nal_prms->u4_first_mb_addr = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(ps_nal_prms->u4_first_mb_addr >= (MAX_MBS_LEVEL_51))
    {
        return ERROR_CORRUPTED_SLICE;
    }
    /* Parse slice type */
    u4_slice_type = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

    if(u4_slice_type > 9) return ERROR_INV_SLC_TYPE_T;

    /* Check the validity of slice prms */
    switch(u4_slice_type)
    {
        case 0:
        case 5:
            u4_slice_type = P_SLICE;
            /* P slice */
            break;
        case 1:
        case 6:
            u4_slice_type = B_SLICE;
            /* B slice */
            break;
        case 2:
        case 7:
            /* I slice */
            u4_slice_type = I_SLICE;
            break;
        default:
            break;
    }

    /* Parse the pps id */
    ps_nal_prms->u1_pps_id = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(ps_nal_prms->u1_pps_id & MASK_ERR_PIC_SET_ID) return ERROR_INV_SLICE_HDR_T;

    /* validate pps id */
    ps_pps += ps_nal_prms->u1_pps_id;
    if(0 == ps_pps->u1_is_valid)
    {
        return NOT_OK;
    }
    /* Derive sps id */
    ps_sps = ps_pps->ps_sps;

    ps_nal_prms->u1_sps_id = ps_sps->u1_seq_parameter_set_id;
    if(CODED_SLICE_EXTENSION_NAL == ps_nal_prms->i4_nal_unit_type)
    {
        ps_sps += MAX_NUM_SEQ_PARAMS;
        ps_nal_prms->u1_sps_id = ps_sps->u1_seq_parameter_set_id;
        ps_nal_prms->u1_sps_id += MAX_NUM_SEQ_PARAMS;
    }

    if(NULL == ps_sps)
    {
        return NOT_OK;
    }
    if(FALSE == ps_sps->u1_is_valid)
    {
        return ERROR_INV_SLICE_HDR_T;
    }
    if(ps_nal_prms->u4_first_mb_addr > (ps_sps->u2_frm_ht_in_mbs * ps_sps->u2_frm_wd_in_mbs))
    {
        return ERROR_CORRUPTED_SLICE;
    }
    *pi4_sps_pps_status = 0;

    /* Parse frame number */
    ps_nal_prms->u2_frm_num = ih264d_get_bits_h264(ps_stream_ctxt, ps_sps->u1_bits_in_frm_num);

    /* IDR picture number */
    if(SVCD_TRUE == ps_nal_prms->i4_idr_pic_flag)
    {
        ps_nal_prms->i4_idr_pic_num = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(ps_nal_prms->i4_idr_pic_num > 65535) return ERROR_INV_SLICE_HDR_T;
    }

    /* Poc lsb */
    if(0 == ps_sps->u1_pic_order_cnt_type)
    {
        ps_nal_prms->i4_poc_lsb =
            ih264d_get_bits_h264(ps_stream_ctxt, ps_sps->u1_log2_max_pic_order_cnt_lsb_minus);

        if(ps_nal_prms->i4_poc_lsb < 0 ||
           ps_nal_prms->i4_poc_lsb >= ps_sps->i4_max_pic_order_cntLsb)
            return ERROR_INV_SLICE_HDR_T;
        if(SVCD_TRUE == ps_pps->u1_pic_order_present_flag)
        {
            ps_nal_prms->i4_delta_poc_bot = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }
    }
    else if((1 == ps_sps->u1_pic_order_cnt_type) && (!ps_sps->u1_delta_pic_order_always_zero_flag))
    {
        ps_nal_prms->ai4_delta_poc[0] = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);

        if(SVCD_TRUE == ps_pps->u1_pic_order_present_flag)
        {
            ps_nal_prms->ai4_delta_poc[1] = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        }
    }

    *pu4_err_code = 0;
    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_int_tgt_lyr_attr                                */
/*                                                                           */
/*  Description   : This routine returns the target layer attributes         */
/*                  (dependency id, temporal id and quality id)              */
/*                                                                           */
/*  Inputs        : 1. Application attributes                                */
/*                  2. Internal attributes (input and output)                */
/*                  3. Nal prms structure                                    */
/*  Globals       : None                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : Updated internal target layer attributes                 */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_int_tgt_lyr_attr(target_lyr_attr_t *ps_app_attr, target_lyr_attr_t *ps_int_attr,
                                  nal_prms_t *ps_nal_prms)
{
    WORD32 i4_dep_id;
    WORD32 i4_quality_id;
    WORD32 i4_temp_id;
    WORD32 i4_prior_id;

    /* sanity checks */
    if((NULL == ps_app_attr) || (NULL == ps_int_attr) || (NULL == ps_nal_prms))
    {
        return NOT_OK;
    }

    i4_dep_id = ps_int_attr->i4_dependency_id;
    i4_quality_id = ps_int_attr->i4_quality_id;
    i4_temp_id = ps_int_attr->i4_temporal_id;
    i4_prior_id = ps_int_attr->i4_priority_id;

    /* check for idr pic flag                                  */
    /* dependency & temporal id is updated only for IDR picture */
    if(SVCD_TRUE == ps_nal_prms->i4_idr_pic_flag)
    {
        if(ps_int_attr->i4_dependency_id < ps_app_attr->i4_dependency_id)
        {
            /* update the internal attributes only if             */
            /* current dep_id -1 == highest dep id decoded so far */
            /* and quality id is equal to 0                       */
            if((ps_nal_prms->i4_dependency_id - 1 == ps_int_attr->i4_dependency_id) &&
               (0 == ps_nal_prms->i4_quality_id))
            {
                /* Set revised target dependency id */
                i4_dep_id = ps_nal_prms->i4_dependency_id;
                i4_temp_id = ps_app_attr->i4_temporal_id;
                i4_prior_id = ps_app_attr->i4_priority_id;
            }
        }
        else
        {
            /* cases when the curr dep is greater than or equal to app dep */
            i4_dep_id = ps_app_attr->i4_dependency_id;
            i4_temp_id = ps_app_attr->i4_temporal_id;
            i4_prior_id = ps_app_attr->i4_priority_id;
        }
    }

    /* Set quality id */
    if(i4_dep_id == ps_app_attr->i4_dependency_id)
    {
        i4_quality_id = ps_app_attr->i4_quality_id;
    }
    else
    {
        i4_quality_id = MAX_QUALITY_ID;
    }

    /* Update the internal attributes */
    ps_int_attr->i4_dependency_id = i4_dep_id;
    ps_int_attr->i4_quality_id = i4_quality_id;
    ps_int_attr->i4_temporal_id = i4_temp_id;
    ps_int_attr->i4_priority_id = i4_prior_id;

    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_discard_nal                                         */
/*                                                                           */
/*  Description   : Determines whether current NAL unit has to be discarded  */
/*                  or not                                                   */
/*                                                                           */
/*  Inputs        : pv_nal_prms - Pointer to NAL header prms                 */
/*                      structure                                            */
/*                  pv_app_lyr_attr - Pointer to application target layer    */
/*                      attributes  structure                                */
/*                  pv_app_lyr_attr - Pointer to internal target layer       */
/*                      attributes  structure                                */
/*                  i4_update_flag - This flag indicates whether the internal*/
/*                      target attrbutes should be updated or not            */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : None                                                     */
/*                                                                           */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay           Draft                               */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_discard_nal(void *pv_nal_prms, void *pv_app_attr, void *pv_int_attr,
                         WORD32 i4_update_flag)
{
    WORD32 i4_discard_nal_flag;
    nal_prms_t *ps_nal_prms;
    target_lyr_attr_t *ps_app_attr;
    target_lyr_attr_t *ps_int_attr;
    WORD32 i4_status;

    ps_nal_prms = (nal_prms_t *) pv_nal_prms;
    ps_app_attr = (target_lyr_attr_t *) pv_app_attr;
    ps_int_attr = (target_lyr_attr_t *) pv_int_attr;

    /* Get the updated target layer attributes */
    if(SVCD_TRUE == i4_update_flag)
    {
        i4_status = isvcd_get_int_tgt_lyr_attr(ps_app_attr, ps_int_attr, ps_nal_prms);
        if(OK != i4_status)
        {
            return NOT_OK;
        }
    }

    i4_discard_nal_flag = SVCD_FALSE;

    if(VCL_NAL == ps_nal_prms->i4_derived_nal_type)
    {
        /*-------------------------------------------------------------------*/
        /*!Discard VCL NAL if any of following is true                       */
        /*! - Dependency id is greater than target dependency id             */
        /*! - Dependency id is equal to target dependency id but quality id  */
        /*!   is greater than target quality id                              */
        /*! - priority id is greater than target priority id                 */
        /*! - Temporal id is greater than target temporal id                 */
        /*! - If dependency id is greater than a NAL unit for which discard  */
        /*!   flag of the NAL header is set                                  */
        /*-------------------------------------------------------------------*/
        if(PREFIX_UNIT_NAL != ps_nal_prms->i4_nal_unit_type)
        {
            if(ps_nal_prms->i4_dependency_id > ps_int_attr->i4_dependency_id)
            {
                i4_discard_nal_flag = SVCD_TRUE;
            }

            if(ps_nal_prms->i4_dependency_id == ps_int_attr->i4_dependency_id &&
               ps_nal_prms->i4_quality_id > ps_int_attr->i4_quality_id)
            {
                i4_discard_nal_flag = SVCD_TRUE;
            }

            if(ps_nal_prms->i4_temporal_id > ps_int_attr->i4_temporal_id)
            {
                i4_discard_nal_flag = SVCD_TRUE;
            }

            if(ps_nal_prms->i4_priority_id > ps_int_attr->i4_priority_id)
            {
                i4_discard_nal_flag = SVCD_TRUE;
            }
        }
        else
        {
            if(0 == ps_int_attr->i4_quality_id && 0 == ps_int_attr->i4_dependency_id)
            {
                i4_discard_nal_flag = SVCD_TRUE;
            }
        }
    }

    return (i4_discard_nal_flag);
}
