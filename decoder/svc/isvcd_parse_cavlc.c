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
 *  isvcd_parse_cavlc.c
 *
 * @brief
 *  This file contains UVLC related functions.
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_parse_bmb_ref_index_cavlc_range1()
 *  - isvcd_parse_bmb_ref_index_cavlc()
 *  - isvcd_parse_pmb_ref_index_cavlc()
 *  - isvcd_parse_pmb_ref_index_cavlc_range1()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <string.h>
#include <stdio.h>

#include "ih264d_bitstrm.h"
#include "isvcd_parse_cavlc.h"
#include "ih264d_error_handler.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "ih264d_cabac.h"
#include "isvcd_structs.h"
#include "ih264d_tables.h"
#include "ih264d_tables.h"
#include "ih264d_mb_utils.h"
#include "ih264d_parse_cavlc.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_bmb_ref_index_cavlc_range1                   */
/*                                                                           */
/*  Description   : This function does the Cavlc  TEV range > 1 parsing of   */
/*                  reference index  for a B MB.                             */
/*                  Range > 1 when num_ref_idx_active_minus1 > 0             */
/*                                                                           */
/*  Inputs        : <What inputs does the function take?>                    */
/*  Globals       : <Does it use any global variables?>                      */
/*  Processing    : <Describe how the function operates - include algorithm  */
/*                  description>                                             */
/*  Outputs       : <What does the function produce?>                        */
/*  Returns       : <What does the function return?>                         */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         19 09 2008   Jay          Draft                                   */
/*                                                                           */
/*****************************************************************************/

void isvcd_parse_bmb_ref_index_cavlc_range1(
    UWORD32 u4_num_part,                  /* Number of partitions in MB      */
    dec_bit_stream_t *ps_bitstrm,         /* Pointer to bitstream Structure. */
    WORD8 *pi1_ref_idx,                   /* pointer to reference index array */
    UWORD32 u4_num_ref_idx_active_minus1, /* Not used for range 1    */
    UWORD8 *pu1_motion_prediction_flag    /*motion_pred_flag */

)
{
    UWORD32 u4_i;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstream_off = &ps_bitstrm->u4_ofst;
    UNUSED(u4_num_ref_idx_active_minus1);
    for(u4_i = 0; u4_i < u4_num_part; u4_i++)
    {
        if(pi1_ref_idx[u4_i] > -1 && (((*pu1_motion_prediction_flag >> u4_i) & 0x01) == 0))
        {
            UWORD32 u4_ref_idx;
            u4_ref_idx = ih264d_tev_range1(pu4_bitstream_off, pu4_bitstrm_buf);

            /* Storing Reference Idx Information */
            pi1_ref_idx[u4_i] = (WORD8) u4_ref_idx;
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_bmb_ref_index_cavlc                          */
/*                                                                           */
/*  Description   : This function does the Cavlc  TEV range > 1 parsing of   */
/*                  reference index  for a B MB.                             */
/*                  Range > 1 when num_ref_idx_active_minus1 > 0             */
/*                                                                           */
/*  Inputs        : <What inputs does the function take?>                    */
/*  Globals       : <Does it use any global variables?>                      */
/*  Processing    : <Describe how the function operates - include algorithm  */
/*                  description>                                             */
/*  Outputs       : <What does the function produce?>                        */
/*  Returns       : <What does the function return?>                         */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         19 09 2008   Jay          Draft                                   */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_parse_bmb_ref_index_cavlc(
    UWORD32 u4_num_part,                  /* Number of partitions in MB      */
    dec_bit_stream_t *ps_bitstrm,         /* Pointer to bitstream Structure. */
    WORD8 *pi1_ref_idx,                   /* pointer to reference index array */
    UWORD32 u4_num_ref_idx_active_minus1, /* Number of active references - 1  */
    UWORD8 *pu1_motion_prediction_flag    /*motion_pred_flag */
)
{
    UWORD32 u4_i;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstream_off = &ps_bitstrm->u4_ofst;

    for(u4_i = 0; u4_i < u4_num_part; u4_i++)
    {
        if(pi1_ref_idx[u4_i] > -1 && (((*pu1_motion_prediction_flag >> u4_i) & 0x01) == 0))
        {
            UWORD32 u4_ref_idx;
            // inlining ih264d_uev
            UWORD32 u4_bitstream_offset = *pu4_bitstream_off;
            UWORD32 u4_word, u4_ldz;

            /***************************************************************/
            /* Find leading zeros in next 32 bits                          */
            /***************************************************************/
            NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
            u4_ldz = CLZ(u4_word);
            /* Flush the ps_bitstrm */
            u4_bitstream_offset += (u4_ldz + 1);
            /* Read the suffix from the ps_bitstrm */
            u4_word = 0;
            if(u4_ldz) GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
            *pu4_bitstream_off = u4_bitstream_offset;
            u4_ref_idx = ((1 << u4_ldz) + u4_word - 1);
            // inlining ih264d_uev
            if(u4_ref_idx > u4_num_ref_idx_active_minus1) return ERROR_REF_IDX;

            /* Storing Reference Idx Information */
            pi1_ref_idx[u4_i] = (WORD8) u4_ref_idx;
        }
    }
    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_pmb_ref_index_cavlc                          */
/*                                                                           */
/*  Description   : This function does the Cavlc  TEV range > 1 parsing of   */
/*                  reference index  for a P MB.                             */
/*                  Range > 1 when num_ref_idx_active_minus1 > 0             */
/*                                                                           */
/*  Inputs        : <What inputs does the function take?>                    */
/*  Globals       : <Does it use any global variables?>                      */
/*  Processing    : <Describe how the function operates - include algorithm  */
/*                  description>                                             */
/*  Outputs       : <What does the function produce?>                        */
/*  Returns       : <What does the function return?>                         */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         19 09 2008   Jay          Draft                                   */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_parse_pmb_ref_index_cavlc(
    UWORD32 u4_num_part,                  /* Number of partitions in MB      */
    dec_bit_stream_t *ps_bitstrm,         /* Pointer to bitstream Structure. */
    WORD8 *pi1_ref_idx,                   /* pointer to reference index array */
    UWORD32 u4_num_ref_idx_active_minus1, /* Number of active references - 1  */
    UWORD8 *pu1_motion_prediction_flag    /*motion_pred_flag_l0 */
)
{
    UWORD32 u4_i;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstream_off = &ps_bitstrm->u4_ofst;

    for(u4_i = 0; u4_i < u4_num_part; u4_i++)
    {
        if(((*pu1_motion_prediction_flag >> u4_i) & 0x01) == 0)
        {
            UWORD32 u4_ref_idx;
            // Inlined ih264d_uev
            UWORD32 u4_bitstream_offset = *pu4_bitstream_off;
            UWORD32 u4_word, u4_ldz;

            /***************************************************************/
            /* Find leading zeros in next 32 bits                          */
            /***************************************************************/
            NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
            u4_ldz = CLZ(u4_word);
            /* Flush the ps_bitstrm */
            u4_bitstream_offset += (u4_ldz + 1);
            /* Read the suffix from the ps_bitstrm */
            u4_word = 0;
            if(u4_ldz) GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
            *pu4_bitstream_off = u4_bitstream_offset;
            u4_ref_idx = ((1 << u4_ldz) + u4_word - 1);

            // Inlined ih264d_uev
            if(u4_ref_idx > u4_num_ref_idx_active_minus1) return ERROR_REF_IDX;

            /* Storing Reference Idx Information */
            pi1_ref_idx[u4_i] = (WORD8) u4_ref_idx;
        }
    }
    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_pmb_ref_index_cavlc_range1                   */
/*                                                                           */
/*  Description   : This function does the Cavlc  TEV range =1 parsing of    */
/*                  reference index  for a P MB. Range is 1 when             */
/*                  num_ref_idx_active_minus1 is 0                           */
/*                                                                           */
/*  Inputs        : <What inputs does the function take?>                    */
/*  Globals       : <Does it use any global variables?>                      */
/*  Processing    : <Describe how the function operates - include algorithm  */
/*                  description>                                             */
/*  Outputs       : <What does the function produce?>                        */
/*  Returns       : <What does the function return?>                         */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         19 09 2008   Jay          Draft                                   */
/*                                                                           */
/*****************************************************************************/
void isvcd_parse_pmb_ref_index_cavlc_range1(
    UWORD32 u4_num_part,                  /* Number of partitions in MB      */
    dec_bit_stream_t *ps_bitstrm,         /* Pointer to bitstream Structure. */
    WORD8 *pi1_ref_idx,                   /* pointer to reference index array */
    UWORD32 u4_num_ref_idx_active_minus1, /* Not used for range 1    */
    UWORD8 *pu1_motion_prediction_flag    /*motion_pred_flag_l0 */
)
{
    UWORD32 u4_i;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstream_off = &ps_bitstrm->u4_ofst;
    UNUSED(u4_num_ref_idx_active_minus1);
    for(u4_i = 0; u4_i < u4_num_part; u4_i++)
    {
        if(((*pu1_motion_prediction_flag >> u4_i) & 0x01) == 0)
        {
            UWORD32 u4_ref_idx;
            u4_ref_idx = ih264d_tev_range1(pu4_bitstream_off, pu4_bitstrm_buf);

            /* Storing Reference Idx Information */
            pi1_ref_idx[u4_i] = (WORD8) u4_ref_idx;
        }
    }
}
