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
 * \file isvcd_nal_parse.c
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
#include "ih264d_debug.h"
#include "ih264_debug.h"
#include "ih264d_inter_pred.h"
#include "isvcd_structs.h"
#include "ih264d_nal.h"
#include "ih264d_error_handler.h"

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
/*  Function Name : isvcd_get_nal_buf                                         */
/*                                                                           */
/*  Description   : This routine returns the NAL buffer structure to use for */
/*                  current NAL unit. This will also perform the initializa -*/
/*                  tion of the structure                                    */
/*  Inputs        : 1. NAL parse structure                                   */
/*                  2. Place holder for nal buffer structure                 */
/*  Globals       : None                                                     */
/*  Processing    : If current NAL unit prefix NAL unit then                 */
/*                      - Resets the prefix nal buffer structure             */
/*                      - Assigns the buffer pointer                         */
/*                  Otherwise                                                */
/*                      - Assigns the buffer pointer                         */
/*  Outputs       :  - Updated NAL buffer strcuture                          */
/*                   - Updates the place holder with correct NAL buffer      */
/*                  structure                                                */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/
void isvcd_get_nal_buf(nal_parse_ctxt_t *ps_nal_parse_ctxt, nal_buf_t **pps_nal_buf)
{
    nal_prms_t *ps_nal_prms;
    nal_buf_t *ps_nal_buf;

    ps_nal_prms = &ps_nal_parse_ctxt->s_nal_prms;

    /* Get the NAL buffer structure */
    if(PREFIX_UNIT_NAL == ps_nal_prms->i4_nal_unit_type)
    {
        ps_nal_buf = &ps_nal_parse_ctxt->s_prefix_nal_buf;

        /* Note: This reset will cause a prefix NAL unit */
        /* which is followed by another prefix NAL unit  */
        /* to be ignored by the module. This is indeed   */
        /* a desired behaviour                           */
        isvcd_nal_buf_reset(ps_nal_buf);
    }
    else
    {
        ps_nal_buf = &ps_nal_parse_ctxt->s_nal_buf;
    }

    /* Initialize the buffer structure */
    ps_nal_buf->i4_valid_flag = SVCD_TRUE;
    if(VCL_NAL == ps_nal_prms->i4_derived_nal_type)
    {
        ps_nal_buf->pu1_buf = ps_nal_parse_ctxt->pu1_vcl_nal_buf;
    }
    else if(NON_VCL_NAL == ps_nal_prms->i4_derived_nal_type)
    {
        ps_nal_buf->pu1_buf = ps_nal_parse_ctxt->pu1_non_vcl_nal_buf;
    }
    else
    {
        ps_nal_buf->pu1_buf = NULL;
        return;
    }

    *pps_nal_buf = ps_nal_buf;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_dqid_ctxt_reset                                     */
/*                                                                           */
/*  Description   : This routine resets the DQID context. This routine shall */
/*                  be invoked once in a picture                             */
/*                                                                           */
/*  Inputs        : DQID context structure                                   */
/*  Globals       : None                                                     */
/*  Processing    : Invalidates all the DQID nodes                           */
/*                                                                           */
/*  Outputs       : Updated DQID context                                     */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_dqid_ctxt_reset(dqid_ctxt_t *ps_dqid_ctxt)
{
    WORD32 i4_lyr_idx;
    WORD32 i4_max_num_lyrs;
    dqid_node_t *ps_dqid_node;

    /* sanity checks */
    if(NULL == ps_dqid_ctxt)
    {
        return NOT_OK;
    }

    i4_max_num_lyrs = ps_dqid_ctxt->i4_max_num_lyrs;
    ps_dqid_node = ps_dqid_ctxt->ps_dqid_node;

    /* Loop over all the layers */
    for(i4_lyr_idx = 0; i4_lyr_idx < i4_max_num_lyrs; i4_lyr_idx++)
    {
        /* Reset the valid flag */
        ps_dqid_node->u1_valid_flag = SVCD_FALSE;

        /* Loop updates */
        ps_dqid_node += 1;
    } /* loop over all the layers */

    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_dqid_node                                       */
/*                                                                           */
/*  Description   : This routine gets a DQID node corresponding to a DQID    */
/*                                                                           */
/*  Inputs        : 1. DQID context                                          */
/*                  2. DQID                                                  */
/*                  3. Place holder for DQID node (output)                   */
/*  Globals       : None                                                     */
/*  Processing    : It performs the following                                */
/*                  - Searches for all elements untill it gets element having*/
/*                    DQID equal to input DQID.                              */
/*                  - If not found it finds a free (un-occupied) node        */
/*                                                                           */
/*  Outputs       : 1. Updated DQID node                                     */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_get_dqid_node(dqid_ctxt_t *ps_dqid_ctxt, UWORD8 u1_dqid, dqid_node_t **pps_dqid_node)
{
    WORD32 i4_lyr_idx;
    WORD32 i4_max_num_lyrs;
    dqid_node_t *ps_dqid_node;
    dqid_node_t *ps_rqrd_dqid_node;

    /* sanity checks */
    if((NULL == ps_dqid_ctxt) || (NULL == pps_dqid_node))
    {
        return NOT_OK;
    }

    i4_max_num_lyrs = ps_dqid_ctxt->i4_max_num_lyrs;
    ps_dqid_node = ps_dqid_ctxt->ps_dqid_node;

    /*Initialization */
    ps_rqrd_dqid_node = NULL;

    /* Loop over all the buffer nodes */
    for(i4_lyr_idx = 0; i4_lyr_idx < i4_max_num_lyrs; i4_lyr_idx++)
    {
        if((SVCD_TRUE == ps_dqid_node->u1_valid_flag) && (u1_dqid == ps_dqid_node->u1_dqid))
        {
            ps_rqrd_dqid_node = ps_dqid_node;
            break;
        }
        /* Loop updates */
        ps_dqid_node += 1;
    } /* Loop over all the buffer nodes */

    if(NULL == ps_rqrd_dqid_node)
    {
        /* If vcl node is not allocated for the requested DQID then allocate buffer */
        ps_dqid_node = ps_dqid_ctxt->ps_dqid_node;
        for(i4_lyr_idx = 0; i4_lyr_idx < i4_max_num_lyrs; i4_lyr_idx++)
        {
            if(SVCD_FALSE == ps_dqid_node->u1_valid_flag)
            {
                break;
            }
            /* Loop updates */
            ps_dqid_node += 1;
        } /* Loop over all the nodes */
        /* Update the node structure */
        ps_rqrd_dqid_node = ps_dqid_node;
    }

    /* sanity checks */
    if(NULL == ps_rqrd_dqid_node)
    {
        return NOT_OK;
    }
    *pps_dqid_node = ps_rqrd_dqid_node;

    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_reset_ctxt                                     */
/*                                                                           */
/*  Description   : This routine performs NAL unit level initialization      */
/*                  This routine shall be called before parsing a NAL unit   */
/*                                                                           */
/*  Inputs        : 1. NAL parse context structure                           */
/*  Globals       : None                                                     */
/*  Processing    : This does initializaiton of NAL unit level tracking      */
/*                  varaibles                                                */
/*                                                                           */
/*  Outputs       : Updated context structure                                */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_nal_reset_ctxt(nal_parse_ctxt_t *ps_nal_parse_ctxt)
{
    nal_unit_t *ps_nal_unit;

    if(NULL == ps_nal_parse_ctxt)
    {
        return NOT_OK;
    }

    /* Reset the NAL boundary detetction */
    ps_nal_parse_ctxt->i4_find_nal_state = NAL_START;
    ps_nal_parse_ctxt->i4_zero_byte_cnt = 0;
    ps_nal_unit = ps_nal_parse_ctxt->pv_nal_unit;
    ps_nal_unit->i4_num_bufs = 0;

    /*Reset emulation prevention */
    isvcd_reset_emulation_ctxt(&ps_nal_parse_ctxt->s_emulation_ctxt);

    /*Reset the NAL header prms */
    isvcd_set_default_nal_prms(&ps_nal_parse_ctxt->s_nal_prms);

    /* Reset other NAL level tracking variables */
    ps_nal_parse_ctxt->i4_discard_nal_flag = SVCD_FALSE;

    /*Reset NAL buffer structure*/
    isvcd_nal_buf_reset(&ps_nal_parse_ctxt->s_nal_buf);

    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_pic_reset_ctxt                                      */
/*                                                                           */
/*  Description   : This routine performs the picture level initialization.  */
/*                  This routine shall be called before parsing a access unit*/
/*                                                                           */
/*  Inputs        : pv_nal_parse_ctxt - Pointer to context structure         */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : 1. Resets the varaibles                                  */
/*                                                                           */
/*  Outputs       : Updated context structure                                */
/*                                                                           */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay           Draft                               */
/*                                                                           */
/*****************************************************************************/
void isvcd_pic_reset_ctxt(nal_parse_ctxt_t *ps_nal_parse_ctxt)
{
    WORD32 i4_status;

    /*-----------------------------------------------------------------------*/
    /*! Reset NAL boundary detetction logic                                  */
    /*-----------------------------------------------------------------------*/
    i4_status = isvcd_nal_reset_ctxt(ps_nal_parse_ctxt);

    UNUSED(i4_status);

    /*-----------------------------------------------------------------------*/
    /*! Reset picture boundary detctetion logic                              */
    /*-----------------------------------------------------------------------*/
    ps_nal_parse_ctxt->i4_is_frst_vcl_nal_in_au = SVCD_TRUE;

    /*-----------------------------------------------------------------------*/
    /*! Reset VCL and non VCL NAL buffer tracking variables                  */
    /*-----------------------------------------------------------------------*/
    ps_nal_parse_ctxt->pu1_non_vcl_nal_buf = ps_nal_parse_ctxt->pv_non_vcl_nal_buf;
    ps_nal_parse_ctxt->pu1_vcl_nal_buf = ps_nal_parse_ctxt->pv_vcl_nal_buf;

    /* reset the bytes left to buffer size */
    ps_nal_parse_ctxt->u4_bytes_left_vcl = MAX_VCL_NAL_BUFF_SIZE;

    /* 85% of the buffer is used. 15% is used to handle error cases*/
    ps_nal_parse_ctxt->u4_bytes_left_non_vcl = (MAX_NON_VCL_NAL_BUFF_SIZE * 0.85);

    /* Offset the buffer to start of vcl data */
    UPDATE_NAL_BUF_PTR(&ps_nal_parse_ctxt->pu1_non_vcl_nal_buf, NON_VCL_NAL,
                       &ps_nal_parse_ctxt->u4_bytes_left_non_vcl);

    UPDATE_NAL_BUF_PTR(&ps_nal_parse_ctxt->pu1_vcl_nal_buf, VCL_NAL,
                       &ps_nal_parse_ctxt->u4_bytes_left_vcl);

    /* Reset previous field */
    ps_nal_parse_ctxt->ps_prev_non_vcl_buf = NULL;
    ps_nal_parse_ctxt->i4_idr_pic_err_flag = 0;

    /*-----------------------------------------------------------------------*/
    /*! Reset other NAL related tracking variables                           */
    /*-----------------------------------------------------------------------*/
    ps_nal_parse_ctxt->i4_num_non_vcl_nals = 0;

    /* Reset the vcl nal node buffer context */
    i4_status = isvcd_dqid_ctxt_reset(&ps_nal_parse_ctxt->s_dqid_ctxt);

    /* Reset target layer update flag */
    ps_nal_parse_ctxt->i4_tgt_lyr_update = SVCD_TRUE;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_get_nal_prms                                        */
/*                                                                           */
/*  Description   : This routine will update the nal prms                    */
/*  Inputs        : 1. Start of bitstream buffer containing NAL header       */
/*                  2. Size of the buffer                                    */
/*                  3. NAL prms structure                                    */
/*                  4. Place holder for error code                           */
/*                  5. Place holder for nal discard flag                     */
/*                  6. NAL parse context structure                           */
/*  Globals       : None                                                     */
/*  Processing    : 1. Parses the NAL header                                 */
/*                  2. Sets the discard flag                                 */
/*                  3. If NAL is not discarded and nal is VCL NAL unit then  */
/*                     decodes the slice prms (prefix nal units are excluded)*/
/*  Outputs       : Updated NAL prms structure                               */
/*                  Updated NAL discrd flag                                  */
/*                  Updates the error code if encountered with error         */
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
WORD32 isvcd_get_nal_prms(UWORD8 *pu1_buf, WORD32 i4_buf_size, nal_prms_t *ps_nal_prms,
                          nal_prms_t *ps_prefix_nal_prms, nal_buf_t *ps_prefix_nal_buf,
                          UWORD32 *pu4_err_code, WORD32 *pi4_sps_pps_status,
                          WORD32 *pi4_nal_discard_flag, nal_parse_ctxt_t *ps_nal_parse_ctxt)
{
    UWORD8 *pu1_input_buf;
    WORD32 i4_status;
    dec_seq_params_t *ps_sps;
    dec_pic_params_t *ps_pps;

    ps_sps = ps_nal_parse_ctxt->pv_seq_prms;
    ps_pps = ps_nal_parse_ctxt->pv_pic_prms;

    *pu4_err_code = 0;
    *pi4_sps_pps_status = NAL_CORRUPT_DATA;

    /* Decode the NAL header */
    isvcd_dec_nal_hdr(pu1_buf, i4_buf_size, ps_nal_parse_ctxt->pv_nal_header_buf, ps_nal_prms,
                      ps_prefix_nal_buf, ps_prefix_nal_prms, pu4_err_code);

    /* If encountered with error return fail */
    if(0 != *pu4_err_code)
    {
        return (NOT_OK);
    }

    if(ACCESS_UNIT_DELIMITER_RBSP == ps_nal_prms->i4_nal_unit_type)
    {
        *pi4_nal_discard_flag = 1;
        return OK;
    }

    /* Set the discard flag */
    *pi4_nal_discard_flag = isvcd_discard_nal(
        (void *) ps_nal_prms, (void *) &ps_nal_parse_ctxt->s_app_attr,
        (void *) &ps_nal_parse_ctxt->s_int_attr, ps_nal_parse_ctxt->i4_tgt_lyr_update);

    /* Parse the slice header if all the following */
    /* conditions are true                         */
    /* 1. NAL is a VCL NAL unit                    */
    /* 2. NAL is not a prefix NAL unit             */
    /* 3. NAL is not discarded                     */
    if((NON_VCL_NAL == ps_nal_prms->i4_derived_nal_type) ||
       (PREFIX_UNIT_NAL == ps_nal_prms->i4_nal_unit_type) || (SVCD_TRUE == *pi4_nal_discard_flag))
    {
        return (OK);
    }

    pu1_input_buf = pu1_buf;
    pu1_input_buf += ps_nal_prms->i4_nal_header_len;
    i4_buf_size -= ps_nal_prms->i4_nal_header_len;

    i4_status =
        isvcd_parse_part_slice_hdr(pu1_input_buf, i4_buf_size, ps_nal_parse_ctxt->pv_nal_header_buf,
                                   ps_sps, ps_pps, ps_nal_prms, pu4_err_code, pi4_sps_pps_status);

    return (i4_status);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_compare_nal_prms                                    */
/*                                                                           */
/*  Description   : Detects the picture boundary for annex B based input     */
/*                  bitstream                                                */
/*                                                                           */
/*  Inputs        : 1. Pointer to NAL prms                                   */
/*                  2. Pass (first pass or second pass (verification)        */
/*                  3. Place holder for picture boundary type                */
/*                  4. Place holder for picture boundary status              */
/*                  4. pointer to bitstream extract context structure        */
/*  Globals       :                                                          */
/*  Processing    : Detects the picture bounadry as described in G.7.4.1.2.4 */
/*                                                                           */
/*  Outputs       : Detects the picture boundary                             */
/*                  Updates the first NAL in AU field                        */
/*                  Updates the picture boundary type if picture boundary is */
/*                      detetcetd otherwise it's value shall be ignored      */
/*                  Updates the picture boundary status with either          */
/*                  PIC_BOUNDARY_TRUE if picture boundary is detetcted or    */
/*                  PIC_BOUNDARY_FALSE otherwise                             */
/*                  Updates the error code                                   */
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
WORD32 isvcd_compare_nal_prms(nal_prms_t *ps_nal_prms, WORD32 i4_pass, WORD32 i4_prev_dqid,
                              WORD32 *pi4_pic_bound_type, WORD32 *pi4_pic_bound_status,
                              nal_parse_ctxt_t *ps_nal_parse_ctxt)
{
    dqid_node_t *ps_dqid_node;
    vcl_node_t *ps_vcl_node;
    WORD32 i4_status;

    /* If DQID is lesser than the DQID of the previous */
    /* NAL then declare the picture boundary           */
    *pi4_pic_bound_type = PIC_BOUND_DQID;
    if(i4_prev_dqid > ps_nal_prms->i4_dqid)
    {
        *pi4_pic_bound_status = PIC_BOUNDARY_TRUE;
        return (OK);
    }

    /* Perform the picture boundary detection only for */
    /* the layers with quality id equal to 0           */
    if((FIRST_PASS == i4_pass) && (0 != (ps_nal_prms->i4_dqid & 0x0F)))
    {
        *pi4_pic_bound_status = PIC_BOUNDARY_FALSE;
        return (OK);
    }

    /* Get the DQID node */
    i4_status =
        isvcd_get_dqid_node(&ps_nal_parse_ctxt->s_dqid_ctxt, (UWORD8) i4_prev_dqid, &ps_dqid_node);
    if((OK != i4_status) || (NULL == ps_dqid_node))
    {
        return NOT_OK;
    }
    /* If the current slice is first slice in the layer */
    /* then do not compare                              */
    if(SVCD_FALSE == ps_dqid_node->u1_valid_flag)
    {
        *pi4_pic_bound_status = PIC_BOUNDARY_FALSE;
        return (OK);
    }

    *pi4_pic_bound_type = PIC_BOUND_SLICE_PRMS;
    *pi4_pic_bound_status = PIC_BOUNDARY_TRUE;
    ps_vcl_node = ps_dqid_node->ps_vcl_node;

    /* Compare NAL ref idc */
    {
        WORD32 i4_prev_ref_pic_flag;
        WORD32 i4_cur_ref_pic_flag;

        i4_prev_ref_pic_flag = (0 != ps_vcl_node->i4_nal_ref_idc);
        i4_cur_ref_pic_flag = (0 != ps_nal_prms->i4_nal_ref_idc);

        if(i4_prev_ref_pic_flag != i4_cur_ref_pic_flag)
        {
            return (OK);
        }
    }

    /* Compare IDR picture flag */
    if(ps_vcl_node->i4_idr_pic_flag != ps_nal_prms->i4_idr_pic_flag)
    {
        return (OK);
    }

    /* Compare PPS id */
    if(ps_vcl_node->u1_pps_id != ps_nal_prms->u1_pps_id)
    {
        return (OK);
    }

    /* Compare idr pic num */
    if((SVCD_TRUE == ps_nal_prms->i4_idr_pic_flag) &&
       (ps_vcl_node->i4_idr_pic_num != ps_nal_prms->i4_idr_pic_num))
    {
        return (OK);
    }

    /* Compare frame number */
    if(ps_vcl_node->u2_frm_num != ps_nal_prms->u2_frm_num)
    {
        return (OK);
    }

    /* Compare poc lsb */
    if(ps_dqid_node->i4_poc_lsb != ps_nal_prms->i4_poc_lsb)
    {
        return (OK);
    }

    /* Compare delta poc bottom */
    if(ps_dqid_node->i4_delta_poc_bot != ps_nal_prms->i4_delta_poc_bot)
    {
        return (OK);
    }

    /* Compare delta poc [0] */
    if(ps_dqid_node->ai4_delta_poc[0] != ps_nal_prms->ai4_delta_poc[0])
    {
        return (OK);
    }

    /* Compare delta poc [0] */
    if(ps_dqid_node->ai4_delta_poc[1] != ps_nal_prms->ai4_delta_poc[1])
    {
        return (OK);
    }

    *pi4_pic_bound_status = PIC_BOUNDARY_FALSE;
    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_detetct_pic_boundary_annex_b                        */
/*                                                                           */
/*  Description   : Detects the picture boundary for annex B based input     */
/*                  bitstream                                                */
/*                                                                           */
/*                                                                           */
/*  Inputs        : 1. Pointer to NAL prms                                   */
/*                  2. Input bitstream structure                             */
/*                  3. Current position of the bitstream pointer             */
/*                  4. Place holder for picture boundary status              */
/*                  5. pointer to bitstream extract context structure        */
/*  Globals       :                                                          */
/*  Processing    : It does the following                                    */
/*                  1. Look for next NAL.                                    */
/*                      If not found then declare picture boundary           */
/*                      Otherwsie goto next step                             */
/*                  2. Parse the NAL header                                  */
/*                      If encountered with error then declare picture       */
/*                      boundary                                             */
/*                      Otherwise goto next step                             */
/*                  3. If picture boundary type is                           */
/*                      DQID change and DQID is not equal previous DQID then */
/*                      declare picture boundary. Otherwise, the comapre the */
/*                      rest of parameters. If during comparison, if there is*/
/*                  4. If picture boundary type is                           */
/*                      SLICE PRMS CHANGE and Dependency id is not equal then*/
/*                      declare picture boundary. Otherwise compre rest of   */
/*                      parameters and goto step 5                           */
/*                  5. If during comparison, if there is                     */
/*                       * an error - then declare picture boundary          */
/*                       * Otherwsie if picture  boundary is not detetcted   */
/*                         then discard the second slice and proceed.        */
/*                                                                           */
/*  Outputs       : Detects the picture boundary                             */
/*                  Updates the first NAL in AU field                        */
/*                  Updates the picture boundary type if picture boundary is */
/*                      detetcetd otherwise it's value shall be ignored      */
/*                  Updates the error code                                   */
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
WORD32 isvcd_detect_pic_boundary_annex_b(nal_prms_t *ps_nal_prms, UWORD8 *pu1_stream_buffer,
                                         WORD32 i4_cur_pos, WORD32 *pi4_pic_bound_status,
                                         nal_parse_ctxt_t *ps_nal_parse_ctxt,
                                         UWORD32 *pu4_num_bytes)
{
    UWORD32 u4_err_code;
    WORD32 i4_zero_cnt;
    WORD32 i4_status;
    nal_prms_t s_nal_prms = {0};
    nal_prms_t s_prefix_nal_prms = {0};
    nal_buf_t s_prefix_nal_buf = {0};
    WORD32 i4_pic_bound_type;
    WORD32 i4_pic_bound_status;
    UWORD8 *pu1_buf;
    WORD32 i4_buf_size;
    WORD32 i4_more_data_flag;
    WORD32 i4_new_lyr_flag;
    WORD32 i4_prev_dqid;
    WORD32 i4_nal_discard_flag;

    /* Initializations */
    i4_zero_cnt = 0;
    s_prefix_nal_buf.i4_valid_flag = SVCD_FALSE;
    *pi4_pic_bound_status = PIC_BOUNDARY_FALSE;
    i4_new_lyr_flag = SVCD_TRUE;

    /* Get the previous layer's DQID                    */
    if(SVCD_TRUE == ps_nal_parse_ctxt->i4_is_frst_vcl_nal_in_au)
    {
        ps_nal_parse_ctxt->i4_prev_dq_id = ps_nal_prms->i4_dqid;
        ps_nal_parse_ctxt->i4_is_frst_vcl_nal_in_au = SVCD_FALSE;
    }
    i4_prev_dqid = ps_nal_parse_ctxt->i4_prev_dq_id;
    ps_nal_parse_ctxt->i4_prev_dq_id = ps_nal_prms->i4_dqid;

    /* Detect the picture boundary */
    if(ps_nal_prms->i4_dqid <= i4_prev_dqid)
    {
        i4_status =
            isvcd_compare_nal_prms(ps_nal_prms, FIRST_PASS, i4_prev_dqid, &i4_pic_bound_type,
                                   &i4_pic_bound_status, ps_nal_parse_ctxt);
        if(OK != i4_status)
        {
            return NOT_OK;
        }
        i4_new_lyr_flag = SVCD_FALSE;

        /* Check whether the picture boundary is detected */
        /* or not */
        if(PIC_BOUNDARY_FALSE == i4_pic_bound_status)
        {
            return (OK);
        }

        /* Otherwise look for next nal and compare again */
        *pi4_pic_bound_status = PIC_BOUNDARY_TRUE;
    }

    do
    {
        WORD32 i4_sps_pps_corrupt_status;
        WORD32 i4_tgt_lyr_bckup;
        /* If following conditions are true then there */
        /* is no data left to decode next NAL and hence*/
        /* no further processing is required           */
        if((NAL_END != ps_nal_parse_ctxt->i4_find_nal_state) ||
           ((WORD64) i4_cur_pos >= (WORD64) *pu4_num_bytes))
        {
            return (OK);
        }

        /* Otherwise fill the parameters */
        pu1_buf = pu1_stream_buffer;
        pu1_buf += i4_cur_pos;
        i4_buf_size = *pu4_num_bytes - i4_cur_pos;

        /* Get the NAL prms. This involves the following things*/
        /* 1. Decode the NAL header                            */
        /* 2. Set the discard flag                             */
        /* 3. Decode the slice header if needed                */
        isvcd_set_default_nal_prms(&s_nal_prms);

        /* take a back up of tgt lyr update flag */
        i4_tgt_lyr_bckup = ps_nal_parse_ctxt->i4_tgt_lyr_update;

        /* the tgt attributes should not be  updaetd while pic boundary det*/
        ps_nal_parse_ctxt->i4_tgt_lyr_update = SVCD_FALSE;

        i4_status = isvcd_get_nal_prms(pu1_buf, i4_buf_size, &s_nal_prms, &s_prefix_nal_prms,
                                       &s_prefix_nal_buf, &u4_err_code, &i4_sps_pps_corrupt_status,
                                       &i4_nal_discard_flag, ps_nal_parse_ctxt);
        /* restore back the tgt lyr update flag */
        ps_nal_parse_ctxt->i4_tgt_lyr_update = i4_tgt_lyr_bckup;
        /* If the error code by the nal prms decoder then declare*/
        /* picture boundary                                     */
        if(0 != u4_err_code)
        {
            return (OK);
        }

        i4_more_data_flag = SVCD_FALSE;

        /* If prefix NAL unit comes then save the nal prms*/
        if(PREFIX_UNIT_NAL == s_nal_prms.i4_nal_unit_type)
        {
            UWORD32 u4_bytes_consumed;
            WORD32 i4_status;

            /* If prefix NAL is not discarded then set the varaibles */
            /* appropriatly */
            if(SVCD_FALSE == i4_nal_discard_flag)
            {
                s_prefix_nal_buf.i4_valid_flag = SVCD_TRUE;
                memcpy(&s_prefix_nal_prms, &s_nal_prms, sizeof(nal_prms_t));
            }

            /* Go to next start code */
            i4_zero_cnt = 0;
            u4_bytes_consumed = 0;
            i4_status = isvcd_nal_find_start_code(pu1_stream_buffer, i4_cur_pos, *pu4_num_bytes,
                                                  &i4_zero_cnt, &u4_bytes_consumed);
            /* If associated NAL unit is  not present then */
            if(SC_FOUND != i4_status)
            {
                return (OK);
            }
            i4_cur_pos += u4_bytes_consumed;
            i4_more_data_flag = SVCD_TRUE;
        }
    } while(SVCD_TRUE == i4_more_data_flag);

    /* Do further picture boundary detection only for */
    /* VCL NAL unit (excliding prefix NAL unit)       */
    if((NON_VCL_NAL == s_nal_prms.i4_derived_nal_type) ||
       (PREFIX_UNIT_NAL == s_nal_prms.i4_nal_unit_type) || (SVCD_TRUE == i4_nal_discard_flag))
    {
        return (OK);
    }

    if(SVCD_FALSE == i4_new_lyr_flag)
    {
        if(PIC_BOUND_DQID == i4_pic_bound_type)
        {
            /* If picture boundary was detetcted based on change*/
            /* in DQID then declare picture boundary if DQID of the third slice is different */
            if(i4_prev_dqid != s_nal_prms.i4_dqid)
            {
                return (OK);
            }
        }
        else
        {
            /* If picture boundary was detetcted based on change in DQID */
            /* then declare picture boundary if dependency id of third slice is different */
            if(PIC_BOUND_SLICE_PRMS != i4_pic_bound_type)
            {
                return NOT_OK;
            }

            if((i4_prev_dqid & 0xF) != (s_nal_prms.i4_dqid & 0xF))
            {
                return (OK);
            }
        }

        isvcd_compare_nal_prms(&s_nal_prms, SECOND_PASS, i4_prev_dqid, &i4_pic_bound_type,
                               &i4_pic_bound_status, ps_nal_parse_ctxt);
        *pi4_pic_bound_status = i4_pic_bound_status;

        if(PIC_BOUNDARY_FALSE == i4_pic_bound_status)
        {
            ps_nal_parse_ctxt->i4_prev_dq_id = i4_prev_dqid;
        }
    }
    else
    {
        if(SVCD_TRUE != i4_new_lyr_flag)
        {
            return NOT_OK;
        }
        /* The NAL header is not corrupted only if any of the following conditions are true */
        /* 1. The DQID of the first slice differs with DQID of the third slice */
        /* 2. Picture boundary is detected between first slice and third slice */
        if(i4_prev_dqid == s_nal_prms.i4_dqid)
        {
            isvcd_compare_nal_prms(&s_nal_prms, SECOND_PASS, i4_prev_dqid, &i4_pic_bound_type,
                                   &i4_pic_bound_status, ps_nal_parse_ctxt);
            /* NAL header is corrupted and hence correct it  */
            if(PIC_BOUNDARY_FALSE == i4_pic_bound_status)
            {
                ps_nal_prms->i4_dqid = s_nal_prms.i4_dqid;
                ps_nal_prms->i4_dependency_id = s_nal_prms.i4_dependency_id;
                ps_nal_prms->i4_quality_id = s_nal_prms.i4_quality_id;
                ps_nal_parse_ctxt->i4_prev_dq_id = ps_nal_prms->i4_dqid;
            }
        }
        *pi4_pic_bound_status = PIC_BOUNDARY_FALSE;
    }
    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_insert_vcl_node                                    */
/*                                                                           */
/*  Description   : This routine inserts a DQID layer into DQID list         */
/*                  (this will add a VCL NAL node into VCL NAL structure     */
/*                                                                           */
/*  Inputs        : 1. vcl nal structure                                     */
/*                  2. VCL node to be inserted                               */
/*  Globals       : None                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : Updated vcl nal structure                                */
/*  Returns       : status                                                   */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay      Draft                                     */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_insert_vcl_node(vcl_nal_t *ps_vcl_nal, vcl_node_t *ps_vcl_node)
{
    vcl_node_t *ps_bot_node;
    vcl_node_t *ps_top_node;
    vcl_node_t *ps_node;
    WORD32 i4_rqrd_dqid;

    /* sanity checks */
    if((NULL == ps_vcl_nal) || (NULL == ps_vcl_node))
    {
        return NOT_OK;
    }

    i4_rqrd_dqid = (ps_vcl_node->i4_dependency_id << 4);
    i4_rqrd_dqid += ps_vcl_node->i4_quality_id;
    ps_node = ps_vcl_nal->ps_bot_node;

    /* Search for node which has a DQID which is */
    /* lesser than taht of the node to inserted  */
    while(NULL != ps_node)
    {
        WORD32 i4_dqid;

        i4_dqid = (ps_node->i4_dependency_id << 4);
        i4_dqid += ps_node->i4_quality_id;

        /* If we get a DQID which is greater than*/
        /* the DQID of the  node to be inserted  */
        /* then break out of the loop and update */
        if(i4_dqid > i4_rqrd_dqid)
        {
            ps_bot_node = ps_node->ps_bot_node;
            break;
        }

        ps_node = ps_node->ps_top_node;
    }

    /* If none of the nodes in the list have DQId */
    /* greater than the node to be inserted then  */
    /* bottom node will be top most node          */
    if(NULL == ps_node)
    {
        ps_bot_node = ps_vcl_nal->ps_top_node;
    }

    /* Insert the node into DQID list */
    if(NULL != ps_bot_node)
    {
        ps_top_node = ps_bot_node->ps_top_node;
    }
    else
    {
        ps_top_node = ps_vcl_nal->ps_bot_node;
    }

    /* Join previous node and specified node */
    if(NULL != ps_bot_node)
    {
        ps_bot_node->ps_top_node = ps_vcl_node;
    }
    else
    {
        ps_vcl_nal->ps_bot_node = ps_vcl_node;
    }
    ps_vcl_node->ps_bot_node = ps_bot_node;

    /* Join next node and specified node */
    if(NULL != ps_top_node)
    {
        ps_top_node->ps_bot_node = ps_vcl_node;
    }
    else
    {
        ps_vcl_nal->ps_top_node = ps_vcl_node;
    }
    ps_vcl_node->ps_top_node = ps_top_node;

    return (OK);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_update_nal_ctxt                                     */
/*                                                                           */
/*  Description   : Updates the vcl nal or non vcl structures.               */
/*                                                                           */
/*  Inputs        : ps_nal_parse_ctxt - Bitstream extract context structure  */
/*                  vcl nal structure pointer                                */
/*                  NON vcl nal structure                                    */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : If VCL NAL then adds a node to DQID list                 */
/*                  otherwise adds information to non vcl structure          */
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
void isvcd_update_nal_ctxt(nal_parse_ctxt_t *ps_nal_parse_ctxt, vcl_nal_t *ps_vcl_nal,
                           non_vcl_nal_t *ps_non_vcl_nal)
{
    /*! If current NAL is VCL NAL then
          - Insert a VCL node into DQID list if neccessery
          - update the information part of NAL unit */
    /*! Otherwise, populate the buffer parameters into non vcl output
    structure */
    nal_prms_t *ps_nal_prms;
    nal_buf_t *ps_nal_buf, *ps_prefix_nal_buf;

    ps_nal_prms = &ps_nal_parse_ctxt->s_nal_prms;
    ps_nal_prms = &ps_nal_parse_ctxt->s_nal_prms;
    ps_nal_buf = &ps_nal_parse_ctxt->s_nal_buf;
    ps_prefix_nal_buf = &ps_nal_parse_ctxt->s_prefix_nal_buf;

    /* If prefix NAL unit then          */
    /* - calculate the SODB length      */
    if(PREFIX_UNIT_NAL == ps_nal_prms->i4_nal_unit_type)
    {
        /* Since we consume the zeroes in start code also */
        /* size has to reduced                            */
        if(NAL_END == ps_nal_parse_ctxt->i4_find_nal_state)
        {
            ps_prefix_nal_buf->i4_buf_size -= 2;
        }

        ps_prefix_nal_buf->u4_max_bits =
            isvcd_nal_rbsp_to_sodb(ps_prefix_nal_buf->pu1_buf, ps_prefix_nal_buf->i4_buf_size, 0);
        memcpy(&ps_nal_parse_ctxt->s_prefix_nal_prms, &ps_nal_parse_ctxt->s_nal_prms,
               sizeof(nal_prms_t));
        return;
    }

    if(ANNEX_B == ps_nal_parse_ctxt->i4_input_bitstream_mode)
    {
        /* Since we consume the zeroes in start code also */
        /* size has to reduced                            */
        if(NAL_END == ps_nal_parse_ctxt->i4_find_nal_state)
        {
            ps_nal_buf->i4_buf_size -= 2;
        }
    }

    if(VCL_NAL == ps_nal_prms->i4_derived_nal_type)
    {
        dqid_node_t *ps_dqid_node;
        vcl_node_t *ps_node;
        WORD32 i4_status;
        dec_pic_params_t *ps_pps;
        dec_seq_params_t *ps_sps;
        vcl_buf_hdr_t *ps_vcl_hdr;
        vcl_buf_hdr_t *ps_prev_vcl_hdr;
        WORD32 i4_slice_offset;

        ps_sps = ps_nal_parse_ctxt->pv_seq_prms;
        ps_sps += ps_nal_prms->u1_sps_id;
        ps_pps = ps_nal_parse_ctxt->pv_pic_prms;
        ps_pps += ps_nal_prms->u1_pps_id;

        /* Get the VCL NAL node */
        i4_status = isvcd_get_dqid_node(&ps_nal_parse_ctxt->s_dqid_ctxt,
                                        (UWORD8) ps_nal_parse_ctxt->i4_prev_dq_id, &ps_dqid_node);

        ps_node = ps_dqid_node->ps_vcl_node;

        if(NULL == ps_node)
        {
            /* no active node has been acquired */
            return;
        }

        /*-------------------------------------------------------------------*/
        /* The DQID list updation should happen only once in a               */
        /* layer. Hence a flag used to determine whether the                 */
        /* layer is already initialized or not.                              */
        /*-------------------------------------------------------------------*/
        if(SVCD_FALSE == ps_dqid_node->u1_valid_flag)
        {
            /* Update the DQID node */
            ps_dqid_node->u1_valid_flag = SVCD_TRUE;
            ps_dqid_node->u1_dqid = (ps_nal_prms->i4_dependency_id << 4);
            ps_dqid_node->u1_dqid += ps_nal_prms->i4_quality_id;
            ps_dqid_node->i4_poc_lsb = ps_nal_prms->i4_poc_lsb;
            ps_dqid_node->i4_delta_poc_bot = ps_nal_prms->i4_delta_poc_bot;
            ps_dqid_node->ai4_delta_poc[0] = ps_nal_prms->ai4_delta_poc[0];
            ps_dqid_node->ai4_delta_poc[1] = ps_nal_prms->ai4_delta_poc[1];

            /* Update the VCL node */
            ps_node->i4_quality_id = ps_nal_prms->i4_quality_id;
            ps_node->i4_dependency_id = ps_nal_prms->i4_dependency_id;
            ps_node->i4_temporal_id = ps_nal_prms->i4_temporal_id;
            ps_node->i4_priority_id = ps_nal_prms->i4_priority_id;
            ps_node->i4_idr_pic_flag = ps_nal_prms->i4_idr_pic_flag;
            ps_node->i4_nal_ref_idc = ps_nal_prms->i4_nal_ref_idc;
            ps_node->i4_nal_unit_type = ps_nal_prms->i4_nal_unit_type;
            ps_node->i4_use_ref_base = ps_nal_prms->i4_use_ref_base_pic_flag;
            ps_node->i4_nal_ref_idc = ps_nal_prms->i4_nal_ref_idc;
            ps_node->u1_sps_id = ps_nal_prms->u1_sps_id;
            ps_node->u1_pps_id = ps_nal_prms->u1_pps_id;
            ps_node->u2_frm_num = ps_nal_prms->u2_frm_num;
            ps_node->i4_idr_pic_num = ps_nal_prms->i4_idr_pic_num;
            ps_node->i4_num_slices = 0;
            ps_node->u1_acc_no_int_pred = 1;
            if(0 == ps_sps->u1_pic_order_cnt_type)
            {
                ps_node->i4_poc_syntax = ps_nal_prms->i4_poc_lsb;
            }
            else
            {
                ps_node->i4_poc_syntax = ps_nal_prms->ai4_delta_poc[0];
            }

            /* Insert the node into DQID list */
            i4_status = isvcd_insert_vcl_node(ps_vcl_nal, ps_node);
            if(OK != i4_status)
            {
                return;
            }

            /* Reset the previous field */
            ps_nal_parse_ctxt->ps_prev_vcl_buf = NULL;
            ps_node->ps_first_vcl_nal = NULL;
        }

        /* Update accumulated no inter layer prediction */
        ps_node->u1_acc_no_int_pred &= (UWORD8) ps_nal_prms->i4_no_int_lyr_pred;

        /****************** Fill VCL BUF header ************/

        /* If prefix NAL unit is present then update  */
        /* the following                              */
        /* - Start of buffer header will be present in*/
        /*   before the start of prefix NAL unit's SODB*/
        /*   data.                                    */
        /*   Note: If memeory left for buffer header  */
        /*   of the prefix NAL unit will have junk    */
        /*   values                                   */

        if(NULL == ps_nal_buf->pu1_buf)
        {
            /* no nal needs to be added into the list hence return */
            return;
        }
        else
        {
            ps_vcl_hdr = (vcl_buf_hdr_t *) (ps_nal_buf->pu1_buf - GET_NAL_BUF_INC(VCL_NAL));
        }

        i4_slice_offset = 0;
        if(SVCD_TRUE == ps_prefix_nal_buf->i4_valid_flag)
        {
            ps_vcl_hdr = (vcl_buf_hdr_t *) (ps_prefix_nal_buf->pu1_buf - GET_NAL_BUF_INC(VCL_NAL));
            i4_slice_offset = ps_nal_buf->pu1_buf - ps_prefix_nal_buf->pu1_buf;
        }

        /* Update the next field of the previous nal  */
        /* unit or if it is the first NAL then update */
        /* VCL node information                       */
        ps_prev_vcl_hdr = ps_nal_parse_ctxt->ps_prev_vcl_buf;
        if(NULL != ps_prev_vcl_hdr)
        {
            ps_prev_vcl_hdr->ps_next = ps_vcl_hdr;
        }
        else
        {
            ps_node->ps_first_vcl_nal = ps_vcl_hdr;
        }

        /* Fill the VCL buffer header */
        ps_vcl_hdr->ps_next = NULL;
        ps_vcl_hdr->i4_no_int_lyr_pred = ps_nal_prms->i4_no_int_lyr_pred;
        ps_vcl_hdr->i4_first_mb_addr = ps_nal_prms->u4_first_mb_addr;
        ps_vcl_hdr->u4_prefix_nal_bits = ps_prefix_nal_buf->u4_max_bits;
        ps_vcl_hdr->i4_slice_offset = 0;
        ps_vcl_hdr->i4_buf_offset = GET_NAL_BUF_INC(VCL_NAL);
        ps_vcl_hdr->i4_slice_offset = i4_slice_offset;

        /* Determine max num bits */
        ps_nal_buf->u4_max_bits = isvcd_nal_rbsp_to_sodb(
            ps_nal_buf->pu1_buf, ps_nal_buf->i4_buf_size, ps_pps->u1_entropy_coding_mode);
        ps_vcl_hdr->u4_max_bits = ps_nal_buf->u4_max_bits;

        /* Updates */
        ps_nal_parse_ctxt->ps_prev_vcl_buf = ps_vcl_hdr;
        ps_node->i4_num_slices += 1;
    }
    /*-----------------------------------------------------------------------*/
    /* If start of NAL and if its a NON VCL NAL then update the              */
    /* start address of the NON VCL NAL                                      */
    /*-----------------------------------------------------------------------*/
    else
    {
        non_vcl_buf_hdr_t *ps_non_vcl_buf_hdr;
        non_vcl_buf_hdr_t *ps_prev_non_vcl_buf_hdr;

        ps_non_vcl_buf_hdr =
            (non_vcl_buf_hdr_t *) (ps_nal_buf->pu1_buf - GET_NAL_BUF_INC(NON_VCL_NAL));

        /* Update NON VCL structure */
        ps_non_vcl_buf_hdr->i4_nal_unit_type = ps_nal_prms->i4_nal_unit_type;
        ps_non_vcl_buf_hdr->ps_next = NULL;
        ps_non_vcl_buf_hdr->i4_buf_offset = GET_NAL_BUF_INC(NON_VCL_NAL);
        ps_non_vcl_buf_hdr->i4_buf_size = ps_nal_buf->i4_buf_size;

        /* Update the next field and first non vcl fields of */
        /* non vcl buffer header structure and non vcl       */
        /* structure respectively                            */
        ps_prev_non_vcl_buf_hdr = ps_nal_parse_ctxt->ps_prev_non_vcl_buf;
        if(NULL != ps_prev_non_vcl_buf_hdr)
        {
            ps_prev_non_vcl_buf_hdr->ps_next = ps_non_vcl_buf_hdr;
        }
        else
        {
            ps_non_vcl_nal->ps_first_non_vcl_nal = ps_non_vcl_buf_hdr;
        }

        /* Updates */
        ps_nal_parse_ctxt->i4_num_non_vcl_nals += 1;
        ps_non_vcl_nal->i4_num_non_vcl_nals = ps_nal_parse_ctxt->i4_num_non_vcl_nals;
        ps_nal_parse_ctxt->ps_prev_non_vcl_buf = ps_non_vcl_buf_hdr;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_idr_err_hdlr                                        */
/*                                                                           */
/*  Description   : This routine shall be invoked to handle a case when a    */
/*                  slice is an IDR picture and it is referring to corrupted */
/*                  SPS or PPS                                               */
/*                                                                           */
/*  Inputs        : 1. VCL NAL structure                                     */
/*                  2. NAL paramters                                         */
/*                  3. NAL parse context structure                           */
/*  Globals       : None                                                     */
/*  Processing    : It will set the highest available dependency id below the*/
/*                  current dependency id as the target layer. Also sets the */
/*                  update target layer flag to FALSE as target layer need not*/
/*                  adopt to the application's target layer in the current   */
/*                  picture                                                  */
/*                                                                           */
/*  Outputs       : Updated vcl nal structure                                */
/*                  Updated internal target layer attributes                 */
/*                  Updated target layer update flag                         */
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
WORD32 isvcd_idr_err_hdlr(vcl_nal_t *ps_vcl_nal, nal_prms_t *ps_nal_prms,
                          nal_parse_ctxt_t *ps_nal_parse_ctxt)
{
    vcl_node_t *ps_vcl_node;
    target_lyr_attr_t *ps_int_attr;

    /* sanity checks */
    if((NULL == ps_vcl_nal) || (NULL == ps_nal_prms) || (NULL == ps_nal_parse_ctxt))
    {
        return NOT_OK;
    }
    UNUSED(ps_nal_prms);

    /* Initializations */
    ps_vcl_node = ps_vcl_nal->ps_top_node;
    ps_int_attr = &ps_nal_parse_ctxt->s_int_attr;

    /* the highest node present in the depedency list will be         */
    /* considered as targte layer and appropriate params will be used */

    /* If not found then delete all the layers in the AU */
    if(NULL == ps_vcl_node)
    {
        ps_int_attr->i4_dependency_id = -1;
        ps_int_attr->i4_quality_id = MAX_QUALITY_ID;
    }
    else
    {
        /* Set the target layer */
        ps_int_attr->i4_dependency_id = ps_vcl_node->i4_dependency_id;
        ps_int_attr->i4_quality_id = ps_vcl_node->i4_quality_id;
    }

    return (OK);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name :svcd_refine_dqid_list                                     */
/*                                                                           */
/*  Description   : Inserts the dummy nodes for each dependency id which     */
/*                  have not come in the bitstream                           */
/*                                                                           */
/*  Inputs        :VCL NAL structure                                         */
/*                  NAL parse context structure                              */
/*  Globals       : None                                                     */
/*  Processing    : For each dependency id till the target dependency id     */
/*                  - If layer already exists (came in the bitstream) then   */
/*                    do nothing                                             */
/*                  - Otherwsie insert a dummy node                          */
/*                                                                           */
/*  Outputs       : Updated VCL NAL structure                                */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay           Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_refine_dqid_list(vcl_nal_t *ps_vcl_nal, nal_parse_ctxt_t *ps_nal_parse_ctxt)
{
    vcl_node_t *ps_node;
    target_lyr_attr_t *ps_int_attr;
    dqid_ctxt_t *ps_dqid_ctxt;
    UWORD8 u1_dep_id;
    WORD32 i4_status;
    WORD32 i4_dep_id;

    ps_int_attr = &ps_nal_parse_ctxt->s_int_attr;
    ps_dqid_ctxt = &ps_nal_parse_ctxt->s_dqid_ctxt;
    i4_dep_id = -1;

    for(u1_dep_id = 0; u1_dep_id <= ps_int_attr->i4_dependency_id; u1_dep_id++)
    {
        dqid_node_t *ps_dqid_node;

        /* Get a DQID node */
        i4_status = isvcd_get_dqid_node(ps_dqid_ctxt, (UWORD8) (u1_dep_id << 4), &ps_dqid_node);
        if(OK != i4_status)
        {
            return NOT_OK;
        }

        /* If node does not exist already then insert a dummy node */
        if(SVCD_FALSE == ps_dqid_node->u1_valid_flag)
        {
            if(1 == ps_nal_parse_ctxt->i4_idr_pic_err_flag)
            {
                ps_int_attr->i4_dependency_id = i4_dep_id;
                ps_int_attr->i4_quality_id = MAX_QUALITY_ID;

                /* remove all the nodes from dependency list */
                /* which are at higher dependency than the   */
                /* value set in init attributes              */
                while(NULL != ps_vcl_nal->ps_top_node)
                {
                    /* if higher dependency */
                    if(ps_vcl_nal->ps_top_node->i4_dependency_id > i4_dep_id)
                    {
                        ps_vcl_nal->ps_top_node = ps_vcl_nal->ps_top_node->ps_bot_node;
                    }
                    else
                    {
                        break;
                    }
                }

                /* if no node exists in the dependency list */
                if(NULL == ps_vcl_nal->ps_top_node)
                {
                    ps_vcl_nal->ps_bot_node = NULL;
                }
                else if(ps_vcl_nal->ps_top_node == ps_vcl_nal->ps_bot_node)
                {
                    /* if a single node exists */
                    ps_vcl_nal->ps_top_node->ps_bot_node = NULL;
                    ps_vcl_nal->ps_bot_node->ps_top_node = NULL;
                }

                return (NOT_OK);
            }
            else
            {
                ps_dqid_node->u1_valid_flag = SVCD_TRUE;
                ps_dqid_node->u1_dqid = (u1_dep_id << 4);

                /* Fill VCL node information */
                ps_node = ps_dqid_node->ps_vcl_node;
                ps_node->i4_dependency_id = u1_dep_id;
                ps_node->i4_quality_id = 0;
                ps_node->ps_first_vcl_nal = NULL;
            }

            /* Insert node into DQID list */
            i4_status = isvcd_insert_vcl_node(ps_vcl_nal, ps_node);
            if(OK != i4_status)
            {
                return (NOT_OK);
            }
        }

        i4_dep_id++;
    } /* End of loop over all the dependency id */
    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_parse_set_target_attr                          */
/*                                                                           */
/*  Description   : Sets the target layer attributes                         */
/*                                                                           */
/*  Inputs        : i4_target_quality_id - Target layer quality id           */
/*                  i4_target_dependency_id - Target layer dependency id     */
/*                  i4_target_temporal_id - Target layer temporal id         */
/*                  i4_target_priority_id - Target layer priority id         */
/*                  pv_nal_parse_ctxt - Pointer module handle                */
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
/*          06 09 2021   Vijay           Draft                               */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_nal_parse_set_target_attr(WORD32 i4_target_quality_id, WORD32 i4_target_dependency_id,
                                       WORD32 i4_target_temporal_id, WORD32 i4_target_priority_id,
                                       void *pv_nal_parse_ctxt)
{
    nal_parse_ctxt_t *ps_nal_parse_ctxt;
    target_lyr_attr_t *ps_app_attr;

    if((i4_target_quality_id > MAX_QUALITY_ID) || (i4_target_dependency_id > MAX_DEPENDENCY_ID))
    {
        return IV_FAIL;
    }

    ps_nal_parse_ctxt = (nal_parse_ctxt_t *) pv_nal_parse_ctxt;
    ps_app_attr = &ps_nal_parse_ctxt->s_app_attr;

    /*-----------------------------------------------------------------------*/
    /*! Register the target information into context structure               */
    /*-----------------------------------------------------------------------*/
    ps_app_attr->i4_quality_id = i4_target_quality_id;
    ps_app_attr->i4_dependency_id = i4_target_dependency_id;
    ps_app_attr->i4_temporal_id = i4_target_temporal_id;
    ps_app_attr->i4_priority_id = i4_target_priority_id;
    return IV_SUCCESS;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_parse_reset_ctxt                               */
/*                                                                           */
/*  Description   : Initializes the bitstream extraction module. Should be   */
/*                  called once in a sequence                                */
/*                                                                           */
/*  Inputs        : i4_input_bitstream_mode - Input bitstream mode RFC or    */
/*                      Annex B                                              */
/*                  i4_input_mode - Input mode - Full input mode or partial  */
/*                      input mode                                           */
/*                  pv_nal_parse_ctxt - Module handle                        */
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
/*          06 09 2021   Vijay           Draft                               */
/*                                                                           */
/*****************************************************************************/

void isvcd_nal_parse_reset_ctxt(WORD32 i4_input_bitstream_mode, WORD32 i4_input_mode,
                                void *pv_nal_parse_ctxt)
{
    nal_parse_ctxt_t *ps_nal_parse_ctxt = (nal_parse_ctxt_t *) pv_nal_parse_ctxt;
    UNUSED(i4_input_mode);

    /*-----------------------------------------------------------------------*/
    /*! Set the input bitstream mode of context structure                    */
    /*-----------------------------------------------------------------------*/
    switch(i4_input_bitstream_mode)
    {
        case ANNEX_B:
        case NON_ANNEX_B:
            break;
        default:
            break;
    }

    ps_nal_parse_ctxt->i4_input_bitstream_mode = i4_input_bitstream_mode;

    /*-----------------------------------------------------------------------*/
    /*! Perform the picture level initialization                             */
    /*-----------------------------------------------------------------------*/
    isvcd_pic_reset_ctxt(pv_nal_parse_ctxt);

    /* Reset the prefix nal unit buffer structure */
    isvcd_nal_buf_reset(&ps_nal_parse_ctxt->s_prefix_nal_buf);

    /*-----------------------------------------------------------------------*/
    /*! Reset other varaibles                                                */
    /*-----------------------------------------------------------------------*/
    ps_nal_parse_ctxt->i4_dec_frst_sc_flag = SVCD_TRUE;
    ps_nal_parse_ctxt->i4_eos_flag = SVCD_FALSE;
    ps_nal_parse_ctxt->u1_pic_boundary_aud_flag = 0;
    ps_nal_parse_ctxt->u4_bytes_left = 0;

    /* Reset target layer attributes */
    {
        target_lyr_attr_t *ps_app_attr;
        target_lyr_attr_t *ps_int_attr;

        ps_app_attr = &ps_nal_parse_ctxt->s_app_attr;
        ps_int_attr = &ps_nal_parse_ctxt->s_int_attr;

        ps_app_attr->i4_dependency_id = MAX_DEPENDENCY_ID;
        ps_app_attr->i4_quality_id = MAX_QUALITY_ID;
        ps_app_attr->i4_temporal_id = MAX_TEMPORAL_ID;
        ps_app_attr->i4_priority_id = MAX_PRIORITY_ID;

        ps_int_attr->i4_dependency_id = -1;
        ps_int_attr->i4_quality_id = MAX_QUALITY_ID;
        ps_int_attr->i4_temporal_id = 0;
        ps_int_attr->i4_priority_id = MAX_PRIORITY_ID;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_parse_partial_signal_eos                       */
/*                                                                           */
/*  Description   : Does processing when end of stream occurs for partial    */
/*                  input mode of operation.                                 */
/*                                                                           */
/*  Inputs        : pv_nal_parse_ctxt - bitstream extract context structure  */
/*                  pv_out_vcl_nal - vcl nal structure                       */
/*                  pv_out_non_vcl_nal - non vcl nal structure               */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : None                                                     */
/*                                                                           */
/*  Returns       : Picture boundary detetcted or not                        */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay           Draft                               */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_nal_parse_partial_signal_eos(void *pv_nal_parse_ctxt, void *pv_out_vcl_nal,
                                          void *pv_out_non_vcl_nal)
{
    nal_parse_ctxt_t *ps_nal_parse_ctxt;
    vcl_nal_t *ps_vcl_nal;

    ps_nal_parse_ctxt = (nal_parse_ctxt_t *) pv_nal_parse_ctxt;
    ps_vcl_nal = (vcl_nal_t *) pv_out_vcl_nal;

    /* for RFC mode */
    if(NON_ANNEX_B == ps_nal_parse_ctxt->i4_input_bitstream_mode)
    {
        /* Reset the end of stream flag so that in    */
        ps_nal_parse_ctxt->i4_eos_flag = SVCD_TRUE;
    }

    if(1 == ps_nal_parse_ctxt->u1_pic_boundary_aud_flag)
    {
        ps_nal_parse_ctxt->i4_eos_flag = SVCD_TRUE;
    }
    /* Update VCL node if it is first call in the */
    /* flush mode                                 */
    if(SVCD_FALSE == ps_nal_parse_ctxt->i4_eos_flag)
    {
        WORD32 i4_status;

        /* Update the unfinished NAL into VCL node if */
        /* all the following conditions are true      */
        /* 1. We have not found the start code and    */
        /*    NAL boundary is not detected yet        */
        /* 2. NAL is not discarded                    */
        if((FIND_NAL_END == ps_nal_parse_ctxt->i4_find_nal_state) &&
           (SVCD_FALSE == ps_nal_parse_ctxt->i4_discard_nal_flag))
        {
            isvcd_update_nal_ctxt(ps_nal_parse_ctxt, pv_out_vcl_nal, pv_out_non_vcl_nal);
        }

        ps_nal_parse_ctxt->i4_idr_pic_err_flag = 0;
        /* Refine based on the no inter layer pred flag*/
        i4_status = isvcd_refine_dqid_list(ps_vcl_nal, ps_nal_parse_ctxt);

        if(!(OK == i4_status))
        {
            return i4_status;
        }
        UNUSED(i4_status);

        /* Reset the context structure variables */
        isvcd_nal_reset_ctxt(ps_nal_parse_ctxt);

        /* Reset the end of stream flag so that in    */
        /* the next flush call the above steps need   */
        /* not be performed                           */
        ps_nal_parse_ctxt->i4_eos_flag = SVCD_TRUE;

        return (PIC_BOUNDARY_TRUE);
    }
    else
    {
        return (FLUSH_DECODED_PICTURE);
    }
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_parse_pic_bound_proc                           */
/*                                                                           */
/*  Description   : Function does the picture end processign and resets      */
/*                                                                           */
/*                                                                           */
/*  Inputs        : ps_nal_parse_ctxt, ps_vcl_nal                            */
/*  Globals       : none                                                     */
/*  Processing    : DQid list refiniment and resets                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   vijayakumar          creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_nal_parse_pic_bound_proc(nal_parse_ctxt_t *ps_nal_parse_ctxt, vcl_nal_t *ps_vcl_nal,
                                    nal_prms_t *ps_nal_prms)
{
    WORD32 i4_status;

    i4_status = isvcd_refine_dqid_list(ps_vcl_nal, ps_nal_parse_ctxt);

    /* in case of IDR pictures if the node     */
    /* which has to be added into dependency   */
    /* list is not valied then the layer below */
    /* that node is set as target layer        */

    if(NOT_OK == i4_status)
    {
        ps_nal_parse_ctxt->i4_discard_nal_flag = SVCD_TRUE;
        ps_vcl_nal->i1_nal_ref_id_next = -1;
    }
    else
    {
        /* update the next access unit params */
        /* will be used by lower level decoder*/
        /* for concealment of frame number    */
        /* applicable for single layer cases  */
        ps_vcl_nal->i1_nal_ref_id_next = ps_nal_prms->i4_nal_ref_idc;

        ps_vcl_nal->u2_frm_num_next = ps_nal_prms->u2_frm_num;
    }

    /* -------- reset few variables in context structure ----*/
    isvcd_pic_reset_ctxt(ps_nal_parse_ctxt);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_parse_vcl_nal_partial                          */
/*                                                                           */
/*  Description   : None                                                     */
/*                                                                           */
/*  Inputs        : pv_nal_parse_ctxt - bitstream extract ctxt               */
/*                      structure                                            */
/*                  pv_input_bitstream_ctxt - bitstream context              */
/*                  pv_out_non_vcl_nal - non vcl nal structure (output)      */
/*                  pv_out_vcl_nal - vcl nal structure (output)              */
/*                  pu4_bytes_consumed - bytes consumed variable(output)     */
/*                  pi4_num_packets_consumed - packets consumed (output/RFC) */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : Updates bytes consumed variable, packets consumed,       */
/*                  output structures (vcl nal , non vcl nal)                */
/*                                                                           */
/*  Returns       : If picture bounadry is detetcted then PIC_BOUNDARY_TRUE  */
/*                  otherwise PIC_BOUNDARY_FALSE                             */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay           Draft                               */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_nal_parse_vcl_nal_partial(void *pv_nal_parse_ctxt, UWORD8 *pu1_stream_buffer,
                                       void *pv_out_non_vcl_nal, void *pv_out_vcl_nal,
                                       UWORD32 *pu4_bytes_consumed, UWORD32 *pu4_num_bytes)
{
    /*! - Search for the NAL boundary
        - If NAL boundary is not found and bytes consumed is lesser than
          minimum buffer size then break out of the loop
        - if it is start of NAL then read the NAL header
        - If it is a VCL NAL then invoke picture boundary detection logic and
          picture boundary is detected then break out of the loop without
          updating the bytes consumed variable
        - NAL discard logic determines whther the current NAL has to be
          discarded or not
        - If NAL is not discarded then populate the vcl or non vcl output
          structures
    */
    nal_parse_ctxt_t *ps_nal_parse_ctxt;
    vcl_nal_t *ps_vcl_nal;
    non_vcl_nal_t *ps_non_vcl_nal;
    nal_unit_t *ps_nal_unit;
    WORD32 i4_nal_start_flag, i4_cur_pos, i4_status;
    WORD32 i4_nal_header_len, i4_more_data_flag;
    UWORD32 u4_bytes_consumed_temp = 0;
    UWORD8 **ppu1_out_buf;
    nal_prms_t *ps_nal_prms;
    WORD32 i4_pic_bound_status;

    ps_nal_parse_ctxt = (nal_parse_ctxt_t *) pv_nal_parse_ctxt;
    ps_vcl_nal = (vcl_nal_t *) pv_out_vcl_nal;
    ps_non_vcl_nal = (non_vcl_nal_t *) pv_out_non_vcl_nal;
    ps_nal_unit = (nal_unit_t *) ps_nal_parse_ctxt->pv_nal_unit;
    ps_nal_prms = &ps_nal_parse_ctxt->s_nal_prms;

    /* Initialization */
    i4_cur_pos = 0;
    *pu4_bytes_consumed = 0;
    i4_nal_header_len = 0;
    i4_nal_start_flag = SVCD_FALSE;
    i4_more_data_flag = SVCD_TRUE;
    i4_pic_bound_status = PIC_BOUNDARY_FALSE;

    ps_non_vcl_nal->i4_num_non_vcl_nals = ps_nal_parse_ctxt->i4_num_non_vcl_nals;

    /* Since we do not perform the picture boundary detection */
    /* on the prefix NAL unit, the current picture's prefix   */
    /* NAL unit will be at the bottom of the buffer. Hence    */
    /* it should be copied to top of the buffer               */
    if(SVCD_TRUE == ps_nal_parse_ctxt->i4_is_frst_vcl_nal_in_au)
    {
        nal_buf_t *ps_prefix_nal_buf;

        ps_prefix_nal_buf = &ps_nal_parse_ctxt->s_prefix_nal_buf;
        if(SVCD_TRUE == ps_prefix_nal_buf->i4_valid_flag)
        {
            WORD32 i4_buf_size;
            UWORD8 *pu1_vcl_nal;

            if(ps_prefix_nal_buf->i4_buf_size > 0)
            {
                i4_buf_size = ps_prefix_nal_buf->i4_buf_size;
                i4_buf_size = UP_ALIGN_8(i4_buf_size + BUFFER_ALIGN_4);
            }
            else
            {
                i4_buf_size = 0;
            }

            pu1_vcl_nal = ps_nal_parse_ctxt->pu1_vcl_nal_buf + i4_buf_size;

            memmove(ps_nal_parse_ctxt->pu1_vcl_nal_buf, ps_prefix_nal_buf->pu1_buf, i4_buf_size);
            ps_prefix_nal_buf->pu1_buf = ps_nal_parse_ctxt->pu1_vcl_nal_buf;
            ps_nal_parse_ctxt->pu1_vcl_nal_buf = pu1_vcl_nal;

            /* subtract the buffer size left */
            ps_nal_parse_ctxt->u4_bytes_left_vcl -= i4_buf_size;
        }
        /* Reset the top and bottom node */
        ps_vcl_nal->ps_top_node = NULL;
        ps_vcl_nal->ps_bot_node = NULL;
        ps_vcl_nal->i1_nal_ref_id_next = -1;
        ps_vcl_nal->u2_frm_num_next = 0;
    }

    /* If number of bytes left in the previous process call  */
    /* is is greater or equal to number of bytes in input    */
    /* buffer of the current process call then declare that  */
    /* end of bitstream has occurred and consume the bytes   */
    /* but do not decode                                     */
    if(ps_nal_parse_ctxt->u4_bytes_left >= (UWORD32) *pu4_num_bytes)
    {
        ps_nal_parse_ctxt->i4_discard_nal_flag = SVCD_TRUE;
        *pu4_bytes_consumed = *pu4_num_bytes;

        i4_status =
            isvcd_nal_parse_partial_signal_eos(ps_nal_parse_ctxt, ps_vcl_nal, ps_non_vcl_nal);
        /* set the next AU params to default values */
        ps_vcl_nal->i1_nal_ref_id_next = -1;
        ps_vcl_nal->u2_frm_num_next = 0;

        return (i4_status);
    }
    ps_nal_parse_ctxt->u4_bytes_left = 0;

    /*************************************************************************/
    /*                      LOOP OVER NALs                                   */
    /*************************************************************************/
    do
    {
        nal_buf_t *ps_nal_buf;
        UWORD32 *pu4_bytes_left;

        /* Find NAL boundary                */
        if(ANNEX_B == ps_nal_parse_ctxt->i4_input_bitstream_mode)
        {
            i4_nal_start_flag = isvcd_get_annex_b_nal_unit(
                pu1_stream_buffer, i4_cur_pos, *pu4_num_bytes,
                &ps_nal_parse_ctxt->i4_find_nal_state, &ps_nal_parse_ctxt->i4_zero_byte_cnt,
                &u4_bytes_consumed_temp, ps_nal_parse_ctxt->pv_nal_unit, &i4_more_data_flag);

            i4_cur_pos += u4_bytes_consumed_temp;
        }

        /*********************************************************************/
        /*          READ NAL HEADER AND NAL DISCARD LOGIC                    */
        /*********************************************************************/

        /* If it is the start of NAL header perform the following */
        /* 1. Decode NAL header                                   */
        /* 2. Determine whether the NAL has to be discarded or not*/
        /* 3. Detect the picture boundary                         */
        if(SVCD_TRUE == i4_nal_start_flag)
        {
            UWORD32 u4_err_code;
            WORD32 i4_sps_pps_corrupt_status;
            WORD32 i4_internal_dep_id_prev;

            /* Get the NAL prms. This involves the following things*/
            /* 1. Decode the NAL header                            */
            /* 2. Set the discard flag                             */
            /* 3. Decode the slice header if needed                */

            /* get the dependency id at which the NAl parse is currently */
            /* present */
            i4_internal_dep_id_prev = ps_nal_parse_ctxt->s_int_attr.i4_dependency_id;

            i4_status = isvcd_get_nal_prms(
                ps_nal_unit->pu1_bufs, ps_nal_unit->i4_buf_sizes, ps_nal_prms,
                &ps_nal_parse_ctxt->s_prefix_nal_prms, &ps_nal_parse_ctxt->s_prefix_nal_buf,
                &u4_err_code, &i4_sps_pps_corrupt_status, &ps_nal_parse_ctxt->i4_discard_nal_flag,
                ps_nal_parse_ctxt);

            if(NON_ANNEX_B == ps_nal_parse_ctxt->i4_input_bitstream_mode)
            {
                ps_nal_parse_ctxt->i4_prev_dq_id = ps_nal_prms->i4_dqid;
            }

            /* If the error code returned by the "picture boundary" */
            /* detetction is                                        */
            /* 1. Insufficient bitstream size: then store the bytes */
            /*    left and break out of the loop                    */
            /* 2. Corrupted slice: then discard the slice           */
            if((NAL_INSUFFICIENT_DATA == (WORD32) u4_err_code) &&
               (NAL_END != ps_nal_parse_ctxt->i4_find_nal_state))
            {
                ps_nal_parse_ctxt->u4_bytes_left = *pu4_num_bytes - *pu4_bytes_consumed;

                /* Reset the NAL level tracking variables */
                isvcd_nal_reset_ctxt(ps_nal_parse_ctxt);
                break;
            }
            else if(0 != u4_err_code)
            {
                ps_nal_parse_ctxt->i4_discard_nal_flag = SVCD_TRUE;

                if(SVCD_TRUE == ps_nal_prms->i4_idr_pic_flag)
                {
                    /* IDR Error handler is called       */
                    /* only if for a given layer the NAL */
                    /* haeder and partial slice decode   */
                    /* routine comes out as no SPS PPS   */
                    /* error. But for Lowest layer in    */
                    /* access unit it is doen always     */
                    if(i4_internal_dep_id_prev != ps_nal_parse_ctxt->s_int_attr.i4_dependency_id)
                    {
                        /* if the target depedency id has been */
                        /* changed while decoding currnet NAL  */

                        if((0 != i4_sps_pps_corrupt_status) ||
                           (-1 == ps_nal_parse_ctxt->i4_prev_dq_id))
                        {
                            i4_status =
                                isvcd_idr_err_hdlr(ps_vcl_nal, ps_nal_prms, ps_nal_parse_ctxt);
                            if(OK != i4_status)
                            {
                                return i4_status;
                            }
                            UNUSED(i4_status);

                            ps_nal_parse_ctxt->i4_tgt_lyr_update = SVCD_FALSE;
                        }
                        else
                        {
                            if(0 == ps_nal_prms->i4_quality_id)
                            {
                                /* over write the frame number */
                                ps_nal_parse_ctxt->s_nal_prms.u2_frm_num = 0;

                                /* Get the previous layer's DQID */
                                if(ps_nal_parse_ctxt->i4_prev_dq_id < ps_nal_prms->i4_dqid)
                                {
                                    ps_nal_parse_ctxt->i4_prev_dq_id = ps_nal_prms->i4_dqid;
                                    ps_nal_parse_ctxt->i4_is_frst_vcl_nal_in_au = SVCD_FALSE;
                                }

                                /* update the nal context with the nal */
                                /* header params */
                                isvcd_update_nal_ctxt(ps_nal_parse_ctxt, ps_vcl_nal,
                                                      ps_non_vcl_nal);
                            }
                        }
                    }
                }
            }

            /* Populate the derived nal type into bitstream extract*/
            /* context structure                                   */
            i4_nal_header_len = ps_nal_prms->i4_nal_header_len;
            ps_nal_parse_ctxt->i4_nal_type = ps_nal_prms->i4_derived_nal_type;

            /* get the accumulated idr pic error flag */
            ps_nal_parse_ctxt->i4_idr_pic_err_flag |=
                ((SVCD_TRUE == ps_nal_prms->i4_idr_pic_flag) &&
                 (SVCD_FALSE == ps_nal_parse_ctxt->i4_discard_nal_flag) &&
                 (i4_internal_dep_id_prev != ps_nal_parse_ctxt->s_int_attr.i4_dependency_id));

            if(ACCESS_UNIT_DELIMITER_RBSP == ps_nal_prms->i4_nal_unit_type)
            {
                i4_pic_bound_status = PIC_BOUNDARY_TRUE;
                ps_nal_parse_ctxt->u1_pic_boundary_aud_flag = 1;
                /* If picture boundary is detected then come out of  */
                /* the loop                                          */
                if(PIC_BOUNDARY_TRUE == i4_pic_bound_status)
                {
                    isvcd_nal_parse_pic_bound_proc(ps_nal_parse_ctxt, ps_vcl_nal, ps_nal_prms);
                    break;
                }
            }
            /* Perform the picture boundary detetction if all the  */
            /* following conditions are TRUE                       */
            /*  1. VCL NAL                                         */
            /*  2. Not a prefix NAL                                */
            /*  3. Not a discardable NAL                           */
            if((VCL_NAL == ps_nal_prms->i4_derived_nal_type) &&
               (PREFIX_UNIT_NAL != ps_nal_prms->i4_nal_unit_type) &&
               (SVCD_FALSE == ps_nal_parse_ctxt->i4_discard_nal_flag))
            {
                if(ANNEX_B == ps_nal_parse_ctxt->i4_input_bitstream_mode)
                {
                    ps_nal_parse_ctxt->u1_pic_boundary_aud_flag = 0;

                    i4_status = isvcd_detect_pic_boundary_annex_b(ps_nal_prms, pu1_stream_buffer,
                                                                  i4_cur_pos, &i4_pic_bound_status,
                                                                  ps_nal_parse_ctxt, pu4_num_bytes);
                }

                /* If picture boundary is detected then come out of  */
                /* the loop                                          */
                if(PIC_BOUNDARY_TRUE == i4_pic_bound_status)
                {
                    isvcd_nal_parse_pic_bound_proc(ps_nal_parse_ctxt, ps_vcl_nal, ps_nal_prms);
                    break;
                }
            }

            if(SVCD_FALSE == ps_nal_parse_ctxt->i4_discard_nal_flag)
            {
                /* Set the active NAL buffer structure and initialize */
                /* the nal buffer structure                           */
                isvcd_get_nal_buf(ps_nal_parse_ctxt, &ps_nal_buf);
                ps_nal_parse_ctxt->ps_nal_buf = ps_nal_buf;
            }
            else
            {
                ps_nal_parse_ctxt->ps_nal_buf = NULL;
            }
        }

        /*-------------------------------------------------------------------*/
        /* In RFC based bitstreams, this is a dummy update (in this mode, the*/
        /* bytes consumed updation is done by picture boundary dectection    */
        /* But for Annex B based streams this is valid update                */
        /*-------------------------------------------------------------------*/
        *pu4_bytes_consumed += u4_bytes_consumed_temp;

        /*********************************************************************/
        /*          EMULATION PREVENTION AND BYTE SWAPPING                   */
        /*********************************************************************/

        /* Determine output buffer */
        ps_nal_buf = ps_nal_parse_ctxt->ps_nal_buf;

        if(VCL_NAL == ps_nal_parse_ctxt->i4_nal_type)
        {
            ppu1_out_buf = &ps_nal_parse_ctxt->pu1_vcl_nal_buf;
            pu4_bytes_left = &ps_nal_parse_ctxt->u4_bytes_left_vcl;
        }
        else
        {
            ppu1_out_buf = &ps_nal_parse_ctxt->pu1_non_vcl_nal_buf;
            pu4_bytes_left = &ps_nal_parse_ctxt->u4_bytes_left_non_vcl;
        }

        /* if 0 bytes left then discard the current NAL */
        if(0 >= (WORD32) *pu4_bytes_left)
        {
            ps_nal_parse_ctxt->i4_discard_nal_flag = SVCD_TRUE;
        }

        /* Perform the emulation prevention and byte swap */
        if(SVCD_FALSE == ps_nal_parse_ctxt->i4_discard_nal_flag)
        {
            UWORD32 u4_output_bytes, u4_buf_inc;

            /* Do emulation prevention and byte swapping on all the packets  */
            /* of RFC or current partial or full Annex B NAL unit            */
            {
                UWORD32 u4_buf_size;

                /* clip the size before emulation prevention */
                u4_buf_size = (UWORD32) CLIP3(0, (WORD32) *pu4_bytes_left,
                                              (ps_nal_unit->i4_buf_sizes - i4_nal_header_len));

                u4_buf_inc = isvcd_nal_byte_swap_emulation(
                    (UWORD32 *) *ppu1_out_buf, &u4_output_bytes,
                    ps_nal_unit->pu1_bufs + i4_nal_header_len, u4_buf_size,
                    NUM_OF_ZERO_BYTES_BEFORE_START_CODE, &ps_nal_parse_ctxt->s_emulation_ctxt);

                i4_nal_header_len = 0;
                u4_buf_inc = UP_ALIGN_8(u4_buf_inc);
                *ppu1_out_buf += u4_buf_inc;
                *pu4_bytes_left -= u4_buf_inc;
                ps_nal_buf->i4_buf_size += u4_output_bytes;
            }
        }

        /*********************************************************************/
        /*                UPDATE VARIABLES                                   */
        /*********************************************************************/
        if(NAL_END == ps_nal_parse_ctxt->i4_find_nal_state)
        {
            if(SVCD_FALSE == ps_nal_parse_ctxt->i4_discard_nal_flag)
            {
                /* This fucntions updates output nal ctxt - vcl nal structure*/
                /* and non vcl nal structure depending upon the current NAL  */
                /* type.                                                     */
                /* This will only update parameters which are available at   */
                /* end of NAL unit like nal unit's total size                */
                isvcd_update_nal_ctxt(ps_nal_parse_ctxt, ps_vcl_nal, ps_non_vcl_nal);

                UPDATE_NAL_BUF_PTR(ppu1_out_buf, ps_nal_prms->i4_derived_nal_type, pu4_bytes_left);
            }

            /* If the prefix NAL unit is not immediatly followed by */
            /* a AVC NAL unit it shall be discarded and hence reset */
            /* is done                                              */
            /* Also if prefix NAL unit is discarded then we should  */
            /* not associate the prefix NAL unit with AVC NAL unit  */
            /* and hence a reset is required                        */
            if((PREFIX_UNIT_NAL != ps_nal_prms->i4_nal_unit_type) ||
               (SVCD_TRUE == ps_nal_parse_ctxt->i4_discard_nal_flag))
            {
                isvcd_nal_buf_reset(&ps_nal_parse_ctxt->s_prefix_nal_buf);
            }

            /* Reset the nal level tracking variables */
            isvcd_nal_reset_ctxt(ps_nal_parse_ctxt);
        }

        /*------------- while loop ends here --------------------------------*/
    } while(SVCD_TRUE == i4_more_data_flag);

    return (i4_pic_bound_status);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_nal_parse_non_vcl_nal                              */
/*                                                                           */
/*  Description   : None                                                     */
/*                                                                           */
/*  Inputs        : pv_nal_parse_ctxt - bitstream extract ctxt               */
/*                      structure                                            */
/*                  pv_input_bitstream_ctxt - bitstream context              */
/*                  pv_out_non_vcl_nal - non vcl nal structure (output)      */
/*                  pu4_bytes_consumed - bytes consumed variable(output)     */
/*                  pi4_num_packets_consumed - packets consumed (output/RFC) */
/*                                                                           */
/*  Globals       : None                                                     */
/*                                                                           */
/*  Processing    : None                                                     */
/*                                                                           */
/*  Outputs       : Updates bytes consumed variable, packets consumed,       */
/*                  output structures (non vcl nal)                          */
/*                                                                           */
/*  Returns       : If vcl nal is found then VCL_NAL_FOUND_TRUE otherwise    */
/*                  VCL_NAL_FOUND_FALSE                                      */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*          DD MM YYYY   Author(s)       Changes                             */
/*          06 09 2021   Vijay      Draft                                    */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_nal_parse_non_vcl_nal(void *pv_nal_parse_ctxt, UWORD8 *pu1_stream_buffer,
                                   void *pv_out_non_vcl_nal, UWORD32 *pu4_bytes_consumed,
                                   UWORD32 *pu4_num_bytes)
{
    /*! - Search for the NAL boundary
        - If NAL boundary is not found and bytes consumed is lesser than
          minimum buffer size then break out of the loop
        - if it is start of NAL then read the NAL header
        - If it is a VCL NAL then return from this fucntion saying that
          VCL NAL found
        - NAL discard logic determines whther the current NAL has to be
          discarded or not
        - If NAL is not discarded then populate the vcl or non vcl output
          structures
    */

    nal_parse_ctxt_t *ps_nal_parse_ctxt;
    non_vcl_nal_t *ps_non_vcl_nal;
    nal_unit_t *ps_nal_unit;
    WORD32 i4_nal_start_flag, i4_cur_pos, i4_status;
    WORD32 i4_nal_header_len, i4_more_data_flag;
    UWORD32 u4_bytes_consumed_temp = 0;
    UWORD8 **ppu1_out_buf;
    nal_prms_t *ps_nal_prms;

    ps_nal_parse_ctxt = (nal_parse_ctxt_t *) pv_nal_parse_ctxt;
    ps_non_vcl_nal = (non_vcl_nal_t *) pv_out_non_vcl_nal;
    ps_nal_unit = (nal_unit_t *) ps_nal_parse_ctxt->pv_nal_unit;
    ps_nal_prms = &ps_nal_parse_ctxt->s_nal_prms;

    /* Initialization */
    i4_cur_pos = 0;
    *pu4_bytes_consumed = 0;
    i4_nal_header_len = 0;
    i4_nal_start_flag = SVCD_FALSE;
    i4_more_data_flag = SVCD_TRUE;
    i4_status = PIC_BOUNDARY_FALSE;

    /* reset the target layer update flag */
    ps_nal_parse_ctxt->i4_tgt_lyr_update = SVCD_FALSE;
    /*************************************************************************/
    /*              SEARCHING FOR THE START OF BITSTREAM                     */
    /*************************************************************************/

    /*-----------------------------------------------------------------------*/
    /* For Annex B based bitstreams the first start code has to decoded      */
    /* The first start code can come after multiple process call also. This  */
    /* has to be carefully handled                                           */
    /*-----------------------------------------------------------------------*/

    if(ANNEX_B == ps_nal_parse_ctxt->i4_input_bitstream_mode &&
       SVCD_TRUE == ps_nal_parse_ctxt->i4_dec_frst_sc_flag)
    {
        WORD32 i4_status;

        i4_status =
            isvcd_get_first_start_code(pu1_stream_buffer, pu4_bytes_consumed, pu4_num_bytes);

        /*-------------------------------------------------------------------*/
        /* If start code found then proceed with bitstream extraction        */
        /*-------------------------------------------------------------------*/

        if(i4_status == SC_NOT_FOUND)
        {
            return (VCL_NAL_FOUND_FALSE);
        }

        i4_cur_pos = *pu4_bytes_consumed;
        ps_nal_parse_ctxt->i4_dec_frst_sc_flag = SVCD_FALSE;
    }

    /* If number of bytes left in the previous process call  */
    /* is is greater or equal to number of bytes in input    */
    /* buffer of the current process call then declare that  */
    /* end of bitstream has occurred and consume the bytes   */
    /* but do not decode                                     */
    if(ps_nal_parse_ctxt->u4_bytes_left >= (UWORD32) *pu4_num_bytes)
    {
        ps_nal_parse_ctxt->i4_discard_nal_flag = SVCD_TRUE;
        *pu4_bytes_consumed = *pu4_num_bytes;

        i4_status = isvcd_nal_parse_partial_signal_eos(ps_nal_parse_ctxt, NULL, ps_non_vcl_nal);
        return (i4_status);
    }

    do
    {
        nal_buf_t *ps_nal_buf;
        UWORD32 *pu4_bytes_left;

        /*********************************************************************/
        /*                  NAL BOUNDARY DETECTION                           */
        /*********************************************************************/
        /*-------------------------------------------------------------------*/
        /* Detect NAL boundary                                               */
        /* After return,  this NAL boundary detetction logic might be in     */
        /* one of following states:                                          */
        /*  - NAL_START                                                      */
        /*  - FIND_NAL_END                                                   */
        /*  - NAL_END                                                        */
        /*-------------------------------------------------------------------*/
        if(ANNEX_B == ps_nal_parse_ctxt->i4_input_bitstream_mode)
        {
            i4_nal_start_flag = isvcd_get_annex_b_nal_unit(
                pu1_stream_buffer, i4_cur_pos, *pu4_num_bytes,
                &ps_nal_parse_ctxt->i4_find_nal_state, &ps_nal_parse_ctxt->i4_zero_byte_cnt,
                &u4_bytes_consumed_temp, ps_nal_parse_ctxt->pv_nal_unit, &i4_more_data_flag);

            i4_cur_pos += u4_bytes_consumed_temp;
        }

        /* If current NAL unit is start of new NAL unit then parse the NAL
            header. If the current NAL unit type is VCL NAL then return from
            this function. otherwise apply NAL discard logic and discard the
            NAL if discard NAL flag is true                                  */

        if(SVCD_TRUE == i4_nal_start_flag)
        {
            UWORD32 u4_err_code;
            WORD32 i4_sps_pps_corrupt_status;

            /* Get the NAL prms. This involves the following things*/
            /* 1. Decode the NAL header                            */
            /* 2. Set the discard flag                             */
            /* 3. Decode the slice header if needed                */
            isvcd_get_nal_prms(ps_nal_unit->pu1_bufs, ps_nal_unit->i4_buf_sizes, ps_nal_prms,
                               &ps_nal_parse_ctxt->s_prefix_nal_prms,
                               &ps_nal_parse_ctxt->s_prefix_nal_buf, &u4_err_code,
                               &i4_sps_pps_corrupt_status, &ps_nal_parse_ctxt->i4_discard_nal_flag,
                               ps_nal_parse_ctxt);
            /* If the error code returned by the "picture boundary" */
            /* detetction is                                        */
            /* 1. Insufficient bitstream size: then store the bytes */
            /*    left and break out of the loop                    */
            /* 2. Corrupted slice: then discard the slice           */
            if((NAL_INSUFFICIENT_DATA == (WORD32) u4_err_code) &&
               (NAL_END != ps_nal_parse_ctxt->i4_find_nal_state))
            {
                ps_nal_parse_ctxt->u4_bytes_left = *pu4_num_bytes - *pu4_bytes_consumed;

                /* Reset the NAL level tracking variables */
                isvcd_nal_reset_ctxt(ps_nal_parse_ctxt);
                break;
            }
            else if(0 != u4_err_code)
            {
                ps_nal_parse_ctxt->i4_discard_nal_flag = SVCD_TRUE;
            }

            /* Populate other paramters based on the nal prms */
            ps_nal_parse_ctxt->i4_nal_type = ps_nal_prms->i4_derived_nal_type;
            i4_nal_header_len = ps_nal_prms->i4_nal_header_len;

            /* If derived NAL unit is VCL_NAL then return from this function */
            if(VCL_NAL == ps_nal_prms->i4_derived_nal_type &&
               PREFIX_UNIT_NAL != ps_nal_prms->i4_nal_unit_type)
            {
                isvcd_pic_reset_ctxt(ps_nal_parse_ctxt);

                return (VCL_NAL_FOUND_TRUE);
            }

            /* Set the active NAL buffer structure and initialize */
            /* the nal buffer structure                           */
            isvcd_get_nal_buf(ps_nal_parse_ctxt, &ps_nal_buf);

            ps_nal_parse_ctxt->ps_nal_buf = ps_nal_buf;
        }

        /* Update the bytes consumed variable */

        *pu4_bytes_consumed += u4_bytes_consumed_temp;

        ps_nal_buf = ps_nal_parse_ctxt->ps_nal_buf;
        if(VCL_NAL == ps_nal_parse_ctxt->i4_nal_type)
        {
            ppu1_out_buf = &ps_nal_parse_ctxt->pu1_vcl_nal_buf;
            pu4_bytes_left = &ps_nal_parse_ctxt->u4_bytes_left_vcl;
        }
        else
        {
            ppu1_out_buf = &ps_nal_parse_ctxt->pu1_non_vcl_nal_buf;
            pu4_bytes_left = &ps_nal_parse_ctxt->u4_bytes_left_non_vcl;
        }

        /* if 0 bytes left then discard the current NAL */
        if(0 >= (WORD32) *pu4_bytes_left)
        {
            ps_nal_parse_ctxt->i4_discard_nal_flag = SVCD_TRUE;
        }

        /* If NAL is not discarded then :
            1) Perform emulation prevention and byte swapping on the RBSP data
            2) Update the NAL unit ctxt:
                a) If VCL NAL then update DQID list
                b) If NON VCL NAL then update the non vcl output structure   */

        if(SVCD_FALSE == ps_nal_parse_ctxt->i4_discard_nal_flag)
        {
            UWORD32 u4_output_bytes, u4_buf_inc;

            {
                UWORD32 u4_buf_size;

                /* clip the size before emulation prevention */
                u4_buf_size = (UWORD32) CLIP3(0, (WORD32) *pu4_bytes_left,
                                              (ps_nal_unit->i4_buf_sizes - i4_nal_header_len));

                u4_buf_inc = isvcd_nal_byte_swap_emulation(
                    (UWORD32 *) *ppu1_out_buf, &u4_output_bytes,
                    ps_nal_unit->pu1_bufs + i4_nal_header_len, u4_buf_size,
                    NUM_OF_ZERO_BYTES_BEFORE_START_CODE, &ps_nal_parse_ctxt->s_emulation_ctxt);
                i4_nal_header_len = 0;

                u4_buf_inc = UP_ALIGN_8(u4_buf_inc);
                *ppu1_out_buf += u4_buf_inc;
                *pu4_bytes_left -= u4_buf_inc;
                ps_nal_buf->i4_buf_size += u4_output_bytes;
            }
        }

        /*********************************************************************/
        /*                UPDATE VARIABLES                                   */
        /*********************************************************************/

        if(NAL_END == ps_nal_parse_ctxt->i4_find_nal_state)
        {
            /*---------------------------------------------------------------*/
            /* - Update the total bits in the NAL. While doing so bits       */
            /* calculated so far should be converted to SODB length          */
            /*---------------------------------------------------------------*/
            if(SVCD_FALSE == ps_nal_parse_ctxt->i4_discard_nal_flag)
            {
                isvcd_update_nal_ctxt(ps_nal_parse_ctxt, NULL, ps_non_vcl_nal);

                UPDATE_NAL_BUF_PTR(ppu1_out_buf, ps_nal_prms->i4_derived_nal_type, pu4_bytes_left);
            }

            /* If the prefix NAL unit is not immediatly followed by */
            /* a AVC NAL unit it shall be discarded and hence reset */
            /* is done                                              */
            /* Also if prefix NAL unit is discarded then we should  */
            /* not associate the prefix NAL unit with AVC NAL unit  */
            /* and hence a reset is required                        */
            if((PREFIX_UNIT_NAL != ps_nal_prms->i4_nal_unit_type) ||
               (SVCD_TRUE == ps_nal_parse_ctxt->i4_discard_nal_flag))
            {
                isvcd_nal_buf_reset(&ps_nal_parse_ctxt->s_prefix_nal_buf);
            }

            /* Reset NAL level tracking variables */
            isvcd_nal_reset_ctxt(ps_nal_parse_ctxt);
        }

        i4_nal_header_len = 0;
        /*------------- while loop ends here --------------------------------*/
    } while(SVCD_TRUE == i4_more_data_flag);

    if(i4_more_data_flag == 0)
    {
        isvcd_pic_reset_ctxt(ps_nal_parse_ctxt);
        return (VCL_NAL_FOUND_TRUE);
    }

    return (VCL_NAL_FOUND_FALSE);
}