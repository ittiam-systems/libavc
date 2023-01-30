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
/* Includes */
/*****************************************************************************/

/* System include files */
#include "stdio.h"

/* User include files */
#include "irc_datatypes.h"
#include "irc_common.h"
#include "irc_cntrl_param.h"
#include "irc_mem_req_and_acq.h"
#include "irc_rd_model.h"
#include "irc_est_sad.h"
#include "irc_fixed_point_error_bits.h"
#include "irc_vbr_storage_vbv.h"
#include "irc_picture_type.h"
#include "irc_bit_allocation.h"
#include "irc_mb_model_based.h"
#include "irc_cbr_buffer_control.h"
#include "irc_vbr_str_prms.h"
#include "irc_rate_control_api.h"
#include "irc_rate_control_api_structs.h"
#include "irc_trace_support.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define DEV_Q 4        /*Q format(Shift) for Deviation range factor */
#define HI_DEV_FCTR 22 /* 1.4*16 */
#define LO_DEV_FCTR 12 /* 0.75*16 */
#define GET_HI_DEV_QP(Qprev) ((((WORD32) Qprev) * HI_DEV_FCTR + (1 << (DEV_Q - 1))) >> DEV_Q)
#define GET_LO_DEV_QP(Qprev) ((((WORD32) Qprev) * LO_DEV_FCTR + (1 << (DEV_Q - 1))) >> DEV_Q)
#define CLIP_QP(Qc, hi_d, lo_d) (((Qc) < (lo_d)) ? ((lo_d)) : (((Qc) > (hi_d)) ? (hi_d) : (Qc)))

/*******************************************************************************
 *  Description   : Gets the frame level qp for the given picture type
 *                  based on bits per pixel and gradient per pixel
 ******************************************************************************/
/* Get frame level QP based on BPP and GPP */
UWORD8 irc_get_frame_level_init_qp(rate_control_handle *ps_rate_control_api, rc_type_e e_rc_type,
                                   picture_type_e e_pic_type, DOUBLE d_bpp, DOUBLE d_gpp)
{
    DOUBLE d_frame_qp;

    UWORD8 u1_min_qp =
        ((rate_control_api_t *) (ps_rate_control_api))->au1_min_max_avc_qp[(e_pic_type << 1)];
    UWORD8 u1_max_qp =
        ((rate_control_api_t *) (ps_rate_control_api))->au1_min_max_avc_qp[(e_pic_type << 1) + 1];

    if((e_rc_type != VBR_STORAGE) && (e_rc_type != VBR_STORAGE_DVD_COMP) &&
       (e_rc_type != CBR_NLDRC) && (e_rc_type != CONST_QP) && (e_rc_type != VBR_STREAMING))
    {
        trace_printf(
            (const WORD8 *) (const WORD8 *) " Only VBR,NLDRC and CONST QP supported for now \n");
        return (0);
    }

    if(d_bpp <= 0.18)
    {
        d_frame_qp = 43.49 + (0.59 * d_gpp) - (106.45 * d_bpp);
    }
    else if(d_bpp <= 0.6)
    {
        d_frame_qp = 25.12 + (0.69 * d_gpp) - (29.23 * (d_bpp - 0.18));
    }
    else
    {
        d_frame_qp = 13.93 + (0.74 * d_gpp) - (18.4 * (d_bpp - 0.6));
    }

    /* Truncating the QP to the Max and Min Qp values possible */
    if(d_frame_qp < u1_min_qp) d_frame_qp = u1_min_qp;
    if(d_frame_qp > u1_max_qp) d_frame_qp = u1_max_qp;

    return ((UWORD8) (d_frame_qp + 0.5));
}

void irc_change_qp_constraints(rate_control_api_t *ps_rate_control_api, UWORD8 *pu1_min_max_qp,
                               UWORD8 *pu1_min_max_avc_qp)
{
    WORD32 i;

    for(i = 0; i < MAX_PIC_TYPE; i++)
    {
        ps_rate_control_api->au1_min_max_qp[(i << 1)] = pu1_min_max_qp[(i << 1)];
        ps_rate_control_api->au1_min_max_qp[(i << 1) + 1] = pu1_min_max_qp[(i << 1) + 1];
        ps_rate_control_api->au1_min_max_avc_qp[(i << 1)] = pu1_min_max_avc_qp[(i << 1)];
        ps_rate_control_api->au1_min_max_avc_qp[(i << 1) + 1] = pu1_min_max_avc_qp[(i << 1) + 1];
    }
}

UWORD8 irc_is_scenecut(rate_control_api_t *ps_rate_control_api)
{
    return ((rate_control_api_t *) (ps_rate_control_api))->u1_scd_detected;
}
