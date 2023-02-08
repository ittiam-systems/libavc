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
 *  isvcd_vui.c
 *
 * @brief
 *  This file contains routines to parse VUI NAL's
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_vui.h"
#include "ih264d_bitstrm.h"
#include "ih264d_parse_cavlc.h"
#include "isvcd_structs.h"
#include "ih264d_error_handler.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_vui_ext_parametres                           */
/*                                                                           */
/*  Description   : This function parses VUI NALs.                           */
/*  Inputs        : ps_vu4          pointer to VUI params                    */
/*                  ps_bitstrm   Bitstream                                   */
/*  Globals       : None                                                     */
/*  Processing    : Parses VUI NAL's units and stores the info               */
/*  Outputs       : None                                                     */
/*  Returns       : None                                                     */
/*                                                                           */
/*  Issues        : None                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 05 2002   Kishore             Draft                            */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_parse_vui_ext_parametres(svc_vui_ext_t *ps_svc_vui_ext, dec_bit_stream_t *ps_bitstrm)
{
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    WORD32 ret;
    UWORD32 u4_i;

    ps_svc_vui_ext->u4_vui_ext_num_entries_minus1 = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
    if(ps_svc_vui_ext->u4_vui_ext_num_entries_minus1 > 1023)
    {
        return ERROR_INV_SPS_PPS_T;
    }

    for(u4_i = 0; u4_i <= ps_svc_vui_ext->u4_vui_ext_num_entries_minus1; u4_i++)
    {
        ps_svc_vui_ext->u1_vui_ext_dependency_id[u4_i] = ih264d_get_bits_h264(ps_bitstrm, 3);
        ps_svc_vui_ext->u1_vui_ext_quality_id[u4_i] = ih264d_get_bits_h264(ps_bitstrm, 4);
        ps_svc_vui_ext->u1_vui_ext_temporal_id[u4_i] = ih264d_get_bits_h264(ps_bitstrm, 3);
        ps_svc_vui_ext->u1_vui_ext_timing_info_present_flag[u4_i] = ih264d_get_bit_h264(ps_bitstrm);

        if(1 == ps_svc_vui_ext->u1_vui_ext_timing_info_present_flag[u4_i])
        {
            ps_svc_vui_ext->u4_vui_ext_num_units_in_tick[u4_i] =
                ih264d_get_bits_h264(ps_bitstrm, 32);
            ps_svc_vui_ext->u4_vui_ext_time_scale[u4_i] = ih264d_get_bits_h264(ps_bitstrm, 32);
            ps_svc_vui_ext->u1_vui_ext_fixed_frame_rate_flag[u4_i] =
                ih264d_get_bit_h264(ps_bitstrm);
        }

        ps_svc_vui_ext->u1_vui_ext_nal_hrd_params_present_flag[u4_i] =
            ih264d_get_bit_h264(ps_bitstrm);
        if(ps_svc_vui_ext->u1_vui_ext_nal_hrd_params_present_flag[u4_i])
        {
            ret = ih264d_parse_hrd_parametres(&ps_svc_vui_ext->s_nal_hrd[u4_i], ps_bitstrm);
            if(ret != OK) return ret;
        }
        ps_svc_vui_ext->u1_vui_ext_vcl_hrd_params_present_flag[u4_i] =
            ih264d_get_bit_h264(ps_bitstrm);
        if(ps_svc_vui_ext->u1_vui_ext_vcl_hrd_params_present_flag[u4_i])
        {
            ret = ih264d_parse_hrd_parametres(&ps_svc_vui_ext->s_vcl_hrd[u4_i], ps_bitstrm);
            if(ret != OK) return ret;
        }
        if(ps_svc_vui_ext->u1_vui_ext_nal_hrd_params_present_flag[u4_i] ||
           ps_svc_vui_ext->u1_vui_ext_vcl_hrd_params_present_flag[u4_i])
        {
            ps_svc_vui_ext->u1_vui_ext_low_delay_hrd_flag[u4_i] = ih264d_get_bit_h264(ps_bitstrm);
        }
        ps_svc_vui_ext->u1_vui_ext_pic_struct_present_flag[u4_i] = ih264d_get_bit_h264(ps_bitstrm);
    }
    return OK;
}
