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
*  isvce_rc_utils.c
*
* @brief
*  Contains get gpp function required by the SVC encoder
*
* @author
*  ittiam
*
* @par List of Functions:
*  - isvce_get_gpp()
*  - isvce_rc_utils_init()
*  - isvce_get_rc_utils_data_size()
*  - isvce_compute_gpp()
*  - isvce_get_gpp_function_selector()
*
* @remarks
*  None
*
*******************************************************************************
*/

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "isvc_structs.h"
#include "isvce_rc_utils.h"
#include "isvce_rc_utils_private_defs.h"

/**
*******************************************************************************
*
* @brief
*   get gpp function
*
* @par Description:
*   computes gradient per pixel value for a given frame
*
* @param[in] ps_input_buf
*  pointer to yuv buffer properties
*
* @returns
*  calculated gpp value
*
* @remarks
*  none
*
*******************************************************************************
*/

static DOUBLE isvce_get_gpp(yuv_buf_props_t *ps_input_buf)
{
    UWORD32 i, j;

    DOUBLE d_gpp_y = 0;
    DOUBLE d_gpp_u = 0;
    DOUBLE d_gpp_v = 0;

    DOUBLE d_gpp = 0;

    UWORD32 u4_width = ps_input_buf->u4_width;
    UWORD32 u4_height = ps_input_buf->u4_height;

    UWORD8 *pu1_input_buf = (UWORD8 *) ps_input_buf->as_component_bufs[0].pv_data;
    WORD32 i4_input_stride = ps_input_buf->as_component_bufs[0].i4_data_stride;

    for(i = 0; i < u4_height - 1; i++)
    {
        for(j = 0; j < u4_width - 1; j++)
        {
            UWORD8 u1_cur_pix = pu1_input_buf[j];
            UWORD8 u1_bot_pix = pu1_input_buf[i4_input_stride + j];
            UWORD8 u1_right_pix = pu1_input_buf[j + 1];

            d_gpp_y += (ABS(u1_cur_pix - u1_bot_pix) + ABS(u1_cur_pix - u1_right_pix));
        }
        pu1_input_buf += i4_input_stride;
    }

    pu1_input_buf = (UWORD8 *) ps_input_buf->as_component_bufs[1].pv_data;
    i4_input_stride = ps_input_buf->as_component_bufs[1].i4_data_stride;

    for(i = 0; i < (u4_height >> 1) - 1; i++)
    {
        for(j = 0; j < u4_width - 2; j += 2)
        {
            UWORD8 u1_cur_pix = pu1_input_buf[j];
            UWORD8 u1_bot_pix = pu1_input_buf[i4_input_stride + j];
            UWORD8 u1_right_pix = pu1_input_buf[j + 2];

            d_gpp_u += (ABS(u1_cur_pix - u1_bot_pix) + ABS(u1_cur_pix - u1_right_pix));

            u1_cur_pix = pu1_input_buf[j + 1];
            u1_bot_pix = pu1_input_buf[i4_input_stride + j + 1];
            u1_right_pix = pu1_input_buf[j + 2 + 1];

            d_gpp_v += (ABS(u1_cur_pix - u1_bot_pix) + ABS(u1_cur_pix - u1_right_pix));
        }
        pu1_input_buf += i4_input_stride;
    }

    d_gpp_y /= (u4_width * u4_height);
    d_gpp_u /= ((u4_width >> 1) * (u4_height >> 1));
    d_gpp_v /= ((u4_width >> 1) * (u4_height >> 1));

    d_gpp = (DOUBLE) ((4 * d_gpp_y) + d_gpp_u + d_gpp_v) / 6;

    return d_gpp;
}

/**
*******************************************************************************
*
* @brief
*   gets the memory size required for compute gpp
*
* @par Description:
*   returns the memory required by the rc utils context and state structs
*   for allocation.
*
* @returns
*
* @remarks
*
*
*******************************************************************************
*/

UWORD32 isvce_get_rc_utils_data_size() { return sizeof(svc_rc_utils_state_t); }

/**
*******************************************************************************
*
* @brief
*   compute gpp process
*
* @par Description:
*   calls the function to compute gpp
*
* @param[in] ps_svc_rc_utils_ctxt
*  pointer to svc rc utils context
*
* @param[in] ps_input_buf
*  pointer to yuv buffer properties
*
* @returns
*  calculated gpp value
*
* @remarks
*  none
*
*******************************************************************************
*/

DOUBLE isvce_compute_gpp(svc_rc_utils_ctxt_t *ps_svc_rc_utils_ctxt, yuv_buf_props_t *ps_input_buf)
{
    svc_rc_utils_state_t *ps_rc_utils_state =
        (svc_rc_utils_state_t *) ps_svc_rc_utils_ctxt->pv_rc_utils_state;

    return ps_rc_utils_state->pf_get_gpp(ps_input_buf);
}

/**
*******************************************************************************
*
* @brief
*   selects which function to call for get gpp based on e_arch
*
* @par Description:
*
* @param[in] ps_rc_utils_state
*  pointer to svc rc utils state
*
* @param[in] e_arch
*  architecure type
*
* @returns
*
* @remarks
*
*******************************************************************************
*/

static void isvce_get_gpp_function_selector(svc_rc_utils_state_t *ps_rc_utils_state,
                                            IV_ARCH_T e_arch)
{
    switch(e_arch)
    {
#if defined(X86)
        case ARCH_X86_SSE42:
        {
            ps_rc_utils_state->pf_get_gpp = isvce_get_gpp_sse42;

            break;
        }
#elif defined(ARMV8)
        case ARCH_ARM_A53:
        case ARCH_ARM_A57:
        case ARCH_ARM_V8_NEON:
        {
            ps_rc_utils_state->pf_get_gpp = isvce_get_gpp_neon;

            break;
        }
#elif !defined(DISABLE_NEON)
        case ARCH_ARM_A9Q:
        case ARCH_ARM_A9A:
        case ARCH_ARM_A9:
        case ARCH_ARM_A7:
        case ARCH_ARM_A5:
        case ARCH_ARM_A15:
        {
            ps_rc_utils_state->pf_get_gpp = isvce_get_gpp_neon;

            break;
        }
#endif
        default:
        {
            ps_rc_utils_state->pf_get_gpp = isvce_get_gpp;

            break;
        }
    }
}

/**
*******************************************************************************
*
* @brief
*   initializes the rc utils context
*
* @par Description:
*   initializes the rc utils context
*
* @param[in] ps_svc_rc_utils_ctxt
*   pointer to svc rc utils context
*
* @param[in] ps_mem_rec
*   pointer to memory allocated to compute gpp process
*
* @param[in] e_arch
*   architecure type
*
* @returns
*
* @remarks
*  none
*
*******************************************************************************
*/

void isvce_rc_utils_init(svc_rc_utils_ctxt_t *ps_svc_rc_utils_ctxt, iv_mem_rec_t *ps_mem_rec,
                         IV_ARCH_T e_arch)
{
    svc_rc_utils_state_t *ps_rc_utils_state;

    UWORD8 *pu1_buf = (UWORD8 *) ps_mem_rec->pv_base;

    ps_rc_utils_state = (svc_rc_utils_state_t *) pu1_buf;

    ps_svc_rc_utils_ctxt->pv_rc_utils_state = ps_rc_utils_state;

    isvce_get_gpp_function_selector(ps_rc_utils_state, e_arch);
}
