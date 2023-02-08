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
 *  isvcd_ii_pred.c
 *
 * @brief
 *  Contains routines that resample for SVC resampling
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_ii_pred_res_init()
 *  - isvcd_ii_get_ref_mb_mode()
 *  - isvcd_ii_get_ref_projections()
 *  - isvcd_ii_pred_compute_flags_mb()
 *  - isvcd_ii_pred_mb()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <assert.h>
#include <string.h>
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264_defs.h"
#include "ih264d_bitstrm.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "isvcd_structs.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_mb_utils.h"
#include "ih264d_deblocking.h"
#include "ih264d_dpb_manager.h"
#include "ih264d_mvpred.h"
#include "ih264d_inter_pred.h"
#include "ih264d_process_pslice.h"
#include "ih264d_error_handler.h"
#include "ih264d_cabac.h"
#include "ih264d_tables.h"
#include "ih264d_parse_slice.h"
#include "ih264d_utils.h"
#include "ih264d_parse_islice.h"
#include "ih264d_process_bslice.h"
#include "ih264d_process_intra_mb.h"
#include "isvcd_mode_mv_resamp.h"
#include "isvcd_ii_pred.h"
#include "ih264_debug.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_ii_pred_res_init                                    */
/*                                                                           */
/*  Description   : this function initialises the resolution level params    */
/*                  into the context structure                               */
/*                                                                           */
/*  Inputs        : pv_ii_pred_ctxt: Intra inter pred  handle                */
/*                  pi2_ref_loc_x             : pointer to buffer having the */
/*                                              projected locations horz     */
/*                  pi2_ref_loc_y             : pointer to buffer having the */
/*                                              projected location vertical  */
/*  Globals       : none                                                     */
/*  Processing    :                                                          */
/*                                                                           */
/*  Outputs       : none                                                     */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_ii_pred_res_init(void *pv_svc_dec)
{
    /* local vaiables */
    intra_inter_pred_ctxt_t *ps_ii_pred_ctxt;
    mode_motion_ctxt_t *ps_ctxt;
    mode_motion_lyr_ctxt *ps_lyr_mem;
    WORD32 i4_base_res_flag;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) pv_svc_dec;

    res_prms_t *ps_res_prms = &ps_svc_lyr_dec->s_res_prms;
    ps_ii_pred_ctxt = (intra_inter_pred_ctxt_t *) ps_svc_lyr_dec->pv_ii_pred_ctxt;
    ps_ctxt = (mode_motion_ctxt_t *) ps_svc_lyr_dec->pv_mode_mv_sample_ctxt;
    i4_base_res_flag = ps_svc_lyr_dec->u1_base_res_flag;

    if((0 != ps_svc_lyr_dec->u1_layer_id) && (SVCD_FALSE == i4_base_res_flag))
    {
        /* if not first resolution layer */
        ps_ii_pred_ctxt->i4_ref_res_lyr_wd = ps_ii_pred_ctxt->i4_cur_res_lyr_wd;
        ps_ii_pred_ctxt->i4_ref_res_lyr_ht = ps_ii_pred_ctxt->i4_cur_res_lyr_ht;
    }

    ps_lyr_mem = &ps_ctxt->as_res_lyr_mem[ps_ctxt->i4_res_id];

    ps_ii_pred_ctxt->pi2_ref_loc_x = ps_lyr_mem->pi2_ref_loc_x;
    ps_ii_pred_ctxt->pi2_ref_loc_y = ps_lyr_mem->pi2_ref_loc_y;

    /* Store the dimensions */
    ps_ii_pred_ctxt->i4_cur_res_lyr_wd = ps_res_prms->i4_res_width;
    ps_ii_pred_ctxt->i4_cur_res_lyr_ht = ps_res_prms->i4_res_height;

    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_ii_get_ref_mb_mode                                  */
/*                                                                           */
/*  Description   : This function is used to find the mb type of the         */
/*                  corresponding MB in the reference layer is INTER or      */
/*                  INTRA                                                    */
/*  Inputs        : pu1_ref_mb_modes : ref mb modes buffer pointer           */
/*                  i4_ref_mode_stride : mb mode buffer stride               */
/*                  i4_x_ref : reference location X                          */
/*                  i4_y_ref : reference location Y                          */
/*  Globals       : none                                                     */
/*  Processing    : it derives the byte corresponding to reference MB and    */
/*                  and gets the mb type                                     */
/*  Outputs       : none                                                     */
/*  Returns       : SVCD_TRUE if INTRA MB else SVCD_FALSE                    */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_ii_get_ref_mb_mode(WORD8 *pi1_ref_mb_modes, WORD32 i4_ref_mode_stride,
                                WORD32 i4_ref_mode_size, WORD32 i4_x_ref, WORD32 i4_y_ref)
{
    WORD32 i4_mb_x, i4_mb_y;
    inter_lyr_mb_prms_t *ps_inter_lyr_mb_prms;
    WORD8 i1_mb_mode;

    i4_mb_x = (i4_x_ref >> MB_WIDTH_SHIFT);
    i4_mb_y = (i4_y_ref >> MB_HEIGHT_SHIFT);

    /* get the location of the byte which has the current mb mode */
    pi1_ref_mb_modes += (i4_mb_y * i4_ref_mode_stride * i4_ref_mode_size);
    pi1_ref_mb_modes += (i4_mb_x * i4_ref_mode_size);
    ps_inter_lyr_mb_prms = (inter_lyr_mb_prms_t *) pi1_ref_mb_modes;
    i1_mb_mode = ps_inter_lyr_mb_prms->i1_mb_mode;

    if(i1_mb_mode <= SVC_INTER_MB)
    {
        /* INTER */
        return (SVCD_FALSE);
    }
    else
    {
        /* INTRA */
        return (SVCD_TRUE);
    }
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_ii_get_ref_projections                              */
/*                                                                           */
/*  Description   : this function projects the corners of current MB and     */
/*                  finds out if any point is falling into an INTRA MB in    */
/*                  reference layer. it also calculates the intersection     */
/*                  point of MB boundaries in the projected region           */
/*  Inputs        : ps_ctxt : Intra Inter context pointer                    */
/*                  ps_ii_mb_ctxt : Curretn MB context pointer               */
/*                  ps_ref_mb_mode : reference MB mode buffer descriptor     */
/*                  i4_mb_x : MB_X of current MB                             */
/*                  i4_mb_y : MB_Y of current MB                             */
/*  Globals       : none                                                     */
/*  Processing    : it derives the intra status of the corners and calculates*/
/*                  the intersection point                                   */
/*  Outputs       : non                                                      */
/*  Returns       : SVCD_TRUE or SVCD_FALSE                                  */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_ii_get_ref_projections(intra_inter_pred_ctxt_t *ps_ctxt,
                                    intra_inter_mb_t *ps_ii_mb_ctxt, mem_element_t *ps_ref_mb_mode,
                                    WORD32 i4_mb_x, WORD32 i4_mb_y)
{
    WORD16 *pi2_ref_loc_x;
    WORD16 *pi2_ref_loc_y;
    WORD8 *pi1_ref_mb_mode;
    WORD32 i4_ref_mode_stride;
    WORD32 i4_element_size;
    WORD32 i4_ref_x, i4_ref_y;
    WORD32 i4_frame_x, i4_frame_y;
    WORD32 i4_flag;

    pi2_ref_loc_x = ps_ctxt->pi2_ref_loc_x;
    pi2_ref_loc_y = ps_ctxt->pi2_ref_loc_y;

    pi1_ref_mb_mode = (WORD8 *) ps_ref_mb_mode->pv_buffer;
    i4_ref_mode_stride = ps_ref_mb_mode->i4_num_element_stride;
    i4_element_size = ps_ref_mb_mode->i4_element_size;

    /* get the current MB frame positions */
    i4_frame_x = i4_mb_x << 4;
    i4_frame_y = i4_mb_y << 4;

    /* reset the flag */
    i4_flag = SVCD_FALSE;

    /* project the (0,0) of current MB and get the ref MB mode */
    i4_ref_x = pi2_ref_loc_x[i4_frame_x];
    i4_ref_y = pi2_ref_loc_y[i4_frame_y];

    if((i4_ref_x < ps_ctxt->i4_ref_res_lyr_wd) && (i4_ref_y < ps_ctxt->i4_ref_res_lyr_ht))
    {
        ps_ii_mb_ctxt->u1_top_left_intra_flag = isvcd_ii_get_ref_mb_mode(
            pi1_ref_mb_mode, i4_ref_mode_stride, i4_element_size, i4_ref_x, i4_ref_y);
    }
    else
    {
        /* If projection is outside the picture boundary */
        ps_ii_mb_ctxt->u1_top_left_intra_flag = SVCD_FALSE;
    }
    /* project the (15,0) of current MB and get the ref MB mode */
    i4_ref_x = pi2_ref_loc_x[i4_frame_x + 15];
    i4_ref_y = pi2_ref_loc_y[i4_frame_y];

    if((i4_ref_x < ps_ctxt->i4_ref_res_lyr_wd) && (i4_ref_y < ps_ctxt->i4_ref_res_lyr_ht))
    {
        ps_ii_mb_ctxt->u1_top_rt_intra_flag = isvcd_ii_get_ref_mb_mode(
            pi1_ref_mb_mode, i4_ref_mode_stride, i4_element_size, i4_ref_x, i4_ref_y);
    }
    else
    {
        ps_ii_mb_ctxt->u1_top_rt_intra_flag = SVCD_FALSE;
    }

    /* project the (0,15) of current MB and get the ref MB mode */
    i4_ref_x = pi2_ref_loc_x[i4_frame_x];
    i4_ref_y = pi2_ref_loc_y[i4_frame_y + 15];

    if((i4_ref_x < ps_ctxt->i4_ref_res_lyr_wd) && (i4_ref_y < ps_ctxt->i4_ref_res_lyr_ht))
    {
        ps_ii_mb_ctxt->u1_bot_left_intra_flag = isvcd_ii_get_ref_mb_mode(
            pi1_ref_mb_mode, i4_ref_mode_stride, i4_element_size, i4_ref_x, i4_ref_y);
    }
    else
    {
        ps_ii_mb_ctxt->u1_bot_left_intra_flag = SVCD_FALSE;
    }

    /* project the (15,15) of current MB and get the ref MB mode */
    i4_ref_x = pi2_ref_loc_x[i4_frame_x + 15];
    i4_ref_y = pi2_ref_loc_y[i4_frame_y + 15];

    if((i4_ref_x < ps_ctxt->i4_ref_res_lyr_wd) && (i4_ref_y < ps_ctxt->i4_ref_res_lyr_ht))
    {
        ps_ii_mb_ctxt->u1_bot_rt_intra_flag = isvcd_ii_get_ref_mb_mode(
            pi1_ref_mb_mode, i4_ref_mode_stride, i4_element_size, i4_ref_x, i4_ref_y);
    }
    else
    {
        ps_ii_mb_ctxt->u1_bot_rt_intra_flag = SVCD_FALSE;
    }

    /* if any of the 4 cormers are falling into intra region
      set the INTRA INTER Flag */
    if((SVCD_TRUE == ps_ii_mb_ctxt->u1_top_left_intra_flag) ||
       (SVCD_TRUE == ps_ii_mb_ctxt->u1_top_rt_intra_flag) ||
       (SVCD_TRUE == ps_ii_mb_ctxt->u1_bot_left_intra_flag) ||
       (SVCD_TRUE == ps_ii_mb_ctxt->u1_bot_rt_intra_flag))
    {
        i4_flag = SVCD_TRUE;
    }

    /* derive the intersection point of MB boundaries */
    if(SVCD_TRUE == i4_flag)
    {
        WORD32 i4_intr_x, i4_intr_y;
        WORD32 i4_ref_mb_init_x, i4_ref_mb_init_y;
        WORD32 i4_ctr;

        /* set the variables to initial values */
        i4_intr_x = 0;
        i4_intr_y = 0;
        i4_ref_mb_init_x = pi2_ref_loc_x[i4_frame_x] >> MB_WIDTH_SHIFT;
        i4_ref_mb_init_y = pi2_ref_loc_y[i4_frame_y] >> MB_HEIGHT_SHIFT;

        /* loop until an Mb boundary is found in horizontal direction */
        for(i4_ctr = 0; i4_ctr < MB_WIDTH; i4_ctr++)
        {
            i4_ref_x = pi2_ref_loc_x[i4_frame_x + i4_ctr];
            i4_ref_x >>= MB_WIDTH_SHIFT;

            /* check if the locations are falling into same MB */
            if(i4_ref_x != i4_ref_mb_init_x)
            {
                break;
            }
            /* increment the position */
            i4_intr_x++;
        }

        /* loop until an Mb boundary is found in vertical direction */
        for(i4_ctr = 0; i4_ctr < MB_HEIGHT; i4_ctr++)
        {
            i4_ref_y = pi2_ref_loc_y[i4_frame_y + i4_ctr];
            i4_ref_y >>= MB_HEIGHT_SHIFT;

            /* check if the locations are falling into same MB */
            if(i4_ref_y != i4_ref_mb_init_y)
            {
                break;
            }
            /* increment the position */
            i4_intr_y++;
        }
        /* store the intersection points */
        ps_ii_mb_ctxt->u1_intersection_x = i4_intr_x;
        ps_ii_mb_ctxt->u1_intersection_y = i4_intr_y;
    }
    else
    {
        /* set to default value */
        ps_ii_mb_ctxt->u1_intersection_x = 0;
        ps_ii_mb_ctxt->u1_intersection_y = 0;
    }

    return (i4_flag);
}
/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_ii_pred_compute_flags_mb                            */
/*                                                                           */
/*  Description   : this function checks all the criteria for an MB to       */
/*                  under go Inter-Intra prediction  and stores the MB mode  */
/*                   as INTER_INTRA for appropriate MBs                      */
/*  Inputs        : refer to comments below                                  */
/*  Globals       : none                                                     */
/*  Processing    : it checks the criteria for anMB to undergo Inter-Intra   */
/*                  pred process and updates the MB mode                     */
/*  Outputs       : MB mode set for each MB with INTRA-INTER status          */
/*  Returns       : SVCD_EOK or SVCD_EFAIL                                   */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_ii_pred_compute_flags_mb(void *pv_ii_pred_ctxt, mem_element_t *ps_ref_mb_mode,
                                      mb_coord_t *ps_coord, void *pv_mb_prms, void *pv_svc_mb_prms,
                                      UWORD8 *pu1_ii_mb_mode)
{
    intra_inter_pred_ctxt_t *ps_ctxt;
    WORD32 i4_mb_x, i4_mb_y;
    dec_svc_mb_info_t *ps_svc_mb_prms;
    UNUSED(pv_mb_prms);

    if((NULL == pv_ii_pred_ctxt) || (NULL == ps_ref_mb_mode) || (NULL == ps_coord) ||
       (NULL == pu1_ii_mb_mode))
    {
        return NOT_OK;
    }

    ps_ctxt = (intra_inter_pred_ctxt_t *) pv_ii_pred_ctxt;
    ps_svc_mb_prms = (dec_svc_mb_info_t *) pv_svc_mb_prms;

    /* get mb co-ordinates */
    i4_mb_x = ps_coord->u2_mb_x;
    i4_mb_y = ps_coord->u2_mb_y;

    {
        intra_inter_mb_t *ps_ii_mb_ctxt;
        WORD32 i4_ii_flag;

        /* get the current MB strcuture pointer */
        ps_ii_mb_ctxt = &ps_ctxt->s_intra_inter_mb_prms;

        /* reset the Intra Inter qualified flag for current MB */
        i4_ii_flag = SVCD_FALSE;

        /* check for base mode flag and Inter MB status */
        if(1 == ps_svc_mb_prms->u1_base_mode_flag)
        {
            /* call the function which calculates the projections
               and returns whether current MB has to under go
               Inter Intra Prediction */
            i4_ii_flag = isvcd_ii_get_ref_projections(ps_ctxt, ps_ii_mb_ctxt, ps_ref_mb_mode,
                                                      i4_mb_x, i4_mb_y);
        }

        /* If the current MB requires Intra Inter prediction */
        if(SVCD_TRUE == i4_ii_flag)
        {
            /* set the mb mode */
            *pu1_ii_mb_mode = SVC_INTRA_INTER_MB;
        }
        else
        {
            /* set all MB params to default values */
            ps_ii_mb_ctxt->u1_bot_left_intra_flag = SVCD_FALSE;
            ps_ii_mb_ctxt->u1_bot_rt_intra_flag = SVCD_FALSE;
            ps_ii_mb_ctxt->u1_top_left_intra_flag = SVCD_FALSE;
            ps_ii_mb_ctxt->u1_top_rt_intra_flag = SVCD_FALSE;
            ps_ii_mb_ctxt->u1_intersection_x = 0;
            ps_ii_mb_ctxt->u1_intersection_y = 0;

            /* set the mb mode to 0 (which has no interpretation) */
            *pu1_ii_mb_mode = 0;
        }
    }
    return (OK);
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_ii_pred_mb                                          */
/*                                                                           */
/*  Description   : This function performs the Intra-Inter Preduction of the */
/*                  given MB                                                 */
/*                                                                           */
/*  Inputs        : ps_mb_ctxt : Intra Inter mb context strcuture            */
/*                  ps_mb_buf : current MB buffers strcuture pointer         */
/*  Globals       : none                                                     */
/*  Processing    : it processes all partitions based on the Intra flag      */
/*                                                                           */
/*  Outputs       : Intra Inter Predecited and reconstructed MB              */
/*  Returns       : none                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Vijay                creation                        */
/*                                                                           */
/*****************************************************************************/
void isvcd_ii_pred_mb(void *pv_svc_dec, dec_mb_info_t *ps_cur_mb_info)
{
    intra_inter_mb_t *ps_mb_ctxt;
    UWORD8 *pu1_rec_y, *pu1_rec_uv;
    UWORD8 *pu1_recon_luma;
    WORD32 i4_recon_luma_stride;
    UWORD8 *pu1_recon_chroma;
    WORD32 i4_recon_chroma_stride;
    UWORD8 *pu1_pred_luma;
    UWORD8 *pu1_pred_chroma;
    WORD32 i4_pred_luma_stride;
    WORD32 i4_pred_chroma_stride;
    WORD32 i4_intr_x, i4_intr_y;
    intra_inter_pred_ctxt_t *ps_ctxt;
    pic_buffer_t *ps_frame_buf;
    svc_dec_lyr_struct_t *ps_svc_lyr_dec = (svc_dec_lyr_struct_t *) pv_svc_dec;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;

    ps_ctxt = (intra_inter_pred_ctxt_t *) ps_svc_lyr_dec->pv_ii_pred_ctxt;
    ps_mb_ctxt = &ps_ctxt->s_intra_inter_mb_prms;
    ps_frame_buf = ps_dec->ps_cur_pic;
    i4_recon_luma_stride = ps_dec->u2_frm_wd_y;
    i4_recon_chroma_stride = ps_dec->u2_frm_wd_uv;

    /* derive the intersection point */
    i4_intr_x = ps_mb_ctxt->u1_intersection_x;
    i4_intr_y = ps_mb_ctxt->u1_intersection_y;

    pu1_rec_y = ps_frame_buf->pu1_buf1 + (ps_cur_mb_info->u2_mbx << 4) +
                (i4_recon_luma_stride * (ps_cur_mb_info->u2_mby << 4));

    pu1_rec_uv = ps_frame_buf->pu1_buf2 + (ps_cur_mb_info->u2_mbx << 3) * YUV420SP_FACTOR +
                 (i4_recon_chroma_stride * (ps_cur_mb_info->u2_mby << 3));

    pu1_pred_luma = ps_svc_lyr_dec->pu1_ii_resamp_buffer_luma;
    pu1_pred_chroma = ps_svc_lyr_dec->pu1_ii_resamp_buffer_chroma;
    i4_pred_luma_stride = MB_SIZE;
    i4_pred_chroma_stride = MB_SIZE;

    /* get the recon and residual buffer pointer */
    pu1_recon_luma = pu1_rec_y;
    pu1_recon_chroma = pu1_rec_uv;

    /*-----------------------------------------------------------------------*/
    /* Reconstruct TOP_LEFT Partition                                        */
    /*-----------------------------------------------------------------------*/
    {
        WORD32 i4_width, i4_height;

        /* assign the appropriate buffer params based on Intra status */
        if(SVCD_TRUE == ps_mb_ctxt->u1_top_left_intra_flag)
        {
            /* Luma Processing */
            isvcd_copy_data(pu1_pred_luma, i4_pred_luma_stride, pu1_recon_luma,
                            i4_recon_luma_stride, i4_intr_x, i4_intr_y);

            /* assign appropriate width and height for chroma */
            i4_width = (((i4_intr_x + 1) >> 1) << 1);
            i4_height = ((i4_intr_y + 1) & ~1);
            i4_height >>= 1;
            /* Chroma Processing (cb and cr interleaved) */
            isvcd_copy_data(pu1_pred_chroma, i4_pred_chroma_stride, pu1_recon_chroma,
                            i4_recon_chroma_stride, i4_width, i4_height);
        }
    }

    /*-----------------------------------------------------------------------*/
    /* Reconstruct TOP_RIGHT Partition                                       */
    /*-----------------------------------------------------------------------*/
    {
        WORD32 i4_width, i4_height;

        /* assign the appropriate buffer params based on Intra status */
        if(SVCD_TRUE == ps_mb_ctxt->u1_top_rt_intra_flag)
        {
            pu1_pred_luma += i4_intr_x;
            pu1_pred_chroma += (((i4_intr_x + 1) >> 1) << 1);

            /* ----------------------- Luma ------------------------ */
            /* get the recon and residual buffer pointer */
            pu1_recon_luma = pu1_rec_y + i4_intr_x;

            /* assign appropriate width and height for luma */
            i4_width = MB_WIDTH - i4_intr_x;
            i4_height = i4_intr_y;

            /* Luma Processing */
            /* Luma Processing */
            isvcd_copy_data(pu1_pred_luma, i4_pred_luma_stride, pu1_recon_luma,
                            i4_recon_luma_stride, i4_width, i4_height);

            /* ----------------------- Chroma ----------------------- */
            /* assign appropriate width and height for luma */
            i4_width = (BLOCK_WIDTH - ((i4_intr_x + 1) >> 1)) << 1;

            /* Height includes for both Cb & Cr */
            i4_height = ((i4_intr_y + 1) & ~1);
            i4_height >>= 1;
            /* get the recon and residual buffer pointer */
            pu1_recon_chroma = pu1_rec_uv;
            {
                WORD32 i4_temp;
                i4_temp = (((i4_intr_x + 1) >> 1) << 1);
                pu1_recon_chroma += i4_temp;
            }

            /* Chroma Processing (cb and cr  interleaved) */
            isvcd_copy_data(pu1_pred_chroma, i4_pred_chroma_stride, pu1_recon_chroma,
                            i4_recon_chroma_stride, i4_width, i4_height);
        }
    }

    /*-----------------------------------------------------------------------*/
    /* Reconstruct BOTTOM_LEFT Partition                                     */
    /*-----------------------------------------------------------------------*/
    {
        WORD32 i4_width, i4_height;

        /* assign the appropriate buffer params based on Intra status */
        if(SVCD_TRUE == ps_mb_ctxt->u1_bot_left_intra_flag)
        {
            pu1_pred_luma = ps_svc_lyr_dec->pu1_ii_resamp_buffer_luma;
            pu1_pred_chroma = ps_svc_lyr_dec->pu1_ii_resamp_buffer_chroma;

            /* increment to current vertical offset */
            pu1_pred_luma += i4_intr_y * i4_pred_luma_stride;
            pu1_pred_chroma += (((i4_intr_y + 1) & ~1) >> 1) * i4_pred_chroma_stride;

            /* ----------------------- Luma ----------------------- */
            /* get the recon and residual buffer pointer */
            pu1_recon_luma = pu1_rec_y;
            pu1_recon_luma += i4_intr_y * i4_recon_luma_stride;

            /* assign appropriate width and height */
            i4_width = i4_intr_x;
            i4_height = MB_HEIGHT - i4_intr_y;

            /* Luma Processing */
            isvcd_copy_data(pu1_pred_luma, i4_pred_luma_stride, pu1_recon_luma,
                            i4_recon_luma_stride, i4_width, i4_height);

            /* ----------------------- Chroma ----------------------- */
            pu1_recon_chroma = pu1_rec_uv;
            {
                WORD32 i4_temp;
                i4_temp = ((i4_intr_y + 1) & ~1) >> 1;
                pu1_recon_chroma += (i4_temp * i4_recon_chroma_stride);
            }
            /* assign appropriate width and height */
            i4_width = ((i4_intr_x + 1) >> 1) << 1;
            i4_height = MB_HEIGHT - (i4_intr_y & ~1);
            i4_height >>= 1;
            /* Chroma Processing (cb and cr interleaved) */
            isvcd_copy_data(pu1_pred_chroma, i4_pred_chroma_stride, pu1_recon_chroma,
                            i4_recon_chroma_stride, i4_width, i4_height);
        }
    }

    /*-----------------------------------------------------------------------*/
    /* Reconstruct BOTTOM_RIGHT Partition                                    */
    /*-----------------------------------------------------------------------*/
    {
        WORD32 i4_width, i4_height;

        /* assign the appropriate buffer params based on Intra status */
        if(SVCD_TRUE == ps_mb_ctxt->u1_bot_rt_intra_flag)
        {
            pu1_pred_luma = ps_svc_lyr_dec->pu1_ii_resamp_buffer_luma;
            pu1_pred_chroma = ps_svc_lyr_dec->pu1_ii_resamp_buffer_chroma;

            /* increment to current vertical offset */
            pu1_pred_luma += i4_intr_x;
            pu1_pred_luma += i4_intr_y * i4_pred_luma_stride;
            pu1_pred_chroma += (((i4_intr_y + 1) & ~1) >> 1) * i4_pred_chroma_stride;
            pu1_pred_chroma += ((i4_intr_x + 1) >> 1) << 1;

            /* ----------------------- Luma ----------------------- */
            /* get the recon and residual buffer pointer horz */
            pu1_recon_luma = pu1_rec_y + i4_intr_x;

            /* get the recon and residual buffer pointer vertical */
            pu1_recon_luma += (i4_intr_y * i4_recon_luma_stride);

            /* assign appropriate width and height */
            i4_width = MB_WIDTH - i4_intr_x;
            i4_height = MB_HEIGHT - i4_intr_y;

            /* Luma Processing */
            isvcd_copy_data(pu1_pred_luma, i4_pred_luma_stride, pu1_recon_luma,
                            i4_recon_luma_stride, i4_width, i4_height);

            /* ----------------------- Chroma ----------------------- */
            /* get the recon and residual buffer pointer horz */
            pu1_recon_chroma = pu1_rec_uv;
            {
                WORD32 i4_temp;
                i4_temp = ((i4_intr_y + 1) & ~1) >> 1;
                i4_temp *= i4_recon_chroma_stride;
                i4_temp += (((i4_intr_x + 1) >> 1) << 1);
                pu1_recon_chroma += i4_temp;
            }

            /* assign appropriate width and height */
            i4_width = (BLOCK_WIDTH - ((i4_intr_x + 1) >> 1)) << 1;
            i4_height = MB_HEIGHT - (i4_intr_y & ~1);
            i4_height >>= 1;
            /* Chroma Processing (cb and cr interleaved) */
            isvcd_copy_data(pu1_pred_chroma, i4_pred_chroma_stride, pu1_recon_chroma,
                            i4_recon_chroma_stride, i4_width, i4_height);
        }
    }
    return;
}