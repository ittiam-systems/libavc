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
 *  isvcd_parse_eislice.c
 *
 * @brief
 *  Contains routines that decode a EI slice type
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_parse_islice_data_cabac()
 *  - isvcd_parse_islice_data_cavlc()
 *  - isvcd_parse_imb_cavlc()
 *  - isvcd_parse_eislice()
 *  - isvcd_parse_eislice_data_cabac()
 *  - isvcd_parse_eislice_data_cavlc()
 *  - isvcd_parse_imb_cabac()
 *  - isvcd_parse_islice()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <string.h>
#include "ih264_defs.h"
#include "ih264d_error_handler.h"
#include "ih264d_debug.h"
#include "ih264d_bitstrm.h"
#include "ih264d_defs.h"
#include "ih264d_tables.h"
#include "isvcd_structs.h"
#include "ih264d_parse_cavlc.h"
#include "ih264d_mb_utils.h"
#include "ih264d_deblocking.h"
#include "ih264d_cabac.h"
#include "ih264d_parse_cabac.h"
#include "ih264d_parse_mb_header.h"
#include "ih264d_parse_slice.h"
#include "ih264d_process_pslice.h"
#include "isvcd_process_epslice.h"
#include "ih264d_process_intra_mb.h"
#include "ih264d_parse_islice.h"
#include "ih264d_error_handler.h"
#include "ih264d_mvpred.h"
#include "ih264d_thread_parse_decode.h"
#include "ithread.h"
#include "ih264d_parse_mb_header.h"
#include "assert.h"
#include "ih264d_utils.h"
#include "ih264d_format_conv.h"
#include "ih264d_parse_headers.h"
#include "isvcd_parse_headers.h"
#include "isvcd_parse_slice.h"
#include "isvcd_mb_utils.h"

void ih264d_init_cabac_contexts(UWORD8 u1_slice_type, dec_struct_t *ps_dec);
void isvcd_init_cabac_contexts(UWORD8 u1_slice_type, dec_struct_t *ps_dec);

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_islice_data_cabac                            */
/*                                                                           */
/*  Description   : This function parses cabac syntax of a inter slice on    */
/*                  N MB basis.                                              */
/*                                                                           */
/*  Inputs        : ps_dec                                                   */
/*                  sliceparams                                              */
/*                  firstMbInSlice                                           */
/*                                                                           */
/*  Processing    : 1. After parsing syntax for N MBs those N MBs are        */
/*                     decoded till the end of slice.                        */
/*                                                                           */
/*  Returns       : 0                                                        */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         24 06 2005   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_parse_islice_data_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                     dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice)
{
    UWORD8 uc_more_data_flag;
    UWORD8 u1_num_mbs, u1_mb_idx;
    dec_mb_info_t *ps_cur_mb_info;
    deblk_mb_t *ps_cur_deblk_mb;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_bit_stream_t *const ps_bitstrm = ps_dec->ps_bitstrm;
    UWORD16 i2_pic_wdin_mbs = ps_dec->u2_frm_wd_in_mbs;
    WORD16 i2_cur_mb_addr;
    UWORD8 u1_mbaff;
    UWORD8 u1_num_mbs_next, u1_end_of_row, u1_tfr_n_mb;
    WORD32 ret = OK;

    ps_dec->u1_qp = ps_slice->u1_slice_qp;
    ih264d_update_qp(ps_dec, 0);
    u1_mbaff = ps_slice->u1_mbaff_frame_flag;

    if(ps_bitstrm->u4_ofst & 0x07)
    {
        ps_bitstrm->u4_ofst += 8;
        ps_bitstrm->u4_ofst &= 0xFFFFFFF8;
    }
    ret = ih264d_init_cabac_dec_envirnoment(&(ps_dec->s_cab_dec_env), ps_bitstrm);
    if(ret != OK) return ret;
    ih264d_init_cabac_contexts(I_SLICE, ps_dec);

    ps_dec->i1_prev_mb_qp_delta = 0;

    /* initializations */
    u1_mb_idx = ps_dec->u1_mb_idx;
    u1_num_mbs = u1_mb_idx;
    uc_more_data_flag = 1;
    i2_cur_mb_addr = u2_first_mb_in_slice << u1_mbaff;
    do
    {
        UWORD16 u2_mbx;
        ps_dec->pv_prev_mb_parse_tu_coeff_data = ps_dec->pv_parse_tu_coeff_data;

        if(i2_cur_mb_addr > ps_dec->ps_cur_sps->u2_max_mb_addr)
        {
            break;
        }

        {
            UWORD8 u1_mb_type;
            ps_cur_mb_info = ps_dec->ps_nmb_info + u1_num_mbs;
            ps_dec->u4_num_mbs_cur_nmb = u1_num_mbs;
            ps_dec->u4_num_pmbair = (u1_num_mbs >> u1_mbaff);
            ps_cur_mb_info->u1_end_of_slice = 0;

            /***************************************************************/
            /* Get the required information for decoding of MB                  */
            /* mb_x, mb_y , neighbour availablity,                              */
            /***************************************************************/
            ps_dec->pf_get_mb_info(ps_dec, i2_cur_mb_addr, ps_cur_mb_info, 0);
            u2_mbx = ps_dec->u2_mbx;

            /*********************************************************************/
            /* initialize u1_tran_form8x8 to zero to aviod uninitialized accesses */
            /*********************************************************************/
            ps_cur_mb_info->u1_tran_form8x8 = 0;
            ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;

            /***************************************************************/
            /* Set the deblocking parameters for this MB                   */
            /***************************************************************/
            ps_cur_deblk_mb = ps_dec->ps_deblk_mbn + u1_num_mbs;
            if(ps_dec->u4_app_disable_deblk_frm == 0)
                ih264d_set_deblocking_parameters(ps_cur_deblk_mb, ps_slice,
                                                 ps_dec->u1_mb_ngbr_availablity,
                                                 ps_dec->u1_cur_mb_fld_dec_flag);

            ps_cur_deblk_mb->u1_mb_type = ps_cur_deblk_mb->u1_mb_type | D_INTRA_MB;

            /* Macroblock Layer Begins */
            /* Decode the u1_mb_type */
            u1_mb_type = ih264d_parse_mb_type_intra_cabac(0, ps_dec);
            if(u1_mb_type > 25) return ERROR_MB_TYPE;
            ps_cur_mb_info->u1_mb_type = u1_mb_type;
            COPYTHECONTEXT("u1_mb_type", u1_mb_type);

            /* Parse Macroblock Data */
            if(25 == u1_mb_type)
            {
                /* I_PCM_MB */
                ps_cur_mb_info->ps_curmb->u1_mb_type = I_PCM_MB;
                ret = ih264d_parse_ipcm_mb(ps_dec, ps_cur_mb_info, u1_num_mbs);
                if(ret != OK) return ret;
                ps_cur_deblk_mb->u1_mb_qp = 0;
            }
            else
            {
                ret = ih264d_parse_imb_cabac(ps_dec, ps_cur_mb_info, u1_mb_type);
                if(ret != OK) return ret;
                ps_cur_deblk_mb->u1_mb_qp = ps_dec->u1_qp;
            }

            if(ps_dec->u1_enable_mb_info)
            {
                ih264d_populate_mb_info_map(ps_dec, ps_cur_mb_info, ps_cur_mb_info->u2_mbx << 1,
                                            ps_cur_mb_info->u2_mby << 1, ps_cur_deblk_mb->u1_mb_qp);
            }
            if(u1_mbaff)
            {
                ih264d_update_mbaff_left_nnz(ps_dec, ps_cur_mb_info);
            }

            if(ps_cur_mb_info->u1_topmb && u1_mbaff)
                uc_more_data_flag = 1;
            else
            {
                uc_more_data_flag = ih264d_decode_terminate(&ps_dec->s_cab_dec_env, ps_bitstrm);
                uc_more_data_flag = !uc_more_data_flag;
                COPYTHECONTEXT("Decode Sliceterm", !uc_more_data_flag);
            }

            if(u1_mbaff)
            {
                if(!uc_more_data_flag && (0 == (i2_cur_mb_addr & 1)))
                {
                    return ERROR_EOB_FLUSHBITS_T;
                }
            }
            /* Next macroblock information */
            i2_cur_mb_addr++;
            /* Store the colocated information */
            {
                mv_pred_t *ps_mv_nmb_start = ps_dec->ps_mv_cur + (u1_num_mbs << 4);
                mv_pred_t s_mvPred = {{0, 0, 0, 0}, {-1, -1}, 0, 0};
                ih264d_rep_mv_colz(ps_dec, &s_mvPred, ps_mv_nmb_start, 0,
                                   (UWORD8) (ps_dec->u1_cur_mb_fld_dec_flag << 1), 4, 4);
            }
            /*if num _cores is set to 3,compute bs will be done in another thread*/
            if(ps_dec->u4_num_cores < 3)
            {
                if(ps_dec->u4_app_disable_deblk_frm == 0)
                    ps_svc_lyr_dec->pf_svc_compute_bs(ps_svc_lyr_dec, ps_cur_mb_info,
                                                      (UWORD16) (u1_num_mbs >> u1_mbaff));
            }
            u1_num_mbs++;
        }

        /****************************************************************/
        /* Check for End Of Row                                         */
        /****************************************************************/
        u1_num_mbs_next = i2_pic_wdin_mbs - u2_mbx - 1;
        u1_end_of_row = (!u1_num_mbs_next) && (!(u1_mbaff && (u1_num_mbs & 0x01)));
        u1_tfr_n_mb =
            (u1_num_mbs == ps_dec->u1_recon_mb_grp) || u1_end_of_row || (!uc_more_data_flag);
        ps_cur_mb_info->u1_end_of_slice = (!uc_more_data_flag);

        if(u1_tfr_n_mb || (!uc_more_data_flag))
        {
            if(ps_dec->u1_separate_parse)
            {
                ih264d_parse_tfr_nmb(ps_dec, u1_mb_idx, u1_num_mbs, u1_num_mbs_next, u1_tfr_n_mb,
                                     u1_end_of_row);
                ps_dec->ps_nmb_info += u1_num_mbs;
                ps_svc_lyr_dec->ps_svc_nmb_info += u1_num_mbs;
            }
            else
            {
                if(ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER)
                {
                    ih264d_decode_recon_tfr_nmb(ps_dec, u1_mb_idx, u1_num_mbs, u1_num_mbs_next,
                                                u1_tfr_n_mb, u1_end_of_row);
                }
                else
                {
                    isvcd_decode_recon_tfr_nmb_base_lyr(ps_svc_lyr_dec, u1_mb_idx, u1_num_mbs,
                                                        u1_num_mbs_next, u1_tfr_n_mb,
                                                        u1_end_of_row);
                }
            }
            ps_dec->u2_total_mbs_coded += u1_num_mbs;
            if(u1_tfr_n_mb) u1_num_mbs = 0;
            u1_mb_idx = u1_num_mbs;
            ps_dec->u1_mb_idx = u1_num_mbs;
        }
    } while(uc_more_data_flag);

    ps_dec->u4_num_mbs_cur_nmb = 0;
    ps_dec->ps_cur_slice->u4_mbs_in_slice = i2_cur_mb_addr - (u2_first_mb_in_slice << u1_mbaff);

    return ret;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_islice_data_cavlc                            */
/*                                                                           */
/*  Description   : This function parses cabac syntax of a inter slice on    */
/*                  N MB basis.                                              */
/*                                                                           */
/*  Inputs        : ps_dec                                                   */
/*                  sliceparams                                              */
/*                  firstMbInSlice                                           */
/*                                                                           */
/*  Processing    : 1. After parsing syntax for N MBs those N MBs are        */
/*                     decoded till the end of slice.                        */
/*                                                                           */
/*  Returns       : 0                                                        */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         24 06 2005   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/

WORD32 isvcd_parse_islice_data_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                     dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice)
{
    UWORD8 uc_more_data_flag;
    UWORD8 u1_num_mbs, u1_mb_idx;
    dec_mb_info_t *ps_cur_mb_info;
    deblk_mb_t *ps_cur_deblk_mb;
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_bit_stream_t *const ps_bitstrm = ps_dec->ps_bitstrm;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD16 i2_pic_wdin_mbs = ps_dec->u2_frm_wd_in_mbs;
    WORD16 i2_cur_mb_addr;
    UWORD8 u1_mbaff;
    UWORD8 u1_num_mbs_next, u1_end_of_row, u1_tfr_n_mb;
    WORD32 ret = OK;

    ps_dec->u1_qp = ps_slice->u1_slice_qp;
    ih264d_update_qp(ps_dec, 0);
    u1_mbaff = ps_slice->u1_mbaff_frame_flag;

    /* initializations */
    u1_mb_idx = ps_dec->u1_mb_idx;
    u1_num_mbs = u1_mb_idx;

    uc_more_data_flag = 1;
    i2_cur_mb_addr = u2_first_mb_in_slice << u1_mbaff;
    do
    {
        UWORD8 u1_mb_type;
        ps_dec->pv_prev_mb_parse_tu_coeff_data = ps_dec->pv_parse_tu_coeff_data;
        if(i2_cur_mb_addr > ps_dec->ps_cur_sps->u2_max_mb_addr)
        {
            break;
        }

        ps_cur_mb_info = ps_dec->ps_nmb_info + u1_num_mbs;
        ps_dec->u4_num_mbs_cur_nmb = u1_num_mbs;
        ps_dec->u4_num_pmbair = (u1_num_mbs >> u1_mbaff);
        ps_cur_mb_info->u1_end_of_slice = 0;

        /***************************************************************/
        /* Get the required information for decoding of MB             */
        /* mb_x, mb_y , neighbour availablity,                         */
        /***************************************************************/
        ps_dec->pf_get_mb_info(ps_dec, i2_cur_mb_addr, ps_cur_mb_info, 0);

        /***************************************************************/
        /* Set the deblocking parameters for this MB                   */
        /***************************************************************/
        ps_cur_deblk_mb = ps_dec->ps_deblk_mbn + u1_num_mbs;
        if(ps_dec->u4_app_disable_deblk_frm == 0)
            ih264d_set_deblocking_parameters(ps_cur_deblk_mb, ps_slice,
                                             ps_dec->u1_mb_ngbr_availablity,
                                             ps_dec->u1_cur_mb_fld_dec_flag);

        ps_cur_deblk_mb->u1_mb_type = ps_cur_deblk_mb->u1_mb_type | D_INTRA_MB;

        /**************************************************************/
        /* Macroblock Layer Begins, Decode the u1_mb_type             */
        /**************************************************************/
        /* Inlined ih264d_uev */
        {
            UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
            UWORD32 u4_word, u4_ldz, u4_temp;

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
            *pu4_bitstrm_ofst = u4_bitstream_offset;
            u4_temp = ((1 << u4_ldz) + u4_word - 1);
            if(u4_temp > 25) return ERROR_MB_TYPE;
            u1_mb_type = u4_temp;
        }
        /* Inlined ih264d_uev */
        ps_cur_mb_info->u1_mb_type = u1_mb_type;
        COPYTHECONTEXT("u1_mb_type", u1_mb_type);

        /**************************************************************/
        /* Parse Macroblock data                                      */
        /**************************************************************/
        if(25 == u1_mb_type)
        {
            /* I_PCM_MB */
            ps_cur_mb_info->ps_curmb->u1_mb_type = I_PCM_MB;
            ret = ih264d_parse_ipcm_mb(ps_dec, ps_cur_mb_info, u1_num_mbs);
            if(ret != OK) return ret;
            ps_cur_deblk_mb->u1_mb_qp = 0;
        }
        else
        {
            ret = ih264d_parse_imb_cavlc(ps_dec, ps_cur_mb_info, u1_num_mbs, u1_mb_type);
            if(ret != OK) return ret;
            ps_cur_deblk_mb->u1_mb_qp = ps_dec->u1_qp;
        }
        if(ps_dec->u1_enable_mb_info)
        {
            ih264d_populate_mb_info_map(ps_dec, ps_cur_mb_info, ps_cur_mb_info->u2_mbx << 1,
                                        ps_cur_mb_info->u2_mby << 1, ps_cur_deblk_mb->u1_mb_qp);
        }
        uc_more_data_flag = MORE_RBSP_DATA(ps_bitstrm);
        if(u1_mbaff)
        {
            ih264d_update_mbaff_left_nnz(ps_dec, ps_cur_mb_info);
            if(!uc_more_data_flag && (0 == (i2_cur_mb_addr & 1)))
            {
                return ERROR_EOB_FLUSHBITS_T;
            }
        }
        /**************************************************************/
        /* Get next Macroblock address                                */
        /**************************************************************/
        i2_cur_mb_addr++;
        /* Store the colocated information */
        {
            mv_pred_t *ps_mv_nmb_start = ps_dec->ps_mv_cur + (u1_num_mbs << 4);
            mv_pred_t s_mvPred = {{0, 0, 0, 0}, {-1, -1}, 0, 0};
            ih264d_rep_mv_colz(ps_dec, &s_mvPred, ps_mv_nmb_start, 0,
                               (UWORD8) (ps_dec->u1_cur_mb_fld_dec_flag << 1), 4, 4);
        }

        /*if num _cores is set to 3,compute bs will be done in another thread*/
        if(ps_dec->u4_num_cores < 3)
        {
            if(ps_dec->u4_app_disable_deblk_frm == 0)
                ps_svc_lyr_dec->pf_svc_compute_bs(ps_svc_lyr_dec, ps_cur_mb_info,
                                                  (UWORD16) (u1_num_mbs >> u1_mbaff));
        }
        u1_num_mbs++;

        /****************************************************************/
        /* Check for End Of Row                                         */
        /****************************************************************/
        u1_num_mbs_next = i2_pic_wdin_mbs - ps_dec->u2_mbx - 1;
        u1_end_of_row = (!u1_num_mbs_next) && (!(u1_mbaff && (u1_num_mbs & 0x01)));
        u1_tfr_n_mb =
            (u1_num_mbs == ps_dec->u1_recon_mb_grp) || u1_end_of_row || (!uc_more_data_flag);
        ps_cur_mb_info->u1_end_of_slice = (!uc_more_data_flag);

        if(u1_tfr_n_mb || (!uc_more_data_flag))
        {
            if(ps_dec->u1_separate_parse)
            {
                ih264d_parse_tfr_nmb(ps_dec, u1_mb_idx, u1_num_mbs, u1_num_mbs_next, u1_tfr_n_mb,
                                     u1_end_of_row);
                ps_dec->ps_nmb_info += u1_num_mbs;
                ps_svc_lyr_dec->ps_svc_nmb_info += u1_num_mbs;
            }
            else
            {
                if(ps_svc_lyr_dec->u1_layer_identifier == TARGET_LAYER)
                {
                    ih264d_decode_recon_tfr_nmb(ps_dec, u1_mb_idx, u1_num_mbs, u1_num_mbs_next,
                                                u1_tfr_n_mb, u1_end_of_row);
                }
                else
                {
                    isvcd_decode_recon_tfr_nmb_base_lyr(ps_svc_lyr_dec, u1_mb_idx, u1_num_mbs,
                                                        u1_num_mbs_next, u1_tfr_n_mb,
                                                        u1_end_of_row);
                }
            }
            ps_dec->u2_total_mbs_coded += u1_num_mbs;
            if(u1_tfr_n_mb) u1_num_mbs = 0;
            u1_mb_idx = u1_num_mbs;
            ps_dec->u1_mb_idx = u1_num_mbs;
        }
    } while(uc_more_data_flag);

    ps_dec->u4_num_mbs_cur_nmb = 0;
    ps_dec->ps_cur_slice->u4_mbs_in_slice = i2_cur_mb_addr - (u2_first_mb_in_slice << u1_mbaff);

    return ret;
}
/*!
**************************************************************************
* \if Function name : ParseIMb \endif
*
* \brief
*    This function parses CAVLC syntax of a I MB. If 16x16 Luma DC transform
*    is also done here. Transformed Luma DC values are copied in their
*    0th pixel location of corrosponding CoeffBlock.
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_imb_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_num,
                             UWORD8 u1_mb_type)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    WORD32 i4_delta_qp;
    UWORD32 u4_temp;
    UWORD32 ui_is_top_mb_available;
    UWORD32 ui_is_left_mb_available;
    UWORD32 u4_cbp;
    UWORD32 u4_offset;
    UWORD32 *pu4_bitstrm_buf;
    WORD32 ret;
    dec_bit_stream_t *const ps_bitstrm = ps_dec->ps_bitstrm;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UNUSED(u1_mb_num);

    ps_cur_mb_info->u1_tran_form8x8 = 0;
    ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;
    ps_cur_mb_info->u1_yuv_dc_block_flag = 0;
    u4_temp = ps_dec->u1_mb_ngbr_availablity;
    ui_is_top_mb_available = BOOLEAN(u4_temp & TOP_MB_AVAILABLE_MASK);
    ui_is_left_mb_available = BOOLEAN(u4_temp & LEFT_MB_AVAILABLE_MASK);
    pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;

    ps_cur_mb_info->ps_curmb->u1_mb_type = P_MB;
    if(!ps_svc_cur_mb_info->u1_base_mode_flag)
    {
        if(u1_mb_type == I_4x4_MB)
        {
            ps_cur_mb_info->ps_curmb->u1_mb_type = I_4x4_MB;
            u4_offset = 0;

            /*--------------------------------------------------------------------*/
            /* Read transform_size_8x8_flag if present                            */
            /*--------------------------------------------------------------------*/
            if(ps_dec->ps_cur_pps->i4_transform_8x8_mode_flag)
            {
                ps_cur_mb_info->u1_tran_form8x8 = ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("transform_size_8x8_flag", ps_cur_mb_info->u1_tran_form8x8);
                ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = ps_cur_mb_info->u1_tran_form8x8;
            }

            /*--------------------------------------------------------------------*/
            /* Read the IntraPrediction modes for LUMA                            */
            /*--------------------------------------------------------------------*/
            if(!ps_cur_mb_info->u1_tran_form8x8)
            {
                UWORD8 *pu1_temp;
                ih264d_read_intra_pred_modes(ps_dec, ((UWORD8 *) ps_dec->pv_parse_tu_coeff_data),
                                             ((UWORD8 *) ps_dec->pv_parse_tu_coeff_data + 16),
                                             ps_cur_mb_info->u1_tran_form8x8);
                pu1_temp = (UWORD8 *) ps_dec->pv_parse_tu_coeff_data;
                pu1_temp += 32;
                ps_dec->pv_parse_tu_coeff_data = (void *) pu1_temp;
            }
            else
            {
                UWORD8 *pu1_temp;
                ih264d_read_intra_pred_modes(ps_dec, ((UWORD8 *) ps_dec->pv_parse_tu_coeff_data),
                                             ((UWORD8 *) ps_dec->pv_parse_tu_coeff_data + 4),
                                             ps_cur_mb_info->u1_tran_form8x8);
                pu1_temp = (UWORD8 *) ps_dec->pv_parse_tu_coeff_data;
                pu1_temp += 8;
                ps_dec->pv_parse_tu_coeff_data = (void *) pu1_temp;
            }
            /*--------------------------------------------------------------------*/
            /* Read the IntraPrediction mode for CHROMA                           */
            /*--------------------------------------------------------------------*/
            /* Inlined ih264d_uev */
            {
                UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
                UWORD32 u4_word, u4_ldz, u4_temp;

                /***************************************************************/
                /* Find leading zeros in next 32 bits                          */
                /***************************************************************/
                NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
                u4_ldz = CLZ(u4_word);
                /* Flush the ps_bitstrm */
                u4_bitstream_offset += (u4_ldz + 1);
                /* Read the suffix from the ps_bitstrm */
                u4_word = 0;
                if(u4_ldz)
                {
                    GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
                }
                *pu4_bitstrm_ofst = u4_bitstream_offset;
                u4_temp = ((1 << u4_ldz) + u4_word - 1);
                if(u4_temp > 3)
                {
                    return ERROR_CHROMA_PRED_MODE;
                }
                ps_cur_mb_info->u1_chroma_pred_mode = u4_temp;
                COPYTHECONTEXT("intra_chroma_pred_mode", ps_cur_mb_info->u1_chroma_pred_mode);
            }
            /*--------------------------------------------------------------------*/
            /* Read the Coded block pattern                                       */
            /*--------------------------------------------------------------------*/
            {
                UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
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
                if(u4_ldz)
                {
                    GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
                }
                *pu4_bitstrm_ofst = u4_bitstream_offset;
                u4_cbp = ((1 << u4_ldz) + u4_word - 1);
            }
            if(u4_cbp > 47)
            {
                return ERROR_CBP;
            }

            u4_cbp = gau1_ih264d_cbp_table[u4_cbp][0];
            COPYTHECONTEXT("coded_block_pattern", u4_cbp);
            ps_cur_mb_info->u1_cbp = u4_cbp;

            /*--------------------------------------------------------------------*/
            /* Read mb_qp_delta                                                   */
            /*--------------------------------------------------------------------*/
            if(ps_cur_mb_info->u1_cbp)
            {
                UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
                UWORD32 u4_word, u4_ldz, u4_abs_val;

                /***************************************************************/
                /* Find leading zeros in next 32 bits                          */
                /***************************************************************/
                NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
                u4_ldz = CLZ(u4_word);

                /* Flush the ps_bitstrm */
                u4_bitstream_offset += (u4_ldz + 1);

                /* Read the suffix from the ps_bitstrm */
                u4_word = 0;
                if(u4_ldz)
                {
                    GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
                }

                *pu4_bitstrm_ofst = u4_bitstream_offset;
                u4_abs_val = ((1 << u4_ldz) + u4_word) >> 1;

                if(u4_word & 0x1)
                {
                    i4_delta_qp = (-(WORD32) u4_abs_val);
                }
                else
                {
                    i4_delta_qp = (u4_abs_val);
                }

                if((i4_delta_qp < -26) || (i4_delta_qp > 25))
                {
                    return ERROR_INV_RANGE_QP_T;
                }

                COPYTHECONTEXT("mb_qp_delta", i1_delta_qp);
                if(i4_delta_qp != 0)
                {
                    ret = ih264d_update_qp(ps_dec, (WORD8) i4_delta_qp);
                    if(ret != OK) return ret;
                }
            }
        }
        else
        {
            u4_offset = 1;
            ps_cur_mb_info->ps_curmb->u1_mb_type = I_16x16_MB;
            /*-------------------------------------------------------------------*/
            /* Read the IntraPrediction mode for CHROMA                          */
            /*-------------------------------------------------------------------*/
            {
                UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
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
                if(u4_ldz)
                {
                    GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
                }
                *pu4_bitstrm_ofst = u4_bitstream_offset;
                u4_temp = ((1 << u4_ldz) + u4_word - 1);

                /* Inlined ih264d_uev */
                if(u4_temp > 3)
                {
                    return ERROR_CHROMA_PRED_MODE;
                }
                ps_cur_mb_info->u1_chroma_pred_mode = u4_temp;
                COPYTHECONTEXT("intra_chroma_pred_mode", ps_cur_mb_info->u1_chroma_pred_mode);
            }
            /*-------------------------------------------------------------------*/
            /* Read the Coded block pattern                                      */
            /*-------------------------------------------------------------------*/
            u4_cbp = gau1_ih264d_cbp_tab[(u1_mb_type - 1) >> 2];
            ps_cur_mb_info->u1_cbp = u4_cbp;

            /*-------------------------------------------------------------------*/
            /* Read mb_qp_delta                                                  */
            /*-------------------------------------------------------------------*/
            {
                UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
                UWORD32 u4_word, u4_ldz, u4_abs_val;

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

                *pu4_bitstrm_ofst = u4_bitstream_offset;
                u4_abs_val = ((1 << u4_ldz) + u4_word) >> 1;

                if(u4_word & 0x1)
                    i4_delta_qp = (-(WORD32) u4_abs_val);
                else
                    i4_delta_qp = (u4_abs_val);

                if((i4_delta_qp < -26) || (i4_delta_qp > 25)) return ERROR_INV_RANGE_QP_T;
            }
            /* inlinined ih264d_sev */
            COPYTHECONTEXT("Delta quant", i1_delta_qp);

            if(i4_delta_qp != 0)
            {
                ret = ih264d_update_qp(ps_dec, (WORD8) i4_delta_qp);
                if(ret != OK) return ret;
            }

            {
                WORD16 i_scaleFactor;
                UWORD32 ui_N = 0;
                WORD16 *pi2_scale_matrix_ptr;
                /*******************************************************************/
                /* for luma DC coefficients the scaling is done during the parsing */
                /* to preserve the precision                                       */
                /*******************************************************************/
                if(ps_dec->s_high_profile.u1_scaling_present)
                {
                    pi2_scale_matrix_ptr = ps_dec->s_high_profile.i2_scalinglist4x4[0];
                }
                else
                {
                    i_scaleFactor = 16;
                    pi2_scale_matrix_ptr = &i_scaleFactor;
                }

                /*---------------------------------------------------------------*/
                /* Decode DC coefficients                                        */
                /*---------------------------------------------------------------*/
                /*---------------------------------------------------------------*/
                /* Calculation of N                                              */
                /*---------------------------------------------------------------*/
                if(ui_is_left_mb_available)
                {
                    if(ui_is_top_mb_available)
                    {
                        ui_N = ((ps_cur_mb_info->ps_top_mb->pu1_nnz_y[0] +
                                 ps_dec->pu1_left_nnz_y[0] + 1) >>
                                1);
                    }
                    else
                    {
                        ui_N = ps_dec->pu1_left_nnz_y[0];
                    }
                }
                else if(ui_is_top_mb_available)
                {
                    ui_N = ps_cur_mb_info->ps_top_mb->pu1_nnz_y[0];
                }

                {
                    WORD16 pi2_dc_coef[16] = {0};
                    WORD32 pi4_tmp[16] = {0};
                    tu_sblk4x4_coeff_data_t *ps_tu_4x4 =
                        (tu_sblk4x4_coeff_data_t *) ps_dec->pv_parse_tu_coeff_data;
                    WORD16 *pi2_coeff_block = (WORD16 *) ps_dec->pv_parse_tu_coeff_data;
                    UWORD32 u4_num_coeff;
                    ps_tu_4x4->u2_sig_coeff_map = 0;
                    ret = ps_dec->pf_cavlc_parse4x4coeff[(ui_N > 7)](pi2_dc_coef, 0, ui_N, ps_dec,
                                                                     &u4_num_coeff);
                    if(ret != OK) return ret;

                    if(EXCEED_OFFSET(ps_bitstrm)) return ERROR_EOB_TERMINATE_T;
                    if(ps_tu_4x4->u2_sig_coeff_map)
                    {
                        memset(pi2_dc_coef, 0, sizeof(pi2_dc_coef));
                        ih264d_unpack_coeff4x4_dc_4x4blk(ps_tu_4x4, pi2_dc_coef,
                                                         ps_dec->pu1_inv_scan);

                        PROFILE_DISABLE_IQ_IT_RECON()
                        ps_dec->pf_ihadamard_scaling_4x4(
                            pi2_dc_coef, pi2_coeff_block, ps_dec->pu2_quant_scale_y,
                            (UWORD16 *) pi2_scale_matrix_ptr, ps_dec->u1_qp_y_div6, pi4_tmp);
                        pi2_coeff_block += 16;
                        ps_dec->pv_parse_tu_coeff_data = (void *) pi2_coeff_block;
                        SET_BIT(ps_cur_mb_info->u1_yuv_dc_block_flag, 0);
                    }
                }
            }
        }
    }
    else
    {
        ps_cur_mb_info->u1_tran_form8x8 = 0;
        ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;
        u4_offset = 0;
        /*--------------------------------------------------------------------*/
        /* Read the Coded block pattern                                       */
        /*--------------------------------------------------------------------*/
        {
            UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
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
            if(u4_ldz)
            {
                GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
            }
            *pu4_bitstrm_ofst = u4_bitstream_offset;
            u4_cbp = ((1 << u4_ldz) + u4_word - 1);
        }
        if(u4_cbp > 47)
        {
            return ERROR_CBP;
        }

        /*  inter cbp table to be used for base mode flag*/
        u4_cbp = gau1_ih264d_cbp_table[u4_cbp][1];
        COPYTHECONTEXT("coded_block_pattern", u4_cbp);
        ps_cur_mb_info->u1_cbp = u4_cbp;
        ps_cur_mb_info->u1_tran_form8x8 = 0;
        ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;

        /*--------------------------------------------------------------------*/
        /* Read transform_size_8x8_flag if present                            */
        /*--------------------------------------------------------------------*/
        if((ps_dec->ps_cur_pps->i4_transform_8x8_mode_flag) && (ps_cur_mb_info->u1_cbp & 0xf))
        {
            ps_cur_mb_info->u1_tran_form8x8 = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("transform_size_8x8_flag", ps_cur_mb_info->u1_tran_form8x8);
            ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = ps_cur_mb_info->u1_tran_form8x8;
        }
        /*--------------------------------------------------------------------*/
        /* Read mb_qp_delta                                                   */
        /*--------------------------------------------------------------------*/
        if(ps_cur_mb_info->u1_cbp)
        {
            UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
            UWORD32 u4_word, u4_ldz, u4_abs_val;

            /***************************************************************/
            /* Find leading zeros in next 32 bits                          */
            /***************************************************************/
            NEXTBITS_32(u4_word, u4_bitstream_offset, pu4_bitstrm_buf);
            u4_ldz = CLZ(u4_word);

            /* Flush the ps_bitstrm */
            u4_bitstream_offset += (u4_ldz + 1);

            /* Read the suffix from the ps_bitstrm */
            u4_word = 0;
            if(u4_ldz)
            {
                GETBITS(u4_word, u4_bitstream_offset, pu4_bitstrm_buf, u4_ldz);
            }

            *pu4_bitstrm_ofst = u4_bitstream_offset;
            u4_abs_val = ((1 << u4_ldz) + u4_word) >> 1;

            if(u4_word & 0x1)
            {
                i4_delta_qp = (-(WORD32) u4_abs_val);
            }
            else
            {
                i4_delta_qp = (u4_abs_val);
            }

            if((i4_delta_qp < -26) || (i4_delta_qp > 25))
            {
                return ERROR_INV_RANGE_QP_T;
            }

            COPYTHECONTEXT("mb_qp_delta", i4_delta_qp);
            if(i4_delta_qp != 0)
            {
                ret = ih264d_update_qp(ps_dec, (WORD8) i4_delta_qp);
                if(ret != OK) return ret;
            }
        }
    }
    if(u4_cbp)
    {
        ret = ih264d_parse_residual4x4_cavlc(ps_dec, ps_cur_mb_info, (UWORD8) u4_offset);
        if(ret != OK) return ret;
        if(EXCEED_OFFSET(ps_bitstrm)) return ERROR_EOB_TERMINATE_T;

        /* Store Left Mb NNZ and TOP chroma NNZ */
    }
    else
    {
        ps_cur_mb_info->u1_qp_div6 = ps_dec->u1_qp_y_div6;
        ps_cur_mb_info->u1_qpc_div6 = ps_dec->u1_qp_u_div6;
        ps_cur_mb_info->u1_qpcr_div6 = ps_dec->u1_qp_v_div6;
        ps_cur_mb_info->u1_qp_rem6 = ps_dec->u1_qp_y_rem6;
        ps_cur_mb_info->u1_qpc_rem6 = ps_dec->u1_qp_u_rem6;
        ps_cur_mb_info->u1_qpcr_rem6 = ps_dec->u1_qp_v_rem6;
        ih264d_update_nnz_for_skipmb(ps_dec, ps_cur_mb_info, CAVLC);
    }

    return OK;
}

/*!
**************************************************************************
* \if Function name : ih264d_decode_eislice \endif
*
* \brief
*    Decodes an EI Slice
*
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_eislice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    WORD32 i_status = OK;
    dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
    dec_slice_params_t *ps_slice = ps_dec->ps_cur_slice;
    dec_bit_stream_t *ps_bitstrm = ps_dec->ps_bitstrm;
    dec_seq_params_t *ps_seq;
    dec_svc_seq_params_t *ps_subset_seq;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;
    dec_subset_seq_params_t *ps_sps_svc_ext = NULL;
    dec_nal_unit_svc_ext_params_t *ps_nal_svc_ext = NULL;
    UWORD32 u4_temp;
    WORD32 i_temp;
    WORD32 ret;
    UWORD32 *pu4_bitstrm_buf = ps_dec->ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_dec->ps_bitstrm->u4_ofst;

    ps_nal_svc_ext = ps_svc_lyr_dec->ps_nal_svc_ext;
    ps_seq = ps_dec->ps_cur_sps;
    ps_subset_seq =
        &ps_svc_lyr_dec->ps_subset_sps[MAX_NUM_SEQ_PARAMS + ps_seq->u1_seq_parameter_set_id];
    ps_sps_svc_ext = &ps_subset_seq->s_sps_svc_ext;
    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;

    if(0 == ps_svc_lyr_dec->ps_nal_svc_ext->u1_quality_id)
    {
        /*--------------------------------------------------------------------*/
        /* Read remaining contents of the slice header                        */
        /*--------------------------------------------------------------------*/
        /* dec_ref_pic_marking function */
        /* G050 */
        if(ps_slice->u1_nal_ref_idc != 0)
        {
            if(!ps_dec->ps_dpb_cmds->u1_dpb_commands_read)
            {
                dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
                dec_seq_params_t *ps_sps_tmp = ps_pps->ps_sps;
                UWORD8 u1_nal_unit_type_tmp = ps_dec->u1_nal_unit_type;

                ps_pps->ps_sps = ps_dec->ps_cur_sps;
                if(ps_svc_lyr_dec->ps_nal_svc_ext->u1_idr_flag)
                    ps_dec->u1_nal_unit_type = IDR_SLICE_NAL;

                i_temp = ih264d_read_mmco_commands(ps_dec);

                ps_pps->ps_sps = ps_sps_tmp;
                ps_dec->u1_nal_unit_type = u1_nal_unit_type_tmp;
                if(i_temp < 0)
                {
                    return ERROR_DBP_MANAGER_T;
                }
                ps_dec->u4_bitoffset = i_temp;
            }
            else
                ps_dec->ps_bitstrm->u4_ofst += ps_dec->u4_bitoffset;

            if(!ps_sps_svc_ext->u1_slice_header_restriction_flag)
            {
                ps_svc_slice_params->u1_store_ref_base_pic_flag = ih264d_get_bit_h264(ps_bitstrm);
                COPYTHECONTEXT("SPS_SVC_EXT: u1_store_ref_base_pic_flag",
                               ps_svc_slice_params->u1_store_ref_base_pic_flag);

                if(0 != ps_svc_slice_params->u1_store_ref_base_pic_flag)
                {
                    return NOT_OK;
                }
                if(((1 == ps_nal_svc_ext->u1_use_ref_base_pic_flag) ||
                    (1 == ps_svc_slice_params->u1_store_ref_base_pic_flag)) &&
                   (!ps_nal_svc_ext->u1_idr_flag))
                {
                    i_status = isvcd_dec_ref_base_pic_marking(
                        &ps_svc_slice_params->s_ref_base_pic_marking_svc_ext, ps_bitstrm);
                    if(i_status != OK)
                    {
                        return i_status;
                    }
                }
            }
        }
    }

    {
        /* G050 */
        /* Read slice_qp_delta */
        WORD64 i8_temp =
            (WORD64) ps_pps->u1_pic_init_qp + ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if((i8_temp < MIN_H264_QP) || (i8_temp > MAX_H264_QP)) return ERROR_INV_RANGE_QP_T;
        ps_slice->u1_slice_qp = (UWORD8) i8_temp;
        COPYTHECONTEXT("Slice Header SVC ext: slice_qp_delta",
                       ps_slice->u1_slice_qp - ps_pps->u1_pic_init_qp);
    }

    if(ps_pps->u1_deblocking_filter_parameters_present_flag == 1)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("Slice Header SVC ext: disable_deblocking_filter_idc", u4_temp);

        if(u4_temp > SLICE_BOUNDARY_DBLK_DISABLED)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        ps_slice->u1_disable_dblk_filter_idc = u4_temp;
        if(u4_temp != 1)
        {
            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf) << 1;
            if((MIN_DBLK_FIL_OFF > i_temp) || (i_temp > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice->i1_slice_alpha_c0_offset = i_temp;
            COPYTHECONTEXT("Slice Header SVC ext: slice_alpha_c0_offset_div2",
                           ps_slice->i1_slice_alpha_c0_offset >> 1);

            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf) << 1;
            if((MIN_DBLK_FIL_OFF > i_temp) || (i_temp > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice->i1_slice_beta_offset = i_temp;
            COPYTHECONTEXT("Slice Header SVC ext: slice_beta_offset_div2",
                           ps_slice->i1_slice_beta_offset >> 1);
        }
        else
        {
            ps_slice->i1_slice_alpha_c0_offset = 0;
            ps_slice->i1_slice_beta_offset = 0;
        }
    }
    else
    {
        ps_slice->u1_disable_dblk_filter_idc = 0;
        ps_slice->i1_slice_alpha_c0_offset = 0;
        ps_slice->i1_slice_beta_offset = 0;
    }

    /* add the remaining part of the code for svc extension from reference */
    ret = isvcd_set_default_slice_header_ext(ps_svc_lyr_dec);
    if(ret != OK)
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    ret = isvcd_parse_slice_header(ps_svc_lyr_dec);
    if(ret != OK)
    {
        return ERROR_INV_SLICE_HDR_T;
    }

    /* Initialization to check if number of motion vector per 2 Mbs */
    /* are exceeding the range or not */
    ps_dec->u2_mv_2mb[0] = 0;
    ps_dec->u2_mv_2mb[1] = 0;

    ret = isvcd_parse_interlayer_resamp_func_init(ps_svc_lyr_dec, u2_first_mb_in_slice);
    if(ret != OK)
    {
        return ERROR_CORRUPTED_SLICE;
    }
    /*set slice header cone to 2 ,to indicate  correct header*/
    ps_dec->u1_slice_header_done = 2;

    if(!ps_svc_slice_params->u1_slice_skip_flag)
    {
        if(ps_pps->u1_entropy_coding_mode)
        {
            SWITCHOFFTRACE;
            SWITCHONTRACECABAC;
            if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
            {
                ps_dec->pf_get_mb_info = ih264d_get_mb_info_cabac_mbaff;
            }
            else
                ps_dec->pf_get_mb_info = isvcd_get_mb_info_cabac_nonmbaff;

            ret = isvcd_parse_eislice_data_cabac(ps_svc_lyr_dec, ps_slice, u2_first_mb_in_slice);
            if(ret != OK) return ret;
            SWITCHONTRACE;
            SWITCHOFFTRACECABAC;
        }
        else
        {
            if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
            {
                ps_dec->pf_get_mb_info = ih264d_get_mb_info_cavlc_mbaff;
            }
            else
                ps_dec->pf_get_mb_info = isvcd_get_mb_info_cavlc_nonmbaff;
            ret = isvcd_parse_eislice_data_cavlc(ps_svc_lyr_dec, ps_slice, u2_first_mb_in_slice);
            if(ret != OK) return ret;
        }
    }
    else
    {
        return ERROR_FEATURE_UNAVAIL;
    }

    return OK;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_eislice_data_cabac                          */
/*                                                                           */
/*  Description   : This function parses cabac syntax of a inter slice on    */
/*                  N MB basis.                                              */
/*                                                                           */
/*  Inputs        : ps_dec                                                   */
/*                  sliceparams                                              */
/*                  firstMbInSlice                                           */
/*                                                                           */
/*  Processing    : 1. After parsing syntax for N MBs those N MBs are        */
/*                     decoded till the end of slice.                        */
/*                                                                           */
/*  Returns       : 0                                                        */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         24 06 2005   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_parse_eislice_data_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                      dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD8 uc_more_data_flag;
    UWORD8 u1_num_mbs, u1_mb_idx;
    dec_mb_info_t *ps_cur_mb_info;
    dec_svc_mb_info_t *ps_svc_cur_mb_info;
    deblk_mb_t *ps_cur_deblk_mb;
    dec_bit_stream_t *const ps_bitstrm = ps_dec->ps_bitstrm;
    UWORD16 i2_pic_wdin_mbs = ps_dec->u2_frm_wd_in_mbs;
    WORD16 i2_cur_mb_addr;
    UWORD8 u1_mbaff;
    UWORD8 u1_num_mbs_next, u1_end_of_row, u1_tfr_n_mb;
    WORD32 ret = OK;
    decoding_envirnoment_t *ps_cab_env;
    UWORD8 *pu1_cur_svc_base_mode_flag;
    UWORD8 u1_left_svc_base_mode_flag;
    UWORD8 u1_top_svc_base_mode_flag;
    UWORD32 u4_a, u4_b, u4_ctxt_inc;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;

    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;
    ps_dec->u1_qp = ps_slice->u1_slice_qp;
    ih264d_update_qp(ps_dec, 0);
    u1_mbaff = ps_slice->u1_mbaff_frame_flag;

    if(ps_bitstrm->u4_ofst & 0x07)
    {
        ps_bitstrm->u4_ofst += 8;
        ps_bitstrm->u4_ofst &= 0xFFFFFFF8;
    }
    ret = ih264d_init_cabac_dec_envirnoment(&(ps_dec->s_cab_dec_env), ps_bitstrm);
    if(ret != OK) return ret;
    isvcd_init_cabac_contexts(I_SLICE, ps_dec);
    ps_dec->i1_prev_mb_qp_delta = 0;
    ps_cab_env = &ps_dec->s_cab_dec_env;
    /* initializations */
    u1_mb_idx = ps_dec->u1_mb_idx;
    u1_num_mbs = u1_mb_idx;
    uc_more_data_flag = 1;
    i2_cur_mb_addr = u2_first_mb_in_slice << u1_mbaff;
    do
    {
        UWORD16 u2_mbx;
        ps_dec->pv_prev_mb_parse_tu_coeff_data = ps_dec->pv_parse_tu_coeff_data;
        if(i2_cur_mb_addr > ps_dec->ps_cur_sps->u2_max_mb_addr)
        {
            break;
        }

        {
            UWORD8 u1_mb_type;
            ps_cur_mb_info = ps_dec->ps_nmb_info + u1_num_mbs;
            ps_svc_cur_mb_info = ps_svc_lyr_dec->ps_svc_nmb_info + u1_num_mbs;
            ps_dec->u4_num_mbs_cur_nmb = u1_num_mbs;
            ps_dec->u4_num_pmbair = (u1_num_mbs >> u1_mbaff);
            ps_cur_mb_info->u1_end_of_slice = 0;
            /***************************************************************/
            /* Get the required information for decoding of MB                  */
            /* mb_x, mb_y , neighbour availablity,                              */
            /***************************************************************/
            ps_dec->pf_get_mb_info(ps_dec, i2_cur_mb_addr, ps_cur_mb_info, 0);
            u2_mbx = ps_dec->u2_mbx;
            ps_svc_cur_mb_info->u1_crop_window_flag =
                *(ps_svc_lyr_dec->pu1_crop_wnd_flag + ps_cur_mb_info->u2_mbx +
                  (ps_cur_mb_info->u2_mby * ps_dec->u2_frm_wd_in_mbs));
            /*********************************************************************/
            /* initialize u1_tran_form8x8 to zero to aviod uninitialized accesses */
            /*********************************************************************/
            ps_cur_mb_info->u1_tran_form8x8 = 0;
            ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;
            /***************************************************************/
            /* Set the deblocking parameters for this MB                   */
            /***************************************************************/
            ps_cur_deblk_mb = ps_dec->ps_deblk_mbn + u1_num_mbs;
            if(ps_dec->u4_app_disable_deblk_frm == 0)
                ih264d_set_deblocking_parameters(ps_cur_deblk_mb, ps_slice,
                                                 ps_dec->u1_mb_ngbr_availablity,
                                                 ps_dec->u1_cur_mb_fld_dec_flag);

            ps_svc_cur_mb_info->u1_base_mode_flag = 0;
            /* Macroblock Layer Begins */
            if(ps_svc_cur_mb_info->u1_crop_window_flag &&
               ps_svc_slice_params->u1_adaptive_base_mode_flag)
            {
                pu1_cur_svc_base_mode_flag =
                    ps_svc_lyr_dec->pu1_svc_base_mode_flag + ps_cur_mb_info->u2_mbx;
                pu1_cur_svc_base_mode_flag +=
                    ps_cur_mb_info->u2_mby * ps_svc_lyr_dec->i4_frm_svc_base_mode_cabac_stride;

                u1_left_svc_base_mode_flag = 0;
                if(ps_dec->u1_mb_ngbr_availablity & LEFT_MB_AVAILABLE_MASK)
                    u1_left_svc_base_mode_flag = *(pu1_cur_svc_base_mode_flag - 1);

                u1_top_svc_base_mode_flag = 0;
                if(ps_dec->u1_mb_ngbr_availablity & TOP_MB_AVAILABLE_MASK)
                    u1_top_svc_base_mode_flag =
                        *(pu1_cur_svc_base_mode_flag -
                          ps_svc_lyr_dec->i4_frm_svc_base_mode_cabac_stride);

                u4_a = 1;
                u4_b = 1;

                if(u1_top_svc_base_mode_flag)
                {
                    u4_a = 0;
                }

                if(u1_left_svc_base_mode_flag)
                {
                    u4_b = 0;
                }

                u4_ctxt_inc = u4_a + u4_b;
                ps_svc_cur_mb_info->u1_base_mode_flag = ih264d_decode_bin(
                    u4_ctxt_inc, ps_svc_lyr_dec->ps_base_mode_flag, ps_bitstrm, ps_cab_env);
                COPYTHECONTEXT("SVC ext: u1_base_mode_flag", ps_cur_mb_info->u1_base_mode_flag);
                *pu1_cur_svc_base_mode_flag = ps_svc_cur_mb_info->u1_base_mode_flag;
            }
            else if(ps_svc_cur_mb_info->u1_crop_window_flag)
            {
                ps_svc_cur_mb_info->u1_base_mode_flag =
                    ps_svc_slice_params->u1_default_base_mode_flag;
            }

            if(!ps_svc_cur_mb_info->u1_base_mode_flag)
            {
                /* Decode the u1_mb_type */
                u1_mb_type = ih264d_parse_mb_type_intra_cabac(0, ps_dec);
                if(u1_mb_type > 25) return ERROR_MB_TYPE;
                ps_cur_mb_info->u1_mb_type = u1_mb_type;
                COPYTHECONTEXT("u1_mb_type", u1_mb_type);
                ps_cur_deblk_mb->u1_mb_type = ps_cur_deblk_mb->u1_mb_type | D_INTRA_MB;
            }
            else
            {
                ps_cur_mb_info->u1_mb_type = MB_INFER;
                ps_cur_deblk_mb->u1_mb_type = ps_cur_deblk_mb->u1_mb_type | D_INTRA_IBL;
            }
            /* Parse Macroblock Data */
            u1_mb_type = ps_cur_mb_info->u1_mb_type;
            if(25 == u1_mb_type)
            {
                /* I_PCM_MB */
                ps_cur_mb_info->ps_curmb->u1_mb_type = I_PCM_MB;
                ret = ih264d_parse_ipcm_mb(ps_dec, ps_cur_mb_info, u1_num_mbs);
                if(ret != OK) return ret;
                ps_cur_deblk_mb->u1_mb_qp = 0;
            }
            else
            {
                ret = isvcd_parse_imb_cabac(ps_svc_lyr_dec, ps_cur_mb_info, ps_svc_cur_mb_info,
                                            u1_mb_type);
                if(ret != OK) return ret;
                ps_cur_deblk_mb->u1_mb_qp = ps_dec->u1_qp;
            }

            if(ps_dec->u1_enable_mb_info)
            {
                ih264d_populate_mb_info_map(ps_dec, ps_cur_mb_info, ps_cur_mb_info->u2_mbx << 1,
                                            ps_cur_mb_info->u2_mby << 1, ps_cur_deblk_mb->u1_mb_qp);
            }
            if(u1_mbaff)
            {
                ih264d_update_mbaff_left_nnz(ps_dec, ps_cur_mb_info);
            }

            if(ps_cur_mb_info->u1_topmb && u1_mbaff)
                uc_more_data_flag = 1;
            else
            {
                uc_more_data_flag = ih264d_decode_terminate(&ps_dec->s_cab_dec_env, ps_bitstrm);
                uc_more_data_flag = !uc_more_data_flag;
                COPYTHECONTEXT("Decode Sliceterm", !uc_more_data_flag);
            }

            if(u1_mbaff)
            {
                if(!uc_more_data_flag && (0 == (i2_cur_mb_addr & 1)))
                {
                    return ERROR_EOB_FLUSHBITS_T;
                }
            }
            /* Next macroblock information */
            i2_cur_mb_addr++;
            /* Store the colocated information */
            {
                mv_pred_t *ps_mv_nmb_start = ps_dec->ps_mv_cur + (u1_num_mbs << 4);
                mv_pred_t s_mvPred = {{0, 0, 0, 0}, {-1, -1}, 0, 0};
                if(ps_mv_nmb_start)
                {
                    ih264d_rep_mv_colz(ps_dec, &s_mvPred, ps_mv_nmb_start, 0,
                                       (UWORD8) (ps_dec->u1_cur_mb_fld_dec_flag << 1), 4, 4);
                }
                else
                {
                    return NOT_OK;
                }
            }

            u1_num_mbs++;
        }

        /****************************************************************/
        /* Check for End Of Row                                         */
        /****************************************************************/
        u1_num_mbs_next = i2_pic_wdin_mbs - u2_mbx - 1;
        u1_end_of_row = (!u1_num_mbs_next) && (!(u1_mbaff && (u1_num_mbs & 0x01)));
        u1_tfr_n_mb =
            (u1_num_mbs == ps_dec->u1_recon_mb_grp) || u1_end_of_row || (!uc_more_data_flag);
        ps_cur_mb_info->u1_end_of_slice = (!uc_more_data_flag);

        if(u1_tfr_n_mb || (!uc_more_data_flag))
        {
            if(ps_dec->u1_separate_parse)
            {
                ih264d_parse_tfr_nmb(ps_dec, u1_mb_idx, u1_num_mbs, u1_num_mbs_next, u1_tfr_n_mb,
                                     u1_end_of_row);
                ps_dec->ps_nmb_info += u1_num_mbs;
                ps_svc_lyr_dec->ps_svc_nmb_info += u1_num_mbs;
            }
            else
            {
                ret = isvcd_decode_recon_tfr_nmb_non_base_lyr(ps_svc_lyr_dec, u1_mb_idx, u1_num_mbs,
                                                              u1_num_mbs_next, u1_tfr_n_mb,
                                                              u1_end_of_row);
                if(ret != OK) return ret;
            }
            ps_dec->u2_total_mbs_coded += u1_num_mbs;
            if(u1_tfr_n_mb) u1_num_mbs = 0;
            u1_mb_idx = u1_num_mbs;
            ps_dec->u1_mb_idx = u1_num_mbs;
        }
    } while(uc_more_data_flag);

    ps_dec->u4_num_mbs_cur_nmb = 0;
    ps_dec->ps_cur_slice->u4_mbs_in_slice = i2_cur_mb_addr - (u2_first_mb_in_slice << u1_mbaff);

    return ret;
}

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_parse_eislice_data_cavlc */
/*                                                                           */
/*  Description   : This function parses cabac syntax of a inter slice on    */
/*                  N MB basis.                                              */
/*                                                                           */
/*  Inputs        : ps_dec                                                   */
/*                  sliceparams                                              */
/*                  firstMbInSlice                                           */
/*                                                                           */
/*  Processing    : 1. After parsing syntax for N MBs those N MBs are        */
/*                     decoded till the end of slice.                        */
/*                                                                           */
/*  Returns       : 0                                                        */
/*                                                                           */
/*  Issues        : <List any issues or problems with this function>         */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         24 06 2005   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/
WORD32 isvcd_parse_eislice_data_cavlc(svc_dec_lyr_struct_t *ps_svc_lyr_dec,
                                      dec_slice_params_t *ps_slice, UWORD16 u2_first_mb_in_slice)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    UWORD8 uc_more_data_flag;
    UWORD8 u1_num_mbs, u1_mb_idx;
    dec_mb_info_t *ps_cur_mb_info;
    dec_svc_mb_info_t *ps_svc_cur_mb_info;
    deblk_mb_t *ps_cur_deblk_mb;
    dec_bit_stream_t *const ps_bitstrm = ps_dec->ps_bitstrm;
    UWORD32 *pu4_bitstrm_ofst = &ps_bitstrm->u4_ofst;
    UWORD32 *pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
    UWORD16 i2_pic_wdin_mbs = ps_dec->u2_frm_wd_in_mbs;
    WORD16 i2_cur_mb_addr;
    UWORD8 u1_mbaff;
    UWORD8 u1_num_mbs_next, u1_end_of_row, u1_tfr_n_mb;
    WORD32 ret = OK;
    dec_slice_svc_ext_params_t *ps_svc_slice_params = NULL;

    ps_svc_slice_params = &ps_svc_lyr_dec->s_svc_slice_params;
    ps_dec->u1_qp = ps_slice->u1_slice_qp;
    ih264d_update_qp(ps_dec, 0);
    u1_mbaff = ps_slice->u1_mbaff_frame_flag;

    u1_mb_idx = ps_dec->u1_mb_idx;
    u1_num_mbs = u1_mb_idx;
    uc_more_data_flag = 1;
    i2_cur_mb_addr = u2_first_mb_in_slice << u1_mbaff;

    do
    {
        UWORD8 u1_mb_type;
        ps_dec->pv_prev_mb_parse_tu_coeff_data = ps_dec->pv_parse_tu_coeff_data;
        if(i2_cur_mb_addr > ps_dec->ps_cur_sps->u2_max_mb_addr)
        {
            break;
        }

        ps_cur_mb_info = ps_dec->ps_nmb_info + u1_num_mbs;
        ps_svc_cur_mb_info = ps_svc_lyr_dec->ps_svc_nmb_info + u1_num_mbs;
        ps_dec->u4_num_mbs_cur_nmb = u1_num_mbs;
        ps_dec->u4_num_pmbair = (u1_num_mbs >> u1_mbaff);
        ps_cur_mb_info->u1_end_of_slice = 0;
        /***************************************************************/
        /* Get the required information for decoding of MB             */
        /* mb_x, mb_y , neighbour availablity,                         */
        /***************************************************************/
        ps_dec->pf_get_mb_info(ps_dec, i2_cur_mb_addr, ps_cur_mb_info, 0);
        ps_svc_cur_mb_info->u1_crop_window_flag =
            *(ps_svc_lyr_dec->pu1_crop_wnd_flag + ps_cur_mb_info->u2_mbx +
              (ps_cur_mb_info->u2_mby * ps_dec->u2_frm_wd_in_mbs));

        /***************************************************************/
        /* Set the deblocking parameters for this MB                   */
        /***************************************************************/
        ps_cur_deblk_mb = ps_dec->ps_deblk_mbn + u1_num_mbs;
        if(ps_dec->u4_app_disable_deblk_frm == 0)
            ih264d_set_deblocking_parameters(ps_cur_deblk_mb, ps_slice,
                                             ps_dec->u1_mb_ngbr_availablity,
                                             ps_dec->u1_cur_mb_fld_dec_flag);

        ps_svc_cur_mb_info->u1_base_mode_flag = 0;
        if(ps_svc_cur_mb_info->u1_crop_window_flag &&
           ps_svc_slice_params->u1_adaptive_base_mode_flag)
        {
            ps_svc_cur_mb_info->u1_base_mode_flag = ih264d_get_bit_h264(ps_bitstrm);
            COPYTHECONTEXT("SVC ext: u1_base_mode_flag", ps_cur_mb_info->u1_base_mode_flag);
        }
        else if(ps_svc_cur_mb_info->u1_crop_window_flag)
        {
            ps_svc_cur_mb_info->u1_base_mode_flag = ps_svc_slice_params->u1_default_base_mode_flag;
        }

        /**************************************************************/
        /* Macroblock Layer Begins, Decode the u1_mb_type                */
        /**************************************************************/
        if(!ps_svc_cur_mb_info->u1_base_mode_flag) /* Inlined ih264d_uev */
        {
            UWORD32 u4_bitstream_offset = *pu4_bitstrm_ofst;
            UWORD32 u4_word, u4_ldz, u4_temp;

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
            *pu4_bitstrm_ofst = u4_bitstream_offset;
            u4_temp = ((1 << u4_ldz) + u4_word - 1);
            if(u4_temp > 25) return ERROR_MB_TYPE;
            u1_mb_type = u4_temp;

            /* Inlined ih264d_uev */
            ps_cur_mb_info->u1_mb_type = u1_mb_type;
            COPYTHECONTEXT("u1_mb_type", u1_mb_type);

            ps_cur_deblk_mb->u1_mb_type = ps_cur_deblk_mb->u1_mb_type | D_INTRA_MB;
        }
        else
        {
            ps_cur_mb_info->u1_mb_type = MB_INFER;
            ps_cur_deblk_mb->u1_mb_type = ps_cur_deblk_mb->u1_mb_type | D_INTRA_IBL;
        }

        /**************************************************************/
        /* Parse Macroblock data                                      */
        /**************************************************************/
        u1_mb_type = ps_cur_mb_info->u1_mb_type;
        if(25 == u1_mb_type)
        {
            /* I_PCM_MB */
            ps_cur_mb_info->ps_curmb->u1_mb_type = I_PCM_MB;
            ret = ih264d_parse_ipcm_mb(ps_dec, ps_cur_mb_info, u1_num_mbs);
            if(ret != OK) return ret;
            ps_cur_deblk_mb->u1_mb_qp = 0;
        }
        else
        {
            ret = isvcd_parse_imb_cavlc(ps_svc_lyr_dec, ps_cur_mb_info, ps_svc_cur_mb_info,
                                        u1_num_mbs, u1_mb_type);
            if(ret != OK) return ret;
            ps_cur_deblk_mb->u1_mb_qp = ps_dec->u1_qp;
        }

        if(ps_dec->u1_enable_mb_info)
        {
            ih264d_populate_mb_info_map(ps_dec, ps_cur_mb_info, ps_cur_mb_info->u2_mbx << 1,
                                        ps_cur_mb_info->u2_mby << 1, ps_cur_deblk_mb->u1_mb_qp);
        }
        uc_more_data_flag = MORE_RBSP_DATA(ps_bitstrm);

        if(u1_mbaff)
        {
            ih264d_update_mbaff_left_nnz(ps_dec, ps_cur_mb_info);
            if(!uc_more_data_flag && (0 == (i2_cur_mb_addr & 1)))
            {
                return ERROR_EOB_FLUSHBITS_T;
            }
        }
        /**************************************************************/
        /* Get next Macroblock address                                */
        /**************************************************************/
        i2_cur_mb_addr++;
        /* Store the colocated information */
        {
            mv_pred_t *ps_mv_nmb_start = ps_dec->ps_mv_cur + (u1_num_mbs << 4);
            mv_pred_t s_mvPred = {{0, 0, 0, 0}, {-1, -1}, 0, 0};
            if(ps_mv_nmb_start)
            {
                ih264d_rep_mv_colz(ps_dec, &s_mvPred, ps_mv_nmb_start, 0,
                                   (UWORD8) (ps_dec->u1_cur_mb_fld_dec_flag << 1), 4, 4);
            }
            else
            {
                return NOT_OK;
            }
        }
        u1_num_mbs++;

        /****************************************************************/
        /* Check for End Of Row                                         */
        /****************************************************************/
        u1_num_mbs_next = i2_pic_wdin_mbs - ps_dec->u2_mbx - 1;
        u1_end_of_row = (!u1_num_mbs_next) && (!(u1_mbaff && (u1_num_mbs & 0x01)));
        u1_tfr_n_mb =
            (u1_num_mbs == ps_dec->u1_recon_mb_grp) || u1_end_of_row || (!uc_more_data_flag);
        ps_cur_mb_info->u1_end_of_slice = (!uc_more_data_flag);

        if(u1_tfr_n_mb || (!uc_more_data_flag))
        {
            if(ps_dec->u1_separate_parse)
            {
                ih264d_parse_tfr_nmb(ps_dec, u1_mb_idx, u1_num_mbs, u1_num_mbs_next, u1_tfr_n_mb,
                                     u1_end_of_row);
                ps_dec->ps_nmb_info += u1_num_mbs;
                ps_svc_lyr_dec->ps_svc_nmb_info += u1_num_mbs;
            }
            else
            {
                ret = isvcd_decode_recon_tfr_nmb_non_base_lyr(ps_svc_lyr_dec, u1_mb_idx, u1_num_mbs,
                                                              u1_num_mbs_next, u1_tfr_n_mb,
                                                              u1_end_of_row);
                if(ret != OK) return ret;
            }
            ps_dec->u2_total_mbs_coded += u1_num_mbs;
            if(u1_tfr_n_mb) u1_num_mbs = 0;
            u1_mb_idx = u1_num_mbs;
            ps_dec->u1_mb_idx = u1_num_mbs;
        }
    } while(uc_more_data_flag);

    ps_dec->u4_num_mbs_cur_nmb = 0;
    ps_dec->ps_cur_slice->u4_mbs_in_slice = i2_cur_mb_addr - (u2_first_mb_in_slice << u1_mbaff);

    return ret;
}

/*!
**************************************************************************
* \if Function name : ParseIMbCab \endif
*
* \brief
*    This function parses CABAC syntax of a I MB. If 16x16 Luma DC transform
*    is also done here. Transformed Luma DC values are copied in their
*    0th pixel location of corrosponding CoeffBlock.
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_imb_cabac(svc_dec_lyr_struct_t *ps_svc_lyr_dec, dec_mb_info_t *ps_cur_mb_info,
                             dec_svc_mb_info_t *ps_svc_cur_mb_info, UWORD8 u1_mb_type)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    WORD8 i1_delta_qp;
    UWORD8 u1_cbp;
    UWORD8 u1_offset;
    /* Variables for handling Cabac contexts */
    ctxt_inc_mb_info_t *p_curr_ctxt = ps_dec->ps_curr_ctxt_mb_info;
    ctxt_inc_mb_info_t *ps_left_ctxt = ps_dec->p_left_ctxt_mb_info;
    dec_bit_stream_t *const ps_bitstrm = ps_dec->ps_bitstrm;
    bin_ctxt_model_t *p_bin_ctxt;
    UWORD8 u1_intra_chrom_pred_mode;
    UWORD8 u1_dc_block_flag = 0;
    WORD32 ret;

    ps_cur_mb_info->u1_yuv_dc_block_flag = 0;

    if(ps_left_ctxt == ps_dec->ps_def_ctxt_mb_info)
    {
        ps_dec->pu1_left_yuv_dc_csbp[0] = 0xf;
    }

    if(ps_dec->ps_cur_slice->u1_slice_type != I_SLICE)
    {
        MEMSET_16BYTES(&ps_dec->pu1_left_mv_ctxt_inc[0][0], 0);
        *((UWORD32 *) ps_dec->pi1_left_ref_idx_ctxt_inc) = 0;
        MEMSET_16BYTES(p_curr_ctxt->u1_mv, 0);
        memset(p_curr_ctxt->i1_ref_idx, 0, 4);
    }

    /* default */
    ps_cur_mb_info->ps_curmb->u1_mb_type = P_MB;
    if(!ps_svc_cur_mb_info->u1_base_mode_flag)
    {
        if(u1_mb_type == I_4x4_MB)
        {
            ps_cur_mb_info->ps_curmb->u1_mb_type = I_4x4_MB;
            p_curr_ctxt->u1_mb_type = CAB_I4x4;
            u1_offset = 0;
            ps_cur_mb_info->u1_tran_form8x8 = 0;
            ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;
            /*--------------------------------------------------------------------*/
            /* Read transform_size_8x8_flag if present                            */
            /*--------------------------------------------------------------------*/
            if(ps_dec->ps_cur_pps->i4_transform_8x8_mode_flag)
            {
                ps_cur_mb_info->u1_tran_form8x8 =
                    ih264d_parse_transform8x8flag_cabac(ps_dec, ps_cur_mb_info);
                COPYTHECONTEXT("transform_size_8x8_flag", ps_cur_mb_info->u1_tran_form8x8);
                p_curr_ctxt->u1_transform8x8_ctxt = ps_cur_mb_info->u1_tran_form8x8;
                ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = ps_cur_mb_info->u1_tran_form8x8;
            }
            else
            {
                p_curr_ctxt->u1_transform8x8_ctxt = 0;
            }

            /*--------------------------------------------------------------------*/
            /* Read the IntraPrediction modes for LUMA                            */
            /*--------------------------------------------------------------------*/
            if(!ps_cur_mb_info->u1_tran_form8x8)
            {
                UWORD8 *pu1_temp;
                ih264d_read_intra_pred_modes_cabac(ps_dec,
                                                   ((UWORD8 *) ps_dec->pv_parse_tu_coeff_data),
                                                   ((UWORD8 *) ps_dec->pv_parse_tu_coeff_data + 16),
                                                   ps_cur_mb_info->u1_tran_form8x8);
                pu1_temp = (UWORD8 *) ps_dec->pv_parse_tu_coeff_data;
                pu1_temp += 32;
                ps_dec->pv_parse_tu_coeff_data = (void *) pu1_temp;
            }
            else
            {
                UWORD8 *pu1_temp;
                ih264d_read_intra_pred_modes_cabac(ps_dec,
                                                   ((UWORD8 *) ps_dec->pv_parse_tu_coeff_data),
                                                   ((UWORD8 *) ps_dec->pv_parse_tu_coeff_data + 4),
                                                   ps_cur_mb_info->u1_tran_form8x8);
                pu1_temp = (UWORD8 *) ps_dec->pv_parse_tu_coeff_data;
                pu1_temp += 8;
                ps_dec->pv_parse_tu_coeff_data = (void *) pu1_temp;
            }
            /*--------------------------------------------------------------------*/
            /* Read the IntraPrediction mode for CHROMA                           */
            /*--------------------------------------------------------------------*/
            u1_intra_chrom_pred_mode = ih264d_parse_chroma_pred_mode_cabac(ps_dec);
            COPYTHECONTEXT("intra_chroma_pred_mode", u1_intra_chrom_pred_mode);
            p_curr_ctxt->u1_intra_chroma_pred_mode = ps_cur_mb_info->u1_chroma_pred_mode =
                u1_intra_chrom_pred_mode;

            /*--------------------------------------------------------------------*/
            /* Read the Coded block pattern                                       */
            /*--------------------------------------------------------------------*/
            u1_cbp = ih264d_parse_ctx_cbp_cabac(ps_dec);
            COPYTHECONTEXT("coded_block_pattern", u1_cbp);
            ps_cur_mb_info->u1_cbp = u1_cbp;
            p_curr_ctxt->u1_cbp = u1_cbp;

            /*--------------------------------------------------------------------*/
            /* Read mb_qp_delta                                                   */
            /*--------------------------------------------------------------------*/
            if(ps_cur_mb_info->u1_cbp)
            {
                ret = ih264d_parse_mb_qp_delta_cabac(ps_dec, &i1_delta_qp);
                if(ret != OK) return ret;
                COPYTHECONTEXT("mb_qp_delta", i1_delta_qp);
                if(i1_delta_qp != 0)
                {
                    ret = ih264d_update_qp(ps_dec, i1_delta_qp);
                    if(ret != OK) return ret;
                }
            }
            else
                ps_dec->i1_prev_mb_qp_delta = 0;
            p_curr_ctxt->u1_yuv_dc_csbp &= 0xFE;
        }
        else
        {
            u1_offset = 1;
            ps_cur_mb_info->ps_curmb->u1_mb_type = I_16x16_MB;
            p_curr_ctxt->u1_mb_type = CAB_I16x16;
            ps_cur_mb_info->u1_tran_form8x8 = 0;
            p_curr_ctxt->u1_transform8x8_ctxt = 0;
            ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;
            /*--------------------------------------------------------------------*/
            /* Read the IntraPrediction mode for CHROMA                           */
            /*--------------------------------------------------------------------*/
            u1_intra_chrom_pred_mode = ih264d_parse_chroma_pred_mode_cabac(ps_dec);
            if(u1_intra_chrom_pred_mode > 3) return ERROR_CHROMA_PRED_MODE;

            COPYTHECONTEXT("Chroma intra_chroma_pred_mode pred mode", u1_intra_chrom_pred_mode);
            p_curr_ctxt->u1_intra_chroma_pred_mode = ps_cur_mb_info->u1_chroma_pred_mode =
                u1_intra_chrom_pred_mode;

            /*--------------------------------------------------------------------*/
            /* Read the Coded block pattern                                       */
            /*--------------------------------------------------------------------*/
            u1_cbp = gau1_ih264d_cbp_tab[(u1_mb_type - 1) >> 2];
            ps_cur_mb_info->u1_cbp = u1_cbp;
            p_curr_ctxt->u1_cbp = u1_cbp;

            /*--------------------------------------------------------------------*/
            /* Read mb_qp_delta                                                   */
            /*--------------------------------------------------------------------*/
            ret = ih264d_parse_mb_qp_delta_cabac(ps_dec, &i1_delta_qp);
            if(ret != OK) return ret;
            COPYTHECONTEXT("mb_qp_delta", i1_delta_qp);
            if(i1_delta_qp != 0)
            {
                ret = ih264d_update_qp(ps_dec, i1_delta_qp);
                if(ret != OK) return ret;
            }

            {
                WORD16 i_scaleFactor;
                WORD16 *pi2_scale_matrix_ptr;
                /*******************************************************************/
                /* for luma DC coefficients the scaling is done during the parsing */
                /* to preserve the precision                                       */
                /*******************************************************************/
                if(ps_dec->s_high_profile.u1_scaling_present)
                {
                    pi2_scale_matrix_ptr = ps_dec->s_high_profile.i2_scalinglist4x4[0];
                }
                else
                {
                    i_scaleFactor = 16;
                    pi2_scale_matrix_ptr = &i_scaleFactor;
                }
                {
                    ctxt_inc_mb_info_t *ps_top_ctxt = ps_dec->p_top_ctxt_mb_info;
                    UWORD8 uc_a, uc_b;
                    UWORD32 u4_ctx_inc;
                    INC_SYM_COUNT(&(ps_dec->s_cab_dec_env));

                    /* if MbAddrN not available then CondTermN = 1 */
                    uc_b = ((ps_top_ctxt->u1_yuv_dc_csbp) & 0x01);

                    /* if MbAddrN not available then CondTermN = 1 */
                    uc_a = ((ps_dec->pu1_left_yuv_dc_csbp[0]) & 0x01);

                    u4_ctx_inc = (uc_a + (uc_b << 1));
                    {
                        WORD16 pi2_dc_coef[16] = {0};
                        tu_sblk4x4_coeff_data_t *ps_tu_4x4 =
                            (tu_sblk4x4_coeff_data_t *) ps_dec->pv_parse_tu_coeff_data;
                        WORD16 *pi2_coeff_block = (WORD16 *) ps_dec->pv_parse_tu_coeff_data;

                        p_bin_ctxt = (ps_dec->p_cbf_t[LUMA_DC_CTXCAT]) + u4_ctx_inc;

                        u1_dc_block_flag = ih264d_read_coeff4x4_cabac(
                            ps_bitstrm, LUMA_DC_CTXCAT,
                            ps_dec->p_significant_coeff_flag_t[LUMA_DC_CTXCAT], ps_dec, p_bin_ctxt);

                        /* Store coded_block_flag */
                        p_curr_ctxt->u1_yuv_dc_csbp &= 0xFE;
                        p_curr_ctxt->u1_yuv_dc_csbp |= u1_dc_block_flag;
                        if(u1_dc_block_flag)
                        {
                            WORD32 pi4_tmp[16] = {0};
                            memset(pi2_dc_coef, 0, sizeof(pi2_dc_coef));
                            ih264d_unpack_coeff4x4_dc_4x4blk(ps_tu_4x4, pi2_dc_coef,
                                                             ps_dec->pu1_inv_scan);

                            PROFILE_DISABLE_IQ_IT_RECON()
                            ps_dec->pf_ihadamard_scaling_4x4(
                                pi2_dc_coef, pi2_coeff_block, ps_dec->pu2_quant_scale_y,
                                (UWORD16 *) pi2_scale_matrix_ptr, ps_dec->u1_qp_y_div6, pi4_tmp);
                            pi2_coeff_block += 16;
                            ps_dec->pv_parse_tu_coeff_data = (void *) pi2_coeff_block;
                            SET_BIT(ps_cur_mb_info->u1_yuv_dc_block_flag, 0);
                        }
                    }
                }
            }
        }
        ps_dec->pu1_left_yuv_dc_csbp[0] &= 0x6;
        ps_dec->pu1_left_yuv_dc_csbp[0] |= u1_dc_block_flag;
    }
    else
    {
        u1_offset = 0;
        /*--------------------------------------------------------------------*/
        /* Read the Coded block pattern                                       */
        /*--------------------------------------------------------------------*/
        u1_cbp = ih264d_parse_ctx_cbp_cabac(ps_dec);
        COPYTHECONTEXT("coded_block_pattern", u1_cbp);

        ps_cur_mb_info->u1_cbp = u1_cbp;
        p_curr_ctxt->u1_cbp = u1_cbp;
        p_curr_ctxt->u1_yuv_dc_csbp &= 0xFE;
        ps_dec->pu1_left_yuv_dc_csbp[0] &= 0x6;
        p_curr_ctxt->u1_mb_type = CAB_INFERRED;
        p_curr_ctxt->u1_intra_chroma_pred_mode = 0;
        ps_cur_mb_info->u1_tran_form8x8 = 0;
        ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = 0;

        /*--------------------------------------------------------------------*/
        /* Read transform_size_8x8_flag if present                            */
        /*--------------------------------------------------------------------*/
        if((ps_dec->ps_cur_pps->i4_transform_8x8_mode_flag) && (u1_cbp & 0xf))
        {
            ps_cur_mb_info->u1_tran_form8x8 =
                ih264d_parse_transform8x8flag_cabac(ps_dec, ps_cur_mb_info);
            COPYTHECONTEXT("transform_size_8x8_flag", ps_cur_mb_info->u1_tran_form8x8);
            p_curr_ctxt->u1_transform8x8_ctxt = ps_cur_mb_info->u1_tran_form8x8;
            ps_cur_mb_info->ps_curmb->u1_tran_form8x8 = ps_cur_mb_info->u1_tran_form8x8;
        }
        else
        {
            p_curr_ctxt->u1_transform8x8_ctxt = 0;
        }

        /*--------------------------------------------------------------------*/
        /* Read mb_qp_delta                                                   */
        /*--------------------------------------------------------------------*/
        if(ps_cur_mb_info->u1_cbp)
        {
            ret = ih264d_parse_mb_qp_delta_cabac(ps_dec, &i1_delta_qp);
            if(ret != OK) return ret;
            COPYTHECONTEXT("mb_qp_delta", i1_delta_qp);
            if(i1_delta_qp != 0)
            {
                ret = ih264d_update_qp(ps_dec, i1_delta_qp);
                if(ret != OK) return ret;
            }
        }
        else
            ps_dec->i1_prev_mb_qp_delta = 0;
    }

    ih264d_parse_residual4x4_cabac(ps_dec, ps_cur_mb_info, u1_offset);
    if(EXCEED_OFFSET(ps_bitstrm)) return ERROR_EOB_TERMINATE_T;
    return OK;
}

/*!
**************************************************************************
* \if Function name : isvcd_decode_islice \endif
*
* \brief
*    Decodes an I Slice
*
*
* \return
*    0 on Success and Error code otherwise
**************************************************************************
*/
WORD32 isvcd_parse_islice(svc_dec_lyr_struct_t *ps_svc_lyr_dec, UWORD16 u2_first_mb_in_slice)
{
    dec_struct_t *ps_dec = &ps_svc_lyr_dec->s_dec;
    dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
    dec_slice_params_t *ps_slice = ps_dec->ps_cur_slice;
    UWORD32 *pu4_bitstrm_buf = ps_dec->ps_bitstrm->pu4_buffer;
    UWORD32 *pu4_bitstrm_ofst = &ps_dec->ps_bitstrm->u4_ofst;
    UWORD32 u4_temp;
    WORD32 i_temp;
    WORD32 ret;

    /*--------------------------------------------------------------------*/
    /* Read remaining contents of the slice header                        */
    /*--------------------------------------------------------------------*/
    /* dec_ref_pic_marking function */
    /* G050 */
    if(ps_slice->u1_nal_ref_idc != 0)
    {
        if(!ps_dec->ps_dpb_cmds->u1_dpb_commands_read)
        {
            dec_pic_params_t *ps_pps = ps_dec->ps_cur_pps;
            dec_seq_params_t *ps_sps_tmp = ps_pps->ps_sps;
            UWORD8 u1_nal_unit_type_tmp = ps_dec->u1_nal_unit_type;

            ps_pps->ps_sps = ps_dec->ps_cur_sps;
            if(ps_svc_lyr_dec->ps_nal_svc_ext->u1_idr_flag)
                ps_dec->u1_nal_unit_type = IDR_SLICE_NAL;

            i_temp = ih264d_read_mmco_commands(ps_dec);
            ps_pps->ps_sps = ps_sps_tmp;
            ps_dec->u1_nal_unit_type = u1_nal_unit_type_tmp;

            if(i_temp < 0)
            {
                return ERROR_DBP_MANAGER_T;
            }
            ps_dec->u4_bitoffset = i_temp;
        }
        else
            ps_dec->ps_bitstrm->u4_ofst += ps_dec->u4_bitoffset;
    }

    {
        /* G050 */
        /* Read slice_qp_delta */
        WORD64 i8_temp =
            (WORD64) ps_pps->u1_pic_init_qp + ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        if((i8_temp < MIN_H264_QP) || (i8_temp > MAX_H264_QP)) return ERROR_INV_RANGE_QP_T;
        ps_slice->u1_slice_qp = (UWORD8) i8_temp;
        COPYTHECONTEXT("SH: slice_qp_delta", ps_slice->u1_slice_qp - ps_pps->u1_pic_init_qp);
    }
    if(ps_pps->u1_deblocking_filter_parameters_present_flag == 1)
    {
        u4_temp = ih264d_uev(pu4_bitstrm_ofst, pu4_bitstrm_buf);
        COPYTHECONTEXT("SH: disable_deblocking_filter_idc", u4_temp);

        if(u4_temp > SLICE_BOUNDARY_DBLK_DISABLED)
        {
            return ERROR_INV_SLICE_HDR_T;
        }
        ps_slice->u1_disable_dblk_filter_idc = u4_temp;
        if(u4_temp != 1)
        {
            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf) << 1;
            if((MIN_DBLK_FIL_OFF > i_temp) || (i_temp > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice->i1_slice_alpha_c0_offset = i_temp;
            COPYTHECONTEXT("SH: slice_alpha_c0_offset_div2",
                           ps_slice->i1_slice_alpha_c0_offset >> 1);

            i_temp = ih264d_sev(pu4_bitstrm_ofst, pu4_bitstrm_buf) << 1;
            if((MIN_DBLK_FIL_OFF > i_temp) || (i_temp > MAX_DBLK_FIL_OFF))
            {
                return ERROR_INV_SLICE_HDR_T;
            }
            ps_slice->i1_slice_beta_offset = i_temp;
            COPYTHECONTEXT("SH: slice_beta_offset_div2", ps_slice->i1_slice_beta_offset >> 1);
        }
        else
        {
            ps_slice->i1_slice_alpha_c0_offset = 0;
            ps_slice->i1_slice_beta_offset = 0;
        }
    }
    else
    {
        ps_slice->u1_disable_dblk_filter_idc = 0;
        ps_slice->i1_slice_alpha_c0_offset = 0;
        ps_slice->i1_slice_beta_offset = 0;
    }

    /* Initialization to check if number of motion vector per 2 Mbs */
    /* are exceeding the range or not */
    ps_dec->u2_mv_2mb[0] = 0;
    ps_dec->u2_mv_2mb[1] = 0;

    ret = isvcd_parse_interlayer_resamp_func_init(ps_svc_lyr_dec, u2_first_mb_in_slice);
    if(ret != OK)
    {
        return ERROR_CORRUPTED_SLICE;
    }
    /*set slice header cone to 2 ,to indicate  correct header*/
    ps_dec->u1_slice_header_done = 2;
    if(ps_pps->u1_entropy_coding_mode)
    {
        SWITCHOFFTRACE;
        SWITCHONTRACECABAC;
        if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
        {
            ps_dec->pf_get_mb_info = ih264d_get_mb_info_cabac_mbaff;
        }
        else
            ps_dec->pf_get_mb_info = isvcd_get_mb_info_cabac_nonmbaff;

        ret = isvcd_parse_islice_data_cabac(ps_svc_lyr_dec, ps_slice, u2_first_mb_in_slice);
        if(ret != OK) return ret;
        SWITCHONTRACE;
        SWITCHOFFTRACECABAC;
    }
    else
    {
        if(ps_dec->ps_cur_slice->u1_mbaff_frame_flag)
        {
            ps_dec->pf_get_mb_info = ih264d_get_mb_info_cavlc_mbaff;
        }
        else
            ps_dec->pf_get_mb_info = isvcd_get_mb_info_cavlc_nonmbaff;
        ret = isvcd_parse_islice_data_cavlc(ps_svc_lyr_dec, ps_slice, u2_first_mb_in_slice);
        if(ret != OK) return ret;
    }

    return OK;
}
