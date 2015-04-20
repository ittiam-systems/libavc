/******************************************************************************
 *
 * Copyright (C) 2015 The Android Open Source Project
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
 * \file ih264d_thread_parse_decode.c
 *
 * \brief
 *    Contains routines that for multi-thread decoder
 *
 * Detailed_description
 *
 * \date
 *    20/02/2012
 *
 * \author  ZR
 **************************************************************************
 */

#include "ih264d_error_handler.h"
#include "ih264d_debug.h"
#include "ithread.h"
#include <string.h>
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "ih264d_tables.h"
#include "ih264d_structs.h"
#include "ih264d_defs.h"
#include "ih264d_mb_utils.h"
#include "ih264d_thread_parse_decode.h"
#include "ih264d_inter_pred.h"

#include "ih264d_process_pslice.h"
#include "ih264d_process_intra_mb.h"
#include "ih264d_deblocking.h"
#include "ih264d_format_conv.h"

void ih264d_deblock_mb_level(dec_struct_t *ps_dec,
                             dec_mb_info_t *ps_cur_mb_info,
                             UWORD32 nmb_index);

void ih264d_copy_intra_pred_line(dec_struct_t *ps_dec,
                                 dec_mb_info_t *ps_cur_mb_info,
                                 UWORD32 nmb_index);

void ih264d_parse_tfr_nmb(dec_struct_t * ps_dec,
                          UWORD8 u1_mb_idx,
                          UWORD8 u1_num_mbs,
                          UWORD8 u1_num_mbs_next,
                          UWORD8 u1_tfr_n_mb,
                          UWORD8 u1_end_of_row)
{
    WORD32 i, u4_mb_num;

    const UWORD32 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    UWORD32 u4_n_mb_start;

    UNUSED(u1_mb_idx);
    UNUSED(u1_num_mbs_next);
    if(u1_tfr_n_mb)
    {


        u4_n_mb_start = (ps_dec->u2_cur_mb_addr + 1) - u1_num_mbs;

        // copy into s_frmMbInfo

        u4_mb_num = u4_n_mb_start;
        ps_dec->ps_parse_cur_slice->u4_num_mbs_done_in_slice += u1_num_mbs;
        u4_mb_num = (ps_dec->u2_cur_mb_addr + 1) - u1_num_mbs;

        for(i = 0; i < u1_num_mbs; i++)
        {
            DATA_SYNC();
            UPDATE_SLICE_NUM_MAP(ps_dec->pu2_slice_num_map, u4_mb_num,
                                 ps_dec->u2_cur_slice_num);
            UPDATE_MB_MAP_MBNUM_BYTE(ps_dec->pu1_dec_mb_map, u4_mb_num);

            u4_mb_num++;
        }

        DATA_SYNC();
        /****************************************************************/
        /* Check for End Of Row in Next iteration                       */
        /****************************************************************/

        /****************************************************************/
        /* Transfer the Following things                                */
        /* N-Mb DeblkParams Data    ( To Ext DeblkParams Buffer )       */
        /* N-Mb Recon Data          ( To Ext Frame Buffer )             */
        /* N-Mb Intrapredline Data  ( Updated Internally)               */
        /* N-Mb MV Data             ( To Ext MV Buffer )                */
        /* N-Mb MVTop/TopRight Data ( To Int MV Top Scratch Buffers)    */
        /****************************************************************/

        /* Swap top and current pointers */

        ps_dec->s_tran_addrecon_parse.pu1_dest_y +=
                        ps_dec->s_tran_addrecon_parse.u4_inc_y[u1_end_of_row];
        ps_dec->s_tran_addrecon_parse.pu1_dest_u +=
                        ps_dec->s_tran_addrecon_parse.u4_inc_uv[u1_end_of_row];
        ps_dec->s_tran_addrecon_parse.pu1_dest_v +=
                        ps_dec->s_tran_addrecon_parse.u4_inc_uv[u1_end_of_row];

        if(u1_end_of_row)
        {
            UWORD16 u2_mb_y;
            UWORD32 u4_frame_stride, y_offset;

            ps_dec->ps_top_mb_row = ps_dec->ps_cur_mb_row;
            ps_dec->ps_cur_mb_row += ((ps_dec->u2_frm_wd_in_mbs) << u1_mbaff);

            u2_mb_y = ps_dec->u2_mby + (1 + u1_mbaff);
            u4_frame_stride = ps_dec->u2_frm_wd_y
                            << ps_dec->ps_cur_slice->u1_field_pic_flag;
            y_offset = (u2_mb_y * u4_frame_stride) << 4;
            ps_dec->s_tran_addrecon_parse.pu1_dest_y =
                            ps_dec->s_cur_pic.pu1_buf1 + y_offset;

            u4_frame_stride = ps_dec->u2_frm_wd_uv
                            << ps_dec->ps_cur_slice->u1_field_pic_flag;
            y_offset = (u2_mb_y * u4_frame_stride) << 3;
            ps_dec->s_tran_addrecon_parse.pu1_dest_u =
                            ps_dec->s_cur_pic.pu1_buf2 + y_offset;
            ps_dec->s_tran_addrecon_parse.pu1_dest_v =
                            ps_dec->s_cur_pic.pu1_buf3 + y_offset;

        }

        ps_dec->ps_deblk_mbn += u1_num_mbs;

        /*
         * The Slice boundary is also a valid condition to transfer. So recalculate
         * the Left increment, in case the number of MBs is lesser than the
         * N MB value. c_numMbs will be equal to N of N MB if the entire N Mb is
         * decoded.
         */
        ps_dec->s_tran_addrecon.u2_mv_left_inc = ((u1_num_mbs >> u1_mbaff) - 1)
                        << (4 + u1_mbaff);
        ps_dec->s_tran_addrecon.u2_mv_top_left_inc = (u1_num_mbs << 2) - 1
                        - (u1_mbaff << 2);

        /* reassign left MV and cur MV pointers */
        ps_dec->ps_mv_left = ps_dec->ps_mv_cur
                        + ps_dec->s_tran_addrecon.u2_mv_left_inc;





        ps_dec->ps_mv_cur += (u1_num_mbs << 4);
        ps_dec->u4_num_mbs_prev_nmb = u1_num_mbs;


        ps_dec->u4_dma_buf_idx = 0;

    }
}

void ih264d_decode_tfr_nmb(dec_struct_t * ps_dec,
                           UWORD8 u1_num_mbs,
                           UWORD8 u1_num_mbs_next,
                           UWORD8 u1_end_of_row)
{

    UWORD32 u1_end_of_row_next;

    const UWORD32 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;

    /****************************************************************/
    /* Check for End Of Row in Next iteration                       */
    /****************************************************************/
    u1_end_of_row_next =
                    u1_num_mbs_next
                                    && ((u1_num_mbs_next)
                                                    <= (ps_dec->u1_recon_mb_grp
                                                                    >> u1_mbaff));

    /****************************************************************/
    /* Transfer the Following things                                */
    /* N-Mb DeblkParams Data    ( To Ext DeblkParams Buffer )       */
    /* N-Mb Recon Data          ( To Ext Frame Buffer )             */
    /* N-Mb Intrapredline Data  ( Updated Internally)               */
    /* N-Mb MV Data             ( To Ext MV Buffer )                */
    /* N-Mb MVTop/TopRight Data ( To Int MV Top Scratch Buffers)    */
    /****************************************************************/
    if(u1_end_of_row)
    {
        ps_dec->i2_dec_thread_mb_y += (1 << u1_mbaff);
    }
    ih264d_transfer_mb_group_data(ps_dec, u1_num_mbs, u1_end_of_row,
                                  u1_end_of_row_next);

    if(u1_end_of_row)
    {
        /* Reset the N-Mb Recon Buf Index to default Values */
        ps_dec->u2_mb_group_cols_y1 = ps_dec->u2_mb_group_cols_y;
        ps_dec->u2_mb_group_cols_cr1 = ps_dec->u2_mb_group_cols_cr;
    }
    /* If next N-Mb Group is the EndOfRow, set the N-Mb Recon Buf Index */
    else if(u1_end_of_row_next)
    {
        ps_dec->u2_mb_group_cols_y1 = (u1_num_mbs_next << 4) + 8;
        ps_dec->u2_mb_group_cols_cr1 = (u1_num_mbs_next << 3) + 8;
    }
}

WORD32 ih264d_decode_recon_tfr_nmb_thread(dec_struct_t * ps_dec, UWORD8 u1_num_mbs, // number of MBs loop should run
                                        UWORD8 u1_num_mbs_next,
                                        UWORD8 u1_end_of_row)
{
    WORD32 i,j;
    dec_mb_info_t * ps_cur_mb_info;
    UWORD32 u4_update_mbaff = 0;
    const UWORD32 u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;
    UWORD32 u1_slice_type, u1_B;
    WORD32 u1_skip_th;
    UWORD32 u1_ipcm_th;
    UWORD32 u4_cond;
    UWORD16 u2_slice_num,u2_cur_dec_mb_num;
    WORD32 ret;

    u1_slice_type = ps_dec->ps_decode_cur_slice->slice_type;

    u1_B = (u1_slice_type == B_SLICE);

    u1_skip_th =
                    ((u1_slice_type != I_SLICE) ?
                                    (u1_B ? B_8x8 : PRED_8x8R0) : -1);

    u1_ipcm_th = ((u1_slice_type != I_SLICE) ? (u1_B ? 23 : 5) : 0);

    u2_cur_dec_mb_num = ps_dec->cur_dec_mb_num;

    /* N Mb MC Loop */
    for(i = 0; i < u1_num_mbs; i++)
    {
        DATA_SYNC();

        // check dec_mb_map
        UWORD32 yield_cnt = 0, u4_max_addr;

        u4_max_addr = ps_dec->ps_cur_sps->u2_max_mb_addr;
        while(1)
        {
            UWORD32 u4_mb_num = u2_cur_dec_mb_num;

            /*introducing 1 MB delay*/
            if(u4_mb_num < u4_max_addr)
                u4_mb_num = u4_mb_num + 1;

            CHECK_MB_MAP_BYTE(u4_mb_num, ps_dec->pu1_dec_mb_map, u4_cond);
            if(u4_cond)
            {
                break;
            }
            else
            {

                {
                    NOP(128);

                }

                DEBUG_THREADS_PRINTF("waiting for mb mapcur_dec_mb_num = %d,ps_dec->u2_cur_mb_addr  = %d\n",u2_cur_dec_mb_num,
                                ps_dec->u2_cur_mb_addr);

            }
        }

        GET_SLICE_NUM_MAP(ps_dec->pu2_slice_num_map, u2_cur_dec_mb_num,
                          u2_slice_num);

        if(u2_slice_num != ps_dec->u2_cur_slice_num_dec_thread)
        {
            ps_dec->u4_cur_slice_decode_done = 1;
            break;
        }

        ps_cur_mb_info = &ps_dec->ps_frm_mb_info[u2_cur_dec_mb_num
                        & PD_MB_BUF_SIZE_MOD];

        ps_dec->u4_dma_buf_idx = 0;
        ps_dec->u4_pred_info_idx = 0;

        if(ps_cur_mb_info->u1_mb_type <= u1_skip_th)
        {

            {
                WORD32 pred_cnt = 0;
                pred_info_pkd_t *ps_pred_pkd;
                UWORD32 u4_pred_info_pkd_idx;
                WORD8 i1_pred;

                u4_pred_info_pkd_idx = ps_cur_mb_info->u4_pred_info_pkd_idx;

                while(pred_cnt < ps_cur_mb_info->u1_num_pred_parts)
                {

                    ps_pred_pkd = ps_dec->ps_pred_pkd + u4_pred_info_pkd_idx;


                    ps_dec->p_form_mb_part_info_thread(ps_pred_pkd,ps_dec,
                                         ps_cur_mb_info->u2_mbx,ps_cur_mb_info->u2_mby,(i >> u1_mbaff),
                                         ps_cur_mb_info);

                    u4_pred_info_pkd_idx++;
                    pred_cnt++;

                }
            }
            ps_dec->p_mc_dec_thread(ps_dec, ps_cur_mb_info);
        }
        else if(ps_cur_mb_info->u1_mb_type == MB_SKIP)
        {
            {
                WORD32 pred_cnt = 0;
                pred_info_pkd_t *ps_pred_pkd;
                UWORD32 u4_pred_info_pkd_idx;
                WORD8 i1_pred;

                u4_pred_info_pkd_idx = ps_cur_mb_info->u4_pred_info_pkd_idx;



                while(pred_cnt < ps_cur_mb_info->u1_num_pred_parts)
                {

                    ps_pred_pkd = ps_dec->ps_pred_pkd + u4_pred_info_pkd_idx;


                    ps_dec->p_form_mb_part_info_thread(ps_pred_pkd,ps_dec,
                                               ps_cur_mb_info->u2_mbx,ps_cur_mb_info->u2_mby,(i >> u1_mbaff),
                                         ps_cur_mb_info);


                    u4_pred_info_pkd_idx++;
                    pred_cnt++;
                }
            }
            /* Decode MB skip */
            ps_dec->p_mc_dec_thread(ps_dec, ps_cur_mb_info);
        }

        u2_cur_dec_mb_num++;
    }

    /* N Mb IQ IT RECON  Loop */
    for(j = 0; j < i; j++)
     {
         DATA_SYNC();


         ps_cur_mb_info = &ps_dec->ps_frm_mb_info[ps_dec->cur_dec_mb_num
                         & PD_MB_BUF_SIZE_MOD];


         if(ps_cur_mb_info->u1_mb_type <= u1_skip_th)
         {
             ih264d_process_inter_mb(ps_dec, ps_cur_mb_info, j);
         }
         else if(ps_cur_mb_info->u1_mb_type != MB_SKIP)
         {
             if((u1_ipcm_th + 25) != ps_cur_mb_info->u1_mb_type)
             {
                 ps_cur_mb_info->u1_mb_type -= (u1_skip_th + 1);
                 ret = ih264d_process_intra_mb(ps_dec, ps_cur_mb_info, j);
                 if(ret != OK)
                     return ret;
             }
         }

         if(ps_dec->u4_mb_level_deblk == 1)
         {

             ih264d_deblock_mb_level(ps_dec, ps_cur_mb_info, j);
         }

         if((ps_dec->u4_num_cores >= 3) && (u1_mbaff == 0))
             ih264d_copy_intra_pred_line(ps_dec, ps_cur_mb_info, j);
         if(u1_mbaff)
         {
             if(u4_update_mbaff)
             {
                 UWORD32 u4_mb_num = ps_cur_mb_info->u2_mbx
                                 + ps_dec->u2_frm_wd_in_mbs
                                                 * (ps_cur_mb_info->u2_mby >> 1);
                 UPDATE_MB_MAP_MBNUM_BYTE(ps_dec->pu1_recon_mb_map, u4_mb_num);
                 u4_update_mbaff = 0;
             }
             else
             {
                 u4_update_mbaff = 1;
             }
         }
         else
         {
             UWORD32 u4_mb_num = ps_cur_mb_info->u2_mbx
                             + ps_dec->u2_frm_wd_in_mbs * ps_cur_mb_info->u2_mby;
             UPDATE_MB_MAP_MBNUM_BYTE(ps_dec->pu1_recon_mb_map, u4_mb_num);
         }
         ps_dec->cur_dec_mb_num++;
     }


    /*handle the last mb in picture case*/
    if(ps_dec->cur_dec_mb_num > ps_dec->ps_cur_sps->u2_max_mb_addr)
        ps_dec->u4_cur_slice_decode_done = 1;

    if(i != u1_num_mbs)
    {
        u1_end_of_row = 0;
        /*Number of MB's left in row*/
        u1_num_mbs_next = u1_num_mbs_next + ((u1_num_mbs - i) >> u1_mbaff);
    }

    ih264d_decode_tfr_nmb(ps_dec, (i), u1_num_mbs_next, u1_end_of_row);

    return OK;
}

WORD32 ih264d_decode_slice_thread(dec_struct_t *ps_dec /* Decoder parameters */
)
{
    UWORD8 u1_num_mbs_next, u1_num_mbsleft, u1_end_of_row = 0; //, u1_slice_end, u1_tfr_n_mb, u1_decode_nmb;
    const UWORD32 i2_pic_wdin_mbs = ps_dec->u2_frm_wd_in_mbs;
    UWORD8 u1_mbaff, u1_num_mbs; //,uc_more_data_flag,u1_mb_idx;

    UWORD16 u2_first_mb_in_slice;

    /*dec_bit_stream_t  *const  ps_bitstrm = ps_dec->ps_bitstrm;
     UWORD32 * pu4_bitstrm_buf = ps_bitstrm->pu4_buffer;
     UWORD32 *pu4_bitstrm_ofst  = &ps_bitstrm->u4_ofst;*/

    UWORD16 i16_mb_x, i16_mb_y;
    UWORD8 u1_field_pic;
    UWORD32 u4_frame_stride, x_offset, y_offset;
    WORD32 ret;

    tfr_ctxt_t *ps_trns_addr;

    if(ps_dec->ps_decode_cur_slice->slice_header_done != 2)
        return ERROR_INV_SLICE_HDR_T;



    u1_mbaff = ps_dec->ps_cur_slice->u1_mbaff_frame_flag;

    u2_first_mb_in_slice = ps_dec->ps_decode_cur_slice->u4_first_mb_in_slice;

    i16_mb_x = MOD(u2_first_mb_in_slice, i2_pic_wdin_mbs);
    i16_mb_y = DIV(u2_first_mb_in_slice, i2_pic_wdin_mbs);
    i16_mb_y <<= u1_mbaff;
    ps_dec->i2_dec_thread_mb_y = i16_mb_y;

    /*if((i16_mb_x > (i2_pic_wdin_mbs - 1))
                    || (i16_mb_y > ps_dec->u2_frm_ht_in_mbs - 1))
    {
    }*/
    if(ps_dec->cur_dec_mb_num == u2_first_mb_in_slice << u1_mbaff)
    {
        ps_dec->u2_mb_skip_error = 0;
    }
    else
    {
        ps_dec->u2_mb_skip_error = 1;
    }
    ps_dec->cur_dec_mb_num = u2_first_mb_in_slice << u1_mbaff;

    // recalculate recon pointers
    u1_field_pic = ps_dec->ps_cur_slice->u1_field_pic_flag;
    u4_frame_stride = ps_dec->u2_frm_wd_y << u1_field_pic;
    x_offset = i16_mb_x << 4;
    y_offset = (i16_mb_y * u4_frame_stride) << 4;

    ps_trns_addr = &(ps_dec->s_tran_addrecon);

    ps_trns_addr->pu1_dest_y = ps_dec->s_cur_pic.pu1_buf1 + x_offset + y_offset;

    u4_frame_stride = ps_dec->u2_frm_wd_uv << u1_field_pic;
    x_offset >>= 1;
    y_offset = (i16_mb_y * u4_frame_stride) << 3;

    x_offset *= YUV420SP_FACTOR;

    ps_trns_addr->pu1_dest_u = ps_dec->s_cur_pic.pu1_buf2 + x_offset + y_offset;
    ps_trns_addr->pu1_dest_v = ps_dec->s_cur_pic.pu1_buf3 + x_offset + y_offset;

    ps_trns_addr->pu1_mb_y = ps_trns_addr->pu1_dest_y;
    ps_trns_addr->pu1_mb_u = ps_trns_addr->pu1_dest_u;
    ps_trns_addr->pu1_mb_v = ps_trns_addr->pu1_dest_v;

    if(ps_dec->u4_mb_level_deblk == 1)
    {
        /*If it is not the first mb in row,the previous MB which needs to be deblocked
         * as there is delay of 1 MB*/
        if(i16_mb_x != 0)
        {
            ps_trns_addr->pu1_mb_y -= MB_SIZE;
            ps_trns_addr->pu1_mb_u -= BLK8x8SIZE * YUV420SP_FACTOR;
            ps_trns_addr->pu1_mb_v -= BLK8x8SIZE;
        }
    }

    /**********Number of Mbs in Slice**********/

    ps_dec->ps_deblk_mbn_dec_thrd = ps_dec->ps_deblk_pic
                    + (u2_first_mb_in_slice << u1_mbaff);

    /* Initialise MC and formMbPartInfo fn ptrs one time based on profile_idc */

    {
        ps_dec->p_mc_dec_thread = ih264d_motion_compensate_bp;
        ps_dec->p_form_mb_part_info_thread = ih264d_form_mb_part_info_bp;
    }
    {
        UWORD8 uc_nofield_nombaff;
        uc_nofield_nombaff = ((ps_dec->ps_cur_slice->u1_field_pic_flag == 0)
                        && (ps_dec->ps_cur_slice->u1_mbaff_frame_flag == 0)
                        && (ps_dec->ps_decode_cur_slice->slice_type != B_SLICE)
                        && (ps_dec->ps_cur_pps->u1_wted_pred_flag == 0));

        if(uc_nofield_nombaff == 0)
        {
            ps_dec->p_mc_dec_thread = ih264d_motion_compensate_mp;
            ps_dec->p_form_mb_part_info_thread = ih264d_form_mb_part_info_mp;
        }

    }

    ps_dec->u4_cur_slice_decode_done = 0;


    while(ps_dec->u4_cur_slice_decode_done != 1)
    {

        u1_num_mbsleft = ((i2_pic_wdin_mbs - i16_mb_x) << u1_mbaff);

        if(u1_num_mbsleft <= ps_dec->u1_recon_mb_grp)
        {
            u1_num_mbs = u1_num_mbsleft;

            /*Indicate number of mb's left in a row*/
            u1_num_mbs_next = 0;
            u1_end_of_row = 1;
            i16_mb_x = 0;
        }
        else
        {
            u1_num_mbs = ps_dec->u1_recon_mb_grp;

            /*Indicate number of mb's left in a row*/
            u1_num_mbs_next = i2_pic_wdin_mbs - i16_mb_x
                            - (ps_dec->u1_recon_mb_grp >> u1_mbaff);
            i16_mb_x += (u1_num_mbs >> u1_mbaff);
            u1_end_of_row = 0;

        }
        ret = ih264d_decode_recon_tfr_nmb_thread(ps_dec, u1_num_mbs, u1_num_mbs_next,
                                           u1_end_of_row);
        if(ret != OK)
            return ret;
    }
    return OK;
}

void ih264d_decode_picture_thread(dec_struct_t *ps_dec )
{

    ithread_set_name("ih264d_decode_picture_thread");

    // run the loop till all slices are decoded

    while(1)
    {
        if(ps_dec->u4_start_frame_decode)
        {
            break;
        }
        else
        {
            NOP(32);

        }
    }

    DEBUG_THREADS_PRINTF("Got start of frame u4_flag\n");

    if(ps_dec->u4_start_frame_decode == 1)
    {
        while(1)
        {
            /*Complete all writes before processing next slice*/
            DATA_SYNC();
            /*wait untill all the slice params have been populated*/
            while(ps_dec->ps_decode_cur_slice->slice_header_done == 0)
            {
                NOP(32); DEBUG_THREADS_PRINTF(" waiting for slice header \n");
            }

            DEBUG_THREADS_PRINTF(" Entering decode slice\n");

            ih264d_decode_slice_thread(ps_dec);
            DEBUG_THREADS_PRINTF(" Exit  ih264d_decode_slice_thread \n");

            /*Complete all writes before processing next slice*/
            DATA_SYNC();

            while(1)
            {
                volatile void * parse_addr, *dec_addr;
                volatile UWORD32 last_slice;

                parse_addr = (volatile void *)ps_dec->ps_parse_cur_slice;
                dec_addr = (volatile void *)ps_dec->ps_decode_cur_slice;
                last_slice = ps_dec->ps_decode_cur_slice->last_slice_in_frame;

                if(last_slice == 1)
                    break;

                if(parse_addr != dec_addr)
                    break;

                DEBUG_THREADS_PRINTF("Waiting for next slice or end of frame\n");

                NOP(32);
            }

            DEBUG_THREADS_PRINTF("Got next slice/end of frame signal \n ");

            if((void *)ps_dec->ps_parse_cur_slice
                            > (void *)ps_dec->ps_decode_cur_slice)
            {
                ps_dec->ps_decode_cur_slice++;
                ps_dec->u2_cur_slice_num_dec_thread++;
            }
            else
            {
                /*Last slice in frame*/
                break;
            }

        }
    }

    if(ps_dec->u4_output_present)
    {
        while(1)
        {
            volatile UWORD32 *u4_flag = &(ps_dec->as_fmt_conv_part[1].u4_flag);

            DEBUG_THREADS_PRINTF(" Format conversion loop in decode *u4_flag = %d\n",*u4_flag);
            if(2 == *u4_flag)
            {
                if(ps_dec->as_fmt_conv_part[1].u4_num_rows_y)
                    ih264d_format_convert(
                                    ps_dec, &(ps_dec->s_disp_op),
                                    ps_dec->as_fmt_conv_part[1].u4_start_y,
                                    ps_dec->as_fmt_conv_part[1].u4_num_rows_y);

                break;
            }
            else if(1 == *u4_flag)
            {
                NOP(32);

            }
            else
                break;

        }
    }

    ithread_exit(0);

}

void ih264d_signal_decode_thread(dec_struct_t *ps_dec)
{
    if(ps_dec->u4_dec_thread_created == 1)
    {

        if(ps_dec->u4_start_frame_decode == 1)
            ps_dec->ps_parse_cur_slice->last_slice_in_frame = 1;
        else
            /*to indicate frame in error*/
            ps_dec->u4_start_frame_decode = 2;

        ithread_join(ps_dec->pv_dec_thread_handle, NULL);
        ps_dec->u4_dec_thread_created = 0;
    }
}
void ih264d_signal_bs_deblk_thread(dec_struct_t *ps_dec)
{
    if(ps_dec->u4_bs_deblk_thread_created)
    {
        /*signal error*/
        if(ps_dec->u4_start_bs_deblk == 0)
            ps_dec->u4_start_bs_deblk = 2;

        ithread_join(ps_dec->pv_bs_deblk_thread_handle, NULL);
        ps_dec->u4_bs_deblk_thread_created = 0;
    }

}
