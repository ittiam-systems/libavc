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
 *  isvcd_cabac.c
 *
 * @brief
 *  This file contains Binary decoding routines.
 *
 * @author
 *  Kishore
 *
 * @par List of Functions:
 *  - isvcd_init_cabac_contexts()
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#include <string.h>
#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "isvcd_structs.h"
#include "ih264d_cabac.h"
#include "isvcd_cabac.h"
#include "ih264d_bitstrm.h"
#include "ih264d_error_handler.h"
#include "ih264d_defs.h"
#include "ih264d_debug.h"
#include "ih264d_tables.h"
#include "isvcd_tables.h"
#include "ih264d_parse_cabac.h"
#include "ih264d_tables.h"

/*****************************************************************************/
/*                                                                           */
/*  Function Name : isvcd_init_cabac_contexts                                */
/*                                                                           */
/*  Description   : This function initializes the cabac contexts             */
/*                  depending upon slice type and Init_Idc value.            */
/*  Inputs        : ps_dec, slice type                                       */
/*  Globals       : <Does it use any global variables?>                      */
/*  Outputs       :                                                          */
/*  Returns       : void                                                     */
/*                                                                           */
/*  Issues        : none                                                     */
/*                                                                           */
/*  Revision History:                                                        */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         06 09 2021   Kishore         Draft                                */
/*                                                                           */
/*****************************************************************************/

void isvcd_init_cabac_contexts(UWORD8 u1_slice_type, dec_struct_t *ps_dec)
{
    bin_ctxt_model_t *p_cabac_ctxt_table_t = ps_dec->p_cabac_ctxt_table_t;
    UWORD8 u1_qp_y = ps_dec->ps_cur_slice->u1_slice_qp;
    UWORD8 u1_cabac_init_Idc = 0;

    if(I_SLICE != u1_slice_type)
    {
        u1_cabac_init_Idc = ps_dec->ps_cur_slice->u1_cabac_init_idc;
    }

    {
        /* MAKING ps_dec->p_ctxt_inc_mb_map a scratch buffer */
        /* 0th entry of CtxtIncMbMap will be always be containing default values
         for CABAC context representing MB not available */
        ctxt_inc_mb_info_t *p_DefCtxt = ps_dec->p_ctxt_inc_mb_map - 1;
        UWORD8 *pu1_temp;
        WORD8 i;
        p_DefCtxt->u1_mb_type = CAB_SKIP;
        p_DefCtxt->u1_cbp = 0x0f;
        p_DefCtxt->u1_intra_chroma_pred_mode = 0;
        p_DefCtxt->u1_yuv_dc_csbp = 0x7;
        p_DefCtxt->u1_transform8x8_ctxt = 0;

        pu1_temp = (UWORD8 *) p_DefCtxt->i1_ref_idx;
        for(i = 0; i < 4; i++, pu1_temp++) (*pu1_temp) = 0;
        pu1_temp = (UWORD8 *) p_DefCtxt->u1_mv;
        for(i = 0; i < 16; i++, pu1_temp++) (*pu1_temp) = 0;
        ps_dec->ps_def_ctxt_mb_info = p_DefCtxt;
    }

    if(u1_slice_type == I_SLICE)
    {
        u1_cabac_init_Idc = 3;
        ps_dec->p_mb_type_t = p_cabac_ctxt_table_t + MB_TYPE_I_SLICE;
    }
    else if(u1_slice_type == P_SLICE)
    {
        ps_dec->p_mb_type_t = p_cabac_ctxt_table_t + MB_TYPE_P_SLICE;
        ps_dec->p_mb_skip_flag_t = p_cabac_ctxt_table_t + MB_SKIP_FLAG_P_SLICE;
        ps_dec->p_sub_mb_type_t = p_cabac_ctxt_table_t + SUB_MB_TYPE_P_SLICE;
    }
    else if(u1_slice_type == B_SLICE)
    {
        ps_dec->p_mb_type_t = p_cabac_ctxt_table_t + MB_TYPE_B_SLICE;
        ps_dec->p_mb_skip_flag_t = p_cabac_ctxt_table_t + MB_SKIP_FLAG_B_SLICE;
        ps_dec->p_sub_mb_type_t = p_cabac_ctxt_table_t + SUB_MB_TYPE_B_SLICE;
    }
    {
        bin_ctxt_model_t *p_cabac_ctxt_table_t_tmp = p_cabac_ctxt_table_t;
        if(ps_dec->ps_cur_slice->u1_field_pic_flag)
        {
            p_cabac_ctxt_table_t_tmp += SIGNIFICANT_COEFF_FLAG_FLD;
        }
        else
        {
            p_cabac_ctxt_table_t_tmp += SIGNIFICANT_COEFF_FLAG_FRAME;
        }
        {
            bin_ctxt_model_t **p_significant_coeff_flag_t = ps_dec->p_significant_coeff_flag_t;
            p_significant_coeff_flag_t[0] = p_cabac_ctxt_table_t_tmp + SIG_COEFF_CTXT_CAT_0_OFFSET;
            p_significant_coeff_flag_t[1] = p_cabac_ctxt_table_t_tmp + SIG_COEFF_CTXT_CAT_1_OFFSET;
            p_significant_coeff_flag_t[2] = p_cabac_ctxt_table_t_tmp + SIG_COEFF_CTXT_CAT_2_OFFSET;
            p_significant_coeff_flag_t[3] = p_cabac_ctxt_table_t_tmp + SIG_COEFF_CTXT_CAT_3_OFFSET;
            p_significant_coeff_flag_t[4] = p_cabac_ctxt_table_t_tmp + SIG_COEFF_CTXT_CAT_4_OFFSET;

            p_significant_coeff_flag_t[5] = p_cabac_ctxt_table_t_tmp + SIG_COEFF_CTXT_CAT_5_OFFSET;
        }
    }

    memcpy(p_cabac_ctxt_table_t, gau1_isvcd_cabac_ctxt_init_table[u1_cabac_init_Idc][u1_qp_y],
           NUM_CABAC_CTXTS_SVC * sizeof(bin_ctxt_model_t));
}
