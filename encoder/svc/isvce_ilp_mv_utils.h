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
*  isvce_ilp_mv_utils.h
*
* @brief
*  Defs to perform experiments in ilp mv
*
*
* @remarks
*  None
*
*******************************************************************************
*/
#ifndef _ISVCE_ILP_MV_UTILS_H_
#define _ISVCE_ILP_MV_UTILS_H_

#include <stdbool.h>

#include "ih264_typedefs.h"
#include "isvc_defs.h"
#include "isvc_macros.h"
#include "isvce_pred_structs.h"
#include "isvce_structs.h"

#define MAX_CAND_IF_NUM_ILP_MV_LT_2 8
#define MAX_CAND_IF_NUM_ILP_MV_GTEQ_2 6

/* nbr_mb.x, nbr_mb.y, pu_pos.x, pu_pos.y */
#define NBR_PU_AND_MB_POS 4

static const WORD8 gai1_nbr_ilp_mv_map[MAX_ILP_MV_IN_NBR_RGN][NBR_PU_AND_MB_POS] = {
    {-1, 0, 3, 0},
    {0, -1, 0, 3},
    {1, 0, 0, 0},
    {0, 1, 0, 0},
};

/**
*******************************************************************************
*
* @brief
*  This function checks if the max difference between ILP MVs is less than four
* or not if number of ILP MVs is greater than or equal to two
*
* @param[in] ps_me
*  Pointer to ilp_me_cands
*
* @returns  One if number of ILP MVs is greater than equal to two and max
* difference between them is less than 4 otherwise returns zero
*
* @remarks none
*
*******************************************************************************
*/
static FORCEINLINE bool isvce_check_max_mv_diff_lt_4(ilp_me_cands_t *ps_ilp_me_cands,
                                                     WORD32 i4_reflist)
{
    UWORD32 i, j;
    UWORD32 u4_mv_diff_x, u4_mv_diff_y;

    for(i = 1; i < ps_ilp_me_cands->u4_num_ilp_mvs; i++)
    {
        for(j = 0; j < i; j++)
        {
            if(((ps_ilp_me_cands->ae_pred_mode[i] == ((PRED_MODE_T) i4_reflist)) ||
                ((ps_ilp_me_cands->ae_pred_mode[i] == BI))) &&
               ((ps_ilp_me_cands->ae_pred_mode[j] == ((PRED_MODE_T) i4_reflist)) ||
                ((ps_ilp_me_cands->ae_pred_mode[j] == BI))))
            {
                u4_mv_diff_x = ABS(ps_ilp_me_cands->as_mv[i][i4_reflist].s_mv.i2_mvx -
                                   ps_ilp_me_cands->as_mv[j][i4_reflist].s_mv.i2_mvx);

                u4_mv_diff_y = ABS(ps_ilp_me_cands->as_mv[i][i4_reflist].s_mv.i2_mvy -
                                   ps_ilp_me_cands->as_mv[j][i4_reflist].s_mv.i2_mvy);

                if(u4_mv_diff_x >= 4 || u4_mv_diff_y >= 4)
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }

    return true;
}

#endif
