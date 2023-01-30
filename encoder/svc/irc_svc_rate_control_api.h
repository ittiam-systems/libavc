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

#ifndef _IRC_SVC_RATE_CONTROL_API_H_
#define _IRC_SVC_RATE_CONTROL_API_H_

/* Dependencies of 'irc_rate_control_api_structs' */
#include "irc_picture_type.h"
#include "irc_rd_model.h"
#include "irc_vbr_storage_vbv.h"
#include "irc_est_sad.h"
#include "irc_bit_allocation.h"
#include "irc_mb_model_based.h"
#include "irc_cbr_buffer_control.h"
#include "irc_vbr_str_prms.h"
#include "irc_common.h"

#include "irc_rate_control_api_structs.h"

/* Get frame level QP based on BPP and GPP */
UWORD8 irc_get_frame_level_init_qp(rate_control_api_t *ps_rate_control_api, rc_type_e e_rc_type,
                                   picture_type_e e_pic_type, DOUBLE d_bpp, DOUBLE d_gpp);

void irc_change_qp_constraints(rate_control_api_t *ps_rate_control_api, UWORD8 *pu1_min_max_qp,
                               UWORD8 *pu1_min_max_avc_qp);

extern UWORD8 irc_is_scenecut(rate_control_api_t *ps_rate_control_api);

#endif
