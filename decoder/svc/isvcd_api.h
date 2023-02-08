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
 *  isvcd_api.h
 *
 * @brief
 *  This file contains all the necessary structure and enumeration definitions
 *    needed for the Application
 *  Program Interface(API) of the Ittiam H264 svc decoder
 *
 * @author
 *  Kishore
 *
 * @remarks
 *  None
 *
 *******************************************************************************
 */

#ifndef _ISVCD_API_H_
#define _ISVCD_API_H_

#include "iv.h"
#include "ivd_ext.h"

extern IV_API_CALL_STATUS_T isvcd_create(ivd_avc_ext_create_ip_t *ps_ip,
                                         ivd_avc_ext_create_op_t *ps_op);

extern IV_API_CALL_STATUS_T isvcd_delete(iv_obj_t *ps_dec_hdl);

extern IV_API_CALL_STATUS_T isvcd_decode_au(iv_obj_t *ps_dec_hdl,
                                            ivd_avc_ext_video_decode_ip_t *ps_ip,
                                            ivd_avc_ext_video_decode_op_t *ps_op);

extern IV_API_CALL_STATUS_T isvcd_ctl_cmd_handler(iv_obj_t *ps_dec_hdl, ivd_avc_ext_ctl_ip_t *ps_ip,
                                                  ivd_avc_ext_ctl_op_t *ps_op);

#endif /*_ISVCD_API_H_*/