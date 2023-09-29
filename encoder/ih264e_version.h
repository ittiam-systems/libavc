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

/**
*******************************************************************************
* @file
*  ih264e_version.h
*
* @brief
*  Contains declarations of miscellaneous utility functions used by the encoder
*
* @author
*  ittiam
*
* @remarks
*  none
*
*******************************************************************************
*/

#ifndef _IH264E_VERSION_H_
#define _IH264E_VERSION_H_


/*****************************************************************************/
/* Function Declarations                                                     */
/*****************************************************************************/

IV_STATUS_T ih264e_get_version(CHAR *pc_version, UWORD32 u4_version_bufsize);

#endif /* _IH264E_VERSION_H_ */
