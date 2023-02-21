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

/*****************************************************************************/
/*                                                                           */
/*  File Name         : ih264d_sei.h                                                */
/*                                                                           */
/*  Description       : This file contains routines to parse SEI NAL's       */
/*                                                                           */
/*  List of Functions : <List the functions defined in this file>            */
/*                                                                           */
/*  Issues / Problems : None                                                 */
/*                                                                           */
/*  Revision History  :                                                      */
/*                                                                           */
/*         DD MM YYYY   Author(s)       Changes (Describe the changes made)  */
/*         25 05 2005   NS              Draft                                */
/*                                                                           */
/*****************************************************************************/

#ifndef _IH264D_SEI_H_
#define _IH264D_SEI_H_

#include "ih264_typedefs.h"
#include "ih264_macros.h"
#include "ih264_platform_macros.h"
#include "ih264d_bitstrm.h"
#include "ih264d_structs.h"
#include "ih264d.h"

#define SEI_BUF_PERIOD      0
#define SEI_PIC_TIMING      1
#define SEI_PAN_SCAN_RECT   2
#define SEI_FILLER          3
#define SEI_UD_REG_T35      4
#define SEI_UD_UN_REG       5
#define SEI_RECOVERY_PT     6
#define SEI_DEC_REF_MARK    7
#define SEI_SPARE_PIC       8
#define SEI_SCENE_INFO      9
#define SEI_SUB_SEQN_INFO   10
#define SEI_SUB_SEQN_LAY_CHAR       11
#define SEI_SUB_SEQN_CHAR   12
#define SEI_FULL_FRAME_FREEZE       13
#define SEI_FULL_FRAME_FREEZE_REL   14
#define SEI_FULL_FRAME_SNAP_SHOT    15
#define SEI_PROG_REF_SEGMENT_START  16
#define SEI_PROG_REF_SEGMENT_END    17
#define SEI_MOT_CON_SLICE_GRP_SET   18
#define SEI_FILM_GRAIN_CHARACTERISTICS 19
#define SEI_MASTERING_DISP_COL_VOL       137
#define SEI_CONTENT_LIGHT_LEVEL_DATA     144
#define SEI_AMBIENT_VIEWING_ENVIRONMENT  148
#define SEI_CONTENT_COLOR_VOLUME         149
#define SEI_SHUTTER_INTERVAL_INFO        205

/* Declaration of dec_struct_t to avoid CCS compilation Error */
struct _DecStruct;
WORD32 ih264d_parse_sei_message(struct _DecStruct *ps_dec,
                                dec_bit_stream_t *ps_bitstrm);
typedef struct
{
    UWORD8 u1_seq_parameter_set_id;
    UWORD32 u4_initial_cpb_removal_delay;
    UWORD32 u4_nitial_cpb_removal_delay_offset;

} buf_period_t;

/**
 * Structure to hold Mastering Display Color Volume SEI
 */
typedef struct
{
    /**
     * Array to store the display_primaries_x values
     */
    UWORD16 au2_display_primaries_x[NUM_SEI_MDCV_PRIMARIES];

    /**
     * Array to store the display_primaries_y values
     */
    UWORD16 au2_display_primaries_y[NUM_SEI_MDCV_PRIMARIES];

    /**
     * Variable to store the white point x value
     */
    UWORD16 u2_white_point_x;

    /**
     * Variable to store the white point y value
     */
    UWORD16 u2_white_point_y;

    /**
     * Variable to store the max display mastering luminance value
     */
    UWORD32 u4_max_display_mastering_luminance;

    /**
     * Variable to store the min display mastering luminance value
     */
    UWORD32 u4_min_display_mastering_luminance;

}sei_mdcv_params_t;


/**
 * Structure for Content Light Level Info
 *
 */
typedef struct
{
    /**
     * The maximum pixel intensity of all samples
     */
    UWORD16 u2_max_content_light_level;

    /**
     * The average pixel intensity of all samples
     */
    UWORD16 u2_max_pic_average_light_level;

}sei_cll_params_t;


/**
 * Structure to hold Ambient viewing environment SEI
 */
typedef struct
{

    /**
     * specifies the environmental illuminance of the ambient viewing environment
     */
    UWORD32 u4_ambient_illuminance;

    /*
     * specify the normalized x chromaticity coordinates of the
     * environmental ambient light in the nominal viewing environment
     */
    UWORD16 u2_ambient_light_x;

    /*
    * specify the normalized y chromaticity coordinates of the
    * environmental ambient light in the nominal viewing environment
    */
    UWORD16 u2_ambient_light_y;

}sei_ave_params_t;


/**
 * Structure to hold Content color volume SEI
 */
typedef struct
{
    /*
     * Flag used to control persistence of CCV SEI messages
     */
    UWORD8 u1_ccv_cancel_flag;

    /*
     * specifies the persistence of the CCV SEI message for the current layer
     */
    UWORD8 u1_ccv_persistence_flag;

    /*
     * specifies the presence of syntax elements ccv_primaries_x and ccv_primaries_y
     */
    UWORD8 u1_ccv_primaries_present_flag;

    /*
     * specifies that the syntax element ccv_min_luminance_value is present
     */
    UWORD8 u1_ccv_min_luminance_value_present_flag;

    /*
     * specifies that the syntax element ccv_max_luminance_value is present
     */
    UWORD8 u1_ccv_max_luminance_value_present_flag;

    /*
     * specifies that the syntax element ccv_avg_luminance_value is present
     */
    UWORD8 u1_ccv_avg_luminance_value_present_flag;

    /*
     * shall be equal to 0 in bitstreams conforming to this version. Other values
     * for reserved_zero_2bits are reserved for future use
     */
    UWORD8 u1_ccv_reserved_zero_2bits;

    /*
     * specify the normalized x chromaticity coordinates of the colour
     * primary component c of the nominal content colour volume
     */
    WORD32 ai4_ccv_primaries_x[NUM_SEI_CCV_PRIMARIES];

    /*
     * specify the normalized y chromaticity coordinates of the colour
     * primary component c of the nominal content colour volume
     */
    WORD32 ai4_ccv_primaries_y[NUM_SEI_CCV_PRIMARIES];

    /*
     * specifies the normalized minimum luminance value
     */
    UWORD32 u4_ccv_min_luminance_value;

    /*
     * specifies the normalized maximum luminance value
     */
    UWORD32 u4_ccv_max_luminance_value;

    /*
     * specifies the normalized average luminance value
     */
    UWORD32 u4_ccv_avg_luminance_value;

}sei_ccv_params_t;

/**
 * Structure to hold Shutter Interval Info SEI
 */
typedef struct
{
    /**
     * specifies if the sei sii is enabled
     */
    UWORD8 u1_shutter_interval_info_present_flag;

    /**
     * specifies the shutter interval temporal sub-layer index
     * of the current picture
     */
    UWORD32 u4_sii_sub_layer_idx;

    /**
     * specify the number of time units that pass in one second
     */
    UWORD32 u4_sii_time_scale;

    /**
     * specifies that the indicated shutter interval is the same for all
     * pictures in the coded video sequence
     */
    UWORD8 u1_fixed_shutter_interval_within_cvs_flag;

    /**
     * specifies the the number of time units of a clock operating at the
     * frequency sii_time_scale Hz that corresponds to the indicated shutter
     * interval of each picture in the coded video sequence
     */
    UWORD32 u4_sii_num_units_in_shutter_interval;

    /**
     * sii_max_sub_layers_minus1 plus 1 specifies the maximum number of
     * shutter interval temporal sub-layers indexes that may be present
     * in the coded video sequence
     */
    UWORD8 u1_sii_max_sub_layers_minus1;

    /**
     * specifies the number of time units of a clock operating at the
     * frequency sii_time_scale Hz that corresponds to the shutter
     * interval of each picture in the coded video sequence
     */
    UWORD32 au4_sub_layer_num_units_in_shutter_interval[SII_MAX_SUB_LAYERS];

} sei_sii_params_t;

typedef struct
{
    /**
     * Flag to control the presence of FGC SEI params
     */
    UWORD8 u1_film_grain_characteristics_cancel_flag;

    /**
     * Specifies the pic order count
     */
    WORD32 i4_poc;

    /**
     * Specifies IDR pic ID
     */
    UWORD32 u4_idr_pic_id;

    /**
     * Specifies film grain model for simulation
     */
    UWORD8 u1_film_grain_model_id;

    /**
     * Specifies separate color format for decoded samples and grain
     */
    UWORD8 u1_separate_colour_description_present_flag;

    /**
     * Specifies the bit depth used for the luma component
     */
    UWORD8 u1_film_grain_bit_depth_luma_minus8;

    /**
     * Specifies the bit depth used for the Cb and Cr components
     */
    UWORD8 u1_film_grain_bit_depth_chroma_minus8;

    /**
     * Specifies the colour space of the FGC in SEI
     */
    UWORD8 u1_film_grain_full_range_flag;

    /**
     * Specifies the colour space of the FGC in SEI
     */
    UWORD8 u1_film_grain_colour_primaries;

    /**
     * Specifies the colour space of the FGC in SEI
     */
    UWORD8 u1_film_grain_transfer_characteristics;

    /**
     * Specifies the colour space of the FGC in SEI
     */
    UWORD8 u1_film_grain_matrix_coefficients;

    /**
     * identifies the blending mode used to blend the simulated film grain with the decoded images
     */
    UWORD8 u1_blending_mode_id;

    /**
     * Specifies a scale factor used in the film grain characterization equations
     */
    UWORD8 u1_log2_scale_factor;

    /**
     * Indicates whether film grain is modelled or not on the colour component
     */
    UWORD8 au1_comp_model_present_flag[SEI_FGC_NUM_COLOUR_COMPONENTS];

    /**
     * Specifies the number of intensity intervals for which
     * a specific set of model values has been estimated
     */
    UWORD8 au1_num_intensity_intervals_minus1[SEI_FGC_NUM_COLOUR_COMPONENTS];

    /**
     * Specifies the number of model values present for each intensity interval in which
     * the film grain has been modelled
     */
    UWORD8 au1_num_model_values_minus1[SEI_FGC_NUM_COLOUR_COMPONENTS];

    /**
     * Specifies the lower bound of the interval of intensity levels for which
     * the set of model values applies
     */
    UWORD8 au1_intensity_interval_lower_bound[SEI_FGC_NUM_COLOUR_COMPONENTS]
                                             [SEI_FGC_MAX_NUM_INTENSITY_INTERVALS];

    /**
     * Specifies the upper bound of the interval of intensity levels for which
     * the set of model values applies
     */
    UWORD8 au1_intensity_interval_upper_bound[SEI_FGC_NUM_COLOUR_COMPONENTS]
                                             [SEI_FGC_MAX_NUM_INTENSITY_INTERVALS];

    /**
     * Represents each one of the model values present for
     * the colour component and intensity interval
     */
    WORD32 ai4_comp_model_value[SEI_FGC_NUM_COLOUR_COMPONENTS][SEI_FGC_MAX_NUM_INTENSITY_INTERVALS]
                               [SEI_FGC_MAX_NUM_MODEL_VALUES];

    /**
     * Specifies the persistence of the film grain characteristics SEI message
     */
    UWORD32 u4_film_grain_characteristics_repetition_period;

} sei_fgc_params_t;

struct _sei
{
    UWORD8 u1_seq_param_set_id;
    buf_period_t s_buf_period;
    UWORD8 u1_pic_struct;
    UWORD16 u2_recovery_frame_cnt;
    UWORD8 u1_exact_match_flag;
    UWORD8 u1_broken_link_flag;
    UWORD8 u1_changing_slice_grp_idc;
    UWORD8 u1_is_valid;

    /**
     *  mastering display color volume info present flag
     */
    UWORD8 u1_sei_mdcv_params_present_flag;

    /*
     * MDCV parameters
     */
    sei_mdcv_params_t s_sei_mdcv_params;

    /**
     * content light level info present flag
     */
    UWORD8 u1_sei_cll_params_present_flag;

    /*
     * CLL parameters
     */
    sei_cll_params_t s_sei_cll_params;

    /**
     * ambient viewing environment info present flag
     */
    UWORD8 u1_sei_ave_params_present_flag;

    /*
     * AVE parameters
     */
    sei_ave_params_t s_sei_ave_params;

    /**
     * content color volume info present flag
     */
    UWORD8 u1_sei_ccv_params_present_flag;

    /*
     * CCV parameters
     */
    sei_ccv_params_t s_sei_ccv_params;

    /**
     * shutter interval info present flag
     */
    UWORD8 u1_sei_sii_params_present_flag;

    /*
     * SII parameters
     */
    sei_sii_params_t s_sei_sii_params;

    /**
     * film grain params info present flag
     */
    UWORD8 u1_sei_fgc_params_present_flag;

    /*
     * film grain characteristics parameters
     */
    sei_fgc_params_t s_sei_fgc_params;
};
typedef struct _sei sei;

WORD32 ih264d_export_sei_mdcv_params(ivd_sei_decode_op_t *ps_sei_decode_op,
                                     sei *ps_sei, sei *ps_sei_export);

WORD32 ih264d_export_sei_cll_params(ivd_sei_decode_op_t *ps_sei_decode_op,
                                    sei *ps_sei, sei *ps_sei_export);

WORD32 ih264d_export_sei_ave_params(ivd_sei_decode_op_t *ps_sei_decode_op,
                                    sei *ps_sei, sei *ps_sei_export);

WORD32 ih264d_export_sei_ccv_params(ivd_sei_decode_op_t *ps_sei_decode_op,
                                    sei *ps_sei, sei *ps_sei_export);

WORD32 ih264d_export_sei_sii_params(ivd_sei_decode_op_t *ps_sei_decode_op, sei *ps_sei,
                                    sei *ps_sei_export);

WORD32 ih264d_export_sei_fgc_params(ivd_sei_decode_op_t *ps_sei_decode_op, sei *ps_sei,
                                    sei *ps_sei_export);

#endif /* _IH264D_SEI_H_ */

