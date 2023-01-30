list(
  APPEND
  LIBSVCENC_SRCS
  "${AVC_ROOT}/encoder/ih264e_bitstream.c"
  "${AVC_ROOT}/encoder/ih264e_cabac_init.c"
  "${AVC_ROOT}/encoder/ih264e_core_coding.c"
  "${AVC_ROOT}/encoder/ih264e_encode_header.c"
  "${AVC_ROOT}/encoder/ih264e_fmt_conv.c"
  "${AVC_ROOT}/encoder/ih264e_globals.c"
  "${AVC_ROOT}/encoder/ih264e_half_pel.c"
  "${AVC_ROOT}/encoder/ih264e_intra_modes_eval.c"
  "${AVC_ROOT}/encoder/ih264e_mc.c"
  "${AVC_ROOT}/encoder/ih264e_me.c"
  "${AVC_ROOT}/encoder/ih264e_modify_frm_rate.c"
  "${AVC_ROOT}/encoder/ih264e_rate_control.c"
  "${AVC_ROOT}/encoder/ih264e_rc_mem_interface.c"
  "${AVC_ROOT}/encoder/ih264e_sei.c"
  "${AVC_ROOT}/encoder/ih264e_time_stamp.c"
  "${AVC_ROOT}/encoder/ih264e_utils.c"
  "${AVC_ROOT}/encoder/ih264e_version.c"
  "${AVC_ROOT}/encoder/ime.c"
  "${AVC_ROOT}/encoder/ime_distortion_metrics.c"
  "${AVC_ROOT}/encoder/irc_bit_allocation.c"
  "${AVC_ROOT}/encoder/irc_cbr_buffer_control.c"
  "${AVC_ROOT}/encoder/irc_est_sad.c"
  "${AVC_ROOT}/encoder/irc_fixed_point_error_bits.c"
  "${AVC_ROOT}/encoder/irc_frame_info_collector.c"
  "${AVC_ROOT}/encoder/irc_mb_model_based.c"
  "${AVC_ROOT}/encoder/irc_picture_type.c"
  "${AVC_ROOT}/encoder/irc_rate_control_api.c"
  "${AVC_ROOT}/encoder/irc_rd_model.c"
  "${AVC_ROOT}/encoder/irc_vbr_storage_vbv.c"
  "${AVC_ROOT}/encoder/irc_vbr_str_prms.c"
  "${AVC_ROOT}/encoder/svc/irc_svc_rate_control_api.c"
  "${AVC_ROOT}/encoder/svc/isvce_api.c"
  "${AVC_ROOT}/encoder/svc/isvce_cabac.c"
  "${AVC_ROOT}/encoder/svc/isvce_cabac_encode.c"
  "${AVC_ROOT}/encoder/svc/isvce_cabac_init.c"
  "${AVC_ROOT}/encoder/svc/isvce_cavlc.c"
  "${AVC_ROOT}/encoder/svc/isvce_core_coding.c"
  "${AVC_ROOT}/encoder/svc/isvce_deblk.c"
  "${AVC_ROOT}/encoder/svc/isvce_downscaler.c"
  "${AVC_ROOT}/encoder/svc/isvce_encode.c"
  "${AVC_ROOT}/encoder/svc/isvce_encode_header.c"
  "${AVC_ROOT}/encoder/svc/isvce_fmt_conv.c"
  "${AVC_ROOT}/encoder/svc/isvce_function_selector_generic.c"
  "${AVC_ROOT}/encoder/svc/isvce_globals.c"
  "${AVC_ROOT}/encoder/svc/isvce_ibl_eval.c"
  "${AVC_ROOT}/encoder/svc/isvce_ilp_mv.c"
  "${AVC_ROOT}/encoder/svc/isvce_intra_modes_eval.c"
  "${AVC_ROOT}/encoder/svc/isvce_mc.c"
  "${AVC_ROOT}/encoder/svc/isvce_me.c"
  "${AVC_ROOT}/encoder/svc/isvce_mode_stat_visualiser.c"
  "${AVC_ROOT}/encoder/svc/isvce_nalu_stat_aggregator.c"
  "${AVC_ROOT}/encoder/svc/isvce_process.c"
  "${AVC_ROOT}/encoder/svc/isvce_rate_control.c"
  "${AVC_ROOT}/encoder/svc/isvce_rc_mem_interface.c"
  "${AVC_ROOT}/encoder/svc/isvce_rc_utils.c"
  "${AVC_ROOT}/encoder/svc/isvce_residual_pred.c"
  "${AVC_ROOT}/encoder/svc/isvce_sub_pic_rc.c"
  "${AVC_ROOT}/encoder/svc/isvce_utils.c"
  "${AVC_ROOT}/encoder/psnr.c")

include_directories(${AVC_ROOT}/encoder)
include_directories(${AVC_ROOT}/encoder/svc)

if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" OR
   "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch32")
  list(
    APPEND
    LIBSVCENC_ASMS
    "${AVC_ROOT}/encoder/arm/svc/isvce_function_selector.c"
    "${AVC_ROOT}/encoder/arm/svc/isvce_function_selector_a9q.c"
    "${AVC_ROOT}/encoder/arm/svc/isvce_function_selector_av8.c"
    "${AVC_ROOT}/encoder/arm/svc/isvce_downscaler_neon.c"
    "${AVC_ROOT}/encoder/arm/svc/isvce_rc_utils_neon.c"
    "${AVC_ROOT}/encoder/arm/svc/isvce_residual_pred_neon.c")

  include_directories(${AVC_ROOT}/encoder/arm)
  include_directories(${AVC_ROOT}/encoder/arm/svc)
endif()

if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "aarch64")
  list(
    APPEND
    LIBSVCENC_ASMS
    "${AVC_ROOT}/encoder/armv8/ih264e_evaluate_intra16x16_modes_av8.s"
    "${AVC_ROOT}/encoder/armv8/ih264e_evaluate_intra_chroma_modes_av8.s"
    "${AVC_ROOT}/encoder/armv8/ih264e_half_pel_av8.s"
    "${AVC_ROOT}/encoder/armv8/ime_distortion_metrics_av8.s")

  include_directories(${AVC_ROOT}/encoder/armv8)
elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "aarch32")
  list(
    APPEND
    LIBSVCENC_ASMS
    "${AVC_ROOT}/encoder/arm/ih264e_evaluate_intra16x16_modes_a9q.s"
    "${AVC_ROOT}/encoder/arm/ih264e_evaluate_intra4x4_modes_a9q.s"
    "${AVC_ROOT}/encoder/arm/ih264e_evaluate_intra_chroma_modes_a9q.s"
    "${AVC_ROOT}/encoder/arm/ih264e_fmt_conv.s"
    "${AVC_ROOT}/encoder/arm/ih264e_half_pel.s"
    "${AVC_ROOT}/encoder/arm/ime_distortion_metrics_a9q.s")
else()
  list(
    APPEND
    LIBSVCENC_SRCS
    "${AVC_ROOT}/encoder/x86/ih264e_function_selector.c"
    "${AVC_ROOT}/encoder/x86/ih264e_function_selector_sse42.c"
    "${AVC_ROOT}/encoder/x86/ih264e_function_selector_ssse3.c"
    "${AVC_ROOT}/encoder/x86/ih264e_half_pel_ssse3.c"
    "${AVC_ROOT}/encoder/x86/ih264e_intra_modes_eval_ssse3.c"
    "${AVC_ROOT}/encoder/x86/ime_distortion_metrics_sse42.c"
    "${AVC_ROOT}/encoder/x86/svc/isvce_downscaler_sse42.c"
    "${AVC_ROOT}/encoder/x86/svc/isvce_function_selector.c"
    "${AVC_ROOT}/encoder/x86/svc/isvce_function_selector_sse42.c"
    "${AVC_ROOT}/encoder/x86/svc/isvce_function_selector_ssse3.c"
    "${AVC_ROOT}/encoder/x86/svc/isvce_rc_utils_sse42.c"
    "${AVC_ROOT}/encoder/x86/svc/isvce_residual_pred_sse42.c")

  include_directories(${AVC_ROOT}/encoder/x86)
  include_directories(${AVC_ROOT}/encoder/x86/svc)
endif()

add_library(libsvcenc STATIC ${LIBAVC_COMMON_SRCS} ${LIBAVC_COMMON_ASMS}
                             ${LIBSVCENC_SRCS} ${LIBSVCENC_ASMS})

target_compile_definitions(libsvcenc PRIVATE N_MB_ENABLE)
