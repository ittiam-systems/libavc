# src files
list(
  APPEND
  LIBAVCENC_SRCS
  "${AVC_ROOT}/encoder/ih264e_api.c"
  "${AVC_ROOT}/encoder/ih264e_bitstream.c"
  "${AVC_ROOT}/encoder/ih264e_cabac.c"
  "${AVC_ROOT}/encoder/ih264e_cabac_encode.c"
  "${AVC_ROOT}/encoder/ih264e_cabac_init.c"
  "${AVC_ROOT}/encoder/ih264e_cavlc.c"
  "${AVC_ROOT}/encoder/ih264e_core_coding.c"
  "${AVC_ROOT}/encoder/ih264e_deblk.c"
  "${AVC_ROOT}/encoder/ih264e_encode.c"
  "${AVC_ROOT}/encoder/ih264e_encode_header.c"
  "${AVC_ROOT}/encoder/ih264e_fmt_conv.c"
  "${AVC_ROOT}/encoder/ih264e_function_selector_generic.c"
  "${AVC_ROOT}/encoder/ih264e_globals.c"
  "${AVC_ROOT}/encoder/ih264e_half_pel.c"
  "${AVC_ROOT}/encoder/ih264e_intra_modes_eval.c"
  "${AVC_ROOT}/encoder/ih264e_mc.c"
  "${AVC_ROOT}/encoder/ih264e_me.c"
  "${AVC_ROOT}/encoder/ih264e_modify_frm_rate.c"
  "${AVC_ROOT}/encoder/ih264e_process.c"
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
  "${AVC_ROOT}/encoder/irc_vbr_str_prms.c")

include_directories(${AVC_ROOT}/encoder)

if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "aarch64")
  list(
    APPEND
    LIBAVCENC_ASMS
    "${AVC_ROOT}/encoder/arm/ih264e_function_selector.c"
    "${AVC_ROOT}/encoder/arm/ih264e_function_selector_a9q.c"
    "${AVC_ROOT}/encoder/arm/ih264e_function_selector_av8.c"
    "${AVC_ROOT}/encoder/armv8/ih264e_evaluate_intra16x16_modes_av8.s"
    "${AVC_ROOT}/encoder/armv8/ih264e_evaluate_intra_chroma_modes_av8.s"
    "${AVC_ROOT}/encoder/armv8/ih264e_half_pel_av8.s"
    "${AVC_ROOT}/encoder/armv8/ime_distortion_metrics_av8.s")

  include_directories(${AVC_ROOT}/encoder/armv8)
elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "aarch32")
  list(
    APPEND
    LIBAVCENC_ASMS
    "${AVC_ROOT}/encoder/arm/ih264e_function_selector.c"
    "${AVC_ROOT}/encoder/arm/ih264e_function_selector_a9q.c"
    "${AVC_ROOT}/encoder/arm/ih264e_function_selector_av8.c"
    "${AVC_ROOT}/encoder/arm/ih264e_evaluate_intra16x16_modes_a9q.s"
    "${AVC_ROOT}/encoder/arm/ih264e_evaluate_intra4x4_modes_a9q.s"
    "${AVC_ROOT}/encoder/arm/ih264e_evaluate_intra_chroma_modes_a9q.s"
    "${AVC_ROOT}/encoder/arm/ih264e_fmt_conv.s"
    "${AVC_ROOT}/encoder/arm/ih264e_half_pel.s"
    "${AVC_ROOT}/encoder/arm/ime_distortion_metrics_a9q.s")

  include_directories(${AVC_ROOT}/encoder/armv8)
else()
  list(
    APPEND
    LIBAVCENC_SRCS
    "${AVC_ROOT}/encoder/x86/ih264e_function_selector.c"
    "${AVC_ROOT}/encoder/x86/ih264e_function_selector_sse42.c"
    "${AVC_ROOT}/encoder/x86/ih264e_function_selector_ssse3.c"
    "${AVC_ROOT}/encoder/x86/ih264e_half_pel_ssse3.c"
    "${AVC_ROOT}/encoder/x86/ih264e_intra_modes_eval_ssse3.c"
    "${AVC_ROOT}/encoder/x86/ime_distortion_metrics_sse42.c")

  include_directories(${AVC_ROOT}/encoder/x86)
endif()

add_library(libavcenc STATIC ${LIBAVC_COMMON_SRCS} ${LIBAVC_COMMON_ASMS}
                             ${LIBAVCENC_SRCS} ${LIBAVCENC_ASMS})

target_compile_definitions(libavcenc PRIVATE N_MB_ENABLE)
