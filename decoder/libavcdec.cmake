# src files
list(
  APPEND
  LIBAVCDEC_SRCS
  "${AVC_ROOT}/decoder/ih264d_api.c"
  "${AVC_ROOT}/decoder/ih264d_bitstrm.c"
  "${AVC_ROOT}/decoder/ih264d_cabac.c"
  "${AVC_ROOT}/decoder/ih264d_cabac_init_tables.c"
  "${AVC_ROOT}/decoder/ih264d_compute_bs.c"
  "${AVC_ROOT}/decoder/ih264d_deblocking.c"
  "${AVC_ROOT}/decoder/ih264d_dpb_mgr.c"
  "${AVC_ROOT}/decoder/ih264d_format_conv.c"
  "${AVC_ROOT}/decoder/ih264d_function_selector_generic.c"
  "${AVC_ROOT}/decoder/ih264d_inter_pred.c"
  "${AVC_ROOT}/decoder/ih264d_mb_utils.c"
  "${AVC_ROOT}/decoder/ih264d_mvpred.c"
  "${AVC_ROOT}/decoder/ih264d_nal.c"
  "${AVC_ROOT}/decoder/ih264d_parse_bslice.c"
  "${AVC_ROOT}/decoder/ih264d_parse_cabac.c"
  "${AVC_ROOT}/decoder/ih264d_parse_cavlc.c"
  "${AVC_ROOT}/decoder/ih264d_parse_headers.c"
  "${AVC_ROOT}/decoder/ih264d_parse_islice.c"
  "${AVC_ROOT}/decoder/ih264d_parse_mb_header.c"
  "${AVC_ROOT}/decoder/ih264d_parse_pslice.c"
  "${AVC_ROOT}/decoder/ih264d_parse_slice.c"
  "${AVC_ROOT}/decoder/ih264d_process_bslice.c"
  "${AVC_ROOT}/decoder/ih264d_process_intra_mb.c"
  "${AVC_ROOT}/decoder/ih264d_process_pslice.c"
  "${AVC_ROOT}/decoder/ih264d_quant_scaling.c"
  "${AVC_ROOT}/decoder/ih264d_sei.c"
  "${AVC_ROOT}/decoder/ih264d_tables.c"
  "${AVC_ROOT}/decoder/ih264d_thread_compute_bs.c"
  "${AVC_ROOT}/decoder/ih264d_thread_parse_decode.c"
  "${AVC_ROOT}/decoder/ih264d_utils.c"
  "${AVC_ROOT}/decoder/ih264d_vui.c")

include_directories(${AVC_ROOT}/decoder)

if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" OR "${CMAKE_SYSTEM_PROCESSOR}"
                                                     STREQUAL "aarch32")
  list(
    APPEND LIBAVCDEC_ASMS "${AVC_ROOT}/decoder/arm/ih264d_function_selector.c"
    "${AVC_ROOT}/decoder/arm/ih264d_function_selector_a9q.c"
    "${AVC_ROOT}/decoder/arm/ih264d_function_selector_av8.c")
else()
  list(
    APPEND LIBAVCDEC_SRCS "${AVC_ROOT}/decoder/x86/ih264d_function_selector.c"
    "${AVC_ROOT}/decoder/x86/ih264d_function_selector_sse42.c"
    "${AVC_ROOT}/decoder/x86/ih264d_function_selector_ssse3.c")
endif()

add_library(libavcdec STATIC ${LIBAVC_COMMON_SRCS} ${LIBAVC_COMMON_ASMS}
                             ${LIBAVCDEC_SRCS} ${LIBAVCDEC_ASMS})
