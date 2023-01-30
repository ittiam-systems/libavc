# src files
list(
  APPEND
  LIBAVC_COMMON_SRCS
  "${AVC_ROOT}/common/svc/isvc_common_tables.c"
  "${AVC_ROOT}/common/svc/isvc_cabac_tables.c"
  "${AVC_ROOT}/common/svc/isvc_intra_resample.c"
  "${AVC_ROOT}/common/svc/isvc_iquant_itrans_recon.c"
  "${AVC_ROOT}/common/svc/isvc_mem_fns.c"
  "${AVC_ROOT}/common/svc/isvc_resi_trans_quant.c")

include_directories(${AVC_ROOT}/common/svc)

# arm/x86 sources
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" OR
   "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch32")
  list(
    APPEND
    LIBAVC_COMMON_ASMS
    "${AVC_ROOT}/common/arm/svc/isvc_intra_sampling_neon.c"
    "${AVC_ROOT}/common/arm/svc/isvc_iquant_itrans_recon_neon.c"
    "${AVC_ROOT}/common/arm/svc/isvc_mem_fns_neon.c"
    "${AVC_ROOT}/common/arm/svc/isvc_resi_trans_quant_neon.c")
  include_directories(${AVC_ROOT}/common/arm/svc)
else()
  list(
    APPEND
    LIBAVC_COMMON_SRCS
    "${AVC_ROOT}/common/x86/svc/isvc_iquant_itrans_recon_dc_ssse3.c"
    "${AVC_ROOT}/common/x86/svc/isvc_iquant_itrans_recon_sse42.c"
    "${AVC_ROOT}/common/x86/svc/isvc_iquant_itrans_recon_ssse3.c"
    "${AVC_ROOT}/common/x86/svc/isvc_mem_fns_sse42.c"
    "${AVC_ROOT}/common/x86/svc/isvc_mem_fns_ssse3.c"
    "${AVC_ROOT}/common/x86/svc/isvc_padding_ssse3.c"
    "${AVC_ROOT}/common/x86/svc/isvc_resi_trans_quant_sse42.c"
    "${AVC_ROOT}/common/x86/svc/isvc_intra_resample_sse42.c")

  include_directories(${AVC_ROOT}/common/x86/svc)
endif()
