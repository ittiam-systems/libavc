list(
  APPEND
  SVCENC_SRCS
  "${AVC_ROOT}/examples/svcenc/main.c"
  "${AVC_ROOT}/examples/svcenc/input.c"
  "${AVC_ROOT}/examples/svcenc/output.c"
  "${AVC_ROOT}/examples/svcenc/psnr.c"
  "${AVC_ROOT}/examples/svcenc/recon.c")

libavc_add_executable(svcenc libsvcenc SOURCES ${SVCENC_SRCS} INCLUDES
                      "${AVC_ROOT}/examples/svcenc/")
target_compile_definitions(svcenc PRIVATE PROFILE_ENABLE MD5_DISABLE)
