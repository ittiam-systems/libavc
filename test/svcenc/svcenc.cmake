list(
  APPEND
  SVCENC_SRCS
  "${AVC_ROOT}/test/svcenc/main.c"
  "${AVC_ROOT}/test/svcenc/input.c"
  "${AVC_ROOT}/test/svcenc/output.c"
  "${AVC_ROOT}/test/svcenc/psnr.c"
  "${AVC_ROOT}/test/svcenc/recon.c")

libavc_add_executable(svcenc libsvcenc SOURCES ${SVCENC_SRCS} INCLUDES
                      "${AVC_ROOT}/test/svcenc/")
target_compile_definitions(svcenc PRIVATE PROFILE_ENABLE MD5_DISABLE)
