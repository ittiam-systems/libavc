list(
  APPEND
  AVCENC_SRCS
  "${AVC_ROOT}/examples/avcenc/input.c"
  "${AVC_ROOT}/examples/avcenc/main.c"
  "${AVC_ROOT}/examples/avcenc/output.c"
  "${AVC_ROOT}/examples/avcenc/psnr.c"
  "${AVC_ROOT}/examples/avcenc/recon.c")

libavc_add_executable(avcenc libavcenc SOURCES ${AVCENC_SRCS})
target_compile_definitions(avcenc PRIVATE PROFILE_ENABLE MD5_DISABLE)
