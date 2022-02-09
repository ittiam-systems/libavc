list(
  APPEND
  AVCENC_SRCS
  "${AVC_ROOT}/test/encoder/input.c"
  "${AVC_ROOT}/test/encoder/main.c"
  "${AVC_ROOT}/test/encoder/output.c"
  "${AVC_ROOT}/test/encoder/psnr.c"
  "${AVC_ROOT}/test/encoder/recon.c")

libavc_add_executable(avcenc libavcenc SOURCES ${AVCENC_SRCS})
