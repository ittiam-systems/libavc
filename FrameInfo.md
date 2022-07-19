## Frame Info exported from libAVC

### Introduction
QP and block type maps for H264 are defined for each 8x8 MB sub-block.
The QP values defined as unsigned 8-bit numbers can range from <0, 51> and the block type can
be INTER/INTRA/SKIP. Set the “u4_frame_info_enable” flag to enable encoder/decoder to populate
and return the qp values and block type data in their output structures ih264e_video_encode_op_t
and ih264d_video_decode_op_t respectively via pu1_8x8_blk_qp_map and pu1_8x8_blk_type_map.

### Mapping to the frame
Let’s say, a frame has a total of ‘n’ MBs (each 16x16). Since the QP and block type are defined
for each 8x8 block, hence each MB will have 4 entries in the maps. Thus, a total of n x 4 entries
for each frame. Qp and block type values for each 8x8 block are stored in raster scan order. Refer
to ih264d.h for details.

### Plugin/Application
The encoder/decoder keeps the QP and block type map as a part of its output handle. The plugins can
access these data through the output structure.
