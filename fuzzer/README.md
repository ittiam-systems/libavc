# Fuzzer for libavc decoder and encoder

This describes steps to build avc_dec_fuzzer and avc_enc_fuzzer.

## Linux x86/x64

###  Requirements
- cmake (3.9.1 or above)
- make
- clang (6.0 or above)
  needs to support -fsanitize=fuzzer, -fsanitize=fuzzer-no-link

### Steps to build
Clone libavc repository
```
$ git clone https://android.googlesource.com/platform/external/libavc
```
Create a directory inside libavc and change directory
```
 $ cd libavc
 $ mkdir build
 $ cd build
```
Build fuzzer with required sanitizers (-DSANITIZE=fuzzer-no-link is mandatory
  to enable fuzzers)
```
 $ cmake .. -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
   -DCMAKE_BUILD_TYPE=Debug -DSANITIZE=fuzzer-no-link,address
 $ make
 ```

### Steps to run
Create a directory CORPUS_DIR and copy some elementary h264 files
(for avc_dec_fuzzer) or yuv files (for avc_enc_fuzzer) to that directory

To run the fuzzers
```
$ ./avc_dec_fuzzer CORPUS_DIR
$ ./avc_enc_fuzzer CORPUS_DIR
```

## Android

### Steps to build
Build the fuzzers
```
  $ mm -j$(nproc) avc_dec_fuzzer
  $ mm -j$(nproc) avc_enc_fuzzer
```

### Steps to run
Create a directory CORPUS_DIR and copy some elementary h264 files
(for avc_dec_fuzzer) or yuv files (for avc_enc_fuzzer) to that folder
Push this directory to device

To run avc_dec_fuzzer on device
```
  $ adb sync data
  $ adb shell /data/fuzz/arm64/avc_dec_fuzzer/avc_dec_fuzzer CORPUS_DIR
```
To run avc_enc_fuzzer on device
```
  $ adb sync data
  $ adb shell /data/fuzz/arm64/avc_enc_fuzzer/avc_enc_fuzzer CORPUS_DIR
```

To run avc_dec_fuzzer on host
```
  $ $ANDROID_HOST_OUT/fuzz/x86_64/avc_dec_fuzzer/avc_dec_fuzzer CORPUS_DIR
```

To run avc_enc_fuzzer on host
```
  $ $ANDROID_HOST_OUT/fuzz/x86_64/avc_enc_fuzzer/avc_enc_fuzzer CORPUS_DIR
```


# Appendix
## libavc encoder fuzzer

## Plugin Design Considerations
The fuzzer plugin for AVC is designed based on the understanding of the
codec and tries to achieve the following:

##### Maximize code coverage
The configuration parameters are not hardcoded, but instead selected based on
incoming data. This ensures more code paths are reached by the fuzzer.

AVC supports the following parameters:
1. Frame Width (parameter name: `u4_wd`)
2. Frame Height (parameter name: `u4_ht`)
3. Input color format  (parameter name: `e_inp_color_fmt`)
4. Architecture type (parameter name: `e_arch`)
5. Rate control mode (parameter name: `e_rc_mode`)
6. Number of cores (parameter name: `u4_num_cores`)
7. Maximum B frames (parameter name: `u4_num_bframes`)
8. Encoder speed preset (parameter name: `u4_enc_speed_preset`)
9. enable constrained intra prediction (parameter name: `u4_constrained_intra_pred`)
10. enable intra 4x4 (parameter name: `u4_enable_intra_4x4`)
11. Qp for I frames (parameter name: `u4_i_qp`)
12. Qp for P frames (parameter name: `u4_p_qp`)
13. Qp for B frames (parameter name: `u4_b_qp`)
14. Target Bitrate (parameter name: `u4_target_bitrate`)
15. Intra refresh period in frames (parameter name: `u4_air_refresh_period`)
16. Enable half pel ME (parameter name: `u4_enable_hpel`)
17. Enable quarter pel ME (parameter name: `u4_enable_qpel`)
18. ME speed preset (parameter name: `u4_me_speed_preset`)
19. Adaptive intra refresh mode (parameter name: `e_air_mode`)
20. Disable deblock level (parameter name: `u4_disable_deblock_level`)
21. Max search range in X direction (parameter name: `u4_srch_rng_x`)
22. Max search range in Y direction (parameter name: `u4_srch_rng_y`)
23. I frame interval (parameter name: `u4_i_frm_interval`)
24. IDR frame interval (parameter name: `u4_idr_frm_interval`)
25. Enable mastering display color volume info (parameter name: `u1_sei_mdcv_params_present_flag`)
26. Enable content light level info (parameter name: `u1_sei_cll_params_present_flag`)
27. Enable ambient viewing environment info (parameter name: `u1_sei_ave_params_present_flag`)
28. Enable content color volume info (parameter name: `u1_sei_ccv_params_present_flag`)
29. Profile (parameter name: `e_profile`)
30. Enable aspect_ratio info (parameter name: `u1_aspect_ratio_info_present_flag`)
31. Enable NAL HRD parameters presence (parameter name: `u1_nal_hrd_parameters_present_flag`)
32. Enable VCL HRD parameters presence (parameter name: `u1_vcl_hrd_parameters_present_flag`)
33. Enable force IDR frame (parameter name: `mIsForceIdrEnabled`)
34. Enable dynamic bitrate change (parameter name: `mIsDynamicBitRateChangeEnabled`)
35. Enable dynamic framerate change (parameter name: `mIsDynamicFrameRateChangeEnabled`)
36. Force IDR frame number (parameter name: `mForceIdrInterval`)
37. Frame number for dynamic bitrate (parameter name: `mDynamicBitRateInterval`)
38. Frame number for dynamic framerate (parameter name: `mDynamicFrameRateInterval`)

| Parameter| Valid Values| Configured Value|
|------------- |-------------| ----- |
| `u4_wd` | In the range `0 to 10239` | All the bits of 1st and 2nd byte of data |
| `u4_ht` | In the range `0 to 10239` | All the bits of 3rd and 4th byte of data |
| `e_inp_color_fmt` | 0. `IV_YUV_420P` 1. `IV_YUV_420SP_UV` 2. `IV_YUV_422ILE` 3. `IV_YUV_420SP_VU` | All the bits of 5th byte of data |
| `e_arch` | 0. `ARCH_ARM_NONEON` 1. `ARCH_NA` | bit 0 and 1 of 6th byte of data |
| `e_rc_mode` | 0. `IVE_RC_NONE` 1. `IVE_RC_STORAGE` 2. `IVE_RC_CBR_NON_LOW_DELAY` 3. `IVE_RC_CBR_LOW_DELAY` | All the bits of 7th byte of data modulus 4 |
| `u4_num_cores` | 0. `0` 1. `1` 2. `2` 3. `3`| bit 0 and 1 of 8th byte of data |
| `u4_num_bframes` | In the range `0 to 7` | bit 0, 1 and 2 of 9th byte of data |
| `u4_enc_speed_preset` | 0. `IVE_CONFIG` 1. `IVE_SLOWEST` 2. `IVE_NORMAL` 3. `IVE_FAST` 4. `IVE_HIGH_SPEED` 5. `IVE_FASTEST` | All the bits of 10th byte of data modulus 6 |
| `u4_constrained_intra_pred` | 0. `0` 1. `1` | bit 0 of 11th byte of data |
| `u4_enable_intra_4x4` | 0. `0` 1. `1` | bit 0 of 12th byte of data |
| `u4_i_qp` | In the range `4 to 51` | All the bits of 13th byte of data |
| `u4_p_qp` | In the range `4 to 51` | All the bits of 14th byte of data |
| `u4_b_qp` | In the range `4 to 51` | All the bits of 15th byte of data |
| `u4_target_bitrate` | In the range `0 to 500000000` | All the bits of 16th and 17th byte of data |
| `u4_target_bitrate` | In the range `0 to 255` | All the bits of 18th byte of data |
| `u4_air_refresh_period` | In the range `1 to 256` | All the bits of 19th byte of data |
| `u4_air_refresh_period` | In the range `1 to 256` | All the bits of 19th byte of data |
| `u4_enable_hpel` | 0. `0` 1. `1` | bit 0 of 20th byte of data |
| `u4_enable_qpel` | 0. `0` 1. `1` | bit 0 of 21st byte of data |
| `u4_me_speed_preset` | 0. `0` 1. `50` 2. `75` 3. `100` | All the bits of 22nd byte of data modulus 4 |
| `e_air_mode` | 0. `IVE_AIR_MODE_NONE` 1. `IVE_AIR_MODE_CYCLIC` 2. `IVE_AIR_MODE_RANDOM` | All the bits of 23rd byte of data modulus 3 |
| `u4_disable_deblock_level` | 0. `0` 1. `1` 2. `2` 3. `3` | bit 0 and 1 of 24th byte of data |
| `u4_srch_rng_x` | In the range `0 to 255` | All the bits of 25th byte of data |
| `u4_srch_rng_y` | In the range `0 to 255`| All the bits of 26th byte of data |
| `u4_i_frm_interval` | In the range `1 to 256` | All the bits of 27th byte of data |
| `u4_idr_frm_interval` | In the range `1 to 256` | All the bits of 28th byte of data |
| `u1_sei_mdcv_params_present_flag` | 0. `0` 1. `1` | bit 0 of 29th byte of data |
| `u1_sei_cll_params_present_flag` | 0. `0` 1. `1` | bit 0 of 30th byte of data |
| `u1_sei_ave_params_present_flag` | 0. `0` 1. `1` | bit 0 of 31st byte of data |
| `u1_sei_ccv_params_present_flag` | 0. `0` 1. `1` | bit 0 of 32nd byte of data |
| `e_profile` | 0. `IV_PROFILE_BASE` 1. `IV_PROFILE_MAIN` | bit 0 and 1 of 33th byte of data modulus 2 |
| `u1_aspect_ratio_info_present_flag` | 0. `0` 1. `1` | bit 0 of 34th byte of data |
| `u1_nal_hrd_parameters_present_flag` | 0. `0` 1. `1` | bit 0 of 35th byte of data |
| `u1_vcl_hrd_parameters_present_flag` | 0. `0` 1. `1` | bit 0 of 36th byte of data |
| `mIsForceIdrEnabled` | 0. `0` 1. `1` | bit 0 of 37th byte of data |
| `mIsDynamicBitRateChangeEnabled` | 0. `0` 1. `1` | bit 0 of 38th byte of data |
| `mIsDynamicFrameRateChangeEnabled` | 0. `0` 1. `1` | bit 0 of 39th byte of data |
| `mForceIdrInterval` | In the range `0 to 7` | bit 0, 1 and 2 of 40th byte of data |
| `mDynamicBitRateInterval` | In the range `0 to 7` | bit 0, 1 and 2 of 41st byte of data |
| `mDynamicFrameRateInterval` | In the range `0 to 7` | bit 0, 1 and 2 of 42nd byte of data |

This also ensures that the plugin is always deterministic for any given input.

##### Maximize utilization of input data
The plugin feeds the entire input data to the codec using a loop.
If the encode operation was successful, the input is advanced by the frame size.
If the encode operation was un-successful, the input is still advanced by frame size so
that the fuzzer can proceed to feed the next frame.

This ensures that the plugin tolerates any kind of input (empty, huge,
malformed, etc) and doesnt `exit()` on any input and thereby increasing the
chance of identifying vulnerabilities.


## References:
 * http://llvm.org/docs/LibFuzzer.html
 * https://github.com/google/oss-fuzz
