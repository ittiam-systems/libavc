## Media Testing ##
---

#### AvcEncoder
The AvcEncoder Test Suite validates the Avc encoder library available in external/libavc.
Run the following steps to build the test suite:
```
m AvcEncoderTest
```

The 32-bit binaries will be created in the following path : ${OUT}/data/nativetest/
The 64-bit binaries will be created in the following path : ${OUT}/data/nativetest64/

To test 64-bit binary push binaries from nativetest64.
```
adb push ${OUT}/data/nativetest64/AvcEncoderTest/AvcEncoderTest /data/local/tmp/
```

To test 32-bit binary push binaries from nativetest.
```
adb push ${OUT}/data/nativetest/AvcEncoderTest/AvcEncoderTest /data/local/tmp/
```

The resource file for the tests is taken from [here](https://storage.googleapis.com/android_media/frameworks/av/media/libstagefright/codecs/m4v_h263/enc/test/Mpeg4H263Encoder.zip ) Download, unzip and push these files into device for testing.

```
adb push AvcEncoder/. /data/local/tmp/
```

usage: AvcEncoderTest -P \<path_to_folder\>
```
adb shell /data/local/tmp/AvcEncoderTest -P /data/local/tmp/
```
Alternatively, the test can also be run using atest command.

```
atest AvcEncoderTest -- --enable-module-dynamic-download=true
```
