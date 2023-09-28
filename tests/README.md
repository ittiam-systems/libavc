# AvcEncTest
The AvcEncoder Test Suite validates the Avc encoder.

## Linux x86/x64

###  Requirements
- cmake (3.9.1 or above)
- make
- clang (12.0 or above)

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

Build with -DENABLE_TESTS=1.
```
 $ cmake .. -DENABLE_TESTS=1 -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
   -DCMAKE_BUILD_TYPE=Debug
 $ make
```

Optionally, enable sanitizers by passing -DSANITIZE
```
 $ cmake .. -DENABLE_TESTS=1 -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
   -DCMAKE_BUILD_TYPE=Debug -DSANITIZE=fuzzer-no-link,address,\
   signed-integer-overflow,unsigned-integer-overflow
 $ make
```

The media files for the tests are present [at](https://storage.googleapis.com/android_media/external/libavc/tests/AvcEncoder.zip).
Download and extract these the current folder.

usage: AvcEncTest -P \<path_to_the local folder\>

```
$./AvcEncTest -P ./
```

## Android

Run the following steps to build the test suite:
```
m AvcEncTest
```

To test 64-bit binary push binaries from nativetest64.
```
adb push ${OUT}/data/nativetest64/AvcEncTest/AvcEncTest /data/local/tmp/
```

To test 32-bit binary push binaries from nativetest.
```
adb push ${OUT}/data/nativetest/AvcEncTest/AvcEncTest /data/local/tmp/
```

The resource file for the tests is taken from [here](https://storage.googleapis.com/android_media/external/libavc/tests/AvcEncoder.zip)

Download, unzip and push these files into device for testing.

```
adb push AvcEncoder/. /data/local/tmp/
```

usage: AvcEncTest -P \<path_to_folder\>
```
adb shell /data/local/tmp/AvcEncTest -P /data/local/tmp/
```
Alternatively, the test can also be run using atest command.

```
atest AvcEncTest -- --enable-module-dynamic-download=true
```
