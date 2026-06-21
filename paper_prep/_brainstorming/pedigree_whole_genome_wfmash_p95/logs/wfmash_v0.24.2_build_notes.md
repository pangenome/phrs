# wfmash v0.24.2 build notes

Task: `fig5-whole-genome-wfmash-p95`

The local Guix profile had an older wfmash (`0.12.5-1+0222f7c`) and is treated
only as a legacy diagnostic source. The primary evidence run uses upstream
wfmash `v0.24.2`.

Source metadata:

- GitHub release: <https://github.com/waveygang/wfmash/releases/tag/v0.24.2>
- Source URL: <https://github.com/waveygang/wfmash/archive/refs/tags/v0.24.2.tar.gz>
- Tag commit: `774c01ffb9df010d6b529520033ed7dce0cb95d5`
- Source tarball SHA256: `db86e4a82538ff019a86fd9d7924805e062b3fb8aafa472781907bf0bd557620`
- GitHub linux binary SHA256: `c1703602931ab0570d8fe043afb7f4397e890fd0e6e17eeed0cc46a6231b61fa`

The GitHub linux binary was downloaded but rejected because it requires newer
system glibc versions than the login node provides. The exact runtime failure is
captured in `logs/wfmash_v0.24.2_release_binary_failure.txt`.

The working binary was built from source with `scripts/build_wfmash_v0242.sh`.
The successful local build used:

```bash
export PATH=/gnu/store/10krix03rl5hqjv2c0qmj44ic9bgd8rc-gcc-toolchain-13.3.0/bin:/gnu/store/is7b1l1rirbgkah1fhh26943p00ycvhf-make-4.4.1/bin:/gnu/store/2ndf2fc56d82iy8xyynkxf1vf5lfnfi6-cmake-minimal-3.31.10/bin:/gnu/store/a3lsdsalcmg5wnk67869af7wljprkbam-pkg-config-0.29.2/bin:$PATH
export BZIP2_PREFIX=/gnu/store/pl09vk5g3cl8fxfln2hjk996pyahqk8m-bzip2-1.0.8
export XZ_PREFIX=/gnu/store/nv5q3a8wf16arzgvgqc3125xbglqg5z2-xz-5.2.8
export CPATH=${BZIP2_PREFIX}/include:${XZ_PREFIX}/include${CPATH:+:$CPATH}
export LIBRARY_PATH=${BZIP2_PREFIX}/lib:${XZ_PREFIX}/lib:/gnu/store/10krix03rl5hqjv2c0qmj44ic9bgd8rc-gcc-toolchain-13.3.0/lib${LIBRARY_PATH:+:$LIBRARY_PATH}
export CPPFLAGS="-I${BZIP2_PREFIX}/include -I${XZ_PREFIX}/include ${CPPFLAGS:-}"
export LDFLAGS="-L${BZIP2_PREFIX}/lib -L${XZ_PREFIX}/lib ${LDFLAGS:-}"
CC=$(command -v gcc) CXX=$(command -v g++) cmake \
  -DPKG_CONFIG_EXECUTABLE=$(command -v pkg-config) \
  -DCMAKE_BUILD_TYPE=Generic \
  -DDISABLE_LTO=ON \
  -DVENDOR_EVERYTHING=ON \
  ../wfmash-v0.24.2-src
make -j 8
mkdir -p libdeflate/lib64
ln -sf ../lib/libdeflate.a libdeflate/lib64/libdeflate.a
make -j 8
```

The `libdeflate/lib64` symlink is a local build-tree compatibility workaround:
v0.24.2 CMake expected `libdeflate/lib64/libdeflate.a`, while the vendored
libdeflate install on this host produced `libdeflate/lib/libdeflate.a`.

Captured build artifacts:

- `logs/wfmash_v0.24.2.version.txt`
- `logs/wfmash_v0.24.2.help.txt`
- `logs/wfmash_v0.24.2.built.sha256`
- `logs/wfmash_v0.24.2.ldd.txt`
- `logs/wfmash_v0.24.2_build_retry*.log`
