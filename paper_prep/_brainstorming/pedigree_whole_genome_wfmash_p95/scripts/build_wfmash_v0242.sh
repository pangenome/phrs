#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$PACKAGE_DIR/build"
TOOLS_DIR="$PACKAGE_DIR/tools/wfmash-v0.24.2-built"
SRC_TARBALL="$BUILD_DIR/wfmash-v0.24.2.tar.gz"
SRC_DIR="$BUILD_DIR/wfmash-v0.24.2-src"
CMAKE_BUILD="$BUILD_DIR/wfmash-v0.24.2-build-vendor-gcc13"

WFMASH_SRC_URL="${WFMASH_SRC_URL:-https://github.com/waveygang/wfmash/archive/refs/tags/v0.24.2.tar.gz}"
WFMASH_SRC_SHA256="${WFMASH_SRC_SHA256:-db86e4a82538ff019a86fd9d7924805e062b3fb8aafa472781907bf0bd557620}"
BZIP2_PREFIX="${BZIP2_PREFIX:-/gnu/store/pl09vk5g3cl8fxfln2hjk996pyahqk8m-bzip2-1.0.8}"
XZ_PREFIX="${XZ_PREFIX:-/gnu/store/nv5q3a8wf16arzgvgqc3125xbglqg5z2-xz-5.2.8}"

export PATH="/gnu/store/10krix03rl5hqjv2c0qmj44ic9bgd8rc-gcc-toolchain-13.3.0/bin:/gnu/store/is7b1l1rirbgkah1fhh26943p00ycvhf-make-4.4.1/bin:/gnu/store/2ndf2fc56d82iy8xyynkxf1vf5lfnfi6-cmake-minimal-3.31.10/bin:/gnu/store/a3lsdsalcmg5wnk67869af7wljprkbam-pkg-config-0.29.2/bin:$PATH"
export CPATH="${BZIP2_PREFIX}/include:${XZ_PREFIX}/include${CPATH:+:$CPATH}"
export LIBRARY_PATH="${BZIP2_PREFIX}/lib:${XZ_PREFIX}/lib:/gnu/store/10krix03rl5hqjv2c0qmj44ic9bgd8rc-gcc-toolchain-13.3.0/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export CPPFLAGS="-I${BZIP2_PREFIX}/include -I${XZ_PREFIX}/include ${CPPFLAGS:-}"
export LDFLAGS="-L${BZIP2_PREFIX}/lib -L${XZ_PREFIX}/lib ${LDFLAGS:-}"

mkdir -p "$BUILD_DIR" "$TOOLS_DIR" "$PACKAGE_DIR/logs"
if [[ ! -s "$SRC_TARBALL" ]]; then
    curl -L --fail "$WFMASH_SRC_URL" -o "$SRC_TARBALL"
fi
echo "${WFMASH_SRC_SHA256}  ${SRC_TARBALL}" | sha256sum -c -

rm -rf "$SRC_DIR"
mkdir -p "$SRC_DIR"
tar -xzf "$SRC_TARBALL" --strip-components=1 -C "$SRC_DIR"
rm -rf "$CMAKE_BUILD"
mkdir -p "$CMAKE_BUILD"
cd "$CMAKE_BUILD"

CC="$(command -v gcc)" CXX="$(command -v g++)" cmake \
    -DPKG_CONFIG_EXECUTABLE="$(command -v pkg-config)" \
    -DCMAKE_BUILD_TYPE=Generic \
    -DDISABLE_LTO=ON \
    -DVENDOR_EVERYTHING=ON \
    "$SRC_DIR"
make -j "${WFMASH_BUILD_THREADS:-8}" || true

# wfmash v0.24.2 CMake currently expects libdeflate/lib64, while the vendored
# libdeflate install on this system uses libdeflate/lib. Provide the expected
# path and continue the incremental build.
mkdir -p libdeflate/lib64
ln -sf ../lib/libdeflate.a libdeflate/lib64/libdeflate.a
make -j "${WFMASH_BUILD_THREADS:-8}"

cp bin/wfmash "$TOOLS_DIR/wfmash"
"$TOOLS_DIR/wfmash" --version | tee "$PACKAGE_DIR/logs/wfmash_v0.24.2.version.txt"
"$TOOLS_DIR/wfmash" --help > "$PACKAGE_DIR/logs/wfmash_v0.24.2.help.txt"
sha256sum "$TOOLS_DIR/wfmash" | tee "$PACKAGE_DIR/logs/wfmash_v0.24.2.built.sha256"
