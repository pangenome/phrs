#!/bin/bash
# Build sweepga main and fix with guix toolchain, /scratch
set -uo pipefail
V=${1:-main}; G=/gnu/store
GT=$G/ndnvicqyk3v45iayahf153w5cpf639iw-gcc-toolchain-14.3.0
CM=$G/c8zkkq4dwblbvkg2589sllnzpnk2nfav-cmake-4.1.3; PK=$G/y9yw118w55d6b9dfh9y2zp6q1hawxmqq-pkg-config-0.29.2
CL=$G/xcshv2szgls27mzcria3glh2fl7f4f03-clang-18.1.8; RS=$G/9yq5hs32py0myfl1k7y5hvw9ms5xi5ck-rust-1.93.0
LIBS="$G/62jjawrmpbsvm7zxzkm862jp4d7b5rl4-htslib-1.21 $G/wk7ns37jhbyj6i406mixd7gsyzb4gmzl-htscodecs-1.6.1 $G/zycg3161779f4ngajivwpik8mfqkp61m-curl-8.6.0 $G/xwwisbnpbg73wzh9rfmghm51i1k32030-bzip2-1.0.8 $G/s93ifrwg45zk7p7vip2xp9syr0527z65-xz-5.4.5 $G/xkka35zrb34gygb4hxj6hb4lijag291l-gsl-2.8 $G/1vn63xy1cq4kc74w3l93b70wb385jbbz-jemalloc-5.3.0 $G/55xifarldd67da6mdnkapx2as1kng1kp-zlib-1.3.1 $G/r447gvm93x5773nz4jm9dbsyrnqqgbqr-libdeflate-1.19 $G/zvv832w40kr2yn5abs5ixkrh1iv6jpkw-openssl-3.0.8"
inc=""; lib=""; pc=""; rp=""
for d in $LIBS $GT; do inc="$inc:$d/include"; lib="$lib:$d/lib"; pc="$pc:$d/lib/pkgconfig"; rp="$rp -Wl,-rpath,$d/lib"; done
W=/scratch/$USER.sweepgabuild; cd $W/$V || exit 1
export PATH=$RS/bin:$GT/bin:$CM/bin:$PK/bin:/usr/bin:/bin CARGO_HOME=$W/cargo_home_guix
export CC=gcc CXX=g++ SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt CARGO_HTTP_CAINFO=/etc/ssl/certs/ca-certificates.crt GIT_SSL_CAINFO=/etc/ssl/certs/ca-certificates.crt
export C_INCLUDE_PATH=${inc#:} CPLUS_INCLUDE_PATH=${inc#:} LIBRARY_PATH=${lib#:} PKG_CONFIG_PATH=${pc#:} LDFLAGS="$rp" CFLAGS="-I$GT/include" PORTABLE=1
export LIBCLANG_PATH=$CL/lib BINDGEN_EXTRA_CLANG_ARGS="-I$CL/lib/clang/18/include $(echo ${inc#:} | sed 's/:/ -I/g; s/^/-I/')"
export CARGO_TARGET_DIR=$W/target_$V
echo "[$(date -Is)] build $V with $(rustc --version), gcc $(gcc -dumpversion), cmake $(cmake --version | head -1 | cut -d' ' -f3)"
RUST_TEST_THREADS=4 nice -n 19 cargo test --release --no-fail-fast -j 8 > $W/test_$V.log 2>&1; grep -E "^test .* \.\.\. " $W/test_$V.log | sort > $W/tests_$V.txt; echo "tests $V: $(grep -c " ok$" $W/tests_$V.txt) ok, $(grep -c "FAILED$" $W/tests_$V.txt) failed, $(grep -c "ignored$" $W/tests_$V.txt) ignored"
