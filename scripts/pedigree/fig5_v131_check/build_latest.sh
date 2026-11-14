#!/bin/bash
# Build latest sweepga and wfmash, /scratch
# run on octopus01: nohup bash build_latest.sh SWEEPGA_SHA WFMASH_SHA > log 2>&1 &
# sweepga via build_sweepga_guix.sh (same guix toolchain); wfmash via cmake with the same toolchain and libraries.
# out: $D/scripts_tmp/latest/{sweepga,FastGA,wfmash,...} + versions.txt
set -uo pipefail
SS=$1; WS=$2; D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; O=$D/scripts_tmp/latest; mkdir -p $O
W=/scratch/$USER.sweepgabuild; mkdir -p $W; export GIT_SSL_CAINFO=/etc/ssl/certs/ca-certificates.crt
log () { echo "[$(date -Is)] $*"; }
cd $W; rm -rf latest_src 2>/dev/null; git clone -q https://github.com/pangenome/sweepga.git latest && cd latest && git checkout -q $SS && log "sweepga $(git log -1 --format='%h %cs')"
bash $D/scripts_tmp/build_sweepga_guix.sh latest
B=$W/target_latest/release; ls $B/sweepga >/dev/null && cp $B/sweepga $O/ && log "sweepga copied: $($O/sweepga --version 2>&1 | head -1)"
for x in FastGA FAtoGDB GIXmake ALNtoPAF PAFtoALN GDBtoFA; do f=$(find $W/target_latest -type f -name $x -perm -u+x | head -1); [ -n "$f" ] && cp $f $O/ && log "bundled $x from $f"; done
G=/gnu/store; GT=$G/ndnvicqyk3v45iayahf153w5cpf639iw-gcc-toolchain-14.3.0; CM=$G/c8zkkq4dwblbvkg2589sllnzpnk2nfav-cmake-4.1.3
LIBS="$G/62jjawrmpbsvm7zxzkm862jp4d7b5rl4-htslib-1.21 $G/xkka35zrb34gygb4hxj6hb4lijag291l-gsl-2.8 $G/1vn63xy1cq4kc74w3l93b70wb385jbbz-jemalloc-5.3.0 $G/55xifarldd67da6mdnkapx2as1kng1kp-zlib-1.3.1 $G/r447gvm93x5773nz4jm9dbsyrnqqgbqr-libdeflate-1.19 $G/xwwisbnpbg73wzh9rfmghm51i1k32030-bzip2-1.0.8 $G/s93ifrwg45zk7p7vip2xp9syr0527z65-xz-5.4.5 $G/zycg3161779f4ngajivwpik8mfqkp61m-curl-8.6.0"
pre=""; rp=""; for d in $LIBS $GT; do pre="$pre;$d"; rp="$rp -Wl,-rpath,$d/lib"; done
cd $W; rm -rf wfmash_latest 2>/dev/null; git clone -q --recursive https://github.com/waveygang/wfmash.git wfmash_latest && cd wfmash_latest && git checkout -q $WS && git submodule update -q --init --recursive && log "wfmash $(git log -1 --format='%h %cs')"
export PATH=$GT/bin:$CM/bin:/usr/bin:/bin CC=gcc CXX=g++
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="${pre#;}" -DCMAKE_EXE_LINKER_FLAGS="$rp" > $W/wfmash_cmake.log 2>&1 || { tail -30 $W/wfmash_cmake.log; exit 1; }
nice -n 19 cmake --build build -j 16 > $W/wfmash_build.log 2>&1 || { grep -E -B2 -A8 "error" $W/wfmash_build.log | head -60; exit 1; }
cp build/bin/wfmash $O/ && log "wfmash copied: $($O/wfmash --version 2>&1 | head -1)"
{ echo "sweepga $SS"; echo "wfmash $WS"; $O/sweepga --version 2>&1 | head -1; $O/wfmash --version 2>&1 | head -1; } > $O/versions.txt
log done
