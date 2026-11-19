#!/bin/bash
# Build minimap2 v2.31 in /scratch, copy to moosefs
set -euo pipefail
G=/gnu/store; GT=$G/ndnvicqyk3v45iayahf153w5cpf639iw-gcc-toolchain-14.3.0; Z=$G/55xifarldd67da6mdnkapx2as1kng1kp-zlib-1.3.1
O=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/scripts_tmp/minimap2-2.31; W=/scratch/$USER.mm2build; mkdir -p $W $O; cd $W
export PATH=$GT/bin:/usr/bin:/bin C_INCLUDE_PATH=$Z/include:$GT/include LIBRARY_PATH=$Z/lib:$GT/lib GIT_SSL_CAINFO=/etc/ssl/certs/ca-certificates.crt
[ -d minimap2 ] || git clone -q --branch v2.31 --depth 1 https://github.com/lh3/minimap2.git
cd minimap2 && make -j 8 CC=gcc LIBS="-lm -lz -lpthread -Wl,-rpath,$Z/lib" > $W/make.log 2>&1 || { tail -20 $W/make.log; exit 1; }
cp minimap2 $O/ && $O/minimap2 --version && echo "built $(git log -1 --format=%h) -> $O/minimap2"
