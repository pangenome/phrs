#!/bin/bash
# Build hifiasm 0.25.0 in /scratch, copy to moosefs
set -euo pipefail
G=/gnu/store; GT=$G/ndnvicqyk3v45iayahf153w5cpf639iw-gcc-toolchain-14.3.0; Z=$G/55xifarldd67da6mdnkapx2as1kng1kp-zlib-1.3.1
O=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/scripts_tmp/hifiasm-0.25.0; W=/scratch/$USER.hifiasmbuild; mkdir -p $W $O; cd $W
export PATH=$GT/bin:/usr/bin:/bin C_INCLUDE_PATH=$Z/include:$G/90qwxdmdw7jdf0dkkgdgxln6b5a57kzx-linux-libre-headers-6.12.17/include CPLUS_INCLUDE_PATH=$Z/include:$G/90qwxdmdw7jdf0dkkgdgxln6b5a57kzx-linux-libre-headers-6.12.17/include LIBRARY_PATH=$Z/lib:$GT/lib GIT_SSL_CAINFO=/etc/ssl/certs/ca-certificates.crt
[ -d hifiasm ] || git clone -q --branch 0.25.0 --depth 1 https://github.com/chhylp123/hifiasm.git
cd hifiasm && make clean > /dev/null 2>&1; make -j 8 CC=gcc CXX=g++ CPPFLAGS="-I$Z/include" LDFLAGS="-L$Z/lib -Wl,-rpath,$Z/lib" > $W/make.log 2>&1 || { tail -20 $W/make.log; exit 1; }
cp hifiasm $O/ && $O/hifiasm --version && echo "built $(git log -1 --format=%h) -> $O/hifiasm"
