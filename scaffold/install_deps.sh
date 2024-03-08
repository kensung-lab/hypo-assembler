#!/bin/bash
cd libs/htslib-1.9;
autoreconf;
./configure --prefix=$(pwd) --disable-bz2 --disable-lzma;
make -j20; make install;
echo "HTSlib installed."

