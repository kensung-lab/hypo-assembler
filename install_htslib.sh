#!/bin/bash
cd external/htslib;
autoreconf;
./configure --prefix=$(pwd) --disable-bz2 --disable-lzma;
make; make install;
echo "HTSlib installed."

