#!/bin/bash
threads=10
if [ "$1" ]; then
    threads="$1"
fi

cd libs/htslib-1.9;
make clean;
autoreconf -i;
./configure --prefix=$(pwd) --disable-bz2 --disable-lzma;
make -j $threads; make install;
echo "HTSlib installed."
