#!/bin/bash

threads=10
if [ "$1" ]; then
    threads="$1"
fi

cd external/install/htslib;
autoreconf;
./configure --prefix=$(pwd) --disable-bz2 --disable-lzma;
make -j $threads; make install;
echo "HTSlib installed."
