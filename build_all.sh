set -e

echo "Building suk"
cd suk
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 10
make
cd ../..

echo "Building overlap"
cd overlap
./install_deps.sh
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 10
make
cd ../..

echo "Building hypo polisher"
cd polisher
./install_deps.sh
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 10
make
cd ../..

echo "Building scaffolder"
cd scaffold
./install_deps.sh
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 10
make
cd ../..

mkdir -p run_all
cp run_all.sh run_all/
cp scan_misjoin.py run_all/
cp suk/build/bin/suk run_all/
cp overlap/run_overlap.sh run_all/
cp overlap/join_overlap.py run_all/
cp overlap/filter_overlap.py run_all/
cp overlap/filter_alignments.py run_all/
cp overlap/build/find_overlap run_all/
cp polisher/build/bin/hypo run_all/
cp scaffold/run_scaffold.sh run_all/
cp scaffold/join_scaffold.py run_all/
cp scaffold/filter_scaffold.py run_all/
cp scaffold/build/find_scaffold run_all/

if ! [ -x "$(command -v minimap2)" ]; then
    echo 'Warning: minimap2 is not installed.'
fi

if ! [ -x "$(command -v samtools)" ]; then
    echo 'Warning: samtools is not installed.'
fi

if ! [ -x "$(command -v kmc)" ]; then
    echo 'Warning: KMC is not installed.'
fi

if ! [ -x "$(command -v flye)" ]; then
    echo 'Warning: flye is not installed.'
fi
