set -e
set -x

usage (){
    echo "  -t <threads>            Number of threads to run make [Default: 10]"
    echo "  -o <directory>          Directory for outputs [Default: run_all]"
    echo "  -n                      Build with march=native"
    echo "  -h                      display this help and exit"
exit 1
}

native=false
threads=10
builddir="run_all"
while getopts "nt:o:h" opt; do
  case $opt in
    n)
        native=true
        ;;
    t)
        threads="$OPTARG"
        ;;
    o)
        builddir="$OPTARG"
        ;;
    h)
        usage
        ;;
    \?) 
        echo "Invalid option -$OPTARG" >&2
        ;;
  esac
done

echo "Building suk"
cd suk
mkdir -p build
cd build
if [ "$native" = true ] ; then
    cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
else
    cmake -DCMAKE_BUILD_TYPE=Release ..
fi    
make -j $threads
make
cd ../..

echo "Building overlap"
cd overlap
./install_deps.sh $threads
mkdir -p build
cd build
if [ "$native" = true ] ; then
    cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
else
    cmake -DCMAKE_BUILD_TYPE=Release ..
fi    
make -j $threads
make
cd ../..

echo "Building hypo polisher"
cd polisher
./install_deps.sh $threads
mkdir -p build
cd build
if [ "$native" = true ] ; then
    cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
else
    cmake -DCMAKE_BUILD_TYPE=Release ..
fi    
make -j $threads
make
cd ../..

echo "Building scaffolder"
cd scaffold
./install_deps.sh $threads
mkdir -p build
cd build
if [ "$native" = true ] ; then
    cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
else
    cmake -DCMAKE_BUILD_TYPE=Release ..
fi    
make -j $threads
make
cd ../..

mkdir -p $builddir
cp run_all.sh $builddir/
cp scan_misjoin.py $builddir/
cp suk/build/bin/suk $builddir/
cp overlap/run_overlap.sh $builddir/
cp overlap/join_overlap.py $builddir/
cp overlap/filter_overlap.py $builddir/
cp overlap/filter_alignments.py $builddir/
cp overlap/build/find_overlap $builddir/
cp polisher/build/bin/hypo $builddir/
cp scaffold/run_scaffold.sh $builddir/
cp scaffold/join_scaffold.py $builddir/
cp scaffold/filter_scaffold.py $builddir/
cp scaffold/get_necessary_reads.py $builddir/
cp scaffold/remove_duplicates.py $builddir/
cp scaffold/build/find_scaffold $builddir/

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
