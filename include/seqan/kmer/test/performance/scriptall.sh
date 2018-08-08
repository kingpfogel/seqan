#!/bin/bash
cd /storage/mi/kingpfogel/buildall/
cmake /home/mi/kingpfogel/dev/seqan/include/seqan/kmer/test/ -DCMAKE_CXX_COMPILER=/home/mi/kingpfogel/miniconda3/bin/g++ -DCMAKE_PREFIX_PATH="$HOME/dev/seqan/util/cmake" -DSEQAN_INCLUDE_PATH="$HOME/dev/seqan/include" -DCMAKE_BUILD_TYPE=RELEASE -DSDSL_INCLUDE_DIRS="../../sdsl-lite/include/" -DCMAKE_CXX_FLAGS="-O3 -march=native -std=c++17"
cd /storage/mi/kingpfogel/buildall/performance/

