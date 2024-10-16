#!/bin/bash
#SBATCH -N1 -n10
#SBATCH -J SmithWaterman_zzb

export OMP_NUM_THREADS=10

./build/SmithWaterman ./data/tiny_q.fna ./data/tiny_t.fna ./data/tiny_ref.txt
./build/SmithWaterman ./data/small_q.fna ./data/small_t.fna ./data/small_ref.txt
./build/SmithWaterman ./data/medium_q.fna ./data/medium_t.fna ./data/medium_ref.txt