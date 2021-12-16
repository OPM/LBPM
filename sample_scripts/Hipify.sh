#!/bin/bash

echo "Hipify cuda file $1"

sed -i 's/cudaError_t/hipError_t/g' $1
sed -i 's/cudaGetLastError/hipGetLastError/g' $1
sed -i 's/cudaSuccess/hipSuccess/g' $1
sed -i 's/cudaGetErrorString/hipGetErrorString/g' $1
sed -i 's/CUDA error/hip error/g' $1
