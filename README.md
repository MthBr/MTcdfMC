MTcdfMC_project
==============================


# The CFD Code Project
install on conda the cdfARE-environment.yml



# WellPlate Code Project

Code Project for the [cumulative distribution function solution for advection-reaction equations project].



# Code Project
The project is named Multi-dimensional delta-Dirac Approximations for Uncertainty Quantification of Flows through Randomly Heterogeneous Porous Media



## Code Development

you can use it as it basic version by
installing on conda the cdfARE-environment.yml

Then set up the kind of test in setup_testSet
Next execute each main script in every folder of the sequence:
data_etl, data_exp, feature_eng, model

Then you can investigate the various model 


## Code Development with AMG

For scalabitly install CUDA, then AMG and then pyamgx


AMG for windows:
https://github.com/NVIDIA/AMGX/releases/tag/v2.3.0

cmake -DCMAKE_C_COMPILER="C:/MinGW/bin/gcc.exe" -DCMAKE_CXX_COMPILER="C:/MinGW/bin/g++.exe" -DCMAKE_CUDA_COMPILER="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.6/bin/nvcc.exe" -DCUDA_ARCH="80"  -DCMAKE_NO_MPI="TRUE" ../


or for linux
cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_BUILD_TYPE=Release \
    -DCUDA_ARCH="80" -DCMAKE_NO_MPI="TRUE" ../
    
    
    
    
