cmake \
-DCMAKE_C_COMPILER=gcc \
-DCMAKE_CXX_COMPILER=g++ \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_BUILD_TYPE=Debug \
-DPHASTA_INCOMPRESSIBLE=OFF \
-DPHASTA_COMPRESSIBLE=ON \
-DLESLIB=/path/to/libles.a \
-DCASES=/path/to/phastaCases/ \
-DPHASTA_TESTING=ON \
-DPHASTA_USE_LESLIB=ON \
..

make

ctest