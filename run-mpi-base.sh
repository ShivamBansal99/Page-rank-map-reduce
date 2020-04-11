mpic++ -std=c++11 -m64 -g -O -I./src  -c mr-pr-mpi-base.cpp -o mr-pr-mpi-base.tem
mpic++ -std=c++11 -g -O mr-pr-mpi-base.tem ./src/libmrmpi_mpicc.a  -o mr-pr-mpi-base.o
