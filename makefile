g++ mr-pr-cpp.cpp /usr/lib/x86_64-linux-gnu/libboost_system.a /usr/lib/x86_64-linux-gnu/libboost_iostreams.a /usr/lib/x86_64-linux-gnu/libboost_filesystem.a -pthread -o PageRank.o -std=c++11-```



mpic++ -std=c++11 -m64 -g -O -I./src  -c mr-pr-mpi-base.cpp -o mr-pr-mpi-base.tem
mpic++ -g -O mr-pr-mpi-base.tem ./src/libmrmpi_mpicc.a  -o mr-pr-mpi-base.o