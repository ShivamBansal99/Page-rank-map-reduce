#include<openmpi.h>
using namespace std;
#ifndef MAPREDUCEMPI_H
#define MAPREDUCEMPI_H
namespace mymapreduce{
	class mapreduce{
		public:
			mapreduce();
			void map();
			void scatt();
			void reduce();
		private:
			
	};
};

#endif