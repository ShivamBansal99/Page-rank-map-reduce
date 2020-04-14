#include<mpi.h>
#include <vector>
#include<bits/stdc++.h>
using namespace std;
#ifndef MAPREDUCEMPI_H
#define MAPREDUCEMPI_H
template<class keytype,class valtype,class key2type,class val2type>
	class mapreduce{
		public:
			mapreduce();
			void mapall(int nprocs,int me, void (*f)(int, map<keytype,vector<valtype> > *) );
			void scatt(int me,int nprocs);
			void reduce(int nprocs,int me, void (* f)(int,keytype, vector<valtype>, map<key2type,vector<val2type> > *) );
			map<key2type, vector<val2type> > get_maps(int me,int nprocs);
		private:
				map<keytype,vector<valtype> > map1;
		
				map<keytype,vector<valtype> > mymaps;
				
				map<keytype,vector<valtype> > allmaps;
				unordered_set<keytype> allkeys;

				map<key2type, vector<valtype> > map2;

			
	};
#include "mapreducempi.cpp"
#endif