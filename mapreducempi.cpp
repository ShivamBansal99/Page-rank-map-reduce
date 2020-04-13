#include "mapreducempi.h"
using namespace std;

namespace mymapreduce{
	template<keytype,valtype,key2type,val2type>
	class mapreduce{
		
		map<keytype,vector<valtype> > map1;
		vector<key2type> allkeys;
		map<key2type, vector<valtype>> map2;
		
		map<keytype,vector<valtype> > mymaps;
		vector<keytype> mykeys;
		
		mapreduce(){
			
		}
		
		
		void map(int nprocs,int me, (void *) f(int, map<keytype,vector<valtype> > *) ){
			f(me,&map1);
		}
		void scatt(){
			
		}
		void reduce(int nprocs,int me, (void *) f(int,keytype, vector<valtype> *) ){
			for(int i=0;i<allkeys.size();i++){
				f(allkeys[i]);
			}
		}
	}
}