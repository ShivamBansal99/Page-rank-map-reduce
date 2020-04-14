#include "mapreducempi.h"
using namespace std;

		template<class keytype,class valtype,class key2type,class val2type>
		mapreduce<keytype,valtype,key2type,val2type>::mapreduce(){
		}




		template<class keytype,class valtype,class key2type,class val2type>
		void mapreduce<keytype,valtype,key2type,val2type>::mapall(int nprocs,int me, void (* f)(int, map<keytype,vector<valtype> > *) ){
			f(me,&map1);
		}
		
		
		
		template<class keytype,class valtype,class key2type,class val2type>
		void mapreduce<keytype,valtype,key2type,val2type>::scatt(int me,int nprocs){
			vector<int> num_keys(nprocs);
			if(me==0){
				
				for(auto t:map1){
					int valsize = t.second.size();
					keytype key = t.first;
					MPI_Send(&(t.first),sizeof(keytype),MPI_CHAR,0,0,MPI_COMM_WORLD);
					for(int j=0;j<valsize;j++){
						valtype val = t.second[j]; 
						if(allmaps.find(key) == map1.end()){
							vector<double> v;
							allmaps[key] =v;
						}
						allmaps[key].push_back(val);
						allkeys.insert(key);
					}
				}
				
				for(int i=1;i<nprocs;i++){
					MPI_Recv(&num_keys[i-1],1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				for(int i=1;i<nprocs;i++){
					int num_vals;
					keytype key;
					for(int t=0;t<num_keys[i-1];t++){
						MPI_Recv(&key,sizeof(keytype),MPI_CHAR,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						MPI_Recv(&num_vals,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						for(int j=0;j<num_vals;j++){
							valtype val; 
							MPI_Recv(&val,sizeof(valtype),MPI_CHAR,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
							if(allmaps.find(key) == allmaps.end()){
								vector<double> v;
								allmaps[key] =v;
							}
							allmaps[key].push_back(val);
							allkeys.insert(key);
						}
					}
				}
				int each = allkeys.size()/(nprocs);
				int each2=each+1;
				int rem = allkeys.size()%(nprocs);
				auto itr=allkeys.begin();
				for(int r=0;r<each;r++){	
					keytype key = *itr;
					
					if(mymaps.find(key) == mymaps.end()){
						vector<double> v;
						mymaps[key] =v;
					}
					
					
					int nval = allmaps[*itr].size();
					for(int j=0;j<nval;j++){
						mymaps[key].push_back(allmaps[*itr][j]);
					}
					itr++;
				}
				if(rem>0){
					keytype key = *itr;
					
					if(mymaps.find(key) == mymaps.end()){
						vector<double> v;
						mymaps[key] =v;
					}
					
					
					int nval = allmaps[*itr].size();
					for(int j=0;j<nval;j++){
						mymaps[key].push_back(allmaps[*itr][j]);
					}
					itr++;
				}
				for(int i=1;i<nprocs && itr!=allkeys.end();i++){
					if(i>rem) MPI_Send(&each,1,MPI_INT,i,0,MPI_COMM_WORLD);
					else MPI_Send(&each2,1,MPI_INT,i,0,MPI_COMM_WORLD);
					for(int r=0;r<each;r++){	
						keytype k = *itr;
						MPI_Send(&k,sizeof(keytype),MPI_CHAR,i,0,MPI_COMM_WORLD);
						int nval = allmaps[*itr].size();
						MPI_Send(&nval,1,MPI_INT,i,0,MPI_COMM_WORLD);
						for(int j=0;j<nval;j++){
							MPI_Send(&(allmaps[*itr][j]),sizeof(valtype),MPI_CHAR,i,0,MPI_COMM_WORLD);
						}
						itr++;
					}
					if(i<rem){
						keytype k = *itr;
						MPI_Send(&k,sizeof(keytype),MPI_CHAR,i,0,MPI_COMM_WORLD);
						int nval = allmaps[*itr].size();
						MPI_Send(&nval,1,MPI_INT,i,0,MPI_COMM_WORLD);
						for(int j=0;j<nval;j++){
							MPI_Send(&(allmaps[*itr][j]),sizeof(valtype),MPI_CHAR,i,0,MPI_COMM_WORLD);
						}
						itr++;
					}
				}
				
			}else{
				int num_keys = map1.size();
				MPI_Send(&num_keys,1,MPI_INT,0,0,MPI_COMM_WORLD);
				
				for(auto t:map1){
					int valsize = t.second.size();
					MPI_Send(&(t.first),sizeof(keytype),MPI_CHAR,0,0,MPI_COMM_WORLD);
					MPI_Send(&(valsize),1,MPI_INT,0,0,MPI_COMM_WORLD);
					for(int j=0;j<valsize;j++){
						MPI_Send(&(t.second[j]),sizeof(valtype),MPI_CHAR,0,0,MPI_COMM_WORLD);
					}
				}
				int each;
				MPI_Recv(&each,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				for(int r=0;r<each;r++){
					keytype key;
					MPI_Recv(&key,sizeof(keytype),MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					int nval;
					if(mymaps.find(key) == mymaps.end()){
						vector<double> v;
						mymaps[key] =v;
					}
					MPI_Recv(&nval,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					for(int j=0;j<nval;j++){
						valtype v;
						MPI_Recv(&v,sizeof(valtype),MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						mymaps[key].push_back(v);
					}
				}
			}
		}
		
		
		
		template<class keytype,class valtype,class key2type,class val2type>
		void mapreduce<keytype,valtype,key2type,val2type>::reduce(int nprocs,int me, void (* f)(int,keytype, vector<valtype>, map<key2type,vector<val2type> > *) ){
			for(auto i:mymaps){
				f(i.second.size(),i.first,i.second, &map2);
			}
		}
		
		
		
		
		template<class keytype,class valtype,class key2type,class val2type>
		map<key2type, vector<val2type> > mapreduce<keytype,valtype,key2type,val2type>::get_maps(){
			return map2;
		}