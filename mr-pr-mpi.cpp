#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <math.h>
#include <string>
#include <cstring>
#include <limits>
#include "mapreducempi.h"
using namespace std;

const char *OUTPUT_ARG = "-o";

double alpha =0.85; // the pagerank damping factor
double convergence=0.001;
unsigned long max_iterations=200;
vector<size_t> num_outgoing; // number of outgoing links per column
vector< vector<size_t> > rows; // the rowns of the hyperlink matrix
vector< vector<size_t> > columns; // the columns of the hyperlink matrix
vector<double> pr; // the pagerank table
vector<double> old_pr;
int me,nprocs;

void reset() {
    num_outgoing.clear();
    rows.clear();
    pr.clear();
}
	

const void error(const char *p,const char *p2) {
    cerr << p <<  ' ' << p2 <<  '\n';
    exit(1);
}

template <class Vector, class T>
bool insert_into_vector(Vector& v, const T& t) {
    typename Vector::iterator i = lower_bound(v.begin(), v.end(), t);
    if (i == v.end() || t < *i) {
        v.insert(i, t);
        return true;
    } else {
        return false;
    }
}


bool add_arc(size_t from, size_t to) {

    bool ret = false;
    size_t max_dim = max(from, to);
    if (rows.size() <= max_dim) {
        max_dim = max_dim + 1;
        rows.resize(max_dim);
		columns.resize(max_dim);
        if (num_outgoing.size() <= max_dim) {
            num_outgoing.resize(max_dim);
        }
    }

    ret = insert_into_vector(rows[to], from);
	insert_into_vector(columns[from], to);

    if (ret) {
        num_outgoing[from]++;
    }

    return ret;
}

int read_file(const string &filename) {
    reset();
    ifstream infile(filename);
	if (!infile) {
	  error("Cannot open file", filename.c_str());
	}
	int from_idx, to_idx; // indices of from and to nodes
    while (infile>>from_idx) {
		infile>>to_idx;
		add_arc(from_idx, to_idx);
	}
    return 0;
}



void mymap(int itask, map<int,vector<double> > * map1)
{
	int num_rows=rows.size();
	int step = num_rows/nprocs;
	for (int m = (itask)*step; m <= (itask)*step + step && m<num_rows; m++) {
		for(int j=0;j<columns[m].size();j++){
			double val = (num_outgoing[m])
				? (1.0*old_pr[m]) / (1.0*num_outgoing[m])
				: 0.0;
			// cerr<<num_outgoing[m]<<' '<<val<<endl;
			// if(columns[m][j]>20000) cerr<<"errorssss"<<endl;
			if(map1->find(columns[m][j]) == map1->end()){
				vector<double> v;
				(*map1)[columns[m][j]] =v;
			}
			(*map1)[columns[m][j]].push_back(val);	
		}
	}
}

void myreduce(int nvalues,int key, vector<double> vals,map<int,vector<double> > * map2) {
	// cerr<<"lol"<<endl;
	double h=0.0;
	assert(nvalues==vals.size());
	for (int i = 0; i <nvalues; i++) {
		/* The current element of the H vector */
		// cerr<<"here "<<(*((double *)(multivalue)))<<endl;
		h += vals[i];
	}
	h*=alpha;
	// cerr<<"key+h "<<(*((size_t *)(key)))<<' '<<h<<endl;
	pr[key] = h;
	if(map2->find(key) == map2->end()){
		vector<double> v;
		(*map2)[key] =v;
	}
	(*map2)[key].push_back(h);	
}

void pagerank() {

    vector<size_t>::iterator ci; // current incoming
    double diff = 1;
    size_t i;
    double sum_pr; // sum of current pagerank vector elements
    double dangling_pr; // sum of current pagerank vector elements for dangling
    			// nodes
    unsigned long num_iterations = 0;

    size_t num_rows = rows.size();
    
    if (num_rows == 0) {
        return;
    }
    
    pr.resize(num_rows);

    pr[0] = 1;

    
    while (diff > convergence && num_iterations < max_iterations) {
		
		// cout<<num_iterations<<' ' << max_iterations<<endl;
        sum_pr = 0;
        dangling_pr = 0;
        
        for (size_t k = 0; k < pr.size(); k++) {
            double cpr = pr[k];
            sum_pr += cpr;
            if (num_outgoing[k] == 0) {
                dangling_pr += cpr;
            }
        }

        if (num_iterations == 0) {
            old_pr = pr;
        } else {
            /* Normalize so that we start with sum equal to one */
            for (i = 0; i < pr.size(); i++) {
                old_pr[i] = pr[i] / sum_pr;
            }
        }

        /*
         * After normalisation the elements of the pagerank vector sum
         * to one
         */
        sum_pr = 1;
        
        /* An element of the A x I vector; all elements are identical */
        double one_Av = alpha * dangling_pr / num_rows;

        /* An element of the 1 x I vector; all elements are identical */
        double one_Iv = (1 - alpha) * sum_pr / num_rows;

        /* The difference to be checked for convergence */
        diff = 0;
		
	

		for (i = 0; i < num_rows; i++) {
			pr[i] = 0;
		}
		
		mapreduce<int,double,int, double> *mr = new mapreduce<int,double,int, double>();
		MPI_Barrier(MPI_COMM_WORLD);
		
		// cerr<<"here1 "<<nprocs<<endl;
		mr->mapall(nprocs,me,&mymap);
		// cerr<<"here2"<<endl;
		// MPI_Barrier(MPI_COMM_WORLD);
		mr->scatt(me,nprocs);
		mr->reduce(nprocs,me,&myreduce);
		
		MPI_Barrier(MPI_COMM_WORLD);
		map<int, vector<double> >mp = mr->get_maps(me,nprocs);
        
		for (i = 0; i < num_rows; i++) {
			if(mp[i].size()>0)
				pr[i]=mp[i][0];
			pr[i] += one_Av + one_Iv;
			diff += fabs(pr[i] - old_pr[i]);
        }
        //diff/=num_rows;
		num_iterations++;
    }
    
}

const void print_params(ostream& out) {
    out << "alpha = " << alpha << " convergence = " << convergence
        << " max_iterations = " << max_iterations
        << endl;
}

const void print_pagerank_v(ofstream& out) {

    size_t i;
    size_t num_rows = pr.size();
    double sum = 0;
    
    out.precision(numeric_limits<double>::digits10);

    for (i = 0; i < num_rows; i++) {
		out << i << " = " << pr[i] << endl;
        sum += pr[i];
    }
    cerr << "s = " << sum << " " << endl;
}

int check_inc(int i, int max) {
    if (i == max) {
        exit(1);
    }
    return i + 1;
}

int main(int argc, char **argv) {

	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&me);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    string input="";
	string output="";
	int i = 1;
	while (i < argc) {
		if (!strcmp(argv[i], OUTPUT_ARG)) {
			i = check_inc(i, argc);
			output = argv[i];
		} else{
			input = argv[i];
		}
		i++;
	}
    if(me==0){
		cerr<<nprocs<<endl;
		if(output.empty() || input.empty()){
			cerr<<"To run this file(assumed to be \"a\") use \n \t a <space> input_file <space> -o <space> output_file";
			MPI_Abort(MPI_COMM_WORLD,1);
			return 1;
		}
		print_params(cerr);
		cerr << "Reading input from " << input << "..." << endl;
	}
    if (!strcmp(input.c_str(), "stdin")) {
            read_file("");
    } else {
        read_file(input);
    }
    if(me==0) cerr << "Calculating pagerank..." << endl;
    pagerank();
    if(me==0){
		cerr << "Done calculating!" << endl;
		ofstream outfile(output);
		if (!outfile) {
		  error("Cannot open file", output.c_str());
		}
		print_pagerank_v(outfile);
	}
	MPI_Finalize();
}

