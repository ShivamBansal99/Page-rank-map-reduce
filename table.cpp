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
using namespace std;

const char *OUTPUT_ARG = "-o";

double alpha =0.85; // the pagerank damping factor
double convergence=0.000001;
unsigned long max_iterations=5;
vector<size_t> num_outgoing; // number of outgoing links per column
vector< vector<size_t> > rows; // the rowns of the hyperlink matrix
vector<double> pr; // the pagerank table


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
        if (num_outgoing.size() <= max_dim) {
            num_outgoing.resize(max_dim);
        }
    }

    ret = insert_into_vector(rows[to], from);

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

void pagerank() {

    vector<size_t>::iterator ci; // current incoming
    double diff = 1;
    size_t i;
    double sum_pr; // sum of current pagerank vector elements
    double dangling_pr; // sum of current pagerank vector elements for dangling
    			// nodes
    unsigned long num_iterations = 0;
    vector<double> old_pr;

    size_t num_rows = rows.size();
    
    if (num_rows == 0) {
        return;
    }
    
    pr.resize(num_rows);

    pr[0] = 1;

    
    while (diff > convergence && num_iterations < max_iterations) {

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
            /* The corresponding element of the H multiplication */
            double h = 0.0;
            for (ci = rows[i].begin(); ci != rows[i].end(); ci++) {
                /* The current element of the H vector */
                double h_v = (num_outgoing[*ci])
                    ? 1.0 / num_outgoing[*ci]
                    : 0.0;
                h += h_v * old_pr[*ci];
            }
            h *= alpha;
            pr[i] = h + one_Av + one_Iv;
            diff += fabs(pr[i] - old_pr[i]);
        }
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
	if(output.empty() || input.empty()){
		cerr<<"To run this file(assumed to be \"a\") use \n \t a <space> input_file <space> -o <space> output_file";
		return 1;
	}
    print_params(cerr);
    cerr << "Reading input from " << input << "..." << endl;
    if (!strcmp(input.c_str(), "stdin")) {
            read_file("");
    } else {
        read_file(input);
    }
    cerr << "Calculating pagerank..." << endl;
    pagerank();
    cerr << "Done calculating!" << endl;
	ofstream outfile(output); 
	if (!outfile) {
	  error("Cannot open file", output.c_str());
	}
    print_pagerank_v(outfile);
}

