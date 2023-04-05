/*
Multi-trait assisted Polygenic Scores (mtPGS)
*/

#ifndef __MVLMM_H__
#define __MVLMM_H__

#include <vector>
#include <string>
#include <iostream>

#include "mtPGS_function.hpp"

class PARAM {
public:
	//parameters
	vector <string> summstat_int;
	vector <string> summstat_ext;
	vector <int> n_ext;
	vector <string> output;
	vector <double> h2;
	vector <double> gencov;
	vector <double> envcov;
	int n_s;
	string ref;
	string vg;
	string ve;
	double mafMax; 
	string block;
	int t;
	int target;
};

class mtPGS {
public:
	//parameters
	string version;
	string date;
	string year;

	//constructor
	mtPGS(void);

	//functions
	void printHeader(void);
	void printHelp(void);
	void Assign(int argc, char ** argv, PARAM &cPar);
	void BatchRun(PARAM &cPar);
};


#endif