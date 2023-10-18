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
	vector <string> summstat;
	vector <int> n;
	vector <string> output;
	string ref;
	string vg;
	string ve; 
	string block;
	string plink;
	string c_t;
	double mafMax;
	int target;
	string r2;
	string pval;
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