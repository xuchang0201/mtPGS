/*
Multi-trait assisted Polygenic Scores (mtPGS)
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "mtPGS_int_ext.hpp"

using namespace std;

int main(int argc, char *argv[])
{
	mtPGS cDB;
	PARAM cPar;

	if (argc <= 1) {
		cDB.printHeader();
		return EXIT_SUCCESS;
	}
	if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') {
		cDB.printHelp();
		return EXIT_SUCCESS;
	}
	cDB.printHeader();
	cDB.Assign(argc, argv, cPar);
	cDB.BatchRun(cPar);
	return EXIT_SUCCESS;
}
