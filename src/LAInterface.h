#ifndef LAINTERFACE
#define LAINTERFACE

#include <vector>
#include <iostream>
#include <string>

extern "C" {
#include "DB.h"
#include "align.h"
}

class Read {
public:
	int id;
	std::string name;
	std::string bases;
	std::string qv;
	Read(int id, std::string name,std::string bases):id(id),bases(bases),name(name) {};
	void showRead();
};




class LAInterface {
public:
	
	HITS_DB   _db1, *db1 = &_db1; 
	HITS_DB   _db2, *db2 = &_db2; 
	Overlap   _ovl, *ovl = &_ovl;
	Alignment _aln, *aln = &_aln;
	char ** flist = NULL;
	int * findx = NULL; 
	int nfiles = 0;
	LAInterface() {};
	int OpenDB2(std::string filename);
	int OpenDB(std::string filename);
	void showRead(int from,int to);
	Read * getRead(int number);
	int CloseDB();
	int CloseDB2();
	
	
};

#endif