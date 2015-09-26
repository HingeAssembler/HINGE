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


class LAlignment { // because class Alignment is taken
	
public:
	LAlignment() {};
	std::string aseq;
	std::string bseq;
	int alen;
	int blen;
	void     *trace;
	int       tlen;
	int       diffs;
	int       abpos, bbpos;
	int       aepos, bepos;
	int flags;	
	
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
	
    FILE   *input;
    int64   novl;
    int     tspace, tbytes, small;
    int     reps, *pts;
    int     input_pts;
	

	LAInterface() {};
	int OpenDB2(std::string filename);
	int OpenDB(std::string filename);
	int OpenAlignment(std::string filename);
	void showRead(int from,int to);
	void showAlignment(int from, int to);
	Read * getRead(int number);
	void getRead(std::vector<Read *> reads, int from, int to);	
	void getAlignment(std::vector<int> &, int n);
	
	int CloseDB();
	int CloseDB2();
	
};

#endif