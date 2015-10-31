#ifndef LAINTERFACE
#define LAINTERFACE

#include <vector>
#include <iostream>
#include <string>

extern "C" {
#include "DB.h"
#include "align.h"
}
typedef std::pair<int,int> Interval;

class Read { // read class
public:
    int id; // id, start from 0
    std::string name; // read name
    std::string bases; // read bases
    std::string qv; // qv currently not available
    std::vector<Interval> intervals;
    int effective_start,effective_end;
    int len;
    Read(int id, int length, std::string name, std::string bases) : id(id), bases(bases), name(name), len(length) { };
    Read(int id, std::string name, std::string bases) : id(id), bases(bases), name(name) { };

    bool active = true;
    void showRead();
};

enum aligntype {
    FORWARD, BACKWARD, COVERING, COVERED, MISMATCH_LEFT, MISMATCH_RIGHT, UNDIFINED, MIDDLE // different type of alignment
/**
 * FORWARD: Alignment and extend to the right
 * BACKWARD: extend to the left
 * COVERING: read a covering read b
 * COVERED: read a covered by read b
 * MISMATCH_LEFT: read a has a chimeric section on the left, and read b align with the rest of read a and extend it to the left
 * MISMATCH_RIGHT: read a has a chimeric section on the right, read b align with the rest of read a and extend it to the right
 * UNDIFIEND: any other exceptions
**/

} ;

class LAlignment { // because class Alignment is taken

public:
    LAlignment() { };
    //std::string aseq;
    //std::string bseq;
    char * aseq;
    char * bseq;


	void show() {printf("%d %d %d [%d...%d] x [%d...%d] %d diffs\n",aid,bid,flags,abpos,aepos,bbpos,bepos,diffs); };
    int aid; // id of read a
    int bid; // id of read b
    int alen; // length of read a
    int blen; // length of read b
    int *trace; // trace
    uint16 *trace_pts;
    int trace_pts_len;
    int tlen;
    int diffs;
    int abpos, bbpos; // begin position of read a and b
    int aepos, bepos; // end position of read a and b
    int flags; // flag = 1 : 'c', flag = 0 : 'n'
    int tps;
    aligntype aln_type;
	bool active = true;
};

class LOverlap { // LOverlap is a simplified version of LAlignment, no trace
public:
    LOverlap() { };
	void show() {printf("%d %d %d [%d...%d]/%d x [%d...%d]/%d %d diffs, %d type\n",aid,bid,flags,abpos,aepos,alen,bbpos,bepos,blen,diffs,aln_type); };
    int aid, bid;
    int alen; // length of read a
    int blen; // length of read b
    int tlen;
    int diffs; //differences
    int abpos, bbpos; // starting position and ending position of alignment in read a
    int aepos, bepos; // starting position and ending position of alignment in read b
    int tps;
    int flags; //flags, reverse complement = 1, same direction = 0
    int aes,aee,bes,bee;
    aligntype aln_type = UNDIFINED;
    void addtype(int);
    static const int CHI_THRESHOLD = 500; // threshold for chimeric/adaptor at the begining
	bool active = true;

};


class LAInterface {
public:

    HITS_DB _db1, *db1 = &_db1; // data base 1
    HITS_DB _db2, *db2 = &_db2; // data base 2
    Overlap _ovl, *ovl = &_ovl; // overlaps
    Alignment _aln, *aln = &_aln; // alignments, those are data structures required to read the data base

    char **flist = NULL;
    int *findx = NULL;
    int nfiles = 0; // n blocks of the read database

    char ** flist2 = NULL;
    int *findx2 = NULL;
    int nfiles2 = 0; // n blocks of read database 2


    FILE *input;
    int64 novl;
    int tspace, tbytes, small;
    int reps, *pts;
    int input_pts;


    LAInterface() { };

    int OpenDB2(std::string filename, std::string filename2); // open 2 databases

    int OpenDB(std::string filename); // open database

    int OpenAlignment(std::string filename); // open .las Alignment file

    void showRead(int from, int to); // show reads in a range
	
    void showRead2(int from, int to); // show reads in a range

    void showAlignment(int from, int to); // show alignment with 'A read' in a range

    void showOverlap(int from, int to); // show alignment with 'A read' in a range

    void resetAlignment(); // rewind the file, need to be called every time before obtaining alignments

    Read *getRead(int number); //get one read

    Read *getRead2(int number); //get one read

    void getRead(std::vector<Read *> &reads, int from, int to); // get reads within a range

    void getRead2(std::vector<Read *> &reads, int from, int to); // get reads within a range


    void getAlignmentB(std::vector<int> &, int n); // get all b reads aligned with a read

    void getOverlap(std::vector<LOverlap *> &, int from, int to); // get overlap(simplified version of alignment) with a read in a range

    void getOverlap(std::vector<LOverlap *> &, int n);

    void getAlignment(std::vector<LAlignment *> &, int from, int to); // get alignment with 'A read' in a range

    void getAlignment(std::vector<LAlignment *> &, int n);

    int CloseDB(); // close database

    int getReadNumber(); // get total number of reads

    int getReadNumber2(); // get total number of reads

    int64 getAlignmentNumber(); // get total number of alignments

    int CloseDB2();

    int LPrint_Alignment(FILE *file, Alignment *align, Work_Data *ework,
                        int indent, int width, int border, int upper, int coord);

    int LPrint_Alignment_exp(FILE *file, LAlignment *align, Work_Data *ework,
                         int indent, int width, int border, int upper, int coord);


    int LCompute_Trace_PTS(Alignment *align, Work_Data *ework, int trace_spacing);


    int Lshow_Alignment_tgs(LAlignment *);

    int generate_consensus(std::vector<LAlignment *> &);


    std::pair<std::string, std::string> Lget_Alignment_tgs(LAlignment *alignment);
};

#endif