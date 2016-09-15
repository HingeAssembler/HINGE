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

enum MatchType {
    FORWARD, BACKWARD, ACOVERB, BCOVERA, UNDEFINED, INTERNAL, NOT_ACTIVE, COVERING,
	COVERED, MIDDLE, MISMATCH_LEFT, MISMATCH_RIGHT, FORWARD_INTERNAL, BACKWARD_INTERNAL // different type of alignment
/**
 * FORWARD: Alignment and extend to the right
 * BACKWARD: extend to the left
 * COVERING: read a covering read b
 * COVERED: read a covered by read b
 * MISMATCH_LEFT: read a has a chimeric section on the left, and read b align with the rest of read a and extend it to the left
 * MISMATCH_RIGHT: read a has a chimeric section on the right, read b align with the rest of read a and extend it to the right
 * UNDEFINED: any other exceptions
 * FORWARD_INTERNAL : forward on read A internal on B
 * BACKWARD_INTERNAL : reverse on read A internal on B
**/

} ;

class LAlignment { // because class Alignment is taken

public:
    LAlignment() { };
    //std::string aseq;
    //std::string bseq;
    char * aseq;
    char * bseq;

    bool recovered = false;

	void show() {printf("%d %d %d [%d...%d] x [%d...%d] %d diffs\n", read_A_id_, read_B_id_,flags,abpos,aepos,bbpos,bepos,diffs); };
    int read_A_id_; // id of read a
    int read_B_id_; // id of read b
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
    MatchType aln_type;
	bool active = true;
};

class LOverlap { // LOverlap is a simplified version of LAlignment, no trace
public:
    LOverlap() { };
	void show() {printf("%d %d %d [%d...%d]/%d x [%d...%d]/%d %d diffs, %d type\n", read_A_id_, read_B_id_,
						reverse_complement_match_,
						read_A_match_start_, read_A_match_end_, alen, read_B_match_start_, read_B_match_end_, blen, diffs,
						match_type_); };
    int read_A_id_, read_B_id_;
    int alen; // length of read a
    int blen; // length of read b
    int tlen;
    int diffs; //differences
    int read_A_match_start_, read_B_match_start_; // starting position and ending position of alignment in read a
    int read_A_match_end_, read_B_match_end_; // starting position and ending position of alignment in read b
    int eff_read_A_match_start_, eff_read_B_match_start_, eff_read_A_match_end_, eff_read_B_match_end_;
    int tps;
    int reverse_complement_match_; //reverse_complement_match_, reverse complement = 1, same direction = 0
    int eff_read_A_read_start_, eff_read_A_read_end_, eff_read_B_read_start_, eff_read_B_read_end_;
    MatchType match_type_ = UNDEFINED;
    void addtype(int max_overhang); //classify overlaps
    void AddTypesAsymmetric(int max_overhang, int min_overhang);
	int GetMatchingPosition(int pos_A);
    static const int CHI_THRESHOLD = 500; // threshold for chimeric/adaptor at the begining
	bool active = true;
    uint16 *trace_pts;
    int trace_pts_len;
    void trim_overlap();
	void TrimOverlapNaive();
    int eff_start_trace_point_index_, eff_end_trace_point_index_;
    int weight;
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

    int openDB2(std::string filename, std::string filename2); // open 2 databases

    int openDB(std::string filename); // open database

    int openAlignmentFile(std::string filename); // open .las Alignment file

    void showRead(int from, int to); // show reads in a range
	
    void showRead2(int from, int to); // show reads in a range

    void showAlignment(int from, int to); // show alignment with 'A read' in a range

    void showOverlap(int from, int to); // show alignment with 'A read' in a range

    void resetAlignment(); // rewind the file, need to be called every time before obtaining alignments

    Read *getRead(int number); //get one read

    Read *getRead2(int number); //get one read

    void getRead(std::vector<Read *> &reads, int from, int to); // get reads within a range
	
	void getQV(std::vector<std::vector<int> > & QV, int from, int to);

    void getRead2(std::vector<Read *> &reads, int from, int to); // get reads within a range


    void getAlignmentB(std::vector<int> &, int n); // get all b reads aligned with a read

    void getOverlap(std::vector<LOverlap *> &, int from, int64 to); // get overlap(simplified version of alignment) with a read in a range

    void getOverlapw(std::vector<LOverlap *> &, int from, int to); // get overlap(simplified version of alignment) with a read in a range

    void getOverlap(std::vector<LOverlap *> &, int n);

    void getAlignment(std::vector<LAlignment *> &, int from, int to); // get alignment with 'A read' in a range

    void getAlignment(std::vector<LAlignment *> &result_vec, std::vector<int> &range);

    void getAlignment(std::vector<LAlignment *> &, int n);

    int closeDB(); // close database

    int getReadNumber(); // get total number of reads

    int getReadNumber2(); // get total number of reads from database 2

    int64 getAlignmentNumber(); // get total number of alignments

    int closeDB2();

    int printAlignment(FILE *file, Alignment *align, Work_Data *ework,
                       int indent, int width, int border, int upper, int coord);

    int printAlignment_exp(FILE *file, LAlignment *align, Work_Data *ework,
                           int indent, int width, int border, int upper, int coord);


    int computeTracePTS(Alignment *align, Work_Data *ework, int trace_spacing);


    int showAlignmentTags(LAlignment *);

    int generateConsensus(std::vector<LAlignment *> &);

    int recoverAlignment(LAlignment *);

    std::vector<int> * getCoverage(std::vector<LOverlap *> alns);

    std::vector<int> * getCoverage(std::vector<LAlignment *> alns);

    std::pair<std::string, std::string> getAlignmentTags(LAlignment *alignment);

    std::vector<std::pair<int, int> > * lowCoverageRegions(std::vector<int> & cov, int min_cov);

    void profileCoverage(std::vector<LOverlap *> &alignments, std::vector<std::pair<int, int> > & coverage,int reso, int cutoff);

    void profileCoveragefine(std::vector<LOverlap *> &alignments, std::vector<std::pair<int, int> > & coverage,int reso, int cutoff, int est_coverage);

    void repeatDetect(std::vector<std::pair<int, int> > & coverage, std::vector<std::pair<int, int> > & repeat);

	int loadPAF(std::string filename, std::vector<LOverlap *> &);

    int loadFASTA(std::string filename, std::vector<Read *> & reads);

};

class Node {
public:
    int id;
    int strand;
    bool pseudo = false;
    Node(int id, int strand): id(id), strand(strand) {};
    Node() {};
    void show() { std::cout<<id; if (strand == 1) std::cout<<"\'";}
};



#endif