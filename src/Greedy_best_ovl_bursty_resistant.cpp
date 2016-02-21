#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <unordered_map>
#include "DB.h"
#include "align.h"
#include "LAInterface.h"
#include "OverlapGraph.h"
#include <algorithm>
#include <fstream>

#include <iostream>
#include <set>
#include <omp.h>
#include "INIReader.h"
#include <tuple>


#define LAST_READ_SYMBOL  '$'


typedef std::tuple<Node, Node, int> Edge_w;

typedef std::pair<Node, Node> Edge_nw;


static int ORDER(const void *l, const void *r) {
    int x = *((int32 *) l);
    int y = *((int32 *) r);
    return (x - y);
}


std::ostream& operator<<(std::ostream& out, const MatchType value){
    static std::map<MatchType, std::string> strings;
    if (strings.size() == 0){
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(FORWARD);
        INSERT_ELEMENT(BACKWARD);
        INSERT_ELEMENT(MISMATCH_LEFT);
        INSERT_ELEMENT(MISMATCH_RIGHT);
        INSERT_ELEMENT(COVERED);
        INSERT_ELEMENT(COVERING);
        INSERT_ELEMENT(UNDEFINED);
        INSERT_ELEMENT(MIDDLE);
#undef INSERT_ELEMENT
    }

    return out << strings[value];
}



bool compare_overlap(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->read_A_match_end_ - ovl1->read_A_match_start_ + ovl1->read_B_match_end_ - ovl1->read_B_match_start_) > (ovl2->read_A_match_end_ - ovl2->read_A_match_start_ + ovl2->read_B_match_end_ - ovl2->read_B_match_start_));
}

bool compare_sum_overlaps(const std::vector<LOverlap * > * ovl1, const std::vector<LOverlap *> * ovl2) {
    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < ovl1->size(); i++) sum1 += (*ovl1)[i]->read_A_match_end_ - (*ovl1)[i]->read_A_match_start_ + (*ovl1)[i]->read_B_match_end_ - (*ovl1)[i]->read_B_match_start_;
    for (int i = 0; i < ovl2->size(); i++) sum2 += (*ovl2)[i]->read_A_match_end_ - (*ovl2)[i]->read_A_match_start_ + (*ovl2)[i]->read_B_match_end_ - (*ovl2)[i]->read_B_match_start_;
    return sum1 > sum2;
}

bool compare_pos(LOverlap * ovl1, LOverlap * ovl2) {
    return (ovl1->read_A_match_start_) > (ovl2->read_A_match_start_);
}

bool compare_overlap_abpos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->read_A_match_start_ < ovl2->read_A_match_start_;
}

bool compare_overlap_aepos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->read_A_match_start_ > ovl2->read_A_match_start_;
}

std::vector<std::pair<int,int>> Merge(std::vector<LOverlap *> & intervals, int cutoff)
{
    //std::cout<<"Merge"<<std::endl;
    std::vector<std::pair<int, int > > ret;
    int n = intervals.size();
    if (n == 0) return ret;

    if(n == 1) {
        ret.push_back(std::pair<int,int>(intervals[0]->read_A_match_start_, intervals[0]->read_A_match_end_));
        return ret;
    }

    sort(intervals.begin(),intervals.end(),compare_overlap_abpos); //sort according to left

    int left= intervals[0]->read_A_match_start_ + cutoff, right = intervals[0]->read_A_match_end_ - cutoff; //left, right means maximal possible interval now

    for(int i = 1; i < n; i++)
    {
        if(intervals[i]->read_A_match_start_ + cutoff <= right)
        {
            right=std::max(right, intervals[i]->read_A_match_end_ - cutoff);
        }
        else
        {
            ret.push_back(std::pair<int, int>(left,right));
            left = intervals[i]->read_A_match_start_ + cutoff;
            right = intervals[i]->read_A_match_end_ - cutoff;
        }
    }
    ret.push_back(std::pair<int, int>(left,right));
    return ret;
}

Interval Effective_length(std::vector<LOverlap *> & intervals, int min_cov) {
    Interval ret;
    sort(intervals.begin(),intervals.end(),compare_overlap_abpos); //sort according to left

    if (intervals.size() > min_cov) {
        ret.first = intervals[min_cov]->read_A_match_start_;
    } else
        ret.first = 0;
    sort(intervals.begin(),intervals.end(),compare_overlap_aepos); //sort according to left
    if (intervals.size() > min_cov) {
        ret.second = intervals[min_cov]->read_A_match_end_;
    } else
        ret.second = 0;
    return ret;
}


int main(int argc, char *argv[]) {


    LAInterface la;
	char * name_db = argv[1];
	char * name_las = argv[2];
    char * name_config = argv[4];
	printf("name of db: %s, name of .las file %s\n", name_db, name_las);
    la.openDB(name_db);
    std::cout<<"# Reads:" << la.getReadNumber() << std::endl;
    la.openAlignmentFile(name_las);
    std::cout<<"# Alignments:" << la.getAlignmentNumber() << std::endl;


	int n_aln = la.getAlignmentNumber();
	int n_read = la.getReadNumber();
    std::vector<LOverlap *> aln;
	la.resetAlignment();
    la.getOverlap(aln,0,n_aln);

    std::vector<Read *> reads;
    la.getRead(reads,0,n_read);

    std::cout << "input data finished" <<std::endl;


    /*
     * load config
     */

    INIReader reader(name_config);

    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<<name_config<<std::endl;
        return 1;
    }

    int LENGTH_THRESHOLD = reader.GetInteger("filter", "length_threshold", -1);
    double QUALITY_THRESHOLD = reader.GetReal("filter", "quality_threshold", 0.0);
    int N_ITER = reader.GetInteger("filter", "n_iter", -1);
    int ALN_THRESHOLD = reader.GetInteger("filter", "aln_threshold", -1);
    int MIN_COV = reader.GetInteger("filter", "min_cov", -1);
    int CUT_OFF = reader.GetInteger("filter", "cut_off", -1);
    int THETA = reader.GetInteger("filter", "theta", -1);
	int N_PROC = reader.GetInteger("running", "n_proc", 4);

    omp_set_num_threads(N_PROC);

    std::vector< std::vector<std::vector<LOverlap*>* > > idx2(n_read); //unordered_map from (aid) to alignments in a vector
    std::vector<Edge_w> edgelist; // save output to edgelist
    std::unordered_map<int, std::vector <LOverlap * > >idx3; // this is the pileup
    std::vector<std::set<int> > has_overlap(n_read);
    std::unordered_map<int, std::unordered_map<int, std::vector<LOverlap *> > > idx; //unordered_map from (aid, bid) to alignments in a vector
/*
 * Index alignments by:
 * 1) read A - idx3
 * 2) read A-read B - idx
 */

    for (int i = 0; i< n_read; i++) {
        has_overlap[i] = std::set<int>();
        idx3[i] = std::vector<LOverlap *>();
    }

    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx[aln[i]->read_A_id_][aln[i]->read_B_id_] = std::vector<LOverlap *>();
        }
    }


    for (int i = 0; i < aln.size(); i++) {
    if (aln[i]->active) {
	    has_overlap[aln[i]->read_A_id_].insert(aln[i]->read_B_id_);
        }
    }

    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx3[aln[i]->read_A_id_].push_back(aln[i]);
        }
    }


    std::cout<<"index data"<<std::endl;
    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx[aln[i]->read_A_id_][aln[i]->read_B_id_].push_back(aln[i]);
        }
    }
    std::cout<<"index data"<<std::endl;


    for (int i = 0; i < n_read; i++)
        for (std::set<int>::iterator j = has_overlap[i].begin(); j != has_overlap[i].end(); j++) {
            idx2[i].push_back(&(idx[i][*j]));
    }

    std::cout<<"index data finished"<<std::endl;
    std::unordered_map<int,std::vector<Interval> > covered_region;


    //for (int i = 0; i < n_read; i++) {
    //    printf("read %d [%d %d]/%d\n", i, reads[i]->effective_start, reads[i]->effective_end, reads[i]->len);
    //}

    /** Reads and alignments filtering
     *
     **/

    for (int n_iter = 0; n_iter < N_ITER; n_iter ++) {

#pragma omp parallel for
        for (int i = 0; i < n_read; i++) {
            if (reads[i]->active) {
                reads[i]->intervals = Merge(idx3[i], CUT_OFF);

                Interval cov = Effective_length(idx3[i], MIN_COV);
                reads[i]->effective_start = cov.first;
                reads[i]->effective_end = cov.second;

            }
        } // find all covered regions, could help remove adaptors

        std::cout<<"covered region"<<std::endl;
# pragma omp parallel for
        for (int i = 0; i < n_read; i++) {
            if (reads[i]->active)
            if ((reads[i]->effective_end - reads[i]->effective_start <
                 LENGTH_THRESHOLD) or (reads[i]->intervals.size() != 1))
                reads[i]->active = false;
        } // filter according to effective length, and interval size

        std::cout<<"filter data"<<std::endl;

        int num_active = 0;

        for (int i = 0; i < n_read; i++) {
            if (reads[i]->active)
                num_active++;
        }
        std::cout << "num active reads " << num_active << std::endl;
        for (int i = 0; i < n_aln; i++) {
            if (aln[i]->active)
            if ((not reads[aln[i]->read_A_id_]->active) or (not reads[aln[i]->read_B_id_]->active) or
                (aln[i]->diffs * 2 /  float(aln[i]->read_B_match_end_ - aln[i]->read_B_match_start_ + aln[i]->read_A_match_end_ - aln[i]->read_A_match_start_) >
                 QUALITY_THRESHOLD) or (aln[i]->read_A_match_end_ - aln[i]->read_A_match_start_ < ALN_THRESHOLD))
                aln[i]->active = false;
        }

# pragma omp parallel for
        for (int i = 0; i < n_aln; i++) {
            if (aln[i]->active) {
                aln[i]->eff_read_A_start_ = reads[aln[i]->read_A_id_]->effective_start;
                aln[i]->eff_read_A_end_ = reads[aln[i]->read_A_id_]->effective_end;
                
				if (aln[i]->reverse_complement_match_ == 0) {
					aln[i]->eff_read_B_start_ = reads[aln[i]->read_B_id_]->effective_start;
                	aln[i]->eff_read_B_end_ = reads[aln[i]->read_B_id_]->effective_end;
				}
				else {
					aln[i]->eff_read_B_start_ = aln[i]->blen - reads[aln[i]->read_B_id_]->effective_end;
                	aln[i]->eff_read_B_end_ = aln[i]->blen - reads[aln[i]->read_B_id_]->effective_start;
				}
            }
        }

        int num_active_aln = 0;
# pragma omp parallel for reduction(+:num_active_aln)
        for (int i = 0; i < n_aln; i++) {
            if (aln[i]->active) {
                num_active_aln++;
                aln[i]->addtype(THETA);
            }
        }
        std::cout << "num active alignments " << num_active_aln << std::endl;

        idx3.clear();
        for (int i = 0; i < aln.size(); i++) {
            if (aln[i]->active) {
                idx3[aln[i]->read_A_id_].push_back(aln[i]);
            }
        }
    }

    //sort the reads by sum of alignment length
# pragma omp parallel for
    for (int i = 0; i < n_read; i++ ) {
        std::sort( idx2[i].begin(), idx2[i].end(), compare_sum_overlaps );
    }

    /**
	 **    get the best overlap for each read and form a graph
	 **/

    for (int n_iter = 0; n_iter < 3; n_iter ++) {
        int no_edge = 0;
        edgelist.clear();
        for (int i = 0; i < n_read; i++) {
            int cf = 0;
            int cb = 0;
            /*
             * For each read, if there is exact right match (type FORWARD), choose the one with longest alignment with read A
             * same for BACKWARD,
             */

            if (reads[i]->active)
                for (int j = 0; j < idx2[i].size(); j++) {

                    if ((*idx2[i][j])[0]->active) {
                        for (int kk = 0; kk < idx2[i][j]->size(); kk++) {
                            if ((reads[(*idx2[i][j])[kk]->read_B_id_]->active) and ((*idx2[i][j])[kk]->active))
                            if (((*idx2[i][j])[kk]->match_type_ == FORWARD) and (cf < 1)) {
                                cf += 1;
                                //add edge
                                if ((*idx2[i][j])[kk]->reverse_complement_match_ == 1) { // n = 0, c = 1 n:this read itself, c:reverse complement
                                    edgelist.push_back(
                                            std::make_tuple(Node((*idx2[i][j])[kk]->read_A_id_, 0),
                                                                  Node((*idx2[i][j])[kk]->read_B_id_, 1), (*idx2[i][j])[kk]->read_A_match_end_ -(*idx2[i][j])[kk]->read_A_match_start_));
                                } else {
                                    edgelist.push_back(
                                            std::make_tuple(Node((*idx2[i][j])[kk]->read_A_id_, 0),
                                                                  Node((*idx2[i][j])[kk]->read_B_id_, 0), (*idx2[i][j])[kk]->read_A_match_end_ -(*idx2[i][j])[kk]->read_A_match_start_));
                                }
                            }
                            if ((reads[(*idx2[i][j])[kk]->read_B_id_]->active) and ((*idx2[i][j])[kk]->active))
                            if (((*idx2[i][j])[kk]->match_type_ == BACKWARD) and (cb < 1)) {
                                cb += 1;
                                //add edge
                                if ((*idx2[i][j])[kk]->reverse_complement_match_ == 1) {
                                    edgelist.push_back(
                                            std::make_tuple(Node((*idx2[i][j])[kk]->read_A_id_, 1),
                                                                  Node((*idx2[i][j])[kk]->read_B_id_, 0), (*idx2[i][j])[kk]->read_A_match_end_ -(*idx2[i][j])[kk]->read_A_match_start_));
                                } else {
                                    edgelist.push_back(
                                            std::make_tuple(Node((*idx2[i][j])[kk]->read_A_id_, 1),
                                                                  Node((*idx2[i][j])[kk]->read_B_id_, 1),(*idx2[i][j])[kk]->read_A_match_end_ -(*idx2[i][j])[kk]->read_A_match_start_));
                                }
                            }
                            if ((cf >= 1) and (cb >= 1)) break;
                        }
                    }
                    if ((cf >= 1) and (cb >= 1)) break;
                }

                if (reads[i]->active) {
                    if (cf == 0) printf("%d has no out-going edges\n", i);
                    if (cb == 0) printf("%d' has no out-going edges\n", i);
                    /*if ((cf == 0) or (cb == 0))
                        for (int j = 0; j < idx2[i].size(); j++) {
                            printf("%d,%d,%d,%d\n", i, idx2[i].size(), j, idx2[i][j].size());
                            for (int k = 0; k < idx2[i][j].size(); k++) {
                                std::cout << " " << idx2[i][j][k]->active << " " << "[" << idx2[i][j][k]->read_A_match_start_ << "," <<
                                idx2[i][j][k]->read_A_match_end_ << "]/" << idx2[i][j][k]->alen << " " << "[" << idx2[i][j][k]->read_B_match_start_ <<
                                "," << idx2[i][j][k]->read_B_match_end_ << "]/" << idx2[i][j][k]->blen << " " << idx2[i][j][k]->match_type_ <<
                                std::endl;
                            }
                        }
                        */
                    if ((cf == 0) or (cb ==0)) { // if no right extension or no left extension, throw away this read
                        no_edge++;
                        reads[i]->active = false;
                    }
                }

            }

            printf("no_edge nodes:%d\n",no_edge); // show number of reads with no extension

# pragma omp parallel for
            for (int ii = 0; ii < n_aln; ii++) {
            if (aln[ii]->active)
            if ((not reads[aln[ii]->read_A_id_]->active) or (not reads[aln[ii]->read_B_id_]->active))
                aln[ii]->active = false;
            }

            int num_active_aln = 0;
# pragma omp parallel for reduction(+:num_active_aln)
            for (int i = 0; i < n_aln; i++) {
                if (aln[i]->active) {
                    num_active_aln++;
                    aln[i]->addtype(THETA);
                }
            }
            std::cout << "num active alignments " << num_active_aln << std::endl;

    }

    std::ofstream out(argv[3], std::ofstream::out);
    //print edges list to file, in node, out node, weight (overlap size)
    for (int i = 0; i < edgelist.size(); i++){
        Node n1,n2;
        int w;

        n1 = std::get<0>(edgelist[i]);
        n2 = std::get<1>(edgelist[i]);
        w = std::get<2>(edgelist[i]);
        n1.show(out);
        out<<"->";
        n2.show(out);
        out << ",";
        out << w;
        out<<std::endl;
    }

    la.closeDB(); //close database
    return 0;
}