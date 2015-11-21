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


std::ostream& operator<<(std::ostream& out, const aligntype value){
    static std::map<aligntype, std::string> strings;
    if (strings.size() == 0){
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(FORWARD);
        INSERT_ELEMENT(BACKWARD);
        INSERT_ELEMENT(MISMATCH_LEFT);
        INSERT_ELEMENT(MISMATCH_RIGHT);
        INSERT_ELEMENT(COVERED);
        INSERT_ELEMENT(COVERING);
        INSERT_ELEMENT(UNDIFINED);
        INSERT_ELEMENT(MIDDLE);
#undef INSERT_ELEMENT
    }

    return out << strings[value];
}



bool compare_overlap(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->aepos - ovl1->abpos + ovl1->bepos - ovl1->bbpos) > (ovl2->aepos - ovl2->abpos + ovl2->bepos - ovl2->bbpos));
}

bool compare_sum_overlaps(const std::vector<LOverlap * > * ovl1, const std::vector<LOverlap *> * ovl2) {
    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < ovl1->size(); i++) sum1 += (*ovl1)[i]->aepos - (*ovl1)[i]->abpos + (*ovl1)[i]->bepos - (*ovl1)[i]->bbpos;
    for (int i = 0; i < ovl2->size(); i++) sum2 += (*ovl2)[i]->aepos - (*ovl2)[i]->abpos + (*ovl2)[i]->bepos - (*ovl2)[i]->bbpos;
    return sum1 > sum2;
}

bool compare_pos(LOverlap * ovl1, LOverlap * ovl2) {
    return (ovl1->abpos) > (ovl2->abpos);
}

bool compare_overlap_abpos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->abpos < ovl2->abpos;
}

bool compare_overlap_aepos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->abpos > ovl2->abpos;
}

std::vector<std::pair<int,int>> Merge(std::vector<LOverlap *> & intervals, int cutoff)
{
    //std::cout<<"Merge"<<std::endl;
    std::vector<std::pair<int, int > > ret;
    int n = intervals.size();
    if (n == 0) return ret;

    if(n == 1) {
        ret.push_back(std::pair<int,int>(intervals[0]->abpos,intervals[0]->aepos));
        return ret;
    }

    sort(intervals.begin(),intervals.end(),compare_overlap_abpos); //sort according to left

    int left=intervals[0]->abpos + cutoff, right = intervals[0]->aepos - cutoff; //left, right means maximal possible interval now

    for(int i = 1; i < n; i++)
    {
        if(intervals[i]->abpos + cutoff <= right)
        {
            right=std::max(right,intervals[i]->aepos - cutoff);
        }
        else
        {
            ret.push_back(std::pair<int, int>(left,right));
            left = intervals[i]->abpos + cutoff;
            right = intervals[i]->aepos - cutoff;
        }
    }
    ret.push_back(std::pair<int, int>(left,right));
    return ret;
}

Interval Effective_length(std::vector<LOverlap *> & intervals, int min_cov) {
    Interval ret;
    sort(intervals.begin(),intervals.end(),compare_overlap_abpos); //sort according to left

    if (intervals.size() > min_cov) {
        ret.first = intervals[min_cov]->abpos;
    } else
        ret.first = 0;
    sort(intervals.begin(),intervals.end(),compare_overlap_aepos); //sort according to left
    if (intervals.size() > min_cov) {
        ret.second = intervals[min_cov]->aepos;
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


	std::vector<std::vector<int>>  QV;
	la.getQV(QV,0,n_read);

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
    std::vector<Edge_w> edgelist, edgelist_ms; // save output to edgelist
    std::unordered_map<int, std::vector <LOverlap * > >idx3,idx4; // this is the pileup
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
        idx4[i] = std::vector<LOverlap *>();
    }

    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx[aln[i]->aid][aln[i]->bid] = std::vector<LOverlap *>();
        }
    }


    for (int i = 0; i < aln.size(); i++) {
    if (aln[i]->active) {
	    has_overlap[aln[i]->aid].insert(aln[i]->bid);
        }
    }

    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx3[aln[i]->aid].push_back(aln[i]);
			idx4[aln[i]->aid].push_back(aln[i]);
        }
    }

	std::cout<<"profile coverage" << std::endl;
    std::ofstream cov("coverage.txt");
for (int i = 0; i < n_read; i ++) {
    std::vector<std::pair<int, int> > coverage, repeat;
    la.profileCoverage(idx3[i],coverage,40);
    cov << "read " << i <<" ";
    for (int j = 0; j < coverage.size(); j++)
        cov << coverage[j].first << ","  << coverage[j].second << " ";
    cov << std::endl;
}

    std::cout<<"index data"<<std::endl;
    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx[aln[i]->aid][aln[i]->bid].push_back(aln[i]);
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

    /**Filtering
     *
     **/

    for (int n_iter = 0; n_iter < N_ITER; n_iter ++) {

		//#pragma omp parallel for
        for (int i = 0; i < n_read; i++) {
            if (reads[i]->active) {
                reads[i]->intervals = Merge(idx3[i], CUT_OFF);

                Interval cov = Effective_length(idx3[i], MIN_COV);
                reads[i]->effective_start = cov.first;
                reads[i]->effective_end = cov.second;
            }
        } // find all covered regions, could help remove adaptors

        std::cout<<"covered region"<<std::endl;
		//# pragma omp parallel for
        for (int i = 0; i < n_read; i++) {
            if (reads[i]->active)
			{
				//printf("read %d, el %d, is%d\n", i,reads[i]->effective_end - reads[i]->effective_start,reads[i]->intervals.size()  );
            	if ((reads[i]->effective_end - reads[i]->effective_start <
            	     LENGTH_THRESHOLD) or (reads[i]->intervals.size() != 1))
            	    reads[i]->active = false;
			}
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
            if ((not reads[aln[i]->aid]->active) or (not reads[aln[i]->bid]->active) or
                (aln[i]->diffs * 2 /  float(aln[i]->bepos - aln[i]->bbpos + aln[i]->aepos - aln[i]->abpos) >
                 QUALITY_THRESHOLD) or (aln[i]->aepos - aln[i]->abpos < ALN_THRESHOLD))
                aln[i]->active = false;
        }

		//# pragma omp parallel for
        for (int i = 0; i < n_aln; i++) {
            if (aln[i]->active) {
                aln[i]->aes = reads[aln[i]->aid]->effective_start;
                aln[i]->aee = reads[aln[i]->aid]->effective_end;

				if (aln[i]->flags == 0) {
					aln[i]->bes = reads[aln[i]->bid]->effective_start;
                	aln[i]->bee = reads[aln[i]->bid]->effective_end;
				}
				else {
					aln[i]->bes = aln[i]->blen - reads[aln[i]->bid]->effective_end;
                	aln[i]->bee = aln[i]->blen - reads[aln[i]->bid]->effective_start;
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
                idx3[aln[i]->aid].push_back(aln[i]);
            }
        }
    }





    //sort the reads
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
                            if ((reads[(*idx2[i][j])[kk]->bid]->active) and ((*idx2[i][j])[kk]->active))
                            if (((*idx2[i][j])[kk]->aln_type == FORWARD) and (cf < 1)) {
                                cf += 1;
                                //add edge
                                if ((*idx2[i][j])[kk]->flags == 1) { // n = 0, c = 1
                                    edgelist.push_back(
                                            std::make_tuple(Node((*idx2[i][j])[kk]->aid, 0),
                                                                  Node((*idx2[i][j])[kk]->bid, 1), (*idx2[i][j])[kk]->aepos -(*idx2[i][j])[kk]->abpos));
                                } else {
                                    edgelist.push_back(
                                            std::make_tuple(Node((*idx2[i][j])[kk]->aid, 0),
                                                                  Node((*idx2[i][j])[kk]->bid, 0), (*idx2[i][j])[kk]->aepos -(*idx2[i][j])[kk]->abpos));
                                }
                            }
                            if ((reads[(*idx2[i][j])[kk]->bid]->active) and ((*idx2[i][j])[kk]->active))
                            if (((*idx2[i][j])[kk]->aln_type == BACKWARD) and (cb < 1)) {
                                cb += 1;
                                //add edge
                                if ((*idx2[i][j])[kk]->flags == 1) {
                                    edgelist.push_back(
                                            std::make_tuple(Node((*idx2[i][j])[kk]->aid, 1),
                                                                  Node((*idx2[i][j])[kk]->bid, 0), (*idx2[i][j])[kk]->aepos -(*idx2[i][j])[kk]->abpos));
                                } else {
                                    edgelist.push_back(
                                            std::make_tuple(Node((*idx2[i][j])[kk]->aid, 1),
                                                                  Node((*idx2[i][j])[kk]->bid, 1),(*idx2[i][j])[kk]->aepos -(*idx2[i][j])[kk]->abpos));
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
                                std::cout << " " << idx2[i][j][k]->active << " " << "[" << idx2[i][j][k]->abpos << "," <<
                                idx2[i][j][k]->aepos << "]/" << idx2[i][j][k]->alen << " " << "[" << idx2[i][j][k]->bbpos <<
                                "," << idx2[i][j][k]->bepos << "]/" << idx2[i][j][k]->blen << " " << idx2[i][j][k]->aln_type <<
                                std::endl;
                            }
                        }
                        */
                    if ((cf == 0) or (cb ==0)) {
                        no_edge++;
                        reads[i]->active = false;
                    }
                }

            }

            printf("no_edge nodes:%d\n",no_edge);

# pragma omp parallel for
            for (int ii = 0; ii < n_aln; ii++) {
            if (aln[ii]->active)
            if ((not reads[aln[ii]->aid]->active) or (not reads[aln[ii]->bid]->active))
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

/** look for missing edges **/


    for (int i = 0; i < edgelist.size(); i++){
        Node n1,n2;
        int w;
        n1 = std::get<0>(edgelist[i]);
        n2 = std::get<1>(edgelist[i]);
        w = std::get<2>(edgelist[i]);
	    std::vector<std::pair<int, int> > coverage, coverage1;
	    la.profileCoverage(idx4[n2.id],coverage,40);
		int es = reads[n2.id]->effective_start;
		int ee = reads[n2.id]->effective_end;
		int es1 = reads[n1.id]->effective_start;
		int ee1 = reads[n1.id]->effective_end;
		int cov_start = 0;
		int cov_end = 0;
		int cov_start1 = 0;
		int cov_end1 = 0;

	    la.profileCoverage(idx4[n1.id],coverage1,40); // profile coverage and look for unusual pattern
		for (int j = 10; j < 20; j++) {
			cov_start1 	+= coverage1[es1/40+j].second;
			cov_end1 	+= coverage1[ee1/40-j].second;
		}

		cov_start1/=10;
		cov_end1/=10;


		for (int j = 10; j < 20; j++) {
			cov_start 	+= coverage[es/40+j].second;
			cov_end 	+= coverage[ee/40-j].second;
		}
		cov_start/=10;
		cov_end/=10;
		if (((cov_start/cov_end>=2) and (cov_start>100) and (n2.strand == 0))
			or ((cov_end /cov_start >= 2) and (cov_end > 100) and (n2.strand == 1)))
				{
					printf("read %d %d cov_start %d, cov_end %d\n", n2.id, n2.strand, cov_start, cov_end);
					//printf("previous %d %d cov_start %d, cov_end %d\n",n1.id, n1.strand, cov_start1, cov_end1);

                    int repeat_end = 0;
                    std::set<int> candidate; // look for candidates
					bool rep = false;
                    
					if (n2.strand == 0) {
						
						for (int j = 1; j < coverage.size(); j++) {
                	        if (coverage[j].second < 0.5 * coverage[j - 1].second) {
                	            printf("  -sharp drop at %d\n", coverage[j].first); // debug output
                	            repeat_end = coverage[j].first;
								rep = true;
                	            break;
                	        }
                	    }
                		if (rep)
                	    for (int k = 0; k < idx3[n2.id].size(); k++) {
                	        LOverlap * ovl_short = idx3[n2.id][k];
                	        if (ovl_short->aepos < repeat_end + 300) {
                	            candidate.insert(ovl_short->bid);
                	        }
                	    }
					} else {
						
						for (int j = coverage.size() - 1; j > 0 ; j--) {
                	        if (coverage[j].second < 0.5 * coverage[j + 1].second) {
                	            printf("  -sharp increase at %d\n", coverage[j].first); // debug output,
                	            repeat_end = coverage[j].first;
								rep = true;
                	            break;
                	        }
                	    }
                	
						if (rep)
                	    for (int k = 0; k < idx3[n2.id].size(); k++) {
                	        LOverlap * ovl_short = idx3[n2.id][k];
                	        if (ovl_short->abpos > repeat_end - 300) {
                	            candidate.insert(ovl_short->bid);
                	        }
                	    }
						
					}
                    /**
                     * Repeat the previous edge looking procedure, look for longest one, except that this time look for read in the candidate set
                     */

                    int cb = 0, cf = 0;
                        int current = n1.id;

                        for (int jj = 0; jj < idx2[current].size(); jj++) {
                            if ((*idx2[current][jj])[0]->active) {
                                for (int kk = 0; kk < idx2[current][jj]->size(); kk++) {
                                    if ((reads[(*idx2[current][jj])[kk]->bid]->active) and ((*idx2[current][jj])[kk]->active))
                                    if (n1.strand == 0)
                                    if (((*idx2[current][jj])[kk]->aln_type == FORWARD) and (cf < 1) and (candidate.find((*idx2[current][jj])[kk]->bid) != candidate.end() ))  {
                                        cf += 1;
                                        printf("add edge\n");
                                        //add edge
                                        if ((*idx2[current][jj])[kk]->flags == 1) { // n = 0, c = 1
                                            edgelist_ms.push_back(
                                                    std::make_tuple(Node((*idx2[current][jj])[kk]->aid, 0),
                                                                    Node((*idx2[current][jj])[kk]->bid, 1), (*idx2[current][jj])[kk]->aepos -(*idx2[current][jj])[kk]->abpos));
                                        } else {
                                            edgelist_ms.push_back(
                                                    std::make_tuple(Node((*idx2[current][jj])[kk]->aid, 0),
                                                                    Node((*idx2[current][jj])[kk]->bid, 0), (*idx2[current][jj])[kk]->aepos -(*idx2[current][jj])[kk]->abpos));
                                        }
                                    }
                                    if ((reads[(*idx2[current][jj])[kk]->bid]->active) and ((*idx2[current][jj])[kk]->active))
                                    if (n1.id == 1)
                                    if (((*idx2[current][jj])[kk]->aln_type == BACKWARD) and (cb < 1)  and (candidate.find((*idx2[current][jj])[kk]->bid) != candidate.end() )) {
                                        cb += 1;
                                        //add edge
                                        printf("add edge\n");
                                        if ((*idx2[current][jj])[kk]->flags == 1) {
                                            edgelist_ms.push_back(
                                                    std::make_tuple(Node((*idx2[current][jj])[kk]->aid, 1),
                                                                    Node((*idx2[current][jj])[kk]->bid, 0), (*idx2[current][jj])[kk]->aepos -(*idx2[current][jj])[kk]->abpos));
                                        } else {
                                            edgelist_ms.push_back(
                                                    std::make_tuple(Node((*idx2[current][jj])[kk]->aid, 1),
                                                                    Node((*idx2[current][jj])[kk]->bid, 1),(*idx2[current][jj])[kk]->aepos -(*idx2[current][jj])[kk]->abpos));
                                        }
                                    }
                                    if ((cf >= 1) or (cb >= 1)) break;
                                }
                            }
                            if ((cf >= 1) or (cb >= 1)) break;
                        }

				}
	}

	printf("missing edge #, %d\n", edgelist_ms.size());


    for (int j = 0; j < edgelist_ms.size(); j++) {
        edgelist.push_back(edgelist_ms[j]);
    }


    std::ofstream out(argv[3], std::ofstream::out);
    //print edges list to file
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
