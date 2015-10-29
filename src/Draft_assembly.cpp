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
#include <sstream>
#include <iostream>
#include <set>
#include <omp.h>
#include "INIReader.h"



#define LAST_READ_SYMBOL  '$'

typedef std::tuple<Node, Node, int> Edge_w;

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

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


std::string reverse_complement(std::string seq) {
    static std::map<char, char> m = {{'a','t'}, {'c','g'}, {'g','c'}, {'t','a'}, {'A','T'}, {'C','G'}, {'T','A'}, {'G','C'}, {'n','n'}, {'N', 'N'}};
    std::reverse(seq.begin(), seq.end());
    for (int i = 0; i < seq.size(); i++) {
        seq[i] = m[seq[i]];
    }
    return seq;
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
    char * name_input = argv[3];
	char * name_output = argv[4];
	char * name_config = argv[5];
	
	printf("name of db: %s, name of .las file %s\n", name_db, name_las);
    la.OpenDB(name_db);
	//la.OpenDB2(name_db); // changed this on Oct 12, may case problem, xf1280@gmail.com
    std::cout<<"# Reads:" << la.getReadNumber() << std::endl;
    la.OpenAlignment(name_las);
    std::cout<<"# Alignments:" << la.getAlignmentNumber() << std::endl;
	//la.resetAlignment();
	//la.showOverlap(0,1);

	int n_aln = la.getAlignmentNumber();
	int n_read = la.getReadNumber();
    std::vector<LOverlap *> aln;
	la.resetAlignment();
    la.getOverlap(aln,0,n_aln);

    std::vector<Read *> reads;
    la.getRead(reads,0,n_read);

    std::cout << "input data finished" <<std::endl;
	/**
	filter reads
	**/

	/**
	 * Remove reads shorter than length threshold
	 */

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
    //int CHI_THRESHOLD = 500; // threshold for chimeric/adaptor at the begining
    int N_ITER = reader.GetInteger("filter", "n_iter", -1);
    int ALN_THRESHOLD = reader.GetInteger("filter", "aln_threshold", -1);
    int MIN_COV = reader.GetInteger("filter", "min_cov", -1);
    int CUT_OFF = reader.GetInteger("filter", "cut_off", -1);
    int THETA = reader.GetInteger("filter", "theta", -1);

	int N_PROC = reader.GetInteger("running", "n_proc", 4);

    omp_set_num_threads(N_PROC);

    //std::unordered_map<std::pair<int,int>, std::vector<LOverlap *> > idx; //unordered_map from (aid, bid) to alignments in a vector
    std::vector< std::vector<std::vector<LOverlap*>* > > idx2(n_read); //unordered_map from (aid) to alignments in a vector
    std::vector< Edge_w > edgelist; // save output to edgelist
    std::unordered_map<int, std::vector <LOverlap * > >idx3; // this is the pileup
    std::vector<std::set<int> > has_overlap(n_read);
    std::unordered_map<int, std::unordered_map<int, std::vector<LOverlap *> > > idx;



    for (int i = 0; i< n_read; i++) {
        //has_overlap[i] = std::set<int>();
        idx3[i] = std::vector<LOverlap *>();
    }

    //for (int i = 0; i < aln.size(); i++)
    //    if (aln[i]->active)
    //        idx[std::pair<int, int>(aln[i]->aid, aln[i]->bid)] = std::vector<LOverlap *>();
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
        }
    }


    std::cout<<"add data"<<std::endl;
    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx[aln[i]->aid][aln[i]->bid].push_back(aln[i]);
        }
    }
    std::cout<<"add data"<<std::endl;


    //sort each a,b pair according to abpos:
    /*for (int i = 0; i < n_read; i++)
        for (std::set<int>::iterator j = has_overlap[i].begin(); j != has_overlap[i].end(); j++) {
            std::sort(idx[std::pair<int,int>(i, *j)].begin(), idx[std::pair<int,int>(i, *j)].end(), compare_pos);
        }
    */
    for (int i = 0; i < n_read; i++)
        for (std::set<int>::iterator j = has_overlap[i].begin(); j != has_overlap[i].end(); j++) {
            idx2[i].push_back(&(idx[i][*j]));
    }

    std::cout<<"add data"<<std::endl;
    std::unordered_map<int,std::vector<Interval> > covered_region;


    //for (int i = 0; i < n_read; i++) {
    //    printf("read %d [%d %d]/%d\n", i, reads[i]->effective_start, reads[i]->effective_end, reads[i]->len);
    //}

    /**Filtering
     *
     **/

    for (int n_iter = 0; n_iter < N_ITER; n_iter ++) {

#pragma omp parallel for
        for (int i = 0; i < n_read; i++) {
            if (reads[i]->active) {
                //covered_region[i] = Merge(idx3[i]);
                reads[i]->intervals = Merge(idx3[i], CUT_OFF);
                /*if (reads[i]->intervals.empty()) {
                    reads[i]->effective_start = 0;
                    reads[i]->effective_end = 0;
                } else {
                    reads[i]->effective_start = reads[i]->intervals.front().first;
                    reads[i]->effective_end = reads[i]->intervals.back().second;
                }*/
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
                 LENGTH_THRESHOLD) /*or (reads[i]->len - reads[i]->effective_end > CHI_THRESHOLD) or (reads[i]->effective_start > CHI_THRESHOLD)*/
                or (reads[i]->intervals.size() != 1))
                reads[i]->active = false;
        } // filter according to effective length, and interval size

        std::cout<<"filter data"<<std::endl;

        int num_active = 0;

        for (int i = 0; i < n_read; i++) {
            if (reads[i]->active)
                num_active++;
        }
        std::cout << "num active reads " << num_active << std::endl;
# pragma omp parallel for
        for (int i = 0; i < n_aln; i++) {
            if (aln[i]->active)
            if ((not reads[aln[i]->aid]->active) or (not reads[aln[i]->bid]->active) or
                (aln[i]->diffs / float(aln[i]->bepos - aln[i]->bbpos + aln[i]->aepos - aln[i]->abpos) >
                 QUALITY_THRESHOLD) or (aln[i]->aepos - aln[i]->abpos < ALN_THRESHOLD))
                aln[i]->active = false;
        }

# pragma omp parallel for
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


# pragma omp parallel for
    for (int i = 0; i < n_read; i++ ) {
        std::sort( idx2[i].begin(), idx2[i].end(), compare_sum_overlaps );
    }


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
                    //idx2[i][j]->show();
                    if ((*idx2[i][j])[0]->active) {
                        for (int kk = 0; kk < idx2[i][j]->size(); kk++) {
                            if (reads[(*idx2[i][j])[kk]->bid]->active)
                            if (((*idx2[i][j])[kk]->aln_type == FORWARD) and (cf < 1)) {
                                cf += 1;
                                //add edge
                                /*if ((*idx2[i][j])[kk]->flags == 1) { // n = 0, c = 1
                                    edgelist.push_back(
                                            std::pair<Node, Node>(Node((*idx2[i][j])[kk]->aid, 0),
                                                                  Node((*idx2[i][j])[kk]->bid, 1)));
                                } else {
                                    edgelist.push_back(
                                            std::pair<Node, Node>(Node((*idx2[i][j])[kk]->aid, 0),
                                                                  Node((*idx2[i][j])[kk]->bid, 0)));
                                }*/
                            }
                            if (reads[(*idx2[i][j])[kk]->bid]->active)
                            if (((*idx2[i][j])[kk]->aln_type == BACKWARD) and (cb < 1)) {
                                cb += 1;
                                //add edge
                                /*if ((*idx2[i][j])[kk]->flags == 1) {
                                    edgelist.push_back(
                                            std::pair<Node, Node>(Node((*idx2[i][j])[kk]->aid, 1),
                                                                  Node((*idx2[i][j])[kk]->bid, 0)));
                                } else {
                                    edgelist.push_back(
                                            std::pair<Node, Node>(Node((*idx2[i][j])[kk]->aid, 1),
                                                                  Node((*idx2[i][j])[kk]->bid, 1)));
                                }*/
                            }
                            if ((cf >= 1) and (cb >= 1)) break;
                        }
                    }
                    if ((cf >= 1) and (cb >= 1)) break;
                }

                if (reads[i]->active) {
                    if (cf == 0) printf("%d has no out-going edges\n", i);
                    if (cb == 0) printf("%d has no in-coming edges\n", i);
                    if ((cf == 0) or (cb ==0)) {
                        no_edge++;
                        reads[i]->active = false;
                    }
                }

            }

            printf("no_edge nodes:%d\n",no_edge);
		}
		
		
		
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


	edgelist.clear();
    std::string edge_line;
    std::ifstream edges_file(name_input);
    while (!edges_file.eof()) {
        std::getline(edges_file,edge_line);
        std::vector<std::string> tokens = split(edge_line, ' ');
        Node node0;
        Node node1;
        int w;
        if (tokens.size() == 3) {
            if (tokens[0].back() == '\'') {
                node0.id = std::stoi(tokens[0].substr(0, tokens[0].size() - 1));
                node0.strand = 1;
            } else {
                node0.id = std::stoi(tokens[0]);
                node0.strand = 0;
            }

            if (tokens[1].back() == '\'') {
                node1.id = std::stoi(tokens[1].substr(0, tokens[1].size() - 1));
                node1.strand = 1;
            } else {
                node1.id = std::stoi(tokens[1]);
                node1.strand = 0;


            }
            w = std::stoi(tokens[2]);
            edgelist.push_back(std::make_tuple(node0,node1,w));
        }
    }

    //edgelist.push_back(std::pair<Node, Node>(edgelist.back().second, edgelist.front().first));

	for (int i = 0; i < edgelist.size(); i++) {
		printf("%d->%d\n", std::get<0>(edgelist[i]), std::get<1>(edgelist[i]));
	}
	
	std::string sequence = "";
	for (int i = 0; i < edgelist.size(); i++){
		if (i == 0) {
            if (std::get<0>(edgelist[i]).strand == 0)
                sequence.append(reads[std::get<0>(edgelist[i]).id]->bases);
            else
                sequence.append(reverse_complement(reads[std::get<0>(edgelist[i]).id]->bases));
		}
		
		std::vector<LOverlap *> currentalns = idx[std::get<0>(edgelist[i]).id][std::get<1>(edgelist[i]).id];
    	LOverlap * currentaln = NULL;
		for (int j = 0; j < currentalns.size(); j++) {
    		std::cout << std::get<0>(edgelist[i]).id << " " << std::get<1>(edgelist[i]).id << " " << currentalns[j]->aln_type << std::endl;
    		if (currentalns[j]->aepos - currentalns[j]->abpos == std::get<2>(edgelist[i]) ) currentaln = currentalns[j];
		}

        std::string current_seq;

        if (std::get<1>(edgelist[i]).strand == 0)
            current_seq = reads[std::get<1>(edgelist[i]).id]->bases;
        else
            current_seq = reverse_complement(reads[std::get<1>(edgelist[i]).id]->bases);

        if (std::get<0>(edgelist[i]).strand == 0) {
            sequence.erase(sequence.end() - (currentaln->alen - currentaln->aepos), sequence.end());
            current_seq.erase(current_seq.begin(), current_seq.begin() + currentaln->bepos);
            sequence.append(current_seq);
        }
        else {
            sequence.erase(sequence.end() - (currentaln->abpos), sequence.end());
            current_seq.erase(current_seq.begin(), current_seq.begin() + currentaln->blen - currentaln->bbpos);
            sequence.append(current_seq);
        }

	}

    
    //need to trim the end


	
	std::cout<<sequence.size()<<std::endl;
	std::ofstream out(name_output);
	
	out << ">draftassembly" << std::endl;
	out << sequence << std::endl;

    la.CloseDB(); //close database
    return 0;
}

