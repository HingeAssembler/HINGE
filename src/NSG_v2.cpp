#include <stdio.h>
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
        INSERT_ELEMENT(ACOVERB);
        INSERT_ELEMENT(BCOVEREA);
        INSERT_ELEMENT(INTERNAL);
        INSERT_ELEMENT(UNDIFINED);
        INSERT_ELEMENT(NOT_ACTIVE);
#undef INSERT_ELEMENT
    }
    return out << strings[value];
}



bool compare_overlap(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->aepos - ovl1->abpos + ovl1->bepos - ovl1->bbpos) > (ovl2->aepos - ovl2->abpos + ovl2->bepos - ovl2->bbpos));
}

bool compare_overlap_effective(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->eaepos - ovl1->eabpos + ovl1->ebepos - ovl1->ebbpos) > (ovl2->eaepos - ovl2->eabpos + ovl2->ebepos - ovl2->ebbpos));
}

bool compare_overlap_weight(LOverlap * ovl1, LOverlap * ovl2) {
    return (ovl1->weight > ovl2->weight);
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

    std::string name_mask = std::string(name_db) + ".mas";


	printf("name of db: %s, name of .las file %s\n", name_db, name_las);
    la.openDB(name_db);
    std::cout<<"# Reads:" << la.getReadNumber() << std::endl;
    la.openAlignmentFile(name_las);
    std::cout<<"# Alignments:" << la.getAlignmentNumber() << std::endl;


	int64 n_aln = la.getAlignmentNumber();
	int n_read = la.getReadNumber();
    std::vector<LOverlap *> aln;
	la.resetAlignment();
    la.getOverlap(aln,0,n_aln);

    std::vector<Read *> reads;
    la.getRead(reads,0,n_read);

	std::vector<std::vector<int>>  QV;
	la.getQV(QV,0,n_read);

    std::cout << "input data finished" <<std::endl;


    INIReader reader(name_config);

    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<<name_config<<std::endl;
        return 1;
    }

    int LENGTH_THRESHOLD = int(reader.GetInteger("filter", "length_threshold", -1));
    double QUALITY_THRESHOLD = reader.GetReal("filter", "quality_threshold", 0.0);
    int N_ITER = (int)reader.GetInteger("filter", "n_iter", -1);
    int ALN_THRESHOLD = (int)reader.GetInteger("filter", "aln_threshold", -1);
    int MIN_COV = (int)reader.GetInteger("filter", "min_cov", -1);
    int CUT_OFF = (int)reader.GetInteger("filter", "cut_off", -1);
    int THETA = (int)reader.GetInteger("filter", "theta", -1);
	int N_PROC = (int)reader.GetInteger("running", "n_proc", 4);

    omp_set_num_threads(N_PROC);

    //std::vector< std::vector<std::vector<LOverlap*>* > > idx2(n_read); //unordered_map from (aid) to alignments in a vector
    std::vector<Edge_w> edgelist, edgelist_ms; // save output to edgelist
    //std::unordered_map<int, std::vector <LOverlap * > >idx3,idx4; // this is the pileup
    std::vector<std::unordered_map<int, std::vector<LOverlap *> > > idx; //unordered_map from (aid, bid) to alignments in a vector
    std::vector<std::vector<LOverlap *>> idx2;
    std::vector<std::unordered_map<int, LOverlap *>> idx3;

/*
 * Index alignments by:
 * 1) read A - idx3
 * 2) read A-read B - idx
 */

    FILE * mask_file;
    mask_file = fopen(name_mask.c_str(), "r");
    int read, rs, re;
    while (fscanf(mask_file,"%d %d %d",&read, &rs, &re) != EOF) {
        reads[read]->effective_start = rs;
        reads[read]->effective_end = re;
    }
    std::cout<<"read mask finished" << std::endl;
    int num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) num_active_read ++;
    }
    std::cout<<"active reads:" << num_active_read<< std::endl;

    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->effective_end - reads[i]->effective_start < LENGTH_THRESHOLD) reads[i]->active = false;
        else num_active_read ++;
    }
    std::cout<<"active reads:" << num_active_read<< std::endl;

    for (int i = 0; i < n_read; i++) {
        idx.push_back(std::unordered_map<int, std::vector<LOverlap *> >());
        idx2.push_back(std::vector<LOverlap *>());
    }
//int num_finished = 0;

# pragma omp parallel for
    for (int i = 0; i < aln.size(); i++) {
        idx[aln[i]->aid][aln[i]->bid] = std::vector<LOverlap *>();
    }

    for (int i = 0; i < aln.size(); i++) {
        idx[aln[i]->aid][aln[i]->bid].push_back(aln[i]);
    }
    std::cout<<"index finished" <<std::endl;

    for (int i = 0; i < n_read; i++) {
        for (std::unordered_map<int, std::vector<LOverlap *> >::iterator it = idx[i].begin(); it!=idx[i].end(); it++) {
            std::sort(it->second.begin(), it->second.end(), compare_overlap);
            if (it->second.size() > 0)
                idx2[i].push_back(it->second[0]);
        }
    }

    int num_overlaps = 0;
    for (int i = 0; i < n_read; i++) {
        num_overlaps += idx2[i].size();
    }
    std::cout<<num_overlaps << " overlaps" << std::endl;

# pragma omp parallel for
    for (int i = 0; i < n_read; i++) {
        bool contained = false;
        for (int j = 0; j < idx2[i].size(); j++){
            idx2[i][j]->aes = reads[idx2[i][j]->aid]->effective_start;
            idx2[i][j]->aee = reads[idx2[i][j]->aid]->effective_end;

            if (idx2[i][j]->flags == 0) {
                idx2[i][j]->bes = reads[idx2[i][j]->bid]->effective_start;
                idx2[i][j]->bee = reads[idx2[i][j]->bid]->effective_end;
            } else {
                idx2[i][j]->bes = idx2[i][j]->blen - reads[idx2[i][j]->bid]->effective_end;
                idx2[i][j]->bee = idx2[i][j]->blen - reads[idx2[i][j]->bid]->effective_start;
            }


            /*for (int j = 0; j < aln[i]->trace_pts_len; j++) {
                std::cout << aln[i]->trace_pts[j] << " ";
            }
            std::cout<<std::endl;*/
            //printf("before %d %d %d %d ov %d %d %d %d\n", aln[i]->aes, aln[i]->aee, aln[i]->bes, aln[i]->bee, aln[i]->abpos, aln[i]->aepos, aln[i]->bbpos, aln[i]->bepos);
            idx2[i][j]->trim_overlap();
            //printf("after %d %d %d %d ov %d %d %d %d\n", aln[i]->aes, aln[i]->aee, aln[i]->bes, aln[i]->bee, aln[i]->eabpos, aln[i]->eaepos, aln[i]->ebbpos, aln[i]->ebepos);
            //num_finished ++;
            //if (num_finished%100000 == 0) printf("%lf percent finished\n", num_finished/double(aln.size())*100);
            if (((idx2[i][j]->ebepos - idx2[i][j]->ebbpos) < ALN_THRESHOLD) or ((idx2[i][j]->eaepos - idx2[i][j]->eabpos) < ALN_THRESHOLD))
            {
                idx2[i][j]->active = false;
                idx2[i][j]->aln_type = NOT_ACTIVE;
            } else {
                idx2[i][j]->addtype2(THETA);
                if (idx2[i][j]->aln_type == BCOVEREA)
                    contained = true;
            }
        }
        if (contained) reads[i]->active = false;
    }

    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) num_active_read ++;
    }
    std::cout<<"remove contained reads, active reads:" << num_active_read<< std::endl;

    for (int iter = 0; iter < N_ITER; iter++) {
        int remove = 0;
# pragma omp parallel for reduction(+:remove)
        for (int i = 0; i < n_read; i++) {
            if (reads[i]->active) {
                int forward = 0;
                int backward = 0;
                for (int j = 0; j < idx2[i].size(); j++)
                    if (idx2[i][j]->active) {
                        //idx2[i][j]->show();
                        //std::cout<<idx2[i][j]->aln_type << std::endl;
                        if ((idx2[i][j]->aln_type == FORWARD) and (reads[idx2[i][j]->bid]->active)) forward++;
                        else if ((idx2[i][j]->aln_type == BACKWARD) and (reads[idx2[i][j]->bid]->active)) backward++;
                    }
                if ((forward == 0) or (backward == 0)) {
                    reads[i]->active = false;
                    remove ++;
                }
                //printf("read %d forward %d backward %d\n", i, forward, backward);
            }
        }
        printf("remove %d reads\n", remove);
    }

    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) num_active_read ++;
    }
    std::cout<<"active reads:" << num_active_read<< std::endl;

    num_overlaps = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
        for (int j = 0; j < idx2[i].size(); j++)
            if (reads[idx2[i][j]->bid]->active) num_overlaps++;
    }
    std::cout<<num_overlaps << " overlaps" << std::endl;

    for (int i = 0; i < n_read; i++) {
        idx3.push_back(std::unordered_map<int, LOverlap*>() );
        if (reads[i]->active)
        for (int j = 0; j < idx2[i].size(); j++) {
            if (reads[idx2[i][j]->bid]->active)
                idx3[i][idx2[i][j]->bid] = idx2[i][j];
        }
    }

    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
            for (std::unordered_map<int, LOverlap*>::iterator it = idx3[i].begin(); it!=idx3[i].end(); it++) {
                int aid = i;
                int bid = it->second->bid;
                idx3[aid][bid]->weight = idx3[aid][bid]->eaepos - idx3[aid][bid]->eabpos + idx3[bid][aid]->eaepos - idx3[bid][aid]->eabpos;
                /*if (idx3[aid][bid]->aln_type == FORWARD) {
                    if (idx3[aid][bid]->flags == 0) idx3[bid][aid]->aln_type = BACKWARD;
                    else idx3[bid][aid]->aln_type = FORWARD;
                }
                if (idx3[aid][bid]->aln_type == BACKWARD) {
                    if (idx3[aid][bid]->flags == 0) idx3[bid][aid]->aln_type = FORWARD;
                    else idx3[bid][aid]->aln_type = BACKWARD;
                }*/
            }
    }


# pragma omp parallel for
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
            std::sort(idx2[i].begin(), idx2[i].end(), compare_overlap_weight);
    }

    FILE * out3;
    out3 = fopen("edges.backup.txt","w");
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
            for (int j = 0; j < idx2[i].size(); j++) {
                if (reads[idx2[i][j]->bid]->active)
                    fprintf(out3,"%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->flags, idx2[i][j]->aln_type,  idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
            }
    }


    FILE * out;
    FILE * out2;
    out = fopen((std::string(argv[3]) + ".1").c_str(), "w");
    out2 = fopen((std::string(argv[3]) + ".2").c_str(), "w");

    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            int forward = 0;
            int backward = 0;
            for (int j = 0; j < idx2[i].size(); j++)
                if (idx2[i][j]->active) {
                    if ((idx2[i][j]->aln_type == FORWARD) and (reads[idx2[i][j]->bid]->active)) {
                        if (forward < 1) {
                            //if (idx2[i][j]->flags == 0) fprintf(out, "%d eb %d %d\n%d be %d %d\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->aepos - idx2[i][j]->abpos, idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->aepos - idx2[i][j]->abpos);
                            //else fprintf(out, "%d ee %d %d\n%d ee %d %d\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->aepos - idx2[i][j]->abpos, idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->aepos - idx2[i][j]->abpos);
                            //if (idx2[i][j]->flags == 0) fprintf(out, "%d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            //else fprintf(out, "%d %d' [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            if (idx2[i][j]->flags == 0) fprintf(out, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            else fprintf(out, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);

                            if (idx2[i][j]->flags == 0) fprintf(out2, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->weight, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            else fprintf(out2, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->weight, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);

                        }
                        forward++;
                    }
                    else if ((idx2[i][j]->aln_type == BACKWARD) and (reads[idx2[i][j]->bid]->active)) {
                        if (backward < 1) {
                            //if (idx2[i][j]->flags == 0) fprintf(out, "%d eb %d %d\n%d be %d %d\n", idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->aepos - idx2[i][j]->abpos, idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->aepos - idx2[i][j]->abpos);
                            //else fprintf(out, "%d bb %d %d\n%d bb %d %d\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->aepos - idx2[i][j]->abpos, idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->aepos - idx2[i][j]->abpos);
                            //if (idx2[i][j]->flags == 0) fprintf(out, "%d' %d' [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            //else fprintf(out, "%d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);

                            if (idx2[i][j]->flags == 0) fprintf(out, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            else fprintf(out, "%d' %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);

                            if (idx2[i][j]->flags == 0) fprintf(out2, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->weight, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            else fprintf(out2, "%d' %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->weight, idx2[i][j]->eabpos, idx2[i][j]->eaepos, idx2[i][j]->ebbpos, idx2[i][j]->ebepos, idx2[i][j]->aes, idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);

                        }
                        backward++;
                    }
                }
        }
    }

    std::cout<<"sort and output finished" <<std::endl;

    /*for (int i = 0; i < n_read; i++) {
        std::cout<<i <<" "<<idx2[i].size() << std::endl;
    }*/

    la.closeDB(); //close database
    return 0;
}