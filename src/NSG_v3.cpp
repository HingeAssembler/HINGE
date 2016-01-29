#include <stdio.h>
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
        INSERT_ELEMENT(BCOVERA);
        INSERT_ELEMENT(INTERNAL);
        INSERT_ELEMENT(UNDEFINED);
        INSERT_ELEMENT(NOT_ACTIVE);
#undef INSERT_ELEMENT
    }
    return out << strings[value];
}



bool compare_overlap(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->read_A_match_end_ - ovl1->read_A_match_start_ + ovl1->read_B_match_end_ - ovl1->read_B_match_start_) >
            (ovl2->read_A_match_end_ - ovl2->read_A_match_start_ + ovl2->read_B_match_end_ - ovl2->read_B_match_start_));
}

bool compare_overlap_effective(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->eff_read_A_match_end_ - ovl1->eff_read_A_match_start_ + ovl1->eff_read_B_match_end_ - ovl1->eff_read_B_match_start_) >
            (ovl2->eff_read_A_match_end_ - ovl2->eff_read_A_match_start_ + ovl2->eff_read_B_match_end_ - ovl2->eff_read_B_match_start_));
}

bool compare_overlap_weight(LOverlap * ovl1, LOverlap * ovl2) {
    return (ovl1->weight > ovl2->weight);
}

bool compare_sum_overlaps(const std::vector<LOverlap * > * ovl1, const std::vector<LOverlap *> * ovl2) {
    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < ovl1->size(); i++)
        sum1 += (*ovl1)[i]->read_A_match_end_ - (*ovl1)[i]->read_A_match_start_ + (*ovl1)[i]->read_B_match_end_ - (*ovl1)[i]->read_B_match_start_;
    for (int i = 0; i < ovl2->size(); i++)
        sum2 += (*ovl2)[i]->read_A_match_end_ - (*ovl2)[i]->read_A_match_start_ + (*ovl2)[i]->read_B_match_end_ - (*ovl2)[i]->read_B_match_start_;
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

    int left= intervals[0]->read_A_match_start_ + cutoff, right = intervals[0]->read_A_match_end_ - cutoff;
    //left, right means maximal possible interval now

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

    std::string name_mask = std::string(name_db) + ".mas";
    std::string name_max = std::string(name_db) + ".max";
    std::string name_homo = std::string(name_db) + ".homologous.txt";
    std::string name_rep = std::string(name_db) + ".repeat.txt";
    std::string name_hg = std::string(name_db) + ".hinges.txt";
    std::string name_cov = std::string(name_db) + ".coverage.txt";
    std::ofstream maximal_reads(name_max);

    std::ifstream homo(name_homo);
    std::vector<int> homo_reads;
    int read_id;
    while (homo >> read_id) homo_reads.push_back(read_id);

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

    //std::vector< std::vector<std::vector<LOverlap*>* > > idx2(n_read);
    // unordered_map from (aid) to alignments in a vector
    std::vector<Edge_w> edgelist, edgelist_ms; // save output to edgelist
    //std::unordered_map<int, std::vector <LOverlap * > >idx3,idx4;
    // this is the pileup
    std::vector<std::unordered_map<int, std::vector<LOverlap *> > > idx; 
    /*
    	idx is a vector of length n_read, each element idx3[read A id] is a map, 
    	from read B id to a vector of overlaps
    */
    std::vector<std::vector<LOverlap *>> idx2;
    /*
    	idx2 is a vector of length n_read, each element idx2[read A id] is a vector, 
    	for each read B, we put the best overlap into that vector
    */
    std::vector<std::unordered_map<int, LOverlap *>> idx3;
	/*
		idx3 is a vector of length n_read, each element idx3[read A id] is a map, 
		from read read B id to the best overlap of read A and read B
	*/

    FILE * mask_file;
    mask_file = fopen(name_mask.c_str(), "r");
    int read, rs, re;
    while (fscanf(mask_file,"%d %d %d",&read, &rs, &re) != EOF) {
        reads[read]->effective_start = rs;
        reads[read]->effective_end = re;
    }
    std::cout<<"read mask finished" << std::endl;

    FILE * repeat_file;
    repeat_file = fopen(name_rep.c_str(), "r");
    FILE * hinge_file;
    hinge_file = fopen(name_hg.c_str(), "r");
    char * line = NULL;
    size_t len = 0;
    std::unordered_map<int, std::vector<std::pair<int, int>>> marked_repeats;

    while (getline(&line, &len, repeat_file) != -1) {
        std::stringstream ss;
        ss.clear();
        ss << line;
        int num;
        ss >> num;
        //printf("%d\n",num);
        marked_repeats[num] = std::vector<std::pair<int, int>>();
        int r1 = 0, r2 = 0;
        while (!ss.eof()) {
            r1 = 0;
            r2 = 0;
            ss >> r1 >> r2;
            if ((r1!=0) and (r2!=0)) {
                //printf("[%d %d]\n", r1, r2);
                marked_repeats[num].push_back(std::pair<int, int>(r1, r2));
            }
        }
        ss.clear();
    }
    fclose(repeat_file);
    std::cout<<"read marked repeats" << std::endl;
    std::unordered_map<int, std::vector<std::pair<int, int>>> marked_hinges;
    while (getline(&line, &len, hinge_file) != -1) {
        std::stringstream ss;
        ss << line;
        int num;
        ss >> num;
        //printf("%d\n",num);
        marked_hinges[num] = std::vector<std::pair<int, int>>();
        int r1 = 0, r2 = 0;
        while (!ss.eof()) {
            r1 = 0; r2 = 0;
            ss >> r1 >> r2;
            if ((r1 != 0) and (r2 != 0)) {
                //printf("[%d %d]\n", r1, r2);
                marked_hinges[num].push_back(std::pair<int, int>(r1, r2));
            }
        }
        ss.clear();
    }
    fclose(hinge_file);

    std::cout<<"read marked hinges" << std::endl;

    if (line)
        free(line);

    int num_active_read = 0;

    //This seems to be an unnecessary stub
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

    //for (int i = 0; i < homo_reads.size(); i++) // filter reads with homologous recombinations
    //    reads[homo_reads[i]]->active = false;

    for (int i = 0; i < n_read; i++) {
        idx.push_back(std::unordered_map<int, std::vector<LOverlap *> >());
        idx2.push_back(std::vector<LOverlap *>());
    }

//int num_finished = 0;

# pragma omp parallel for
    for (int i = 0; i < aln.size(); i++) {
        idx[aln[i]->read_A_id][aln[i]->read_B_id] = std::vector<LOverlap *>();
    }

    for (int i = 0; i < aln.size(); i++) {
        idx[aln[i]->read_A_id][aln[i]->read_B_id].push_back(aln[i]);
    }
    std::cout<<"index finished" <<std::endl;

    for (int i = 0; i < n_read; i++) {
        for (std::unordered_map<int, std::vector<LOverlap *> >::iterator it = idx[i].begin(); it!=idx[i].end(); it++) {
            std::sort(it->second.begin(), it->second.end(), compare_overlap);//Sort overlaps by lengths
            if (it->second.size() > 0)//Is this not just max? Why sort?
                idx2[i].push_back(it->second[0]);
        }
    }

    int num_overlaps = 0;
    for (int i = 0; i < n_read; i++) {//Isn't this just 0 or 1?
        num_overlaps += idx2[i].size();
    }
    std::cout<<num_overlaps << " overlaps" << std::endl;
//Figure out contained read
# pragma omp parallel for
    for (int i = 0; i < n_read; i++) {
        bool contained = false;
        for (int j = 0; j < idx2[i].size(); j++){
            idx2[i][j]->eff_read_A_start_ = reads[idx2[i][j]->read_A_id]->effective_start;
            idx2[i][j]->eff_read_A_end_ = reads[idx2[i][j]->read_A_id]->effective_end;

            if (idx2[i][j]->flags == 0) {//Where are flags set?
                idx2[i][j]->eff_read_B_start_ = reads[idx2[i][j]->read_B_id]->effective_start;
                idx2[i][j]->eff_read_B_end_ = reads[idx2[i][j]->read_B_id]->effective_end;
            } else {//looks like an iverted match
                idx2[i][j]->eff_read_B_start_ = idx2[i][j]->blen - reads[idx2[i][j]->read_B_id]->effective_end;
                idx2[i][j]->eff_read_B_end_ = idx2[i][j]->blen - reads[idx2[i][j]->read_B_id]->effective_start;
            }


            /*for (int j = 0; j < aln[i]->trace_pts_len; j++) {
                std::cout << aln[i]->trace_pts[j] << " ";
            }
            std::cout<<std::endl;*/
            //printf("before %d %d %d %d ov %d %d %d %d\n", aln[i]->aes, aln[i]->aee,
            // aln[i]->bes, aln[i]->bee, aln[i]->read_A_match_start_, aln[i]->read_A_match_end_, aln[i]->read_B_match_start_, aln[i]->read_B_match_end_);
            idx2[i][j]->trim_overlap();//What does this do?
            //printf("after %d %d %d %d ov %d %d %d %d\n", aln[i]->aes, aln[i]->aee,
            // aln[i]->bes, aln[i]->bee, aln[i]->eff_read_A_match_start_, aln[i]->eff_read_A_match_end_, aln[i]->eff_read_B_match_start_, aln[i]->eff_read_B_match_end_);
            //num_finished ++;
            //if (num_finished%100000 == 0) printf("%lf percent finished\n", num_finished/double(aln.size())*100);
            if (((idx2[i][j]->eff_read_B_match_end_ - idx2[i][j]->eff_read_B_match_start_) < ALN_THRESHOLD)
                or ((idx2[i][j]->eff_read_A_match_end_ - idx2[i][j]->eff_read_A_match_start_) < ALN_THRESHOLD))
            {
                idx2[i][j]->active = false;
                idx2[i][j]->aln_type = NOT_ACTIVE;
            } else {
                idx2[i][j]->addtype2(THETA);
                if (idx2[i][j]->aln_type == BCOVERA)
                    contained = true;
            }
        }
        if (contained) reads[i]->active = false;
    }

    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            num_active_read ++;
            maximal_reads << i << std::endl;
        }
    }
    std::cout<<"removed contained reads, active reads:" << num_active_read<< std::endl;

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
                        if ((idx2[i][j]->aln_type == FORWARD) and (reads[idx2[i][j]->read_B_id]->active)) forward++;
                        else if ((idx2[i][j]->aln_type == BACKWARD) and (reads[idx2[i][j]->read_B_id]->active)) backward++;
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
            if (reads[idx2[i][j]->read_B_id]->active) num_overlaps++;
    }
    std::cout<<num_overlaps << " overlaps" << std::endl;

    for (int i = 0; i < n_read; i++) {
        idx3.push_back(std::unordered_map<int, LOverlap*>() );
        if (reads[i]->active)
        for (int j = 0; j < idx2[i].size(); j++) {
            if (reads[idx2[i][j]->read_B_id]->active)
                idx3[i][idx2[i][j]->read_B_id] = idx2[i][j];
        }
    }

    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
            for (std::unordered_map<int, LOverlap*>::iterator it = idx3[i].begin(); it!=idx3[i].end(); it++) {
                int aid = i;
                int bid = it->second->bid;
                idx3[aid][bid]->weight =
                        idx3[aid][bid]->eff_aepos - idx3[aid][bid]->eff_abpos
                        + idx3[bid][aid]->eff_aepos - idx3[bid][aid]->eff_abpos;
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
                if (reads[idx2[i][j]->read_B_id]->active)
                    fprintf(out3, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n",
                            idx2[i][j]->read_A_id, idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->flags,
                            idx2[i][j]->aln_type, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_,
                            idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_, idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
            }
    }



    FILE * out;
    FILE * out2;
    out = fopen((std::string(argv[3]) + ".1").c_str(), "w");
    out2 = fopen((std::string(argv[3]) + ".2").c_str(), "w");
    out3 = fopen((std::string(argv[3]) + ".hinges").c_str(), "w");




    class Hinge {
    public:
        int pos;
        int type; // 1, -1
        bool active;
        bool active2;
        Hinge(int pos, int t, bool active):pos(pos),type(t), active(active), active2(false) {};
        Hinge():pos(0),type(1), active(true) {};
    };

    std::unordered_map<int, std::vector<Hinge> > hinges_vec;



    int n = 0;
    for (int i = 0; i < n_read; i++) {
        hinges_vec[i] = std::vector<Hinge>();
        for (int j = 0; j < marked_hinges[i].size(); j++) {
            hinges_vec[i].push_back(Hinge(marked_hinges[i][j].first, marked_hinges[i][j].second , true));
            if (reads[i]->active) {
                n ++;
                //printf("%d %d %d\n", i, marked_hinges[i][j].first, marked_hinges[i][j].second);
            }
        }
    }

    printf("%d hinges\n", n);


    n = 0;
    for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < hinges_vec[i].size(); j++) {
            if ((reads[i]->active) and (hinges_vec[i][j].active)) n++;
        }
    }
    printf("%d active hinges\n", n);





    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            int forward = 0;
            int backward = 0;
            for (int j = 0; j < idx2[i].size(); j++)
                if (idx2[i][j]->active) {
                    if ((idx2[i][j]->aln_type == FORWARD) and (reads[idx2[i][j]->read_B_id]->active)) {
                        if (forward < 1) {

                            /*if (idx2[i][j]->flags == 0)
                             fprintf(out, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                             idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_,
                             idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->aes,
                             idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            else
                             fprintf(out, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                             idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_,
                             idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->aes,
                             idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);

                            if (idx2[i][j]->flags == 0)
                             fprintf(out2, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                             idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_,
                             idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->aes,
                             idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            else
                             fprintf(out2, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                             idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_,
                             idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->aes,
                             idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            */
                            //remove certain hinges

                            for (int k = 0; k < hinges_vec[i].size(); k++) {
                                if ((idx2[i][j]->eff_read_A_match_start_ < hinges_vec[i][k].pos - 400) and (hinges_vec[i][k].type == 1))
                                    hinges_vec[i][k].active = false;
                            }

                            /*for (int k = 0; k < hinges_vec[idx2[i][j]->bid].size(); k++) {
                                if ((hinges_vec[idx2[i][j]->bid][k].type == 1)
                                and (idx2[i][j]->flags == 1) and
                                (idx2[i][j]->eff_read_B_match_end_ > idx2[i][j]->blen - hinges_vec[idx2[i][j]->bid][k].pos + 300))
                                    hinges_vec[idx2[i][j]->bid][k].active2 = false;

                                if ((hinges_vec[idx2[i][j]->bid][k].type == -1)
                                and (idx2[i][j]->flags == 0)
                                and (idx2[i][j]->eff_read_B_match_end_ > hinges_vec[idx2[i][j]->bid][k].pos + 300))
                                    hinges_vec[idx2[i][j]->bid][k].active2 = false;
                            }*/

                            //if ((repeat_status_back[i])
                            // and (idx2[i][j]->eff_read_A_match_start_ < marked_repeats[i].front().first - 200)) {
                            //    repeat_status_back[i] = false;
                                //printf("remove %d\n",i);
                            //}


                        }
                        forward++;
                    }
                    else if ((idx2[i][j]->aln_type == BACKWARD) and (reads[idx2[i][j]->read_B_id]->active)) {
                        if (backward < 1) {

                            /*if (idx2[i][j]->flags == 0)
                             fprintf(out, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                             idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_,
                             idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->aes,
                             idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            else
                             fprintf(out, "%d' %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                             idx2[i][j]->aid, idx2[i][j]->bid, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_,
                             idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->aes,
                             idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);

                            if (idx2[i][j]->flags == 0)
                             fprintf(out2, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                             idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_,
                             idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->aes,
                             idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            else
                             fprintf(out2, "%d' %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                             idx2[i][j]->bid, idx2[i][j]->aid, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_,
                             idx2[i][j]->eff_read_A_match_end_, idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->aes,
                             idx2[i][j]->aee, idx2[i][j]->bes, idx2[i][j]->bee);
                            */
                            // remove certain hinges

                            //if ((repeat_status_front[i])
                            // and (idx2[i][j]->eff_read_A_match_end_ > marked_repeats[i].back().second + 200)
                            // and (idx2[i][j]->eff_read_A_match_end_ > marked_repeats[i].front().second + 200)) {
                            //    repeat_status_front[i] = false;
                                //printf("remove %d\n",i);
                            //}

                            for (int k = 0; k < hinges_vec[i].size(); k++) {
                                if ((idx2[i][j]->eff_read_A_match_end_ > hinges_vec[i][k].pos + 400) and (hinges_vec[i][k].type == -1))
                                    hinges_vec[i][k].active = false;
                            }


                            /*for (int k = 0; k < hinges_vec[idx2[i][j]->bid].size(); k++) {
                                if ((hinges_vec[idx2[i][j]->bid][k].type == 1)
                                and (idx2[i][j]->flags == 0)
                                and (idx2[i][j]->eff_read_B_match_start_ < hinges_vec[idx2[i][j]->bid][k].pos - 300))
                                    hinges_vec[idx2[i][j]->bid][k].active2 = false;

                                if ((hinges_vec[idx2[i][j]->bid][k].type == -1)
                                and (idx2[i][j]->flags == 1)
                                and (idx2[i][j]->eff_read_B_match_start_ < idx2[i][j]->blen - hinges_vec[idx2[i][j]->bid][k].pos - 300))
                                    hinges_vec[idx2[i][j]->bid][k].active2 = false;
                            }*/

                        }
                        backward++;
                    }
                }
        }
    }


    n = 0;
    for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < hinges_vec[i].size(); j++) {
            if ((reads[i]->active) and ((hinges_vec[i][j].active) or hinges_vec[i][j].active2)) {
                printf("%d %d %d\n", i, marked_hinges[i][j].first, marked_hinges[i][j].second);
                n++;
            }
        }
    }
    printf("after filter %d active hinges\n", n);


    // filter hinges
    //
    std::vector<bool> repeat_status_front;
    std::vector<bool> repeat_status_back;


    for (int i = 0; i < n_read; i++) {
        bool in = false;
        bool out = false;
        for (int j = 0; j < hinges_vec[i].size(); j++) {
            if (((hinges_vec[i][j].active) or (hinges_vec[i][j].active2)) and (hinges_vec[i][j].type == 1)) in = true;
            if (((hinges_vec[i][j].active) or (hinges_vec[i][j].active2)) and (hinges_vec[i][j].type == -1)) out = true;
        }
        repeat_status_front.push_back(out);
        repeat_status_back.push_back(in);
    }




    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            int forward = 0;
            int backward = 0;
            for (int j = 0; j < idx2[i].size(); j++)
                if (idx2[i][j]->active) {
                    if ((idx2[i][j]->aln_type == FORWARD) and (reads[idx2[i][j]->read_B_id]->active)) {
                        /*if (not repeat_status_back[i])*/ {
                            if (forward < 1) {

                                if (idx2[i][j]->flags == 0)
                                    fprintf(out, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_A_id,
                                            idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                            idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                            idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                                else
                                    fprintf(out, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_A_id,
                                            idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                            idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                            idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);

                                if (idx2[i][j]->flags == 0)
                                    fprintf(out2, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_B_id,
                                            idx2[i][j]->read_A_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                            idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                            idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                                else
                                    fprintf(out2, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_B_id,
                                            idx2[i][j]->read_A_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                            idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                            idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);

                            }
                        }
                        forward++;
                    }
                    else if ((idx2[i][j]->aln_type == BACKWARD) and (reads[idx2[i][j]->read_B_id]->active)) {
                        if (backward < 1) {

                            /*if (not repeat_status_front[i])*/ {
                                if (idx2[i][j]->flags == 0)
                                    fprintf(out, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_A_id,
                                            idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                            idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                            idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                                else
                                    fprintf(out, "%d' %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_A_id,
                                            idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                            idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                            idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);

                                if (idx2[i][j]->flags == 0)
                                    fprintf(out2, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_B_id,
                                            idx2[i][j]->read_A_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                            idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                            idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                                else
                                    fprintf(out2, "%d' %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_B_id,
                                            idx2[i][j]->read_A_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                            idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                            idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                            }
                        }
                        backward++;
                    }
                }
        }
    }




    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {

            for (int j = 0; j < idx2[i].size(); j++)
                if (idx2[i][j]->active) {
                    if ((idx2[i][j]->aln_type == FORWARD) and (reads[idx2[i][j]->read_B_id]->active)) {
                        if (repeat_status_back[i]) // choose the correct repeat status (with hinges)
                        if (((idx2[i][j]->flags == 0) and repeat_status_front[idx2[i][j]->read_B_id])
                            or ((idx2[i][j]->flags == 1) and repeat_status_back[idx2[i][j]->read_B_id])) {
                            if (idx2[i][j]->flags == 0)
                                fprintf(out3, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_A_id,
                                        idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                        idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                        idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                            else
                                fprintf(out3, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_A_id,
                                        idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                        idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                        idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                            if (idx2[i][j]->flags == 0)
                                fprintf(out3, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_B_id,
                                        idx2[i][j]->read_A_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                        idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                        idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                            else
                                fprintf(out3, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_B_id,
                                        idx2[i][j]->read_A_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                        idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                        idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                        }



                    }
                    else if  ((idx2[i][j]->aln_type == BACKWARD) and (reads[idx2[i][j]->read_B_id]->active)) {
                        if (repeat_status_front[i])
                        if (((idx2[i][j]->flags == 0) and repeat_status_back[idx2[i][j]->read_B_id])
                            or ((idx2[i][j]->flags == 1) and repeat_status_front[idx2[i][j]->read_B_id])) {
                            if (idx2[i][j]->flags == 0)
                                fprintf(out3, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_A_id,
                                        idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                        idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                        idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                            else
                                fprintf(out3, "%d' %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_A_id,
                                        idx2[i][j]->read_B_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                        idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                        idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);

                            if (idx2[i][j]->flags == 0)
                                fprintf(out3, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_B_id,
                                        idx2[i][j]->read_A_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                        idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                        idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                            else
                                fprintf(out3, "%d' %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", idx2[i][j]->read_B_id,
                                        idx2[i][j]->read_A_id, idx2[i][j]->weight, idx2[i][j]->eff_read_A_match_start_, idx2[i][j]->eff_read_A_match_end_,
                                        idx2[i][j]->eff_read_B_match_start_, idx2[i][j]->eff_read_B_match_end_, idx2[i][j]->eff_read_A_start_, idx2[i][j]->eff_read_A_end_,
                                        idx2[i][j]->eff_read_B_start_, idx2[i][j]->eff_read_B_end_);
                        }

                    }

                }
        }
    }




    //printf("%d %d %d\n", n_read, repeat_status_front.size(), repeat_status_back.size());

    //
    //second pass, for those with repeats
    //

    std::cout<<"sort and output finished" <<std::endl;

    /*for (int i = 0; i < n_read; i++) {
        std::cout<<i <<" "<<idx2[i].size() << std::endl;
    }*/

    la.closeDB(); //close database
    return 0;
}
