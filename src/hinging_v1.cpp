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
#include <iomanip>


#define LAST_READ_SYMBOL  '$'

#define HINGED_EDGE 1
#define UNHINGED_EDGE -1
#define REVERSE_COMPLEMENT_MATCH 1
#define SAME_DIRECTION_MATCH 0

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
    return ((ovl1->read_A_match_end_ - ovl1->read_A_match_start_
             + ovl1->read_B_match_end_ - ovl1->read_B_match_start_) >
            (ovl2->read_A_match_end_ - ovl2->read_A_match_start_
             + ovl2->read_B_match_end_ - ovl2->read_B_match_start_));
}

bool compare_overlap_effective(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->eff_read_A_match_end_ - ovl1->eff_read_A_match_start_
             + ovl1->eff_read_B_match_end_ - ovl1->eff_read_B_match_start_) >
            (ovl2->eff_read_A_match_end_ - ovl2->eff_read_A_match_start_
             + ovl2->eff_read_B_match_end_ - ovl2->eff_read_B_match_start_));
}

bool compare_overlap_weight(LOverlap * ovl1, LOverlap * ovl2) {
    return (ovl1->weight > ovl2->weight);
}

bool compare_sum_overlaps(const std::vector<LOverlap * > * ovl1, const std::vector<LOverlap *> * ovl2) {
    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < ovl1->size(); i++)
        sum1 += (*ovl1)[i]->read_A_match_end_ - (*ovl1)[i]->read_A_match_start_
                + (*ovl1)[i]->read_B_match_end_ - (*ovl1)[i]->read_B_match_start_;
    for (int i = 0; i < ovl2->size(); i++)
        sum2 += (*ovl2)[i]->read_A_match_end_ - (*ovl2)[i]->read_A_match_start_
                + (*ovl2)[i]->read_B_match_end_ - (*ovl2)[i]->read_B_match_start_;
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

bool ProcessAlignment(LOverlap * match, Read * read_A, Read * read_B, int ALN_THRESHOLD, int THETA){
    //Function takes as input pointers to a match, and the read_A and read_B of that match, set constants
    //ALN_THRESHOLD and THETA
    //It inputs the effective read start and end into the match class object
    //Next it trims match
    //Finally it figures out the type of match we have here by calling AddTypesAsymmetric() on the
    //class object

    //std::cout<<" In ProcessAlignment"<<std::endl;
    bool contained=false;
    match->eff_read_A_start_ = read_A->effective_start;
    match->eff_read_A_end_ = read_A->effective_end;


    match->eff_read_B_start_ = read_B->effective_start;
    match->eff_read_B_end_ = read_B->effective_end;


    match->trim_overlap();
    //std::cout<< contained<<std::endl;
    if (((match->eff_read_B_match_end_ - match->eff_read_B_match_start_) < ALN_THRESHOLD)
        or ((match->eff_read_A_match_end_ - match->eff_read_A_match_start_) < ALN_THRESHOLD))
    {   std::ofstream ofs ("Hinging.if.txt", std::ofstream::app);
        ofs << "In deactivating match "
        << "Read A id "<< match->read_A_id_ << " "<< read_A->id
        << " Read A eff start " << read_A->effective_start
        << " Read A eff end " << read_A->effective_end
        << " Read A length " << read_A->len
        << " Read A match start "<<  match->eff_read_A_match_start_
        << " Read A match end " << match->eff_read_A_match_end_
        << "Read B id "<< match->read_B_id_ << " "<< read_B->id
        << " Read B start " << read_B->effective_start
        << " Read B end " << read_B->effective_end
        << " Read B length " << read_B->len
        << " Read B match start "<<  match->eff_read_B_match_start_
        << " Read B match end " << match->eff_read_B_match_end_
        << " Match direction "<< match->reverse_complement_match_
        << "\n" << std::endl;
        ofs.close();
        match->active = false;
        match->match_type_ = NOT_ACTIVE;
    } else {
        match->AddTypesAsymmetric(THETA);
        if (match->match_type_ == BCOVERA) {
            std::ofstream ofs ("Hinging.else.txt", std::ofstream::app);
            ofs <<  "===============================================\n"
                << "Read A id "<< std::setfill('0') << std::setw(5) <<match->read_A_id_
                << " " << std::setfill('0') << std::setw(5) << read_A->id
                << "\nRead B id "  << std::setfill('0') << std::setw(5) << match->read_B_id_
                << " "<< std::setfill('0') << std::setw(5) <<  read_B->id
                << "\nRead A eff start "<< std::setfill('0') << std::setw(5)  << read_A->effective_start
                << " Read A eff end "<< std::setfill('0') << std::setw(5)  << read_A->effective_end
                << " Read A length " << std::setfill('0') << std::setw(5)  << read_A->len
                << " Read A match start "<< std::setfill('0') << std::setw(5) <<  match->read_A_match_start_
                << " Read A eff match start " << std::setfill('0') << std::setw(5) <<  match->eff_read_A_match_start_
                << " Read A match end " << std::setfill('0') << std::setw(5)  << match->read_A_match_end_
                << " Read A eff match end " << std::setfill('0') << std::setw(5)  << match->eff_read_A_match_end_
                << "\nRead B eff start "  << std::setfill('0') << std::setw(5) << read_B->effective_start
                << " Read B eff end " << std::setfill('0') << std::setw(5)  << read_B->effective_end
                << " Read B length " << std::setfill('0') << std::setw(5)  << read_B->len
                << " Read B match start "<< std::setfill('0') << std::setw(5) <<  match->read_B_match_start_
                << " Read B eff match start " << std::setfill('0') << std::setw(5) <<  match->eff_read_B_match_start_
                << " Read B match end " << std::setfill('0') << std::setw(5)  << match->read_B_match_end_
                << " Read B eff match end "  << std::setfill('0') << std::setw(5)  << match->eff_read_B_match_end_
                << "\nMatch direction "  << std::setfill('0') << std::setw(5)  << match->reverse_complement_match_
                << "\n" << std::endl;
            contained = true;
            ofs.close();
        }
        //std::cout<< contained<< std::endl;
    }

    match->weight =
            match->eff_read_A_match_end_ - match->eff_read_A_match_start_
            + match->eff_read_B_match_end_ - match->eff_read_B_match_start_;
    
    return contained;
}

class Hinge {
public:
    int pos;
    int type; // 1, -1
    bool active;
    bool active2;
    Hinge(int pos, int t, bool active):pos(pos),type(t), active(active), active2(false) {};
    Hinge():pos(0),type(1), active(true) {};
};

bool isValidHinge(LOverlap *match, std::vector<Hinge> &read_hinges){
    //Returns true if read_hinges (a vector of all hinges corresponding to a read )
    // has a hinge of appropriate type within tolerance from positions of start of the
    // overlap on read_B of the overlap given.
    int tolerance=100;//TODO put as #define
    int position=match->eff_read_B_match_start_;
    int type; //TODO : Make enum
    if (match->match_type_==FORWARD_INTERNAL)
        type=1;
    else if (match->match_type_==BACKWARD_INTERNAL)
        type=-1;

    if (match->reverse_complement_match_==1){
        type=-type;
        position=match->eff_read_B_match_end_;
    }

    bool valid=false;
    for (int index=0; index < read_hinges.size(); index++) {
        if ((abs(position - read_hinges[index].pos) < tolerance) and (type == read_hinges[index].type))
            valid = true;
        return valid;
    }
}

void PrintOverlapToFile(FILE * file_pointer, LOverlap * match) {

    int direction = match->reverse_complement_match_;
    int hinged;

    if ((match->match_type_ == FORWARD) or (match->match_type_ == BACKWARD))
        hinged = UNHINGED_EDGE;

    else if ((match->match_type_ == FORWARD_INTERNAL) or (match->match_type_ == BACKWARD_INTERNAL))
        hinged = HINGED_EDGE;

    if ((match->match_type_ == FORWARD_INTERNAL) or (match->match_type_ == FORWARD)) {
        fprintf(file_pointer, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                match->read_A_id_,
                match->read_B_id_, match->weight,
                direction, hinged,
                match->eff_read_A_match_start_,
                match->eff_read_A_match_end_,
                match->eff_read_B_match_start_,
                match->eff_read_B_match_end_,
                match->eff_read_A_start_,
                match->eff_read_A_end_,
                match->eff_read_B_start_,
                match->eff_read_B_end_);
    }
    else if ((match->match_type_ == BACKWARD_INTERNAL) or (match->match_type_ == BACKWARD)){
        fprintf(file_pointer, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                match->read_B_id_,
                match->read_A_id_,
                match->weight,
                direction, hinged,
                match->eff_read_B_match_start_,
                match->eff_read_B_match_end_,
                match->eff_read_A_match_start_,
                match->eff_read_A_match_end_,
                match->eff_read_B_start_,
                match->eff_read_B_end_,
                match->eff_read_A_start_,
                match->eff_read_A_end_);
    }
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
    std::cout << "In 308"<< std::endl;
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
    //std::vector<std::vector<LOverlap *>> idx2;
    /*
    	idx2 is a vector of length n_read, each element idx2[read A id] is a vector, 
    	for each read B, we put the best overlap into that vector
    */
    //std::vector<std::unordered_map<int, LOverlap *>> idx3;
	/*
		idx3 is a vector of length n_read, each element idx3[read A id] is a map, 
		from read read B id to the best overlap of read A and read B
	*/
    std::vector<std::vector<LOverlap *>> matches_forward, matches_backward;
    //matches_forward is the vector of vectors where matches_forward[read_id] is a vector of matches of read_id
    //of type FORWARD, and FORWARD_INTERNAL
    //matches_backward is the vector of vectors where matches_backward[read_id] is a vector of matches of read_id
    //of type BACKWARD, and BACKWARD_INTERNAL


    std::vector<std::vector<LOverlap *>> edges_forward, edges_backward; 
    // edges_forward is a "filtered" version of matches_forward, where every (active) read has at exactly one outgoing match
    // edges_backward is a "filtered" version of matches_backward, where every (active) read has at exactly one incoming match


    std::vector<std::vector<LOverlap *>> intersection_edges_forward, intersection_edges_backward; 
    //Stores the intersection of edges constructing the intersection list of edges


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

    for (int i = 0; i < n_read; i++) {
        //An initialisation for loop
        //TODO Preallocate memory. Much more efficient.
        idx.push_back(std::unordered_map<int, std::vector<LOverlap *> >());
        //idx2.push_back(std::vector<LOverlap *>());
        matches_forward.push_back(std::vector<LOverlap *>());
        matches_backward.push_back(std::vector<LOverlap *>());
        edges_forward.push_back(std::vector<LOverlap *>());
        edges_backward.push_back(std::vector<LOverlap *>());
        intersection_edges_forward.push_back(std::vector<LOverlap *>());
        intersection_edges_backward.push_back(std::vector<LOverlap *>());
    }

//int num_finished = 0;
    int num_overlaps = 0;
    int num_forward_overlaps(0),num_forward_internal_overlaps(0), num_reverse_overlaps(0),
            num_reverse_internal_overlaps(0), rev_complemented_matches(0);
//# pragma omp parallel for
    for (int i = 0; i < aln.size(); i++) {
        idx[aln[i]->read_A_id_][aln[i]->read_B_id_] = std::vector<LOverlap *>();
    }

    for (int i = 0; i < aln.size(); i++) {
        idx[aln[i]->read_A_id_][aln[i]->read_B_id_].push_back(aln[i]);
    }
    std::cout<<"index finished. " << "Number reads "<< n_read <<std::endl;

    for (int i = 0; i < n_read; i++) {
        bool contained=false;
        //std::cout<< "Testing opt " << i << std::endl;
        for (std::unordered_map<int, std::vector<LOverlap *> >::iterator it = idx[i].begin(); it!=idx[i].end(); it++) {
            std::sort(it->second.begin(), it->second.end(), compare_overlap);//Sort overlaps by lengths
            //std::cout<<"Giving input to ProcessAlignment "<<it->second.size() <<std::endl;
            if (it->second.size() > 0) {
                //Figure out if read is contained
                bool contained_alignment = ProcessAlignment(it->second[0],reads[it->second[0]->read_A_id_],
                                                         reads[it->second[0]->read_B_id_], ALN_THRESHOLD, THETA);
                contained=contained or contained_alignment;
                //Filter matches that matter.
                //TODO Figure out a way to do this more efficiently
                if ((it->second[0]->match_type_== FORWARD) or (it->second[0]->match_type_== FORWARD_INTERNAL))
                    matches_forward[i].push_back(it->second[0]);
                else if ((it->second[0]->match_type_== BACKWARD) or (it->second[0]->match_type_== BACKWARD_INTERNAL))
                    matches_backward[i].push_back(it->second[0]);
            }
        }
        if (contained) reads[i]->active = false;
    }



    for (int i = 0; i < n_read; i++) {//Isn't this just 0 or 1?
        num_overlaps += matches_forward[i].size()+ matches_backward[i].size();
        for (int j=0; j < matches_forward[i].size(); j++)
            rev_complemented_matches+= matches_forward[i][j]->reverse_complement_match_;
        for (int j=0; j < matches_backward[i].size(); j++)
            rev_complemented_matches+= matches_backward[i][j]->reverse_complement_match_;
    }
    std::cout<<num_overlaps << " overlaps" << std::endl;

    std::cout<<num_overlaps << " overlaps " << num_forward_overlaps << " fwd overlaps "
    << num_forward_internal_overlaps << " fwd internal overlaps "<< num_reverse_overlaps
    << " backward overlaps " << num_reverse_internal_overlaps
    << " backward internal overlaps "<< rev_complemented_matches << " reverse complement overlaps" << std::endl;

    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            num_active_read ++;
            maximal_reads << i << std::endl;
        }
    }
    std::cout<<"removed contained reads, active reads:" << num_active_read<< std::endl;

    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) num_active_read ++;
    }
    std::cout<<"active reads:" << num_active_read<< std::endl;

    num_overlaps = 0;
    num_forward_overlaps=0;
    num_forward_internal_overlaps=0;
    num_reverse_overlaps=0;
    num_reverse_internal_overlaps=0;
    rev_complemented_matches=0;
    int rev_complemented_fwd_matches(0),rev_complemented_bck_matches(0), rev_complemented_fwd_int_matches(0),
        rev_complemented_bck_int_matches(0);

    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (reads[matches_forward[i][j]->read_B_id_]->active) {
                    num_overlaps++;
                    if (matches_forward[i][j]->match_type_==FORWARD) {
                        num_forward_overlaps++;
                        rev_complemented_fwd_matches+=matches_forward[i][j]->reverse_complement_match_;
                    }
                    else if (matches_forward[i][j]->match_type_==FORWARD_INTERNAL) {
                        num_forward_internal_overlaps++;
                        rev_complemented_fwd_int_matches+=matches_forward[i][j]->reverse_complement_match_;
                    }
                    if (matches_forward[i][j]->reverse_complement_match_==1)
                        rev_complemented_matches++;

                }
            }
            //std::cout <<"First for done "<<std::endl;
            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (reads[matches_backward[i][j]->read_B_id_]->active) {
                    num_overlaps++;
                    if (matches_backward[i][j]->match_type_==BACKWARD) {
                        num_reverse_overlaps++;
                        rev_complemented_bck_matches+=matches_backward[i][j]->reverse_complement_match_;
                    }
                    else if (matches_backward[i][j]->match_type_==BACKWARD_INTERNAL) {
                        num_reverse_internal_overlaps++;
                        rev_complemented_bck_int_matches+=matches_backward[i][j]->reverse_complement_match_;
                    }
                    if (matches_backward[i][j]->reverse_complement_match_==1)
                        rev_complemented_matches++;
                }
            }
        }
    }
    std::cout<<num_overlaps << " overlaps " << num_forward_overlaps << " fwd overlaps "
    << num_forward_internal_overlaps << " fwd internal overlaps "<< num_reverse_overlaps
    << " backward overlaps " << num_reverse_internal_overlaps
    << " backward internal overlaps "<< rev_complemented_matches << " reverse complement overlaps\n"
    << rev_complemented_fwd_matches <<" rev cmplment fwd matches "
    << rev_complemented_fwd_int_matches << " rev cmplement fwd int matches "
    << rev_complemented_bck_matches << " rev cmplment bck matches "
    << rev_complemented_bck_int_matches << " rev cmplement bck int matches " << std::endl;



# pragma omp parallel for
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            std::sort(matches_forward[i].begin(), matches_forward[i].end(), compare_overlap_weight);
            std::sort(matches_backward[i].begin(), matches_backward[i].end(), compare_overlap_weight);
        }
    }

    FILE * out3;
    out3 = fopen("edges.fwd.backup.txt","w");
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (reads[matches_forward[i][j]->read_B_id_]->active)
                    fprintf(out3, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n",
                            matches_forward[i][j]->read_A_id_, matches_forward[i][j]->read_B_id_, 
                            matches_forward[i][j]->weight, matches_forward[i][j]->reverse_complement_match_, 
                            matches_forward[i][j]->match_type_, matches_forward[i][j]->eff_read_A_match_start_,
                            matches_forward[i][j]->eff_read_A_match_end_,
                            matches_forward[i][j]->eff_read_B_match_start_, matches_forward[i][j]->eff_read_B_match_end_,
                            matches_forward[i][j]->eff_read_A_start_, matches_forward[i][j]->eff_read_A_end_,
                            matches_forward[i][j]->eff_read_B_start_, matches_forward[i][j]->eff_read_B_end_);
            }
    }
    fclose(out3);
    out3 = fopen("edges.bkw.backup.txt","w");
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (reads[matches_backward[i][j]->read_B_id_]->active)
                    fprintf(out3, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n",
                            matches_backward[i][j]->read_A_id_, matches_backward[i][j]->read_B_id_,
                            matches_backward[i][j]->weight, matches_backward[i][j]->reverse_complement_match_,
                            matches_backward[i][j]->match_type_, matches_backward[i][j]->eff_read_A_match_start_,
                            matches_backward[i][j]->eff_read_A_match_end_,
                            matches_backward[i][j]->eff_read_B_match_start_, matches_backward[i][j]->eff_read_B_match_end_,
                            matches_backward[i][j]->eff_read_A_start_, matches_backward[i][j]->eff_read_A_end_,
                            matches_backward[i][j]->eff_read_B_start_, matches_backward[i][j]->eff_read_B_end_);
            }
    }
    fclose(out3);



    FILE * out;
    FILE * out2;
    FILE * out4;
    out = fopen((std::string(argv[3]) + ".1").c_str(), "w");
    out2 = fopen((std::string(argv[3]) + ".2").c_str(), "w");

    // Output file for matches 
    out3 = fopen((std::string(argv[3]) + ".hinges").c_str(), "w");

    // Output file for edges
    out4 = fopen((std::string(argv[3]) + ".hinges.edges").c_str(),"w");



    std::unordered_map<int, std::vector<Hinge> > hinges_vec;



    int n = 0;
    for (int i = 0; i < n_read; i++) {
        hinges_vec[i] = std::vector<Hinge>();
        for (int j = 0; j < marked_hinges[i].size(); j++) {
            hinges_vec[i].push_back(Hinge(marked_hinges[i][j].first, marked_hinges[i][j].second , true));
            if (reads[i]->active) {
                n ++;
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
        //This is in essence the filtering step
        //For each read find the best forward match, and remove all incoming hinges starting after the start
        //of the match corresponding to this.
        if (reads[i]->active) {
            int forward = 0;
            int backward = 0;
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (matches_forward[i][j]->active) {
                    if ((matches_forward[i][j]->match_type_ == FORWARD) and
                            (reads[matches_forward[i][j]->read_B_id_]->active)) {
                        if (forward < 1) {
                            //remove certain hinges
                            for (int k = 0; k < hinges_vec[i].size(); k++) {
                                if ((matches_forward[i][j]->eff_read_A_match_start_ < hinges_vec[i][k].pos - 400)
                                    and (hinges_vec[i][k].type == 1))
                                    hinges_vec[i][k].active = false;
                            }
                        }
                        forward++;
                    }
                }
            }

            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (matches_backward[i][j]->active) {
                    if ((matches_backward[i][j]->match_type_ == BACKWARD) and
                            (reads[matches_backward[i][j]->read_B_id_]->active)) {
                        if (backward < 1) {
                            //remove certain hinges
                            for (int k = 0; k < hinges_vec[i].size(); k++) {
                                if ((matches_backward[i][j]->eff_read_A_match_start_ < hinges_vec[i][k].pos - 400)
                                    and (hinges_vec[i][k].type == -1))
                                    hinges_vec[i][k].active = false;
                            }
                        }
                        backward++;
                    }
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



    //Perform greedy graph construction and write outputs out and out2
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            int forward = 0;
            int backward = 0;
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (matches_forward[i][j]->active) {
                    if ((matches_forward[i][j]->match_type_ == FORWARD) and
                        (reads[matches_forward[i][j]->read_B_id_]->active)) {
                        /*if (not repeat_status_back[i])*/
                        {
                            if (forward < 1) {

                                if (matches_forward[i][j]->reverse_complement_match_ == 0)
                                    fprintf(out, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                                            matches_forward[i][j]->read_A_id_,
                                            matches_forward[i][j]->read_B_id_, matches_forward[i][j]->weight,
                                            matches_forward[i][j]->eff_read_A_match_start_,
                                            matches_forward[i][j]->eff_read_A_match_end_,
                                            matches_forward[i][j]->eff_read_B_match_start_,
                                            matches_forward[i][j]->eff_read_B_match_end_,
                                            matches_forward[i][j]->eff_read_A_start_,
                                            matches_forward[i][j]->eff_read_A_end_,
                                            matches_forward[i][j]->eff_read_B_start_,
                                            matches_forward[i][j]->eff_read_B_end_);
                                else
                                    fprintf(out, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                                            matches_forward[i][j]->read_A_id_,
                                            matches_forward[i][j]->read_B_id_, matches_forward[i][j]->weight,
                                            matches_forward[i][j]->eff_read_A_match_start_,
                                            matches_forward[i][j]->eff_read_A_match_end_,
                                            matches_forward[i][j]->eff_read_B_match_start_,
                                            matches_forward[i][j]->eff_read_B_match_end_,
                                            matches_forward[i][j]->eff_read_A_start_,
                                            matches_forward[i][j]->eff_read_A_end_,
                                            matches_forward[i][j]->eff_read_B_start_,
                                            matches_forward[i][j]->eff_read_B_end_);

                                if (matches_forward[i][j]->reverse_complement_match_ == 0)
                                    fprintf(out2, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                                            matches_forward[i][j]->read_B_id_,
                                            matches_forward[i][j]->read_A_id_, matches_forward[i][j]->weight,
                                            matches_forward[i][j]->eff_read_A_match_start_,
                                            matches_forward[i][j]->eff_read_A_match_end_,
                                            matches_forward[i][j]->eff_read_B_match_start_,
                                            matches_forward[i][j]->eff_read_B_match_end_,
                                            matches_forward[i][j]->eff_read_A_start_,
                                            matches_forward[i][j]->eff_read_A_end_,
                                            matches_forward[i][j]->eff_read_B_start_,
                                            matches_forward[i][j]->eff_read_B_end_);
                                else
                                    fprintf(out2, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                                            matches_forward[i][j]->read_B_id_,
                                            matches_forward[i][j]->read_A_id_, matches_forward[i][j]->weight,
                                            matches_forward[i][j]->eff_read_A_match_start_,
                                            matches_forward[i][j]->eff_read_A_match_end_,
                                            matches_forward[i][j]->eff_read_B_match_start_,
                                            matches_forward[i][j]->eff_read_B_match_end_,
                                            matches_forward[i][j]->eff_read_A_start_,
                                            matches_forward[i][j]->eff_read_A_end_,
                                            matches_forward[i][j]->eff_read_B_start_,
                                            matches_forward[i][j]->eff_read_B_end_);

                            }
                        }
                        forward++;
                    }
                }
            }
            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (matches_backward[i][j]->active) {
                    if ((matches_backward[i][j]->match_type_ == BACKWARD) and
                        (reads[matches_backward[i][j]->read_B_id_]->active)) {
                        /*if (not repeat_status_back[i])*/
                        {
                            if (backward < 1) {

                                if (matches_backward[i][j]->reverse_complement_match_ == 0)
                                    fprintf(out, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                                            matches_backward[i][j]->read_A_id_,
                                            matches_backward[i][j]->read_B_id_, matches_backward[i][j]->weight,
                                            matches_backward[i][j]->eff_read_A_match_start_,
                                            matches_backward[i][j]->eff_read_A_match_end_,
                                            matches_backward[i][j]->eff_read_B_match_start_,
                                            matches_backward[i][j]->eff_read_B_match_end_,
                                            matches_backward[i][j]->eff_read_A_start_,
                                            matches_backward[i][j]->eff_read_A_end_,
                                            matches_backward[i][j]->eff_read_B_start_,
                                            matches_backward[i][j]->eff_read_B_end_);
                                else
                                    fprintf(out, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                                            matches_backward[i][j]->read_A_id_,
                                            matches_backward[i][j]->read_B_id_, matches_backward[i][j]->weight,
                                            matches_backward[i][j]->eff_read_A_match_start_,
                                            matches_backward[i][j]->eff_read_A_match_end_,
                                            matches_backward[i][j]->eff_read_B_match_start_,
                                            matches_backward[i][j]->eff_read_B_match_end_,
                                            matches_backward[i][j]->eff_read_A_start_,
                                            matches_backward[i][j]->eff_read_A_end_,
                                            matches_backward[i][j]->eff_read_B_start_,
                                            matches_backward[i][j]->eff_read_B_end_);

                                if (matches_backward[i][j]->reverse_complement_match_ == 0)
                                    fprintf(out2, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                                            matches_backward[i][j]->read_B_id_,
                                            matches_backward[i][j]->read_A_id_, matches_backward[i][j]->weight,
                                            matches_backward[i][j]->eff_read_A_match_start_,
                                            matches_backward[i][j]->eff_read_A_match_end_,
                                            matches_backward[i][j]->eff_read_B_match_start_,
                                            matches_backward[i][j]->eff_read_B_match_end_,
                                            matches_backward[i][j]->eff_read_A_start_,
                                            matches_backward[i][j]->eff_read_A_end_,
                                            matches_backward[i][j]->eff_read_B_start_,
                                            matches_backward[i][j]->eff_read_B_end_);
                                else
                                    fprintf(out2, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                                            matches_backward[i][j]->read_B_id_,
                                            matches_backward[i][j]->read_A_id_, matches_backward[i][j]->weight,
                                            matches_backward[i][j]->eff_read_A_match_start_,
                                            matches_backward[i][j]->eff_read_A_match_end_,
                                            matches_backward[i][j]->eff_read_B_match_start_,
                                            matches_backward[i][j]->eff_read_B_match_end_,
                                            matches_backward[i][j]->eff_read_A_start_,
                                            matches_backward[i][j]->eff_read_A_end_,
                                            matches_backward[i][j]->eff_read_B_start_,
                                            matches_backward[i][j]->eff_read_B_end_);

                            }
                        }
                        backward++;
                    }
                }
            }
        }
    }

    num_overlaps = 0;
    num_forward_overlaps=0;
    num_forward_internal_overlaps=0;
    num_reverse_overlaps=0;
    num_reverse_internal_overlaps=0;
    rev_complemented_matches=0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (reads[matches_forward[i][j]->read_B_id_]->active) {
                    num_overlaps++;
                    if (matches_forward[i][j]->match_type_==FORWARD)
                        num_forward_overlaps++;
                    else if (matches_forward[i][j]->match_type_==FORWARD_INTERNAL)
                        num_forward_internal_overlaps++;
                    if (matches_forward[i][j]->reverse_complement_match_==1)
                        rev_complemented_matches++;

                }
            }
            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (reads[matches_backward[i][j]->read_B_id_]->active) {
                    num_overlaps++;
                    if (matches_backward[i][j]->match_type_==BACKWARD)
                        num_reverse_overlaps++;
                    else if (matches_backward[i][j]->match_type_==BACKWARD_INTERNAL)
                        num_reverse_internal_overlaps++;
                    if (matches_backward[i][j]->reverse_complement_match_==1)
                        rev_complemented_matches++;
                }
            }
        }
    }
    std::cout<<num_overlaps << " overlaps " << num_forward_overlaps << " fwd overlaps "
    << num_forward_internal_overlaps << " fwd internal overlaps "<< num_reverse_overlaps
    << " backward overlaps " << num_reverse_internal_overlaps
    << " backward internal overlaps "<< rev_complemented_matches << " reverse complement overlaps" << std::endl;


    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            for (int j = 0; j < matches_forward[i].size(); j++){
                if (matches_forward[i][j]->active) {
                    if ((reads[matches_forward[i][j]->read_B_id_]->active)) {
                        if ((matches_forward[i][j]->match_type_ == FORWARD)){
                            fprintf(out3,"Printed from forward\n");
                            PrintOverlapToFile(out3,matches_forward[i][j]);
                            edges_forward[i].push_back(matches_forward[i][j]);

                            break;

                        }
                        else if ((matches_forward[i][j]->match_type_ == FORWARD_INTERNAL)
                                and isValidHinge(matches_forward[i][j], hinges_vec[matches_forward[i][j]->read_B_id_])){
                            fprintf(out3,"Printed from forward internal\n");
                            PrintOverlapToFile(out3,matches_forward[i][j]);
                            edges_forward[i].push_back(matches_forward[i][j]);
                            break;

                        }
                    }
                }
            }
            for (int j = 0; j < matches_backward[i].size(); j++){
                if (matches_backward[i][j]->active) {
                    if ((reads[matches_backward[i][j]->read_B_id_]->active)) {
                        if ((matches_backward[i][j]->match_type_ == BACKWARD)){
                            fprintf(out3,"Printed from backward\n");
                            PrintOverlapToFile(out3,matches_backward[i][j]);
                            edges_backward[i].push_back(matches_backward[i][j]);
                            break;

                        }
                        else if ((matches_backward[i][j]->match_type_ == BACKWARD_INTERNAL)
                                 and isValidHinge(matches_backward[i][j], hinges_vec[matches_backward[i][j]->read_B_id_])) {
                            fprintf(out3,"Printed from backward internal\n");
                            PrintOverlapToFile(out3,matches_backward[i][j]);
                            edges_backward[i].push_back(matches_backward[i][j]);
                            break;

                        }
                    }
                }
            }
        }
    }


    // Find intersection forward edges
    // For each edge in edges_forward, push it into intersection_edges_forward if:
    // - It ends in a hinge, or
    // - It doesn't end in a hinge, but it is picked as a reverse edges as well
    for (int i=0; i < edges_forward.size(); i++) {
        if (edges_forward[i].size() > 0) {  // size should be either zero or one
            if (edges_forward[i][0]->match_type_ == FORWARD_INTERNAL) {
                intersection_edges_forward[i].push_back(edges_forward[i][0]);
            }
            else { // match_type_ should be FORWARD

                int read_B_id = edges_forward[i][0]->read_B_id_;

                if (edges_forward[i][0]->reverse_complement_match_ != 1) {
                    if (edges_backward[read_B_id].size() > 0) {
                        if ((edges_backward[read_B_id][0]->match_type_ == BACKWARD) and
                                (edges_backward[read_B_id][0]->read_B_id_ == i)) {
                            PrintOverlapToFile(out4,edges_forward[i][0]);
                            intersection_edges_forward[i].push_back(edges_forward[i][0]);
                        }
                    }
                }
                else { // reverse complement match
                    if (edges_forward[read_B_id].size() > 0) {
                        if ((edges_forward[read_B_id][0]->match_type_ == FORWARD) and
                                (edges_forward[read_B_id][0]->read_B_id_ == i)) {
                            PrintOverlapToFile(out4,edges_forward[i][0]);
                            intersection_edges_forward[i].push_back(edges_forward[i][0]);
                        }
                    }
                }
            }
        }
    }




    // Find intersection backward edges
    // For each edge in edges_forward, push it into intersection_edges_forward if:
    // - It ends in a hinge, or
    // - It doesn't end in a hinge, but it is picked as a reverse edges as well
    for (int i=0; i < edges_backward.size(); i++) {
        if (edges_backward[i].size() > 0) {  // size should be either zero or one
            if (edges_backward[i][0]->match_type_ == BACKWARD_INTERNAL) {
                intersection_edges_backward[i].push_back(edges_backward[i][0]);
            }
            else { // match_type_ should be FORWARD

                int read_B_id = edges_backward[i][0]->read_B_id_;

                if (edges_backward[i][0]->reverse_complement_match_ != 1) {
                    if (edges_forward[read_B_id].size() > 0) {
                        if ((edges_forward[read_B_id][0]->match_type_ == FORWARD) and
                                (edges_forward[read_B_id][0]->read_B_id_ == i)) {
                            PrintOverlapToFile(out4,edges_backward[i][0]);
                            intersection_edges_backward[i].push_back(edges_backward[i][0]);
                        }
                    }
                }
                else { // reverse complement match
                    if (edges_backward[read_B_id].size() > 0) {
                        if ((edges_backward[read_B_id][0]->match_type_ == BACKWARD) and
                                (edges_backward[read_B_id][0]->read_B_id_ == i)) {
                            PrintOverlapToFile(out4,edges_backward[i][0]);
                            intersection_edges_backward[i].push_back(edges_backward[i][0]);
                        }
                    }
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
