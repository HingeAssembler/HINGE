#include <stdio.h>
#include <unistd.h>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <omp.h>
#include <tuple>
#include <iomanip>

#include "spdlog/spdlog.h"
#include "cmdline.h"
#include "INIReader.h"
#include "DB.h"
#include "align.h"
#include "LAInterface.h"

#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#define LAST_READ_SYMBOL  '$'

#define HINGED_EDGE 1
#define UNHINGED_EDGE -1
#define REVERSE_COMPLEMENT_MATCH 1
#define SAME_DIRECTION_MATCH 0

using namespace boost;

typedef adjacency_list <vecS, vecS, undirectedS> Graph;
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

bool ProcessAlignment(LOverlap * match, Read * read_A, Read * read_B, int ALN_THRESHOLD,
                      int THETA, int THETA2, bool trim){
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

    // removed the following if, so that things agree with the convention for reverse complement matches

    match->eff_read_B_start_ = read_B->effective_start;
    match->eff_read_B_end_ = read_B->effective_end;

//    if (match->reverse_complement_match_ == 0) {
//        match->eff_read_B_start_ = read_B->effective_start;
//        match->eff_read_B_end_ = read_B->effective_end;
//    } else {
//        match->eff_read_B_start_ = read_B->len - read_B->effective_end;
//        match->eff_read_B_end_ = read_B->len - read_B->effective_start;
//    }

    /*printf("bef %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", match->read_A_id_, match->read_B_id_,
     * match->reverse_complement_match_,
        match->read_A_match_start_, match->read_A_match_end_, match->read_B_match_start_, match->read_B_match_end_,
           match->eff_read_A_start_, match->eff_read_A_end_, match->eff_read_B_start_, match->eff_read_B_end_
    );*/

    if (trim)
        match->trim_overlap();
    else {
        match->eff_read_B_match_start_ = match->read_B_match_start_;
        match->eff_read_B_match_end_ = match->read_B_match_end_;
        match->eff_read_A_match_start_ = match->read_A_match_start_;
        match->eff_read_A_match_end_ = match->read_A_match_end_;
    }
    /*printf("aft %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", match->read_A_id_, match->read_B_id_,
     * match->reverse_complement_match_,
           match->eff_read_A_match_start_, match->eff_read_A_match_end_, match->eff_read_B_match_start_,
           match->eff_read_B_match_end_,
           match->eff_read_A_start_, match->eff_read_A_end_, match->eff_read_B_start_, match->eff_read_B_end_
    );*/
    //std::cout<< contained<<std::endl;
    if (((match->eff_read_B_match_end_ - match->eff_read_B_match_start_) < ALN_THRESHOLD)
        or ((match->eff_read_A_match_end_ - match->eff_read_A_match_start_) < ALN_THRESHOLD) or (!match->active))

    {
        match->active = false;
        match->match_type_ = NOT_ACTIVE;
    } else {
        match->AddTypesAsymmetric(THETA,THETA2);
        if (match->match_type_ == BCOVERA) {
            contained = true;
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
    Hinge(int pos, int t, bool active):pos(pos),type(t), active(active) {};
    Hinge():pos(0),type(1), active(true) {};
};

// if we uncomment this, we need to make sure it works with the new convention of B_match_start and
// B_match_end for reverse complement matches

//bool isValidHinge(LOverlap *match, std::vector<Hinge> &read_hinges){
//    //Returns true if read_hinges (a vector of all hinges corresponding to a read )
//    // has a hinge of appropriate type within tolerance from positions of start of the
//    // overlap on read_B of the overlap given.
//    int tolerance=100;//TODO put as #define
//    int position=match->eff_read_B_match_start_;   // parei aqui
//    int type; //TODO : Make enum
//    if (match->match_type_==FORWARD_INTERNAL)
//        type=1;
//    else if (match->match_type_==BACKWARD_INTERNAL)
//        type=-1;
//
//    if (match->reverse_complement_match_==1){
//        type=-type;
//        position=match->eff_read_B_match_end_;
//    }
//
//    bool valid=false;
//    for (int index=0; index < read_hinges.size(); index++) {
//        if ((abs(position - read_hinges[index].pos) < tolerance) and (type == read_hinges[index].type))
//            valid = true;
//        return valid;
//    }
//}



void PrintOverlapToFile(FILE * file_pointer, LOverlap * match) {

    int direction = match->reverse_complement_match_;
    int hinged;

    if ((match->match_type_ == FORWARD) or (match->match_type_ == BACKWARD))
        hinged = UNHINGED_EDGE;

    else if ((match->match_type_ == FORWARD_INTERNAL) or (match->match_type_ == BACKWARD_INTERNAL))
        hinged = HINGED_EDGE;

    if ((match->match_type_ == FORWARD_INTERNAL) or (match->match_type_ == FORWARD)) {
        fprintf(file_pointer, "%d %d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                match->read_A_id_,
                match->read_B_id_,
                match->weight,
                0,
                direction,
                hinged,
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
        fprintf(file_pointer, "%d %d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                match->read_B_id_,
                match->read_A_id_,
                match->weight,
                direction,
                0,
                hinged,
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




void PrintOverlapToFile2(FILE * file_pointer, LOverlap * match, int hinge_pos) {

    int direction = match->reverse_complement_match_;
    int hinged;

//    if ((match->match_type_ == FORWARD) or (match->match_type_ == BACKWARD))
//        hinged = UNHINGED_EDGE;
//
//    else if ((match->match_type_ == FORWARD_INTERNAL) or (match->match_type_ == BACKWARD_INTERNAL))
//        hinged = HINGED_EDGE;

//    if ((match->match_type_ == FORWARD) or (match->match_type_ == BACKWARD))
//        hinged = 0;
//    else if (match->match_type_ == FORWARD_INTERNAL)
//        hinged = 1;
//    else if (match->match_type_ == BACKWARD_INTERNAL)
//        hinged = -1;

    if (match->match_type_ == FORWARD) {
        fprintf(file_pointer, "%d %d %d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                match->read_A_id_,
                match->read_B_id_,
                match->weight,
                0,
                direction,
                0,
                -1, // hinge pos
                match->eff_read_A_match_start_,
                match->eff_read_A_match_end_,
                match->eff_read_B_match_start_,
                match->eff_read_B_match_end_,
                match->eff_read_A_start_,
                match->eff_read_A_end_,
                match->eff_read_B_start_,
                match->eff_read_B_end_);
    }
    else if (match->match_type_ == BACKWARD) {
        fprintf(file_pointer, "%d %d %d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                match->read_B_id_,
                match->read_A_id_,
                match->weight,
                direction,
                0,
                0,
                -1, // hinge pos
                match->eff_read_B_match_start_,
                match->eff_read_B_match_end_,
                match->eff_read_A_match_start_,
                match->eff_read_A_match_end_,
                match->eff_read_B_start_,
                match->eff_read_B_end_,
                match->eff_read_A_start_,
                match->eff_read_A_end_);
    }
    else if (match->match_type_ == FORWARD_INTERNAL) {

        fprintf(file_pointer, "%d %d %d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                match->read_A_id_,
                match->read_B_id_,
                match->weight,
                0,
                direction,
                1, // hinged forward
                hinge_pos,
                match->eff_read_A_match_start_,
                match->eff_read_A_match_end_,
                match->eff_read_B_match_start_,
                match->eff_read_B_match_end_,
                match->eff_read_A_start_,
                match->eff_read_A_end_,
                match->eff_read_B_start_,
                match->eff_read_B_end_);
    }
    else if (match->match_type_ == BACKWARD_INTERNAL) {
        fprintf(file_pointer, "%d %d %d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
                match->read_B_id_,
                match->read_A_id_,
                match->weight,
                direction,
                0,
                -1, // hinged backward
                hinge_pos,
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

    cmdline::parser cmdp;
    cmdp.add<std::string>("db", 'b', "db file name", false, "");
    cmdp.add<std::string>("las", 'l', "las file name", false, "");
    cmdp.add<std::string>("paf", 'p', "paf file name", false, "");
    cmdp.add<std::string>("config", 'c', "configuration file name", false, "");
    cmdp.add<std::string>("fasta", 'f', "fasta file name", false, "");
    cmdp.add<std::string>("prefix", 'x', "(intermediate output) input file prefix", true, "");
    cmdp.add<std::string>("out", 'o', "final output file name", true, "");
    cmdp.add<std::string>("log", 'g', "log folder name", false, "log");



//    cmdp.add<std::string>("restrictreads",'r',"restrict to reads in the file",false,"");


    cmdp.parse_check(argc, argv);

    LAInterface la;
    const char *name_db = cmdp.get<std::string>("db").c_str(); //.db file of reads to load
    const char *name_las = cmdp.get<std::string>("las").c_str();//.las file of alignments
    const char *name_paf = cmdp.get<std::string>("paf").c_str();
    const char *name_fasta = cmdp.get<std::string>("fasta").c_str();
    const char *name_config = cmdp.get<std::string>("config").c_str();//name of the configuration file, in INI format
    std::string out = cmdp.get<std::string>("prefix");
    std::string out_name = cmdp.get<std::string>("out");
//    const char * name_restrict = cmdp.get<std::string>("restrictreads").c_str();


    std::string name_mask = out + ".mas";
    std::string name_max = out + ".max";
    std::string name_homo = out + ".homologous.txt";
    std::string name_rep = out + ".repeat.txt";
    std::string name_hg = out + ".hinges.txt";
    std::string name_cov = out + ".coverage.txt";
    std::string name_garbage = out + ".garbage.txt";
    std::string name_contained = out + ".contained.txt";
    std::string name_deadend = out_name + ".deadends.txt";


    std::ofstream deadend_out(name_deadend);
    std::ofstream maximal_reads(name_max);
    std::ofstream garbage_out(name_garbage);
    std::ofstream contained_out(name_contained);
    std::ifstream homo(name_homo);
    std::vector<int> homo_reads;


    bool delete_telomere = false;  // TODO: command line option to set this true

    int read_id;
    while (homo >> read_id) homo_reads.push_back(read_id);


    namespace spd = spdlog;

    //auto console = spd::stdout_logger_mt("console");
    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
    sinks.push_back(
            std::make_shared<spdlog::sinks::daily_file_sink_st>(cmdp.get<std::string>("log") + "/log", "txt", 23, 59));
    auto console = std::make_shared<spdlog::logger>("log", std::begin(sinks), std::end(sinks));
    spdlog::register_logger(console);

    console->info("Hinging layout");
    char *buff = (char *) malloc(sizeof(char) * 2000);
    getwd(buff);
    console->info("current user {}, current working directory {}", getlogin(), buff);
    free(buff);
    console->info("name of db: {}, name of .las file {}", name_db, name_las);
    console->info("name of fasta: {}, name of .paf file {}", name_fasta, name_paf);
    console->info("filter files prefix: {}", out);
    console->info("output prefix: {}", out_name);


    std::ifstream ini_file(name_config);
    std::string str((std::istreambuf_iterator<char>(ini_file)),
                    std::istreambuf_iterator<char>());

    console->info("Parameters passed in \n{}", str);

    if (strlen(name_db) > 0)
        la.openDB(name_db);


    if (strlen(name_las) > 0)
        la.openAlignmentFile(name_las);

    int64 n_aln = 0;

    if (strlen(name_las) > 0) {
        n_aln = la.getAlignmentNumber();
        console->info("Load alignments from {}", name_las);
        console->info("# Alignments: {}", n_aln);
    }

    int n_read;
    if (strlen(name_db) > 0)
        n_read = la.getReadNumber();

    std::vector<Read *> reads; //Vector of pointers to all reads

    if (strlen(name_fasta) > 0) {
        n_read = la.loadFASTA(name_fasta, reads);
    }

    console->info("# Reads: {}", n_read); // output some statistics

    std::vector<LOverlap *> aln;//Vector of pointers to all alignments

    if (strlen(name_las) > 0) {
        la.resetAlignment();
        la.getOverlap(aln, 0, n_aln);
    }

    if (strlen(name_paf) > 0) {
        n_aln = la.loadPAF(std::string(name_paf), aln);
        console->info("Load alignments from {}", name_paf);
        console->info("# Alignments: {}", n_aln);
    }

    if (n_aln == 0) {
        console->error("No alignments!");
        return 1;
    }


    if (strlen(name_db) > 0) {
        la.getRead(reads, 0, n_read);
    }

    console->info("Input data finished");

    INIReader reader(name_config);

    if (reader.ParseError() < 0) {
        console->warn("Can't load {}", name_config);
        return 1;
    }

    int LENGTH_THRESHOLD = int(reader.GetInteger("filter", "length_threshold", -1));
    double QUALITY_THRESHOLD = reader.GetReal("filter", "quality_threshold", 0.0);
    int N_ITER = (int) reader.GetInteger("filter", "n_iter", -1);
    int ALN_THRESHOLD = (int) reader.GetInteger("filter", "aln_threshold", -1);
    int MIN_COV = (int) reader.GetInteger("filter", "min_cov", -1);
    int CUT_OFF = (int) reader.GetInteger("filter", "cut_off", -1);
    int THETA = (int) reader.GetInteger("filter", "theta", -1);
    int THETA2 = (int) reader.GetInteger("filter", "theta2", 0);
    int N_PROC = (int) reader.GetInteger("running", "n_proc", 4);
    int HINGE_SLACK = (int) reader.GetInteger("layout", "hinge_slack", 1000);
    //This is the amount by which  a forward overlap
    //must be longer than a forward internal overlap to be preferred while
    //building a graph.
    int HINGE_TOLERANCE = (int) reader.GetInteger("layout", "hinge_tolerance", 150);
    //This is how far an overlap must start from a hinge to be considered an internal
    //overlap.
    int KILL_HINGE_OVERLAP_ALLOWANCE = (int) reader.GetInteger("layout", "kill_hinge_overlap", 300);
    int KILL_HINGE_INTERNAL_ALLOWANCE = (int) reader.GetInteger("layout", "kill_hinge_internal", 40);

    int MATCHING_HINGE_SLACK = (int) reader.GetInteger("layout", "matching_hinge_slack", 200);

    int NUM_EVENTS_TELOMERE = (int) reader.GetInteger("layout", "num_events_telomere", 7);

    int MIN_CONNECTED_COMPONENT_SIZE = (int) reader.GetInteger("layout", "min_connected_component_size", 8);


    console->info("LENGTH_THRESHOLD = {}", LENGTH_THRESHOLD);
    console->info("QUALITY_THRESHOLD = {}", QUALITY_THRESHOLD);
    console->info("ALN_THRESHOLD = {}", ALN_THRESHOLD);
    console->info("MIN_COV = {}", MIN_COV);
    console->info("CUT_OFF = {}", CUT_OFF);
    console->info("THETA = {}", THETA);
    console->info("N_ITER = {}", N_ITER);
    console->info("THETA2 = {}", THETA2);
    console->info("N_PROC = {}", N_PROC);
    console->info("HINGE_SLACK = {}", HINGE_SLACK);
    console->info("HINGE_TOLERANCE = {}", HINGE_TOLERANCE);
    console->info("KILL_HINGE_OVERLAP_ALLOWANCE = {}", KILL_HINGE_OVERLAP_ALLOWANCE);
    console->info("KILL_HINGE_INTERNAL_ALLOWANCE = {}", KILL_HINGE_INTERNAL_ALLOWANCE);
    console->info("MATCHING_HINGE_SLACK = {}", MATCHING_HINGE_SLACK);
    console->info("MIN_CONNECTED_COMPONENT_SIZE = {}", MIN_CONNECTED_COMPONENT_SIZE);






    omp_set_num_threads(N_PROC);
    //std::vector< std::vector<std::vector<LOverlap*>* > > idx2(n_read);
    // unordered_map from (aid) to alignments in a vector
    std::vector<Edge_w> edgelist, edgelist_ms; // save output to edgelist
    //std::unordered_map<int, std::vector <LOverlap * > >idx3,idx4;
    // this is the pileup
    std::vector<std::unordered_map<int, std::vector<LOverlap *> > > idx_ab;
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
    
    if (strlen(name_db) > 0)
    la.closeDB(); //close database
    return 0;
}
