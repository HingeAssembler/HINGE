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
#include "OverlapGraph.h"

//#include <ogdf/basic/Graph.h>
//#include <ogdf/basic/graph_generators.h>
//#include <ogdf/basic/Hashing.h>
//#include <ogdf/basic/simple_graph_alg.h>
//#include <ogdf/basic/HashArray.h>

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
        fprintf(file_pointer, "%d %d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] [%d %d] [%d %d]\n",
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
                match->eff_read_B_end_,

                match->read_A_match_start_,
                match->read_A_match_end_,
                match->read_B_match_start_,
                match->read_B_match_end_


                );
    }
    else if ((match->match_type_ == BACKWARD_INTERNAL) or (match->match_type_ == BACKWARD)){
        fprintf(file_pointer, "%d %d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] [%d %d] [%d %d]\n",
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
                match->eff_read_A_end_,

                match->read_A_match_start_,
                match->read_A_match_end_,
                match->read_B_match_start_,
                match->read_B_match_end_

        );
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
    std::vector<std::vector<LOverlap *>> matches_forward, matches_backward;
    //matches_forward is the vector of vectors where matches_forward[read_id] is a vector of matches of read_id
    //of type FORWARD, and FORWARD_INTERNAL
    //matches_backward is the vector of vectors where matches_backward[read_id] is a vector of matches of read_id
    //of type BACKWARD, and BACKWARD_INTERNAL


    std::vector<std::vector<LOverlap *>> edges_forward, edges_backward;
    // edges_forward is a "filtered" version of matches_forward, where every (active) read has at exactly
    // one outgoing match
    // edges_backward is a "filtered" version of matches_backward, where every (active) read has at exactly
    // one incoming match


    std::vector<std::vector<LOverlap *>> intersection_edges_forward, intersection_edges_backward;
    //Stores the intersection of edges constructing the intersection list of edges


    FILE *mask_file;
    mask_file = fopen(name_mask.c_str(), "r");
    int read, rs, re;

    while (fscanf(mask_file, "%d %d %d", &read, &rs, &re) != EOF) {
        reads[read]->effective_start = rs;
        reads[read]->effective_end = re;
    }
    console->info("read mask finished");

    FILE *repeat_file;
    repeat_file = fopen(name_rep.c_str(), "r");
    FILE *hinge_file;
    hinge_file = fopen(name_hg.c_str(), "r");
    char *line = NULL;
    size_t len = 0;
    std::unordered_map<int, std::vector<std::pair<int, int>>> marked_repeats;

    int telomere_cnt = 0;

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
            if ((r1 != 0) and (r2 != 0)) {
                //printf("[%d %d]\n", r1, r2);
                marked_repeats[num].push_back(std::pair<int, int>(r1, r2));
            }
        }
        ss.clear();

        if ((delete_telomere) and (marked_repeats[num].size() > NUM_EVENTS_TELOMERE)) {
            reads[num]->active = false;
            telomere_cnt++;
        }

    }
    fclose(repeat_file);
    console->info("read marked repeats");
    console->info("killed {} reads with many repeats",telomere_cnt);

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
            r1 = 0;
            r2 = 0;
            ss >> r1 >> r2;
            if ((r1 != 0) and (r2 != 0)) {
                //printf("[%d %d]\n", r1, r2);
                marked_hinges[num].push_back(std::pair<int, int>(r1, r2));
            }
        }
        ss.clear();
    }
    fclose(hinge_file);

    console->info("read marked hinges");

    if (line)
        free(line);

    int num_active_read = 0;

    //This seems to be an unnecessary stub
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) num_active_read++;
    }
    console->info("active reads: {}", num_active_read);


    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->effective_end - reads[i]->effective_start < LENGTH_THRESHOLD) {
            reads[i]->active = false;
            garbage_out << i << std::endl;
        }
        else num_active_read++;
    }
    console->info("active reads: {}", num_active_read);

    for (int i = 0; i < n_read; i++) {
        //An initialisation for loop
        //TODO Preallocate memory. Much more efficient.
        idx_ab.push_back(std::unordered_map<int, std::vector<LOverlap *> >());
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
    int num_forward_overlaps(0), num_forward_internal_overlaps(0), num_reverse_overlaps(0),
            num_reverse_internal_overlaps(0), rev_complemented_matches(0);
//# pragma omp parallel for
    for (int i = 0; i < aln.size(); i++) {
        idx_ab[aln[i]->read_A_id_][aln[i]->read_B_id_] = std::vector<LOverlap *>();
    }

    for (int i = 0; i < aln.size(); i++) {
        idx_ab[aln[i]->read_A_id_][aln[i]->read_B_id_].push_back(aln[i]);
    }

    int n_overlaps = 0;
    int n_rev_overlaps = 0;
    for (int i = 0; i < aln.size(); i++) {
        n_overlaps++;
        n_rev_overlaps += aln[i]->reverse_complement_match_;
    }

    console->info("overlaps {} rev_overlaps {}", n_overlaps, n_rev_overlaps);

    console->info("index finished");
    console->info("Number reads {}", n_read);


    for (int i = 0; i < n_read; i++) {
        bool contained = false;
        //std::cout<< "Testing opt " << i << std::endl;
        if (reads[i]->active == false) {
            continue;
        }

        int containing_read;

        for (std::unordered_map<int, std::vector<LOverlap *> >::iterator it = idx_ab[i].begin();
             it != idx_ab[i].end(); it++) {
            std::sort(it->second.begin(), it->second.end(), compare_overlap);//Sort overlaps by lengths
            //std::cout<<"Giving input to ProcessAlignment "<<it->second.size() <<std::endl;

            if (it->second.size() > 0) {
                //Figure out if read is contained
                LOverlap *ovl = it->second[0];
                bool contained_alignment;

                if (strlen(name_db) > 0)
                    contained_alignment = ProcessAlignment(ovl, reads[ovl->read_A_id_],
                                                           reads[ovl->read_B_id_], ALN_THRESHOLD, THETA, THETA2, true);
                else
                    contained_alignment = ProcessAlignment(ovl, reads[ovl->read_A_id_],
                                                           reads[ovl->read_B_id_], ALN_THRESHOLD, THETA, THETA2, false);
                if (contained_alignment == true) {
                    containing_read = ovl->read_B_id_;
                }

                if (reads[ovl->read_B_id_]->active == true)
                    contained = contained or contained_alignment;

                //Filter matches that matter.
                //TODO Figure out a way to do this more efficiently
                if ((ovl->match_type_ == FORWARD) or (ovl->match_type_ == FORWARD_INTERNAL))
                    matches_forward[i].push_back(it->second[0]);
                else if ((ovl->match_type_ == BACKWARD) or (ovl->match_type_ == BACKWARD_INTERNAL))
                    matches_backward[i].push_back(it->second[0]);

            }

        }
        if (contained) {
            reads[i]->active = false;
            contained_out << i << "\t" << containing_read << std::endl;

        }
    }


    for (int i = 0; i < n_read; i++) {//Isn't this just 0 or 1?
        num_overlaps += matches_forward[i].size() + matches_backward[i].size();
        for (int j = 0; j < matches_forward[i].size(); j++)
            rev_complemented_matches += matches_forward[i][j]->reverse_complement_match_;
        for (int j = 0; j < matches_backward[i].size(); j++)
            rev_complemented_matches += matches_backward[i][j]->reverse_complement_match_;
    }
    console->info("{} overlaps", num_overlaps);
    console->info("{} rev overlaps", rev_complemented_matches);

    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            num_active_read++;
            maximal_reads << i << std::endl;
        }
    }
    console->info("removed contained reads, active reads: {}", num_active_read);

    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) num_active_read++;
    }
    console->info("active reads: {}", num_active_read);

    num_overlaps = 0;
    num_forward_overlaps = 0;
    num_forward_internal_overlaps = 0;
    num_reverse_overlaps = 0;
    num_reverse_internal_overlaps = 0;
    rev_complemented_matches = 0;
    int rev_complemented_fwd_matches(0), rev_complemented_bck_matches(0), rev_complemented_fwd_int_matches(0),
            rev_complemented_bck_int_matches(0);

    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (reads[matches_forward[i][j]->read_B_id_]->active) {
                    num_overlaps++;
                    if (matches_forward[i][j]->match_type_ == FORWARD) {
                        num_forward_overlaps++;
                        rev_complemented_fwd_matches += matches_forward[i][j]->reverse_complement_match_;
                    }
                    else if (matches_forward[i][j]->match_type_ == FORWARD_INTERNAL) {
                        num_forward_internal_overlaps++;
                        rev_complemented_fwd_int_matches += matches_forward[i][j]->reverse_complement_match_;
                    }
                    if (matches_forward[i][j]->reverse_complement_match_ == 1)
                        rev_complemented_matches++;
                }
            }
            //std::cout <<"First for done "<<std::endl;
            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (reads[matches_backward[i][j]->read_B_id_]->active) {
                    num_overlaps++;
                    if (matches_backward[i][j]->match_type_ == BACKWARD) {
                        num_reverse_overlaps++;
                        rev_complemented_bck_matches += matches_backward[i][j]->reverse_complement_match_;
                    }
                    else if (matches_backward[i][j]->match_type_ == BACKWARD_INTERNAL) {
                        num_reverse_internal_overlaps++;
                        rev_complemented_bck_int_matches += matches_backward[i][j]->reverse_complement_match_;
                    }
                    if (matches_backward[i][j]->reverse_complement_match_ == 1)
                        rev_complemented_matches++;
                }
            }
        }
    }
    /*std::cout<<num_overlaps << " overlaps " << num_forward_overlaps << " fwd overlaps "
    << num_forward_internal_overlaps << " fwd internal overlaps "<< num_reverse_overlaps
    << " backward overlaps " << num_reverse_internal_overlaps
    << " backward internal overlaps "<< rev_complemented_matches << " reverse complement overlaps\n"
    << rev_complemented_fwd_matches <<" rev cmplment fwd matches "
    << rev_complemented_fwd_int_matches << " rev cmplement fwd int matches "
    << rev_complemented_bck_matches << " rev cmplment bck matches "
    << rev_complemented_bck_int_matches << " rev cmplement bck int matches " << std::endl;*/

# pragma omp parallel for
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            std::sort(matches_forward[i].begin(), matches_forward[i].end(), compare_overlap_weight);
            std::sort(matches_backward[i].begin(), matches_backward[i].end(), compare_overlap_weight);
        }
    }

    // temporary
    FILE *G_out;
    G_out = fopen("edges.g_out.txt", "w");
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (reads[matches_forward[i][j]->read_B_id_]->active) {
                    fprintf(G_out, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n",
                            matches_forward[i][j]->read_A_id_, matches_forward[i][j]->read_B_id_,
                            matches_forward[i][j]->weight, matches_forward[i][j]->reverse_complement_match_,
                            matches_forward[i][j]->match_type_, matches_forward[i][j]->eff_read_A_match_start_,
                            matches_forward[i][j]->eff_read_A_match_end_,
                            matches_forward[i][j]->eff_read_B_match_start_,
                            matches_forward[i][j]->eff_read_B_match_end_,
                            matches_forward[i][j]->eff_read_A_start_, matches_forward[i][j]->eff_read_A_end_,
                            matches_forward[i][j]->eff_read_B_start_, matches_forward[i][j]->eff_read_B_end_);
                    break;
                }
            }
        }
    }

    fprintf(G_out, "bkw\n");

    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {
            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (reads[matches_backward[i][j]->read_B_id_]->active) {
                    fprintf(G_out, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n",
                            matches_backward[i][j]->read_A_id_, matches_backward[i][j]->read_B_id_,
                            matches_backward[i][j]->weight, matches_backward[i][j]->reverse_complement_match_,
                            matches_backward[i][j]->match_type_, matches_backward[i][j]->eff_read_A_match_start_,
                            matches_backward[i][j]->eff_read_A_match_end_,
                            matches_backward[i][j]->eff_read_B_match_start_,
                            matches_backward[i][j]->eff_read_B_match_end_,
                            matches_backward[i][j]->eff_read_A_start_, matches_backward[i][j]->eff_read_A_end_,
                            matches_backward[i][j]->eff_read_B_start_, matches_backward[i][j]->eff_read_B_end_);
                    break;
                }
            }
        }
    }

    FILE *out_backup;
    out_backup = fopen("edges.fwd.backup.txt", "w");
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (reads[matches_forward[i][j]->read_B_id_]->active)
                    fprintf(out_backup, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n",
                            matches_forward[i][j]->read_A_id_, matches_forward[i][j]->read_B_id_,
                            matches_forward[i][j]->weight, matches_forward[i][j]->reverse_complement_match_,
                            matches_forward[i][j]->match_type_, matches_forward[i][j]->eff_read_A_match_start_,
                            matches_forward[i][j]->eff_read_A_match_end_,
                            matches_forward[i][j]->eff_read_B_match_start_,
                            matches_forward[i][j]->eff_read_B_match_end_,
                            matches_forward[i][j]->eff_read_A_start_, matches_forward[i][j]->eff_read_A_end_,
                            matches_forward[i][j]->eff_read_B_start_, matches_forward[i][j]->eff_read_B_end_);
            }
    }
    fclose(out_backup);
    out_backup = fopen("edges.bkw.backup.txt", "w");
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active)
            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (reads[matches_backward[i][j]->read_B_id_]->active)
                    fprintf(out_backup, "%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n",
                            matches_backward[i][j]->read_A_id_, matches_backward[i][j]->read_B_id_,
                            matches_backward[i][j]->weight, matches_backward[i][j]->reverse_complement_match_,
                            matches_backward[i][j]->match_type_, matches_backward[i][j]->eff_read_A_match_start_,
                            matches_backward[i][j]->eff_read_A_match_end_,
                            matches_backward[i][j]->eff_read_B_match_start_,
                            matches_backward[i][j]->eff_read_B_match_end_,
                            matches_backward[i][j]->eff_read_A_start_, matches_backward[i][j]->eff_read_A_end_,
                            matches_backward[i][j]->eff_read_B_start_, matches_backward[i][j]->eff_read_B_end_);
            }
    }
    fclose(out_backup);


    FILE *out_g1;
    FILE *out_g2;
    FILE *out_hg;
    FILE *out_hg2;
    FILE *out_greedy;
    FILE *out_skipped;

    out_g1 = fopen((std::string(out_name) + ".edges.1").c_str(), "w");
    out_g2 = fopen((std::string(out_name) + ".edges.2").c_str(), "w");

    // Output files for edges
    out_hg = fopen((std::string(out_name) + ".edges.hinges").c_str(), "w");
    out_hg2 = fopen((std::string(out_name) + ".edges.hinges2").c_str(), "w");
    out_greedy = fopen((std::string(out_name) + ".edges.greedy").c_str(), "w");
    out_skipped = fopen((std::string(out_name) + ".edges.skipped").c_str(), "w");

    // All hinges ikmported from the hinges.txt file
    std::unordered_map<int, std::vector<Hinge> > hinges_vec;

    // Hinges that we were previously killed in filter.cpp due to bridging
    std::unordered_map<int, std::vector<Hinge> > killed_hinges_vec;

    // Hinges that will be killed for being matched with a hinge in killed_hinges_vec
    std::unordered_map<int, std::vector<Hinge> > new_killed_hinges_vec;

    int n = 0;
    int kh = 0;
    for (int i = 0; i < n_read; i++) {
        hinges_vec[i] = std::vector<Hinge>();
        std::set<std::pair<int, int> > surviving_hinges(marked_hinges[i].begin(), marked_hinges[i].end());
        for (int j = 0; j < marked_hinges[i].size(); j++) {
            hinges_vec[i].push_back(Hinge(marked_hinges[i][j].first, marked_hinges[i][j].second, true));
            if (reads[i]->active) {
                n++;
            }
        }
        for (int j = 0; j < marked_repeats[i].size(); j++) {
            if (surviving_hinges.find(marked_repeats[i][j]) == surviving_hinges.end()) {
                killed_hinges_vec[i].push_back(Hinge(marked_repeats[i][j].first, marked_repeats[i][j].second, false));
                if (reads[i]->active) {
                    kh++;
                }
            }
        }
    }
    console->info("{} killed hinges", kh);
    console->info("{} hinges", n);

    std::ofstream killed_out(out + ".killed.hinges");
    for (int i = 0; i < n_read; i++) {
        killed_out << i << " ";
        for (int j = 0; j < killed_hinges_vec[i].size(); j++) {
            killed_out << killed_hinges_vec[i][j].type << " " << killed_hinges_vec[i][j].pos << " ";
        }
        killed_out << std::endl;
    }

    n = 0;
    for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < hinges_vec[i].size(); j++) {
            if ((reads[i]->active) and (hinges_vec[i][j].active)) n++;
        }
    }
    console->info("{} active hinges", n);

    /**
     * Switch to naive hinge filtering
     * Keep the hinge only if there are HINGE_READS reads that start near the hinge and continue to the end of the read
     */

    /*int HINGE_READS = 1;

    for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < hinges_vec[i].size(); j++) {
            int num_near_hinge_reads = 0;
            if ((reads[i]->active) and (hinges_vec[i][j].active) and (hinges_vec[i][j].type == 1)) {
                // count reads that start near the hinge and continue to the end of the read
                printf("read %d hinge %d type %d pos %d ", i, j, 1, hinges_vec[i][j].pos);
                num_near_hinge_reads = 0;
                for (int k = 0; k < matches_forward[i].size(); k ++ ) {
                    if ((matches_forward[i][k]->match_type_ == FORWARD) and
                            (reads[matches_forward[i][k]->read_B_id_]->active)
                            and abs((matches_forward[i][k]->eff_read_A_match_start_ -  hinges_vec[i][j].pos ) < 300))
                        num_near_hinge_reads ++;
                }
                printf("num %d\n", num_near_hinge_reads);
            } else if ((reads[i]->active) and (hinges_vec[i][j].active) and (hinges_vec[i][j].type == -1)) {
                printf("read %d hinge %d type %d pos %d ", i, j, -1, hinges_vec[i][j].pos);
                num_near_hinge_reads = 0;
                for (int k = 0; k < matches_backward[i].size(); k ++ ) {
                    if ((matches_backward[i][k]->match_type_ == BACKWARD) and
                            (reads[matches_backward[i][k]->read_B_id_]->active)
                            and (abs(matches_backward[i][k]->eff_read_A_match_end_ -  hinges_vec[i][j].pos ) < 300))
                        num_near_hinge_reads ++;
                }
                printf("num %d\n", num_near_hinge_reads);
            }
            //if (num_near_hinge_reads != HINGE_READS) hinges_vec[i][j].active = false;
        }
    }*/




    // TODO: Technically we dont need this filtering, as we can use the hinge graph
    // construction to do the filtering as well



    for (int i = 0; i < n_read; i++) {
        //This is in essence the filtering step
        //For each read find the best forward match, and remove all incoming hinges starting after the start
        //of the match corresponding to this.
        //Update 2/19: Now, we remove any in-hinge (out-hinge) if there is a FORWARD or FORWARD_INTERNAL match
        // (BACKWARD or BACKWARD_INTERNAL) that starts on or before (after) the hinge. 40 is error margin.

        if (reads[i]->active) {
            int forward = 0;
            int backward = 0;
            for (int j = 0; j < matches_forward[i].size(); j++) {
                if (matches_forward[i][j]->active) {
                    if (((matches_forward[i][j]->match_type_ == FORWARD) or
                         (matches_forward[i][j]->match_type_ == FORWARD_INTERNAL)) and
                        (reads[matches_forward[i][j]->read_B_id_]->active)) {

                        for (int k = 0; k < hinges_vec[i].size(); k++) {
                            if ((((matches_forward[i][j]->eff_read_A_match_start_ <
                                   hinges_vec[i][k].pos + KILL_HINGE_INTERNAL_ALLOWANCE) and
                                  (matches_forward[i][j]->match_type_ == FORWARD_INTERNAL))
                                 or ((matches_forward[i][j]->eff_read_A_match_start_ <
                                      hinges_vec[i][k].pos - KILL_HINGE_OVERLAP_ALLOWANCE) and
                                     (matches_forward[i][j]->match_type_ == FORWARD)))
                                and (hinges_vec[i][k].type == 1)) {
                                hinges_vec[i][k].active = false;

                            }
                        }
                        //}
                        //forward++;
                    }
                }
            }

            for (int j = 0; j < matches_backward[i].size(); j++) {
                if (matches_backward[i][j]->active) {
                    if (((matches_backward[i][j]->match_type_ == BACKWARD) or
                         (matches_backward[i][j]->match_type_ == BACKWARD_INTERNAL)) and
                        (reads[matches_backward[i][j]->read_B_id_]->active)) {
                        // if (backward < 1) {
                        //remove certain hinges
                        for (int k = 0; k < hinges_vec[i].size(); k++) {
                            if ((((matches_backward[i][j]->eff_read_A_match_end_ >
                                   hinges_vec[i][k].pos - KILL_HINGE_INTERNAL_ALLOWANCE) and
                                  (matches_backward[i][j]->match_type_ == BACKWARD_INTERNAL)) or
                                 ((matches_backward[i][j]->eff_read_A_match_end_ >
                                   hinges_vec[i][k].pos + KILL_HINGE_OVERLAP_ALLOWANCE) and
                                  (matches_backward[i][j]->match_type_ == BACKWARD)))
                                and (hinges_vec[i][k].type == -1)) {
                                hinges_vec[i][k].active = false;

                            }
                        }
                        //}
                        //backward++;
                    }
                }
            }
        }
    }



    console->info("Building hinge graph");

    //ogdf::Graph hinge_graph;
    //ogdf::HashArray<int, ogdf::node> hinge_graph_node_list;

    int num_hinges(0);
    for (int i = 0; i < n_read; i++) {
        num_hinges+=hinges_vec[i].size();
    }

    //ogdf::Graph hinge_graph;
    //ogdf::HashArray<int, ogdf::node> hinge_graph_node_list;
    console->info("num hinges {}", num_hinges);
    Graph hinge_graph (num_hinges);
    int hg(0);
    std::map< std::pair <int, int>, int> hinge_graph_node_map;
    std::map<int, std:: pair<int,int> > hinge_graph_node_map_rev;

    for (int i=0; i< hinges_vec.size(); i++){
        for(int j=0; j < hinges_vec[i].size(); j++){
            hinge_graph_node_map[std::make_pair(i,j)]=hg;
            hinge_graph_node_map_rev[hg]= std::make_pair(i,j);
            hg++;
        }
    }

    // Hinge graph construction
    // En passant, we identify the new_killed_hinges

    FILE *out_hgraph;
    out_hgraph = fopen((std::string(out_name) + ".hgraph").c_str(), "w");

    FILE *out_debug;
    out_debug = fopen((std::string(out_name) + ".debug").c_str(), "w");

    FILE * OverlapDebugFile;
    OverlapDebugFile = fopen("overlap_debug.txt", "w");

    int pos_B;

    for (int i = 0; i < n_read; i++) {

        if (reads[i]->active) {

            for (int k = 0; k < hinges_vec[i].size(); k++) {

                for (int j = 0; j < matches_forward[i].size(); j++) {
                    if (matches_forward[i][j]->active) {
                        if (((matches_forward[i][j]->match_type_ == FORWARD) or
                             (matches_forward[i][j]->match_type_ == FORWARD_INTERNAL)) and
                            (reads[matches_forward[i][j]->read_B_id_]->active)) {


                            // Here we check whether read B has a hinge matching hinges_vec[i][k]

                            // Should we also check whether hinges are active?

                            pos_B = matches_forward[i][j]->GetMatchingPosition(hinges_vec[i][k].pos);

//                            console->info("Matching position is {}", pos_B); // for debugging
                            int req_hinge_type;

                            int rev_int = 0;

                            if (matches_forward[i][j]->reverse_complement_match_ == true) {
                                req_hinge_type = -1 * hinges_vec[i][k].type;
                                rev_int = 1;
                            }
                            else {
                                req_hinge_type = hinges_vec[i][k].type;
                            }
//                            std::cout << req_hinge_type << std::endl;


                            int b_id = matches_forward[i][j]->read_B_id_;


                            for (int l = 0; l < hinges_vec[b_id].size(); l++) {

                                if ((hinges_vec[b_id][l].pos < pos_B + MATCHING_HINGE_SLACK) and
                                    (hinges_vec[b_id][l].pos > pos_B - MATCHING_HINGE_SLACK)) {

                                    // found a matching hinge

                                    if (req_hinge_type == hinges_vec[b_id][l].type) {


                                        std::pair <int,int> first_coord, second_coord;

                                        first_coord=std::make_pair(i,k);
                                        second_coord=std::make_pair(b_id,l);


                                        if (hinges_vec[i][k].type == 1) {

                                            add_edge(hinge_graph_node_map[first_coord], hinge_graph_node_map[second_coord], hinge_graph);
                                            fprintf(out_hgraph, "%d %d %d %d %d %d\n",
                                                    i,
                                                    b_id,
                                                    hinges_vec[i][k].pos,
                                                    hinges_vec[b_id][l].pos, 1,
                                                    rev_int);

                                        }
                                        else {

                                            add_edge(hinge_graph_node_map[second_coord], hinge_graph_node_map[first_coord], hinge_graph);
                                            fprintf(out_hgraph, "%d %d %d %d %d %d\n",
                                                    b_id,
                                                    i,
                                                    hinges_vec[b_id][l].pos,
                                                    hinges_vec[i][k].pos, 1,
                                                    rev_int);
                                        }
                                    }
                                }

                            }

                            for (int l = 0; l < killed_hinges_vec[b_id].size(); l++) {
//                                std::cout << i <<"\t" << b_id <<"\t" << k << "\t" << l <<std::endl;

                                if ((killed_hinges_vec[b_id][l].pos < pos_B + MATCHING_HINGE_SLACK) and
                                    (killed_hinges_vec[b_id][l].pos > pos_B - MATCHING_HINGE_SLACK)) {

                                    // found a matching hinge
                                    if (req_hinge_type == killed_hinges_vec[b_id][l].type) {

                                        if (hinges_vec[i][k].type == 1) {

                                            fprintf(out_hgraph, "%d %d %d %d %d %d\n",
                                                    i,
                                                    b_id,
                                                    hinges_vec[i][k].pos,
                                                    killed_hinges_vec[b_id][l].pos, 0,
                                                    rev_int);
                                        }
                                        else {

                                            fprintf(out_hgraph, "%d %d %d %d %d %d\n",
                                                    b_id,
                                                    i,
                                                    killed_hinges_vec[b_id][l].pos,
                                                    hinges_vec[i][k].pos, 0,
                                                    rev_int);
                                        }

                                        if (matches_forward[i][j]->match_type_ == FORWARD) {

                                            new_killed_hinges_vec[i].push_back(Hinge(hinges_vec[i][k].pos,hinges_vec[i][k].type,false));

                                            if (hinges_vec[i][k].type == -1) {
                                                console->info("This should not have happened.");
                                                // If this is a -1 hinge, read i should also bridge the repeat,
                                                // and hinges_vec[i][k] would have been killed in filter

                                                fprintf(out_debug,"%d %d %d %d %d [%d %d] [%d %d] [%d %d] [%d %d] \n",
                                                        matches_forward[i][j]->read_A_id_, matches_forward[i][j]->read_B_id_,
                                                        matches_forward[i][j]->weight, matches_forward[i][j]->reverse_complement_match_,
                                                        matches_forward[i][j]->match_type_, matches_forward[i][j]->eff_read_A_match_start_,
                                                        matches_forward[i][j]->eff_read_A_match_end_,
                                                        matches_forward[i][j]->eff_read_B_match_start_,
                                                        matches_forward[i][j]->eff_read_B_match_end_,
                                                        matches_forward[i][j]->eff_read_A_start_, matches_forward[i][j]->eff_read_A_end_,
                                                        matches_forward[i][j]->eff_read_B_start_, matches_forward[i][j]->eff_read_B_end_);

                                                fprintf(out_debug, "%d %d %d %d\n", hinges_vec[i][k].pos,
                                                        hinges_vec[i][k].type,
                                                        killed_hinges_vec[b_id][l].pos,
                                                        killed_hinges_vec[b_id][l].type);

                                            }

                                        }


                                    }


                                }


                            }


                        }

                    }
                }


                for (int j = 0; j < matches_backward[i].size(); j++) {
                    if (matches_backward[i][j]->active) {
                        if (((matches_backward[i][j]->match_type_ == BACKWARD) or
                             (matches_backward[i][j]->match_type_ == BACKWARD_INTERNAL)) and
                            (reads[matches_backward[i][j]->read_B_id_]->active)) {

                            // Need to check whether read B has a hinge matching hinges_vec[i][k]


                            pos_B = matches_backward[i][j]->GetMatchingPosition(hinges_vec[i][k].pos);

//                            console->info("Matching position is {}", pos_B); // for debugging

                            int req_hinge_type;

                            int rev_int = 0;

                            if (matches_backward[i][j]->reverse_complement_match_ == true) {
                                req_hinge_type = -1 * hinges_vec[i][k].type;
                                rev_int = 1;
                            }
                            else {
                                req_hinge_type = hinges_vec[i][k].type;
                            }
//                            std::cout << req_hinge_type << std::endl;

                            int b_id = matches_backward[i][j]->read_B_id_;
                            for (int l = 0; l < hinges_vec[b_id].size(); l++) {

                                if ((hinges_vec[b_id][l].pos < pos_B + MATCHING_HINGE_SLACK) and
                                    (hinges_vec[b_id][l].pos > pos_B - MATCHING_HINGE_SLACK)) {


                                    // found a matching hinge
                                    std::pair <int,int> first_coord, second_coord;
                                    first_coord=std::make_pair(i,k);
                                    second_coord=std::make_pair(b_id,l);


                                    if (req_hinge_type == hinges_vec[b_id][l].type) {


                                        if (hinges_vec[i][k].type == -1) {

                                            add_edge(hinge_graph_node_map[first_coord], hinge_graph_node_map[second_coord], hinge_graph);
                                            fprintf(out_hgraph, "%d %d %d %d %d %d\n",
                                                    i,
                                                    b_id,
                                                    hinges_vec[i][k].pos,
                                                    hinges_vec[b_id][l].pos, 1,
                                                    rev_int);

                                        }
                                        else {

                                            add_edge(hinge_graph_node_map[second_coord], hinge_graph_node_map[first_coord], hinge_graph);
                                            fprintf(out_hgraph, "%d %d %d %d %d %d\n",
                                                    b_id,
                                                    i,
                                                    hinges_vec[b_id][l].pos,
                                                    hinges_vec[i][k].pos, 1,
                                                    rev_int);
                                        }


                                    }
                                }

                            }
                            for (int l = 0; l < killed_hinges_vec[b_id].size(); l++) {

                                if ((killed_hinges_vec[b_id][l].pos < pos_B + MATCHING_HINGE_SLACK) and
                                    (killed_hinges_vec[b_id][l].pos > pos_B - MATCHING_HINGE_SLACK)) {

                                    // found a matching hinge

                                    if (req_hinge_type == killed_hinges_vec[b_id][l].type) {

                                        if (hinges_vec[i][k].type == -1) {

                                            fprintf(out_hgraph, "%d %d %d %d %d %d\n",
                                                    i,
                                                    b_id,
                                                    hinges_vec[i][k].pos,
                                                    killed_hinges_vec[b_id][l].pos, 0,
                                                    rev_int);

                                        }
                                        else {

                                            fprintf(out_hgraph, "%d %d %d %d %d %d\n",
                                                    b_id,
                                                    i,
                                                    killed_hinges_vec[b_id][l].pos,
                                                    hinges_vec[i][k].pos, 0,
                                                    rev_int);
                                        }

                                    }

                                    if (matches_backward[i][j]->match_type_ == BACKWARD) {

                                        new_killed_hinges_vec[i].push_back(Hinge(hinges_vec[i][k].pos,hinges_vec[i][k].type,false));

                                        if (hinges_vec[i][k].type != -1) {
                                            console->info("This should not have happened 2.");
                                            // If this is a +1 hinge, read i should also bridge the repeat,
                                            // and hinges_vec[i][k] would have been killed in filter
                                        }

                                    }


                                }

                            }


                        }
                    }
                }

            }
        }
    }


    console->info("Hinge graph built");
    std::vector<int> component(num_vertices(hinge_graph));
    int num = connected_components(hinge_graph, &component[0]);

    std::vector<int>::size_type i;
    std::cout << "Total number of components: " << num << std::endl;

    std::map<int,int> component_size;
    for (i = 0; i != component.size(); ++i){   // are we skipping i=0?
        if ( component_size.find(component[i]) == component_size.end() ){
            component_size[component[i]]=1;
        }
        else
            component_size[component[i]]+=1;
    }
//    std::unordered_map<int, std::vector<Hinge> > filtered_hinges_vec;
//    for (int i = 0; i < n_read; i++) {
//        filtered_hinges_vec[i] = std::vector<Hinge>();
//    }




    for (int i = 0; i != component.size(); ++i) {
        if (component_size[component[i]] < MIN_CONNECTED_COMPONENT_SIZE) {
            int ind1, ind2;
            ind1 = hinge_graph_node_map_rev[i].first;
            ind2 = hinge_graph_node_map_rev[i].second;
            hinges_vec[ind1][ind2].active=false;
//            filtered_hinges_vec[ind1].push_back(hinges_vec[ind1][ind2]);
        }

    }


//    std::map< std::pair <int, int>, int> hinge_graph_node_map;


    std::map<int, std::pair <int, int>> component_sink;
    for (i = 0; i != component.size(); ++i){
        int ind1, ind2;
        ind1 = hinge_graph_node_map_rev[i].first;
        ind2 = hinge_graph_node_map_rev[i].second;

        // for now let us just pick an arbitrary active hinge as the component main sink
        if ( hinges_vec[ind1][ind2].active == true )
            component_sink[component[i]]= std::make_pair(ind1,ind2);

    }

    n = 0;
    FILE *out_hglist;
    out_hglist = fopen((std::string(out_name) + ".hinge.list").c_str(), "w");
    for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < hinges_vec[i].size(); j++) {
            if ((reads[i]->active) and ((hinges_vec[i][j].active))) {
                fprintf(out_hglist, "%d %d %d\n", i, marked_hinges[i][j].first, marked_hinges[i][j].second);
                n++;
            }
        }
    }
    fclose(out_hglist);
    console->info("after filter {} active hinges", n);


    // filter hinges
    std::vector<bool> repeat_status_front;
    std::vector<bool> repeat_status_back;

    for (int i = 0; i < n_read; i++) {
        bool in = false;
        bool out = false;
        for (int j = 0; j < hinges_vec[i].size(); j++) {
            if ((hinges_vec[i][j].active)  and (hinges_vec[i][j].type == 1)) in = true;
            if ((hinges_vec[i][j].active)  and (hinges_vec[i][j].type == -1)) out = true;
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

                                PrintOverlapToFile(out_greedy, matches_forward[i][j]);

                                if (matches_forward[i][j]->reverse_complement_match_ == 0)
                                    fprintf(out_g1, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
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
                                    fprintf(out_g1, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
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
                                    fprintf(out_g2, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
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
                                    fprintf(out_g2, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
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

                                PrintOverlapToFile(out_greedy, matches_backward[i][j]);

                                if (matches_backward[i][j]->reverse_complement_match_ == 0)
                                    fprintf(out_g1, "%d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
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
                                    fprintf(out_g1, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
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
                                    fprintf(out_g2, "%d' %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
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
                                    fprintf(out_g2, "%d %d' %d [%d %d] [%d %d] [%d %d] [%d %d]\n",
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

    /*std::cout<<num_overlaps << " overlaps " << num_forward_overlaps << " fwd overlaps "
    << num_forward_internal_overlaps << " fwd internal overlaps "<< num_reverse_overlaps
    << " backward overlaps " << num_reverse_internal_overlaps
    << " backward internal overlaps "<< rev_complemented_matches << " reverse complement overlaps" << std::endl;
     */
/*    int test_read_no(1291);
    for (int j = 0; j < matches_forward[test_read_no].size(); j++){
        if (matches_forward[test_read_no][j]->active){
            std::cout << matches_forward[test_read_no][j]->read_A_id_ << "\t"
                    << matches_forward[test_read_no][j]->read_B_id_ << "\t"
                    << matches_forward[test_read_no][j]->match_type_ << "\t"
                    << matches_forward[test_read_no][j]->weight << std::endl ;
        }
    }*/

    std::ofstream debug_fle("hinge_debug.txt");

    console->info("Starting to build assembly graph.");


//    int eff_b_id;
    int hinge_pos = -1;

    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) {

            int forward = 0;
            int forward_internal = 0;
            int backward = 0;
            int backward_internal = 0;


            LOverlap * chosen_match = NULL;


            for (int j = 0; j < matches_forward[i].size(); j++){
//                if(i==227622) {
//                    debug_fle << i << "\t" << forward << "\t" << matches_forward[i][j]->match_type_
//                            << "\t" << matches_forward[i][j]->weight << "\t"
//                    << matches_forward[i][j]->active << std::endl;
//                }

                if (matches_forward[i][j]->active) {


                    if ((reads[matches_forward[i][j]->read_B_id_]->active)) { // and (forward == 0)) {
                        //printf("hinge size %d\n", hinges_vec[matches_forward[i][j]->read_B_id_].size());


                        if ((matches_forward[i][j]->match_type_ == FORWARD) and (forward == 0)) {

                            // check if read j has new_killed_hinge
                            //TODO: should this be checked for FORWARD_INTERNAL as well?

                            bool poisoned = false;

                            for (int k = 0; k < new_killed_hinges_vec[i].size(); k++) {

                                if ( (matches_forward[i][j]->reverse_complement_match_ != 1) and
                                        (new_killed_hinges_vec[i][k].type == -1) and
                                        (new_killed_hinges_vec[i][k].pos > matches_forward[i][j]->eff_read_B_match_end_) ) {

                                    //TODO: do we need a tolerance in the comparison above?

                                    PrintOverlapToFile(out_skipped, matches_forward[i][j]);
                                    poisoned = true;

                                }
                                else if ( (matches_forward[i][j]->reverse_complement_match_ == 1) and
                                          (new_killed_hinges_vec[i][k].type == 1) and
                                          (new_killed_hinges_vec[i][k].pos < matches_forward[i][j]->eff_read_B_match_start_) ) {


                                    PrintOverlapToFile(out_skipped, matches_forward[i][j]);
                                    poisoned = true;

                                }

                            }

                            if (not poisoned) {
                                chosen_match = matches_forward[i][j];
                                hinge_pos = -1;
                                forward = 1;
                                //break;
                            }

                        }
                        else if ((matches_forward[i][j]->match_type_ == FORWARD_INTERNAL)
                                //and isValidHinge(matches_forward[i][j], hinges_vec[matches_forward[i][j]->read_B_id_])
                                  and (hinges_vec[matches_forward[i][j]->read_B_id_].size() > 0)
                                    and (forward_internal == 0)){

                            // In the case of a forward_internal match we check whether
                            // the hinge on read B is an in-hinge
                            // (or an out-hinge if it's a reverse complement match)

//                            int hinge_index = 0;

                            int read_B_match_start = matches_forward[i][j]->read_B_match_start_;
                            if (matches_forward[i][j]->reverse_complement_match_ == 1) {
                                read_B_match_start = matches_forward[i][j]->read_B_match_end_;
                            }

                            for (int k = 0; k < hinges_vec[matches_forward[i][j]->read_B_id_].size(); k++) {
                                if ( (read_B_match_start >
                                        hinges_vec[matches_forward[i][j]->read_B_id_][k].pos - HINGE_TOLERANCE)
                                            and (read_B_match_start <
                                        hinges_vec[matches_forward[i][j]->read_B_id_][k].pos + HINGE_TOLERANCE)
                                            and (hinges_vec[matches_forward[i][j]->read_B_id_][k].type ==
                                        (1-2*matches_forward[i][j]->reverse_complement_match_))
                                            and (hinges_vec[matches_forward[i][j]->read_B_id_][k].active) ) {

                                    if ((forward == 0) or
                                            (matches_forward[i][j]->weight > chosen_match->weight - 2*HINGE_SLACK)) {

                                        chosen_match = matches_forward[i][j];
                                        forward = 1;
                                        forward_internal = 1;

                                        hinge_pos = hinges_vec[matches_forward[i][j]->read_B_id_][k].pos;




                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            if (chosen_match != NULL) {
                PrintOverlapToFile(out_hg,chosen_match);

                edges_forward[i].push_back(chosen_match);

                PrintOverlapToFile2(out_hg2,chosen_match,hinge_pos);

                chosen_match = NULL;
            }
            else {

                // Deadend debugging

                // Forward dead-end
                deadend_out << i;
//                deadend_out << "\t Active: " << reads[i]->active << std::endl;
                deadend_out << "\t matches_forward size: " << matches_forward[i].size() << std::endl;

            }

            for (int j = 0; j < matches_backward[i].size(); j++){


                if (matches_backward[i][j]->active) {

                    if ((reads[matches_backward[i][j]->read_B_id_]->active)) {


                        if ((matches_backward[i][j]->match_type_ == BACKWARD) and (backward == 0)){

                            // check if read j has new_killed_hinge

                            bool poisoned = false;

                            for (int k = 0; k < new_killed_hinges_vec[i].size(); k++) {

                                if ( (matches_backward[i][j]->reverse_complement_match_ != 1) and
                                     (new_killed_hinges_vec[i][k].type == 1) and
                                     (new_killed_hinges_vec[i][k].pos < matches_backward[i][j]->eff_read_B_match_start_) ) {

                                    //TODO: do we need a tolerance in the comparison above?

                                    PrintOverlapToFile(out_skipped, matches_backward[i][j]);
                                    poisoned = true;

                                }
                                else if ( (matches_backward[i][j]->reverse_complement_match_ == 1) and
                                          (new_killed_hinges_vec[i][k].type == -1) and
                                          (new_killed_hinges_vec[i][k].pos > matches_backward[i][j]->eff_read_B_match_end_) ) {

                                    PrintOverlapToFile(out_skipped, matches_backward[i][j]);
                                    poisoned = true;

                                }

                            }

                            if (not poisoned) {
                                chosen_match = matches_backward[i][j];
                                backward = 1;
                                hinge_pos = -1;
                            }

                        }
                        else if ((matches_backward[i][j]->match_type_ == BACKWARD_INTERNAL)
                                and (hinges_vec[matches_backward[i][j]->read_B_id_].size() > 0)
                                    and (backward_internal == 0)) {

                            // In the case of a backward_internal match
                            // we check whether the hinge on read B is an in-hinge
                            // (or an in-hinge if it's a reverse complement match)


                            int read_B_match_end = matches_backward[i][j]->read_B_match_end_;
                            if (matches_backward[i][j]->reverse_complement_match_ == 1) {
                                read_B_match_end = matches_backward[i][j]->read_B_match_start_;
                            }

                            for (int k = 0; k < hinges_vec[matches_backward[i][j]->read_B_id_].size(); k++) {


                                if ( (read_B_match_end >
                                        hinges_vec[matches_backward[i][j]->read_B_id_][k].pos - HINGE_TOLERANCE)
                                     and (read_B_match_end <
                                        hinges_vec[matches_backward[i][j]->read_B_id_][k].pos + HINGE_TOLERANCE)
                                     and (hinges_vec[matches_backward[i][j]->read_B_id_][k].type ==
                                        (-1+2*matches_backward[i][j]->reverse_complement_match_))
                                     and (hinges_vec[matches_backward[i][j]->read_B_id_][k].active) ) {

                                    if ((backward == 0) or
                                            (matches_backward[i][j]->weight > chosen_match->weight - 2*HINGE_SLACK)) {
                                        chosen_match = matches_backward[i][j];
                                        backward = 1;
                                        backward_internal = 1;


                                        hinge_pos = hinges_vec[matches_backward[i][j]->read_B_id_][k].pos;

//                                        int hinge_graph_id = hinge_graph_node_map[std::make_pair(matches_backward[i][j]->read_B_id_,k)];

//                                        eff_b_id = component_sink[component[hinge_graph_id]].first;

                                    }

                                    break;
                                }
                            }
                        }
                    }
                }
            }

            if (chosen_match != NULL) {
                PrintOverlapToFile(out_hg,chosen_match);

                edges_backward[i].push_back(chosen_match);

                PrintOverlapToFile2(out_hg2,chosen_match,hinge_pos);

            }
            else {
                // Deadend debugging

                // Backward dead-end
                deadend_out << i;
//                deadend_out << "\t Active: " << reads[i]->active << std::endl;
                deadend_out << "\t matches_backward size: " << matches_backward[i].size() << std::endl;

            }
        }
    }

    console->info("sort and output finished");
    console->info("version 0.0.3");

    if (strlen(name_db) > 0)
    la.closeDB(); //close database
    return 0;
}
