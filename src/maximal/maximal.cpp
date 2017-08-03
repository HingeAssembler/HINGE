#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <tuple>
#include <random>
#include <omp.h>
#include <time.h>
#include <glob.h>


#include "INIReader.h"
#include "spdlog/spdlog.h"
#include "DB.h"
#include "align.h"
#include "LAInterface.h"
#include "cmdline.h"

std::string lastN(std::string input, int n)
{
    return input.substr(input.size() - n);
}

inline std::vector<std::string> glob(const std::string& pat){
    using namespace std;
    glob_t glob_result;
    int i = 1;
    std::string search_name;
    search_name = pat + "."+std::to_string(i)+".las";
    std::cout << search_name << endl;
    glob(search_name.c_str(),GLOB_TILDE,NULL,&glob_result);
//    std::cout << "Number of files " << glob_result.gl_pathc << std::endl;

    vector<string> ret;


    while (glob_result.gl_pathc != 0){
        ret.push_back(string(glob_result.gl_pathv[0]));
        i ++;
        search_name = pat + "."+std::to_string(i)+".las";
        glob(search_name.c_str(),GLOB_TILDE,NULL,&glob_result);
//        std::cout << "Number of files " << glob_result.gl_pathc << std::endl;
    }

    std::cout << "-------------------------"<< std::endl;
    std::cout << "Number of files " << i-1 << std::endl;
    std::cout << "Input string " << pat.c_str() << std::endl;
    std::cout << "-------------------------"<< std::endl;

    globfree(&glob_result);
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
    match->eff_read_A_read_start_ = read_A->effective_start;
    match->eff_read_A_read_end_ = read_A->effective_end;

    // removed the following if, so that things agree with the convention for reverse complement matches

    match->eff_read_B_read_start_ = read_B->effective_start;
    match->eff_read_B_read_end_ = read_B->effective_end;

//    if (match->reverse_complement_match_ == 0) {
//        match->eff_read_B_read_start_ = read_B->effective_start;
//        match->eff_read_B_read_end_ = read_B->effective_end;
//    } else {
//        match->eff_read_B_read_start_ = read_B->len - read_B->effective_end;
//        match->eff_read_B_read_end_ = read_B->len - read_B->effective_start;
//    }

    /*printf("bef %d %d %d [%d %d] [%d %d] [%d %d] [%d %d]\n", match->read_A_id_, match->read_B_id_,
     * match->reverse_complement_match_,
        match->read_A_match_start_, match->read_A_match_end_, match->read_B_match_start_, match->read_B_match_end_,
           match->eff_read_A_read_start_, match->eff_read_A_read_end_, match->eff_read_B_read_start_, match->eff_read_B_read_end_
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
           match->eff_read_A_read_start_, match->eff_read_A_read_end_, match->eff_read_B_read_start_, match->eff_read_B_read_end_
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

    match->length = match->read_A_match_end_ - match->read_A_match_start_
                    + match->read_B_match_end_ - match->read_B_match_start_;

    return contained;
}

std::vector<std::pair<int,int>> Merge(std::vector<LOverlap *> & intervals, int cutoff)
//Returns sections of read a which are covered by overlaps. Each overlap is considered as
// <start_pos+cutoff,end_pos-cutoff>.
{
    //std::cout<<"Merge"<<std::endl;
    std::vector<std::pair<int, int > > ret;
    int n = intervals.size(); // Length of the vector intervals
    if (n == 0) return ret;

    if(n == 1) {
        ret.push_back(std::pair<int,int>(intervals[0]->read_A_match_start_, intervals[0]->read_A_match_end_));
        return ret;
    }

    //Where is sort defined ? Is this std::sort?
    sort(intervals.begin(),intervals.end(),compare_overlap_abpos); //sort according to left (start position of
    // overlap beginning on a)

    int left= intervals[0]->read_A_match_start_ + cutoff, right = intervals[0]->read_A_match_end_ - cutoff;
    //left, right means maximal possible interval now

    for(int i = 1; i < n; i++) {
        //Ovl1 ~ Ovl2 if Ovl1 and Ovl2 have a nonzero intersection. (that is both the b read maps
        // to the same position on the a read)
        //This defines a chain of  connected overlaps. This for loop returns a a vector ret which
        // is a pair of <start of connected overlaps, end of connected overlaps>
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

//Interval = pair<int, int>. Defined in LAInterface.h
Interval Effective_length(std::vector<LOverlap *> & intervals, int min_cov) {
//Returns <start_pos, end_pos>
//start_pos : the first position at which Read a of the overlaps have at least min_cov matches on it.
//end_pos : the last position that the  (#overlaps- min_cov)th read (in order of start positions ends).
//Should compare_overlap_aepos actually compare read_A_match_end_? If that is done, then the end_pos
// will be the last position
// on the a read so that all positions beyond have less than min_cov matches on them
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

bool bridge(LOverlap* ovl, int s, int e){
    //Returns True if [s e] on read a is bridged by ovl. False else.
    //Put 500 in a typedef perhaps?
    return ((ovl->read_A_match_start_ < s - 500) and (ovl->read_A_match_end_ > e + 500));
}

float number_of_bridging_reads(std::vector<LOverlap *> ovl_reads, int hinge_location, int hinge_type,int threshold){
    int num_bridging_reads=0;
    //int threshold=100;
    std::vector<int> read_ends;
    if (hinge_type==1){
        for (int i=0; i < ovl_reads.size(); i++){
            if ((ovl_reads[i]->read_A_match_start_ > hinge_location-threshold ) and
                (ovl_reads[i]->read_A_match_start_ < hinge_location+threshold ))
                read_ends.push_back(ovl_reads[i]->read_A_match_end_);
        }
    }
    else if (hinge_type==-1){
        for (int i=0; i < ovl_reads.size(); i++){
            if ((ovl_reads[i]->read_A_match_end_ > hinge_location-threshold ) and
                (ovl_reads[i]->read_A_match_end_ < hinge_location+threshold ))
                read_ends.push_back(ovl_reads[i]->read_A_match_start_);
        }
    }
    std::sort(read_ends.begin(),read_ends.end(), std::greater<int>());
    int start_point=0;
    int num_bins=0;
    for (int i=0; i<read_ends.size(); i++) {
        std::cout << hinge_location <<"\t"<< read_ends[i]<< std::endl;
        if (read_ends[start_point] - read_ends[i] > 2 * threshold) {
            num_bins++;
            start_point = i;
        }
    }
    return num_bins/((float)1);
}


int main(int argc, char *argv[]) {

    mkdir("log",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    cmdline::parser cmdp;
    cmdp.add<std::string>("db", 'b', "db file name", false, "");
    cmdp.add<std::string>("las", 'l', "las file name", false, "");
    cmdp.add<std::string>("paf", 'p', "paf file name", false, "");
    cmdp.add<std::string>("config", 'c', "configuration file name", false, "");
    cmdp.add<std::string>("fasta", 'f', "fasta file name", false, "");
    cmdp.add<std::string>("prefix", 'x', "prefix of (intermediate) output", false, "out");
    cmdp.add<std::string>("restrictreads",'r',"restrict to reads in the file",false,"");
    cmdp.add<std::string>("log", 'g', "log folder name", false, "log");
    cmdp.add("mlas", '\0', "multiple las files");
    cmdp.add("debug", '\0', "debug mode");
    cmdp.parse_check(argc, argv);

    LAInterface la;
    const char * name_db = cmdp.get<std::string>("db").c_str(); //.db file of reads to load
    const char * name_las_base = cmdp.get<std::string>("las").c_str();//.las file of alignments
    const char * name_paf = cmdp.get<std::string>("paf").c_str();
    const char * name_fasta = cmdp.get<std::string>("fasta").c_str();
    const char * name_config = cmdp.get<std::string>("config").c_str();//name of the configuration file, in INI format
    std::string out = cmdp.get<std::string>("prefix");
    bool has_qv = true;
    const char * name_restrict = cmdp.get<std::string>("restrictreads").c_str();
    std::string name_mask = out + ".mas";

    namespace spd = spdlog;

    //auto console = spd::stdout_logger_mt("console",true);

    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
    sinks.push_back(std::make_shared<spdlog::sinks::daily_file_sink_st>(cmdp.get<std::string>("log") + "/log", "txt", 23, 59));
    auto console = std::make_shared<spdlog::logger>("log", begin(sinks), end(sinks));
    spdlog::register_logger(console);
    //auto console = std::make_shared<spdlog::logger>("name", begin(sinks), end(sinks));


    console->info("Getting maximal reads");

    bool db_and_las, db_or_las, fa_and_paf, fa_or_paf;
    db_and_las = (strlen(name_db) > 0) and (strlen(name_las_base) > 0);
    db_or_las = (strlen(name_db) > 0) or (strlen(name_las_base) > 0);
    fa_and_paf = (strlen(name_fasta) > 0) and (strlen(name_paf) > 0);
    fa_or_paf = (strlen(name_fasta) > 0) or (strlen(name_paf) > 0);

    if (db_or_las and fa_or_paf){
        console->error("Pass in either a db and a las or a fasta and a paf");
        return 1;
    }

    if (( not fa_and_paf) and (not db_and_las)){
        console->error("Pass in at least one of the following two combinations: a db and a las or a fasta and a paf");
        return 1;
    }

    std::string name_las_string;
    if (cmdp.exist("mlas")) {
        if (not db_and_las){
            console->error("--mlas works only with db and las");
            return 1;
        }
        name_las_string = std::string(name_las_base);
    }
    else if (strlen(name_las_base) > 0) {
        if (lastN(std::string(name_las_base), 4) == ".las")
            name_las_string = std::string(name_las_base);
        else
            name_las_string = std::string(name_las_base) + ".las";
    }


    const char * name_las = name_las_string.c_str();
    /**
     * There are two sets of input, the first is db+las, which corresponds to daligner as an overlapper,
     * the other is fasta + paf, which corresponds to minimap as an overlapper.
     */




    console->info("name of db: {}, name of .las file {}", name_db, name_las);
    console->info("name of fasta: {}, name of .paf file {}", name_fasta, name_paf);


    std::ifstream ini_file(name_config);
    std::string str((std::istreambuf_iterator<char>(ini_file)),
                    std::istreambuf_iterator<char>());

    console->info("Parameters passed in \n{}", str);

    if (strlen(name_db) > 0)
        la.openDB(name_db);


    std::vector<std::string> name_las_list;
    std::string name_las_str(name_las);
    console->info("Las files: {}", name_las_str);
    if (cmdp.exist("mlas")) {
        console->info("Calling glob.");
        name_las_list = glob(name_las_str);
    }
    else
        name_las_list.push_back(name_las_str);




    int n_read;
    if (strlen(name_db) > 0)
        n_read = la.getReadNumber();

    std::vector<Read *> reads; //Vector of pointers to all reads

    if (strlen(name_fasta) > 0) {
        n_read = la.loadFASTA(name_fasta,reads);
        has_qv = false;
    }


    console->info("# Reads: {}", n_read); // output some statistics




    std::vector<std::vector<int>>  QV;

    if (strlen(name_db) > 0) {
        la.getRead(reads,0,n_read);
        if (la.getQV(QV,0,n_read) != 0) // load QV track from .db file
            has_qv = false;
    }



    if (has_qv)
        for (int i = 0; i < n_read; i++) {
            for (int j = 0; j < QV[i].size(); j++) QV[i][j] = int(QV[i][j] < 40);
        }
    //Binarize QV vector, 40 is the threshold
    std::set<int> reads_to_keep, reads_to_keep_initial;
    char * line = NULL;
    size_t len = 0;
    if (strlen(name_restrict) > 0){
        FILE * restrict_reads;
        restrict_reads = fopen(name_restrict, "r");
        while (getline(&line, &len, restrict_reads) != -1){
            std::stringstream ss;
            ss.clear();
            ss << line;
            int num;
            ss >> num;
            reads_to_keep.insert(num);
        }
        fclose(restrict_reads);
        console->info("Reads to debug loaded from: {}", name_restrict);
        console->info("Number of reads to debug loaded: {}", reads_to_keep.size());
    }
    else
        console->info("No debug restrictions.");



    if (strlen(name_las_list[0].c_str()) > 0)
        la.openAlignmentFile(name_las_list[0]); // get tspace

    std::vector<std::pair<int, int> > QV_mask(n_read);
    // QV_mask is the mask based on QV for reads, for each read, it has one pair [start, end]

    if (has_qv) {
        for (int i = 0; i < n_read; i++) {
            int s = 0, e = 0;
            int max = 0, maxs = s, maxe = e;

            for (int j = 0; j < QV[i].size(); j++) {
                if ((QV[i][j] == 1) and (j<QV[i].size() - 1)) {
                    e ++;
                }
                else {
                    if (e - s > max) {
                        maxe = e ; maxs = s;
                        max = e - s;
                    }

                    s = j+1;
                    e = j+1;
                }
            }
            // get the longest consecutive region that has good QV
            //printf("maxs %d maxe %d size%d\n",maxs, maxe,QV[i].size());

            QV_mask[i] = (std::pair<int, int>(maxs*la.tspace, maxe*la.tspace));
            // tspace the the interval of trace points
            // create mask by QV
        }
    }

    INIReader reader(name_config);
    if (reader.ParseError() < 0) {
        console->warn("Can't load {}", name_config);
        return 1;
    }

    int LENGTH_THRESHOLD = reader.GetInteger("filter", "length_threshold", -1);
    double QUALITY_THRESHOLD = reader.GetReal("filter", "quality_threshold", 0.0);
    int N_ITER = reader.GetInteger("filter", "n_iter", -1);
    int ALN_THRESHOLD = reader.GetInteger("filter", "aln_threshold", -1);
    int MIN_COV = reader.GetInteger("filter", "min_cov", -1);
    int CUT_OFF = reader.GetInteger("filter", "cut_off", -1);
    int THETA = reader.GetInteger("filter", "theta", -1);
    int THETA2 = (int) reader.GetInteger("filter", "theta2", 0);
    int N_PROC = reader.GetInteger("running", "n_proc", 4);
    int EST_COV = reader.GetInteger("filter", "ec", 0); // load the estimated coverage (probably from other programs) from ini file, if it is zero, then estimate it
    int reso = 40; // resolution of masks, repeat annotation, coverage, etc  = 40 basepairs
    bool use_qv_mask = reader.GetBoolean("filter", "use_qv", true);
    bool use_coverage_mask = reader.GetBoolean("filter", "coverage", true);
    int COVERAGE_FRACTION = (int) reader.GetInteger("filter", "coverage_frac_repeat_annotation", 3);
    const int MIN_REPEAT_ANNOTATION_THRESHOLD = (int) reader.GetInteger("filter", "min_repeat_annotation_threshold", 10);
    const int MAX_REPEAT_ANNOTATION_THRESHOLD = (int) reader.GetInteger("filter", "max_repeat_annotation_threshold", 20);
    const int REPEAT_ANNOTATION_GAP_THRESHOLD = (int) reader.GetInteger("filter", "repeat_annotation_gap_threshold",300);
    //How far two hinges of the same type can be
    const int NO_HINGE_REGION = (int) reader.GetInteger("filter", "no_hinge_region",500);
    const int HINGE_MIN_SUPPORT = (int) reader.GetInteger("filter", "hinge_min_support", 7);
    //Minimum number of reads that have to start in a reso length interval to be considered in hinge calling
    const int HINGE_BIN_PILEUP_THRESHOLD = (int) reader.GetInteger("filter", "hinge_min_pileup", 7);
    //Minimum number of reads to have in a pileup to consider a hinge bridged
    const int HINGE_READ_UNBRIDGED_THRESHOLD = (int) reader.GetInteger("filter", "hinge_unbridged", 6);
    //Number of reads that one has to see before a pileup to declare a potential hinge unbridged
    int HINGE_BIN_LENGTH = (int) reader.GetInteger("filter", "hinge_bin", 100);
    //Physical length of the bins considered
    const int HINGE_TOLERANCE_LENGTH = (int) reader.GetInteger("filter", "hinge_tolerance_length", 100);
    bool USE_TWO_MATCHES = (int) reader.GetInteger("layout", "use_two_matches", 1);

    //Reads starting at +/- HINGE_TOLERANCE_LENGTH are considered reads starting at hinges
    HINGE_BIN_LENGTH=2*HINGE_TOLERANCE_LENGTH;

    console->info("use_qv_mask set to {}",use_qv_mask);
    use_qv_mask = use_qv_mask and has_qv;

    console->info("use_qv_mask set to {}",use_qv_mask);

    omp_set_num_threads(N_PROC);
    console->info("number processes set to {}", N_PROC);

    console->info("LENGTH_THRESHOLD = {}",LENGTH_THRESHOLD);
    console->info("QUALITY_THRESHOLD = {}",QUALITY_THRESHOLD);
    console->info("N_ITER = {}",N_ITER);
    console->info("ALN_THRESHOLD = {}",ALN_THRESHOLD);
    console->info("MIN_COV = {}",MIN_COV);
    console->info("CUT_OFF = {}",CUT_OFF);
    console->info("THETA = {}",THETA);
    console->info("EST_COV = {}",EST_COV);
    console->info("reso = {}",reso);
    console->info("use_coverage_mask = {}",use_coverage_mask);
    console->info("COVERAGE_FRACTION = {}",COVERAGE_FRACTION);
    console->info("MIN_REPEAT_ANNOTATION_THRESHOLD = {}",MIN_REPEAT_ANNOTATION_THRESHOLD);
    console->info("MAX_REPEAT_ANNOTATION_THRESHOLD = {}",MAX_REPEAT_ANNOTATION_THRESHOLD);
    console->info("REPEAT_ANNOTATION_GAP_THRESHOLD = {}",REPEAT_ANNOTATION_GAP_THRESHOLD);
    console->info("NO_HINGE_REGION = {}",NO_HINGE_REGION);
    console->info("HINGE_MIN_SUPPORT = {}",HINGE_MIN_SUPPORT);
    console->info("HINGE_BIN_PILEUP_THRESHOLD = {}",HINGE_BIN_PILEUP_THRESHOLD);
    console->info("HINGE_READ_UNBRIDGED_THRESHOLD = {}",HINGE_READ_UNBRIDGED_THRESHOLD);
    console->info("HINGE_BIN_LENGTH = {}",HINGE_BIN_LENGTH);
    console->info("HINGE_TOLERANCE_LENGTH = {}",HINGE_TOLERANCE_LENGTH);




    std::vector<LOverlap *> aln;//Vector of pointers to all alignments
    std::vector< std::vector<std::pair<int, int> > > coverages(n_read);
    std::vector< std::vector<std::pair<int, int> > > cutoff_coverages(n_read);
    std::vector< std::vector<std::pair<int, int> > > cgs(n_read); //coverage gradient;
    std::vector<std::pair<int, int>> maskvec;
    std::vector<std::vector<std::pair<int, int> > > repeat_annotation;
    std::unordered_map<int, std::vector<std::pair<int, int>> > hinges;


    std::ofstream cov(out + ".coverage.txt");
    std::ofstream homo(out + ".homologous.txt");
    std::ofstream filtered(out + ".filtered.fasta");
    std::ofstream contained_out(out + ".contained.txt");
    std::ofstream maximal_reads(out + ".max");


    FILE *mask_file;
    mask_file = fopen(name_mask.c_str(), "r");
    int read, rs, re;

    while (fscanf(mask_file, "%d %d %d", &read, &rs, &re) != EOF) {
        reads[read]->effective_start = rs;
        reads[read]->effective_end = re;
    }
    console->info("read mask finished");

    int num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->active) num_active_read++;
    }
    console->info("active reads at start: {}", num_active_read);


    num_active_read = 0;
    for (int i = 0; i < n_read; i++) {
        if (reads[i]->effective_end - reads[i]->effective_start < LENGTH_THRESHOLD) {
            reads[i]->active = false;
        }
        else num_active_read++;
    }
    console->info("active reads after correcting for read lengths: {}", num_active_read);

    int number_of_parts;
    if (strlen(name_las) > 0)
        number_of_parts = name_las_list.size();
    else if(strlen(name_paf) > 0)
        number_of_parts = 1;
    else {
        console->error("Need to provide either las and db or paf and fasta");
        return 1;
    }

    console->info("number of las files: {}", number_of_parts);

    for (int part = 0; part < name_las_list.size(); part++) {


        console->info("name of las: {}", name_las_list[part]);




        int64 n_aln = 0;

        if(strlen(name_las_base)> 0) {

            if (strlen(name_las_list[part].c_str()) > 0)
                la.openAlignmentFile(name_las_list[part]);
            if (strlen(name_las_list[part].c_str()) > 0) {
                n_aln = la.getAlignmentNumber();
                console->info("Load alignments from {}", name_las_list[part]);
                console->info("# Alignments: {}", n_aln);
            }
            if (strlen(name_las_list[part].c_str()) > 0) {
                la.resetAlignment();
                la.getOverlap(aln, 0, n_read);
            }
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

        console->info("Input data finished, part {}/{}", part + 1, number_of_parts);



        int r_begin = aln.front()->read_A_id_;
        int r_end = aln.back()->read_A_id_;


        std::vector<std::vector <LOverlap * > > idx_pileup; // this is the pileup
        std::vector<std::vector <LOverlap * > > idx_pileup_dedup; // this is the deduplicated pileup
        std::vector<std::unordered_map<int, std::vector<LOverlap *> > > idx_ab; //unordered_map from (aid, bid) to alignments in a vector



        for (int i = 0; i< n_read; i++) {
            idx_pileup.push_back(std::vector<LOverlap *>());
            idx_pileup_dedup.push_back(std::vector<LOverlap *>());
            idx_ab.push_back(std::unordered_map<int, std::vector<LOverlap *>> ());
            repeat_annotation.push_back(std::vector<std::pair<int, int> >());
            maskvec.push_back(std::pair<int, int>());
        }

        for (int i = 0; i < aln.size(); i++) {
            if (aln[i]->read_A_id_ == aln[i]->read_B_id_) {
                aln[i]->active = false;
            }
            if (aln[i]->active) {
                idx_pileup[aln[i]->read_A_id_].push_back(aln[i]);
            }
        }




        for (int i = 0; i < n_read; i++) {// sort overlaps of a reads
            std::sort(idx_pileup[i].begin(), idx_pileup[i].end(), compare_overlap);
        }

        for (int i = 0; i < aln.size(); i++) {
            idx_ab[aln[i]->read_A_id_][aln[i]->read_B_id_] = std::vector<LOverlap *>();
        }

        for (int i = 0; i < aln.size(); i++) {
            idx_ab[aln[i]->read_A_id_][aln[i]->read_B_id_].push_back(aln[i]);
        }


        for (int i = 0; i < n_read; i++) {
            for (std::unordered_map<int, std::vector<LOverlap *> >::iterator it = idx_ab[i].begin(); it!= idx_ab[i].end(); it++) {
                std::sort(it->second.begin(), it->second.end(), compare_overlap);
                if (it->second.size() > 0)
                    idx_pileup_dedup[i].push_back(it->second[0]);
            }
        }

        console->info("profile coverage (with and without CUT_OFF)");

        //std::vector< std::vector<std::pair<int, int> > > his;
        for (int i = r_begin; i <= r_end; i ++) {
            std::vector<std::pair<int, int> > coverage;

            std::vector<std::pair<int, int> > cutoff_coverage;


            //TODO : Implement set based gradient
            std::vector<std::pair<int, int> > cg;
            //profileCoverage: get the coverage based on pile-o-gram
            la.profileCoverage(idx_pileup[i], cutoff_coverage, reso, CUT_OFF);
            la.profileCoverage(idx_pileup[i], coverage, reso, 0);
            cov << "read " << i <<" ";
            for (int j = 0; j < coverage.size(); j++)
                cov << coverage[j].first << ","  << coverage[j].second << " ";
            cov << std::endl;

            //Computes coverage gradients.
            if (coverage.size() >= 2)
                for (int j = 0; j < coverage.size() - 1; j++) {
                    cg.push_back(std::pair<int,int>(coverage[j].first, coverage[j+1].second - coverage[j].second));
                }
            else cg.push_back(std::pair<int, int> (0,0));

            coverages[i] = (coverage);
            cutoff_coverages[i] = (cutoff_coverage);
            cgs[i] = (cg);
        }

        console->info("profile coverage done part {}/{}", part + 1, number_of_parts);


        std::set<int> rand_reads;
        srand(time(NULL));
        rand_reads.insert(0);
        int temp_index(0);
        while (rand_reads.size() < (r_end - r_begin)/500){
            temp_index ++;
            int rd_id=rand()%(r_end - r_begin) + r_begin;
            if (reads[rd_id]->len > 5000)
                rand_reads.insert(rd_id);
            if (temp_index > 20000)
                break;
        }

        int num_slot = 0;
        long int total_cov = 0;

        std::vector<int> read_coverage;
        long int read_cov=0;
        int read_slot =0;
        //Finding the average coverage, probing a small proportion of reads

//    for (std::set<int>::iterator it=rand_reads.begin();it!=rand_reads.end(); ++it) {
        for (int i =r_begin; i <= r_end;  i++){
            if (reads[i]->len < 5000)
                continue;
            read_cov=0;
            read_slot=0;
            for (int j = 0; j < coverages[i].size(); j++) {
                //printf("%d\n", coverages[i][j].second);
                read_cov+=coverages[i][j].second;
                read_slot++;
            }
            total_cov += read_cov;
            num_slot += read_slot;
            int mean_read_cov=read_cov / std::max(1,read_slot);
            read_coverage.push_back(mean_read_cov);
        }




        size_t median_id = read_coverage.size() / 2;
        if (median_id > 0)
            std::nth_element(read_coverage.begin(), read_coverage.begin()+median_id, read_coverage.end());

        int cov_est= read_coverage[median_id];

        int mean_cov_est = total_cov / num_slot;


        //get estimated coverage

        if (EST_COV != 0) cov_est = EST_COV;
        console->info("Estimated mean coverage: {}", mean_cov_est); //if the coverage is specified by ini file, cover the estimated one
        console->info("Estimated median coverage: {}", cov_est);


        // mask vector, same format as mask_QV
        if (MIN_COV < cov_est/3)
            MIN_COV = cov_est/3;

        if (reads_to_keep.size()>0) {
            reads_to_keep_initial = reads_to_keep;
            for (std::set<int>::iterator iter = reads_to_keep_initial.begin();
                 iter != reads_to_keep_initial.end(); ++iter) {
                int i = *iter;
                for (std::unordered_map<int, std::vector<LOverlap *> >::iterator it = idx_ab[i].begin();
                     it != idx_ab[i].end(); it++) {
                    if (it->second.size() > 0) {
                        LOverlap *ovl = it->second[0];
                        reads_to_keep.insert(ovl->read_B_id_);
                    }
                }
            }
            console->info("After accounting for neighbours of reads selected, have {} reads", reads_to_keep.size());
        }

        std::unordered_map<int, std::vector<LOverlap *> > matches_forward, matches_backward;

        for (int i = r_begin; i <= r_end; i ++) {
            //An initialisation for loop
            //TODO Preallocate memory. Much more efficient.
            //idx2.push_back(std::vector<LOverlap *>());
            matches_forward[i] = std::vector<LOverlap *>();
            matches_backward[i] = std::vector<LOverlap *>();
        }




        for (int i = r_begin; i <= r_end; i ++) {
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


                if ((it->second.size() > 1) and (USE_TWO_MATCHES)) {
                    //Figure out if read is contained
                    LOverlap *ovl = it->second[1];
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
                        matches_forward[i].push_back(it->second[1]);
                    else if ((ovl->match_type_ == BACKWARD) or (ovl->match_type_ == BACKWARD_INTERNAL))
                        matches_backward[i].push_back(it->second[1]);

                }




            }
            if (contained) {
                reads[i]->active = false;
                contained_out << i << "\t" << containing_read << std::endl;

            }
        }
        int num_overlaps = 0;
        int num_forward_overlaps(0), num_forward_internal_overlaps(0), num_reverse_overlaps(0),
                num_reverse_internal_overlaps(0), rev_complemented_matches(0);
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
        for (int i = r_begin; i <= r_end; i ++) {
            if (reads[i]->active) {
                num_active_read++;
                maximal_reads << i << std::endl;
            }
        }
        console->info("removed contained reads, active reads: {}", num_active_read);

        num_active_read = 0;
        for (int i = r_begin; i <= r_end; i ++) {
            if (reads[i]->active) num_active_read++;
        }
        console->info("active reads: {}", num_active_read);
        console->info("total reads: {}", r_end-r_begin+1);

        if (strlen(name_las) > 0) {
            for (int i = 0; i < aln.size(); i++) {
                delete aln[i];
            }
            aln.clear();
        }
    }



    if (strlen(name_db)>0)
        la.closeDB(); //close database
    return 0;




}
