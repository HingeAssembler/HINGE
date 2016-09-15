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

#include "INIReader.h"
#include "spdlog/spdlog.h"
#include "DB.h"
#include "align.h"
#include "LAInterface.h"
#include "cmdline.h"




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
    cmdp.add("debug", '\0', "debug mode");
    cmdp.parse_check(argc, argv);

    LAInterface la;
    const char * name_db = cmdp.get<std::string>("db").c_str(); //.db file of reads to load
    const char * name_las = cmdp.get<std::string>("las").c_str();//.las file of alignments
    const char * name_paf = cmdp.get<std::string>("paf").c_str();
    const char * name_fasta = cmdp.get<std::string>("fasta").c_str();
    const char * name_config = cmdp.get<std::string>("config").c_str();//name of the configuration file, in INI format
    std::string out = cmdp.get<std::string>("prefix");
    bool has_qv = true;
    const char * name_restrict = cmdp.get<std::string>("restrictreads").c_str();
    /**
     * There are two sets of input, the first is db+las, which corresponds to daligner as an overlapper,
     * the other is fasta + paf, which corresponds to minimap as an overlapper.
     */

    namespace spd = spdlog;

    //auto console = spd::stdout_logger_mt("console",true);

    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
    sinks.push_back(std::make_shared<spdlog::sinks::daily_file_sink_st>(cmdp.get<std::string>("log") + "/log", "txt", 23, 59));
    auto console = std::make_shared<spdlog::logger>("log", begin(sinks), end(sinks));
    spdlog::register_logger(console);
    //auto console = std::make_shared<spdlog::logger>("name", begin(sinks), end(sinks));


    console->info("Reads filtering");

    if (cmdp.exist("debug")) {
        char *buff = (char *) malloc(sizeof(char) * 2000);
        getwd(buff);
        console->info("current user {}, current working directory {}", getlogin(), buff);
        free(buff);
    }

    console->info("name of db: {}, name of .las file {}", name_db, name_las);
    console->info("name of fasta: {}, name of .paf file {}", name_fasta, name_paf);


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
        n_read = la.loadFASTA(name_fasta,reads);
        has_qv = false;
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


    std::vector<std::vector<int>>  QV;

    if (strlen(name_db) > 0) {
        la.getRead(reads,0,n_read);
        la.getQV(QV,0,n_read); // load QV track from .db file
    }

    console->info("Input data finished");

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
    int N_PROC = reader.GetInteger("running", "n_proc", 4);
    int EST_COV = reader.GetInteger("filter", "ec", 0); // load the estimated coverage (probably from other programs) from ini file, if it is zero, then estimate it
    int reso = 40; // resolution of masks, repeat annotation, coverage, etc  = 40 basepairs
    bool use_qv_mask = reader.GetBoolean("filter", "use_qv", true);
    bool use_coverage_mask = reader.GetBoolean("filter", "coverage", true);
    int COVERAGE_FRACTION = (int) reader.GetInteger("filter", "coverage_frac_repeat_annotation", 3);
    const int MIN_REPEAT_ANNOTATION_THRESHOLD = (int) reader.GetInteger("filter", "min_repeat_annotation_threshold", 10);
    const int MAX_REPEAT_ANNOTATION_THRESHOLD = (int) reader.GetInteger("filter", "max_repeat_annotation_threshold", 20);
    const int REPEAT_ANNOTATION_GAP_THRESHOLD = (int) reader.GetInteger("filter", "repeat_annotation_gap_threshold",300);
    const int NO_HINGE_REGION = (int) reader.GetInteger("filter", "no_hinge_region",500);
    //How far two hinges of the same type can be
    const int HINGE_MIN_SUPPORT = (int) reader.GetInteger("filter", "hinge_min_support", 7);
    //Minimum number of reads that have to start in a reso length interval to be considered in hinge calling
    const int HINGE_BIN_PILEUP_THRESHOLD = (int) reader.GetInteger("filter", "hinge_min_pileup", 7);
    //Minimum number of reads to have in a pileup to consider a hinge bridged
    const int HINGE_READ_UNBRIDGED_THRESHOLD = (int) reader.GetInteger("filter", "hinge_unbridged", 6);
    //Number of reads that one has to see before a pileup to declare a potential hinge unbridged
    int HINGE_BIN_LENGTH = (int) reader.GetInteger("filter", "hinge_bin", 100);
    //Physical length of the bins considered
    const int HINGE_TOLERANCE_LENGTH = (int) reader.GetInteger("filter", "hinge_tolerance_length", 100);
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


    std::vector<std::vector <LOverlap * > > idx_pileup; // this is the pileup
    std::vector<std::vector <LOverlap * > > idx_pileup_dedup; // this is the deduplicated pileup

    std::vector<std::unordered_map<int, std::vector<LOverlap *> > > idx_ab; //unordered_map from (aid, bid) to alignments in a vector



    for (int i = 0; i< n_read; i++) {
        idx_pileup.push_back(std::vector<LOverlap *>());
        idx_pileup_dedup.push_back(std::vector<LOverlap *>());
        idx_ab.push_back(std::unordered_map<int, std::vector<LOverlap *>> ());
    }

    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx_pileup[aln[i]->read_A_id_].push_back(aln[i]);
        }
    }

# pragma omp parallel for
    for (int i = 0; i < n_read; i++) {// sort overlaps of a reads
        std::sort(idx_pileup[i].begin(), idx_pileup[i].end(), compare_overlap);
    }

# pragma omp parallel for
    for (int i = 0; i < aln.size(); i++) {
        idx_ab[aln[i]->read_A_id_][aln[i]->read_B_id_] = std::vector<LOverlap *>();
    }

# pragma omp parallel for
    for (int i = 0; i < aln.size(); i++) {
        idx_ab[aln[i]->read_A_id_][aln[i]->read_B_id_].push_back(aln[i]);
    }


# pragma omp parallel for
    for (int i = 0; i < n_read; i++) {
        for (std::unordered_map<int, std::vector<LOverlap *> >::iterator it = idx_ab[i].begin(); it!= idx_ab[i].end(); it++) {
            std::sort(it->second.begin(), it->second.end(), compare_overlap);
            if (it->second.size() > 0)
                idx_pileup_dedup[i].push_back(it->second[0]);
        }
    }

    console->info("profile coverage (with and without CUT_OFF)");


    std::ofstream cov(out + ".coverage.txt");
    std::ofstream homo(out + ".homologous.txt");
    std::ofstream rep(out + ".repeat.txt");
    std::ofstream filtered(out + ".filtered.fasta");
    std::ofstream hg(out + ".hinges.txt");
    std::ofstream mask(out + ".mas");

    std::vector< std::vector<std::pair<int, int> > > coverages(n_read);
    std::vector< std::vector<std::pair<int, int> > > cutoff_coverages(n_read);
    std::vector< std::vector<std::pair<int, int> > > cgs(n_read); //coverage gradient;
    //std::vector< std::vector<std::pair<int, int> > > his;
     for (int i = 0; i < n_read; i ++) {
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

    console->info("profile coverage done");

    std::set<int> rand_reads;
    srand(time(NULL));
    rand_reads.insert(0);
    while (rand_reads.size() < n_read/500){
        int rd_id=rand()%n_read;
        if (reads[rd_id]->len > 5000)
            rand_reads.insert(rd_id);
    }

    int num_slot = 0;
    long int total_cov = 0;

    std::vector<int> read_coverage;
    long int read_cov=0;
    int read_slot =0;
    //Finding the average coverage, probing a small proportion of reads

//    for (std::set<int>::iterator it=rand_reads.begin();it!=rand_reads.end(); ++it) {
    for (int i =0; i <n_read;  i++){
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

    std::vector<std::pair<int, int>> maskvec;
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

    for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < cutoff_coverages[i].size(); j++) {
            cutoff_coverages[i][j].second -= MIN_COV;
            if (cutoff_coverages[i][j].second < 0) cutoff_coverages[i][j].second = 0;
        }

        //get the longest consecutive region that has decent coverage, decent coverage = estimated coverage / 3
        int start = 0;
        int end = start;
        int maxlen = 0, maxstart = 0, maxend = 0;
        for (int j = 0; j < cutoff_coverages[i].size(); j++) {
            if (cutoff_coverages[i][j].second > 0) {
                end = cutoff_coverages[i][j].first;
            } else {
                if (end > start) {
                    //std::cout<<"read" << i << " "<<start+reso << "->" << end << std::endl;
                    if (end - start - reso > maxlen) {
                        maxlen = end - start - reso;
                        maxstart = start + reso;
                        maxend = end;
                    }
                }
                start = cutoff_coverages[i][j].first;
                end = start;
            }
        }
        //std::cout << i << " " << maxstart << " " << maxend << std::endl;
        //int s = std::max(maxstart, QV_mask[i].first);
        //int l = std::min(maxend, QV_mask[i].second) - std::max(maxstart, QV_mask[i].first);
        //if (l < 0) l = 0;
        //filtered << ">read_" << i << std::endl;
        //filtered << reads[i]->bases.substr(s,l) << std::endl;

        if (reads_to_keep.size()>0) {
            if (reads_to_keep.find(i) == reads_to_keep.end()){
//                std::cout<<"setting masks equal";
                maxend=maxstart;
                QV_mask[i].second=QV_mask[i].first;
            }
        }
        if ((use_qv_mask) and (use_coverage_mask)) {
            maskvec.push_back(
                    std::pair<int, int>(std::max(maxstart, QV_mask[i].first), std::min(maxend, QV_mask[i].second)));
            //get the interestion of two masks
            mask << i << " " << std::max(maxstart, QV_mask[i].first) << " " << std::min(maxend, QV_mask[i].second) << std::endl;
        } else if ((use_coverage_mask) and (not use_qv_mask)) {
            maskvec.push_back(std::pair<int, int>(maxstart, maxend));
            mask << i << " " << maxstart << " " << maxend << std::endl;
        } else {
            maskvec.push_back(std::pair<int, int>(QV_mask[i].first, QV_mask[i].second));
            mask << i << " " << QV_mask[i].first << " " << QV_mask[i].second << std::endl;
        }
    }

    /*FILE* temp_out1;
    FILE* temp_out2;
    temp_out1=fopen("coverage.debug.txt","w");
    temp_out2=fopen("coverage_gradient.debug.txt","w");

    for (int i=0; i< n_read ; i++) {
        fprintf(temp_out1,"%d \t", i);
        for (int j=0; j < coverages[i].size(); j++){
            fprintf(temp_out1,"%d:%d \t", coverages[i][j].first,coverages[i][j].second);
        }
        fprintf(temp_out1,"\n");
    }

    for (int i=0; i< n_read ; i++) {
        fprintf(temp_out2,"%d \t", i);
        for (int j=0; j < cgs[i].size(); j++){
            fprintf(temp_out2,"%d:%d \t", cgs[i][j].first,cgs[i][j].second);
        }
        fprintf(temp_out2,"\n");
    }
    fclose(temp_out1);
    fclose(temp_out2);*/

    /*for (int i = 0; i < maskvec.size(); i++) {
        printf("read %d %d %d\n", i, maskvec[i].first, maskvec[i].second);
        printf("QV: read %d %d %d\n", i, QV_mask[i].first, QV_mask[i].second);
    }*/


    //binarize coverage gradient;


    std::vector<std::vector<std::pair<int, int> > > repeat_annotation;
    //detect repeats based on coverage gradient, mark it has rising (1) or falling (-1)
    for (int i = 0; i < n_read; i++) {
        std::vector<std::pair<int, int> > anno;
        for (int j = 0; j < cgs[i].size()-1; j++) { // changed, remove the last one
            //std::cout<< i << " " << cgs[i][j].first << " " << cgs[i][j].second << std::endl;

            if ((cgs[i][j].first >= maskvec[i].first + NO_HINGE_REGION) and (cgs[i][j].first <= maskvec[i].second - NO_HINGE_REGION)) {
                if (cgs[i][j].second > std::min(
                        std::max((coverages[i][j].second+MIN_COV)/COVERAGE_FRACTION, MIN_REPEAT_ANNOTATION_THRESHOLD),
                        MAX_REPEAT_ANNOTATION_THRESHOLD))
                    anno.push_back(std::pair<int, int>(cgs[i][j].first, 1));
                else if (cgs[i][j].second < - std::min(
                        std::max((coverages[i][j].second+MIN_COV)/COVERAGE_FRACTION, MIN_REPEAT_ANNOTATION_THRESHOLD),
                        MAX_REPEAT_ANNOTATION_THRESHOLD))
                    anno.push_back(std::pair<int, int>(cgs[i][j].first, -1));
            }
        }
        repeat_annotation.push_back(anno);
    }


    // clean it a bit, merge consecutive 1, or consecutive -1, or adjacent 1 and -1 if their position is within gap_threshold (could be bursty error)
    for (int i = 0; i < n_read; i++) {
        for (std::vector<std::pair<int, int> >::iterator iter = repeat_annotation[i].begin(); iter < repeat_annotation[i].end(); ) {
            if (iter+1 < repeat_annotation[i].end()){
                if (((iter->second == 1) and ((iter + 1)->second == 1)) and
                        ((iter+1)->first - iter->first < REPEAT_ANNOTATION_GAP_THRESHOLD)) {
                    repeat_annotation[i].erase((iter + 1));
                } else if (((iter->second == -1) and ((iter + 1)->second == -1)) and
                        ((iter+1)->first - iter->first < REPEAT_ANNOTATION_GAP_THRESHOLD)) {
                    iter = repeat_annotation[i].erase(iter);
                } else iter++;
            } else iter ++;
        }
    }

    //remove gaps
//    for (int i = 0; i < n_read; i++) {
//        for (std::vector<std::pair<int, int> >::iterator iter = repeat_annotation[i].begin(); iter < repeat_annotation[i].end(); ) {
//            if (iter+1 < repeat_annotation[i].end()){
//                if ((iter->second == -1) and ((iter+1)->second == 1) and
//                        ((iter+1)->first - iter->first < REPEAT_ANNOTATION_GAP_THRESHOLD)){
//                    iter = repeat_annotation[i].erase(iter);
//                    iter = repeat_annotation[i].erase(iter); // fill gaps
//                } else if ((iter->second == 1) and ((iter+1)->second == -1) and
//                        ((iter+1)->first - iter->first < REPEAT_ANNOTATION_GAP_THRESHOLD)) {
//                    iter = repeat_annotation[i].erase(iter);
//                    iter = repeat_annotation[i].erase(iter);
//                } else iter++;
//            } else iter ++;
//        }
//    }


    /*temp_out1=fopen("repeat_annotation.debug.txt","w");
    for (int i = 0; i < n_read; i++) {
        fprintf(temp_out1,"%d \t%d\t",i,repeat_annotation[i].size());
        for (std::vector<std::pair<int, int> >::iterator iter = repeat_annotation[i].begin(); iter < repeat_annotation[i].end();iter++) {
            fprintf(temp_out1,"%d:%d\t",iter->first,iter->second);
        }
        fprintf(temp_out1,"\n");
    }
    fclose(temp_out1);*/
    // need a better hinge detection

    // get hinges from repeat annotation information
    std::unordered_map<int, std::vector<std::pair<int, int>> > hinges;
    // n_read pos -1 = in_hinge 1 = out_hinge
    std::ofstream  debug_file("debug.txt");
    for (int i = 0; i < n_read; i++) {
        //std::cout << i <<std::endl;
        hinges[i] = std::vector<std::pair<int, int>>();

        int coverage_at_start(0);
        int num_at_start(0);
        int num_at_end(0);
        int coverage_at_end(0);
        float avg_coverage_at_start;
        float avg_coverage_at_end;
        for (int j = 0; j < coverages[i].size(); j++){
            if ((coverages[i][j].first <= maskvec[i].first + NO_HINGE_REGION) and
                (coverages[i][j].first >= maskvec[i].first )){
                coverage_at_start += coverages[i][j].second;
                num_at_start++;
            }
            if ((coverages[i][j].first <= maskvec[i].second ) and
                (coverages[i][j].first >= maskvec[i].second - NO_HINGE_REGION )){
                coverage_at_end += coverages[i][j].second;
                num_at_end++;
            }
        }

        avg_coverage_at_end = (float)coverage_at_end/num_at_end;
        avg_coverage_at_start = (float)coverage_at_start/num_at_start;
        if (std::abs(avg_coverage_at_end-avg_coverage_at_start) < 10){
            continue;
        }

        for (int j = 0; j < repeat_annotation[i].size(); j++) {

            if (repeat_annotation[i][j].second == -1) { // look for out hinges, negative gradient

                bool bridged = true;
                int support = 0;
                int num_reads_at_end=1;

                std::vector<std::pair<int,int> > read_other_ends;


                for (int k = 0; k < idx_pileup[i].size(); k++) {

                    int left_overhang, right_overhang;
                    int temp_id;
                    temp_id=idx_pileup[i][k]->read_B_id_;

                    if (idx_pileup[i][k]->reverse_complement_match_==0){
                        right_overhang= std::max(maskvec[temp_id].second-idx_pileup[i][k]->read_B_match_end_,0);
                        left_overhang= std::max(idx_pileup[i][k]->read_B_match_start_- maskvec[temp_id].first,0);
                    }
                    else if (idx_pileup[i][k]->reverse_complement_match_==1) {
                        right_overhang= std::max(idx_pileup[i][k]->read_B_match_start_- maskvec[temp_id].first,0);
                        left_overhang= std::max(maskvec[temp_id].second-idx_pileup[i][k]->read_B_match_end_,0);
                    }



                    if (right_overhang > THETA) {
                        if ((idx_pileup[i][k]->read_A_match_end_ >
                             repeat_annotation[i][j].first - HINGE_TOLERANCE_LENGTH)
                            and (idx_pileup[i][k]->read_A_match_end_ <
                                 repeat_annotation[i][j].first + HINGE_TOLERANCE_LENGTH)) {


                            std::pair <int,int> other_end;
                            other_end.first=idx_pileup[i][k]->read_A_match_start_;
                            other_end.second=left_overhang;
                            read_other_ends.push_back(other_end);
                            support++;
                        }
                    }
                }

                if (support < HINGE_MIN_SUPPORT){
                    continue;
                }

                std::sort(read_other_ends.begin(),read_other_ends.end(), pairAscend);

                int num_reads_considered=0;
                int num_reads_extending_to_end=0;
                int num_reads_with_internal_overlaps=0;

                for (int id = 0; id < read_other_ends.size() ; ++id) {
                    if (read_other_ends[id].first -maskvec[i].first < HINGE_BIN_LENGTH){
                        num_reads_considered++;
                        num_reads_extending_to_end++;

                        if ((num_reads_extending_to_end > HINGE_READ_UNBRIDGED_THRESHOLD) or
                            ((num_reads_considered > HINGE_READ_UNBRIDGED_THRESHOLD) and
                            (read_other_ends[id].first - read_other_ends[0].first > HINGE_BIN_LENGTH))) {
                            bridged=false;
                            break;
                        }
                    }
                    else if  (read_other_ends[id].second < THETA){
                        num_reads_considered++;

                        if ((num_reads_extending_to_end > HINGE_READ_UNBRIDGED_THRESHOLD) or
                            ((num_reads_considered > HINGE_READ_UNBRIDGED_THRESHOLD) and
                             (read_other_ends[id].first - read_other_ends[0].first > HINGE_BIN_LENGTH))) {
                            bridged=false;
                            break;
                        }
                    }
                    else if (read_other_ends[id].second > THETA) {
                        num_reads_with_internal_overlaps++;
                        num_reads_considered++;
                        int id1=id+1;
                        int pileup_length=1;

                        while (id1 < read_other_ends.size()){
                            if (read_other_ends[id1].first - read_other_ends[id].first  < HINGE_BIN_LENGTH){
                                pileup_length++;
                                id1++;
                            }
                            else{
                                break;
                            }
                        }

                        if (pileup_length > HINGE_BIN_PILEUP_THRESHOLD){
                            bridged=true;
                            break;
                        }
                    }
                }
                if ((not bridged) and (support > HINGE_MIN_SUPPORT))
                    hinges[i].push_back(std::pair<int, int>(repeat_annotation[i][j].first,-1));


            } else { // look for in_hinges, positive gradient
                bool bridged = true;
                int support = 0;
                int num_reads_at_end=1;

                std::vector<std::pair<int,int> > read_other_ends;


                for (int k = 0; k < idx_pileup[i].size(); k++) {
                    int left_overhang, right_overhang;
                    int temp_id;
                    temp_id=idx_pileup[i][k]->read_B_id_;

                    if (idx_pileup[i][k]->reverse_complement_match_==0){
                        right_overhang= std::max(maskvec[temp_id].second-idx_pileup[i][k]->read_B_match_end_,0);
                        left_overhang= std::max(idx_pileup[i][k]->read_B_match_start_- maskvec[temp_id].first,0);
                    }
                    else if (idx_pileup[i][k]->reverse_complement_match_==1) {
                        right_overhang= std::max(idx_pileup[i][k]->read_B_match_start_- maskvec[temp_id].first,0);
                        left_overhang= std::max(maskvec[temp_id].second-idx_pileup[i][k]->read_B_match_end_,0);
                    }


                    if (left_overhang > THETA) {
                        if ((idx_pileup[i][k]->read_A_match_start_ >
                             repeat_annotation[i][j].first - HINGE_TOLERANCE_LENGTH)
                            and (idx_pileup[i][k]->read_A_match_start_ <
                                 repeat_annotation[i][j].first + HINGE_TOLERANCE_LENGTH)) {

                            std::pair <int,int> other_end;
                            other_end.first=idx_pileup[i][k]->read_A_match_end_;
                            other_end.second=right_overhang;
                            read_other_ends.push_back(other_end);
                            support++;
                        }
                    }
                }
                if (support < HINGE_MIN_SUPPORT){
                    continue;
                }


                std::sort(read_other_ends.begin(),read_other_ends.end(),pairDescend);//Sort in descending order


                int num_reads_considered=0;
                int num_reads_extending_to_end=0;
                int num_reads_with_internal_overlaps=0;



                for (int id = 0; id < read_other_ends.size() ; ++id) {
                    if (maskvec[i].second-read_other_ends[id].first < HINGE_BIN_LENGTH){
                        num_reads_considered++;
                        num_reads_extending_to_end++;

                        if ((num_reads_extending_to_end > HINGE_READ_UNBRIDGED_THRESHOLD) or
                            ((num_reads_considered > HINGE_READ_UNBRIDGED_THRESHOLD) and
                             (read_other_ends[0].first - read_other_ends[id].first > HINGE_BIN_LENGTH))) {
                            bridged=false;
                            break;
                        }
                    }
                    else if  (read_other_ends[id].second < THETA){
                        num_reads_considered++;

                        if ((num_reads_extending_to_end > HINGE_READ_UNBRIDGED_THRESHOLD) or
                            ((num_reads_considered > HINGE_READ_UNBRIDGED_THRESHOLD) and
                             (read_other_ends[0].first - read_other_ends[id].first > HINGE_BIN_LENGTH))) {
                            bridged=false;
                            break;
                        }
                    }
                    else if (read_other_ends[id].second > THETA) {
                        num_reads_with_internal_overlaps++;
                        num_reads_considered++;
                        int id1=id+1;
                        int pileup_length=1;

                        while (id1 < read_other_ends.size()){
                            if (read_other_ends[id].first - read_other_ends[id1].first  < HINGE_BIN_LENGTH){
                                pileup_length++;
                                id1++;
                            }
                            else{
                                break;
                            }
                        }

                        if (pileup_length > HINGE_BIN_PILEUP_THRESHOLD){
                            bridged=true;
                            break;
                        }
                    }
                }

                if ((not bridged) and (support > HINGE_MIN_SUPPORT))
                    hinges[i].push_back(std::pair<int, int>(repeat_annotation[i][j].first, 1));


            }
        }
    }

    //output hinges

    int ra_cnt = 0;

    for (int i = 0; i < n_read; i++) {
        rep << i << " ";
        for (int j = 0; j < repeat_annotation[i].size(); j++) {
            rep << repeat_annotation[i][j].first << " " << repeat_annotation[i][j].second << " ";
        }
        ra_cnt += repeat_annotation[i].size();
        rep << std::endl;
    }
    rep.close();
    console->info("Number of hinges before filtering: {}", ra_cnt);

    int hg_cnt = 0;

    for (int i = 0; i < n_read; i++) {
        hg << i << " ";
        for (int j = 0; j < hinges[i].size(); j++) {
            hg << hinges[i][j].first << " " << hinges[i][j].second << " ";
        }
        hg_cnt += hinges[i].size();
        hg << std::endl;
    }
    hg.close();

    console->info("Number of hinges: {}", hg_cnt);

    if (strlen(name_db)>0)
        la.closeDB(); //close database
    return 0;
}
