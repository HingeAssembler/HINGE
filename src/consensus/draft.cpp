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

extern "C" {
#include "common.h"
}


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



std::vector<int> get_mapping(std::string aln_tag1, std::string aln_tag2) {
    int pos = 0;
    int count = 0;
    int count2 = 0;

    std::vector<int> ret;
    while (pos < aln_tag1.size()) {
        if (aln_tag1[pos] != '-') {
            ret.push_back(count2);
            count ++;
        }
        if (aln_tag2[pos] != '-') {
            count2 ++;
        }
        pos++;
    }
    return ret;
}



std::string reverse_complement(std::string seq) {
    static std::map<char, char> m = {{'a','t'}, {'c','g'}, {'g','c'}, {'t','a'}, {'A','T'}, {'C','G'}, {'T','A'}, {'G','C'}, {'n','n'}, {'N', 'N'}, {'-', '-'}};
    std::reverse(seq.begin(), seq.end());
    for (int i = 0; i < seq.size(); i++) {
        seq[i] = m[seq[i]];
    }
    return seq;
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



bool compare_overlap(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->read_A_match_end_ - ovl1->read_A_match_start_
             + ovl1->read_B_match_end_ - ovl1->read_B_match_start_) >
            (ovl2->read_A_match_end_ - ovl2->read_A_match_start_
             + ovl2->read_B_match_end_ - ovl2->read_B_match_start_));
}


bool compare_overlap_weight(LOverlap * ovl1, LOverlap * ovl2) {
    return (ovl1->weight > ovl2->weight);
}


bool compare_overlap_abpos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->read_A_match_start_ < ovl2->read_A_match_start_;
}

bool compare_overlap_aepos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->read_A_match_start_ > ovl2->read_A_match_start_;
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
    cmdp.add<std::string>("path", 0, "path file name", false, "path");


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
    std::string path_name = cmdp.get<std::string>("path");
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

    console->info("draft consensus");
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


    int MIN_COV2 = reader.GetInteger("draft", "min_cov", -1);
    int EDGE_TRIM = reader.GetInteger("draft", "trim", -1);
    int EDGE_SAFE = reader.GetInteger("draft", "edge_safe", -1);
    int TSPACE = reader.GetInteger("draft", "tspace", -1);
    int STEP = reader.GetInteger("draft", "step", -1);

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
    std::vector<Edge_w> edgelist, edgelist_ms; // save output to edgelist
    std::vector<std::unordered_map<int, std::vector<LOverlap *> > > idx_ab;


    for (int i = 0; i < n_read; i++) {
        //An initialisation for loop
        //TODO Preallocate memory. Much more efficient.
        idx_ab.push_back(std::unordered_map<int, std::vector<LOverlap *> >());
    }

    for (int i = 0; i < aln.size(); i++) {
        idx_ab[aln[i]->read_A_id_][aln[i]->read_B_id_] = std::vector<LOverlap *>();
    }

    for (int i = 0; i < aln.size(); i++) {
        idx_ab[aln[i]->read_A_id_][aln[i]->read_B_id_].push_back(aln[i]);
    }


    std::unordered_map<int, std::vector<LOverlap *> > idx3; // this is the pileup
    std::vector<std::set<int> > has_overlap(n_read);
    std::unordered_map<int, std::unordered_map<int, std::vector<LOverlap *> > > idx;


    for (int i = 0; i < n_read; i++) {
        //has_overlap[i] = std::set<int>();
        idx3[i] = std::vector<LOverlap *>();
    }

    //for (int i = 0; i < aln.size(); i++)
    //    if (aln[i]->active)
    //        idx[std::pair<int, int>(aln[i]->aid, aln[i]->bid)] = std::vector<LOverlap *>();
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


    std::cout << "add data" << std::endl;
    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx[aln[i]->read_A_id_][aln[i]->read_B_id_].push_back(aln[i]);
        }
    }
    std::cout << "add data" << std::endl;


    std::string name_input= out + ".edges.list";
    std::ifstream edges_file(name_input);

    std::string name_output = out_name + ".fasta";
    std::ofstream out_fa(name_output);

    int num_contig = 0;
    int num_one_read_contig = 0;
    while (true) {
        if (edges_file.eof()) break;
        edgelist.clear();
        std::string edge_line;
        while (!edges_file.eof()) {
            std::getline(edges_file, edge_line);
            //std::cout << edge_line << std::endl;

            std::vector<std::string> tokens = split(edge_line, ' ');

            if (tokens.size() == 1) {
                break;
            }
            //std::cout << tokens.size() << std::endl;

            Node node0;
            Node node1;
            int w;
            if (tokens.size() > 5 ) {
                node0.id = std::stoi(tokens[0]);
                node0.strand = std::stoi(tokens[1]);

                node1.id = std::stoi(tokens[2]);
                node1.strand = std::stoi(tokens[3]);;

                w = std::stoi(tokens[4]);
                edgelist.push_back(std::make_tuple(node0, node1, w));
            }


            if (tokens.size() == 4) {
                out_fa << ">OneReadContig" << num_one_read_contig << std::endl;



                int node_id = std::stoi(tokens[0]);
                int node_strand = std::stoi(tokens[1]);
                int from = std::stoi(tokens[2]);
                int to = std::stoi(tokens[3]);


                std::string current_seq;


                if (node_strand == 0) current_seq = reads[node_id]->bases;
                else current_seq = reverse_complement(reads[node_id]->bases);

                out_fa << current_seq.substr(from, to-from) << std::endl;

                num_one_read_contig++;
            }
        }

        std::cout << "list size:" << edgelist.size() << std::endl;
        if (edgelist.size() == 0) continue;


        std::vector<LAlignment *> full_alns;
        std::vector<LAlignment *> selected;
        std::unordered_map<int, std::vector<LAlignment *>> idx_aln;
        la.resetAlignment();
        std::vector<int> range;

        for (int i = 0; i < edgelist.size(); i++) {
            range.push_back(std::get<0>(edgelist[i]).id);
            idx_aln[std::get<0>(edgelist[i]).id] = std::vector<LAlignment *>();
        }

        std::sort(range.begin(), range.end());

        la.getAlignment(full_alns, range);

        for (auto i:full_alns) {
            idx_aln[i->read_A_id_].push_back(i);
        }

        for (int i = 0; i < edgelist.size(); i++) {
            int aid = std::get<0>(edgelist[i]).id;
            int bid = std::get<1>(edgelist[i]).id;
            bool found = false;
            for (int j = 0; j < idx_aln[std::get<0>(edgelist[i]).id].size(); j++) {
                //printf("%d %d %d %d\n",bid, idx_aln[aid][j]->bid, idx_aln[aid][j]->read_A_match_end_ - idx_aln[aid][j]->read_A_match_start_, std::get<2>(edgelist[i]));
                if ((idx_aln[aid][j]->read_B_id_ == bid) and \
            (idx_aln[aid][j]->aepos - idx_aln[aid][j]->abpos + idx_aln[aid][j]->bepos - idx_aln[aid][j]->bbpos == std::get<2>(edgelist[i]))) {
                    selected.push_back(idx_aln[aid][j]);
                    found = true;
                    break;
                }
                if (found) continue;
            }
        }

        std::cout << "selected:" << selected.size() << std::endl;

        std::unordered_map<int, std::unordered_map<int, std::pair<std::string, std::string> > > aln_tags_map;
        std::vector<std::pair<std::string, std::string> > aln_tags_list;
        std::vector<std::pair<std::string, std::string> > aln_tags_list_true_strand;


        for (int i = 0; i < selected.size(); i++) {
            la.recoverAlignment(selected[i]);
            //printf("%d %d\n",selected[i]->tlen, selected[i]->trace_pts_len);
            std::pair<std::string, std::string> res = la.getAlignmentTags(selected[i]);
            aln_tags_map[selected[i]->read_A_id_][selected[i]->read_B_id_] = res;
            aln_tags_list.push_back(res);
        }



        std::string sequence = "";

        std::vector<LOverlap *> bedges;
        std::vector<std::string> breads;

        std::vector<std::vector<std::pair<int, int> > > pitfalls;


        range.clear();
        for (int i = 0; i < edgelist.size(); i++) {
            range.push_back(std::get<0>(edgelist[i]).id);
        }

        std::vector<std::vector<int> *> coverages;

        for (int i = 0; i < range.size(); i++) {
            int aread = range[i];
            if (idx3[aread].size() > 0) {
                std::vector<int> *res = la.getCoverage(idx3[aread]);
                std::vector<std::pair<int, int> > *res2 = la.lowCoverageRegions(*res, MIN_COV2);
                //delete res;
                coverages.push_back(res);
                //printf("%d %d: (%d %d) ", i, aread, 0, idx3[aread][0]->alen);
                //for (int j = 0; j < res2->size(); j++) {
                //    printf("[%d %d] ", res2->at(j).first, res2->at(j).second);
                //}
                //printf("\n");
                pitfalls.push_back(*res2);
                delete res2;
            }
        }


        /***
         * Prepare the data
         */

        for (int i = 0; i < edgelist.size(); i++) {

            std::vector<LOverlap *> currentalns = idx[std::get<0>(edgelist[i]).id][std::get<1>(edgelist[i]).id];

            LOverlap *currentaln = NULL;

            for (int j = 0; j < currentalns.size(); j++) {
                //std::cout << std::get<0>(edgelist[i]).id << " " << std::get<1>(edgelist[i]).id << " " << currentalns[j]->match_type_ << std::endl;
                if (currentalns[j]->read_A_match_end_ - currentalns[j]->read_A_match_start_ + currentalns[j]->read_B_match_end_ - currentalns[j]->read_B_match_start_ ==
                    std::get<2>(edgelist[i]))
                    currentaln = currentalns[j];
            }

            if (currentaln == NULL) exit(1);
            //currentaln->show();

            std::string current_seq;
            std::string next_seq;

            std::string aln_tags1;
            std::string aln_tags2;


            if (std::get<0>(edgelist[i]).strand == 0)
                current_seq = reads[std::get<0>(edgelist[i]).id]->bases;
            else
                current_seq = reverse_complement(reads[std::get<0>(edgelist[i]).id]->bases);

            if (std::get<0>(edgelist[i]).strand == 0) {
                aln_tags1 = aln_tags_list[i].first;
                aln_tags2 = aln_tags_list[i].second;
            } else {
                aln_tags1 = reverse_complement(aln_tags_list[i].first);
                aln_tags2 = reverse_complement(aln_tags_list[i].second);
            }

            aln_tags_list_true_strand.push_back(std::pair<std::string, std::string>(aln_tags1, aln_tags2));

            if (std::get<1>(edgelist[i]).strand == 0)
                next_seq = reads[std::get<1>(edgelist[i]).id]->bases;
            else
                next_seq = reverse_complement(reads[std::get<1>(edgelist[i]).id]->bases);

            int abpos, aepos, alen, bbpos, bepos, blen, aes, aee, bes, bee;

            alen = currentaln->alen;
            blen = currentaln->blen;


            if (std::get<0>(edgelist[i]).strand == 0) {
                abpos = currentaln->read_A_match_start_;
                aepos = currentaln->read_A_match_end_;

                aes = currentaln->eff_read_A_start_;
                aee = currentaln->eff_read_A_end_;

            } else {
                abpos = alen - currentaln->read_A_match_end_;
                aepos = alen - currentaln->read_A_match_start_;

                aes = alen - currentaln->eff_read_A_end_;
                aee = alen - currentaln->eff_read_A_start_;
            }

            if (((std::get<1>(edgelist[i]).strand == 0))) {
                bbpos = currentaln->read_B_match_start_;
                bepos = currentaln->read_B_match_end_;

                bes = currentaln->eff_read_B_start_;
                bee = currentaln->eff_read_B_end_;

            } else {
                bbpos = blen - currentaln->read_B_match_end_;
                bepos = blen - currentaln->read_B_match_start_;

                bes = blen - currentaln->eff_read_B_end_;
                bee = blen - currentaln->eff_read_B_start_;

            }
            aes = 0;
            bes = 0;
            aee = alen;
            bee = blen;

            printf("%d %d [[%d %d] << [%d %d]] x [[%d %d] << [%d %d]]\n", std::get<0>(edgelist[i]).id, std::get<1>(edgelist[i]).id, abpos, aepos, aes, aee, bbpos, bepos, bes, bee);

            LOverlap *new_ovl = new LOverlap();
            new_ovl->read_A_match_start_ = abpos;
            new_ovl->read_A_match_end_ = aepos;
            new_ovl->read_B_match_start_ = bbpos;
            new_ovl->read_B_match_end_ = bepos;
            new_ovl->eff_read_A_end_ = aee;
            new_ovl->eff_read_A_start_ = aes;
            new_ovl->eff_read_B_end_ = bee;
            new_ovl->eff_read_B_start_ = bes;
            new_ovl->alen = currentaln->alen;
            new_ovl->blen = currentaln->blen;
            new_ovl->read_A_id_ = std::get<0>(edgelist[i]).id;
            new_ovl->read_B_id_ = std::get<1>(edgelist[i]).id;


            bedges.push_back(new_ovl);
            breads.push_back(current_seq);


        }
        //need to trim the end



        std::vector<std::vector<int> > mappings;
        for (int i = 0; i < range.size(); i++) {
            mappings.push_back(get_mapping(aln_tags_list_true_strand[i].first, aln_tags_list_true_strand[i].second));
        }

        std::cout << bedges.size() << " " << breads.size() << " " << selected.size() << " "
        << aln_tags_list.size() << " " << pitfalls.size() << " " << aln_tags_list_true_strand.size()
        << " " << mappings.size() << " " << coverages.size() << std::endl;

        /*for (int i = 0; i < bedges.size() - 1; i++) {
            printf("%d %d %d %d %d\n", bedges[i]->read_B_match_start_, bedges[i]->read_B_match_end_, bedges[i+1]->read_A_match_start_, bedges[i+1]->read_A_match_end_, bedges[i]->read_B_match_end_ - bedges[i+1]->read_A_match_start_);
        }*/


        int tspace = TSPACE; // set lane length to be 500
        int nlane = 0;


        //printf("%d %d\n", mappings[0][800], mappings[0][1000]); // debug output
        //printf("%s\n%s\n", breads[0].substr(bedges[0]->read_A_match_start_ + 800, 50).c_str(),
        //       breads[1].substr(bedges[0]->read_B_match_start_ + mappings[0][800], 50).c_str()); //debug output


        std::vector<std::vector<std::pair<int, int>>> lanes;

        std::string draft_assembly = "";


        int currentlane = 0;
        int current_starting_read = 0;
        int current_starting_space = 1;
        int current_starting_offset = 0;
        int n_bb_reads = range.size();
        std::vector<std::vector<int>> trace_pts(n_bb_reads);
        bool revert = false;


        int rmax = -1;
        /**
         * Move forward and put "trace points"
         */
        while (current_starting_read < n_bb_reads - 1) {
            int currentread = current_starting_read;
            int additional_offset = 0;
            while (bedges[current_starting_read]->read_A_match_start_ + current_starting_space * tspace +
                   current_starting_offset + additional_offset <
                   bedges[current_starting_read]->read_A_match_end_ - EDGE_SAFE) {
                int waypoint = bedges[current_starting_read]->read_A_match_start_ + tspace * current_starting_space +
                               current_starting_offset + additional_offset;
                //if ((waypoint - bedges[current_starting_read]->read_A_match_start_) < EDGE_SAFE)
                //    waypoint += EDGE_SAFE;

                //int next_waypoint = mappings[currentread][waypoint - bedges[current_starting_read]->read_A_match_start_] + bedges[current_starting_read]->read_B_match_start_;
                std::vector<std::pair<int, int> > lane;

                while ((waypoint > bedges[currentread]->read_A_match_start_) and
                       (waypoint < bedges[currentread]->read_A_match_end_)) {

                    printf("%d %d\n", currentread, waypoint);
                    trace_pts[currentread].push_back(waypoint);


                    /*if (waypoint > bedges[currentread]->read_A_match_end_ - EDGE_SAFE) {
                        printf("Reaching the end, neglect low coverage\n");
                    }

                    if ((coverages[currentread]->at(waypoint) < MIN_COV2) and (waypoint < bedges[currentread]->read_A_match_end_ - EDGE_SAFE)) {
                        revert = true;
                        printf("Low coverage, revert\n");
                        break;
                    }*/


                    lane.push_back(std::pair<int, int>(currentread, waypoint));
                    if (currentread > rmax) rmax = currentread;
                    //int previous_wp = waypoint;
                    waypoint = mappings[currentread][waypoint - bedges[currentread]->read_A_match_start_] +
                               bedges[currentread]->read_B_match_start_;
                    //printf("%s\n%s\n", breads[currentread].substr(previous_wp,50).c_str(), breads[currentread+1].substr(waypoint,50).c_str());
                    currentread++;
                    if (currentread >= n_bb_reads) break;
                }
                if (currentread < n_bb_reads) if (waypoint < bedges[currentread]->alen) {
                    lane.push_back(std::pair<int, int>(currentread, waypoint));
                    if (currentread > rmax) rmax = currentread;
                }
                /*if (revert) {
                    printf("revert\n");
                    revert = false;
                    while (currentread >= current_starting_read) {
                        trace_pts[currentread].pop_back();
                        currentread --;
                        additional_offset += STEP;
                    }
                    currentread = current_starting_read;
                }
                else*/
                {
                    if (currentread >= rmax)
                        lanes.push_back(lane);
                    current_starting_space++;
                    currentread = current_starting_read;

                }

            }

            current_starting_read++;
            current_starting_space = 1;//get next space;
            if (trace_pts[current_starting_read].size() == 0)
                current_starting_offset = 0;
            else
                current_starting_offset =
                        trace_pts[current_starting_read].back() - bedges[current_starting_read]->read_A_match_start_;
        }


        /**
         * Show trace points on reads
         */
        for (int i = 0; i < n_bb_reads; i++) {
            printf("Read %d:", i);
            for (int j = 0; j < trace_pts[i].size(); j++) {
                printf("%d ", trace_pts[i][j]);
            }
            printf("\n");
        }

        /**
         * Show lanes
         */

        for (int i = 0; i < lanes.size(); i++) {

            printf("Lane %d\n", i);
            for (int j = 0; j < lanes[i].size(); j++) {
                printf("[%d %d] ", lanes[i][j].first, lanes[i][j].second);
            }
            printf("\n");
        }


        printf("In total %d lanes\n", lanes.size());
        if (lanes.size() == 0) {
            draft_assembly = breads[0];
            out_fa << ">DraftAssemblyContig" << num_contig << std::endl;
            out_fa << draft_assembly << std::endl;
            num_contig++;
            continue;
        }



        /**
         * Consequtive lanes form a column (ladder)
         */

        std::vector<std::vector<std::tuple<int, int, int> > > ladders;

        for (int i = 0; i < lanes.size() - 1; i++) {
            std::vector<std::pair<int, int> > lane1 = lanes[i];
            std::vector<std::pair<int, int> > lane2 = lanes[i + 1];
            std::vector<std::tuple<int, int, int> > ladder;
            int pos = 0;
            for (int j = 0; j < lane2.size(); j++) {
                while ((lane1[pos].first != lane2[j].first) and (pos < lane1.size() - 1)) pos++;
                if ((lane1[pos].first == lane2[j].first))
                    ladder.push_back(std::make_tuple(lane2[j].first, lane1[pos].second, lane2[j].second));
            }
            ladders.push_back(ladder);
        }


        /**
         * show ladders
         */
        for (int i = 0; i < ladders.size(); i++) {
            printf("Ladder %d\n", i);
            for (int j = 0; j < ladders[i].size(); j++) {
                //printf("[%d %d-%d] ", std::get<0>(ladders[i][j]), std::get<1>(ladders[i][j]), std::get<2>(ladders[i][j]) );
                //printf("%s\n", breads[std::get<0>(ladders[i][j])].substr(std::get<1>(ladders[i][j]),std::get<2>(ladders[i][j])-std::get<1>(ladders[i][j])).c_str());

            }

            if (ladders[i].size() == 0) {
                printf("low coverage!\n");
                continue;
            }

            if (ladders[i].size() > 1) {


                int mx = 0;
                int maxcoverage = 0;
                for (int j = 0; j < ladders[i].size(); j++) {
                    int mincoverage = 10000;
                    int read = std::get<0>(ladders[i][j]);
                    int start = std::get<1>(ladders[i][j]);
                    int end = std::get<2>(ladders[i][j]);
                    for (int pos = start; pos < end; pos++) {
                        if (coverages[read]->at(pos) < mincoverage) mincoverage = coverages[read]->at(pos);
                    }
                    if (mincoverage > maxcoverage) {
                        maxcoverage = mincoverage;
                        mx = j;
                    }
                }

                std::cout << "ladder " << i << " num reads " << ladders[i].size() << " possibly error here " <<
                maxcoverage << "\n!";


                //if (ladders[i].size() == 2) {
                //    draft_assembly += breads[std::get<0>(ladders[i][mx])].substr(std::get<1>(ladders[i][mx]),
                //                                                                 std::get<2>(ladders[i][mx]) -
                //                                                                 std::get<1>(ladders[i][mx]));
                //    continue;
               // }


                std::string base = breads[std::get<0>(ladders[i][mx])].substr(std::get<1>(ladders[i][mx]),
                                                                              std::get<2>(ladders[i][mx]) -
                                                                              std::get<1>(ladders[i][mx]));;
                int seq_count = ladders[i].size();
                printf("seq_count:%d, max %d\n", seq_count, mx);
                align_tags_t **tags_list;
                tags_list = (align_tags_t **) calloc(seq_count, sizeof(align_tags_t *));
                consensus_data *consensus;

                int alen = (std::get<2>(ladders[i][mx]) - std::get<1>(ladders[i][mx]));
                for (int j = 0; j < ladders[i].size(); j++) {

                    int blen = (std::get<2>(ladders[i][j]) - std::get<1>(ladders[i][j]));
                    char *aseq = (char *) malloc(
                            (20 + (std::get<2>(ladders[i][mx]) - std::get<1>(ladders[i][mx]))) * sizeof(char));
                    char *bseq = (char *) malloc(
                            (20 + (std::get<2>(ladders[i][j]) - std::get<1>(ladders[i][j]))) * sizeof(char));
                    strcpy(aseq, breads[std::get<0>(ladders[i][mx])].substr(std::get<1>(ladders[i][mx]),
                                                                            std::get<2>(ladders[i][mx]) -
                                                                            std::get<1>(ladders[i][mx])).c_str());
                    strcpy(bseq, breads[std::get<0>(ladders[i][j])].substr(std::get<1>(ladders[i][j]),
                                                                           std::get<2>(ladders[i][j]) -
                                                                           std::get<1>(ladders[i][j])).c_str());


                    aln_range *arange = (aln_range *) calloc(1, sizeof(aln_range));
                    arange->s1 = 0;
                    arange->e1 = strlen(bseq);
                    arange->s2 = 0;
                    arange->e2 = strlen(aseq);
                    arange->score = 5;

                    //printf("blen %d alen%d\n",strlen(bseq), strlen(aseq));
                    //printf("before get tags\n");

                    alignment *alng = _align(bseq, blen, aseq, alen, 150, 1);

                    char *q_aln_str = (char *) malloc((5 + strlen(alng->q_aln_str)) * sizeof(char));
                    char *t_aln_str = (char *) malloc((5 + strlen(alng->t_aln_str)) * sizeof(char));


                    strcpy(q_aln_str + 1, alng->q_aln_str);
                    strcpy(t_aln_str + 1, alng->t_aln_str);
                    q_aln_str[0] = 'T';
                    t_aln_str[0] = 'T';


                    for (int pos = 0; pos < strlen(q_aln_str); pos++) q_aln_str[pos] = toupper(q_aln_str[pos]);
                    for (int pos = 0; pos < strlen(t_aln_str); pos++) t_aln_str[pos] = toupper(t_aln_str[pos]);

                    //printf("Q:%s\nT:%s\n", q_aln_str, t_aln_str);

                    tags_list[j] = get_align_tags(q_aln_str,
                                                  t_aln_str,
                                                  strlen(alng->q_aln_str) + 1,
                                                  arange, (unsigned int) j, 0);
                    //free(aseq);
                    //free(bseq);

                    /*for (int k = 0; k < tags_list[j]->len; k++) {
                        printf("%d %d %ld %d %c %c\n",j, k, tags_list[j]->align_tags[k].t_pos,
                               tags_list[j]->align_tags[k].delta,
                                //tags_list[j]->align_tags[k].p_q_base,
                               aseq[tags_list[j]->align_tags[k].t_pos],
                               tags_list[j]->align_tags[k].q_base);
                    }*/
                    free(q_aln_str);
                    free(t_aln_str);
                    free(aseq);
                    free(bseq);
                    free_alignment(alng);

                }

                //printf("%d %d\n%s\n",seq_count, strlen(seq), seq);

                consensus = get_cns_from_align_tags(tags_list, seq_count, alen + 1, 1);
                printf("Consensus len :%d\n",strlen(consensus->sequence));
                draft_assembly += std::string(consensus->sequence);

                free_consensus_data(consensus);
                for (int j = 0; j < seq_count; j++)
                    free_align_tags(tags_list[j]);

            } else {
                draft_assembly += breads[std::get<0>(ladders[i][0])].substr(std::get<1>(ladders[i][0]),
                                                                            std::get<2>(ladders[i][0]) -
                                                                            std::get<1>(ladders[i][0]));
            }

            printf("\n");
        }



        /*for (int i = 0; i < mapping.size(); i++)
            printf("%d %d\n", i, mapping[i]);
        printf("[%d %d], [%d %d]\n", bedges[0]->read_A_match_start_, bedges[0]->read_A_match_end_, bedges[0]->read_B_match_start_, bedges[0]->read_B_match_end_);*/

        std::cout << sequence.size() << std::endl;
        std::cout << draft_assembly.size() << std::endl;


        out_fa << ">Draft_assembly" << num_contig << std::endl;
        out_fa << draft_assembly << std::endl;
        num_contig++;

    }

    if (strlen(name_db) > 0)
    la.closeDB(); //close database
    return 0;
}
