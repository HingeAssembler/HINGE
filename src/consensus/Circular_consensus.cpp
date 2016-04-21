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

extern "C" {
#include "common.h"
}



#define LAST_READ_SYMBOL  '$'

typedef std::tuple<Node, Node, int> Edge_w;

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
        INSERT_ELEMENT(MISMATCH_LEFT);
        INSERT_ELEMENT(MISMATCH_RIGHT);
        INSERT_ELEMENT(COVERED);
        INSERT_ELEMENT(COVERING);
        INSERT_ELEMENT(UNDEFINED);
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
    return ((ovl1->read_A_match_end_ - ovl1->read_A_match_start_ + ovl1->read_B_match_end_ - ovl1->read_B_match_start_) > (ovl2->read_A_match_end_ - ovl2->read_A_match_start_ + ovl2->read_B_match_end_ - ovl2->read_B_match_start_));
}

bool compare_sum_overlaps(const std::vector<LOverlap * > * ovl1, const std::vector<LOverlap *> * ovl2) {
    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < ovl1->size(); i++) sum1 += (*ovl1)[i]->read_A_match_end_ - (*ovl1)[i]->read_A_match_start_ + (*ovl1)[i]->read_B_match_end_ - (*ovl1)[i]->read_B_match_start_;
    for (int i = 0; i < ovl2->size(); i++) sum2 += (*ovl2)[i]->read_A_match_end_ - (*ovl2)[i]->read_A_match_start_ + (*ovl2)[i]->read_B_match_end_ - (*ovl2)[i]->read_B_match_start_;
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

    int left= intervals[0]->read_A_match_start_ + cutoff, right = intervals[0]->read_A_match_end_ - cutoff; //left, right means maximal possible interval now

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

    char * filename = argv[1];
    char * out = argv[2];
    printf("filename: %s\n", filename);

    std::string line;
    std::ifstream myfile (filename);
    std::string str;
    getline (myfile,line);
    std::vector<std::string> reads_list;
    if (myfile.is_open()) {
        while ( getline (myfile,line) ) {
            //std::cout << line << '\n';
            if (line.front() == '>') {
                //std::cout << str << std::endl;
                reads_list.push_back(str);
                str = "";
            }
            else {
                str += line;
            }
        }
        reads_list.push_back(str);
        myfile.close();
    } else {
        std::cout << "Cannot open file!" << std::endl;
    }

    //for (auto i:reads_list)
    //    std::cout<<i<<std::endl;


    for (int i = 0; i < reads_list.size(); i ++ ) {
        std::transform(reads_list[i].begin(),reads_list[i].end(),reads_list[i].begin(), ::toupper);
    }

    std::string seq = reads_list[0];
    printf("%d\n",seq.size());

    for (int i = 0; i < reads_list.size() - 1; i++) {
        std::string first_seq = reads_list[i];
        char *first_seq_cstr = (char *) malloc((first_seq.size()+20) * sizeof(char));
        strcpy(first_seq_cstr, first_seq.c_str());

        std::string second_seq = reads_list[i + 1];
        char *second_seq_cstr = (char *) malloc((second_seq.size()+20) * sizeof(char));
        strcpy(second_seq_cstr, second_seq.c_str());


        first_seq_cstr++; //there is a white space
        second_seq_cstr++;

        //printf("%s\n", first_seq_cstr);
        kmer_lookup *lk_ptr;
        seq_array sa_ptr;
        seq_addr_array sda_ptr;
        kmer_match *kmer_match_ptr;

        int K = 8;
        lk_ptr = allocate_kmer_lookup(1 << (K * 2));
        sa_ptr = allocate_seq((seq_coor_t) first_seq.size()); // seed sequence
        sda_ptr = allocate_seq_addr((seq_coor_t) first_seq.size());

        //printf("fine %d\n", strlen(first_seq_cstr));

        add_sequence(0, K, first_seq_cstr, strlen(first_seq_cstr), sda_ptr, sa_ptr, lk_ptr);


        kmer_match_ptr = find_kmer_pos_for_seq(second_seq_cstr, strlen(second_seq_cstr), K, sda_ptr, lk_ptr);

        //for (int i = 0; i < kmer_match_ptr->count; i++) {
        //    printf("%d %d\n", kmer_match_ptr->query_pos[i], kmer_match_ptr->target_pos[i]);
        //}

#define INDEL_ALLOWENCE_0 6
        aln_range *arange;

        arange = find_best_aln_range(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0,
                                     5);  // narrow band to avoid aligning through big indels
        printf("%d: %d, %d, %d, %d\n", i, arange->s1, arange->e1, arange->s2, arange->e2);

        seq.erase(seq.end()- (reads_list[i].size() - arange->e2) , seq.end());
        std::cout << "erase:" <<(reads_list[i].size() - arange->e2) << std::endl;
        std::string apd = reads_list[i + 1];
        apd.erase(apd.begin(), apd.begin() + arange->e1);
        std::cout<<"grow:" << apd.size() << std::endl;


        std::cout<<seq.substr(seq.size() - 100, std::string::npos) << apd.substr(0,100) << std::endl;

        seq.append(apd);

        free(first_seq_cstr - 1);
        free(second_seq_cstr - 1);
        free_kmer_match(kmer_match_ptr);
        free_aln_range(arange);
        free_seq_addr_array(sda_ptr);
        free_seq_array(sa_ptr);
        free_kmer_lookup(lk_ptr);


        // Now it gives segment error, need to be fixed,
        // Update:fixed

    }

    //std::cout<<seq<<std::endl;
    std::cout<<seq.size()<<std::endl;


    std::ofstream outfile(out);

    outfile << ">Draft_assembly_" <<seq.size() <<std::endl;
    outfile<<seq<<std::endl;

    /*consensus_data * consensus;

    char ** input_seq;
    char ** seq_id;

    int seq_count = reads_list.size() + 1;

    input_seq = (char **)calloc( seq_count, sizeof(char *));
    seq_id = (char **)calloc( seq_count, sizeof(char *));


    input_seq[0] = (char *) malloc((seq.size() + 20) * sizeof(char));
    strcpy(input_seq[0], seq.c_str());
    seq_id[0] = (char *)calloc( 20, sizeof(char));
    strcpy(seq_id[0], "hi");

    for (int i = 0; i < seq_count-1; i++) {
        input_seq[i+1] = (char *) malloc((reads_list[i].size() + 20) * sizeof(char));
        strcpy(input_seq[i+1], reads_list[i].c_str());
        seq_id[i+1] = (char *)calloc( 20, sizeof(char));
        strcpy(seq_id[i+1], "hi");

    }

    printf("finish reading\n");

    consensus = generate_consensus_large(input_seq, seq_count - 1, 8, 8, 12, 1, 0.90);

    if (strlen(consensus->sequence) > 500) {
        printf(">%s\n%s\n", seq_id[0], consensus->sequence);
    }
    */


    return 0;
}

