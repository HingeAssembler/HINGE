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

    char * filename = argv[1];
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

    std::string first_seq = reads_list.front();
    char * first_seq_cstr = (char *) malloc(first_seq.size() * sizeof(char));
    strcpy(first_seq_cstr, first_seq.c_str());

    std::string second_seq = reads_list[1];
    char * second_seq_cstr = (char *) malloc(second_seq.size() * sizeof(char));
    strcpy(second_seq_cstr, second_seq.c_str());


    first_seq_cstr ++; //there is a white space
    second_seq_cstr ++;

    printf("%s\n",first_seq_cstr);
    kmer_lookup * lk_ptr;
    seq_array sa_ptr;
    seq_addr_array sda_ptr;
    kmer_match * kmer_match_ptr;

    int K = 8;
    lk_ptr = allocate_kmer_lookup( 1 << (K * 2) );
    sa_ptr = allocate_seq( (seq_coor_t) first_seq.size() ); // seed sequence
    sda_ptr = allocate_seq_addr( (seq_coor_t) first_seq.size() );

    printf("fine %d\n", strlen(first_seq_cstr));

    add_sequence( 0, K, first_seq_cstr, strlen(first_seq_cstr), sda_ptr, sa_ptr, lk_ptr);


    kmer_match_ptr = find_kmer_pos_for_seq(second_seq_cstr, strlen(second_seq_cstr), K, sda_ptr, lk_ptr);



    for (int i = 0; i <  kmer_match_ptr->count; i++ ) {
        printf("%d %d\n", kmer_match_ptr->query_pos[i], kmer_match_ptr->target_pos[i]);
    }

#define INDEL_ALLOWENCE_0 6
    aln_range * arange;

    arange = find_best_aln_range(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels
    printf("%d, %d, %d, %d\n", arange->s1, arange->e1, arange->s2, arange->e2);


    return 0;
}

