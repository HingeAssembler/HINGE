#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"
#include "LAInterface.h"

#include <iostream>
#include <fstream>
#include <tuple>
#include <string>
#include <algorithm>
#include <map>
#include <unordered_map>
extern "C" {
#include "common.h"
}
#include "INIReader.h"


static char ToU[4] = { 'A', 'C', 'G', 'T' };


char toLower(char c) {

    char base = c;

    switch (c) {
        case 'A': base = 'a'; break;
        case 'C': base = 'c'; break;
        case 'G': base = 'g'; break;
        case 'T': base = 't'; break;
    }

    return base;

}


int main(int argc, char *argv[]) {

    std::string name_db1 = std::string(argv[1]);
    std::string name_db2 = std::string(argv[2]);
    std::string name_las = std::string(argv[3]);
	char * name_out = argv[4];
	char * name_config = argv[5];

    std::ofstream out(name_out);

	INIReader reader(name_config);

	if (reader.ParseError() < 0) {
	    std::cout << "Can't load "<<name_config<<std::endl;
	    return 1;
	}


	int LENGTH_THRESHOLD = reader.GetInteger("consensus", "min_length", -1);
	printf("length threshold:%d\n", LENGTH_THRESHOLD);
    //std::cout<<name_db1 << " " << name_db2 << " " << name_las <<std::endl;

    LAInterface la;
    //std::cout << "hello" << std::endl;
    Read *test_read;

    la.openDB2(name_db1, name_db2);

	int n_reads = la.getReadNumber2();
    int n_contigs = la.getReadNumber();

    std::cout<<"# Contigs:" << n_contigs << std::endl;
    std::cout<<"# Reads:" << n_reads << std::endl;


    la.openAlignmentFile(name_las);

	int n_alns = la.getAlignmentNumber();

    std::cout<<"# Alignments:" << n_alns << std::endl;


    std::vector<LAlignment *> res;
    la.resetAlignment();
	la.getAlignment(res, 0, n_alns); // get all alignments


	std::vector<std::vector<LAlignment *>> idx;

    printf("%lu\n", res.size());

	for (int i = 0; i < n_contigs; i++)
		idx.push_back(std::vector<LAlignment *>());

	for (int i = 0; i < n_alns; i++) {
        idx[res[i]->read_A_id_].push_back(res[i]);
    }

	for (int i = 0; i < n_contigs; i++) {
        std::sort(idx[i].begin(), idx[i].end(), compare_overlap_aln);
        printf("%d %lu\n", i, idx[i].size());
    }

    std::cout << "Getting read lengths" << std::endl;
    std::vector<Read *> reads_vec;
    la.getRead(reads_vec, 0, n_contigs);
    for (int i = 0; i < n_contigs; i++){
        std::cout << i << "\t" << (reads_vec[i]->bases).size() << std::endl;
    }

    std::cout << "Building consensus sequences..." << std::endl;

    for (int i = 0; i < n_contigs; i++) {



        int k = 0;
        for (k = 0; k < idx[i].size(); k++)
            if (idx[i][k]->aepos - idx[i][k]->abpos < LENGTH_THRESHOLD)
                break;

        int seq_count = k;

        std::cout << "Contig " << i << ": " << seq_count << " reads" << std::endl;

        if (seq_count == 0) {
            out << ">Consensus" << i << std::endl;
            out << reads_vec[i]->bases << std::endl;
            continue;
        }

        std::vector<std::vector<int>> contig_base_scores;

        std::vector<int> insertion_score (idx[i][0]->alen,0);
        std::vector<std::vector<int>> insertion_base_scores; // handling single insertions only

        std::vector<int> cov_depth (idx[i][0]->alen,0);

        std::vector<int> zero_scores (5,0); // scores for A,C,G,T,- are initialized at 0
        for (int j = 0; j < idx[i][0]->alen; j++) {
            contig_base_scores.push_back(zero_scores);
            insertion_base_scores.push_back(zero_scores);
        }

        for (int j = 0; j < seq_count ; j ++) {

            la.recoverAlignment(idx[i][j]);
            std::pair<std::string, std::string>  alignment = la.getAlignmentTags(idx[i][j]);

            int pos_in_contig = idx[i][j]->abpos;
            

            for (int m = 0; m < alignment.first.length(); m++) {

                int base = -1;
                switch (alignment.second[m]) {
                    case 'A': base = 0; break;
                    case 'C': base = 1; break;
                    case 'G': base = 2; break;
                    case 'T': base = 3; break;
                    case '-': base = 4; break;
                }

                if (alignment.first[m] != '-') {

                    if (base != -1) {
                        contig_base_scores[pos_in_contig][base]++;
                        cov_depth[pos_in_contig]++;
                    }

                    pos_in_contig++;
                }
                else if (base != -1) {
                    insertion_score[pos_in_contig]++;
                    insertion_base_scores[pos_in_contig][base]++;
                }

            }

        }

        int good_bases = 0;
        int insertions = 0; // insertion here means that a base is inserted in the consensus
        int deletions = 0; // deletion here means that the base from the draft is deleted in the consensus

        int consensus_length = 0;

        int low_coverage_bases = 0;

        long int sum_coverage = 0;

        out << ">Consensus" << i << std::endl;


        for (int j=0; j < idx[i][0]->alen ; j++) {

            sum_coverage += cov_depth[j];

            if (cov_depth[j] < 3) {
//                std::cout << "Low coverage." << std::endl;

                low_coverage_bases++;
                out << toLower(reads_vec[i]->bases[j]);
                continue;
            }

            if (insertion_score[j] > cov_depth[j]/2) {
                int max_insertion_base = 0;
                for (int b=1; b<4; b++) {
                    if (insertion_base_scores[j][b] > insertion_base_scores[j][max_insertion_base]) max_insertion_base = b;
                }
                out << ToU[max_insertion_base];
                consensus_length++;
                insertions++;
            }

            int max_base = 0;

            for (int b=1; b<5; b++) {
                if (contig_base_scores[j][b] > contig_base_scores[j][max_base]) max_base = b;
            }

            if (max_base < 4) {
                out << ToU[max_base];
                good_bases++;
                consensus_length++;
            }
            else {
                deletions++;
            }


        }
        out << std::endl;


        printf("Average coverage: %f\n",(1.0*sum_coverage)/idx[i][0]->alen);
        printf("Good bases: %d/%d\n",good_bases,idx[i][0]->alen);
        printf("Insertions: %d/%d\n",insertions,idx[i][0]->alen);
        printf("Deletions: %d/%d\n",deletions,idx[i][0]->alen);
        printf("Low coverage bases: %d/%d\n",low_coverage_bases,idx[i][0]->alen);
        printf("Consensus length: %d\n",consensus_length);


    }




    la.closeDB(); //close database*/
    return 0;
}
