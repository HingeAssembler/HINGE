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


#define LAST_READ_SYMBOL  '$'

bool compare_overlap(LAlignment * ovl1, LAlignment * ovl2) {
    return ((ovl1->aepos - ovl1->abpos + ovl1->bepos - ovl1->bbpos) > (ovl2->aepos - ovl2->abpos + ovl2->bepos - ovl2->bbpos));
}


static int ORDER(const void *l, const void *r) {
    int x = *((int32 *) l);
    int y = *((int32 *) r);
    return (x - y);
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


	la.showRead2(4,6);

    la.openAlignmentFile(name_las);

	int n_alns = la.getAlignmentNumber();

    std::cout<<"# Alignments:" << n_alns << std::endl;

	//la.showAlignment(0,1);

	std::vector<LOverlap *> res1;
	  	la.resetAlignment();
	  	la.getOverlap(res1, 0, 1); // get alignment(overlap) for reads [3,5)

	for (auto i:res1)
		i->show();
	printf("\n");


    std::vector<LAlignment *> res;
    la.resetAlignment();
	la.getAlignment(res, 0, n_alns); // get all alignments


	std::vector<std::vector<LAlignment *>> idx;

    printf("%d\n", res.size());

	for (int i = 0; i < n_contigs; i++)
		idx.push_back(std::vector<LAlignment *>());

	for (int i = 0; i < n_alns; i++) {
        idx[res[i]->read_A_id_].push_back(res[i]);
    }

	for (int i = 0; i < n_contigs; i++) {
        std::sort(idx[i].begin(), idx[i].end(), compare_overlap);
        printf("%d %d\n", i, idx[i].size());
    }

    for (int i = 0; i < n_contigs; i++) {

        int j = 0;
        for (j = 0;j < idx[i].size(); j++)
            if (idx[i][j]->aepos - idx[i][j]->abpos < 2000)
                break;




        int seq_count = j;
        printf("seq_count:%d\n",seq_count);

        align_tags_t ** tags_list;
        tags_list = (align_tags_t **) calloc( seq_count+1, sizeof(align_tags_t *) );

        test_read = la.getRead(i); //get read 0
        std::string base_structure = test_read->bases;

        std::transform(base_structure.begin(), base_structure.end(),base_structure.begin(), ::toupper);

        aln_range * arange;
        arange = (aln_range*) calloc(1 , sizeof(aln_range));
        arange->s1 = 0;
        arange->s2 = 0;
        arange->e1 = base_structure.size();
        arange->e2 = base_structure.size();

        char * seq = (char *) malloc(base_structure.size()* sizeof(char));
        strcpy(seq, base_structure.c_str());

        tags_list[0] = get_align_tags(  seq,
                                        seq,
                                        strlen(seq),
                                        arange, 0, 0);



        for (int j = 0; j < seq_count ; j ++) {


            la.recoverAlignment(idx[i][j]);
            std::pair<std::string, std::string>  alignment = la.getAlignmentTags(idx[i][j]);


            char * t_aln_str = (char *) malloc((alignment.first.size()+20)* sizeof(char));
            char * q_aln_str = (char *) malloc((alignment.second.size()+20)* sizeof(char));

            strcpy(q_aln_str, alignment.second.c_str());
            strcpy(t_aln_str, alignment.first.c_str());


            //printf("%s\n%s\n",q_aln_str,t_aln_str);


            seq_coor_t aln_str_size = strlen(q_aln_str);
            if (j%1000 == 0) printf("%d/%d, %d\n", j, seq_count, aln_str_size);

            aln_range * arange = (aln_range*) calloc(1 , sizeof(aln_range));
            arange->s1 = idx[i][j]->bbpos;
            arange->e1 = idx[i][j]->bepos;
            arange->s2 = idx[i][j]->abpos;
            arange->e2 = idx[i][j]->aepos;
            arange->score = 5;
            //printf("before get tags\n");
            tags_list[j+1] = get_align_tags( q_aln_str,
                                             t_aln_str,
                                             aln_str_size,
                                             arange, (unsigned int)j + 1, 0);

            //printf("after get tags\n");

            if (q_aln_str!=NULL) free(q_aln_str);
            if (t_aln_str!=NULL) free(t_aln_str);
            if (arange!=NULL) free_aln_range(arange);
            q_aln_str = NULL;
            t_aln_str = NULL;
            arange = NULL;
        }

        consensus_data * consensus;
        consensus = get_cns_from_align_tags_large( tags_list, seq_count+1, strlen(seq), 6 );

        printf("Length:%d\n", strlen(consensus->sequence));
        out << ">Consensus" << i << std::endl << consensus->sequence<<std::endl;
        free_consensus_data(consensus);
        for (int j = 0; j <seq_count + 1; j++)
            free_align_tags(tags_list[j]);

    }


    la.closeDB(); //close database*/
    return 0;
}
