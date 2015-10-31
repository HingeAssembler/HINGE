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
#include <tuple>
#include <string>
#include <algorithm>
extern "C" {
#include "common.h"
}

#define LAST_READ_SYMBOL  '$'

static int ORDER(const void *l, const void *r) {
    int x = *((int32 *) l);
    int y = *((int32 *) r);
    return (x - y);
}

int main(int argc, char *argv[]) {
    LAInterface la;
    std::cout << "hello" << std::endl;
    Read *test_read;

    la.OpenDB("G");
    std::cout<<"# Reads:" << la.getReadNumber() << std::endl;

    la.showRead(1, 3); //show read [1,3)

    test_read = la.getRead(0); //get read 0
    test_read->showRead(); // show read 0

    la.OpenAlignment("G.1.las");
    la.showAlignment(0, 2); // show alignments of read [0,2)

    std::cout<<"# Alignments:" << la.getAlignmentNumber() << std::endl;

    la.resetAlignment();
    std::vector<int> res;
    la.getAlignmentB(res, 1); //get alignment for read 1
    for (auto i:res)
        printf("%d ", i);
    printf("\n");

    std::vector<LOverlap *> res1;
    la.resetAlignment();
    la.getOverlap(res1, 3, 5); // get alignment(overlap) for reads [3,5)
	
	for (auto i:res1)
		i->show();
	printf("\n");

    std::vector<LAlignment *> res2;
    la.resetAlignment();
    //la.getAlignment(res2, 0, 3);// get alignment for reads [0,3)
	
	la.getAlignment(res2, 0, 1);


    int seq_count = res2.size();
    align_tags_t ** tags_list;
    tags_list = (align_tags_t **) calloc( seq_count+1, sizeof(align_tags_t *) );


    test_read = la.getRead(0); //get read 0
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


    for (int i = 0; i < seq_count; i ++) {
        res2[i]->show();

        std::pair<std::string, std::string>  alignment = la.Lget_Alignment_tgs(res2[i]);

        //alignment.first.erase (std::remove(alignment.first.begin(), alignment.first.end(), '-'), alignment.first.end());
        //alignment.second.erase (std::remove(alignment.second.begin(), alignment.second.end(), '-'), alignment.second.end());

        //std::cout << alignment.first.size() <<alignment.first << std::endl;
        //std::cout <<  alignment.second.size() <<alignment.second << std::endl;

        char * t_aln_str = (char *) malloc(alignment.first.size()* sizeof(char));
        char * q_aln_str = (char *) malloc(alignment.second.size()* sizeof(char));

        strcpy(q_aln_str, alignment.second.c_str());
        strcpy(t_aln_str, alignment.first.c_str());

        seq_coor_t aln_str_size = strlen(q_aln_str);
        aln_range * arange;
        arange = (aln_range*) calloc(1 , sizeof(aln_range));
        arange->s1 = res2[i]->bbpos;
        arange->e1 = res2[i]->bepos;
        arange->s2 = res2[i]->abpos;
        arange->e2 = res2[i]->aepos;

        tags_list[i+1] = get_align_tags( q_aln_str,
                                         t_aln_str,
                                         aln_str_size,
                                         arange, (unsigned int)i + 1, 0);

        free(q_aln_str);
        free(t_aln_str);
        free_aln_range(arange);

    }

    //print consensus

    consensus_data * consensus;
    consensus = get_cns_from_align_tags( tags_list, seq_count+1, strlen(seq), 6 );

    printf("Consensus:%s\n", consensus->sequence);

    free_consensus_data(consensus);
    for (int i = 0; i <seq_count + 1; i++)
        free_align_tags(tags_list[i]);

    la.CloseDB(); //close database
    return 0;
}