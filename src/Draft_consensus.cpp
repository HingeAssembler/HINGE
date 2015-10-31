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

    std::string name_db1 = std::string(argv[1]);
    std::string name_db2 = std::string(argv[2]);
    std::string name_las = std::string(argv[3]);

    //std::cout<<name_db1 << " " << name_db2 << " " << name_las <<std::endl;

    LAInterface la;
    //std::cout << "hello" << std::endl;
    Read *test_read;

    la.OpenDB2(name_db1, name_db2);
    std::cout<<"# Contigs:" << la.getReadNumber() << std::endl;
    std::cout<<"# Reads:" << la.getReadNumber2() << std::endl;


	//la.showRead2(4,6);
	
    la.OpenAlignment(name_las);

    std::cout<<"# Alignments:" << la.getAlignmentNumber() << std::endl;
	
	//la.showAlignment(0,1);
	
	//std::vector<LOverlap *> res1;
	//  	la.resetAlignment();
	//  	la.getOverlap(res1, 0, 1); // get alignment(overlap) for reads [3,5)
	    	
	//for (auto i:res1)
	//	i->show();
	//printf("\n");
	
	
    std::vector<LAlignment *> res;
    la.resetAlignment();
	la.getAlignment(res, 0, 1); // get all alignments



    std::vector<LAlignment *> filtered;

    for (auto i:res) {
        if ((i->aepos - i->abpos) > 4000) {
            filtered.push_back(i);
        }
    }




    //for (auto i:res) {
	//	i->show();
	//}
	
	
    int seq_count = filtered.size();
	printf("seq_count:%d\n",seq_count);
	
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



	//printf("%s\n%s\n",seq,seq);




	
	printf("num:%d\n",seq_count);

	for (int i = 0; i < seq_count ; i ++) {


        std::pair<std::string, std::string>  alignment = la.Lget_Alignment_tgs(filtered[i]);

		
        char * t_aln_str = (char *) malloc((alignment.first.size()+20)* sizeof(char));
        char * q_aln_str = (char *) malloc((alignment.second.size()+20)* sizeof(char));

        strcpy(q_aln_str, alignment.second.c_str());
        strcpy(t_aln_str, alignment.first.c_str());


		//printf("%s\n%s\n",q_aln_str,t_aln_str);
		
	
        seq_coor_t aln_str_size = strlen(q_aln_str);
		printf("%d/%d, %d\n", i, seq_count, aln_str_size);
        
        aln_range * arange = (aln_range*) calloc(1 , sizeof(aln_range));
        arange->s1 = filtered[i]->bbpos;
        arange->e1 = filtered[i]->bepos;
        arange->s2 = filtered[i]->abpos;
        arange->e2 = filtered[i]->aepos;
		arange->score = 5;
		//printf("before get tags\n");
        tags_list[i+1] = get_align_tags( q_aln_str,
                                         t_aln_str,
                                         aln_str_size,
                                         arange, (unsigned int)i + 1, 0);
		
										 //printf("after get tags\n");
        
		if (q_aln_str!=NULL) free(q_aln_str);
        if (t_aln_str!=NULL) free(t_aln_str);
        if (arange!=NULL) free_aln_range(arange);
    }

    //print consensus

    consensus_data * consensus;
    consensus = get_cns_from_align_tags_large( tags_list, seq_count+1, strlen(seq), 6 );

    printf(">Consensus\n%s\n", consensus->sequence);

    free_consensus_data(consensus);
    for (int i = 0; i <seq_count + 1; i++)
        free_align_tags(tags_list[i]);
	
    la.CloseDB(); //close database*/
    return 0;
}