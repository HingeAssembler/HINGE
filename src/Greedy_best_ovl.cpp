#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <map>
#include "DB.h"
#include "align.h"
#include "LAInterface.h"
#include "OverlapGraph.h"
#include <algorithm>
#include <fstream>

#include <iostream>


#define LAST_READ_SYMBOL  '$'

static int ORDER(const void *l, const void *r) {
    int x = *((int32 *) l);
    int y = *((int32 *) r);
    return (x - y);
}


bool compare_overlap(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->aepos - ovl1->abpos + ovl1->bepos - ovl1->bbpos) > (ovl2->aepos - ovl2->abpos + ovl2->bepos - ovl2->bbpos));

}



int main(int argc, char *argv[]) {
    LAInterface la;
	char * name_db = argv[1];
	char * name_las = argv[2];
	printf("name of db: %s, name of .las file %s\n", name_db, name_las);
    la.OpenDB(name_db);
    std::cout<<"# Reads:" << la.getReadNumber() << std::endl;

    la.OpenAlignment(name_las);
    std::cout<<"# Alignments:" << la.getAlignmentNumber() << std::endl;

	int n_aln = la.getAlignmentNumber();
	int n_read = la.getReadNumber();
    std::vector<LOverlap *> aln;
    la.getOverlap(aln,0,n_aln);

	std::map<std::pair<int,int>, int> idx; //map from (aid, bid) to alignment id
	std::map<int, std::vector<LOverlap*>> idx2; //map from (aid) to alignment id in a vector

    std::vector< std::pair<Node, Node> > edgelist;

	for (int i = 0; i < n_read; i++ ) {
		idx2[i] = std::vector< LOverlap * >(); // initialize idx2
	}

	for (int i = 0; i < aln.size(); i++) {
		idx[std::pair<int,int>(aln[i]->aid, aln[i]->bid )] = i;
		idx2[aln[i]->aid].push_back(aln[i]);
	}

    //sort the reads

    for (int i = 0; i < n_read; i++ ) {
        std::sort( idx2[i].begin(), idx2[i].end(), compare_overlap );
    }

    /****
	get the best overlap for each read and form a graph
	****/

    for (int i = 0; i < idx2.size(); i++) {
        int cf = 0;
        int cb = 0;
        for (int j = 0; j< idx2[i].size(); j++) {
            //idx2[i][j]->show();
            if ((idx2[i][j]->aln_type == FORWARD) and (cf == 0)) {
                cf = 1;
                //add edge
                if (idx2[i][j]->flags == 1) { // n = 0, c = 1
                    edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j]->aid,0),Node(idx2[i][j]->bid,1)));
                } else {
                    edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j]->aid,0),Node(idx2[i][j]->bid,0)));
                }
            }

            if ((idx2[i][j]->aln_type == BACKWARD) and (cb == 0)) {
                cb = 1;
                //add edge
                if (idx2[i][j]->flags == 1) {
                    edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j]->aid,1),Node(idx2[i][j]->bid,0)));
                } else {
                    edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j]->aid,1),Node(idx2[i][j]->bid,1)));
                }
            }
            if ((cf == 1) and (cb == 1)) break;
        }
    }

    std::ofstream out  (argv[3], std::ofstream::out);


    //print edges list
    for (int i = 0; i < edgelist.size(); i++){
        edgelist[i].first.show(out);
        out<<"->";
        edgelist[i].second.show(out);
        out<<std::endl;
    }



	//for (int i=0; i <idx2[1500].size(); i++) {
	//	printf("%d\n",idx2[1500][i]);
	//}

	/*std::cout << "hello" << std::endl;
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
    la.getAlignment(res2, 0, 3);// get alignment for reads [0,3)

	for (auto i:res2)
		i->show();
	printf("\n");
    */
    la.CloseDB(); //close database
    return 0;
}



