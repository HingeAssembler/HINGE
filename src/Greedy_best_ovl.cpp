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
	la.OpenDB2(name_db);
    std::cout<<"# Reads:" << la.getReadNumber() << std::endl;

    la.OpenAlignment(name_las);
    std::cout<<"# Alignments:" << la.getAlignmentNumber() << std::endl;
	//la.resetAlignment();
	//la.showOverlap(0,1);

	int n_aln = la.getAlignmentNumber();
	int n_read = la.getReadNumber();
    std::vector<LOverlap *> aln;
	la.resetAlignment();
    la.getOverlap(aln,0,n_aln);

	int covered = 0, forward = 0, backward = 0, covering = 0, undefined = 0;
	int pcovered = 0, pforward = 0, pbackward = 0, pcovering = 0, pundefined = 0;
	
	/**
	filter reads
	**/
	
	int LENGTH_THRESHOLD = 12000;
	for (int i = 0; i < aln.size(); i++) {
		if ((aln[i]->alen < LENGTH_THRESHOLD) or (aln[i]->blen < LENGTH_THRESHOLD)) aln[i]->active = false;
	}
	
	
	
	for (int i = 0; i < aln.size(); i++) {
		int CHI_THRESHOLD = 300;
		
		    if ((aln[i]->abpos > 0) and (aln[i]->aepos > aln[i]->alen - CHI_THRESHOLD) ) {
		        forward ++; 
		    }

		    else if ( ( aln[i]->abpos < CHI_THRESHOLD) and (aln[i]->aepos < aln[i]->alen)) {
		        backward ++ ;
		    }
			
			else if ((aln[i]->abpos < CHI_THRESHOLD) and (aln[i]->aepos > aln[i]->alen - CHI_THRESHOLD) ) {
			        if (aln[i]->blen > aln[i]->alen) covered ++ ;
			}
		    
			else if ((aln[i]->bbpos < CHI_THRESHOLD) and (aln[i]->bepos > aln[i]->blen - CHI_THRESHOLD) ) {
		        if (aln[i]->alen >aln[i]-> blen) covering ++; 
		    }
			else undefined ++;
		
		
		if ((aln[i]->abpos > 0) and (aln[i]->aepos == aln[i]->alen) ) {
		        pforward ++;
		    }

		    else if 
			( ( aln[i]->abpos == 0) and (aln[i]->aepos < aln[i]->alen)) {
					        pbackward ++; 
			}

		    else if ((aln[i]->abpos  == 0) and (aln[i]->aepos == aln[i]->alen) ) {
		        if (aln[i]->blen > aln[i]->alen) pcovered ++; 
		    }
			
		    else if ((aln[i]->bbpos == 0) and (aln[i]->bepos == aln[i]->blen) ) {
		        if (aln[i]->alen > aln[i]->blen) pcovering ++;
		    }	
			else pundefined ++;
	}
	
	printf("covered %d forward %d backward %d covering %d undefined %d\ncovered %d forward %d backward %d covering %d undefined %d\n",covered, forward, backward, covering, undefined, pcovered, pforward, pbackward, pcovering, pundefined);


	std::map<std::pair<int,int>, int> idx; //map from (aid, bid) to alignment id
	std::map<int, std::vector<LOverlap*>> idx2; //map from (aid) to alignment id in a vector

    std::vector< std::pair<Node, Node> > edgelist;

	for (int i = 0; i < n_read; i++ ) {
		idx2[i] = std::vector< LOverlap * >(); // initialize idx2
	}

	for (int i = 0; i < aln.size(); i++) {
    if (aln[i]->diffs / float(aln[i]->bepos - aln[i]->bbpos + aln[i]->aepos - aln[i]->abpos) < 0.22 ) {

		idx[std::pair<int,int>(aln[i]->aid, aln[i]->bid )] = i;
		idx2[aln[i]->aid].push_back(aln[i]);

	    }
    }

	for (int j = 0; j< idx2[1].size(); j++)
		idx2[1][j]->show();

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

		/*
		 * For each read, if there is exact right match (type FORWARD), choose the one with longest alignment with read A
		 * same for BACKWARD,
		 */
        for (int j = 0; j< idx2[i].size(); j++) {
        	    //idx2[i][j]->show();
			if (idx2[i][j]->active) {
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
			}
        	if ((cf == 1) and (cb == 1)) break;
        	
		}
		/*
		 * For each read, if there is no exact right or left match, choose that one with a chimeric end, but still choose
		 * the one with longest alignment, same for BACKWARD
		 */

        for (int j = 0; j< idx2[i].size(); j++) {
            //idx2[i][j]->show();
            if (idx2[i][j]->active) {
				if ((idx2[i][j]->aln_type == MISMATCH_RIGHT) and (cf == 0)) {
            	    cf = 1;
            	    //add edge
            	    if (idx2[i][j]->flags == 1) { // n = 0, c = 1
            	        edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j]->aid,0),Node(idx2[i][j]->bid,1)));
            	    } else {
            	        edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j]->aid,0),Node(idx2[i][j]->bid,0)));
            	    }
            	}
            	
            	if ((idx2[i][j]->aln_type == MISMATCH_LEFT) and (cb == 0)) {
            	    cb = 1;
            	    //add edge
            	    if (idx2[i][j]->flags == 1) {
            	        edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j]->aid,1),Node(idx2[i][j]->bid,0)));
            	    } else {
            	        edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j]->aid,1),Node(idx2[i][j]->bid,1)));
            	    }
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


    la.CloseDB(); //close database
    return 0;
}



