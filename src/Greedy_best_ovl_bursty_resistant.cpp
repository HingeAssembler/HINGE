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
#include <set>


#define LAST_READ_SYMBOL  '$'

static int ORDER(const void *l, const void *r) {
    int x = *((int32 *) l);
    int y = *((int32 *) r);
    return (x - y);
}


bool compare_overlap(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->aepos - ovl1->abpos + ovl1->bepos - ovl1->bbpos) > (ovl2->aepos - ovl2->abpos + ovl2->bepos - ovl2->bbpos));
}

bool compare_sum_overlaps(const std::vector<LOverlap * > & ovl1, const std::vector<LOverlap *> & ovl2) {
    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < ovl1.size(); i++) sum1 += ovl1[i]->aepos - ovl1[i]->abpos + ovl1[i]->bepos - ovl1[i]->bbpos;
    for (int i = 0; i < ovl2.size(); i++) sum2 += ovl2[i]->aepos - ovl2[i]->abpos + ovl2[i]->bepos - ovl2[i]->bbpos;
    return sum1 > sum2;
}

bool compare_pos(LOverlap * ovl1, LOverlap * ovl2) {
    return (ovl1->abpos) > (ovl2->abpos);
}

bool compare_overlap_abpos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->abpos < ovl2->abpos;
}

std::vector<std::pair<int,int>> Merge(std::vector<LOverlap *> & intervals)
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

    int left=intervals[0]->abpos, right = intervals[0]->aepos; //left, right means maximal possible interval now

    for(int i = 1; i < n; i++)
    {
        if(intervals[i]->abpos <= right)
        {
            right=std::max(right,intervals[i]->aepos);
        }
        else
        {
            ret.push_back(std::pair<int, int>(left,right));
            left = intervals[i]->abpos;
            right = intervals[i]->aepos;
        }
    }
    ret.push_back(std::pair<int, int>(left,right));
    return ret;
}



int main(int argc, char *argv[]) {
    LAInterface la;
	char * name_db = argv[1];
	char * name_las = argv[2];
	printf("name of db: %s, name of .las file %s\n", name_db, name_las);
    la.OpenDB(name_db);
	//la.OpenDB2(name_db); // changed this on Oct 12, may case problem, xf1280@gmail.com
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

	/**
	 * Remove reads shorter than length threshold
	 */
	int LENGTH_THRESHOLD = 12000;
	for (int i = 0; i < aln.size(); i++) {
		if ((aln[i]->alen < LENGTH_THRESHOLD) or (aln[i]->blen < LENGTH_THRESHOLD)) aln[i]->active = false;
	}





/* try and see the numbers of different types
	
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
*/

	std::map<std::pair<int,int>, std::vector<LOverlap *> > idx; //map from (aid, bid) to alignments in a vector
	std::map<int, std::vector<std::vector<LOverlap*>> > idx2; //map from (aid) to alignments in a vector
    std::vector< std::pair<Node, Node> > edgelist; // save output to edgelist
    std::map<int, std::vector <LOverlap * > >idx3;

    std::map< int, std::set<int> > has_overlap;

    for (int i = 0; i< n_read; i++) {
        has_overlap[i] = std::set<int>();
        idx3[i] = std::vector<LOverlap *>();
    }
    for (int i = 0; i < aln.size(); i++) {
        //printf("%d,%d\n",aln[i]->aid, aln[i]->bid);
        if (aln[i]->diffs / float(aln[i]->bepos - aln[i]->bbpos + aln[i]->aepos - aln[i]->abpos) < 0.23 ) {
            idx[std::pair<int,int>(aln[i]->aid, aln[i]->bid )] = std::vector< LOverlap *> ();
        }
    }

    for (int i = 0; i < aln.size(); i++) {
    if (aln[i]->diffs / float(aln[i]->bepos - aln[i]->bbpos + aln[i]->aepos - aln[i]->abpos) < 0.23 ) {
		idx[std::pair<int,int>(aln[i]->aid, aln[i]->bid )].push_back(aln[i]);
        has_overlap[aln[i]->aid].insert(aln[i]->bid);
	    idx3[aln[i]->aid].push_back(aln[i]);
        }
    }

    //sort each a,b pair according to abpos:
    for (int i = 0; i < n_read; i++)
        for (std::set<int>::iterator j = has_overlap[i].begin(); j != has_overlap[i].end(); j++) {
            std::sort(idx[std::pair<int,int>(i, *j)].begin(), idx[std::pair<int,int>(i, *j)].end(), compare_pos);
        }

    for (int i = 0; i < n_read; i++)
        for (std::set<int>::iterator j = has_overlap[i].begin(); j != has_overlap[i].end(); j++) {
            idx2[i].push_back(idx[std::pair<int,int>(i, *j)]);
        }

    //sort the reads
    for (int i = 0; i < n_read; i++ ) {
        std::sort( idx2[i].begin(), idx2[i].end(), compare_sum_overlaps );
    }

    std::cout<<"here" << std::endl;

    std::map<int,std::vector<std::pair<int,int> > > covered_region;
    for (int i = 0; i < n_read; i++ ) {
        covered_region[i] = Merge(idx3[i]);
    } // find all covered regions, could help remove adaptors

    std::cout<<"here2" << std::endl;

    for (int i = 0; i < n_read; i++) {
        printf("\n read %d:", i);
        for (int j = 0; j < covered_region[i].size(); j++) {
            printf("[%d, %d] ",covered_region[i][j].first, covered_region[i][j].second);
        }
    }


    /*
    for (int i = 0; i < n_read; i++ ) {
        printf("read:%d\n",i);
        for (int j = 0; j<idx2[i].size(); j++) {
            int sum = 0;
            for (int k = 0; k < idx2[i][j].size(); k++) {
                sum +=  idx2[i][j][k]->aepos + idx2[i][j][k]->bepos - idx2[i][j][k]->abpos - idx2[i][j][k]->bbpos;
            }

            sum /= 2;
            printf("sum:%d\n",sum);
        }
    }*/

    /*
     * Debug output
     */
    /*for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < idx2[i].size(); j++) {
            printf("%d,%d,%d\n",i, idx2[i].size(),idx2[i][j].size() );

            for (int k = 0; k<idx2[i][j].size(); k++) {
                printf(" -%d\n",idx2[i][j][k]->abpos);
            }
        }
    }*/


    //from the list, find the first one that can extend read A to the right, this will form a graph


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
			if (idx2[i][j][0]->active) {
        	    if ((idx2[i][j].front()->aln_type == FORWARD) and (cf < 1)) {
        	        cf += 1;
        	        //add edge
        	        if (idx2[i][j][0]->flags == 1) { // n = 0, c = 1
        	            edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j][0]->aid,0),Node(idx2[i][j][0]->bid,1)));
        	        } else {
        	            edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j][0]->aid,0),Node(idx2[i][j][0]->bid,0)));
        	        }
        	    }

        	    if ((idx2[i][j].back()->aln_type == BACKWARD) and (cb < 1)) {
        	        cb += 1;
        	        //add edge
        	        if (idx2[i][j][0]->flags == 1) {
        	            edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j][0]->aid,1),Node(idx2[i][j][0]->bid,0)));
        	        } else {
        	            edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j][0]->aid,1),Node(idx2[i][j][0]->bid,1)));
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
            if (idx2[i][j][0]->active) {
				if ((idx2[i][j].front()->aln_type == MISMATCH_RIGHT) and (cf < 1)) {
            	    cf += 1;
            	    //add edge
            	    if (idx2[i][j][0]->flags == 1) { // n = 0, c = 1
            	        edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j][0]->aid,0),Node(idx2[i][j][0]->bid,1)));
            	    } else {
            	        edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j][0]->aid,0),Node(idx2[i][j][0]->bid,0)));
            	    }
            	}

            	if ((idx2[i][j].back()->aln_type == MISMATCH_LEFT) and (cb < 1)) {
            	    cb += 1;
            	    //add edge
            	    if (idx2[i][j][0]->flags == 1) {
            	        edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j][0]->aid,1),Node(idx2[i][j][0]->bid,0)));
            	    } else {
            	        edgelist.push_back(std::pair<Node, Node> (Node(idx2[i][j][0]->aid,1),Node(idx2[i][j][0]->bid,1)));
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



