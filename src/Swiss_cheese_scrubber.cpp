//
// Created by Fei Xia on 10/8/15.
//

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
#include <sstream>
#include <limits>
#include <iostream>

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
    return ((ovl1->aepos - ovl1->abpos + ovl1->bepos - ovl1->bbpos) < (ovl2->aepos - ovl2->abpos + ovl2->bepos - ovl2->bbpos));
}


//merge intervals
std::vector<std::pair<int,int>> Merge(std::vector<LOverlap *> & intervals)
{
    std::vector<std::pair<int, int > > ret;
    const int n = intervals.size();
    if(n == 1) {
        ret.push_back(std::pair<int,int>(intervals[0]->abpos,intervals[0]->aepos));
    }

    sort(intervals.begin(),intervals.end(),compare_overlap); //首先按照区间的左端点进行排序

    int left=intervals[0]->abpos, right = intervals[0]->aepos; //left, right表示当前区间的可能的最大范围

    for(int i=1;i<n;++i)
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




int main(int argc, char ** argv) {

    std::cout<<"hello world"<<std::endl;

    LAInterface la;
    char * name_db = argv[1];
    char * name_las = argv[2];
    printf("name of db: %s, name of .las file %s\n", name_db, name_las);
    la.OpenDB(name_db);
    //la.OpenDB2(name_db);
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

    std::map<std::pair<int,int>, std::vector<LOverlap *> > idx; //map from (aid, bid) to alignment id
    std::map<int, std::vector<LOverlap*>> idx2; //map from (aid) to alignment id in a vector

    for (int i = 0; i < n_read; i++ ) {
        idx2[i] = std::vector< LOverlap * >(); // initialize idx2
    }

    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->diffs / float(aln[i]->bepos - aln[i]->bbpos + aln[i]->aepos - aln[i]->abpos) < 0.5 ) {
            idx2[aln[i]->aid].push_back(aln[i]);
        }
    }
	
    for (int i = 0; i < n_read; i++ ) {
        std::sort(idx2[i].begin(), idx2[i].end(), compare_overlap);
    }
	
	for (int i = 0; i < n_read; i++) 
		for (int j = 0; j<idx2[i].size(); j++) {
            if (idx.find(std::pair<int,int>(idx2[i][j]->aid, idx2[i][j]->bid )) == idx.end())
                idx[std::pair<int,int>(idx2[i][j]->aid, idx2[i][j]->bid )] = std::vector<LOverlap *> ();

            idx[std::pair<int,int>(idx2[i][j]->aid, idx2[i][j]->bid )].push_back(idx2[i][j]);
		}

    std::vector<Read *> reads;
    la.getRead(reads,0,n_read);


    for (int i = 0; i < n_read; i++ ) {
        std::vector<std::pair<int,int> > covered_region;
        covered_region = Merge(idx2[i]);
        
    }


    /*
    std::vector< std::pair<Node, Node> > edgelist;
    char * name_edges = argv[3];

    std::string edge_line;
    std::ifstream edges_file(name_edges);

    while (!edges_file.eof()) {
        std::getline(edges_file,edge_line);
        std::vector<std::string> tokens = split(edge_line, ' ');
        Node node0;
        Node node1;
        if (tokens.size() == 2) {
            if (tokens[0].back() == '\'') {
                node0.id = std::stoi(tokens[0].substr(0, tokens[0].size() - 1));
                node0.strand = 1;
            } else {
                node0.id = std::stoi(tokens[0]);
                node0.strand = 0;
            }

            if (tokens[1].back() == '\'') {
                node1.id = std::stoi(tokens[1].substr(0, tokens[1].size() - 1));
                node1.strand = 1;
            } else {
                node1.id = std::stoi(tokens[1]);
                node1.strand = 0;
            }
            edgelist.push_back(std::pair<Node,Node>(node0,node1));
        }
    }*/


    la.CloseDB(); //close database
    return 0;
}
