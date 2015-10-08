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


int main(int argc, char ** argv) {

    std::cout<<"hello world"<<std::endl;

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


    std::vector<Read *> reads;
    la.getRead(reads,0,n_read);

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
    }

    /**
     * Use DFS edges to generate assembly
     */

    std::string sequence;
    for (int i = 0; i < edgelist.size(); i++){
        edgelist[i].first.show();
        std::cout<<"->";
        edgelist[i].second.show();
        std::cout<<std::endl;

    }




    la.CloseDB(); //close database
    return 0;
}