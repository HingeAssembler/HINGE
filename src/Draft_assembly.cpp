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


std::string reverse_complement(std::string seq) {
    static std::map<char, char> m = {{'a','t'}, {'c','g'}, {'g','c'}, {'t','a'}, {'A','T'}, {'C','G'}, {'T','A'}, {'G','C'}, {'n','n'}, {'N', 'N'}};
    std::reverse(seq.begin(), seq.end());
    for (int i = 0; i < seq.size(); i++) {
        seq[i] = m[seq[i]];
    }
    return seq;
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

    std::map<std::pair<int,int>, LOverlap *> idx; //map from (aid, bid) to alignment id
    std::map<int, std::vector<LOverlap*>> idx2; //map from (aid) to alignment id in a vector

    for (int i = 0; i < n_read; i++ ) {
        idx2[i] = std::vector< LOverlap * >(); // initialize idx2
    }

    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->diffs / float(aln[i]->bepos - aln[i]->bbpos + aln[i]->aepos - aln[i]->abpos) < 0.5 ) {
            idx[std::pair<int,int>(aln[i]->aid, aln[i]->bid )] = aln[i];
            idx2[aln[i]->aid].push_back(aln[i]);
        }
    }

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

    std::string sequence = "";
    Node lastnode;
    LOverlap * lastoverlap = NULL;

    for (int i = 0; i < edgelist.size(); i++){
        //edgelist[i].first.show();
        //std::cout<<"->";
        //edgelist[i].second.show();
        //std::cout<<std::endl;
        if (edgelist[i].first.id != lastnode.id) {
            //std::cout<<std::endl;
            //edgelist[i].first.show();
            //std::cout<<"->";
            /**
             * start a new "sequence"
             */
            if (edgelist[i].first.strand == 0)
                sequence.append(reads[edgelist[i].first.id]->bases);
            else
                sequence.append(reverse_complement(reads[edgelist[i].first.id]->bases));
        }
        //edgelist[i].second.show();
        //std::cout<<"->";
        LOverlap * currentaln = NULL;
        lastnode = edgelist[i].second;
        if (idx.find(std::pair<int,int>(edgelist[i].first.id, edgelist[i].second.id)) != idx.end()) {
            currentaln = idx[std::pair<int,int>(edgelist[i].first.id, edgelist[i].second.id)];
        }

        if (currentaln == NULL) std::cout<<"Warning, cannot find alignment!\n";
        /***
         * play with current alignment
         */
        std::string current_seq;

        if (edgelist[i].second.strand == 0)
            current_seq = reads[edgelist[i].second.id]->bases;
        else
            current_seq = reverse_complement(reads[edgelist[i].second.id]->bases);


        /*
        std::cout<< "seq " << (currentaln->alen - currentaln->aepos) << std::endl;
        sequence.erase(sequence.end() - (currentaln->alen - currentaln->aepos), sequence.end());


        current_seq = current_seq.substr(currentaln->bepos);

        std::cout<<"next " << current_seq.size() << " " << (currentaln->blen - currentaln->bepos) << std::endl;

        */
        //if ((currentaln->bbpos > 0) and (current_seq.size() > currentaln->bbpos))
        //    current_seq.erase(current_seq.begin(), current_seq.begin() + currentaln->bbpos);

        if (edgelist[i].first.strand == 0) {
            sequence.erase(sequence.end() - (currentaln->alen - currentaln->aepos), sequence.end());
            current_seq.erase(current_seq.begin(), current_seq.begin() + currentaln->bepos);
            sequence.append(current_seq);
        }
        else {
            sequence.erase(sequence.end() - (currentaln->abpos), sequence.end());
            current_seq.erase(current_seq.begin(), current_seq.begin() + currentaln->blen - currentaln->bbpos);
            sequence.append(current_seq);
        }
        //std::cout<<current_seq<<std::endl;
        lastoverlap = currentaln;
    }


    //std::cout<<sequence<<std::endl;
    std::cout<<sequence.size()<<std::endl;

    std::ofstream out(argv[4]);
    out << "> draft assembly" <<std::endl;
    out << sequence;
    la.CloseDB(); //close database
    return 0;
}