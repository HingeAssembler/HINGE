//
// Created by Fei Xia on 10/2/15.
//

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#ifndef PROJECT_OVERLAPGRAPH_H
#define PROJECT_OVERLAPGRAPH_H

class OverlapNode {
    public:
        int nodeid;
        std::vector<int> edgelist;
        int next; // next node id, temporary, for greedy algorithm
        OverlapNode() {};
};

class Node {
public:
    int id;
    int strand;
    bool pseudo = false;
    Node(int id, int strand): id(id), strand(strand) {};
    Node() {};
    void show(std::ofstream & out ) { out<<id; if (strand == 1) out<<"\'";}
    void show() { std::cout<<id; if (strand == 1) std::cout<<"\'";}
};



/*
class OverlapEdge {

};
*/

class OverlapGraph {
    public:
        std::vector< OverlapNode * > nodes;
        std::map<int, OverlapNode *> node_map;
        OverlapGraph() {} ;
        void Graph2FastG() {};

};

#endif //PROJECT_OVERLAPGRAPH_H
