#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <unordered_map>
#include "DB.h"
#include "align.h"
#include "LAInterface.h"
#include "OverlapGraph.h"
#include <algorithm>
#include <fstream>

#include <iostream>
#include <set>
#include <omp.h>
#include "INIReader.h"
#include <tuple>


#define LAST_READ_SYMBOL  '$'


typedef std::tuple<Node, Node, int> Edge_w;

typedef std::pair<Node, Node> Edge_nw;


static int ORDER(const void *l, const void *r) {
    int x = *((int32 *) l);
    int y = *((int32 *) r);
    return (x - y);
}


std::ostream& operator<<(std::ostream& out, const aligntype value){
    static std::map<aligntype, std::string> strings;
    if (strings.size() == 0){
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(FORWARD);
        INSERT_ELEMENT(BACKWARD);
        INSERT_ELEMENT(MISMATCH_LEFT);
        INSERT_ELEMENT(MISMATCH_RIGHT);
        INSERT_ELEMENT(COVERED);
        INSERT_ELEMENT(COVERING);
        INSERT_ELEMENT(UNDIFINED);
        INSERT_ELEMENT(MIDDLE);
#undef INSERT_ELEMENT
    }

    return out << strings[value];
}



bool compare_overlap(LOverlap * ovl1, LOverlap * ovl2) {
    return ((ovl1->aepos - ovl1->abpos + ovl1->bepos - ovl1->bbpos) > (ovl2->aepos - ovl2->abpos + ovl2->bepos - ovl2->bbpos));
}

bool compare_sum_overlaps(const std::vector<LOverlap * > * ovl1, const std::vector<LOverlap *> * ovl2) {
    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0; i < ovl1->size(); i++) sum1 += (*ovl1)[i]->aepos - (*ovl1)[i]->abpos + (*ovl1)[i]->bepos - (*ovl1)[i]->bbpos;
    for (int i = 0; i < ovl2->size(); i++) sum2 += (*ovl2)[i]->aepos - (*ovl2)[i]->abpos + (*ovl2)[i]->bepos - (*ovl2)[i]->bbpos;
    return sum1 > sum2;
}

bool compare_pos(LOverlap * ovl1, LOverlap * ovl2) {
    return (ovl1->abpos) > (ovl2->abpos);
}

bool compare_overlap_abpos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->abpos < ovl2->abpos;
}

bool compare_overlap_aepos(LOverlap * ovl1, LOverlap * ovl2) {
    return ovl1->abpos > ovl2->abpos;
}

std::vector<std::pair<int,int>> Merge(std::vector<LOverlap *> & intervals, int cutoff)
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

    int left=intervals[0]->abpos + cutoff, right = intervals[0]->aepos - cutoff; //left, right means maximal possible interval now

    for(int i = 1; i < n; i++)
    {
        if(intervals[i]->abpos + cutoff <= right)
        {
            right=std::max(right,intervals[i]->aepos - cutoff);
        }
        else
        {
            ret.push_back(std::pair<int, int>(left,right));
            left = intervals[i]->abpos + cutoff;
            right = intervals[i]->aepos - cutoff;
        }
    }
    ret.push_back(std::pair<int, int>(left,right));
    return ret;
}

Interval Effective_length(std::vector<LOverlap *> & intervals, int min_cov) {
    Interval ret;
    sort(intervals.begin(),intervals.end(),compare_overlap_abpos); //sort according to left

    if (intervals.size() > min_cov) {
        ret.first = intervals[min_cov]->abpos;
    } else
        ret.first = 0;
    sort(intervals.begin(),intervals.end(),compare_overlap_aepos); //sort according to left
    if (intervals.size() > min_cov) {
        ret.second = intervals[min_cov]->aepos;
    } else
        ret.second = 0;
    return ret;
}

bool bridge(LOverlap* ovl, int s, int e){
    return ((ovl->abpos < s - 500) and (ovl->aepos > e + 500));
}


int main(int argc, char *argv[]) {

    LAInterface la;
	char * name_db = argv[1];
	char * name_las = argv[2];
    char * name_mask = argv[3];
    char * name_config = argv[4];
	printf("name of db: %s, name of .las file %s\n", name_db, name_las);
    la.openDB(name_db);
    std::cout<<"# Reads:" << la.getReadNumber() << std::endl;
    la.openAlignmentFile(name_las);
    std::cout<<"# Alignments:" << la.getAlignmentNumber() << std::endl;

	int n_aln = la.getAlignmentNumber();
	int n_read = la.getReadNumber();
    std::vector<LOverlap *> aln;
	la.resetAlignment();
    la.getOverlap(aln,0,n_aln);
    std::vector<Read *> reads;
    la.getRead(reads,0,n_read);
	std::vector<std::vector<int>>  QV;
	la.getQV(QV,0,n_read);
    std::cout << "input data finished" <<std::endl;

    for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < QV[i].size(); j++) QV[i][j] = int(QV[i][j] < 40);
    }

    std::vector<std::pair<int, int> > QV_mask;

    for (int i = 0; i < n_read; i++) {
        int s = 0, e = 0;
        int max = 0, maxs = s, maxe = e;

        for (int j = 0; j < QV[i].size(); j++) {
            if (QV[i][j] == 1) {
                e ++;
            }
            else {
                if (e - s > max) {
                    maxe = e ; maxs = s;
                    max = e - s;
                }

                s = j+1;
                e = j+1;
            }
        }

        QV_mask.push_back(std::pair<int, int>(maxs*la.tspace, maxe*la.tspace));
        // create mask by QV
    }


    /*for (int i = 0; i < 200; i++) {
        printf("read %d: ",i+1);

        printf("%d %d ", QV[i].size(), reads[i]->len);

        for (int j = 0; j < QV[i].size(); j++) printf("%d ", QV[i][j]);
        printf("[%d %d]", QV_mask[i].first, QV_mask[i].second);
        printf("\n");
    }
    //for debug
    */


    INIReader reader(name_config);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<<name_config<<std::endl;
        return 1;
    }

    int LENGTH_THRESHOLD = reader.GetInteger("filter", "length_threshold", -1);
    double QUALITY_THRESHOLD = reader.GetReal("filter", "quality_threshold", 0.0);
    int N_ITER = reader.GetInteger("filter", "n_iter", -1);
    int ALN_THRESHOLD = reader.GetInteger("filter", "aln_threshold", -1);
    int MIN_COV = reader.GetInteger("filter", "min_cov", -1);
    int CUT_OFF = reader.GetInteger("filter", "cut_off", -1);
    int THETA = reader.GetInteger("filter", "theta", -1);
	int N_PROC = reader.GetInteger("running", "n_proc", 4);
    int EST_COV = reader.GetInteger("filter", "ec", 0);
    int reso = 40;

    omp_set_num_threads(N_PROC);

    std::vector<std::vector <LOverlap * > >idx3; // this is the pileup

    for (int i = 0; i< n_read; i++) {
        idx3.push_back(std::vector<LOverlap *>());
    }

    for (int i = 0; i < aln.size(); i++) {
        if (aln[i]->active) {
            idx3[aln[i]->aid].push_back(aln[i]);
        }
    }

# pragma omp parallel for
    for (int i = 0; i < n_read; i++) {
        std::sort(idx3[i].begin(), idx3[i].end(), compare_overlap);
    }

    std::cout<<"profile coverage" << std::endl;
    std::ofstream cov(std::string(name_db) + ".coverage.txt");
    std::ofstream homo(std::string(name_db) + ".homologous.txt");
    std::ofstream rep(std::string(name_db) + ".repeat.txt");
    std::ofstream filtered(std::string(name_db) + ".filtered.fasta");
    std::ofstream hg(std::string(name_db) + ".hinges.txt");



    std::vector< std::vector<std::pair<int, int> > > coverages;
    std::vector< std::vector<std::pair<int, int> > > cgs; //coverage gradient;

    for (int i = 0; i < n_read; i ++) {
        std::vector<std::pair<int, int> > coverage;
        std::vector<std::pair<int, int> > cg;
        la.profileCoverage(idx3[i], coverage, reso, CUT_OFF);
        cov << "read " << i <<" ";
        for (int j = 0; j < coverage.size(); j++)
            cov << coverage[j].first << ","  << coverage[j].second << " ";
        cov << std::endl;

        if (coverage.size() >= 2)
            for (int j = 0; j < coverage.size() - 1; j++) {
                cg.push_back(std::pair<int,int>(coverage[j].first, coverage[j+1].second - coverage[j].second));
            }
        else cg.push_back(std::pair<int, int> (0,0));

        coverages.push_back(coverage);
        cgs.push_back(cg);
    }

    int num_slot = 0;
    int total_cov = 0;
    for (int i = 0; i < n_read/500; i++) {
        for (int j = 0; j < coverages[i].size(); j++) {
            printf("%d\n", coverages[i][j].second);
            total_cov += coverages[i][j].second;
            num_slot ++;
        }
    }

    std::cout << "Estimated coverage:" << total_cov / float(num_slot) << std::endl;
    int cov_est = total_cov / num_slot;


    if (EST_COV != 0) cov_est = EST_COV;
    std::cout << "Estimated coverage:" << cov_est << std::endl;

    /*coverages.clear();
    for (int i = 0; i < n_read; i ++) {
        std::vector<std::pair<int, int> > coverage;
        la.profileCoveragefine(idx3[i], coverage, reso, CUT_OFF, total_cov/num_slot);
        cov << "read " << i <<" ";
        for (int j = 0; j < coverage.size(); j++)
            cov << coverage[j].first << ","  << coverage[j].second << " ";
        cov << std::endl;
        coverages.push_back(coverage);
    }*/

    std::vector<std::pair<int, int>> maskvec;
    if (MIN_COV < cov_est/3)
        MIN_COV = cov_est/3;

    std::ofstream mask(name_mask);
    for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < coverages[i].size(); j++) {
            coverages[i][j].second -= MIN_COV;
            if (coverages[i][j].second < 0) coverages[i][j].second = 0;
        }
        int start = 0;
        int end = start;
        int maxlen = 0, maxstart = 0, maxend = 0;
        for (int j = 0; j < coverages[i].size(); j++) {
            if (coverages[i][j].second > 0) {
                end = coverages[i][j].first;
            } else {
                if (end > start) {
                    //std::cout<<"read" << i << " "<<start+reso << "->" << end << std::endl;
                    if (end - start - reso > maxlen) {
                        maxlen = end - start - reso;
                        maxstart = start + reso;
                        maxend = end;
                    }
                }
                start = coverages[i][j].first;
                end = start;
            }
        }
        //mask << i << " " << maxstart << " " << maxend << std::endl;
        mask << i << " " << std::max(maxstart, QV_mask[i].first) << " " << std::min(maxend, QV_mask[i].second) << std::endl;
        int s = std::max(maxstart, QV_mask[i].first);
        int l = std::min(maxend, QV_mask[i].second) - std::max(maxstart, QV_mask[i].first);
        if (l < 0) l = 0;
        filtered << ">read_" << i << std::endl;
        filtered << reads[i]->bases.substr(s,l) << std::endl;

        //maskvec.push_back(std::pair<int, int>(maxstart + 200, maxend - 200));
        maskvec.push_back(std::pair<int, int>(std::max(maxstart + 200, QV_mask[i].first), std::min(maxend - 200, QV_mask[i].second)));
    }







    //binarize coverage gradient;

    std::vector<std::vector<std::pair<int, int> > > repeat_anno;

    for (int i = 0; i < n_read; i++) {
        std::vector<std::pair<int, int> > anno;
        for (int j = 0; j < cgs[i].size(); j++) {
            //std::cout<< i << " " << cgs[i][j].first << " " << cgs[i][j].second << std::endl;
            if ((cgs[i][j].first >= maskvec[i].first) and (cgs[i][j].first <= maskvec[i].second)) {
                if (cgs[i][j].second > cov_est / 2) anno.push_back(std::pair<int, int>(cgs[i][j].first, 1));
                else if (cgs[i][j].second < -cov_est / 2) anno.push_back(std::pair<int, int>(cgs[i][j].first, -1));
            }
        }
        repeat_anno.push_back(anno);
    }

    int gap_thre = 300;

    // clean it a bit
    for (int i = 0; i < n_read; i++) {
        for (std::vector<std::pair<int, int> >::iterator iter = repeat_anno[i].begin(); iter < repeat_anno[i].end(); ) {
            if (iter+1 < repeat_anno[i].end()){
                if ((iter->second == 1) and ((iter+1)->second == -1) and ((iter+1)->first - iter->first < gap_thre)){
                    iter = repeat_anno[i].erase(iter);
                    iter = repeat_anno[i].erase(iter); // fill gaps
                } else if (((iter->second == 1) and ((iter + 1)->second == 1)) and ((iter+1)->first - iter->first < gap_thre)) {
                    repeat_anno[i].erase((iter + 1));
                } else if (((iter->second == -1) and ((iter + 1)->second == -1)) and ((iter+1)->first - iter->first < gap_thre)) {
                    iter = repeat_anno[i].erase(iter);
                } else iter++;
            } else iter ++;
        }
    }

    /*for (int i = 0; i < n_read; i++) {
        for (int j = 0; j < repeat_anno[i].size(); j++){
            std::cout<< i << " " << repeat_anno[i][j].first << " " << repeat_anno[i][j].second << std::endl;
        }
    }*/

    for (int i = 0; i < n_read; i++) {
        rep << i << " ";
        if (repeat_anno[i].size() > 0)
            if (repeat_anno[i].front().second == -1)
                rep << -1 << " "<<repeat_anno[i].front().first<<" ";
        bool active = true;
        for (int j = 0; j < repeat_anno[i].size(); j++) {
            if (j+1<repeat_anno[i].size())
            if ((repeat_anno[i][j].second == 1) and (repeat_anno[i][j+1].second == -1)) {
                rep << repeat_anno[i][j].first << " " << repeat_anno[i][j + 1].first << " ";
                //check if all internal repeats are bridged by pile-B reads.
                int s = repeat_anno[i][j].first;
                int e = repeat_anno[i][j+1].first;

                bool bridged = false;

                //test bridging

                for (int k = 0; k < idx3[i].size(); k++) {
                    if (bridge(idx3[i][k], s, e)) {
                        bridged = true;
                        break;
                    }
                }

                if (not bridged) active = false;
            }


        }
        if (not active) {
            //printf("read %d comes from homologus recombination\n", i);
            homo << i << std::endl;
        }

        if (repeat_anno[i].size() > 0)
        if (repeat_anno[i].back().second == 1)
            rep << repeat_anno[i].back().first<<" " << -1 << " ";

        rep << std::endl;
    }
    // need a better hinge detection

    std::unordered_map<int, std::vector<std::pair<int, int>> > hinges;
    // n_read pos -1 = in_hinge 1 = out_hinge

    for (int i = 0; i < n_read; i++) {
        hinges[i] = std::vector<std::pair<int, int>>();

        for (int j = 0; j < repeat_anno[i].size(); j++) {
            if (repeat_anno[i][j].second == -1) { // look for in hinges, negative gradient
                bool bridged = true;
                int support = 0;
                for (int k = 0; k < idx3[i].size(); k++) {
                    //printf("%d %d %d %d\n", idx3[i][k]->abpos, idx3[i][k]->aepos, maskvec[i].first, repeat_anno[i][j].first);
                    if ((idx3[i][k]->abpos < maskvec[i].first + 200) and (idx3[i][k]->aepos > repeat_anno[i][j].first - 300) and (idx3[i][k]->aepos < repeat_anno[i][j].first + 300)) {
                        support ++;
                        if (support > 7) {
                            bridged = false;
                            break;
                        }
                    }
                }
                if (not bridged) hinges[i].push_back(std::pair<int, int>(repeat_anno[i][j].first,-1));

            } else { // look for out_hinges, positive gradient
                bool bridged = true;
                int support = 0;
                for (int k = 0; k < idx3[i].size(); k++) {
                    if ((idx3[i][k]->abpos > repeat_anno[i][j].first - 300) and (idx3[i][k]->abpos < repeat_anno[i][j].first + 300) and (idx3[i][k]->aepos > maskvec[i].second - 200)) {
                        support ++;
                        if (support > 7) {
                            bridged = false;
                            break;
                        }
                    }
                }
                if (not bridged) hinges[i].push_back(std::pair<int, int>(repeat_anno[i][j].first, 1));

            }

        }

    }

    for (int i = 0; i < n_read; i++) {
        hg << i << " ";
        for (int j = 0; j < hinges[i].size(); j++) {
            hg << hinges[i][j].first << " " << hinges[i][j].second << " ";
        }
        hg << std::endl;
    }

    la.closeDB(); //close database
    return 0;
}