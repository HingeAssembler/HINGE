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
	
	la.getAlignment(res2, 0, 2);


    for (int i = 0; i < 5; i ++) {
        res2[i]->show();
        /*int tlen = res2[i]->tlen;
        int *trace = (int *) res2[i]->trace;
        int u;
        printf(" ");

        for (u = 0; u < tlen; u++) {
            printf("%d,", (int) trace[u]);
        }
        printf("\n");
        la.Lshow_Alignment_tgs(res2[i]);*/

        /*printf("\n");
        for (int j = 0; j < res2[i]->trace_pts_len; j++)
            printf(" %d", res2[i]->trace_pts[j]);
        printf("\n");*/

    }
    la.CloseDB(); //close database
    return 0;
}