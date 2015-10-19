#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <sstream>
#include <iostream>
#include "LAInterface.h"
#include "align.h"
#include "DB.h"


void Read::showRead() {
    std::cout << "read #" << id << std::endl;
    std::cout << ">" << name << std::endl;
    std::cout << bases << std::endl;
}


int LAInterface::OpenDB2(std::string filename) {
    char *fn = new char[filename.length() + 1];
    strcpy(fn, filename.c_str());

    int status = Open_DB(fn, this->db1);
    if (status < 0)
        exit(1);
    if (this->db1->part > 0) {
        fprintf(stderr, "%s: Cannot be called on a block: %s\n", "test", fn);
        exit(1);
    }

    this->db2 = this->db1;
    Trim_DB(db1);

    delete[] fn;
    return 0;
}


int LAInterface::OpenDB(std::string filename) {
    char *fn = new char[filename.length() + 1];
    strcpy(fn, filename.c_str());

    int status = Open_DB(fn, this->db1);
    if (status < 0)
        exit(1);
    if (this->db1->part > 0) {
        fprintf(stderr, "%s: Cannot be called on a block: %s\n", "test", fn);
        exit(1);
    }

    this->db2 = this->db1;
    Trim_DB(db1);

    FILE *dstub;

    char *fn2 = new char[filename.length() + 1 + 3];
    strcpy(fn2, fn);
    strcat(fn2, ".db");

    dstub = Fopen(fn2, "r");
    if (dstub == NULL)
        exit(1);

    if (fscanf(dstub, DB_NFILE, &nfiles) != 1) SYSTEM_ERROR

    printf("%d files\n", nfiles);

    flist = (char **) Malloc(sizeof(char *) * nfiles, "Allocating file list");
    findx = (int *) Malloc(sizeof(int *) * (nfiles + 1), "Allocating file index");

    if (flist == NULL || findx == NULL)
        exit(1);

    findx += 1;
    findx[-1] = 0;

    for (int i = 0; i < nfiles; i++) {
        char prolog[MAX_NAME], fname[MAX_NAME];

        if (fscanf(dstub, DB_FDATA, findx + i, fname, prolog) != 3) SYSTEM_ERROR
        if ((flist[i] = Strdup(prolog, "Adding to file list")) == NULL)
            exit(1);
    }

    fclose(dstub);

    delete[] fn;


    return 0;
}

int LAInterface::CloseDB() {
    Close_DB(db1);
    return 0;
}


int LAInterface::CloseDB2() {
    Close_DB(db1);
    Close_DB(db2);
    return 0;
}

void LAInterface::showRead(int from, int to) {
    if (flist == NULL || findx == NULL)
        exit(1);
    HITS_READ *reads;
    HITS_TRACK *first;
    char *read, **entry;
    int c, b, e, i;
    int hilight, substr;
    int map;
    int       (*iscase)(int);

    read = New_Read_Buffer(db1);
    int UPPER = 1;
    int WIDTH = 80;
    //printf("2");
    {
        entry = NULL;
        first = db1->tracks;
    }


    hilight = 'A' - 'a';
    iscase = islower;

    map = 0;
    reads = db1->reads;
    substr = 0;

    c = 0;

    b = from;
    e = to;

    for (i = b; i < e; i++) {
        int len;
        int fst, lst;
        int flags, qv;
        HITS_READ *r;
        HITS_TRACK *track;

        r = reads + i;
        len = r->rlen;

        flags = r->flags;
        qv = (flags & DB_QV);

        {
            while (i < findx[map - 1])
                map -= 1;
            while (i >= findx[map])
                map += 1;
            printf(">%s/%d/%d_%d", flist[map], r->origin, r->fpulse, r->fpulse + len);
            if (qv > 0)
                printf(" RQ=0.%3d", qv);
        }
        printf("\n");


        Load_Read(db1, i, read, UPPER);

        for (track = first; track != NULL; track = track->next) {
            int64 *anno;
            int *data;
            int64 s, f, j;
            int bd, ed, m;

            anno = (int64 *) track->anno;
            data = (int *) track->data;

            s = (anno[i] >> 2);
            f = (anno[i + 1] >> 2);
            if (s < f) {
                for (j = s; j < f; j += 2) {
                    bd = data[j];
                    ed = data[j + 1];
                    for (m = bd; m < ed; m++)
                        if (iscase(read[m]))
                            read[m] = (char) (read[m] + hilight);
                    if (j == s)
                        printf("> %s:", track->name);
                    printf(" [%d,%d]", bd, ed);
                }
                printf("\n");
            }
        }


        fst = 0;
        lst = len;

        {
            int j;

            for (j = fst; j + WIDTH < lst; j += WIDTH)
                printf("%.*s\n", WIDTH, read + j);
            if (j < lst)
                printf("%.*s\n", lst - j, read + j);
        }

    }

}


Read *LAInterface::getRead(int number) {

    std::stringstream ss;
    std::string read_name;
    std::string read_bases;
    if (flist == NULL || findx == NULL)
        exit(1);
    HITS_READ *reads;
    HITS_TRACK *first;
    char *read, **entry;
    int c, b, e, i;
    int hilight, substr;
    int map;
    int       (*iscase)(int);
    read = New_Read_Buffer(db1);
    int UPPER = 1;
    int WIDTH = 80;
    //printf("2");
    entry = NULL;
    first = db1->tracks;
    hilight = 'A' - 'a';

    map = 0;
    reads = db1->reads;
    substr = 0;

    c = 0;

    b = number;
    e = number + 1;

    for (i = b; i < e; i++) {
        int len;
        int fst, lst;
        int flags, qv;
        HITS_READ *r;
        HITS_TRACK *track;

        r = reads + i;
        len = r->rlen;

        flags = r->flags;
        qv = (flags & DB_QV);

        {
            while (i < findx[map - 1])
                map -= 1;
            while (i >= findx[map])
                map += 1;
            ss << flist[map] << '/' << r->origin << '/' << r->fpulse << '_' << r->fpulse + len;
            if (qv > 0)
                ss << "RQ=" << qv;
        }

        ss >> read_name;

        Load_Read(db1, i, read, UPPER);

        for (track = first; track != NULL; track = track->next) {
            int64 *anno;
            int *data;
            int64 s, f, j;
            int bd, ed, m;

            anno = (int64 *) track->anno;
            data = (int *) track->data;

            s = (anno[i] >> 2);
            f = (anno[i + 1] >> 2);
            if (s < f) {
                for (j = s; j < f; j += 2) {
                    bd = data[j];
                    ed = data[j + 1];
                    for (m = bd; m < ed; m++)
                        if (iscase(read[m]))
                            read[m] = (char) (read[m] + hilight);
                    if (j == s)
                        printf("> %s:", track->name);
                    printf(" [%d,%d]", bd, ed);
                }
                printf("\n");
            }
        }

        read_bases = std::string(read);
        fst = 0;
        lst = len;


    }
    Read *new_r = new Read(number, read_name, read_bases);
    return new_r;
}


int LAInterface::OpenAlignment(std::string filename) {
    db2 = db1;

    char *fn = new char[filename.size() + 1];
    strcpy(fn, filename.c_str());

    input = Fopen(fn, "r");
    if (input == NULL)
        exit(1);

    if (fread(&novl, sizeof(int64), 1, input) != 1) SYSTEM_ERROR
    if (fread(&tspace, sizeof(int), 1, input) != 1) SYSTEM_ERROR

    if (tspace <= TRACE_XOVR) {
        small = 1;
        tbytes = sizeof(uint8);
    }
    else {
        small = 0;
        tbytes = sizeof(uint16);
    }

    printf("\n%s: ", fn);
    Print_Number(novl, 0, stdout);
    printf(" records\n");


    return 0;
}


void LAInterface::showAlignment(int from, int to) {
    int j;
    uint16 *trace;
    Work_Data *work;
    int tmax;
    int in, npt, idx, ar;
    int64 tps;
    char *abuffer, *bbuffer;
    int ar_wide, br_wide;
    int ai_wide, bi_wide;
    int mn_wide, mx_wide;
    int tp_wide;
    int blast, match, seen, lhalf, rhalf;
    bool ALIGN = true;
    bool REFERENCE = false;
    bool CARTOON = false;
    bool OVERLAP = true;
    bool FLIP = false;
    bool UPPERCASE = false;
    bool MAP = false;
    int INDENT = 4;
    int WIDTH = 100;
    int BORDER = 10;

    aln->path = &(ovl->path);
    if (ALIGN || REFERENCE) {
        work = New_Work_Data();
        abuffer = New_Read_Buffer(db1);
        bbuffer = New_Read_Buffer(db2);
    }
    else {
        abuffer = NULL;
        bbuffer = NULL;
        work = NULL;
    }

    tmax = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16) * tmax, "Allocating trace vector");
    if (trace == NULL)
        exit(1);
    in = 0;

    //if (pts!=NULL) free(pts);
    //pts = NULL;
    pts = new int[4];
    pts[0] = from + 1;
    pts[1] = to;
    pts[2] = INT32_MAX;

    npt = pts[0];
    idx = 1;

    ar_wide = Number_Digits((int64) db1->nreads);
    br_wide = Number_Digits((int64) db2->nreads);
    ai_wide = Number_Digits((int64) db1->maxlen);
    bi_wide = Number_Digits((int64) db2->maxlen);
    if (db1->maxlen < db2->maxlen) {
        mn_wide = ai_wide;
        mx_wide = bi_wide;
        tp_wide = Number_Digits((int64) db1->maxlen / tspace + 2);
    }
    else {
        mn_wide = bi_wide;
        mx_wide = ai_wide;
        tp_wide = Number_Digits((int64) db2->maxlen / tspace + 2);
    }
    ar_wide += (ar_wide - 1) / 3;
    br_wide += (br_wide - 1) / 3;
    ai_wide += (ai_wide - 1) / 3;
    bi_wide += (bi_wide - 1) / 3;
    mn_wide += (mn_wide - 1) / 3;
    tp_wide += (tp_wide - 1) / 3;
    if (FLIP) {
        int x;
        x = ar_wide;
        ar_wide = br_wide;
        br_wide = x;
        x = ai_wide;
        ai_wide = bi_wide;
        bi_wide = x;
    }

    //  For each record do

    blast = -1;
    match = 0;
    seen = 0;
    lhalf = rhalf = 0;

    for (j = 0; j < novl; j++)

        //  Read it in

    {
        //printf("j:%d/%d\n",j,novl);
        Read_Overlap(input, ovl);
        if (ovl->path.tlen > tmax) {
            tmax = ((int) 1.2 * ovl->path.tlen) + 100;
            trace = (uint16 *) Realloc(trace, sizeof(uint16) * tmax, "Allocating trace vector");
            if (trace == NULL)
                exit(1);
        }
        ovl->path.trace = (void *) trace;
        Read_Trace(input, ovl, tbytes);
        //  Determine if it should be displayed

        ar = ovl->aread + 1;
        if (in) {
            while (ar > npt) {
                npt = pts[idx++];
                if (ar < npt) {
                    in = 0;
                    break;
                }
                npt = pts[idx++];
            }
        }
        else {
            while (ar >= npt) {
                npt = pts[idx++];
                if (ar <= npt) {
                    in = 1;
                    break;
                }
                npt = pts[idx++];
            }
        }
        if (!in)
            continue;

        //  If -o check display only overlaps

        aln->alen = db1->reads[ovl->aread].rlen;
        aln->blen = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps = ovl->path.tlen / 2;

        if (OVERLAP) {
            if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
                continue;
            if (ovl->path.aepos != aln->alen && ovl->path.bepos != aln->blen)
                continue;
        }

        //  If -M option then check the completeness of the implied mapping

        if (MAP) {
            while (ovl->bread != blast) {
                if (!match && seen && !(lhalf && rhalf)) {
                    printf("Missing ");
                    Print_Number((int64) blast + 1, br_wide + 1, stdout);
                    printf(" %d ->%lld\n", db2->reads[blast].rlen, db2->reads[blast].coff);
                }
                match = 0;
                seen = 0;
                lhalf = rhalf = 0;
                blast += 1;
            }
            seen = 1;
            if (ovl->path.abpos == 0)
                rhalf = 1;
            if (ovl->path.aepos == aln->alen)
                lhalf = 1;
            if (ovl->path.bbpos != 0 || ovl->path.bepos != aln->blen)
                continue;
            match = 1;
        }

        //  Display it

        if (ALIGN || CARTOON || REFERENCE)
            printf("\n");
        if (FLIP) {
            Flip_Alignment(aln, 0);
            Print_Number((int64) ovl->bread + 1, ar_wide + 1, stdout);
            printf("  ");
            Print_Number((int64) ovl->aread + 1, br_wide + 1, stdout);
        }
        else {
            Print_Number((int64) ovl->aread , ar_wide + 1, stdout);
            printf("  ");
            Print_Number((int64) ovl->bread , br_wide + 1, stdout);
        }
        if (COMP(ovl->flags))
            printf(" c");
        else
            printf(" n");
        printf("   [");
        Print_Number((int64) ovl->path.abpos, ai_wide, stdout);
        printf("..");
        Print_Number((int64) ovl->path.aepos, ai_wide, stdout);
        printf("]%d x [",aln->alen);
        Print_Number((int64) ovl->path.bbpos, bi_wide, stdout);
        printf("..");
        Print_Number((int64) ovl->path.bepos, bi_wide, stdout);
        printf("]%d", aln->blen);

        if (ALIGN || CARTOON || REFERENCE) {
            if (ALIGN || REFERENCE) {
                char *aseq, *bseq;
                int amin, amax;
                int bmin, bmax;

                if (FLIP)
                    Flip_Alignment(aln, 0);
                if (small)
                    Decompress_TraceTo16(ovl);

                amin = ovl->path.abpos - BORDER;
                if (amin < 0) amin = 0;
                amax = ovl->path.aepos + BORDER;
                if (amax > aln->alen) amax = aln->alen;
                if (COMP(aln->flags)) {
                    bmin = (aln->blen - ovl->path.bepos) - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = (aln->blen - ovl->path.bbpos) + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }
                else {
                    bmin = ovl->path.bbpos - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = ovl->path.bepos + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }

                aseq = Load_Subread(db1, ovl->aread, amin, amax, abuffer, 0);
                bseq = Load_Subread(db2, ovl->bread, bmin, bmax, bbuffer, 0);

                aln->aseq = aseq - amin;
                if (COMP(aln->flags)) {
                    Complement_Seq(bseq, bmax - bmin);
                    aln->bseq = bseq - (aln->blen - bmax);
                }
                else
                    aln->bseq = bseq - bmin;

                Compute_Trace_PTS(aln, work, tspace);

                if (FLIP) {
                    if (COMP(aln->flags)) {
                        Complement_Seq(aseq, amax - amin);
                        Complement_Seq(bseq, bmax - bmin);
                        aln->aseq = aseq - (aln->alen - amax);
                        aln->bseq = bseq - bmin;
                    }
                    Flip_Alignment(aln, 1);
                }
            }
            if (CARTOON) {
                printf("  (");
                Print_Number(tps, tp_wide, stdout);
                printf(" trace pts)\n\n");
                Alignment_Cartoon(stdout, aln, INDENT, mx_wide);
            }
            else {
                printf(" :   = ");
                Print_Number((int64) ovl->path.diffs, mn_wide, stdout);
                printf(" diffs  (");
                Print_Number(tps, tp_wide, stdout);
                printf(" trace pts)\n");
            }
            if (REFERENCE)
                Print_Reference(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
            if (ALIGN)
                Print_Alignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
        }
        else {
            printf(" :   < ");
            Print_Number((int64) ovl->path.diffs, mn_wide, stdout);
            printf(" diffs  (");
            Print_Number(tps, tp_wide, stdout);
            printf(" trace pts)\n");
        }
    }

    free(trace);
    if (ALIGN) {
        free(bbuffer - 1);
        free(abuffer - 1);
        Free_Work_Data(work);
    }


}


void LAInterface::getAlignmentB(std::vector<int> &result, int from) {

    int j;
    uint16 *trace;
    Work_Data *work;
    int tmax;
    int in, npt, idx, ar;
    int64 tps;
    char *abuffer, *bbuffer;
    int ar_wide, br_wide;
    int ai_wide, bi_wide;
    int mn_wide, mx_wide;
    int tp_wide;
    int blast, match, seen, lhalf, rhalf;
    bool ALIGN = false;
    bool REFERENCE = false;
    bool CARTOON = false;
    bool OVERLAP = true;
    bool FLIP = false;
    bool UPPERCASE = false;
    bool MAP = false;
    int INDENT = 4;
    int WIDTH = 100;
    int BORDER = 10;

    aln->path = &(ovl->path);
    if (ALIGN || REFERENCE) {
        work = New_Work_Data();
        abuffer = New_Read_Buffer(db1);
        bbuffer = New_Read_Buffer(db2);
    }
    else {
        abuffer = NULL;
        bbuffer = NULL;
        work = NULL;
    }

    tmax = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16) * tmax, "Allocating trace vector");
    if (trace == NULL)
        exit(1);
    in = 0;

    //if (pts!=NULL) free(pts);
    //pts = NULL;
    pts = new int[4];
    pts[0] = from + 1;
    pts[1] = from + 1;
    pts[2] = INT32_MAX;

    npt = pts[0];
    idx = 1;

    ar_wide = Number_Digits((int64) db1->nreads);
    br_wide = Number_Digits((int64) db2->nreads);
    ai_wide = Number_Digits((int64) db1->maxlen);
    bi_wide = Number_Digits((int64) db2->maxlen);
    if (db1->maxlen < db2->maxlen) {
        mn_wide = ai_wide;
        mx_wide = bi_wide;
        tp_wide = Number_Digits((int64) db1->maxlen / tspace + 2);
    }
    else {
        mn_wide = bi_wide;
        mx_wide = ai_wide;
        tp_wide = Number_Digits((int64) db2->maxlen / tspace + 2);
    }
    ar_wide += (ar_wide - 1) / 3;
    br_wide += (br_wide - 1) / 3;
    ai_wide += (ai_wide - 1) / 3;
    bi_wide += (bi_wide - 1) / 3;
    mn_wide += (mn_wide - 1) / 3;
    tp_wide += (tp_wide - 1) / 3;
    if (FLIP) {
        int x;
        x = ar_wide;
        ar_wide = br_wide;
        br_wide = x;
        x = ai_wide;
        ai_wide = bi_wide;
        bi_wide = x;
    }

    //  For each record do

    blast = -1;
    match = 0;
    seen = 0;
    lhalf = rhalf = 0;

    for (j = 0; j < novl; j++)

        //  Read it in

    {
        //printf("j:%d/%d\n",j,novl);
        Read_Overlap(input, ovl);
        if (ovl->path.tlen > tmax) {
            tmax = ((int) 1.2 * ovl->path.tlen) + 100;
            trace = (uint16 *) Realloc(trace, sizeof(uint16) * tmax, "Allocating trace vector");
            if (trace == NULL)
                exit(1);
        }
        ovl->path.trace = (void *) trace;
        Read_Trace(input, ovl, tbytes);
        //  Determine if it should be displayed

        ar = ovl->aread + 1;
        if (in) {
            while (ar > npt) {
                npt = pts[idx++];
                if (ar < npt) {
                    in = 0;
                    break;
                }
                npt = pts[idx++];
            }
        }
        else {
            while (ar >= npt) {
                npt = pts[idx++];
                if (ar <= npt) {
                    in = 1;
                    break;
                }
                npt = pts[idx++];
            }
        }
        if (!in)
            continue;

        //  If -o check display only overlaps

        aln->alen = db1->reads[ovl->aread].rlen;
        aln->blen = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps = ovl->path.tlen / 2;

        if (OVERLAP) {
            if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
                continue;
            if (ovl->path.aepos != aln->alen && ovl->path.bepos != aln->blen)
                continue;
        }

        //  If -M option then check the completeness of the implied mapping

        if (MAP) {
            while (ovl->bread != blast) {
                if (!match && seen && !(lhalf && rhalf)) {
                    printf("Missing ");
                    Print_Number((int64) blast + 1, br_wide + 1, stdout);
                    printf(" %d ->%lld\n", db2->reads[blast].rlen, db2->reads[blast].coff);
                }
                match = 0;
                seen = 0;
                lhalf = rhalf = 0;
                blast += 1;
            }
            seen = 1;
            if (ovl->path.abpos == 0)
                rhalf = 1;
            if (ovl->path.aepos == aln->alen)
                lhalf = 1;
            if (ovl->path.bbpos != 0 || ovl->path.bepos != aln->blen)
                continue;
            match = 1;
        }

        //  Display it

        if (ALIGN || CARTOON || REFERENCE)
            printf("\n");
        if (FLIP) {
            Flip_Alignment(aln, 0);
            //Print_Number((int64) ovl->bread+1,ar_wide+1,stdout);
            //printf("  ");
            //Print_Number((int64) ovl->aread+1,br_wide+1,stdout);
        }
        else { //Print_Number((int64) ovl->aread+1,ar_wide+1,stdout);
            //printf("  ");
            //Print_Number((int64) ovl->bread+1,br_wide+1,stdout);
            result.push_back(ovl->bread);
        }
        //if (COMP(ovl->flags))
        //  printf(" c");
        //else
        //  printf(" n");
        //printf("   [");
        //Print_Number((int64) ovl->path.abpos,ai_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.aepos,ai_wide,stdout);
        //printf("] x [");
        //Print_Number((int64) ovl->path.bbpos,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.bepos,bi_wide,stdout);
        //printf("]");

        if ((ALIGN || CARTOON || REFERENCE) && (false)) {
            if (ALIGN || REFERENCE) {
                char *aseq, *bseq;
                int amin, amax;
                int bmin, bmax;

                if (FLIP)
                    Flip_Alignment(aln, 0);
                if (small)
                    Decompress_TraceTo16(ovl);

                amin = ovl->path.abpos - BORDER;
                if (amin < 0) amin = 0;
                amax = ovl->path.aepos + BORDER;
                if (amax > aln->alen) amax = aln->alen;
                if (COMP(aln->flags)) {
                    bmin = (aln->blen - ovl->path.bepos) - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = (aln->blen - ovl->path.bbpos) + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }
                else {
                    bmin = ovl->path.bbpos - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = ovl->path.bepos + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }

                aseq = Load_Subread(db1, ovl->aread, amin, amax, abuffer, 0);
                bseq = Load_Subread(db2, ovl->bread, bmin, bmax, bbuffer, 0);

                aln->aseq = aseq - amin;
                if (COMP(aln->flags)) {
                    Complement_Seq(bseq, bmax - bmin);
                    aln->bseq = bseq - (aln->blen - bmax);
                }
                else
                    aln->bseq = bseq - bmin;

                Compute_Trace_PTS(aln, work, tspace);

                if (FLIP) {
                    if (COMP(aln->flags)) {
                        Complement_Seq(aseq, amax - amin);
                        Complement_Seq(bseq, bmax - bmin);
                        aln->aseq = aseq - (aln->alen - amax);
                        aln->bseq = bseq - bmin;
                    }
                    Flip_Alignment(aln, 1);
                }
            }
            if (CARTOON) {
                printf("  (");
                Print_Number(tps, tp_wide, stdout);
                printf(" trace pts)\n\n");
                Alignment_Cartoon(stdout, aln, INDENT, mx_wide);
            }
            else {
                printf(" :   = ");
                Print_Number((int64) ovl->path.diffs, mn_wide, stdout);
                printf(" diffs  (");
                Print_Number(tps, tp_wide, stdout);
                printf(" trace pts)\n");
            }
            if (REFERENCE)
                Print_Reference(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
            if (ALIGN)
                Print_Alignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
        }
        else {// printf(" :   < ");
            // Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
            // printf(" diffs  (");
            // Print_Number(tps,tp_wide,stdout);
            // printf(" trace pts)\n");
        }
    }

    free(trace);
    if (ALIGN) {
        free(bbuffer - 1);
        free(abuffer - 1);
        Free_Work_Data(work);
    }

}


void LAInterface::getRead(std::vector<Read *> &reads_vec, int from, int to) {

    std::stringstream ss;
    std::string read_name;
    std::string read_bases;
    if (flist == NULL || findx == NULL)
        exit(1);
    HITS_READ *reads;
    HITS_TRACK *first;
    char *read, **entry;
    int c, b, e, i;
    int hilight, substr;
    int map;
    int       (*iscase)(int);
    read = New_Read_Buffer(db1);
    int UPPER = 1;
    int WIDTH = 80;
    entry = NULL;
    first = db1->tracks;
    hilight = 'A' - 'a';

    map = 0;
    reads = db1->reads;
    substr = 0;

    c = 0;

    b = from;
    e = to;

    for (i = b; i < e; i++) {
        int len;
        int fst, lst;
        int flags, qv;
        HITS_READ *r;
        HITS_TRACK *track;

        r = reads + i;
        len = r->rlen;

        flags = r->flags;
        qv = (flags & DB_QV);

        {
            while (i < findx[map - 1])
                map -= 1;
            while (i >= findx[map])
                map += 1;
            ss << flist[map] << '/' << r->origin << '/' << r->fpulse << '_' << r->fpulse + len;
            if (qv > 0)
                ss << "RQ=" << qv;
        }

        ss >> read_name;

        Load_Read(db1, i, read, UPPER);

        for (track = first; track != NULL; track = track->next) {
            int64 *anno;
            int *data;
            int64 s, f, j;
            int bd, ed, m;

            anno = (int64 *) track->anno;
            data = (int *) track->data;

            s = (anno[i] >> 2);
            f = (anno[i + 1] >> 2);
            if (s < f) {
                for (j = s; j < f; j += 2) {
                    bd = data[j];
                    ed = data[j + 1];
                    for (m = bd; m < ed; m++)
                        if (iscase(read[m]))
                            read[m] = (char) (read[m] + hilight);
                    if (j == s)
                        printf("> %s:", track->name);
                    printf(" [%d,%d]", bd, ed);
                }
                printf("\n");
            }
        }

        read_bases = std::string(read);
        fst = 0;
        lst = len;
        Read *new_r = new Read(i, read_name, read_bases);
        reads_vec.push_back(new_r);

    }

}


void LAInterface::resetAlignment() {
    rewind(input);

    if (fread(&novl, sizeof(int64), 1, input) != 1) SYSTEM_ERROR
    if (fread(&tspace, sizeof(int), 1, input) != 1) SYSTEM_ERROR

    if (tspace <= TRACE_XOVR) {
        small = 1;
        tbytes = sizeof(uint8);
    }
    else {
        small = 0;
        tbytes = sizeof(uint16);
    }

    //printf("\n%s: ", "read again");
    //Print_Number(novl, 0, stdout);
    //printf(" records\n");


}


void LAInterface::getOverlap(std::vector<LOverlap *> &result_vec, int from, int to) {

    int j;
    uint16 *trace;
    Work_Data *work;
    int tmax;
    int in, npt, idx, ar;
    int64 tps;
    char *abuffer, *bbuffer;
    int ar_wide, br_wide;
    int ai_wide, bi_wide;
    int mn_wide, mx_wide;
    int tp_wide;
    int blast, match, seen, lhalf, rhalf;
    bool ALIGN = false;
    bool REFERENCE = false;
    bool CARTOON = false;
    bool OVERLAP = false;
    bool FLIP = false;
    bool UPPERCASE = false;
    bool MAP = false;
    int INDENT = 4;
    int WIDTH = 100;
    int BORDER = 10;

    aln->path = &(ovl->path);
    if (ALIGN || REFERENCE) {
        work = New_Work_Data();
        abuffer = New_Read_Buffer(db1);
        bbuffer = New_Read_Buffer(db2);
    }
    else {
        abuffer = NULL;
        bbuffer = NULL;
        work = NULL;
    }

    tmax = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16) * tmax, "Allocating trace vector");
    if (trace == NULL)
        exit(1);
    in = 0;

    //if (pts!=NULL) free(pts);
    //pts = NULL;
    pts = new int[4];
    pts[0] = from + 1;
    pts[1] = to + 0;
    pts[2] = INT32_MAX;

    npt = pts[0];
    idx = 1;

    ar_wide = Number_Digits((int64) db1->nreads);
    br_wide = Number_Digits((int64) db2->nreads);
    ai_wide = Number_Digits((int64) db1->maxlen);
    bi_wide = Number_Digits((int64) db2->maxlen);
    if (db1->maxlen < db2->maxlen) {
        mn_wide = ai_wide;
        mx_wide = bi_wide;
        tp_wide = Number_Digits((int64) db1->maxlen / tspace + 2);
    }
    else {
        mn_wide = bi_wide;
        mx_wide = ai_wide;
        tp_wide = Number_Digits((int64) db2->maxlen / tspace + 2);
    }
    ar_wide += (ar_wide - 1) / 3;
    br_wide += (br_wide - 1) / 3;
    ai_wide += (ai_wide - 1) / 3;
    bi_wide += (bi_wide - 1) / 3;
    mn_wide += (mn_wide - 1) / 3;
    tp_wide += (tp_wide - 1) / 3;
    if (FLIP) {
        int x;
        x = ar_wide;
        ar_wide = br_wide;
        br_wide = x;
        x = ai_wide;
        ai_wide = bi_wide;
        bi_wide = x;
    }

    //  For each record do

    blast = -1;
    match = 0;
    seen = 0;
    lhalf = rhalf = 0;

    for (j = 0; j < novl; j++)

        //  Read it in
    {
        //printf("j:%d/%d\n",j,novl);
        Read_Overlap(input, ovl);
        if (ovl->path.tlen > tmax) {
            tmax = ((int) 1.2 * ovl->path.tlen) + 100;
            trace = (uint16 *) Realloc(trace, sizeof(uint16) * tmax, "Allocating trace vector");
            if (trace == NULL)
                exit(1);
        }
        ovl->path.trace = (void *) trace;
        Read_Trace(input, ovl, tbytes);
        //  Determine if it should be displayed

        ar = ovl->aread + 1;
        if (in) {
            while (ar > npt) {
                npt = pts[idx++];
                if (ar < npt) {
                    in = 0;
                    break;
                }
                npt = pts[idx++];
            }
        }
        else {
            while (ar >= npt) {
                npt = pts[idx++];
                if (ar <= npt) {
                    in = 1;
                    break;
                }
                npt = pts[idx++];
            }
        }
        if (!in)
            continue;

        //  If -o check display only overlaps
        aln->alen = db1->reads[ovl->aread].rlen;
        aln->blen = db1->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps = ovl->path.tlen / 2;
        LOverlap *new_ovl = new LOverlap();
/*
        if (OVERLAP) {
            if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
                continue;
            if (ovl->path.aepos != aln->alen && ovl->path.bepos != aln->blen)
                continue;
        }

        //  If -M option then check the completeness of the implied mapping

        if (MAP) {
            while (ovl->bread != blast) {
                if (!match && seen && !(lhalf && rhalf)) {
                    printf("Missing ");
                    Print_Number((int64) blast + 1, br_wide + 1, stdout);
                    printf(" %d ->%lld\n", db2->reads[blast].rlen, db2->reads[blast].coff);
                }
                match = 0;
                seen = 0;
                lhalf = rhalf = 0;
                blast += 1;
            }
            seen = 1;
            if (ovl->path.abpos == 0)
                rhalf = 1;
            if (ovl->path.aepos == aln->alen)
                lhalf = 1;
            if (ovl->path.bbpos != 0 || ovl->path.bepos != aln->blen)
                continue;
            match = 1;
        }
*/
        //  Display it

        //if (ALIGN || CARTOON || REFERENCE)
        //printf("\n");
        //printf(" %d ",j);
        if (FLIP) {
            Flip_Alignment(aln, 0);
            //Print_Number((int64) ovl->bread+1,ar_wide+1,stdout);
            //printf("  ");
            //Print_Number((int64) ovl->aread+1,br_wide+1,stdout);
        }
        else { 
			//Print_Number((int64) ovl->aread+1,ar_wide+1,stdout);
            //printf("  ");
            //Print_Number((int64) ovl->bread+1,br_wide+1,stdout);
            new_ovl->aid = ovl->aread;
            new_ovl->bid = ovl->bread;
        }
        if (COMP(ovl->flags))
        {   //printf(" c");
            new_ovl->flags = 1;
        }
        else {
            new_ovl->flags = 0;
            //printf(" n");
        }
        //printf("   [");
        //Print_Number((int64) ovl->path.abpos,ai_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.aepos,ai_wide,stdout);
        //printf("] x [");
        //Print_Number((int64) ovl->path.bbpos,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.bepos,bi_wide,stdout);
        //printf("]");
        //printf("%d",aln->blen);

        new_ovl->abpos = ovl->path.abpos;
		new_ovl->aepos = ovl->path.aepos;
        new_ovl->bbpos = ovl->path.bbpos;
        new_ovl->bepos = ovl->path.bepos;
        new_ovl->alen = aln->alen;
        new_ovl->blen = aln->blen;
        new_ovl->diffs = ovl->path.diffs;
        new_ovl->tlen = ovl->path.tlen;
        new_ovl->tps = tps;
        new_ovl->addtype();
        result_vec.push_back(new_ovl);
        if ((ALIGN || CARTOON || REFERENCE) && false) {
            if (ALIGN || REFERENCE) {
                char *aseq, *bseq;
                int amin, amax;
                int bmin, bmax;

                if (FLIP)
                    Flip_Alignment(aln, 0);
                if (small)
                    Decompress_TraceTo16(ovl);

                amin = ovl->path.abpos - BORDER;
                if (amin < 0) amin = 0;
                amax = ovl->path.aepos + BORDER;
                if (amax > aln->alen) amax = aln->alen;
                if (COMP(aln->flags)) {
                    bmin = (aln->blen - ovl->path.bepos) - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = (aln->blen - ovl->path.bbpos) + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }
                else {
                    bmin = ovl->path.bbpos - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = ovl->path.bepos + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }

                aseq = Load_Subread(db1, ovl->aread, amin, amax, abuffer, 0);
                bseq = Load_Subread(db2, ovl->bread, bmin, bmax, bbuffer, 0);

                aln->aseq = aseq - amin;
                if (COMP(aln->flags)) {
                    Complement_Seq(bseq, bmax - bmin);
                    aln->bseq = bseq - (aln->blen - bmax);
                }
                else
                    aln->bseq = bseq - bmin;

                Compute_Trace_PTS(aln, work, tspace);

                if (FLIP) {
                    if (COMP(aln->flags)) {
                        Complement_Seq(aseq, amax - amin);
                        Complement_Seq(bseq, bmax - bmin);
                        aln->aseq = aseq - (aln->alen - amax);
                        aln->bseq = bseq - bmin;
                    }
                    Flip_Alignment(aln, 1);
                }
            }
            if (CARTOON) {
                printf("  (");
                Print_Number(tps, tp_wide, stdout);
                printf(" trace pts)\n\n");
                Alignment_Cartoon(stdout, aln, INDENT, mx_wide);
            }
            else {
                printf(" :   = ");
                Print_Number((int64) ovl->path.diffs, mn_wide, stdout);
                printf(" diffs  (");
                Print_Number(tps, tp_wide, stdout);
                printf(" trace pts)\n");
            }
            if (REFERENCE)
                Print_Reference(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
            if (ALIGN)
                Print_Alignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
        }
        else {// printf(" :   < ");
            // Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
            // printf(" diffs  (");
            // Print_Number(tps,tp_wide,stdout);
            // printf(" trace pts)\n");
        }
    }

    free(trace);
    if (ALIGN) {
        free(bbuffer - 1);
        free(abuffer - 1);
        Free_Work_Data(work);
    }

}


void LAInterface::getOverlap(std::vector<LOverlap *> &result_vec, int n) {

    getOverlap(result_vec, n, n + 1);

}

void LAInterface::getAlignment(std::vector<LAlignment *> &result_vec, int from) {

    getAlignment(result_vec, from, from + 1);

}

void LAInterface::getAlignment(std::vector<LAlignment *> &result_vec, int from, int to) {

    int j;
    uint16 *trace;
    Work_Data *work;
    int tmax;
    int in, npt, idx, ar;
    int64 tps;
    char *abuffer, *bbuffer;
    int ar_wide, br_wide;
    int ai_wide, bi_wide;
    int mn_wide, mx_wide;
    int tp_wide;
    int blast, match, seen, lhalf, rhalf;
    bool ALIGN = false;
    bool REFERENCE = false;
    bool CARTOON = false;
    bool OVERLAP = false;
    bool FLIP = false;
    bool UPPERCASE = false;
    bool MAP = false;
    int INDENT = 4;
    int WIDTH = 100;
    int BORDER = 10;

    aln->path = &(ovl->path);
    if (ALIGN || REFERENCE) {
        work = New_Work_Data();
        abuffer = New_Read_Buffer(db1);
        bbuffer = New_Read_Buffer(db2);
    }
    else {
        abuffer = NULL;
        bbuffer = NULL;
        work = NULL;
    }

    tmax = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16) * tmax, "Allocating trace vector");
    if (trace == NULL)
        exit(1);
    in = 0;

    //if (pts!=NULL) free(pts);
    //pts = NULL;
    pts = new int[4];
    pts[0] = from + 1;
    pts[1] = to ;
    pts[2] = INT32_MAX;

    npt = pts[0];
    idx = 1;

    ar_wide = Number_Digits((int64) db1->nreads);
    br_wide = Number_Digits((int64) db2->nreads);
    ai_wide = Number_Digits((int64) db1->maxlen);
    bi_wide = Number_Digits((int64) db2->maxlen);
    if (db1->maxlen < db2->maxlen) {
        mn_wide = ai_wide;
        mx_wide = bi_wide;
        tp_wide = Number_Digits((int64) db1->maxlen / tspace + 2);
    }
    else {
        mn_wide = bi_wide;
        mx_wide = ai_wide;
        tp_wide = Number_Digits((int64) db2->maxlen / tspace + 2);
    }
    ar_wide += (ar_wide - 1) / 3;
    br_wide += (br_wide - 1) / 3;
    ai_wide += (ai_wide - 1) / 3;
    bi_wide += (bi_wide - 1) / 3;
    mn_wide += (mn_wide - 1) / 3;
    tp_wide += (tp_wide - 1) / 3;
    if (FLIP) {
        int x;
        x = ar_wide;
        ar_wide = br_wide;
        br_wide = x;
        x = ai_wide;
        ai_wide = bi_wide;
        bi_wide = x;
    }

    //  For each record do

    blast = -1;
    match = 0;
    seen = 0;
    lhalf = rhalf = 0;

    for (j = 0; j < novl; j++)

        //  Read it in

    {
        //printf("j:%d/%d\n",j,novl);
        Read_Overlap(input, ovl);
        if (ovl->path.tlen > tmax) {
            tmax = ((int) 1.2 * ovl->path.tlen) + 100;
            trace = (uint16 *) Realloc(trace, sizeof(uint16) * tmax, "Allocating trace vector");
            if (trace == NULL)
                exit(1);
        }
        ovl->path.trace = (void *) trace;
        Read_Trace(input, ovl, tbytes);
        //  Determine if it should be displayed

        ar = ovl->aread + 1;
        if (in) {
            while (ar > npt) {
                npt = pts[idx++];
                if (ar < npt) {
                    in = 0;
                    break;
                }
                npt = pts[idx++];
            }
        }
        else {
            while (ar >= npt) {
                npt = pts[idx++];
                if (ar <= npt) {
                    in = 1;
                    break;
                }
                npt = pts[idx++];
            }
        }
        if (!in)
            continue;

        //  If -o check display only overlaps

        aln->alen = db1->reads[ovl->aread].rlen;
        aln->blen = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps = ovl->path.tlen / 2;
        LAlignment *new_al = new LAlignment();
        new_al->aid = ovl->aread;
        new_al->bid = ovl->bread;

        if (COMP(ovl->flags))
            //printf(" c");
            new_al->flags = 1;
        else
            new_al->flags = 0;
        //printf(" n");
        //printf("   [");
        //Print_Number((int64) ovl->path.abpos,ai_wide,stdout);
        new_al->abpos = ovl->path.abpos;
        //printf("..");
        //Print_Number((int64) ovl->path.aepos,ai_wide,stdout);
        new_al->aepos = ovl->path.aepos;
        //printf("] x [");
        //Print_Number((int64) ovl->path.bbpos,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.bepos,bi_wide,stdout);
        //printf("]");
        new_al->bbpos = ovl->path.bbpos;
        new_al->bepos = ovl->path.bepos;
        new_al->alen = aln->alen;
        new_al->blen = aln->blen;
        new_al->diffs = ovl->path.diffs;
        new_al->tlen = ovl->path.tlen;
        new_al->tps = tps;
        new_al->trace = (uint16 *) Realloc(trace, sizeof(uint16) * tmax, "Allocating trace vector");
        if (new_al->trace == NULL)
            exit(1);
        memcpy(new_al->trace, (void *) trace, sizeof(uint16) * tmax);


        if (OVERLAP) {
            if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
                continue;
            if (ovl->path.aepos != aln->alen && ovl->path.bepos != aln->blen)
                continue;
        }

        //  If -M option then check the completeness of the implied mapping

        if (MAP) {
            while (ovl->bread != blast) {
                if (!match && seen && !(lhalf && rhalf)) {
                    printf("Missing ");
                    Print_Number((int64) blast + 1, br_wide + 1, stdout);
                    printf(" %d ->%lld\n", db2->reads[blast].rlen, db2->reads[blast].coff);
                }
                match = 0;
                seen = 0;
                lhalf = rhalf = 0;
                blast += 1;
            }
            seen = 1;
            if (ovl->path.abpos == 0)
                rhalf = 1;
            if (ovl->path.aepos == aln->alen)
                lhalf = 1;
            if (ovl->path.bbpos != 0 || ovl->path.bepos != aln->blen)
                continue;
            match = 1;
        }

        //  Display it

        //if (ALIGN || CARTOON || REFERENCE)
            //printf("\n");
        if (FLIP) {
            Flip_Alignment(aln, 0);
            //Print_Number((int64) ovl->bread+1,ar_wide+1,stdout);
            //printf("  ");
            //Print_Number((int64) ovl->aread+1,br_wide+1,stdout);
        }
        else { //Print_Number((int64) ovl->aread+1,ar_wide+1,stdout);

            //printf("  ");
            //Print_Number((int64) ovl->bread+1,br_wide+1,stdout);
            //result.push_back(ovl->bread);
        }
        //if (COMP(ovl->flags))
        //  printf(" c");
        //else
        //  printf(" n");
        //printf("   [");
        //Print_Number((int64) ovl->path.abpos,ai_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.aepos,ai_wide,stdout);
        //printf("] x [");
        //Print_Number((int64) ovl->path.bbpos,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.bepos,bi_wide,stdout);
        //printf("]");

        if ((ALIGN || CARTOON || REFERENCE) && (true)) {
            if (ALIGN || REFERENCE) {
                char *aseq, *bseq;
                int amin, amax;
                int bmin, bmax;

                if (FLIP)
                    Flip_Alignment(aln, 0);
                if (small)
                    Decompress_TraceTo16(ovl);

                amin = ovl->path.abpos - BORDER;
                if (amin < 0) amin = 0;
                amax = ovl->path.aepos + BORDER;
                if (amax > aln->alen) amax = aln->alen;
                if (COMP(aln->flags)) {
                    bmin = (aln->blen - ovl->path.bepos) - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = (aln->blen - ovl->path.bbpos) + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }
                else {
                    bmin = ovl->path.bbpos - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = ovl->path.bepos + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }

                aseq = Load_Subread(db1, ovl->aread, amin, amax, abuffer, 0);
                bseq = Load_Subread(db2, ovl->bread, bmin, bmax, bbuffer, 0);

                //new_al->aseq = std::string(aseq);
                //new_al->bseq = std::string(bseq);


                aln->aseq = aseq - amin;
                if (COMP(aln->flags)) {
                    Complement_Seq(bseq, bmax - bmin);
                    aln->bseq = bseq - (aln->blen - bmax);
                }
                else
                    aln->bseq = bseq - bmin;

                Compute_Trace_PTS(aln, work, tspace);

                if (FLIP) {
                    if (COMP(aln->flags)) {
                        Complement_Seq(aseq, amax - amin);
                        Complement_Seq(bseq, bmax - bmin);
                        aln->aseq = aseq - (aln->alen - amax);
                        aln->bseq = bseq - bmin;
                    }
                    Flip_Alignment(aln, 1);
                }
            }
            if (CARTOON) {
                //printf("  (");
                //Print_Number(tps, tp_wide, stdout);
                //printf(" trace pts)\n\n");
                //Alignment_Cartoon(stdout, aln, INDENT, mx_wide);
            }
            else {
                //printf(" :   = ");
                //Print_Number((int64) ovl->path.diffs, mn_wide, stdout);
                //printf(" diffs  (");
                //Print_Number(tps, tp_wide, stdout);
                //printf(" trace pts)\n");
            }
            if (REFERENCE)
                Print_Reference(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
            if (ALIGN)
                Print_Alignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
        }
        else {// printf(" :   < ");
            // Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
            // printf(" diffs  (");
            // Print_Number(tps,tp_wide,stdout);
            // printf(" trace pts)\n");
        }
        result_vec.push_back(new_al);

    }

    free(trace);
    if (ALIGN) {
        free(bbuffer - 1);
        free(abuffer - 1);
        Free_Work_Data(work);
    }

}

int LAInterface::getReadNumber() {
    return db1->nreads;
}

int64 LAInterface::getAlignmentNumber() {
    resetAlignment();
    return novl;

}

void LOverlap::addtype() {
    if ((abpos < CHI_THRESHOLD) and (aepos > alen - CHI_THRESHOLD) ) {
        if (blen > alen) aln_type = COVERED;
    }
	/**
	 A:   ==========>
	 B: ===============>
	**/

    else if ((abpos > 0) and (aepos > alen - CHI_THRESHOLD) ) {
        if (bbpos < CHI_THRESHOLD) aln_type = FORWARD;
    }
	/**
	 A:   ==========>
	 B:        ==========>
	**/

    else if ( ( abpos < CHI_THRESHOLD) and (aepos < alen)) {
        if (bepos >  blen - CHI_THRESHOLD ) aln_type = BACKWARD;
    }
	/**
	 A:     ==========>
	 B:  ==========>
	**/
    else if ((bbpos < CHI_THRESHOLD) and (bepos > blen - CHI_THRESHOLD) ) {
        if (alen > blen) aln_type = COVERING;
    }
	
	else if ((abpos > 0) and (bbpos<CHI_THRESHOLD)) {
		aln_type = MISMATCH_RIGHT;
	}
	/**
	 A:   ======..xxxx>
	 B:      ===..xxxx===>
	**/
	
	else if ((aepos < alen) and (bepos > blen - CHI_THRESHOLD)) {
		aln_type = MISMATCH_LEFT;
	}
	/**
	 A:   		xxxx..===>
	 B:     ====xxxx..=>
	**/
}



void LAInterface::showOverlap(int from, int to) {
	int        j;
	    uint16    *trace;
	    Work_Data *work;
	    int        tmax;
	    int        in, npt, idx, ar;
	    int64      tps;

	    char      *abuffer, *bbuffer;
	    int        ar_wide, br_wide;
	    int        ai_wide, bi_wide;
	    int        mn_wide, mx_wide;
	    int        tp_wide;
	    int        blast, match, seen, lhalf, rhalf;
	    bool ALIGN = false;
	    bool REFERENCE = false;
	    bool CARTOON = false;
	    bool OVERLAP = false;
	    bool FLIP = false;
	    bool UPPERCASE = false;
	    bool MAP = false;
	    int INDENT = 4;
	    int WIDTH = 100;
	    int BORDER = 10;

	    aln->path = &(ovl->path);
	    if (ALIGN || REFERENCE)
	      { work = New_Work_Data();
	        abuffer = New_Read_Buffer(db1);
	        bbuffer = New_Read_Buffer(db2);
	      }
	    else
	      { abuffer = NULL;
	        bbuffer = NULL;
	        work = NULL;
	      }

	    tmax  = 1000;
	    trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
	    if (trace == NULL)
	      exit (1);

	    in  = 0;
	    npt = pts[0];
	    idx = 1;

	    ar_wide = Number_Digits((int64) db1->nreads);
	    br_wide = Number_Digits((int64) db2->nreads);
	    ai_wide = Number_Digits((int64) db1->maxlen);
	    bi_wide = Number_Digits((int64) db2->maxlen);
	    if (db1->maxlen < db2->maxlen)
	      { mn_wide = ai_wide;
	        mx_wide = bi_wide;
	        tp_wide = Number_Digits((int64) db1->maxlen/tspace+2);
	      }
	    else
	      { mn_wide = bi_wide;
	        mx_wide = ai_wide;
	        tp_wide = Number_Digits((int64) db2->maxlen/tspace+2);
	      }
	    ar_wide += (ar_wide-1)/3;
	    br_wide += (br_wide-1)/3;
	    ai_wide += (ai_wide-1)/3;
	    bi_wide += (bi_wide-1)/3;
	    mn_wide += (mn_wide-1)/3;
	    tp_wide += (tp_wide-1)/3;

	    if (FLIP)
	      { int x;
	        x = ar_wide; ar_wide = br_wide; br_wide = x;
	        x = ai_wide; ai_wide = bi_wide; bi_wide = x;
	      }

	    //  For each record do

	    blast = -1;
	    match = 0;
	    seen  = 0;
	    lhalf = rhalf = 0;
		
		
	    pts = new int[4];
	    pts[0] = from + 1;
	    pts[1] = to ;
	    pts[2] = INT32_MAX;

	    npt = pts[0];
	    idx = 1;
		
	    for (j = 0; j < novl; j++)

	       //  Read it in

	      { Read_Overlap(input,ovl);
	        if (ovl->path.tlen > tmax)
	          { tmax = ((int) 1.2*ovl->path.tlen) + 100;
	            trace = (uint16 *) Realloc(trace,sizeof(uint16)*tmax,"Allocating trace vector");
	            if (trace == NULL)
	              exit (1);
	          }
	        ovl->path.trace = (void *) trace;
	        Read_Trace(input,ovl,tbytes);

	        //  Determine if it should be displayed

	        ar = ovl->aread+1;
	        if (in)
	          { while (ar > npt)
	              { npt = pts[idx++];
	                if (ar < npt)
	                  { in = 0;
	                    break;
	                  }
	                npt = pts[idx++];
	              }
	          }
	        else
	          { while (ar >= npt)
	              { npt = pts[idx++];
	                if (ar <= npt)
	                  { in = 1;
	                    break;
	                  }
	                npt = pts[idx++];
	              }
	          }
	        if (!in)
	          continue;

	        //  If -o check display only overlaps

	        aln->alen  = db1->reads[ovl->aread].rlen;
	        aln->blen  = db2->reads[ovl->bread].rlen;
	        aln->flags = ovl->flags;
	        tps        = ovl->path.tlen/2;

	        if (OVERLAP)
	          { if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
	              continue;
	            if (ovl->path.aepos != aln->alen && ovl->path.bepos != aln->blen)
	              continue;
	          }

	        //  If -M option then check the completeness of the implied mapping

	        if (MAP)
	          { while (ovl->bread != blast)
	              { if (!match && seen && !(lhalf && rhalf))
	                  { printf("Missing ");
	                    Print_Number((int64) blast+1,br_wide+1,stdout);
	                    printf(" %d ->%lld\n",db2->reads[blast].rlen,db2->reads[blast].coff);
	                  }
	                match = 0;
	                seen  = 0; 
	                lhalf = rhalf = 0;
	                blast += 1;
	              }
	            seen = 1;
	            if (ovl->path.abpos == 0)
	              rhalf = 1;
	            if (ovl->path.aepos == aln->alen)
	              lhalf = 1;
	            if (ovl->path.bbpos != 0 || ovl->path.bepos != aln->blen)
	              continue;
	            match = 1;
	          }

	        //  Display it
            
	        if (ALIGN || CARTOON || REFERENCE)
	          printf("\n");
	        if (FLIP)
	          { Flip_Alignment(aln,0);
	            Print_Number((int64) ovl->bread+1,ar_wide+1,stdout);
	            printf("  ");
	            Print_Number((int64) ovl->aread+1,br_wide+1,stdout);
	          }
	        else
	          { Print_Number((int64) ovl->aread+1,ar_wide+1,stdout);
	            printf("  ");
	            Print_Number((int64) ovl->bread+1,br_wide+1,stdout);
	          }
	        if (COMP(ovl->flags))
	          printf(" c");
	        else
	          printf(" n");
	        printf("   [");
	        Print_Number((int64) ovl->path.abpos,ai_wide,stdout);
	        printf("..");
	        Print_Number((int64) ovl->path.aepos,ai_wide,stdout);
	        printf("] x [");
	        Print_Number((int64) ovl->path.bbpos,bi_wide,stdout);
	        printf("..");
	        Print_Number((int64) ovl->path.bepos,bi_wide,stdout);
	        printf("]%d",aln->blen);

	        if (ALIGN || CARTOON || REFERENCE)
	          { if (ALIGN || REFERENCE)
	              { char *aseq, *bseq;
	                int   amin,  amax;
	                int   bmin,  bmax;

	                if (FLIP)
	                  Flip_Alignment(aln,0);
	                if (small)
	                  Decompress_TraceTo16(ovl);

	                amin = ovl->path.abpos - BORDER;
	                if (amin < 0) amin = 0;
	                amax = ovl->path.aepos + BORDER;
	                if (amax > aln->alen) amax = aln->alen;
	                if (COMP(aln->flags))
	                  { bmin = (aln->blen-ovl->path.bepos) - BORDER;
	                    if (bmin < 0) bmin = 0;
	                    bmax = (aln->blen-ovl->path.bbpos) + BORDER;
	                    if (bmax > aln->blen) bmax = aln->blen;
	                  }
	                else
	                  { bmin = ovl->path.bbpos - BORDER;
	                    if (bmin < 0) bmin = 0;
	                    bmax = ovl->path.bepos + BORDER;
	                    if (bmax > aln->blen) bmax = aln->blen;
	                  }

	                aseq = Load_Subread(db1,ovl->aread,amin,amax,abuffer,0);
	                bseq = Load_Subread(db2,ovl->bread,bmin,bmax,bbuffer,0);

	                aln->aseq = aseq - amin;
	                if (COMP(aln->flags))
	                  { Complement_Seq(bseq,bmax-bmin);
	                    aln->bseq = bseq - (aln->blen - bmax);
	                  }
	                else
	                  aln->bseq = bseq - bmin;

	                Compute_Trace_PTS(aln,work,tspace);

	                if (FLIP)
	                  { if (COMP(aln->flags))
	                      { Complement_Seq(aseq,amax-amin);
	                        Complement_Seq(bseq,bmax-bmin);
	                        aln->aseq = aseq - (aln->alen - amax);
	                        aln->bseq = bseq - bmin;
	                      }
	                    Flip_Alignment(aln,1);
	                  }
	              }
	            if (CARTOON)
	              { printf("  (");
	                Print_Number(tps,tp_wide,stdout);
	                printf(" trace pts)\n\n");
	                Alignment_Cartoon(stdout,aln,INDENT,mx_wide);
	              }
	            else
	              { printf(" :   = ");
	                Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
	                printf(" diffs  (");
	                Print_Number(tps,tp_wide,stdout);
	                printf(" trace pts)\n");
	              }
	            if (REFERENCE)
	              Print_Reference(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide);
	            if (ALIGN)
	              Print_Alignment(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide);
	          }
	        else
	          { printf(" :   < ");
	            Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
	            printf(" diffs  (");
	            Print_Number(tps,tp_wide,stdout);
	            printf(" trace pts)\n");
	          }
	      }

	    free(trace);
	    if (ALIGN)
	      { free(bbuffer-1);
	        free(abuffer-1);
	        Free_Work_Data(work);
	      }
	  
		  return;
}

