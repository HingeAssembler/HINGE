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
#include <algorithm>
#include "LAInterface.h"
#include "align.h"
#include "DB.h"

#include <unordered_map>
#include <tgmath.h>
#include <iomanip>
#include <fstream>


void Read::showRead() {
    std::cout << "read #" << id << std::endl;
    std::cout << ">" << name << std::endl;
    std::cout << bases << std::endl;
}


int LAInterface::openDB2(std::string filename, std::string filename2) {
    char *fn = new char[filename.length() + 1];
    strcpy(fn, filename.c_str());

    char *fn2 = new char[filename2.length() + 1];
    strcpy(fn2, filename2.c_str());

    int status = Open_DB(fn, this->db1);
    if (status < 0)
        exit(1);
    if (this->db1->part > 0) {
        fprintf(stderr, "%s: Cannot be called on a block: %s\n", "test", fn);
        exit(1);
    }


    status = Open_DB(fn2, this->db2);
    if (status < 0)
        exit(1);
    if (this->db2->part > 0) {
        fprintf(stderr, "%s: Cannot be called on a block: %s\n", "test", fn);
        exit(1);
    }

    Trim_DB(db1);
    Trim_DB(db2);



    char *fn_1 = new char[filename.length() + 1 + 3];
    strcpy(fn_1, fn);
    strcat(fn_1, ".db");

    FILE * dstub = Fopen(fn_1, "r");
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



    char *fn_2 = new char[filename2.length() + 1 + 3];
    strcpy(fn_2, fn2);
    strcat(fn_2, ".db");

    dstub = Fopen(fn_2, "r");
    if (dstub == NULL)
        exit(1);

    if (fscanf(dstub, DB_NFILE, &nfiles2) != 1) SYSTEM_ERROR

    printf("%d files\n", nfiles2);

    flist2 = (char **) Malloc(sizeof(char *) * nfiles2, "Allocating file list");
    findx2 = (int *) Malloc(sizeof(int *) * (nfiles2 + 1), "Allocating file index");

    if (flist2 == NULL || findx2 == NULL)
        exit(1);

    findx2 += 1;
    findx2[-1] = 0;

    for (int i = 0; i < nfiles2; i++) {
        char prolog[MAX_NAME], fname[MAX_NAME];

        if (fscanf(dstub, DB_FDATA, findx2 + i, fname, prolog) != 3) SYSTEM_ERROR
        if ((flist2[i] = Strdup(prolog, "Adding to file list")) == NULL)
            exit(1);
    }

    fclose(dstub);


    delete [] fn;
    delete [] fn2;
    delete [] fn_1;
    delete [] fn_2;
    return 0;
}


int LAInterface::openDB(std::string filename) {
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

int LAInterface::closeDB() {
    Close_DB(db1);
    return 0;
}


int LAInterface::closeDB2() {
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


void LAInterface::showRead2(int from, int to) {
    if (flist2 == NULL || findx2 == NULL)
        exit(1);
    HITS_READ *reads;
    HITS_TRACK *first;
    char *read, **entry;
    int c, b, e, i;
    int hilight, substr;
    int map;
    int       (*iscase)(int);

    read = New_Read_Buffer(db2);
    int UPPER = 1;
    int WIDTH = 80;
    //printf("2");
    {
        entry = NULL;
        first = db2->tracks;
    }


    hilight = 'A' - 'a';
    iscase = islower;

    map = 0;
    reads = db2->reads;
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


        Load_Read(db2, i, read, UPPER);

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

Read *LAInterface::getRead2(int number) {

    std::stringstream ss;
    std::string read_name;
    std::string read_bases;
    if (flist2 == NULL || findx2 == NULL)
        exit(1);
    HITS_READ *reads;
    HITS_TRACK *first;
    char *read, **entry;
    int c, b, e, i;
    int hilight, substr;
    int map;
    int       (*iscase)(int);
    read = New_Read_Buffer(db2);
    int UPPER = 1;
    int WIDTH = 80;
    //printf("2");
    entry = NULL;
    first = db2->tracks;
    hilight = 'A' - 'a';

    map = 0;
    reads = db2->reads;
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
            while (i < findx2[map - 1])
                map -= 1;
            while (i >= findx2[map])
                map += 1;
            ss << flist2[map] << '/' << r->origin << '/' << r->fpulse << '_' << r->fpulse + len;
            if (qv > 0)
                ss << "RQ=" << qv;
        }

        ss >> read_name;

        Load_Read(db2, i, read, UPPER);

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


int LAInterface::openAlignmentFile(std::string filename) {

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

                computeTracePTS(aln, work, tspace);

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
                printAlignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
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
        //if (COMP(ovl->reverse_complement_match_))
        //  printf(" c");
        //else
        //  printf(" n");
        //printf("   [");
        //Print_Number((int64) ovl->path.read_A_match_start_,ai_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.read_A_match_end_,ai_wide,stdout);
        //printf("] x [");
        //Print_Number((int64) ovl->path.read_B_match_start_,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.read_B_match_end_,bi_wide,stdout);
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

                Compute_Trace_PTS(aln, work, tspace,GREEDIEST);

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
        Read *new_r = new Read(i, len, read_name, read_bases);
        reads_vec.push_back(new_r);

    }

}


void LAInterface::getRead2(std::vector<Read *> &reads_vec, int from, int to) {

    std::stringstream ss;
    std::string read_name;
    std::string read_bases;
    if (flist2 == NULL || findx2 == NULL)
        exit(1);
    HITS_READ *reads;
    HITS_TRACK *first;
    char *read, **entry;
    int c, b, e, i;
    int hilight, substr;
    int map;
    int       (*iscase)(int);
    read = New_Read_Buffer(db2);
    int UPPER = 1;
    int WIDTH = 80;
    entry = NULL;
    first = db2->tracks;
    hilight = 'A' - 'a';

    map = 0;
    reads = db2->reads;
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
            while (i < findx2[map - 1])
                map -= 1;
            while (i >= findx2[map])
                map += 1;
            ss << flist2[map] << '/' << r->origin << '/' << r->fpulse << '_' << r->fpulse + len;
            if (qv > 0)
                ss << "RQ=" << qv;
        }

        ss >> read_name;

        Load_Read(db2, i, read, UPPER);

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
        Read *new_r = new Read(i, len, read_name, read_bases);
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
    int tmax;
    int in, npt, idx, ar;
    int64 tps;

    aln->path = &(ovl->path);

    tmax = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16) * tmax, "Allocating trace vector");
    if (trace == NULL)
        exit(1);
    in = 0;

    pts = new int[4];
    pts[0] = from + 1;
    pts[1] = to + 0;
    pts[2] = INT32_MAX;

    npt = pts[0];
    idx = 1;

    //  For each record do

    for (j = 0; j < novl; j++)
        //  Read it in
    {
        if (j % (novl/100) == 0) {
            printf("%d percent finished\n", j/(novl/100));
        }
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

        aln->alen = db1->reads[ovl->aread].rlen;
        aln->blen = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps = ovl->path.tlen / 2;
        LOverlap *new_ovl = new LOverlap();

        if (COMP(ovl->flags))
        {   new_ovl->reverse_complement_match_ = 1;
        }
        else {
            new_ovl->reverse_complement_match_ = 0;
        }

        if (small)
            Decompress_TraceTo16(ovl);

        new_ovl->trace_pts_len = ovl->path.tlen;
        new_ovl->trace_pts = (uint16 *)malloc(ovl->path.tlen * sizeof(uint16));

        memcpy(new_ovl->trace_pts, ovl->path.trace, ovl->path.tlen * sizeof(uint16));

        new_ovl->read_A_id_ = ovl->aread;
        new_ovl->read_B_id_ = ovl->bread;
        new_ovl->read_A_match_start_ = ovl->path.abpos;
		new_ovl->read_A_match_end_ = ovl->path.aepos;
        new_ovl->alen = aln->alen;
        new_ovl->blen = aln->blen;

        if (new_ovl->reverse_complement_match_ == 0) {
            new_ovl->read_B_match_start_ = ovl->path.bbpos;
            new_ovl->read_B_match_end_ = ovl->path.bepos;
        }
        else {
            new_ovl->read_B_match_start_ = new_ovl->blen - ovl->path.bepos;
            new_ovl->read_B_match_end_ = new_ovl->blen - ovl->path.bbpos;
        }

        new_ovl->diffs = ovl->path.diffs;
        new_ovl->tlen = ovl->path.tlen;
        new_ovl->tps = tps;
        result_vec.push_back(new_ovl);
    }
    free(trace);
}

void LAInterface::getOverlapw(std::vector<LOverlap *> &result_vec, int from, int to) {

    int j;
    uint16 *trace;
    int tmax;
    int in, npt, idx, ar;
    int64 tps;

    aln->path = &(ovl->path);

    tmax = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16) * tmax, "Allocating trace vector");
    if (trace == NULL)
        exit(1);
    in = 0;

    pts = new int[4];
    pts[0] = from + 1;
    pts[1] = to + 0;
    pts[2] = INT32_MAX;

    npt = pts[0];
    idx = 1;

    //  For each record do

    for (j = 0; j < novl; j++)
        //  Read it in
    {
        if (j % (novl/100) == 0) {
            printf("%d percent finished\n", j/(novl/100));
        }
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

        aln->alen = db1->reads[ovl->aread].rlen;
        aln->blen = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps = ovl->path.tlen / 2;
        LOverlap *new_ovl = new LOverlap();

        if (COMP(ovl->flags))
        {   new_ovl->reverse_complement_match_ = 1;
        }
        else {
            new_ovl->reverse_complement_match_ = 0;
        }

        if (small)
            Decompress_TraceTo16(ovl);

        new_ovl->trace_pts_len = ovl->path.tlen;
        //new_ovl->trace_pts = (uint16 *)malloc(ovl->path.tlen * sizeof(uint16));
        //memcpy(new_ovl->trace_pts, ovl->path.trace, ovl->path.tlen * sizeof(uint16));

        new_ovl->trace_pts = 0;
        new_ovl->read_A_id_ = ovl->aread;
        new_ovl->read_B_id_ = ovl->bread;
        new_ovl->read_A_match_start_ = ovl->path.abpos;
        new_ovl->read_A_match_end_ = ovl->path.aepos;
        new_ovl->read_B_match_start_ = ovl->path.bbpos;
        new_ovl->read_B_match_end_ = ovl->path.bepos;
        new_ovl->alen = aln->alen;
        new_ovl->blen = aln->blen;
        new_ovl->diffs = ovl->path.diffs;
        new_ovl->tlen = ovl->path.tlen;
        new_ovl->tps = tps;
        result_vec.push_back(new_ovl);
    }
    free(trace);
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
    bool ALIGN = true;
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

		//printf("j:%d/%d\n",j,novl);

        //  If -o check display only overlaps

        aln->alen = db1->reads[ovl->aread].rlen;
        aln->blen = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps = ovl->path.tlen / 2;
        LAlignment *new_al = new LAlignment();
        new_al->read_A_id_ = ovl->aread;
        new_al->read_B_id_ = ovl->bread;

        if (COMP(ovl->flags))
            //printf(" c");
            new_al->flags = 1;
        else
            new_al->flags = 0;
        //printf(" n");
        //printf("   [");
        //Print_Number((int64) ovl->path.read_A_match_start_,ai_wide,stdout);
        new_al->abpos = ovl->path.abpos;
        //printf("..");
        //Print_Number((int64) ovl->path.read_A_match_end_,ai_wide,stdout);
        new_al->aepos = ovl->path.aepos;
        //printf("] x [");
        //Print_Number((int64) ovl->path.read_B_match_start_,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.read_B_match_end_,bi_wide,stdout);
        //printf("]");
        new_al->bbpos = ovl->path.bbpos;
        new_al->bepos = ovl->path.bepos;
        new_al->alen = aln->alen;
        new_al->blen = aln->blen;
        new_al->diffs = ovl->path.diffs;
        new_al->tlen = ovl->path.tlen;
        new_al->tps = tps;

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
        //if (COMP(ovl->reverse_complement_match_))
        //  printf(" c");
        //else
        //  printf(" n");
        //printf("   [");
        //Print_Number((int64) ovl->path.read_A_match_start_,ai_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.read_A_match_end_,ai_wide,stdout);
        //printf("] x [");
        //Print_Number((int64) ovl->path.read_B_match_start_,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.read_B_match_end_,bi_wide,stdout);
        //printf("]");


        if ((ALIGN || CARTOON || REFERENCE)  || true) {
            if (ALIGN || REFERENCE) {
                char *aseq, *bseq;
                int amin, amax;
                int bmin, bmax;

                if (FLIP)
                    Flip_Alignment(aln, 0);
                //if (small)
                //    Decompress_TraceTo16(ovl);


                if (small)
                    Decompress_TraceTo16(ovl);

                new_al->trace_pts_len = ovl->path.tlen;
                new_al->trace_pts = (uint16 *)malloc(ovl->path.tlen * sizeof(uint16));

                memcpy(new_al->trace_pts, ovl->path.trace, ovl->path.tlen * sizeof(uint16));

                /*{
                    printf("\n");

                    uint16 *pp = (uint16 *) ovl->path.trace;
                    for (int uu = 0; uu < ovl->path.tlen; uu++) {
                        printf("%d ", pp[uu]);
                        new_al->trace_pts[uu] = pp[uu];
                    }


                    printf("\n");


                }*/
#ifdef DOALIGN
                amin = ovl->path.read_A_match_start_ - BORDER;
                if (amin < 0) amin = 0;
                amax = ovl->path.read_A_match_end_ + BORDER;
                if (amax > aln->alen) amax = aln->alen;
                if (COMP(aln->reverse_complement_match_)) {
                    bmin = (aln->blen - ovl->path.read_B_match_end_) - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = (aln->blen - ovl->path.read_B_match_start_) + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }
                else {
                    bmin = ovl->path.read_B_match_start_ - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = ovl->path.read_B_match_end_ + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }

                aseq = Load_Subread(db1, ovl->aread, amin, amax, abuffer, 0);
                bseq = Load_Subread(db2, ovl->bread, bmin, bmax, bbuffer, 0);


                aln->aseq = aseq - amin;
                if (COMP(aln->reverse_complement_match_)) {
                    Complement_Seq(bseq, bmax - bmin);
                    aln->bseq = bseq - (aln->blen - bmax);
                }
                else
                    aln->bseq = bseq - bmin;




                Compute_Trace_PTS(aln, work, tspace,GREEDIEST);


                /*new_al->aseq = (char *) malloc(new_al->alen * sizeof(char));
                new_al->bseq = (char *) malloc(new_al->blen * sizeof(char));

                memcpy(new_al->aseq, aln->aseq, new_al->alen* sizeof(char));
                memcpy(new_al->bseq, aln->bseq, new_al->blen* sizeof(char));*/
                new_al->aseq = NULL;
                new_al->bseq = NULL;


                /*{
                    int tlen = aln->path->tlen;
                    int *trace = (int *) aln->path->trace;
                    int u;
                    printf(" ");
                    for (u = 0; u < tlen; u++)
                        printf("%d,", (int) trace[u]);
                    printf("\n");
                }*/

                new_al->tlen =  aln->path->tlen;
                new_al->trace = (int *) malloc(sizeof(int) * aln->path->tlen*2);
                //if (new_al->trace == NULL)
                //    exit(1);
                //memcpy(new_al->trace, (void *) aln->path->trace, sizeof(int) * sizeof(int) * aln->path->tlen);

                //free(trace);

                //printf("after\n");
                {
                    int tlen = aln->path->tlen;
                    int *trace = (int *) aln->path->trace;
                    int u;
                    //printf(" ");
                    for (u = 0; u < tlen; u++) {
                        //printf("%d,", (int) trace[u]);
                        new_al->trace[u] = (int)trace[u];
                    }
                    //printf("\n");
                }

#endif
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
            //if (ALIGN)
                //printAlignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
            //    printAlignment_exp(stdout, new_al, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);



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


void LAInterface::getAlignment(std::vector<LAlignment *> &result_vec, std::vector<int> & range) {

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
    pts = new int[range.size()*2+20];
    for (int k = 0; k < range.size(); k++) {
        pts[k*2] = range[k] + 1;
        pts[k*2+1] = range[k] + 1;
    }
    pts[range.size()*2] = INT32_MAX;

    /*for (int i = 0; i < range.size()*2+2; i++) {
        printf("%d\n",pts[i]);
    }*/

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

        //printf("j:%d/%d\n",j,novl);

        //  If -o check display only overlaps

        aln->alen = db1->reads[ovl->aread].rlen;
        aln->blen = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;
        tps = ovl->path.tlen / 2;
        LAlignment *new_al = new LAlignment();
        new_al->read_A_id_ = ovl->aread;
        new_al->read_B_id_ = ovl->bread;

        if (COMP(ovl->flags))
            //printf(" c");
            new_al->flags = 1;
        else
            new_al->flags = 0;
        //printf(" n");
        //printf("   [");
        //Print_Number((int64) ovl->path.read_A_match_start_,ai_wide,stdout);
        new_al->abpos = ovl->path.abpos;
        //printf("..");
        //Print_Number((int64) ovl->path.read_A_match_end_,ai_wide,stdout);
        new_al->aepos = ovl->path.aepos;
        //printf("] x [");
        //Print_Number((int64) ovl->path.read_B_match_start_,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.read_B_match_end_,bi_wide,stdout);
        //printf("]");
        new_al->bbpos = ovl->path.bbpos;
        new_al->bepos = ovl->path.bepos;
        new_al->alen = aln->alen;
        new_al->blen = aln->blen;
        new_al->diffs = ovl->path.diffs;
        new_al->tlen = ovl->path.tlen;
        new_al->tps = tps;

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
        //if (COMP(ovl->reverse_complement_match_))
        //  printf(" c");
        //else
        //  printf(" n");
        //printf("   [");
        //Print_Number((int64) ovl->path.read_A_match_start_,ai_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.read_A_match_end_,ai_wide,stdout);
        //printf("] x [");
        //Print_Number((int64) ovl->path.read_B_match_start_,bi_wide,stdout);
        //printf("..");
        //Print_Number((int64) ovl->path.read_B_match_end_,bi_wide,stdout);
        //printf("]");


        if ((ALIGN || CARTOON || REFERENCE)  || true) {
            if (ALIGN || REFERENCE) {
                char *aseq, *bseq;
                int amin, amax;
                int bmin, bmax;

                if (FLIP)
                    Flip_Alignment(aln, 0);
                //if (small)
                //    Decompress_TraceTo16(ovl);


                if (small)
                    Decompress_TraceTo16(ovl);

                new_al->trace_pts_len = ovl->path.tlen;
                new_al->trace_pts = (uint16 *)malloc(ovl->path.tlen * sizeof(uint16));
                memcpy(new_al->trace_pts, ovl->path.trace, ovl->path.tlen * sizeof(uint16));

                /*{
                    printf("\n");

                    uint16 *pp = (uint16 *) ovl->path.trace;
                    for (int uu = 0; uu < ovl->path.tlen; uu++) {
                        printf("%d ", pp[uu]);
                        new_al->trace_pts[uu] = pp[uu];
                    }


                    printf("\n");


                }*/


//#define DOALIGN
#ifdef DOALIGN
                amin = ovl->path.read_A_match_start_ - BORDER;
                if (amin < 0) amin = 0;
                amax = ovl->path.read_A_match_end_ + BORDER;
                if (amax > aln->alen) amax = aln->alen;
                if (COMP(aln->reverse_complement_match_)) {
                    bmin = (aln->blen - ovl->path.read_B_match_end_) - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = (aln->blen - ovl->path.read_B_match_start_) + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }
                else {
                    bmin = ovl->path.read_B_match_start_ - BORDER;
                    if (bmin < 0) bmin = 0;
                    bmax = ovl->path.read_B_match_end_ + BORDER;
                    if (bmax > aln->blen) bmax = aln->blen;
                }

                aseq = Load_Subread(db1, ovl->aread, amin, amax, abuffer, 0);
                bseq = Load_Subread(db2, ovl->bread, bmin, bmax, bbuffer, 0);


                aln->aseq = aseq - amin;
                if (COMP(aln->reverse_complement_match_)) {
                    Complement_Seq(bseq, bmax - bmin);
                    aln->bseq = bseq - (aln->blen - bmax);
                }
                else
                    aln->bseq = bseq - bmin;




                Compute_Trace_PTS(aln, work, tspace,GREEDIEST);

#endif
                /*new_al->aseq = (char *) malloc(new_al->alen * sizeof(char));
                new_al->bseq = (char *) malloc(new_al->blen * sizeof(char));

                memcpy(new_al->aseq, aln->aseq, new_al->alen* sizeof(char));
                memcpy(new_al->bseq, aln->bseq, new_al->blen* sizeof(char));*/
                //new_al->aseq = NULL;
                //new_al->bseq = NULL;


                /*{
                    int tlen = aln->path->tlen;
                    int *trace = (int *) aln->path->trace;
                    int u;
                    printf(" ");
                    for (u = 0; u < tlen; u++)
                        printf("%d,", (int) trace[u]);
                    printf("\n");
                }*/

#ifdef DOALIGN
                new_al->tlen =  aln->path->tlen;
                new_al->trace = (int *) malloc(sizeof(int) * aln->path->tlen*2);
                //if (new_al->trace == NULL)
                //    exit(1);
                //memcpy(new_al->trace, (void *) aln->path->trace, sizeof(int) * sizeof(int) * aln->path->tlen);

                //free(trace);

                //printf("after\n");
                {
                    int tlen = aln->path->tlen;
                    int *trace = (int *) aln->path->trace;
                    int u;
                    //printf(" ");
                    for (u = 0; u < tlen; u++) {
                        //printf("%d,", (int) trace[u]);
                        new_al->trace[u] = (int)trace[u];
                    }
                    //printf("\n");
                }

#endif
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
            //if (ALIGN)
            //printAlignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);
            //    printAlignment_exp(stdout, new_al, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide);



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


int LAInterface::getReadNumber2() {
    return db2->nreads;
}

int64 LAInterface::getAlignmentNumber() {
    resetAlignment();
    return novl;

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

	                Compute_Trace_PTS(aln,work,tspace,GREEDIEST);

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

typedef struct            //  Hidden from the user, working space for each thread
{ int     vecmax;
    void   *vector;
    int     celmax;
    void   *cells;
    int     pntmax;
    void   *points;
    int     tramax;
    void   *trace;
} _Work_Data;

static int enlarge_vector(_Work_Data *work, int newmax)
{ void *vec;
    int   max;

    max = ((int) (newmax*1.2)) + 10000;
    vec = Realloc(work->vector,max,"Enlarging DP vector");
    if (vec == NULL)
        EXIT(1);
    work->vecmax = max;
    work->vector = vec;
    return (0);
}

static char ToL[8] = { 'a', 'c', 'g', 't', '.', '[', ']', '-' };
static char ToU[8] = { 'A', 'C', 'G', 'T', '.', '[', ']', '-' };

int LAInterface::printAlignment(FILE *file, Alignment *align, Work_Data *ework,
                                int indent, int width, int border, int upper, int coord)
{ _Work_Data *work  = (_Work_Data *) ework;
    int        *trace = (int *) align->path->trace;
    int         tlen  = align->path->tlen;

    char *Abuf, *Bbuf, *Dbuf;
    int   i, j, o;
    char *a, *b;
    char  mtag, dtag;
    int   prefa, prefb;
    int   aend, bend;
    int   sa, sb;
    int   match, diff;
    char *N2A;

    if (trace == NULL) return (0);

#ifdef SHOW_TRACE
    fprintf(file,"\nTrace:\n");
  for (i = 0; i < tlen; i++)
    fprintf(file,"  %3d\n",trace[i]);
#endif

    o = sizeof(char)*3*(width+1);
    if (o > work->vecmax)
    if (enlarge_vector(work,o))
        EXIT(1);

    if (upper)
        N2A = ToU;
    else
        N2A = ToL;

    Abuf = (char *) work->vector;
    Bbuf = Abuf + (width+1);
    Dbuf = Bbuf + (width+1);

    aend = align->path->aepos;
    bend = align->path->bepos;

    Abuf[width] = Bbuf[width] = Dbuf[width] = '\0';
    /* buffer/output next column */
#define COLUMN(x,y)							\
{ int u, v;								\
  if (o >= width)							\
    { fprintf(file,"\n");						\
      fprintf(file,"%*s",indent,"");					\
      if (coord > 0)							\
        { if (sa <= aend)						\
            fprintf(file," %*d",coord,sa);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %s\n",Abuf);					\
          fprintf(file,"%*s %*s %s\n",indent,"",coord,"",Dbuf);		\
          fprintf(file,"%*s",indent,"");				\
          if (sb <= bend)						\
            fprintf(file," %*d",coord,sb);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %s",Bbuf);					\
        }								\
      else								\
        { fprintf(file," %s\n",Abuf);					\
          fprintf(file,"%*s %s\n",indent,"",Dbuf);			\
          fprintf(file,"%*s %s",indent,"",Bbuf);			\
        }								\
      fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));		\
      o  = 0;								\
      sa = i;								\
      sb = j;								\
      match = diff = 0;							\
    }									\
  u = (x);								\
  v = (y);								\
  if (u == 4 || v == 4)							\
    Dbuf[o] = ' ';							\
  else if (u == v)							\
    Dbuf[o] = mtag;							\
  else									\
    Dbuf[o] = dtag;							\
  Abuf[o] = N2A[u];							\
  Bbuf[o] = N2A[v];							\
  o += 1;								\
}

    a = align->aseq - 1;
    b = align->bseq - 1;

    o  = 0;
    i = j = 1;

    prefa = align->path->abpos;
    prefb = align->path->bbpos;

    if (prefa > border)
    { i = prefa-(border-1);
        prefa = border;
    }
    if (prefb > border)
    { j = prefb-(border-1);
        prefb = border;
    }

    sa   = i;
    sb   = j;
    mtag = ':';
    dtag = ':';

    while (prefa > prefb)
    { COLUMN(a[i],4)
        i += 1;
        prefa -= 1;
    }
    while (prefb > prefa)
    { COLUMN(4,b[j])
        j += 1;
        prefb -= 1;
    }
    while (prefa > 0)
    { COLUMN(a[i],b[j])
        i += 1;
        j += 1;
        prefa -= 1;
    }

    mtag = '[';
    if (prefb > 0)
    COLUMN(5,5)

    mtag  = '|';
    dtag  = '*';

    match = diff = 0;

    { int p, c;      /* Output columns of alignment til reach trace end */

        for (c = 0; c < tlen; c++)
            if ((p = trace[c]) < 0)
            { p = -p;
                while (i != p)
                { COLUMN(a[i],b[j])
                    if (a[i] == b[j])
                        match += 1;
                    else
                        diff += 1;
                    i += 1;
                    j += 1;
                }
                COLUMN(7,b[j])
                j += 1;
                diff += 1;
            }
            else
            { while (j != p)
                { COLUMN(a[i],b[j])
                    if (a[i] == b[j])
                        match += 1;
                    else
                        diff += 1;
                    i += 1;
                    j += 1;
                }
                COLUMN(a[i],7)
                i += 1;
                diff += 1;
            }
        p = align->path->aepos;
        while (i <= p)
        { COLUMN(a[i],b[j])
            if (a[i] == b[j])
                match += 1;
            else
                diff += 1;
            i += 1;
            j += 1;
        }
    }

    { int c;     /* Output remaining column including unaligned suffix */

        mtag = ']';
        if (a[i] != 4 && b[j] != 4 && border > 0)
        COLUMN(6,6)

        mtag = ':';
        dtag = ':';

        c = 0;
        while (c < border && (a[i] != 4 || b[j] != 4))
        { if (a[i] != 4)
            if (b[j] != 4)
            { COLUMN(a[i],b[j])
                i += 1;
                j += 1;
            }
            else
            { COLUMN(a[i],4)
                i += 1;
            }
            else
            { COLUMN(4,b[j])
                j += 1;
            }
            c += 1;
        }
    }

    /* Print remainder of buffered col.s */

    fprintf(file,"\n");
    fprintf(file,"%*s",indent,"");
    if (coord > 0)
    { if (sa <= aend)
            fprintf(file," %*d",coord,sa);
        else
            fprintf(file," %*s",coord,"");
        fprintf(file," %.*s\n",o,Abuf);
        fprintf(file,"%*s %*s %.*s\n",indent,"",coord,"",o,Dbuf);
        fprintf(file,"%*s",indent,"");
        if (sb <= bend)
            fprintf(file," %*d",coord,sb);
        else
            fprintf(file," %*s",coord,"");
        fprintf(file," %.*s",o,Bbuf);
    }
    else
    { fprintf(file," %.*s\n",o,Abuf);
        fprintf(file,"%*s %.*s\n",indent,"",o,Dbuf);
        fprintf(file,"%*s %.*s",indent,"",o,Bbuf);
    }
    if (diff+match > 0)
        fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));
    else
        fprintf(file,"\n");

    //fprintf(file, "Cool!\n");
    fflush(file);
    return (0);
}


typedef void Work_Data;

typedef struct
{ int  *Stop;          //  Ongoing stack of alignment indels
    char *Aabs, *Babs;   //  Absolute base of A and B sequences

    int  **PVF, **PHF;   //  List of waves for iterative np algorithms
    int   mida,  midb;   //  mid point division for mid-point algorithms

    int   *VF,   *VB;    //  Forward/Reverse waves for nd algorithms
    //  (defunct: were used for O(nd) algorithms)
} Trace_Waves;

static int enlarge_trace(_Work_Data *work, int newmax)
{ void *vec;
    int   max;

    max = ((int) (newmax*1.2)) + 10000;
    vec = Realloc(work->trace,max,"Enlarging trace vector");
    if (vec == NULL)
        EXIT(1);
    work->tramax = max;
    work->trace  = vec;
    return (0);
}

static int iter_np(char *A, int M, char *B, int N, Trace_Waves *wave)
{ int  **PVF = wave->PVF;
    int  **PHF = wave->PHF;
    int    D;
    int    del = M-N;

    { int  *F0, *F1, *F2;
        int  *HF;
        int   low, hgh;
        int   posl, posh;

#ifdef DEBUG_ALIGN
        printf("\n%*s BASE %ld,%ld: %d vs %d\n",depth,"",A-wave->Aabs,B-wave->Babs,M,N);
    printf("%*s A = ",depth,"");
    for (D = 0; D < M; D++)
      printf("%c",ToA[(int) A[D]]);
    printf("\n");
    printf("%*s B = ",depth,"");
    for (D = 0; D < N; D++)
      printf("%c",ToA[(int) B[D]]);
    printf("\n");
#endif

        if (del >= 0)
        { low = 0;
            hgh = del;
        }
        else
        { low = del;
            hgh = 0;
        }

        posl = -INT32_MAX;
        posh =  INT32_MAX;
        if (wave->Aabs == wave->Babs)
        { if (B == A)
            { EPRINTF(EPLACE,"Error: self comparison starts on diagonal 0 (Compute_Trace)\n");
                EXIT(-1);
            }
            else if (B < A)
                posl = (B-A)+1;
            else
                posh = (B-A)-1;
        }

        F1 = PVF[-2];
        F0 = PVF[-1];

        for (D = low-1; D <= hgh+1; D++)
            F1[D] = F0[D] = -2;
        F0[0] = -1;

        low += 1;
        hgh -= 1;

        for (D = 0; 1; D += 1)
        { int   k, i, j;
            int   am, ac, ap;
            char *a;

            F2 = F1;
            F1 = F0;
            F0 = PVF[D];
            HF = PHF[D];

            if ((D & 0x1) == 0)
            { if (low > posl)
                    low -= 1;
                if (hgh < posh)
                    hgh += 1;
            }
            F0[hgh+1] = F0[low-1] = -2;

#define FS_MOVE(mdir,pdir)			\
  ac = F1[k]+1;					\
  if (ac < am)					\
    if (ap < am)				\
      { HF[k] = mdir;				\
        j = am;					\
      }						\
    else					\
      { HF[k] = pdir;				\
        j = ap;					\
      }						\
  else						\
    if (ap < ac)				\
      { HF[k] = 0;				\
        j = ac;					\
      }						\
    else					\
      { HF[k] = pdir;				\
        j = ap;					\
      }						\
						\
  if (N < i)					\
    while (j < N && B[j] == a[j])		\
      j += 1;					\
  else						\
    while (j < i && B[j] == a[j])		\
      j += 1;					\
  F0[k] = j;

            j = -2;
            a = A + hgh;
            i = M - hgh;
            for (k = hgh; k > del; k--)
            { ap = j+1;
                am = F2[k-1];
                FS_MOVE(-1,4)
                a -= 1;
                i += 1;
            }

            j = -2;
            a = A + low;
            i = M - low;
            for (k = low; k < del; k++)
            { ap = F2[k+1]+1;
                am = j;
                FS_MOVE(2,1)
                a += 1;
                i -= 1;
            }

            ap = F0[del+1]+1;
            am = j;
            FS_MOVE(2,4)

#ifdef DEBUG_AWAVE
            print_awave(F0,low,hgh);
        print_awave(HF,low,hgh);
#endif

            if (F0[del] >= N)
                break;
        }
    }

    { int   k, h, m, e, c;
        char *a;
        int   ap = (wave->Aabs-A)-1;
        int   bp = (B-wave->Babs)+1;

        PHF[0][0] = 3;

        c = N;
        k = del;
        e = PHF[D][k];
        PHF[D][k] = 3;
        while (e != 3)
        { h = k+e;
            if (e > 1)
                h -= 3;
            else if (e == 0)
                D -= 1;
            else
                D -= 2;
            if (h < k)       // => e = -1 or 2
            { a = A + k;
                if (k < 0)
                    m = -k;
                else
                    m = 0;
                if (PVF[D][h] <= c)
                    c = PVF[D][h]-1;
                while (c >= m && a[c] == B[c])
                    c -= 1;
                if (e < 1)  //  => edge is 2, others are 1, and 0
                { if (c <= PVF[D+2][k+1])
                    { e = 4;
                        h = k+1;
                        D = D+2;
                    }
                    else if (c == PVF[D+1][k])
                    { e = 0;
                        h = k;
                        D = D+1;
                    }
                    else
                        PVF[D][h] = c+1;
                }
                else      //   => edge is 0, others are 1, and 2 (if k != del), 0 (otherwise)
                { if (k == del)
                        m = D;
                    else
                        m = D-2;
                    if (c <= PVF[m][k+1])
                    { if (k == del)
                            e = 4;
                        else
                            e = 1;
                        h = k+1;
                        D = m;
                    }
                    else if (c == PVF[D-1][k])
                    { e = 0;
                        h = k;
                        D = D-1;
                    }
                    else
                        PVF[D][h] = c+1;
                }
            }
            m = PHF[D][h];
            PHF[D][h] = e;
            e = m;

            k = h;
        }

        k = D = 0;
        e = PHF[D][k];
        while (e != 3)
        { h = k-e;
            c = PVF[D][k];
            if (e > 1)
                h += 3;
            else if (e == 0)
                D += 1;
            else
                D += 2;
            if (h > k)
                *wave->Stop++ = bp+c;
            else if (h < k)
                *wave->Stop++ = ap-(c+k);
            k = h;
            e = PHF[D][h];
        }

#ifdef DEBUG_SCRIPT
        k = D = 0;
    e = PHF[D][k];
    while (e != 3)
      { h = k-e;
        c = PVF[D][k];
        if (e > 1)
          h += 3;
        else if (e == 0)
          D += 1;
        else
          D += 2;
        if (h > k)
          printf("%*s  D %d(%d)\n",depth,"",(c-k)-(ap-1),c+bp);
        else if (h < k)
          printf("%*s  I %d(%d)\n",depth,"",c+(bp-1),(c+k)-ap);
        else
          printf("%*s  %d S %d\n",depth,"",(c+k)-(ap+1),c+(bp-1));
        k = h;
        e = PHF[D][h];
      }
#endif
    }

    return (D + abs(del));
}



int LAInterface::computeTracePTS(Alignment *align, Work_Data *ework, int trace_spacing)
{ _Work_Data *work = (_Work_Data *) ework;
    Trace_Waves wave;

    Path   *path;
    char   *aseq, *bseq;
    uint16 *points;
    int     tlen;
    int     ab, bb;
    int     ae, be;
    int     diffs;

    path   = align->path;
    aseq   = align->aseq;
    bseq   = align->bseq;
    tlen   = path->tlen;
    points = (uint16 *) path->trace;

    { int64 s;
        int   d;
        int   M, N;
        int   dmax, nmax;
        int   **PVF, **PHF;

        M = path->aepos-path->abpos;
        N = path->bepos-path->bbpos;
        if (M < N)
            s = N*sizeof(int);
        else
            s = M*sizeof(int);
        if (s > work->tramax)
        if (enlarge_trace(work,s))
            EXIT(1);

        nmax = 0;
        dmax = 0;
        for (d = 1; d < tlen; d += 2)
        { if (points[d-1] > dmax)
                dmax = points[d-1];
            if (points[d] > nmax)
                nmax = points[d];
        }
        if (tlen <= 1)
            nmax = N;
        if (points[d-1] > dmax)
            dmax = points[d-1];

        s = (dmax+3)*2*((trace_spacing+nmax+3)*sizeof(int) + sizeof(int *));

        if (s > work->vecmax)
        if (enlarge_vector(work,s))
            EXIT(1);

        wave.PVF = PVF = ((int **) (work->vector)) + 2;
        wave.PHF = PHF = PVF + (dmax+3);

        s = trace_spacing+nmax+3;
        PVF[-2] = ((int *) (PHF + (dmax+1))) + (nmax+1);
        for (d = -1; d <= dmax; d++)
            PVF[d] = PVF[d-1] + s;
        PHF[-2] = PVF[dmax] + s;
        for (d = -1; d <= dmax; d++)
            PHF[d] = PHF[d-1] + s;
    }

    wave.Stop = (int *) (work->trace);
    wave.Aabs = aseq;
    wave.Babs = bseq;

    { int i, d;

        diffs = 0;
        ab = path->abpos;
        ae = (ab/trace_spacing)*trace_spacing;
        bb = path->bbpos;
        tlen -= 2;
        for (i = 1; i < tlen; i += 2)
        { ae = ae + trace_spacing;
            be = bb + points[i];
            d  = iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave);
            if (d < 0)
                EXIT(1);
            diffs += d;
            ab = ae;
            bb = be;
        }
        ae = path->aepos;
        be = path->bepos;
        d  = iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave);
        if (d < 0)
            EXIT(1);
        diffs += d;
    }

    path->trace = work->trace;
    path->tlen  = wave.Stop - ((int *) path->trace);
    path->diffs = diffs;

    return (0);
}

int LAInterface::showAlignmentTags(LAlignment *alignment) {

    //load aseq and bseq first

    //printf("A:%s\n",alignment->aseq);
    //printf("B:%s\n",alignment->bseq);
    int amin, amax, bmin, bmax;
    const int BORDER = 10;
    amin = alignment->abpos - BORDER;
    if (amin < 0) amin = 0;
    amax = alignment->aepos + BORDER;
    if (amax > alignment->alen) amax = alignment->alen;
    if (alignment->flags == 1) {
        bmin = (alignment->blen - alignment->bepos) - BORDER;
        if (bmin < 0) bmin = 0;
        bmax = (alignment->blen - alignment->bbpos) + BORDER;
        if (bmax > alignment->blen) bmax = alignment->blen;
    }
    else {
        bmin = alignment->bbpos - BORDER;
        if (bmin < 0) bmin = 0;
        bmax = alignment->bepos + BORDER;
        if (bmax > alignment->blen) bmax = alignment->blen;
    }


    char * abuffer = New_Read_Buffer(db1);
    char * bbuffer = New_Read_Buffer(db2);


    char * aseq = Load_Subread(db1, alignment->read_A_id_, amin, amax, abuffer, 0);
    char * bseq = Load_Subread(db2, alignment->read_B_id_, bmin, bmax, bbuffer, 0);


    alignment->aseq = aseq - amin;
    if (alignment->flags == 1) {
        Complement_Seq(bseq, bmax - bmin);
        alignment->bseq = bseq - (alignment->blen - bmax);
    }
    else
        alignment->bseq = bseq - bmin;




    char *Abuf, *Bbuf, *Dbuf;
    int   i, j, o;
    char *a, *b;
    char  mtag, dtag;
    int   prefa, prefb;
    int   aend, bend;
    int   sa, sb;
    int   match, diff;
    char *N2A;
    int border = 10;

    int tlen = alignment->tlen;
    int * trace = alignment->trace;

    a = alignment->aseq - 1;
    b = alignment->bseq - 1;

    i = j = 1;

    prefa = alignment->abpos;
    prefb = alignment->bbpos;

    if (prefa > border)
    { i = prefa-(border-1);
        prefa = border;
    }
    if (prefb > border)
    { j = prefb-(border-1);
        prefb = border;
    }

    sa   = i;
    sb   = j;
    mtag = ':';
    dtag = ':';

#define COLUMN(x,y) \
    {               \
        printf(" %c-%c ",ToU[x],ToU[y]); \
    }               \


    while (prefa > prefb)
    { COLUMN(a[i],4)
        i += 1;
        prefa -= 1;
    }
    while (prefb > prefa)
    { COLUMN(4,b[j])
        j += 1;
        prefb -= 1;
    }
    while (prefa > 0)
    { COLUMN(a[i],b[j])
        i += 1;
        j += 1;
        prefa -= 1;
    }

    mtag = '[';
    if (prefb > 0)
    COLUMN(5,5)

    mtag  = '|';
    dtag  = '*';

    match = diff = 0;

    { int p, c;      /* Output columns of alignment til reach trace end */

        for (c = 0; c < tlen; c++)
            if ((p = trace[c]) < 0)
            { p = -p;
                //printf("%d\n",trace[c]);
                while (i != p)
                { COLUMN(a[i],b[j])
                    if (a[i] == b[j])
                        match += 1;
                    else
                        diff += 1;
                    i += 1;
                    j += 1;
                }
                COLUMN(7,b[j])
                j += 1;
                diff += 1;
            }
            else
            { while (j != p)
                { COLUMN(a[i],b[j])
                    if (a[i] == b[j])
                        match += 1;
                    else
                        diff += 1;
                    i += 1;
                    j += 1;
                }
                COLUMN(a[i],7)
                i += 1;
                diff += 1;
            }
        p = alignment->aepos;
        while (i <= p)
        { COLUMN(a[i],b[j])
            if (a[i] == b[j])
                match += 1;
            else
                diff += 1;
            i += 1;
            j += 1;
        }
    }

    { int c;     /* Output remaining column including unaligned suffix */

        mtag = ']';
        if (a[i] != 4 && b[j] != 4 && border > 0)
        COLUMN(6,6)

        mtag = ':';
        dtag = ':';

        c = 0;
        while (c < border && (a[i] != 4 || b[j] != 4))
        { if (a[i] != 4)
            if (b[j] != 4)
            { COLUMN(a[i],b[j])
                i += 1;
                j += 1;
            }
            else
            { COLUMN(a[i],4)
                i += 1;
            }
            else
            { COLUMN(4,b[j])
                j += 1;
            }
            c += 1;
        }
    }

    free(abuffer - 1);
    free(bbuffer - 1);

    alignment->aseq = NULL;
    alignment->bseq = NULL;


    return 0;
}


std::pair<std::string, std::string> LAInterface::getAlignmentTags(LAlignment *alignment) {

    //load aseq and bseq first

    //printf("A:%s\n",alignment->aseq);
    //printf("B:%s\n",alignment->bseq);
    int amin, amax, bmin, bmax;
    const int BORDER = 10;
    amin = alignment->abpos - BORDER;
    if (amin < 0) amin = 0;
    amax = alignment->aepos + BORDER;
    if (amax > alignment->alen) amax = alignment->alen;
    if (alignment->flags == 1) {
        bmin = (alignment->blen - alignment->bepos) - BORDER;
        if (bmin < 0) bmin = 0;
        bmax = (alignment->blen - alignment->bbpos) + BORDER;
        if (bmax > alignment->blen) bmax = alignment->blen;
    }
    else {
        bmin = alignment->bbpos - BORDER;
        if (bmin < 0) bmin = 0;
        bmax = alignment->bepos + BORDER;
        if (bmax > alignment->blen) bmax = alignment->blen;
    }


    char * abuffer = New_Read_Buffer(db1);
    char * bbuffer = New_Read_Buffer(db2);


    char * aseq = Load_Subread(db1, alignment->read_A_id_, amin, amax, abuffer, 0);
    char * bseq = Load_Subread(db2, alignment->read_B_id_, bmin, bmax, bbuffer, 0);


    alignment->aseq = aseq - amin;
    if (alignment->flags == 1) {
        Complement_Seq(bseq, bmax - bmin);
        alignment->bseq = bseq - (alignment->blen - bmax);
    }
    else
        alignment->bseq = bseq - bmin;




    char *Abuf, *Bbuf, *Dbuf;
    int   i, j, o;
    char *a, *b;
    char  mtag, dtag;
    int   prefa, prefb;
    int   aend, bend;
    int   sa, sb;
    int   match, diff;
    char *N2A;
    int border = 10;

    int tlen = alignment->tlen;
    int * trace = alignment->trace; // get the trace from here

    a = alignment->aseq - 1;
    b = alignment->bseq - 1;

    i = j = 1;

    prefa = alignment->abpos;
    prefb = alignment->bbpos;

    if (prefa > border)
    { i = prefa-(border-1);
        prefa = border;
    }
    if (prefb > border)
    { j = prefb-(border-1);
        prefb = border;
    }

    sa   = i;
    sb   = j;
    mtag = ':';
    dtag = ':';

    std::string aa = "";
    std::string bb = "";

#define COLUMN(x,y) \
    {               \
        aa.append(1,ToU[x]); \
        bb.append(1,ToU[y]); \
    }               \


    while (prefa > prefb)
    { //COLUMN(a[i],4)
        i += 1;
        prefa -= 1;
    }
    while (prefb > prefa)
    { //COLUMN(4,b[j])
        j += 1;
        prefb -= 1;
    }
    while (prefa > 0)
    { //COLUMN(a[i],b[j])
        i += 1;
        j += 1;
        prefa -= 1;
    }

    mtag = '[';
    if (prefb > 0)
    //COLUMN(5,5)

    mtag  = '|';
    dtag  = '*';

    match = diff = 0;

    { int p, c;      /* Output columns of alignment til reach trace end */

        for (c = 0; c < tlen; c++)
            if ((p = trace[c]) < 0)
            { p = -p;
                //printf("%d\n",trace[c]);
                while (i != p)
                { COLUMN(a[i],b[j])
                    if (a[i] == b[j])
                        match += 1;
                    else
                        diff += 1;
                    i += 1;
                    j += 1;
                }
                COLUMN(7,b[j])
                j += 1;
                diff += 1;
            }
            else
            { while (j != p)
                { COLUMN(a[i],b[j])
                    if (a[i] == b[j])
                        match += 1;
                    else
                        diff += 1;
                    i += 1;
                    j += 1;
                }
                COLUMN(a[i],7)
                i += 1;
                diff += 1;
            }
        p = alignment->aepos;
        while (i <= p)
        { COLUMN(a[i],b[j])
            if (a[i] == b[j])
                match += 1;
            else
                diff += 1;
            i += 1;
            j += 1;
        }
    }

 /*   { int c;     // Output remaining column including unaligned suffix

        mtag = ']';
        if (a[i] != 4 && b[j] != 4 && border > 0)
        COLUMN(6,6)

        mtag = ':';
        dtag = ':';

        c = 0;
        while (c < border && (a[i] != 4 || b[j] != 4))
        { if (a[i] != 4)
            if (b[j] != 4)
            { COLUMN(a[i],b[j])
                i += 1;
                j += 1;
            }
            else
            { COLUMN(a[i],4)
                i += 1;
            }
            else
            { COLUMN(4,b[j])
                j += 1;
            }
            c += 1;
        }
    }*/


    //printf("%s\n%s\n", aa.c_str(), bb.c_str());
    free(abuffer - 1);
    free(bbuffer - 1);

    alignment->aseq = NULL;
    alignment->bseq = NULL;


    return std::pair<std::string, std::string>(aa,bb);
}



int LAInterface::printAlignment_exp(FILE *file, LAlignment *align, Work_Data *ework,
                                    int indent, int width, int border, int upper, int coord)
{ _Work_Data *work  = (_Work_Data *) ework;
    int        *trace = (int *) align->trace;
    int         tlen  = align->tlen;

    char *Abuf, *Bbuf, *Dbuf;
    int   i, j, o;
    char *a, *b;
    char  mtag, dtag;
    int   prefa, prefb;
    int   aend, bend;
    int   sa, sb;
    int   match, diff;
    char *N2A;

    if (trace == NULL) return (0);

#ifdef SHOW_TRACE
    fprintf(file,"\nTrace:\n");
  for (i = 0; i < tlen; i++)
    fprintf(file,"  %3d\n",trace[i]);
#endif

    o = sizeof(char)*3*(width+1);
    if (o > work->vecmax)
    if (enlarge_vector(work,o))
        EXIT(1);

    if (upper)
        N2A = ToU;
    else
        N2A = ToL;

    Abuf = (char *) work->vector;
    Bbuf = Abuf + (width+1);
    Dbuf = Bbuf + (width+1);

    aend = align->aepos;
    bend = align->bepos;

    Abuf[width] = Bbuf[width] = Dbuf[width] = '\0';
    /* buffer/output next column */
#define COLUMN(x,y)							\
{ int u, v;								\
  if (o >= width)							\
    { fprintf(file,"\n");						\
      fprintf(file,"%*s",indent,"");					\
      if (coord > 0)							\
        { if (sa <= aend)						\
            fprintf(file," %*d",coord,sa);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %s\n",Abuf);					\
          fprintf(file,"%*s %*s %s\n",indent,"",coord,"",Dbuf);		\
          fprintf(file,"%*s",indent,"");				\
          if (sb <= bend)						\
            fprintf(file," %*d",coord,sb);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %s",Bbuf);					\
        }								\
      else								\
        { fprintf(file," %s\n",Abuf);					\
          fprintf(file,"%*s %s\n",indent,"",Dbuf);			\
          fprintf(file,"%*s %s",indent,"",Bbuf);			\
        }								\
      fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));		\
      o  = 0;								\
      sa = i;								\
      sb = j;								\
      match = diff = 0;							\
    }									\
  u = (x);								\
  v = (y);								\
  if (u == 4 || v == 4)							\
    Dbuf[o] = ' ';							\
  else if (u == v)							\
    Dbuf[o] = mtag;							\
  else									\
    Dbuf[o] = dtag;							\
  Abuf[o] = N2A[u];							\
  Bbuf[o] = N2A[v];							\
  o += 1;								\
}

    a = align->aseq - 1;
    b = align->bseq - 1;

    o  = 0;
    i = j = 1;

    prefa = align->abpos;
    prefb = align->bbpos;

    if (prefa > border)
    { i = prefa-(border-1);
        prefa = border;
    }
    if (prefb > border)
    { j = prefb-(border-1);
        prefb = border;
    }

    sa   = i;
    sb   = j;
    mtag = ':';
    dtag = ':';

    while (prefa > prefb)
    { COLUMN(a[i],4)
        i += 1;
        prefa -= 1;
    }
    while (prefb > prefa)
    { COLUMN(4,b[j])
        j += 1;
        prefb -= 1;
    }
    while (prefa > 0)
    { COLUMN(a[i],b[j])
        i += 1;
        j += 1;
        prefa -= 1;
    }

    mtag = '[';
    if (prefb > 0)
    COLUMN(5,5)

    mtag  = '|';
    dtag  = '*';

    match = diff = 0;

    { int p, c;      /* Output columns of alignment til reach trace end */

        for (c = 0; c < tlen; c++)
            if ((p = trace[c]) < 0)
            { p = -p;
                while (i != p)
                { COLUMN(a[i],b[j])
                    if (a[i] == b[j])
                        match += 1;
                    else
                        diff += 1;
                    i += 1;
                    j += 1;
                }
                COLUMN(7,b[j])
                j += 1;
                diff += 1;
            }
            else
            { while (j != p)
                { COLUMN(a[i],b[j])
                    if (a[i] == b[j])
                        match += 1;
                    else
                        diff += 1;
                    i += 1;
                    j += 1;
                }
                COLUMN(a[i],7)
                i += 1;
                diff += 1;
            }
        p = align->aepos;
        while (i <= p)
        { COLUMN(a[i],b[j])
            if (a[i] == b[j])
                match += 1;
            else
                diff += 1;
            i += 1;
            j += 1;
        }
    }

    { int c;     /* Output remaining column including unaligned suffix */

        mtag = ']';
        if (a[i] != 4 && b[j] != 4 && border > 0)
        COLUMN(6,6)

        mtag = ':';
        dtag = ':';

        c = 0;
        while (c < border && (a[i] != 4 || b[j] != 4))
        { if (a[i] != 4)
            if (b[j] != 4)
            { COLUMN(a[i],b[j])
                i += 1;
                j += 1;
            }
            else
            { COLUMN(a[i],4)
                i += 1;
            }
            else
            { COLUMN(4,b[j])
                j += 1;
            }
            c += 1;
        }
    }

    /* Print remainder of buffered col.s */

    fprintf(file,"\n");
    fprintf(file,"%*s",indent,"");
    if (coord > 0)
    { if (sa <= aend)
            fprintf(file," %*d",coord,sa);
        else
            fprintf(file," %*s",coord,"");
        fprintf(file," %.*s\n",o,Abuf);
        fprintf(file,"%*s %*s %.*s\n",indent,"",coord,"",o,Dbuf);
        fprintf(file,"%*s",indent,"");
        if (sb <= bend)
            fprintf(file," %*d",coord,sb);
        else
            fprintf(file," %*s",coord,"");
        fprintf(file," %.*s",o,Bbuf);
    }
    else
    { fprintf(file," %.*s\n",o,Abuf);
        fprintf(file,"%*s %.*s\n",indent,"",o,Dbuf);
        fprintf(file,"%*s %.*s",indent,"",o,Bbuf);
    }
    if (diff+match > 0)
        fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));
    else
        fprintf(file,"\n");

    //fprintf(file, "Cool!\n");
    fflush(file);
    return (0);
}

int LAInterface::generateConsensus(std::vector<LAlignment *> &alns) {

    int seq_count = alns.size();

    //TBD


    return 0;
}

int LAInterface::recoverAlignment(LAlignment *alignment) {


    int j;
    uint16 *trace;
    Work_Data *work;
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
    bool OVERLAP = false;
    bool FLIP = false;
    bool UPPERCASE = false;
    bool MAP = false;
    int INDENT = 4;
    int WIDTH = 100;
    int BORDER = 10;

    //int tmax = 3000;
    //trace = (uint16 *) malloc(sizeof(uint16) * tmax);
    //if (trace == NULL)
    //    exit(1);

    int amin, amax, bmin, bmax;


    work = New_Work_Data();
    abuffer = New_Read_Buffer(db1);
    bbuffer = New_Read_Buffer(db2);

    Overlap * ovl = (Overlap *) malloc(sizeof(Overlap));
    Alignment * aln = (Alignment *) malloc(sizeof (Alignment));

    aln->path = &(ovl->path);
    Path * path = &(ovl->path);

    path->abpos = alignment->abpos;
    path->aepos = alignment->aepos;
    path->bbpos = alignment->bbpos;
    path->bepos = alignment->bepos;
    path->diffs = alignment->diffs;
    path->tlen = alignment->tlen;
    aln->alen = alignment->alen;
    aln->blen = alignment->blen;
    aln->flags = (uint32)alignment->flags;
    ovl->aread = alignment->read_A_id_;
    ovl->bread = alignment->read_B_id_;

    path->trace = (uint16 *)malloc(path->tlen * sizeof(uint16));
    memcpy(path->trace, alignment->trace_pts, path->tlen * sizeof(uint16));


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

    char * aseq = Load_Subread(db1, ovl->aread, amin, amax, abuffer, 0);
    char * bseq = Load_Subread(db2, ovl->bread, bmin, bmax, bbuffer, 0);


    aln->aseq = aseq - amin;
    if (COMP(aln->flags)) {
        Complement_Seq(bseq, bmax - bmin);
        aln->bseq = bseq - (aln->blen - bmax);
    }
    else
        aln->bseq = bseq - bmin;


    Compute_Trace_PTS(aln, work, tspace,GREEDIEST);


    /*{
        int tlen = aln->path->tlen;
        int *trace = (int *) aln->path->trace;
        int u;
        printf(" ");
        for (u = 0; u < tlen; u++)
            printf("%d,", (int) trace[u]);
        printf("\n");
    }*/



    alignment->tlen =  aln->path->tlen;
    alignment->trace = (int *) malloc(sizeof(int) * aln->path->tlen*2);
    {
        int tlen = aln->path->tlen;
        int *trace = (int *) aln->path->trace;
        int u;
        //printf(" ");
        for (u = 0; u < tlen; u++) {
            //printf("%d,", (int) trace[u]);
            alignment->trace[u] = (int)trace[u];
        }
        //printf("\n");
    }

    free(bbuffer - 1);
    free(abuffer - 1);
    Free_Work_Data(work);
	//free(path->trace);



    return 0;
}

std::vector<int> * LAInterface::getCoverage(std::vector<LOverlap *> alns) {
    std::vector<int> * res = new std::vector<int>( alns[0]->alen, 0 );

    for (int i = 0; i < alns.size(); i++) {
        for (int j = alns[i]->read_A_match_start_; j < alns[i]->read_A_match_end_; j++)
            (*res)[j] ++;
    }

    return res;
}


std::vector<int> *LAInterface::getCoverage(std::vector<LAlignment *> alns) {
    std::vector<int> * res = new std::vector<int>( alns[0]->alen, 0 );

    for (int i = 0; i < alns.size(); i++) {
        for (int j = alns[i]->abpos; j < alns[i]->aepos; j++)
            (*res)[j] ++;
    }

    return res;
}

std::vector<std::pair<int, int> > * LAInterface::lowCoverageRegions(std::vector<int> &cov, int min_cov) {
    std::vector<std::pair < int, int>> * reg = new std::vector<std::pair < int, int>> ();
    int pos = 0;
    while (pos < cov.size()) {
        int start = 0;
        if (cov[pos] < min_cov){
            start = pos;
            while ((cov[pos] < min_cov) and (pos < cov.size()))
                pos ++;
            reg->push_back(std::pair<int, int >(start, pos) ); //low coverage region in [a,b)
            }
        else pos ++;
    }
    return reg;
}

bool compare_event(std::pair<int, int> event1,std::pair<int, int> event2) {
    return event1.first < event2.first;
}


void LAInterface::profileCoverage(std::vector<LOverlap *> &alignments, std::vector<std::pair<int, int> > & coverage,int reso, int cutoff) {
    //Returns coverage, which is a pair of ints <i*reso, coverage at position i*reso of read a>
    std::vector<std::pair<int, int> > events;
    for (int i = 0; i < alignments.size(); i ++) {
        events.push_back(std::pair<int, int>(alignments[i]->read_A_match_start_ + cutoff, 1));
        events.push_back(std::pair<int, int>(alignments[i]->read_A_match_end_ - cutoff, -1));
    }

    std::sort(events.begin(), events.end(), compare_event);

    int pos = 0;
    int i = 0;
    int count = 0;
    while (pos < events.size()) {
        while ((events[pos].first < i*reso) and (pos < events.size())) {
            count += events[pos].second;
            pos++;
        }
        coverage.push_back(std::pair<int, int>(i*reso, count));
        i++;
    }
    return;
}


void LAInterface::profileCoveragefine(std::vector<LOverlap *> &alignments, std::vector<std::pair<int, int> > & coverage,int reso, int cutoff, int est_coverage) {
    std::vector<std::pair<int, int> > events;
    int sz = alignments.size();
    if (sz > est_coverage) sz = est_coverage;

    for (int i = 0; i < sz; i ++) {
        events.push_back(std::pair<int, int>(alignments[i]->read_A_match_start_ + cutoff, 1));
        events.push_back(std::pair<int, int>(alignments[i]->read_A_match_end_ - cutoff, -1));
    }

    std::sort(events.begin(), events.end(), compare_event);

    int pos = 0;
    int i = 0;
    int count = 0;
    while (pos < events.size()) {
        while ((events[pos].first < i*reso) and (pos < events.size())) {
            count += events[pos].second;
            pos++;
        }
        coverage.push_back(std::pair<int, int>(i*reso, count));
        i++;
    }
    return;
}



void LAInterface::repeatDetect(std::vector<std::pair<int, int> > & coverage, std::vector<std::pair<int, int> > & repeat) {
    for (int i = 1; i < coverage.size(); i++) {
        if (coverage[i].second > 2*coverage[i-1].second) repeat.push_back(std::pair<int, int>(coverage[i].first, 1));
        if (coverage[i].second < 0.5*coverage[i-1].second) repeat.push_back(std::pair<int, int>(coverage[i].first, -1));
    }
    return;
}


static int qv_map[51] =
  { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
    'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D',
    'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
    'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
    'Y'
  };
  
void LAInterface::getQV(std::vector<std::vector<int> > & QV, int from, int to) {
	int b,e;
    b = from;
    e = to;
	HITS_READ * reads  = db1->reads;
	bool UPPER = true;
	
    int64      *qv_idx;
    uint8      *qv_val;
	
    //if (DOIQV)
      { int status, kind;
        HITS_TRACK *track;
        status = Check_Track(db1,"qual",&kind);
        if (status == -2)
          { fprintf(stderr,"%s: .qual-track does not exist for this db.\n",Prog_Name);
            exit (1);
          }
        if (status == -1)
          { fprintf(stderr,"%s: .qual-track not sync'd with db.\n",Prog_Name);
            exit (1);
          }
        track = Load_Track(db1,"qual");
        qv_idx = (int64 *) track->anno;
        qv_val = (uint8 *) track->data;
      }
	  
	
    for (int i = b; i < e; i++)
      { 
		int         len;
        int         fst, lst;
        int         flags, qv;
        HITS_READ  *r;

        r   = reads + i;
        len = r->rlen;
        /*if (DORED)
          printf("R %d\n",i+1);*/

        flags = r->flags;
        qv    = (flags & DB_QV);
        /*if (DOHDR)
          { if (DAM)
              { char header[MAX_NAME];

                fseeko(hdrs,r->coff,SEEK_SET);
                fgets(header,MAX_NAME,hdrs);
                header[strlen(header)-1] = '\0';
                printf("H %ld %s\n",strlen(header),header);
                printf("L %d %d %d\n",r->origin,r->fpulse,r->fpulse+len);
              }
            else
              { while (i < findx[map-1])
                  map -= 1;
                while (i >= findx[map])
                  map += 1;
                printf("H %ld %s\n",strlen(flist[map]),flist[map]);
                printf("L %d %d %d\n",r->origin,r->fpulse,r->fpulse+len);
                if (qv > 0)
                  printf("Q: %d\n",qv);
              }
          }*/

        /*if (DOQVS)
          Load_QVentry(db,i,entry,UPPER);*/
        /*if (DOSEQ)
          Load_Read(db,i,read,UPPER);*/

        /*for (m = 0; m < MTOP; m++)
          { int64 *anno;
            int   *data;
            int64  s, f, j;

            anno = (int64 *) MTRACK[m]->anno;
            data = (int *) MTRACK[m]->data;

            s = (anno[i] >> 2);
            f = (anno[i+1] >> 2);
            printf("T%d %lld ",m,(f-s)/2);
            if (s < f)
              { for (j = s; j < f; j += 2)
                  printf(" %d %d",data[j],data[j+1]);
              }
            printf("\n");
          }

        if (substr)
          { fst = iter->beg;
            lst = iter->end;
          }
        else
          { fst = 0;
            lst = len;
          }

        if (DOSEQ)
          { printf("S %d ",lst-fst);
            printf("%.*s\n",lst-fst,read+fst);
          }
		*/
        //if (DOIQV)
          { int64 k, e;
			  std::vector<int> qv;
            k = qv_idx[i];
            e = qv_idx[i+1];
            //printf("I %lld ",e-k);
            while (k < e) {
				qv.push_back(qv_val[k++]);
                //putchar(qv_map[qv_val[k++]]);
				
            }
               //printf("\n");
			QV.push_back(qv);
          }
        /*if (DOQVS)
          { int k;

            for (k = 0; k < 5; k++)
              { printf("%c %d ",qvname[k],lst-fst);
                printf("%.*s\n",lst-fst,entry[k]+fst);
              }
          }*/
	  }
	return;
}


void LOverlap::trim_overlap() {
    /**
     * Trim overlap: the reads are trimmed according to qualities and coverage,
     * To be consistent, the overlap needs to be trimmed.
     * Rather than running DAligner on trimmed reads, this function  trims the overlap according to trace points.
     * It finds the trace point that are not trimmed in both reads.
     */

    //before trimming, the positions are read_A_match_start_, read_B_match_start_, read_A_match_end_
    // and read_B_match_end_, we add a eff_ prefix to it after trimming

    this->eff_read_B_match_start_ = this->read_B_match_start_;
    this->eff_read_B_match_end_ = this->read_B_match_end_;
    this->eff_read_A_match_start_ = this->read_A_match_start_;
    this->eff_read_A_match_end_ = this->read_A_match_end_;


    std::vector<std::pair<int,int> > trace_points;

    if (this->reverse_complement_match_ == 0) {
        trace_points.push_back(std::pair<int,int>(this->read_A_match_start_, this->read_B_match_start_));
    }
    else {
        trace_points.push_back(std::pair<int,int>(this->read_A_match_start_, this->read_B_match_end_));
    }

    int rev_sign = 1 - 2*this->reverse_complement_match_;

    int current_position_read_A = this->read_A_match_start_;
    // this for loop change trace points stored trace_pts[] into coordinate pairs vector: tps
    for (int j = 0; j < this->trace_pts_len/2-1; j++) {
        if (current_position_read_A % 100 != 0)
            current_position_read_A = int(ceil(current_position_read_A / 100.0)) * 100;
        else
            current_position_read_A += 100;
        trace_points.push_back(std::pair<int,int>(current_position_read_A,
                                                  trace_points.back().second + rev_sign * this->trace_pts[2 * j + 1]));
    }
    if (this->reverse_complement_match_ == 0) {
        trace_points.push_back(std::pair<int, int>(this->read_A_match_end_, this->read_B_match_end_));
    }
    else {
        trace_points.push_back(std::pair<int, int>(this->read_A_match_end_, this->read_B_match_start_));
    }


    //printf("[%6d %6d] [%6d %6d]\n", this->eff_read_A_start_, this->eff_read_A_end_, this->eff_read_B_start_, this->eff_read_B_end_);

    //printf("[%6d %6d] [%6d %6d]\n", this->eff_read_A_match_start_, this->eff_read_A_match_end_, this->eff_read_B_match_start_, this->eff_read_B_match_end_);

    /*for (int j = 0; j < trace_points.size(); j++) {
        printf("a%d b%d ", trace_points[j].first, trace_points[j].second);
    }
    printf("\n");

     // for debugging
    */

    this->eff_start_trace_point_index_ = trace_points.size();
    this->eff_end_trace_point_index_ = 0;



    if (this->reverse_complement_match_ == 0) {

        //for trace point pairs, get the first one that is in untrimmed regions for both reads

        for (int i = 0; i < trace_points.size(); i++) {
            if ( (trace_points[i].first >= this->eff_read_A_start_) and
                (trace_points[i].second >= this->eff_read_B_start_) ) {
                this->eff_read_A_match_start_ = trace_points[i].first;
                this->eff_read_B_match_start_ = trace_points[i].second;
                this->eff_start_trace_point_index_ = i;
                break;
            }
        }

        //for trace point pairs, get the last one that is in untrimmed regions for both reads
        for (int i = (int) trace_points.size() - 1; i >= 0; i--) {
            if ((trace_points[i].first <= this->eff_read_A_end_) and
                (trace_points[i].second <= this->eff_read_B_end_)) {
                this->eff_read_A_match_end_ = trace_points[i].first;
                this->eff_read_B_match_end_ = trace_points[i].second;
                this->eff_end_trace_point_index_ = i;
                break;
            }
        }

    }
    else {

        for (int i = 0; i < trace_points.size(); i++) {
            if ( (trace_points[i].first >= this->eff_read_A_start_) and
                 (trace_points[i].second <= this->eff_read_B_end_) ) {
                this->eff_read_A_match_start_ = trace_points[i].first;
                this->eff_read_B_match_end_ = trace_points[i].second;
                this->eff_start_trace_point_index_ = i; // "start" with respect to A
                break;
            }
        }

        for (int i = (int) trace_points.size() - 1; i >= 0; i--) {
            if ((trace_points[i].first <= this->eff_read_A_end_) and
                (trace_points[i].second >= this->eff_read_B_start_)) {
                this->eff_read_A_match_end_ = trace_points[i].first;
                this->eff_read_B_match_start_ = trace_points[i].second;
                this->eff_end_trace_point_index_ = i;
                break;
            }
        }

    }

    if (this->eff_start_trace_point_index_ >= this->eff_end_trace_point_index_)
    {
        this->active = false;
    }

    /*printf("[%6d %6d] [%6d %6d]\n", this->eff_read_A_match_start_, this->eff_read_A_match_end_, this->eff_read_B_match_start_, this->eff_read_B_match_end_);

    int overhang_read_A_left = this->eff_read_A_match_start_ - this->eff_read_A_start_;
    int overhang_read_A_right = this->eff_read_A_end_ - this->eff_read_A_match_end_;
    int overhang_read_B_left = this->eff_read_B_match_start_ - this->eff_read_B_start_;
    int overhang_read_B_right = this->eff_read_B_end_ - this->eff_read_B_match_end_;

    printf("trim A_left %6d, A_right %6d, B_left %6d, B_right %6d\n",
           overhang_read_A_left, overhang_read_A_right,
           overhang_read_B_left, overhang_read_B_right);

    */

}


void LOverlap::TrimOverlapNaive(){
    this->eff_read_B_match_start_ = std::max (this->read_B_match_start_,this->eff_read_B_start_);
    this->eff_read_B_match_end_ = std::min (this->read_B_match_end_,this->eff_read_B_end_);
    this->eff_read_A_match_start_ = std::max (this->read_A_match_start_,this->eff_read_A_start_);
    this->eff_read_A_match_end_ = std::min (this->read_A_match_end_,this->eff_read_A_end_);;
}


// This function is no longer used in hinging_v1.cpp

void LOverlap::addtype(int max_overhang) {
    /**
     * addtype is a function for classifying overlaps, edges are classified into forward, backward, internal match, bcovera and acoverb,
        it is based on effective positions, rather than positions
     */

    int overhang = std::min(this->eff_read_A_match_start_ - this->eff_read_A_start_, this->eff_read_B_match_start_ - this->eff_read_B_start_) + std::min(this->eff_read_A_end_ - this->eff_read_A_match_end_, this->eff_read_B_end_ - this->eff_read_B_match_end_);

    //int tol = 0;
    if (overhang > max_overhang)
        this->match_type_ = INTERNAL;
    else if ((this->eff_read_A_match_start_ - this->eff_read_A_start_ <= this->eff_read_B_match_start_ - this->eff_read_B_start_) and (this->eff_read_A_end_ - this->eff_read_A_match_end_ <= this->eff_read_B_end_ - this->eff_read_B_match_end_))
        this->match_type_ = BCOVERA;
    else if ((this->eff_read_A_match_start_ - this->eff_read_A_start_ >= this->eff_read_B_match_start_ - this->eff_read_B_start_) and (this->eff_read_A_end_ - this->eff_read_A_match_end_ >= this->eff_read_B_end_ - this->eff_read_B_match_end_))
        this->match_type_ = ACOVERB;
    else if (this->eff_read_A_match_start_ - this->eff_read_A_start_ > this->eff_read_B_match_start_ - this->eff_read_B_start_) {
        if ((this->eff_read_B_end_ - this->eff_read_B_match_end_ > 0) and (this->eff_read_A_match_start_ - this->eff_read_A_start_ > 0))
            this->match_type_ = FORWARD;
    }
    else {
        if ((this->eff_read_B_match_start_ - this->eff_read_B_start_ > 0) and (this->eff_read_A_end_ - this->eff_read_A_match_end_ > 0))
            this->match_type_ = BACKWARD;
    }
}

void LOverlap::AddTypesAsymmetric(int max_overhang) {
    //Getting a parameter max_overhang, which is the maximum overlap that one can attribute to bad DAligner ends
    //The function sets the class variable match_type_ according to the relative positions of the reads.
    //Possible things it can set to are:
    // BCOVERA, ACOVERB, INTERNAL, FORWARD, FORWARD_INTERNAL, BACKWARD, BACKWARD_INTERNAL
    int overhang_read_A_left = this->eff_read_A_match_start_ - this->eff_read_A_start_;
    int overhang_read_A_right = this->eff_read_A_end_ - this->eff_read_A_match_end_;
    int overhang_read_B_left = this->eff_read_B_match_start_ - this->eff_read_B_start_;
    int overhang_read_B_right = this->eff_read_B_end_ - this->eff_read_B_match_end_;


    //printf("     A_left %6d, A_right %6d, B_left %6d, B_right %6d\n",
    //       overhang_read_A_left, overhang_read_A_right,
    //       overhang_read_B_left, overhang_read_B_right);

    if (this->reverse_complement_match_ == 1) {
        //Exchange overhang left and right of read B if match is reverse complement
        overhang_read_B_left = this->eff_read_B_end_ - this->eff_read_B_match_end_;
        overhang_read_B_right = this->eff_read_B_match_start_ - this->eff_read_B_start_;
    }


    if (std::max(overhang_read_A_left, overhang_read_A_right) < max_overhang)
       // and ((overhang_read_A_left <= overhang_read_B_left)
       //      and (overhang_read_A_right <= overhang_read_B_right)))
        this->match_type_ = BCOVERA;
    else if (std::max(overhang_read_B_left, overhang_read_B_right) < max_overhang)
         //and (overhang_read_A_left >= overhang_read_B_left)
         //     and (overhang_read_A_right >= overhang_read_B_right))
    this->match_type_ = ACOVERB;
    else if ((std::min(overhang_read_A_left, overhang_read_A_right) > max_overhang))
        this->match_type_ = INTERNAL;
    else if (overhang_read_A_left <= max_overhang) {
        //Check if read B if a left extension. As we've handled internal,
        //we know that this is a BACKWARD or BACKWARD_INTERNAL match
        if ((overhang_read_B_right <= max_overhang) and (overhang_read_B_left >= max_overhang)) {
            //Alignment internal in B. (It may be an overlap or a non extending overlap)
            this->match_type_ = BACKWARD;
        }
        else if ((overhang_read_B_right >= max_overhang) and (overhang_read_B_left >= max_overhang)) {
            //Alignment is a overlap on B.
            this->match_type_ = BACKWARD_INTERNAL;
        }
    }
    else if (overhang_read_A_right <= max_overhang) {
        //Check if read B if a right extension. As we've handled internal,
        //we know that this is a FORWARD or FORWARD_INTERNAL match
        if ((overhang_read_B_left <= max_overhang) and (overhang_read_B_right >= max_overhang)) {
            //Alignment internal in B. (It may be an overlap or a non extending overlap)
            this->match_type_ = FORWARD;
        }
        else if ((overhang_read_B_left >= max_overhang) and (overhang_read_B_right >= max_overhang)) {
            //Alignment is a overlap on B.
            this->match_type_ = FORWARD_INTERNAL;
        }
    else{
            this->match_type_ = UNDEFINED;
        }
    }

    /*std::ofstream ofs ("overlapt.txt", std::ofstream::app);
    ofs <<  "===============================================\n"
    << "Read A id "<< std::setfill('0') << std::setw(5) <<this->read_A_id_
    << "\nRead B id "  << std::setfill('0') << std::setw(5) << this->read_B_id_
    << "\nRead A eff start "<< std::setfill('0') << std::setw(5)  << this->eff_read_A_start_
    << " Read A eff end "<< std::setfill('0') << std::setw(5)  << this->eff_read_A_end_
    << " Read A length " << std::setfill('0') << std::setw(5)  << this->alen
    << " Read A match start "<< std::setfill('0') << std::setw(5) <<  this->read_A_match_start_
    << " Read A eff match start " << std::setfill('0') << std::setw(5) <<  this->eff_read_A_match_start_
    << " Read A match end " << std::setfill('0') << std::setw(5)  << this->read_A_match_end_
    << " Read A eff match end " << std::setfill('0') << std::setw(5)  << this->eff_read_A_match_end_
    << "\nRead B eff start "  << std::setfill('0') << std::setw(5) << this->eff_read_B_start_
    << " Read B eff end " << std::setfill('0') << std::setw(5)  << this->eff_read_B_end_
    << " Read B length " << std::setfill('0') << std::setw(5)  << this->blen
    << " Read B match start "<< std::setfill('0') << std::setw(5) <<  this->read_B_match_start_
    << " Read B eff match start " << std::setfill('0') << std::setw(5) <<  this->eff_read_B_match_start_
    << " Read B match end " << std::setfill('0') << std::setw(5)  << this->read_B_match_end_
    << " Read B eff match end "  << std::setfill('0') << std::setw(5)  << this->eff_read_B_match_end_
    << "\nReverse complement "  << std::setfill('0') << std::setw(5)  << this->reverse_complement_match_
    << "\nMatch type "<<this->match_type_
    << "\n" << std::endl;
    ofs.close();*/
}