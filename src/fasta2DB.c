/************************************************************************************\
*                                                                                    *
* Copyright (c) 2014, Dr. Eugene W. Myers (EWM). All rights reserved.                *
*                                                                                    *
* Redistribution and use in source and binary forms, with or without modification,   *
* are permitted provided that the following conditions are met:                      *
*                                                                                    *
*  · Redistributions of source code must retain the above copyright notice, this     *
*    list of conditions and the following disclaimer.                                *
*                                                                                    *
*  · Redistributions in binary form must reproduce the above copyright notice, this  *
*    list of conditions and the following disclaimer in the documentation and/or     *
*    other materials provided with the distribution.                                 *
*                                                                                    *
*  · The name of EWM may not be used to endorse or promote products derived from     *
*    this software without specific prior written permission.                        *
*                                                                                    *
* THIS SOFTWARE IS PROVIDED BY EWM ”AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND       *
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL EWM BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  *
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY      *
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     *
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN  *
* IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                      *
*                                                                                    *
* For any issues regarding this software and its use, contact EWM at:                *
*                                                                                    *
*   Eugene W. Myers Jr.                                                              *
*   Bautzner Str. 122e                                                               *
*   01099 Dresden                                                                    *
*   GERMANY                                                                          *
*   Email: gene.myers@gmail.com                                                      *
*                                                                                    *
\************************************************************************************/

/*******************************************************************************************
 *
 *  Add .fasta files to a DB:
 *     Adds the given fasta files in the given order to <path>.db.  If the db does not exist
 *     then it is created.  All .fasta files added to a given data base must have the same
 *     header format and follow Pacbio's convention.  A file cannot be added twice and this
 *     is enforced.  The command either builds or appends to the .<path>.idx and .<path>.bps
 *     files, where the index file (.idx) contains information about each read and their offsets
 *     in the base-pair file (.bps) that holds the sequences where each base is compessed
 *     into 2-bits.  The two files are hidden by virtue of their names beginning with a '.'.
 *     <path>.db is effectively a stub file with given name that contains an ASCII listing
 *     of the files added to the DB and possibly the block partitioning for the DB if DBsplit
 *     has been called upon it.
 *
 *  Author:  Gene Myers
 *  Date  :  May 2013
 *  Modify:  DB upgrade: now *add to* or create a DB depending on whether it exists, read
 *             multiple .fasta files (no longer a stdin pipe).
 *  Date  :  April 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "DB.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-v] <path:string> ( -f<file> | <input:fasta> ... )";

static char number[128] =
    { 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };

typedef struct
  { int    argc;
    char **argv;
    FILE  *input;
    int    count;
    char  *name;
  } File_Iterator;

File_Iterator *init_file_iterator(int argc, char **argv, FILE *input, int first)
{ File_Iterator *it;

  it = Malloc(sizeof(File_Iterator),"Allocating file iterator");
  it->argc  = argc;
  it->argv  = argv;
  it->input = input;
  if (input == NULL)
    it->count = first;
  else
    { it->count = 1;
      rewind(input);
    }
  return (it);
}

int next_file(File_Iterator *it)
{ static char nbuffer[MAX_NAME+8];

  if (it->input == NULL)
    { if (it->count >= it->argc)
        return (0);
      it->name = it->argv[it->count++];
    }
  else
    { char *eol;

      if (fgets(nbuffer,MAX_NAME+8,it->input) == NULL)
        { if (feof(it->input))
            return (0);
          SYSTEM_ERROR;
        }
      if ((eol = index(nbuffer,'\n')) == NULL)
        { fprintf(stderr,"%s: Line %d in file list is longer than %d chars!\n",
                         Prog_Name,it->count,MAX_NAME+7);
          it->name = NULL;
        }
      *eol = '\0';
      it->count += 1;
      it->name  = nbuffer;
    }
  return (1);
}


int main(int argc, char *argv[])
{ FILE  *istub, *ostub;
  char  *dbname;
  char  *root, *pwd;

  FILE  *bases, *indx;
  int64  boff, ioff;

  int    ifiles, ofiles;
  char **flist;

  HITS_DB db;
  int     ureads;
  int64   offset;

  FILE   *IFILE;
  int     VERBOSE;

  //   Usage: [-v] <path:string> ( -f<file> | <input:fasta> ... )

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("fasta2DB")

    IFILE = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'f':
            IFILE = fopen(argv[i]+2,"r");
            if (IFILE == NULL)
              { fprintf(stderr,"%s: Cannot open file of inputs '%s'\n",Prog_Name,argv[i]+2);
                exit (1);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if ((IFILE == NULL && argc <= 2) || (IFILE != NULL && argc != 2))
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  //  Try to open DB file, if present then adding to DB, otherwise creating new DB.  Set up
  //  variables as follows:
  //    dbname = full name of db = <pwd>/<root>.db
  //    istub  = open db file (if adding) or NULL (if creating)
  //    ostub  = new image of db file (will overwrite old image at end)
  //    bases  = .bps file positioned for appending
  //    indx   = .idx file positioned for appending
  //    ureads = # of reads currently in db
  //    offset = offset in .bps at which to place next sequence
  //    ioff   = offset in .idx file to truncate to if command fails
  //    boff   = offset in .bps file to truncate to if command fails
  //    ifiles = # of .fasta files to add
  //    ofiles = # of .fasta files already in db
  //    flist  = [0..ifiles+ofiles] list of file names (root only) added to db so far

  { int     i;

    root   = Root(argv[1],".db");
    pwd    = PathTo(argv[1]);
    dbname = Strdup(Catenate(pwd,"/",root,".db"),"Allocating db name");
    if (dbname == NULL)
      exit (1);

    if (IFILE == NULL)
      ifiles = argc-2;
    else
      { File_Iterator *ng;

        ifiles = 0;
        ng = init_file_iterator(argc,argv,IFILE,2);
        while (next_file(ng))
          ifiles += 1;
        free(ng);
      }

    istub = fopen(dbname,"r");
    if (istub == NULL)
      { ofiles = 0;

        bases = Fopen(Catenate(pwd,PATHSEP,root,".bps"),"w+");
        indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"w+");
        if (bases == NULL || indx == NULL)
          exit (1);

        fwrite(&db,sizeof(HITS_DB),1,indx);

        ureads  = 0;
        offset  = 0;
        boff    = 0;
        ioff    = 0;
      }
    else
      { if (fscanf(istub,DB_NFILE,&ofiles) != 1)
          SYSTEM_ERROR

        bases = Fopen(Catenate(pwd,PATHSEP,root,".bps"),"r+");
        indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r+");
        if (bases == NULL || indx == NULL)
          exit (1);

        if (fread(&db,sizeof(HITS_DB),1,indx) != 1)
          SYSTEM_ERROR
        fseeko(bases,0,SEEK_END);
        fseeko(indx, 0,SEEK_END);

        ureads = db.ureads;
        offset = ftello(bases);
        boff   = offset;
        ioff   = ftello(indx);
      }

    flist  = (char **) Malloc(sizeof(char *)*(ofiles+ifiles),"Allocating file list");
    ostub  = Fopen(Catenate(pwd,"/",root,".dbx"),"w+");
    if (ostub == NULL || flist == NULL)
      exit (1);

    fprintf(ostub,DB_NFILE,ofiles+ifiles);
    for (i = 0; i < ofiles; i++)
      { int  last;
        char prolog[MAX_NAME], fname[MAX_NAME];

        if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
          SYSTEM_ERROR
        if ((flist[i] = Strdup(fname,"Adding to file list")) == NULL)
          goto error;
        fprintf(ostub,DB_FDATA,last,fname,prolog);
      }
  }

  { int            maxlen;
    int64          totlen, count[4];
    int            pmax, rmax;
    HITS_READ     *prec;
    char          *read;
    int            c;
    File_Iterator *ng;

    //  Buffer for reads all in the same well

    pmax = 100;
    prec = (HITS_READ *) Malloc(sizeof(HITS_READ)*pmax,"Allocating record buffer");
    if (prec == NULL)
      goto error;

    //  Buffer for accumulating .fasta sequence over multiple lines

    rmax  = MAX_NAME + 60000;
    read  = (char *) Malloc(rmax+1,"Allocating line buffer");
    if (read == NULL)
      goto error;

    totlen = 0;              //  total # of bases in new .fasta files
    maxlen = 0;              //  longest read in new .fasta files
    for (c = 0; c < 4; c++)  //  count of acgt in new .fasta files
      count[c] = 0;

    //  For each new .fasta file do:

    ng = init_file_iterator(argc,argv,IFILE,2);
    while (next_file(ng))
      { FILE *input;
        char *path, *core, *prolog;
        int   nline, eof, rlen, pcnt;
        int   pwell;

        if (ng->name == NULL) goto error;

        //  Open it: <path>/<core>.fasta, check that core is not too long,
        //           and checking that it is not already in flist.

        path  = PathTo(ng->name);
        core  = Root(ng->name,".fasta");
        if ((input = Fopen(Catenate(path,"/",core,".fasta"),"r")) == NULL)
          goto error;
        free(path);
        if (strlen(core) >= MAX_NAME)
          { fprintf(stderr,"%s: File name over %d chars: '%.200s'\n",
                           Prog_Name,MAX_NAME,core);
            goto error;
          }

        { int j;

          for (j = 0; j < ofiles; j++)
            if (strcmp(core,flist[j]) == 0)
              { fprintf(stderr,"%s: File %s.fasta is already in database %s.db\n",
                               Prog_Name,core,Root(argv[1],".db"));
                goto error;
              }
        }

        //  Get the header of the first line.  If the file is empty skip.

        pcnt  = 0;
        rlen  = 0;
        nline = 1;
        eof   = (fgets(read,MAX_NAME,input) == NULL);
        if (eof || strlen(read) < 1)
          { fprintf(stderr,"Skipping '%s', file is empty!\n",core);
            fclose(input);
            free(core);
            continue;
          }

        //   Add the file name to flist

        if (VERBOSE)
          { fprintf(stderr,"Adding '%s' ...\n",core);
            fflush(stderr);
          }
        flist[ofiles++] = core;

        // Check that the first line has PACBIO format, and record prolog in 'prolog'.

        if (read[strlen(read)-1] != '\n')
          { fprintf(stderr,"File %s.fasta, Line 1: Fasta line is too long (> %d chars)\n",
                           core,MAX_NAME-2);
            goto error;
          }
        if (!eof && read[0] != '>')
          { fprintf(stderr,"File %s.fasta, Line 1: First header in fasta file is missing\n",core);
            goto error;
          }

        { char *find;
          int   well, beg, end, qv;

          find = index(read+1,'/');
          if (find != NULL && sscanf(find+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv) >= 3)
            { *find = '\0';
              prolog = Strdup(read+1,"Extracting prolog");
              *find = '/';
              if (prolog == NULL) goto error;
            }
          else
            { fprintf(stderr,"File %s.fasta, Line %d: Pacbio header line format error\n",
                             core,nline);
              goto error;
            }
        }

        //  Read in all the sequences until end-of-file

        { int i, x;

          pwell = -1;
          while (!eof)
            { int   beg, end, clen;
              int   well, qv;
              char *find;

              find = index(read+(rlen+1),'/');
              if (find == NULL)
                { fprintf(stderr,"File %s.fasta, Line %d: Pacbio header line format error\n",
                                 core,nline);
                  goto error;
                }
              *find = '\0';
              if (strcmp(read+(rlen+1),prolog) != 0)
                { fprintf(stderr,"File %s.fasta, Line %d: Pacbio header line name inconsisten\n",
                                 core,nline);
                  goto error;
                }
              *find = '/';
              x = sscanf(find+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv);
              if (x < 3)
                { fprintf(stderr,"File %s.fasta, Line %d: Pacbio header line format error\n",
                                 core,nline);
                  goto error;
                }
              else if (x == 3)
                qv = 0;

              rlen  = 0;
              while (1)
                { eof = (fgets(read+rlen,MAX_NAME,input) == NULL);
                  nline += 1;
                  x = strlen(read+rlen)-1;
                  if (read[rlen+x] != '\n')
                    { fprintf(stderr,"File %s.fasta, Line %d:",core,nline);
                      fprintf(stderr," Fasta line is too long (> %d chars)\n",MAX_NAME-2);
                      goto error;
                    }
                  if (eof || read[rlen] == '>')
                    break;
                  rlen += x;
                  if (rlen + MAX_NAME > rmax)
                    { rmax = ((int) (1.2 * rmax)) + 1000 + MAX_NAME;
                      read = (char *) realloc(read,rmax+1);
                      if (read == NULL)
                        { fprintf(stderr,"File %s.fasta, Line %d:",core,nline);
                          fprintf(stderr," Out of memory (Allocating line buffer)\n");
                          goto error;
                        }
                    }
                }
              read[rlen] = '\0';

              for (i = 0; i < rlen; i++)
                { x = number[(int) read[i]];
                  count[x] += 1;
                  read[i]   = (char) x;
                }
              ureads += 1;
              totlen += rlen;
              if (rlen > maxlen)
                maxlen = rlen;

              prec[pcnt].origin = well;
              prec[pcnt].fpulse = beg;
              prec[pcnt].rlen   = rlen;
              prec[pcnt].boff   = offset;
              prec[pcnt].coff   = -1;
              prec[pcnt].flags  = qv;

              Compress_Read(rlen,read);
              clen = COMPRESSED_LEN(rlen);
              fwrite(read,1,clen,bases);
              offset += clen;

              if (pwell == well)
                { prec[pcnt].flags |= DB_CSS;
                  pcnt += 1;
                  if (pcnt >= pmax)
                    { pmax = ((int) (pcnt*1.2)) + 100;
                      prec = (HITS_READ *) realloc(prec,sizeof(HITS_READ)*pmax);
                      if (prec == NULL)
                        { fprintf(stderr,"File %s.fasta, Line %d: Out of memory",core,nline);
                          fprintf(stderr," (Allocating read records)\n");
                          goto error;
                        }
                    }
                }
              else if (pcnt == 0)
                pcnt += 1;
              else
                { x = 0;
                  for (i = 1; i < pcnt; i++)
                    if (prec[i].rlen > prec[x].rlen)
                      x = i;
                  prec[x].flags |= DB_BEST;
                  fwrite(prec,sizeof(HITS_READ),pcnt,indx);
                  prec[0] = prec[pcnt];
                  pcnt = 1;
                }
              pwell = well;
            }

          //  Complete processing of .fasta file: flush last well group, write file line
          //      in db image, free prolog, and close file

          x = 0;
          for (i = 1; i < pcnt; i++)
            if (prec[i].rlen > prec[x].rlen)
              x = i;
          prec[x].flags |= DB_BEST;
          fwrite(prec,sizeof(HITS_READ),pcnt,indx);

          fprintf(ostub,DB_FDATA,ureads,core,prolog);
        }

        free(prolog);
        fclose(input);
      }

    //  Finished loading all sequences: update relevant fields in db record

    db.ureads = ureads;
    if (istub == NULL)
      { for (c = 0; c < 4; c++)
          db.freq[c] = (float) ((1.*count[c])/totlen);
        db.totlen = totlen;
        db.maxlen = maxlen;
        db.cutoff = -1;
      }
    else
      { for (c = 0; c < 4; c++)
          db.freq[c] = (float) ((db.freq[c]*db.totlen + (1.*count[c]))/(db.totlen + totlen));
        db.totlen += totlen;
        if (maxlen > db.maxlen)
          db.maxlen = maxlen;
      }
  }

  //  If db has been previously partitioned then calculate additional partition points and
  //    write to new db file image

  if (db.cutoff >= 0)
    { int64      totlen, dbpos, size;
      int        nblock, ireads, tfirst, rlen;
      int        ufirst, cutoff, allflag;
      HITS_READ  record;
      int        i;

      if (VERBOSE)
        { fprintf(stderr,"Updating block partition ...\n");
          fflush(stderr);
        }

      //  Read the block portion of the existing db image getting the indices of the first
      //    read in the last block of the exisiting db as well as the partition parameters.
      //    Copy the old image block information to the new block information (except for
      //    the indices of the last partial block)

      if (fscanf(istub,DB_NBLOCK,&nblock) != 1)
        SYSTEM_ERROR
      dbpos = ftello(ostub);
      fprintf(ostub,DB_NBLOCK,0);
      if (fscanf(istub,DB_PARAMS,&size,&cutoff,&allflag) != 3)
        SYSTEM_ERROR
      fprintf(ostub,DB_PARAMS,size,cutoff,allflag); 
      if (allflag)
        allflag = 0;
      else
        allflag = DB_BEST;
      size *= 1000000ll;

      nblock -= 1;
      for (i = 0; i <= nblock; i++)
        { if (fscanf(istub,DB_BDATA,&ufirst,&tfirst) != 2)
            SYSTEM_ERROR
          fprintf(ostub,DB_BDATA,ufirst,tfirst);
        }

      //  Seek the first record of the last block of the existing db in .idx, and then
      //    compute and record partition indices for the rest of the db from this point
      //    forward.

      fseeko(indx,sizeof(HITS_DB)+sizeof(HITS_READ)*ufirst,SEEK_SET);
      totlen = 0;
      ireads = 0;
      for (i = ufirst; i < ureads; i++)
        { if (fread(&record,sizeof(HITS_READ),1,indx) != 1)
            SYSTEM_ERROR
          rlen = record.rlen;
          if (rlen >= cutoff && (record.flags & DB_BEST) >= allflag)
            { ireads += 1;
              tfirst += 1;
              totlen += rlen;
              if (totlen >= size)
                { fprintf(ostub," %9d %9d\n",i+1,tfirst);
                  totlen = 0;
                  ireads = 0;
                  nblock += 1;
                }
            }
        }

      if (ireads > 0)
        { fprintf(ostub,DB_BDATA,ureads,tfirst);
          nblock += 1;
        }

      db.treads = tfirst;

      fseeko(ostub,dbpos,SEEK_SET);
      fprintf(ostub,DB_NBLOCK,nblock);    //  Rewind and record the new number of blocks
    }
  else
    db.treads = ureads;

  rewind(indx);
  fwrite(&db,sizeof(HITS_DB),1,indx);   //  Write the finalized db record into .idx

  rewind(ostub);                        //  Rewrite the number of files actually added
  fprintf(ostub,DB_NFILE,ofiles);

  if (istub != NULL)
    fclose(istub);
  fclose(ostub);
  fclose(indx);
  fclose(bases);

  rename(Catenate(pwd,"/",root,".dbx"),dbname);   //  New image replaces old image

  exit (0);

  //  Error exit:  Either truncate or remove the .idx and .bps files as appropriate.
  //               Remove the new image file <pwd>/<root>.dbx

error:
  if (ioff != 0)
    { fseeko(indx,0,SEEK_SET);
      if (ftruncate(fileno(indx),ioff) < 0)
        SYSTEM_ERROR
    }
  if (boff != 0)
    { fseeko(bases,0,SEEK_SET);
      if (ftruncate(fileno(bases),boff) < 0)
        SYSTEM_ERROR
    }
  fclose(indx);
  fclose(bases);
  if (ioff == 0)
    unlink(Catenate(pwd,PATHSEP,root,".idx"));
  if (boff == 0)
    unlink(Catenate(pwd,PATHSEP,root,".bps"));

  if (istub != NULL)
    fclose(istub);
  fclose(ostub);
  unlink(Catenate(pwd,"/",root,".dbx"));

  exit (1);
}
