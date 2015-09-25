
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



void Read::showRead(){
	std::cout<<"read #"<<id<<std::endl;
	std::cout<<">" << name << std::endl;
	std::cout<< bases << std::endl;
}


int LAInterface::OpenDB2(std::string filename) {
	char * fn = new char[filename.length() + 1];
	strcpy(fn, filename.c_str());
		
    int status = Open_DB(fn,this->db1);
    if (status < 0)
      exit (1);
    if (this->db1->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n", "test" ,fn);
        exit (1);
      }
	  
    this->db2 = this->db1;
    Trim_DB(db1);
	
	delete [] fn;
	return 0;
}


int LAInterface::OpenDB(std::string filename) {
	char * fn = new char[filename.length() + 1];
	strcpy(fn, filename.c_str());
		
    int status = Open_DB(fn,this->db1);
    if (status < 0)
      exit (1);
    if (this->db1->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n", "test" ,fn);
        exit (1);
      }
	  
    //this->db2 = this->db1;
    //Trim_DB(db1);
	
	FILE* dstub;
	
	char * fn2 = new char[filename.length() + 1 + 3];
	strcpy(fn2,fn);
	strcat(fn2,".db");
	
    dstub  = Fopen(fn2,"r");
    if (dstub == NULL)
      exit (1);
    
    if (fscanf(dstub,DB_NFILE,&nfiles) != 1)
      SYSTEM_ERROR
		  
	printf("%d files\n",nfiles);

    flist = (char **) Malloc(sizeof(char *)*nfiles,"Allocating file list");
    findx = (int *) Malloc(sizeof(int *)*(nfiles+1),"Allocating file index");
    
	if (flist == NULL || findx == NULL)
      exit (1);

    findx += 1;
    findx[-1] = 0;

    for (int i = 0; i < nfiles; i++)
      { char prolog[MAX_NAME], fname[MAX_NAME];

        if (fscanf(dstub,DB_FDATA,findx+i,fname,prolog) != 3)
          SYSTEM_ERROR
        if ((flist[i] = Strdup(prolog,"Adding to file list")) == NULL)
          exit (1);
      }

    fclose(dstub);
	
	delete [] fn;
	
	
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

	//printf("1");

    if (flist == NULL || findx == NULL)
      exit (1);
	
	
      HITS_READ  *reads;
      HITS_TRACK *first;
      char       *read, **entry;
      int         c, b, e, i;
      int         hilight, substr;
      int         map;
      int       (*iscase)(int);

      read  = New_Read_Buffer(db1);
	  bool DOQVS = false;
	  bool DOSEQ = true;
	  int UPPER = 1;
	  int WIDTH = 80;
	  //printf("2");
      if (DOQVS)
        { entry = New_QV_Buffer(db1);
          first = db1->tracks->next;
        }
      else
        { entry = NULL;
          first = db1->tracks;
        }
	

      if (UPPER == 1)
        { hilight = 'A'-'a';
          iscase  = islower;
        }
      else
        { hilight = 'a'-'A';
          iscase  = isupper;
        }

      map    = 0;
      reads  = db1->reads;
      substr = 0;

      c = 0;
      //while (1)
      //  { 
            /*{ if (c >= reps)
                break;
              b = pts[c]-1;
              e = pts[c+1];
              if (e > db1->nreads)
                e = db1->nreads;
              c += 2;
            }*/
			
	  	b = from;
	  	e = to;

          for (i = b; i < e; i++)
            { int         len;
              int         fst, lst;
              int         flags, qv;
              HITS_READ  *r;
              HITS_TRACK *track;

              r   = reads + i;
              len = r->rlen;

              flags = r->flags;
              qv    = (flags & DB_QV);
              
                { while (i < findx[map-1])
                    map -= 1;
                  while (i >= findx[map])
                    map += 1;
                  //if (QUIVA)
                  //  printf("@%s/%d/%d_%d",flist[map],r->origin,r->fpulse,r->fpulse+len);
                  //else
                    printf(">%s/%d/%d_%d",flist[map],r->origin,r->fpulse,r->fpulse+len);
                  if (qv > 0)
                    printf(" RQ=0.%3d",qv);
                }
              printf("\n");

              if (DOQVS)
                Load_QVentry(db1,i,entry,UPPER);
              if (DOSEQ)
                Load_Read(db1,i,read,UPPER);

              for (track = first; track != NULL; track = track->next)
                { int64 *anno;
                  int   *data;
                  int64  s, f, j;
                  int    bd, ed, m;

                  anno = (int64 *) track->anno;
                  data = (int *) track->data;

                  s = (anno[i] >> 2);
                  f = (anno[i+1] >> 2);
                  if (s < f)
                    { for (j = s; j < f; j += 2)
                        { bd = data[j];
                          ed = data[j+1];
                          if (DOSEQ)
                            for (m = bd; m < ed; m++)
                              if (iscase(read[m]))
                                read[m] = (char) (read[m] + hilight);
                          if (j == s)
                            printf("> %s:",track->name);
                          printf(" [%d,%d]",bd,ed);
                        }
                      printf("\n");
                    }
                }

              //if (substr)
              //  { fst = iter->beg;
              //    lst = iter->end;
              //  }
              //else
                //{ 
					fst = 0;
                  	lst = len;
					//}

              /*if (QUIVA)
                { int k;

                  for (k = 0; k < 5; k++)
                    printf("%.*s\n",lst-fst,entry[k]+fst);
                }
					*/
              //else
                { if (DOQVS)
                    { int j, k;

                      printf("\n");
                      for (j = fst; j+WIDTH < lst; j += WIDTH)
                        { if (DOSEQ)
                            printf("%.*s\n",WIDTH,read+j);
                          for (k = 0; k < 5; k++)
                            printf("%.*s\n",WIDTH,entry[k]+j);
                          printf("\n");
                        }
                      if (j < lst)
                        { if (DOSEQ)
                            printf("%.*s\n",lst-j,read+j);
                          for (k = 0; k < 5; k++)
                            printf("%.*s\n",lst-j,entry[k]+j);
                          printf("\n");
                        }
                    }
                  else if (DOSEQ)
                    { int j;
    
                      for (j = fst; j+WIDTH < lst; j += WIDTH)
                        printf("%.*s\n",WIDTH,read+j);
                      if (j < lst)
                        printf("%.*s\n",lst-j,read+j);
                    }
                }
            }
			//}
	
}


Read * LAInterface::getRead(int number) {

	//printf("1");
	std::stringstream ss;
	std::string read_name;
	std::string read_bases;
	
    if (flist == NULL || findx == NULL)
      exit (1);
	
	
      HITS_READ  *reads;
      HITS_TRACK *first;
      char       *read, **entry;
      int         c, b, e, i;
      int         hilight, substr;
      int         map;
      int       (*iscase)(int);

      read  = New_Read_Buffer(db1);
	  bool DOQVS = false;
	  bool DOSEQ = true;
	  int UPPER = 1;
	  int WIDTH = 80;
	  //printf("2");
      if (DOQVS)
        { entry = New_QV_Buffer(db1);
          first = db1->tracks->next;
        }
      else
        { entry = NULL;
          first = db1->tracks;
        }
	

      if (UPPER == 1)
        { hilight = 'A'-'a';
          iscase  = islower;
        }
      else
        { hilight = 'a'-'A';
          iscase  = isupper;
        }

      map    = 0;
      reads  = db1->reads;
      substr = 0;

      c = 0;
      //while (1)
      //  { 
            /*{ if (c >= reps)
                break;
              b = pts[c]-1;
              e = pts[c+1];
              if (e > db1->nreads)
                e = db1->nreads;
              c += 2;
            }*/
			
	  	b = number;
	  	e = number + 1;

          for (i = b; i < e; i++)
            { int         len;
              int         fst, lst;
              int         flags, qv;
              HITS_READ  *r;
              HITS_TRACK *track;

              r   = reads + i;
              len = r->rlen;

              flags = r->flags;
              qv    = (flags & DB_QV);
              
                { while (i < findx[map-1])
                    map -= 1;
                  while (i >= findx[map])
                    map += 1;
                  //if (QUIVA)
                  //  printf("@%s/%d/%d_%d",flist[map],r->origin,r->fpulse,r->fpulse+len);
                  //else
                    //printf(">%s/%d/%d_%d",flist[map],r->origin,r->fpulse,r->fpulse+len);
				  ss << flist[map] << '/' << r->origin << '/' << r->fpulse << '_' << r->fpulse+len;
                  if (qv > 0)
                    //printf(" RQ=0.%3d",qv);
					  ss << "RQ=" << qv;
                }
				
				ss >> read_name;
              //printf("\n");

              if (DOQVS)
                Load_QVentry(db1,i,entry,UPPER);
              if (DOSEQ)
                Load_Read(db1,i,read,UPPER);

              for (track = first; track != NULL; track = track->next)
                { int64 *anno;
                  int   *data;
                  int64  s, f, j;
                  int    bd, ed, m;

                  anno = (int64 *) track->anno;
                  data = (int *) track->data;

                  s = (anno[i] >> 2);
                  f = (anno[i+1] >> 2);
                  if (s < f)
                    { for (j = s; j < f; j += 2)
                        { bd = data[j];
                          ed = data[j+1];
                          if (DOSEQ)
                            for (m = bd; m < ed; m++)
                              if (iscase(read[m]))
                                read[m] = (char) (read[m] + hilight);
                          if (j == s)
                            printf("> %s:",track->name);
                          printf(" [%d,%d]",bd,ed);
                        }
                      printf("\n");
                    }
                }

				read_bases = std::string(read);
              //if (substr)
              //  { fst = iter->beg;
              //    lst = iter->end;
              //  }
              //else
                //{ 
					fst = 0;
                  	lst = len;
					//}

              /*if (QUIVA)
                { int k;

                  for (k = 0; k < 5; k++)
                    printf("%.*s\n",lst-fst,entry[k]+fst);
                }
					*/
              //else
                /*{ if (DOQVS)
                    { int j, k;

                      printf("\n");
                      for (j = fst; j+WIDTH < lst; j += WIDTH)
                        { if (DOSEQ)
                            printf("%.*s\n",WIDTH,read+j);
                          for (k = 0; k < 5; k++)
                            printf("%.*s\n",WIDTH,entry[k]+j);
                          printf("\n");
                        }
                      if (j < lst)
                        { if (DOSEQ)
                            printf("%.*s\n",lst-j,read+j);
                          for (k = 0; k < 5; k++)
                            printf("%.*s\n",lst-j,entry[k]+j);
                          printf("\n");
                        }
                    }
                  else if (DOSEQ)
                    { int j;
    
                      for (j = fst; j+WIDTH < lst; j += WIDTH)
                        printf("%.*s\n",WIDTH,read+j);
                      if (j < lst)
                        printf("%.*s\n",lst-j,read+j);
                    }
                }*/
					
					
            }
			//}
			
			Read * new_r = new Read(number, read_name, read_bases);
			return new_r;
}