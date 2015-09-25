
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "LAInterface.h"


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
	  
    this->db2 = this->db1;
    Trim_DB(db1);
	
	delete [] fn;
	return 0;
}
