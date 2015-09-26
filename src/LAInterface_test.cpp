

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

static int ORDER(const void *l, const void *r)
{ int x = *((int32 *) l);
  int y = *((int32 *) r);
  return (x-y);
}

int main(int argc, char *argv[])
{
	LAInterface la;
	std::cout<<"hello"<<std::endl;
	Read * test_read;
	
	la.OpenDB("G");
	la.showRead(1,3); //show read [1,3)
	
	test_read = la.getRead(0); //get read 0
	test_read->showRead(); // show read 0 
	
	la.OpenAlignment("G.1.las");
	la.showAlignment(0,2); // show alignments of read [0,2) 
	la.OpenAlignment("G.1.las");
	
	std::vector<int> res;
	la.getAlignment(res,1);
	for (auto i:res)
		printf("%d ",i);
	//tbd
	
	la.CloseDB();
	return 0;
}
