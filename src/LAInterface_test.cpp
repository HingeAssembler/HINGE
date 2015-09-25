

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
	la.showRead(0,1); //show read [0,1)
	
	test_read = la.getRead(0);
	test_read->showRead();
	
	
	
	la.CloseDB();
	return 0;
}
