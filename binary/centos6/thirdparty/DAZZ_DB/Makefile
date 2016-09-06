CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = fasta2DB DB2fasta quiva2DB DB2quiva DBsplit DBdust Catrack DBshow DBstats DBrm simulator \
      fasta2DAM DAM2fasta DBdump

all: $(ALL)

fasta2DB: fasta2DB.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o fasta2DB fasta2DB.c DB.c QV.c -lm

DB2fasta: DB2fasta.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DB2fasta DB2fasta.c DB.c QV.c -lm

quiva2DB: quiva2DB.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o quiva2DB quiva2DB.c DB.c QV.c -lm

DB2quiva: DB2quiva.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DB2quiva DB2quiva.c DB.c QV.c -lm

DBsplit: DBsplit.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBsplit DBsplit.c DB.c QV.c -lm

DBdust: DBdust.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBdust DBdust.c DB.c QV.c -lm

Catrack: Catrack.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o Catrack Catrack.c DB.c QV.c -lm

DBshow: DBshow.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBshow DBshow.c DB.c QV.c -lm

DBdump: DBdump.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBdump DBdump.c DB.c QV.c -lm

DBstats: DBstats.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBstats DBstats.c DB.c QV.c -lm

DBrm: DBrm.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBrm DBrm.c DB.c QV.c -lm

simulator: simulator.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o simulator simulator.c DB.c QV.c -lm

fasta2DAM: fasta2DAM.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o fasta2DAM fasta2DAM.c DB.c QV.c -lm

DAM2fasta: DAM2fasta.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DAM2fasta DAM2fasta.c DB.c QV.c -lm

DBupgrade.Sep.25.2014: DBupgrade.Sep.25.2014.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBupgrade.Sep.25.2014 DBupgrade.Sep.25.2014.c DB.c QV.c -lm

DBupgrade.Dec.31.2014: DBupgrade.Dec.31.2014.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBupgrade.Dec.31.2014 DBupgrade.Dec.31.2014.c DB.c QV.c -lm

DUSTupgrade.Jan.1.2015: DUSTupgrade.Jan.1.2015.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DUSTupgrade.Jan.1.2015 DUSTupgrade.Jan.1.2015.c DB.c QV.c -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f DBupgrade.Sep.25.2014 DBupgrade.Dec.31.2014 DUSTupgrade.Jan.1.2015
	rm -f dazz.db.tar.gz

install:
	cp $(ALL) ~/bin

package:
	make clean
	tar -zcf dazz.db.tar.gz README Makefile *.h *.c
