# Makefile for KabatMan
# Copyright (c) 1994, Dr. Andrew C.R. Martin
#
# N.B. The ansi program used to gerenate prototype files from C source
# may be obtained from me on request
#
CC     = cc -O3 -L$(HOME)/lib -I$(HOME)/include
COPT   = 
LIBS   = -lm -lbiop -lgen -lm -lxml2
ANSI   = ansi -p
EXE    = kabatman
OFILES = kabatman.o RdKabat.o BuildSelect.o BuildWhere.o \
         ExecSearch.o KabCho.o subgroup.o
PFILES = kabatman.p RdKabat.p BuildSelect.p BuildWhere.p \
         ExecSearch.p KabCho.p subgroup.p


all    : $(EXE) splitkabat


splitkabat : splitkabat.o
	$(CC) $(COPT) -o splitkabat splitkabat.o


$(EXE) : $(OFILES)
	$(CC) $(COPT) -o $(EXE) $(OFILES) $(LIBS)

.c.o   :
	$(CC) $(COPT) -DNOBIOPLIB -o $@ -c $<

protos : $(PFILES)

.c.p   :
	$(ANSI) $< $@

clean  :
	/bin/rm -f $(OFILES) splitkabat.o
