# GRR20176143 Cláudio Torres Júnior
# GRR20171607 Gabriela Stein
# -----------------------------------------------------------------------------

	CC = gcc -std=c11 -g -O3 -mavx -march=native
    CFLAGS = $(LIKWID_FLAGS) -Wall
	LFLAGS = $(LIKWID_LIBS) -lm

    PROG = cgSolver
    OBJS = utils.o \
             sistemarandom.o \
             gradienteconjugado.o \
             $(PROG).o

    MODULOS   = matriz utils

       LIKWID = /home/soft/likwid
 LIKWID_FLAGS = -DLIKWID_PERFMON -I$(LIKWID)/include
  LIKWID_LIBS = -L$(LIKWID)/lib -llikwid

.PHONY: all clean limpa purge faxina distclean debug avx likwid

%.o: %.c %.h utils.h
	$(CC) $(CFLAGS) -c $<

all: $(PROG)

debug:     CFLAGS += -DDEBUG

likwid:    CFLAGS += -DLIKWID_PERFMON
likwid:    LFLAGS += -llikwid

likwid debug: $(PROG)

$(PROG):  $(PROG).o

$(PROG): $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean:
	@echo "Limpando ...."
	@rm -f *~ *.bak *.tmp

purge distclean:   clean
	@echo "Faxina ...."
	@rm -f  $(PROG) *.o core a.out
	@rm -f marker.out *.log