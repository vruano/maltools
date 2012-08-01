.PHONY: all clean mrproper

CPP ?= g++
CFLAGS ?= -L/opt/local/lib
OBJDIR ?= .
BINDIR ?= .
SRCDIR ?= .

PROGS ?= hmm bwa2som bam2snps

all: $(addprefix $(BINDIR)/,$(PROGS))

clean:
	rm -f $(OBJDIR)/hmm.o

mrproper: clean
	rm -f $(BINDIR)/hmm rm -f $(BINDIR)

$(OBJDIR)/hmm.o : $(SRCDIR)/hmm.cpp $(SRCDIR)/hmm.h
	$(CPP) $(CFLAGS) `pkg-config --cflags gsl` -I $(SRCDIR) -c $< -o $@

$(BINDIR)/hmm: $(OBJDIR)/hmm.o
	$(CPP) $(CFLAGS) $< `pkg-config --libs gsl` -o $@   

$(BINDIR)/bwa2som: $(SRCDIR)/bwa2som.cpp
	$(CPP) $(CFLAGS) $< -O3 -o $@ -lpng -L $(SAMTOOLS_DIR)/samtools -I $(SAMTOOLS_DIR)/samtools -lbam -lz

$(BINDIR)/bam2snps: $(SRCDIR)/bam2snps.cpp
	$(CPP) $(CFLAGS) $< -O3 -o $@ 
