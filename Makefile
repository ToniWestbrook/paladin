CC=			gcc
#CC=			clang --analyze
CFLAGS=		-g -Wall -Wno-unused-function -Wno-unused-result -O2 
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)
LOBJS=		utils.o kthread.o kstring.o ksw.o bwt.o bntseq.o bwa.o bwamem.o bwamem_pair.o bwamem_extra.o malloc_wrap.o
AOBJS=		is.o bwtindex.o kopen.o align.o translations.o protein.o uniprot.o bwashm.o
PROG=		paladin
INCLUDES=	
LIBS=		-lm -lz -lpthread
SUBDIRS=	.

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

paladin:libbwa.a $(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) main.o -o $@ -L. -lbwa -lcurl $(LIBS)

libbwa.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.

bntseq.o: bntseq.h utils.h kseq.h malloc_wrap.h khash.h
bwa.o: bntseq.h bwa.h bwt.h ksw.h utils.h kstring.h malloc_wrap.h kseq.h
bwamem.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h ksw.h kvec.h
bwamem.o: ksort.h utils.h kbtree.h
bwamem_extra.o: bwa.h bntseq.h bwt.h bwamem.h kstring.h malloc_wrap.h
bwamem_pair.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h kvec.h
bwamem_pair.o: utils.h ksw.h
bwashm.o: bwa.h bntseq.h bwt.h
bwt.o: utils.h bwt.h kvec.h malloc_wrap.h
bwtindex.o: bntseq.h bwt.h utils.h malloc_wrap.h
align.o: bwa.h bntseq.h bwt.h bwamem.h kvec.h malloc_wrap.h utils.h kseq.h
is.o: malloc_wrap.h
kopen.o: malloc_wrap.h
kstring.o: kstring.h malloc_wrap.h
ksw.o: ksw.h malloc_wrap.h
main.o: main.h kstring.h malloc_wrap.h utils.h
malloc_wrap.o: malloc_wrap.h
protein.o: protein.h translations.h utils.h kseq.h malloc_wrap.h khash.h
translations.o: translations.h 
utils.o: utils.h ksort.h malloc_wrap.h kseq.h
uniprot.o: uniprot.h
