# (C) Makefile for program #1             ATMS 502/CSE 566, Spring 2017
#
# A "makefile" is a series of directions needed to build a program.
# It is used by the Linux command "make".
#
# A makefile typically defines:
#   1. The compiler (that which turns code, like pgm1.c, into an intermediate object file, like pgm1.o)
#   2. Compiler options - any special (not-standard) settings telling the compiler what to do
#   3. The object files (".o" versions of source files ending in ".c" if fortran, or ".c" if C)
#      The compiler makes object files from your source code files, and once all of those are
#      up to date, it combines object files into a running program.
#
# Lines starting with "#" are comments - like this one!
# Type "man make" for more information on make, or look at these sites:
#
#	http://mrbook.org/blog/tutorials/make/
#	https://www.cs.usask.ca/staff/oster/makefiles.html
#	http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
#       http://www.cprogramming.com/tutorial/makefiles_continued.html
#
#  In this case: typing -
#    "make"          gets you an annoying "help" message.
#    "make pgm1"     actually does what you want: builds running program "pgm1"
#    "make clean"    deletes all the object (something.o), graphics and program files 
#    "make listing"  creates a text listing of all your code, with line numbers
#    "make archive"  creates a 'tape archive' (pgm1.tar) file with copies of all your files.
#
#  Note: you Only have to "make clean" if you have changed your OPTIONS.
#  In regular use, just type "make pgm1" to compile/build your program;
#  make will figure out which source code files you have changed, compiling
#  only those files, and then build your running program.

PROGRAM	= pgm6
OBJECTS = pgm6.o advection.o bc.o ic.o  stats.o  sfc.o advect1d.o contr.o error.o diffusion.o pgf.o update.o putfield.o
SOURCE  = pgm6.c advection.c bc.c ic.c  stats.c  sfc.c advect1d.c contr.c error.c diffusion.c pgf.c	update.c putfield.c
ARCHIVE = pgm6.tar
CC	= ncargcc
OPTIONS	= -openmp  -fno-alias  -fno-fnalias  -fargument-noalias  -restrict 
# use following for debugging (put "#" in front of OPTIONS above, remove "#" before OPTIONS below)
#OPTIONS = -g -debug extended -traceback -O0
# use following for extensive debugging
#OPTIONS = -g -debug extended -traceback -O0 -Wuninitialized -Wcheck -check=uninit -Wmissing-prototypes -ftrapuv -fp-stack-check

help:
	@echo Try:
	@echo make $(PROGRAM) .... to build the program named $(PROGRAM)
	@echo make clean .... to clean up, removing object files and program $(PROGRAM)
	@echo make listing .... to make a printable listing ... you may want to adjust the page length
	@echo make archive .... to make an archive tar file you can transfer or submit

$(PROGRAM):	$(OBJECTS)
	$(CC) $(OPTIONS) -o $(PROGRAM) $(OBJECTS)

%.o:	%.c
	$(CC) -c $(OPTIONS) $<

clean:	
	rm -f $(OBJECTS) $(PROGRAM) gmeta gmeta.zip *.gif

listing:
	@echo Creating code listing named listing.txt ...
	pr -F --length=59 --page-width=72 -n    \
           $(SOURCE) Makefile > listing.txt
	@echo listing.txt is `cat listing.txt | wc -l` lines long.

archive:
	@echo Creating code archive tar-file $(ARCHIVE) ...
	tar cf $(ARCHIVE) $(SOURCE) Makefile
	@ls -l $(ARCHIVE)
