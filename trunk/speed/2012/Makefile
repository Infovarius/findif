## The default is to compile source files and build the program.

PROGRAM_NAME = tor
DIR = ./
FILES = $(DIR)main.c $(DIR)pde.c $(DIR)util.c $(DIR)rkad.c
IO = $(DIR)ioutil.c
# IO = $(DIR)iobinutil.c
COMPILER = mpicc
# COMPILER = gcc
KEYS = -lm -O3

all:
	$(MAKE) compile
	$(MAKE) link

compile:
	@echo "Compiling all source files..."
	for fname in $(FILES); do 	\
		$(COMPILER) -c $$f;	\
	done;				\
	$(COMPILER) -c $(IO)

#debug:
#	dtcc $(FILES)

link:
	@echo "Building program..."
	$(COMPILER) $(FILES) $(IO) -o $(PROGRAM_NAME) $(KEYS)

clean:
	@echo "Deleting all compiled files..."
	@rm *.o -f
	@rm *.log -f
	@rm *.err -f
	@rm *.gch -f

uninst un uninstall:	
	@$(MAKE) clean
	@echo "Deleting all results..."
	@rm $(PROGRAM_NAME) -f
	@rm *.dat -f
	@rm *.sta -f
	@rm node -f
	@rm coord -f
	@rm *.cp -f

start:
	@$(MAKE) clean
	@echo "Initialisating calculation..."
	@rm *.cp -f
	@rm *.dat -f
	@rm coord -f
	@rm node -f
	@rm *.sta -f
