#
# @fileoverview Copyright (c) 2024-202x, by Stefano Gualandi, UniPv,
#               via Ferrata, 1, Pavia, Italy, 27100
#
# @author stefano.gualandi@unipv.it (Stefano Gualandi)
#
#

# Directory for my files
MYHOME          = ${PWD}
BIN             = ${MYHOME}/bin
INCLUDE         = ${MYHOME}/include
LIB             = ${MYHOME}/lib
SRC             = ${MYHOME}/src

OPTFLAG = -O3 -m64 -fPIC -ffast-math -DNDEBUG -Wall -DLINUX -Wall -std=c++14
LDFLAGS = -O3 -m64 -fPIC -DNDEBUG -lm -pthread

COMPILER    = clang++ ${OPTFLAG}
LINKER      = clang++ ${LDFLAGS}

# Directory for output files
OUT_DIR=bin lib

GUROBI_INC = /Library/gurobi1200/macos_universal2/include
GUROBI_LIB = -L/Library/gurobi1200/macos_universal2/lib/ -lgurobi120

LEMON_INC = /Users/gualandi/solvers/lemon-dev
LEMON_LIB = /Users/gualandi/solvers/lemon-dev/build/lemon/libemon.a

laptop: ${OUT_DIR}
	${COMPILER} -c -g -pg ${SRC}/ATSP_CLI.cpp -o ${LIB}/ATSP_CLI.o -I${INCLUDE} -I${GUROBI_INC} -I${LEMON_INC}
	${LINKER} -o ${BIN}/subtour ${LIB}/ATSP_CLI.o ${GUROBI_LIB} ${LEMON_LIB}

# Create subdirectory for output files (bin,lib)
MKDIR_P = mkdir -p

lib:
	${MKDIR_P} lib

bin:
	${MKDIR_P} bin

# Be careful!
clean::
	rm -f *.o
	rm -f ${LIB}/*.o
	rm -f *~
	rm -f ${SRC}/*~ ${INCLUDE}/*~
