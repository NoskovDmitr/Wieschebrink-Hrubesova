#
# Galois Field Arithmetic Library
# By Arash Partow - 2000
#
# URL: http://www.partow.net/projects/galois/index.html
#
# Copyright Notice:
# Free use of this library is permitted under the
# guidelines and in accordance with the most
# current version of the Common Public License.
# http://www.opensource.org/licenses/cpl.php
#

COMPILER         = g++ 
OPTIMIZATION_OPT = -O3
OPTIONS          = -g -pedantic -ansi -Wall -lstdc++ $(OPTIMIZATION_OPT) -o 
OPTIONS_LIBS     = -g -pedantic -ansi -Wall $(OPTIMIZATION_OPT) -c 


CPP_SRC = GaloisField.cpp \
	grsmatrix.cpp \
	matrix.cpp \
	permmatrix.cpp \
	randmatrix.cpp \
	regmatrix.cpp \
	attacker.cpp \
	cryptosystem.cpp \
	polynomial.cpp \
	testin.cpp \
	dimtests.cpp \
	timetests.cpp


OBJECTS = $(CPP_SRC:.cpp=.o)


%.o: %.hpp %.cpp
	$(COMPILER) $(OPTIONS_LIBS) $*.cpp


all: $(OBJECTS) dimensionTests crashTests

crashTests: crashTests.cpp $(OBJECTS)
	$(COMPILER) $(OPTIONS) crashTests crashTests.cpp $(OBJECTS)

dimensionTests: dimensionTests.cpp $(OBJECTS)
	$(COMPILER) $(OPTIONS) dimensionTests dimensionTests.cpp $(OBJECTS)

clean:
	rm -f core *.o *.bak *stackdump *#
