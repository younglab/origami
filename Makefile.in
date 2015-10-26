BAMLIBS = -lbamtools
BAMINC = -I/usr/local/include/bamtools

all: bamexes

bamexes:
	cd src && g++ $(BAMINC) -c mapped-reads-merge.cpp
	cd src && g++ -o mapped-reads-merge mapped-reads-merge.o $(BAMLIBS)
