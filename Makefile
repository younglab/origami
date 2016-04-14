BAMLIBS = -lbamtools
BAMINC = -I/usr/local/include/bamtools

all: bamexes

bamexes:
	cd src && g++ $(BAMINC) -c mapped-reads-merge.cpp
	cd src && g++ $(BAMINC) -c bam-read-extension.cpp
	cd src && g++ -o mapped-reads-merge mapped-reads-merge.o $(BAMLIBS)
	cd src && g++ -o bam-read-extension bam-read-extension.o $(BAMLIBS)
