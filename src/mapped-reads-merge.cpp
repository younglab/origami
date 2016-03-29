#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace BamTools;

namespace {
  int calc_rounded_perc(double num,double denom) {
    return(std::floor(num/denom*100));
  }
}

int main(int argc,char **argv) {
  
    if( argc < 5) {
      std::cerr << "Error: not enough arguments" << std::endl;
      std::cerr << "mapped-reads-merge <first read BAM> <second read BAM> <output BAM> <logfile>" << std::endl;
      return EXIT_FAILURE;
    }

  
    const char *leftfile = argv[1];
    const char *rightfile = argv[2];
    const char *outputfile = argv[3];
    const char *reportfile = argv[4];
    BamReader left, right;
    
    left.Open(leftfile);
    right.Open(rightfile);
    
    const SamHeader header = left.GetHeader();
    const RefVector references = left.GetReferenceData();
    
    BamWriter writer;
    writer.Open(outputfile,header,references);
    
    BamAlignment l, r;
    unsigned long totalpairs = 0, leftaligned = 0, rightaligned = 0, leftnotaligned = 0, rightnotaligned = 0,
      pairmapped=0,singleton=0,neither=0;
    while( left.GetNextAlignment(l) && right.GetNextAlignment(r)) {
        totalpairs++;
      
        if(l.IsMapped()) {
          leftaligned++;
        } else {
          leftnotaligned++;
        }
        
        if(r.IsMapped()) {
          rightaligned++;
        } else {
          rightnotaligned++;
        }
        
        if(l.IsMapped() && r.IsMapped()) {
          pairmapped++;
        } else if( l.IsMapped() || r.IsMapped()) {
          singleton++;
        } else {
          neither++;
        }
      
        l.MateRefID = r.RefID;
        l.MatePosition = r.Position;
        l.SetIsPaired(true);
        l.SetIsMateMapped(r.IsMapped());
        l.SetIsMateReverseStrand(r.IsReverseStrand());
        l.SetIsFirstMate(true);
        l.SetIsProperPair(true);
        writer.SaveAlignment(l);
        
        r.MateRefID = l.RefID;
        r.MatePosition = l.Position;
        r.SetIsPaired(true);
        r.SetIsMateMapped(l.IsMapped());
        r.SetIsMateReverseStrand(l.IsReverseStrand());
        r.SetIsSecondMate(true);
        r.SetIsProperPair(true);
        writer.SaveAlignment(r);
    }
    
    writer.Close();
    left.Close();
    right.Close();
    
    std::ofstream logfile;
    logfile.open(reportfile);
    
    logfile << "Alignment statistics: " << std::endl;
    logfile << "Total Read Pairs: " << totalpairs << std::endl;
    logfile << "Total Mapping Reads: " << (leftaligned+rightaligned) << " (" << calc_rounded_perc(leftaligned+rightaligned,2*totalpairs) << "%)" << std::endl;
    logfile << "Both Pairs Mapped: " << (pairmapped) << " (" << calc_rounded_perc(pairmapped,totalpairs) << "%)" << std::endl;
    logfile << "Singletons: " << singleton << " (" << calc_rounded_perc(singleton,totalpairs) << "%)" << std::endl;
    logfile << "Unmapped Pairs: " << neither << " (" << calc_rounded_perc(neither,totalpairs) << "%)" << std::endl;
    
    logfile.close();
    
    return 0;
}