#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <iostream>
#include <cstdlib>

using namespace BamTools;


int main(int argc, const char **argv ) {
  if(argc < 3) {
    std::cerr << "bam-read-extension: <BAM file> <extension> <out BAM file>" << std::endl;
    return EXIT_FAILURE;
  }
  
  const char *bamfile = argv[1];
  long extension = std::atoi(argv[2]);
  const char *outfile = argv[3];
  
  BamReader bam;
  bam.Open(bamfile);
  
  const SamHeader header = bam.GetHeader();
  const RefVector references = bam.GetReferenceData();
  
  BamWriter out;
  out.Open(outfile,header,references);
  
  BamAlignment a;
  
  while(bam.GetNextAlignment(a)) {
    if(a.IsMapped()) {
      a.Length = extension;
    }
    
    out.SaveAlignment(a);
  }
  
  out.Close();
  bam.Close();
  
  return 0;
}