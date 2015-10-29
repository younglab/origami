#include <api/BamReader.h>
#include <api/BamWriter.h>


using namespace BamTools;

int main(int argc,char **argv) {
    const char *leftfile = argv[1];
    const char *rightfile = argv[2];
    const char *outputfile = argv[3];
    BamReader left, right;
    
    left.Open(leftfile);
    right.Open(rightfile);
    
    const SamHeader header = left.GetHeader();
    const RefVector references = left.GetReferenceData();
    
    BamWriter writer;
    writer.Open(outputfile,header,references);
    
    BamAlignment l, r;
    while( left.GetNextAlignment(l) && right.GetNextAlignment(r)) {
        l.MateRefID = r.RefID;
        l.MatePosition = r.Position;
        l.SetIsPaired(true);
        l.SetIsMateMapped(r.IsMapped());
        l.SetIsMateReverseStrand(r.IsReverseStrand());
        l.SetIsFirstMate(true);
        writer.SaveAlignment(l);
        
        r.MateRefID = l.RefID;
        r.MatePosition = l.Position;
        r.SetIsPaired(true);
        r.SetIsMateMapped(l.IsMapped());
        r.SetIsMateReverseStrand(l.IsReverseStrand());
        r.SetIsSecondMate(true);
        writer.SaveAlignment(r);
    }
    
    writer.Close();
    left.Close();
    right.Close();
    
    return 0;
}