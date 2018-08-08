#include <algorithm>
#undef NDEBUG
#include <cassert>
#include <stdio.h>
#include <vector>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/kmer.h>

static const uint32_t filterMetadataSize = 256;
static const uint8_t INT_WIDTH = 0x40;

using namespace seqan;

int main()
{
    uint64_t noBins{3};
    uint64_t kmerSize{3};



    // ==========================================================================
    // Test functionality of the offset parameter
    // ==========================================================================
    std::cout<< "Offset functionality tests" <<std::endl; 
    typedef DirectAddressing       TSpec;
    typedef Uncompressed TFilter;
    typedef SimpleShape TShapeSpec;
    unsigned const offset = 3;
    unsigned const offset2 = 2;
    unsigned const overlapping = 1;

    std::vector<DnaString> dnas{DnaString("AAA"), DnaString("AAT"), DnaString("AAC"), DnaString("AAG"), DnaString("ATA"), DnaString("ATT"), DnaString("ATC"), DnaString("ATG"), DnaString("ACA"), DnaString("ACT"), DnaString("ACC"), DnaString("ACG"), DnaString("AGA"), DnaString("AGT"), DnaString("AGC"), DnaString("AGG"), DnaString("TAA"), DnaString("TAT"), DnaString("TAC"), DnaString("TAG"), DnaString("TTA"), DnaString("TTT"), DnaString("TTC"), DnaString("TTG"), DnaString("TCA"), DnaString("TCT"), DnaString("TCC"), DnaString("TCG"), DnaString("TGA"), DnaString("TGT"), DnaString("TGC"), DnaString("TGG"), DnaString("CAA"), DnaString("CAT"), DnaString("CAC"), DnaString("CAG"), DnaString("CTA"), DnaString("CTT"), DnaString("CTC"), DnaString("CTG"), DnaString("CCA"), DnaString("CCT"), DnaString("CCC"), DnaString("CCG"), DnaString("CGA"), DnaString("CGT"), DnaString("CGC"), DnaString("CGG"), DnaString("GAA"), DnaString("GAT"), DnaString("GAC"), DnaString("GAG"), DnaString("GTA"), DnaString("GTT"), DnaString("GTC"), DnaString("GTG"), DnaString("GCA"), DnaString("GCT"), DnaString("GCC"), DnaString("GCG"), DnaString("GGA"), DnaString("GGT"), DnaString("GGC"), DnaString("GGG")};
    KmerFilter<Dna, TSpec, TFilter, TShapeSpec, offset> nonoverlapping (noBins, kmerSize);
    KmerFilter<Dna, TSpec, TFilter, TShapeSpec, offset2> offset_2 (noBins, kmerSize);
    KmerFilter<Dna, TSpec, TFilter, TShapeSpec, overlapping > standard (noBins, kmerSize);
    // ==========================================================================
    // inserting all k-mers 4^3
    // ==========================================================================
    for (DnaString dna: dnas){
        insertKmer(nonoverlapping, dna, 1,1);
        insertKmer(standard, dna, 1,1);
        insertKmer(offset_2, dna, 1,1);
    }
    
    // ==========================================================================
    // omitting output for 4 different patterns and 3 different offsets.
    // ==========================================================================
    std::vector<DnaString> ps{DnaString("AAAGGGTTTCCC"),DnaString("AAAGGGTTTCC"),DnaString("AAAGGGTTTC"),DnaString("AAAGGGTTT")};
    for(DnaString p: ps){
        std::vector<uint32_t> counts(nonoverlapping.noOfBins, 0);
        std::vector<uint32_t> counts2(standard.noOfBins, 0);
        std::vector<uint32_t> counts3(offset_2.noOfBins, 0);       
        select(nonoverlapping,counts, p);
        select(standard,counts2, p);
        select(offset_2,counts3, p);
        std::cout<< "testing offset for pattern " <<p <<std::endl;  
        uint16_t x = (length(p) - kmerSize) % offset2;
        uint16_t exp= (x!=0)? (1 + (length(p) - kmerSize + offset2 - x) / offset2): ((length(p) - kmerSize + offset2 - x) / offset2);
        
        std::cout<< "k=3, o=1: actual: " <<counts2[1] << " expected: " << std::ceil((double)(length(p)-3+1)/overlapping) <<std::endl;   
        std::cout<< "k=3, o=2: actual: " <<counts3[1] << " expected: " << exp <<std::endl; 
        std::cout<< "k=3, o=3: actual: " <<counts[1] << " expected: " << std::ceil((double)length(p)/offset) <<std::endl;
        

    }
    return 0;
}
#define NDEBUG
