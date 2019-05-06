//
//  scrinvex.h
//  scrinvex
//
//  Created by Aaron Graubert on 2/26/19.
//  Copyright © 2019 Aaron Graubert. All rights reserved.
//

#ifndef scrinvex_h
#define scrinvex_h

#include <GTF.h>
#include <BamReader.h>

using namespace rnaseqc;

namespace scrinvex {
    typedef std::tuple<unsigned long, unsigned long, unsigned long> countTuple;
    
    class InvexCounter {
        // barcode -> counts
        std::unordered_map<std::string,  countTuple> counts;
        
    public:
        InvexCounter() : counts() {
            
        }
        
        countTuple& getCounts(const std::string&);
        std::set<std::string>& getBarcodes(std::set<std::string>&) const;
    };
    
    // gid -> invex counter
    typedef std::unordered_map<std::string, InvexCounter> geneCounters;
    
    // gene id -> (genic aligned length, exonic aligned length)
    typedef std::unordered_map<std::string, std::tuple<unsigned int, unsigned int> > alignmentLengthTracker;
    
    void countRead(geneCounters&, std::list<Feature>&, Alignment&, chrom, const std::unordered_set<std::string>&, InvexCounter* = nullptr);
    void dropFeatures(std::list<Feature>&, geneCounters&, std::ostream&);
//    void dropFeatures(std::list<Feature>&, geneCounters&, std::ostream&);
    void trimFeatures(Alignment&, std::list<Feature>&, geneCounters&, std::ostream&);
//    void trimFeatures(Alignment&, std::list<Feature>&, geneCounters&, std::ostream&);
    chrom getChrom(Alignment&, SeqLib::HeaderSequenceVector&);
    
    const std::size_t GENIC_ALIGNED_LENGTH = 0, EXONIC_ALIGNED_LENGTH = 1, INTRONS = 0, JUNCTIONS = 1, EXONS = 2;
    const std::string BARCODE_TAG = "CB", UMI_TAG = "UB", MISMATCH_TAG = "NM";
    
    extern unsigned int missingBC, missingUMI, skippedBC;
    extern std::unordered_map<std::string, unsigned long> intergenicCounts; //bc -> readcounts for intergenic reads
}

#endif /* scrinvex_h */
