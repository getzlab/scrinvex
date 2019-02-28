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
    class InvexCounter;
    
    // gid -> invex counter
    typedef std::unordered_map<std::string, InvexCounter> geneCounters;
    
    class InvexCounter {
        // barcode -> counts
        std::unordered_map<std::string, std::tuple<unsigned long, unsigned long, unsigned long> > counts;
        
    public:
        InvexCounter() : counts() {
            
        }
        
        std::tuple<unsigned long, unsigned long, unsigned long>& getCounts(const std::string&);
        std::set<std::string>& getBarcodes(std::set<std::string>&) const;
        
//        friend void ::countRead(geneCounters&, std::list<Feature>&, Alignment&, chrom);
//        friend void dropFeatures(std::list<Feature>&, geneCounters&, std::ostream&, std::ostream&, std::ostream&);
//        friend void trimFeatures(Alignment&, std::list<Feature>&, geneCounters&, std::ostream&, std::ostream&, std::ostream&);
    };
    
    // gene id -> (genic aligned length, exonic aligned length)
    typedef std::unordered_map<std::string, std::tuple<unsigned int, unsigned int> > alignmentLengthTracker;
    
    void countRead(geneCounters&, std::list<Feature>&, Alignment&, chrom);
    void dropFeatures(std::list<Feature>&, geneCounters&, std::ostream&, std::ostream&, std::ostream&);
    void trimFeatures(Alignment&, std::list<Feature>&, geneCounters&, std::ostream&, std::ostream&, std::ostream&);
    chrom getChrom(Alignment&, SeqLib::HeaderSequenceVector&);
    
    const std::size_t GENIC_ALIGNED_LENGTH = 0, EXONIC_ALIGNED_LENGTH = 1, INTRONS = 0, JUNCTIONS = 1, EXONS = 2;
    const std::string BARCODE_TAG = "CB", UMI_TAG = "UB", MISMATCH_TAG = "NM";
    
    extern unsigned int missingBC, missingUMI;
}

#endif /* scrinvex_h */