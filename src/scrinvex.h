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

class InvexCounter {
    unsigned long introns, exons, junctions;

public:
    InvexCounter() : introns(0ul), exons(0ul), junctions(0ul) {

    }

    void countRead(std::list<Feature>&, Alignment&, chrom);
    friend std::ofstream& operator<<(std::ofstream&, const InvexCounter&);
};

// gene id -> (genic aligned length, exonic aligned length)
typedef std::unordered_map<std::string, std::tuple<unsigned int, unsigned int> > alignmentLengthTracker;
// barcode -> invex counter
typedef std::unordered_map<std::string, InvexCounter> geneCounters;

void dropFeatures(std::list<Feature>&);
chrom getChrom(Alignment&, SeqLib::HeaderSequenceVector&);
std::ofstream& operator<<(std::ofstream&, const InvexCounter&);
std::ofstream& operator<<(std::ofstream&, const geneCounters&);

const std::size_t GENIC_ALIGNED_LENGTH = 0, EXONIC_ALIGNED_LENGTH = 1;
const std::string BARCODE_TAG = "CB", UMI_TAG = "UB", MISMATCH_TAG = "NM";


#endif /* scrinvex_h */
