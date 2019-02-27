//
//  scrinvex.cpp
//  scrinvex
//
//  Created by Aaron Graubert on 2/26/19.
//  Copyright © 2019 Aaron Graubert. All rights reserved.
//

#include "scrinvex.h"
#include <stdio.h>
#include <memory>
#include <args.hxx>
#include <Metrics.h>
#include <Expression.h>

using namespace args;
using namespace std;

int main(int argc, char* argv[])
{
    ArgumentParser parser("SCRINVEX");
    HelpFlag help(parser, "help", "Display this message and quit", {'h', "help"});
    Positional<string> gtfFile(parser, "gtf", "The input GTF file containing features to check the bam against");
    Positional<string> bamFile(parser, "bam", "The input SAM/BAM file containing reads to process");
    Positional<string> outputDir(parser, "output", "Output directory");
    
    parser.ParseCLI(argc, argv);
    
    Feature line; //current feature being read from the gtf
    ifstream reader(gtfFile.Get());
    if (!reader.is_open())
    {
        cerr << "Unable to open GTF file: " << gtfFile.Get() << endl;
        return 10;
    }
    
    unsigned long featcnt = 0, alignmentCount = 0;
    map<chrom, list<Feature>> features;
    while (reader >> line)
    {
        if (line.type == "gene" || line.type == "exon")
        {
            features[line.chromosome].push_back(line);
            ++featcnt;
        }
    }
    for (auto entry : features)
        entry.second.sort(compIntervalStart);
    cout << featcnt << " features loaded" << endl;
    
    geneCounters counts; //barcode -> counts
    
    {
        SeqlibReader bam;
        if (!bam.open(bamFile.Get()))
        {
            cerr << "Unable to open BAM file: " << bamFile.Get() << endl;
            return 10;
        }
        SeqLib::BamHeader header = bam.getHeader();
        SeqLib::HeaderSequenceVector sequences = header.GetHeaderSequenceVector();
        
        bool hasOverlap = false;
        for(auto sequence = sequences.begin(); sequence != sequences.end(); ++sequence)
        {
            chrom chrom = chromosomeMap(sequence->Name);
            if (features.find(chrom) != features.end())
            {
                hasOverlap = true;
                break;
            }
        }
        if (!hasOverlap)
        {
            cerr << "BAM file shares no contigs with GTF" << endl;
            return 11;
        }
        
        Alignment alignment;
        
        int32_t last_position = 0; // For some reason, htslib has decided that this will be the datatype used for positions
        chrom current_chrom = 0;
        
        while (bam.next(alignment))
        {
            ++alignmentCount;
            if (!(alignment.SecondaryFlag() || alignment.QCFailFlag()) && alignment.MappedFlag())
            {
                chrom chr = getChrom(alignment, sequences); //parse out a chromosome shorthand
                if (chr != current_chrom)
                {
                    dropFeatures(features[current_chrom]);
                    current_chrom = chr;
                }
                else if (last_position > alignment.Position())
                    cerr << "Warning: The input bam does not appear to be sorted. An unsorted bam will yield incorrect results" << endl;
                last_position = alignment.Position();
                trimFeatures(alignment, features[chr]); //drop features that appear before this read
                string barcode;
                alignment.GetZTag(BARCODE_TAG, barcode);
                counts[barcode].countRead(features[chr], alignment, chr);
            }
        }
    }
    
    cout << "Generating Report" << endl;
    
    ofstream report(outputDir.Get() + "/" + "report.tsv");
//    report << "Barcode\tIntrons\tJunctions\tExons" << endl;
//    for (auto entry : counts)
//    {
//        report << entry.first << "\t";
//        report << entry.second << endl;
//    }
    report << counts;
    report.close();
    
    
}

void dropFeatures(std::list<Feature> &features)
{
    for (auto feat = features.begin(); feat != features.end(); ++feat) if (feat->type == "gene") fragmentTracker.erase(feat->feature_id);
    features.clear();
}

void InvexCounter::countRead(std::list<Feature> &features, Alignment &alignment, chrom chromosome)
{
    vector<Feature> alignedSegments;
    extractBlocks(alignment, alignedSegments, chromosome, false);
    alignmentLengthTracker lengths;
    string umi;
    alignment.GetZTag(UMI_TAG, umi);
    for (Feature &segment : alignedSegments)
    {
        //Check memory usage in profiler. This isn't the best syntax, but hopefully shared_ptr understands
        shared_ptr<list<Feature> > intersections = shared_ptr<list<Feature> >(intersectBlock(segment, features));
        for (Feature &genomeFeature : *intersections)
        {
            // strandedness?
            if (genomeFeature.type == "exon" && fragmentTracker[genomeFeature.gene_id].count(umi) == 0)
                get<exonicAlignedLength>(lengths[genomeFeature.gene_id]) += partialIntersect(genomeFeature, segment);
            else if (genomeFeature.type == "gene" && fragmentTracker[genomeFeature.feature_id].count(umi) == 0)
                get<genicAlignedLength>(lengths[genomeFeature.gene_id]) += partialIntersect(genomeFeature, segment);
        }
    }
    for (auto entry : lengths)
    {
        unsigned int genicLength = get<genicAlignedLength>(entry.second), exonicLength = get<exonicAlignedLength>(entry.second);
        if (genicLength > 0)
        {
            if (genicLength > exonicLength)
            {
                if (exonicLength) ++(this->junctions);
                else ++(this->introns);
            }
            else ++(this->exons);
            fragmentTracker[entry.first].insert(umi);
        }
    }
}

chrom getChrom(Alignment &alignment, SeqLib::HeaderSequenceVector &sequences)
{
    return chromosomeMap(sequences[alignment.ChrID()].Name);
}

ofstream& operator<<(ofstream &stream, const InvexCounter &counter)
{
    stream << counter.introns << "\t" << counter.junctions << "\t" << counter.exons;
    return stream;
}

std::ofstream& operator<<(std::ofstream &stream, const geneCounters &counters)
{
    stream << "Barcode\tIntrons\tJunctions\tExons" << endl;
    for (auto entry : counters)
    {
        stream << entry.first << "\t";
        stream << entry.second << endl;
    }
    return stream;
}

