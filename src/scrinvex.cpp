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
#include <boost/filesystem.hpp>
#include <Metrics.h>
#include <Expression.h>

using namespace args;
using namespace std;
using namespace scrinvex;

namespace scrinvex {
    unsigned int missingBC, missingUMI, skippedBC;
}

int main(int argc, char* argv[])
{
    ArgumentParser parser("SCRINVEX");
    HelpFlag help(parser, "help", "Display this message and quit", {'h', "help"});
    Positional<string> gtfFile(parser, "gtf", "The input GTF file containing features to check the bam against");
    Positional<string> bamFile(parser, "bam", "The input SAM/BAM file containing reads to process");
    Positional<string> outputDir(parser, "output", "Output Directory");
    ValueFlag<string> sampleName(parser, "sample", "The name of the current sample.  Default: The bam's filename", {'s', "sample"});
    ValueFlag<string> barcodeFile(parser, "barcodes", "Path to filtered barcodes.tsv file from cellranger. Only barcodes listed in the file will be used. Default: All barcodes present in bam", {'b', "barcodes"});
    ValueFlag<unsigned int> mappingQualityThreshold(parser,"quality", "Set the lower bound on read quality for coverage counting. Reads below this quality are skipped. Default: 255", {'q', "quality"});
    try
    {
        parser.ParseCLI(argc, argv);

        if (!gtfFile) throw ValidationError("No GTF file provided");
        if (!bamFile) throw ValidationError("No BAM file provided");
        if (!outputDir) throw ValidationError("No output directory provided");

        const string SAMPLENAME = sampleName ? sampleName.Get() : boost::filesystem::path(bamFile.Get()).filename().string();
        const unsigned int MAPQ = mappingQualityThreshold ? mappingQualityThreshold.Get() : 255u;

        Feature line; //current feature being read from the gtf
        ifstream reader(gtfFile.Get());
        if (!reader.is_open())
        {
            cerr << "Unable to open GTF file: " << gtfFile.Get() << endl;
            return 10;
        }

        cout << "                                                                           ***                         " << endl;
        cout << "                       ***                                              *****                         " << endl;
        cout << "                     *****                   ***         *            *****                           " << endl;
        cout << "              **   *****  ***********     * ******       **          ****                             " << endl;
        cout << "            ***** ****    ***********    ** ********     ***         ***   ****        ***        " << endl;
        cout << "         ******   ***     ***  *****    *** ***  *****   ***     **  ***   ******    *****        " << endl;
        cout << "       ******     ***     ********      *** ***    ***** ***     *** ******* ***********          " << endl;
        cout << "     ******       ***     *******       *** ***      *** ***     *** *******   *******            " << endl;
        cout << "   ************** ***     *********     *** ***      *** ***     *** ***      *********           " << endl;
        cout << " **************** ***     ***   *****   *** ***      *** ***     *** ***     *****  *****         " << endl;
        cout << "           *****  *****   ***     ***** *** ***      *** ***     *** ***    ***       *****       " << endl;
        cout << "         *****      ***** ***       *** *** ***      *** ***    **** ***** **           *****     " << endl;
        cout << "       *****          *****             *** ***      *** ***  *****    *****              *****   " << endl;
        cout << "     *****              ****                ***          ********        *****              ****  " << endl;
        cout << "    ****                  ***               **           ******            *****              *** " << endl;
        cout << "                            **                           ****                ***                **" << endl;
        cout << "                                                         **                                       " << endl;
        
        unordered_set<string> goodBarcodes;
        if (barcodeFile)
        {
            cout << "Reading barcodes" << endl;
            ifstream barcodeReader(barcodeFile.Get());
            string barcode;
            while (barcodeReader >> barcode) goodBarcodes.insert(barcode);
            cout << "Filtering input using " << goodBarcodes.size() << " barcodes" << endl;
        }

        cout << "Parsing GTF" << endl;

        unsigned long featcnt = 0, alignmentCount = 0;
        map<chrom, list<Feature>> features;
        while (reader >> line)
        {
            // Only record Genes and Exons. Transcripts not important for scrinvex
            if (line.type == FeatureType::Gene || line.type == FeatureType::Exon)
            {
                features[line.chromosome].push_back(line);
                ++featcnt;
            }
        }
        // Sort all features by position
        for (auto entry : features)
            entry.second.sort(compIntervalStart);
        cout << featcnt << " features loaded" << endl;

        geneCounters counts; //barcode -> counts

        SeqlibReader bam;
        if (!bam.open(bamFile.Get()))
        {
            cerr << "Unable to open BAM file: " << bamFile.Get() << endl;
            return 10;
        }
        
        // Get the list of contigs present in bam header
        SeqLib::BamHeader header = bam.getHeader();
        SeqLib::HeaderSequenceVector sequences = header.GetHeaderSequenceVector();
        
        // Intersect bam header with gtf contigs to make sure they share the same naming scheme
        bool hasOverlap = false;
        for(auto sequence : sequences)
        {
            chrom chrom = chromosomeMap(sequence.Name);
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
        
        //use boost to ensure that the output directory exists before the metrics are dumped to it
        if (!boost::filesystem::exists(outputDir.Get()))
        {
            boost::filesystem::create_directories(outputDir.Get());
        }
        
        // Open all output files
        ofstream introns(outputDir.Get() + "/" + SAMPLENAME + ".introns.tsv");
        ofstream junctions(outputDir.Get() + "/" + SAMPLENAME + ".junctions.tsv");
        ofstream exons(outputDir.Get() + "/" + SAMPLENAME + ".exons.tsv");
        introns << "gene_id\tbarcode\tintrons" << endl;
        junctions << "gene_id\tbarcode\tjunctions" << endl;
        exons << "gene_id\tbarcode\texons" << endl;
        
        cout << "Parsing BAM" << endl;
        
        while (bam.next(alignment))
        {
            ++alignmentCount;
            // Only consider uniquely mapped reads
            if (!(alignment.SecondaryFlag() || alignment.QCFailFlag()) && alignment.MappedFlag() && alignment.MapQuality() >= MAPQ)
            {
                chrom chr = getChrom(alignment, sequences); //parse out a chromosome shorthand
                if (chr != current_chrom)
                {
                    // If we've switched chromosomes, drop all features from that chromosome
                    // Saves memory and also writes out the coverage data
                    dropFeatures(features[current_chrom], counts, introns, junctions, exons);
                    current_chrom = chr;
                }
                else if (last_position > alignment.Position())
                    cerr << "Warning: The input bam does not appear to be sorted. An unsorted bam will yield incorrect results" << endl;
                last_position = alignment.Position();
                trimFeatures(alignment, features[chr], counts, introns, junctions, exons); //drop features that appear before this read
                countRead(counts, features[chr], alignment, chr, goodBarcodes);
            }
        }
        
        cout << "Finalizing data" << endl;
        // Drop all remaining genes to ensure their coverage data has been written
        for (auto contig : features) if (contig.second.size()) dropFeatures(contig.second, counts, introns, junctions, exons);
        introns.close();
        junctions.close();
        exons.close();
        
        if (missingUMI + missingBC)
            cerr << "There were " << missingBC << " reads without a barcode (CB) and " << missingUMI << " reads without a UMI (UB)" << endl;
        
        if (skippedBC)
            cerr << "Skipped " << skippedBC << " reads with barcodes not listed in " << barcodeFile.Get() << endl;

        return 0;
    }
    catch (args::Help)
    {
        cout << parser;
        return 4;
    }
    catch (args::ParseError &e)
    {
        cerr << parser << endl;
        cerr << "Argument parsing error: " << e.what() << endl;
        return 5;
    }
    catch (args::ValidationError &e)
    {
        cerr << parser << endl;
        cerr << "Argument validation error: " << e.what() << endl;
        return 6;
    }
    catch (fileException &e)
    {
        cerr << e.error << endl;
        return 10;
    }
    catch (invalidContigException &e)
    {
        cerr << "GTF referenced a contig not present in the FASTA: " << e.error << endl;
        return 11;
    }
    catch (gtfException &e)
    {
        cerr << "Failed to parse the GTF: " << e.error << endl;
        return 11;
    }
    catch (boost::filesystem::filesystem_error &e)
    {
        cerr << "Filesystem error:  " << e.what() << endl;
        return 8;
    }
    catch(ios_base::failure &e)
    {
        cerr << "Encountered an IO failure" << endl;
        cerr << e.what() << endl;
        return 10;
    }
    catch(std::bad_alloc &e)
    {
        cerr << "Memory allocation failure. Out of memory" << endl;
        cerr << e.what() << endl;
        return 10;
    }
    catch (...)
    {
        cerr << parser << endl;
        cerr << "Unknown error" << endl;
        return -1;
    }
}

namespace scrinvex {

    std::tuple<unsigned long, unsigned long, unsigned long>& InvexCounter::getCounts(const std::string &barcode)
    {
        return this->counts[barcode];
    }

    set<string>& InvexCounter::getBarcodes(set<string> &destination) const
    {
        // Insert all barcodes from this InvexCounter into a given set
        for (auto entry : this->counts) destination.insert(entry.first);
        return destination;
    }

    void countRead(geneCounters &counts, std::list<Feature> &features, Alignment &alignment, chrom chromosome, const std::unordered_set<std::string> &goodBarcodes)
    {
        // We don't expect goodBarcodes to change, so just grab its size once
        static const size_t n_barcodes = goodBarcodes.size();
        
        // Now parse the Cigar string to get the set of intervals that this read aligns over
        vector<Feature> alignedSegments;
        extractBlocks(alignment, alignedSegments, chromosome, false);
        
        alignmentLengthTracker lengths;
        
        // Extract the barcode and umi. Check that they're present and barcode is in the set of good barcodes
        string barcode, umi;
        if (!alignment.GetZTag(BARCODE_TAG, barcode))
        {
            ++missingBC;
            return;
        }
        if (!alignment.GetZTag(UMI_TAG, umi))
        {
            ++missingUMI;
            return;
        }
        if (n_barcodes > 0 && goodBarcodes.count(barcode) == 0)
        {
            ++skippedBC;
            return;
        }
        
        // Intersect all aligned segments with the list of features.
        for (Feature &segment : alignedSegments)
        {
            shared_ptr<list<Feature> > intersections = shared_ptr<list<Feature> >(intersectBlock(segment, features));
            for (Feature &genomeFeature : *intersections)
            {
                // Count the total number of read bases which align to genes and exons
                if (genomeFeature.type == FeatureType::Exon && fragmentTracker[genomeFeature.gene_id].count(umi) == 0)
                    get<EXONIC_ALIGNED_LENGTH>(lengths[genomeFeature.gene_id]) += partialIntersect(genomeFeature, segment);
                else if (genomeFeature.type == FeatureType::Gene && fragmentTracker[genomeFeature.feature_id].count(umi) == 0)
                    get<GENIC_ALIGNED_LENGTH>(lengths[genomeFeature.gene_id]) += partialIntersect(genomeFeature, segment);
            }
        }
        // For every gene that this read aligned to
        for (auto entry : lengths)
        {
            unsigned int genicLength = get<GENIC_ALIGNED_LENGTH>(entry.second), exonicLength = get<EXONIC_ALIGNED_LENGTH>(entry.second);
            if (genicLength > 0)
            {
                if (genicLength > exonicLength) // Read did not align entirely to exons
                {
                    if (exonicLength) get<JUNCTIONS>(counts[entry.first].getCounts(barcode)) += 1; // Read aligned a little to exons
                    else get<INTRONS>(counts[entry.first].getCounts(barcode)) += 1; // Read aligned entirely to introns
                }
                else get<EXONS>(counts[entry.first].getCounts(barcode)) += 1; // Read aligned entirely to exons
                
                // Now add the UMI to the tracker so we skip UMI duplicates
                fragmentTracker[entry.first].insert(umi);
            }
        }
    }

    chrom getChrom(Alignment &alignment, SeqLib::HeaderSequenceVector &sequences)
    {
        return chromosomeMap(sequences[alignment.ChrID()].Name);
    }

    void dropFeatures(std::list<Feature> &features, geneCounters &counts, std::ostream &introns, std::ostream &junctions, std::ostream &exons)
    {
        for (Feature &feat : features) if (feat.type == FeatureType::Gene) {
            // For all genes, dump their coverage data
            fragmentTracker.erase(feat.feature_id);
            InvexCounter &invex = counts[feat.feature_id];
            set<string> barcodes;
            for (const string &barcode : invex.getBarcodes(barcodes))
            {
                auto data = invex.getCounts(barcode);
                if (get<INTRONS>(data)) introns << feat.feature_id << "\t" << barcode << "\t" << get<INTRONS>(data) << endl;
                if (get<JUNCTIONS>(data)) junctions << feat.feature_id << "\t" << barcode << "\t" << get<JUNCTIONS>(data) << endl;
                if (get<EXONS>(data)) exons << feat.feature_id << "\t" << barcode << "\t" << get<EXONS>(data) << endl;
            }
        }
        features.clear();
    }
    void trimFeatures(Alignment &alignment, std::list<Feature> &features, geneCounters &counts, std::ostream &introns, std::ostream &junctions, std::ostream &exons)
    {
        auto cursor = features.begin();
        while (cursor != features.end() && cursor->end < alignment.Position())
        {
            if (cursor->type == FeatureType::Gene) {
                // For all genes, dump their coverage data
                fragmentTracker.erase(cursor->feature_id);
                InvexCounter &invex = counts[cursor->feature_id];
                set<string> barcodes;
                for (const string &barcode : invex.getBarcodes(barcodes))
                {
                    auto data = invex.getCounts(barcode);
                    if (get<INTRONS>(data)) introns << cursor->feature_id << "\t" << barcode << "\t" << get<INTRONS>(data) << endl;
                    if (get<JUNCTIONS>(data)) junctions << cursor->feature_id << "\t" << barcode << "\t" << get<JUNCTIONS>(data) << endl;
                    if (get<EXONS>(data)) exons << cursor->feature_id << "\t" << barcode << "\t" << get<EXONS>(data) << endl;
                }
            }
            ++cursor;
        }
        // Also erase all features up to (but not including) the cursor. They are now outside the search window
        // This saves memory and greatly improves runtime
        features.erase(features.begin(), cursor);
    }
}
