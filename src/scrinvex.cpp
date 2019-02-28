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
    try
    {
        parser.ParseCLI(argc, argv);

        if (!gtfFile) throw ValidationError("No GTF file provided");
        if (!bamFile) throw ValidationError("No BAM file provided");
        if (!outputDir) throw ValidationError("No output directory provided");

        const string SAMPLENAME = sampleName ? sampleName.Get() : boost::filesystem::path(bamFile.Get()).filename().string();

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
            if (line.type == FeatureType::Gene || line.type == FeatureType::Exon)
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

            ofstream introns(outputDir.Get() + "/" + SAMPLENAME + ".introns.tsv");
            ofstream junctions(outputDir.Get() + "/" + SAMPLENAME + ".junctions.tsv");
            ofstream exons(outputDir.Get() + "/" + SAMPLENAME + ".exons.tsv");
            introns << "# sparse" << endl << "gene_id\tbarcode\tintrons" << endl;
            junctions << "# sparse" << endl << "gene_id\tbarcode\tjunctions" << endl;
            exons << "# sparse" << endl << "gene_id\tbarcode\texons" << endl;

            cout << "Parsing BAM" << endl;

            while (bam.next(alignment))
            {
                ++alignmentCount;
                if (!(alignment.SecondaryFlag() || alignment.QCFailFlag()) && alignment.MappedFlag())
                {
                    chrom chr = getChrom(alignment, sequences); //parse out a chromosome shorthand
                    if (chr != current_chrom)
                    {
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
            for (auto contig : features) if (contig.second.size()) dropFeatures(contig.second, counts, introns, junctions, exons);
            introns.close();
            junctions.close();
            exons.close();

            if (missingUMI + missingBC)
                cerr << "There were " << missingBC << " reads without a barcode (CB) and " << missingUMI << " reads without a UMI (UB)" << endl;
            
            if (skippedBC)
                cerr << "Skipped " << skippedBC << " reads with barcodes not listed in " << barcodeFile.Get() << endl;
        }

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
        for (auto entry : this->counts) destination.insert(entry.first);
        return destination;
    }

    void countRead(geneCounters &counts, std::list<Feature> &features, Alignment &alignment, chrom chromosome, const std::unordered_set<std::string> &goodBarcodes)
    {
        vector<Feature> alignedSegments;
        extractBlocks(alignment, alignedSegments, chromosome, false);
        alignmentLengthTracker lengths;
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
        if (goodBarcodes.count(barcode) == 0)
        {
            ++skippedBC;
            return;
        }
        for (Feature &segment : alignedSegments)
        {
            shared_ptr<list<Feature> > intersections = shared_ptr<list<Feature> >(intersectBlock(segment, features));
            for (Feature &genomeFeature : *intersections)
            {
                if (genomeFeature.type == FeatureType::Exon && fragmentTracker[genomeFeature.gene_id].count(umi) == 0)
                    get<EXONIC_ALIGNED_LENGTH>(lengths[genomeFeature.gene_id]) += partialIntersect(genomeFeature, segment);
                else if (genomeFeature.type == FeatureType::Gene && fragmentTracker[genomeFeature.feature_id].count(umi) == 0)
                    get<GENIC_ALIGNED_LENGTH>(lengths[genomeFeature.gene_id]) += partialIntersect(genomeFeature, segment);
            }
        }
        for (auto entry : lengths)
        {
            unsigned int genicLength = get<GENIC_ALIGNED_LENGTH>(entry.second), exonicLength = get<EXONIC_ALIGNED_LENGTH>(entry.second);
            if (genicLength > 0)
            {
                if (genicLength > exonicLength)
                {
                    if (exonicLength) get<JUNCTIONS>(counts[entry.first].getCounts(barcode)) += 1;
                    else get<INTRONS>(counts[entry.first].getCounts(barcode)) += 1;
                }
                else get<EXONS>(counts[entry.first].getCounts(barcode)) += 1;
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
        features.erase(features.begin(), cursor);
    }
}
