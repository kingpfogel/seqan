// ==========================================================================
//                                 build_filter.cpp
// ==========================================================================
// Copyright (c) 2017-2022, Temesgen H. Dadi, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Temesgen H. Dadi or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL TEMESGEN H. DADI OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#define BUILD_FILTER
// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------
#include <string>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <future>
#include <thread>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------
#include <seqan/index.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "../yara/store_seqs.h"
#include "../yara/misc_types.h"
#include "../yara/bits_matches.h"
#include "../yara/misc_options.h"
#include "../yara/misc_options_dis.h"
#include "../yara/index_fm.h"


using namespace seqan;

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString      kmersDir;
    CharString      bigIndexDir;
    CharString      indicesDir;

    uint32_t        numberOfBins;
    unsigned        threadsCount;

    uint64_t        contigsSize;
    uint64_t        contigsLength;
    uint64_t        contigsSum;

    Options() :
    numberOfBins(64),
    threadsCount(1)
    {
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "count_unique_kmers");
    setShortDescription(parser, "Count unique kmers");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE FILES DIR \\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE FILE DIR"));

    setHelpText(parser, 0, "A directory containing the kmers for each bin.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "INDICES DIR"));
    setHelpText(parser, 1, "A directory containing fm indices for each bin.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "BIG INDICES DIR"));
    setHelpText(parser, 2, "A directory containing the big indice.");

    addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins (indices)",
                                     ArgParseOption::INTEGER));

    setMinValue(parser, "number-of-bins", "1");
    setMaxValue(parser, "number-of-bins", "10000");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threadsCount);
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.kmersDir, parser, 0);
    getArgumentValue(options.indicesDir, parser, 1);
    getArgumentValue(options.bigIndexDir, parser, 2);

    if (isSet(parser, "number-of-bins")) getOptionValue(options.numberOfBins, parser, "number-of-bins");
    if (isSet(parser, "threads")) getOptionValue(options.threadsCount, parser, "threads");

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function get_unique_kmers()
// ----------------------------------------------------------------------------

template<typename TContigsSize, typename TContigsLen, typename TContigsSum>
void get_unique_kmers(Options & options)
{
    typedef YaraFMConfig<TContigsSize, TContigsLen, TContigsSum>    TIndexConfig;
    typedef FMIndex<void, TIndexConfig>                             TIndexSpec;
    typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;

    std::string comExt = commonExtension(options.kmersDir, options.numberOfBins);
    std::map<CharString, uint32_t> bin_map;
    std::map<uint32_t, uint32_t> big_bin_map;

    std::vector<uint64_t> uniq_counts(options.numberOfBins, 0);
    std::vector<uint64_t> kmer_counts(options.numberOfBins, 0);

    typedef SeqStore<void, YaraContigsConfig<> >                    TContigs;

    // Fill bin_map: bin_map[ref_name] = binNo
    for (uint32_t binNo = 0; binNo < options.numberOfBins; ++binNo)
    {
        CharString fm_index_file;
        appendFileName(fm_index_file, options.indicesDir, binNo);

        TContigs tmpContigs;
        if (!open(tmpContigs, toCString(fm_index_file), OPEN_RDONLY))
            throw RuntimeError("Error while opening tmp_contigs.");

        for (uint32_t i = 0; i < length(tmpContigs.names); ++i)
        {
            CharString s = (CharString)tmpContigs.names[i];
            bin_map[s] = binNo;
        }
        clear(tmpContigs);
    }

    TContigs allContigs;
    if (!open(allContigs, toCString(options.bigIndexDir), OPEN_RDONLY))
        throw RuntimeError("Error while opening allContigs.");

    // Fill big_bin_map: big_bin_map[FM_ID] = binNo
    for (uint32_t i = 0; i < length(allContigs.names); ++i)
    {
        CharString s = (CharString)allContigs.names[i];
        big_bin_map[i] = bin_map[s];
        // if (i == 0)
            // std::cout << "big_bin_map.i = " << i  << "\tBin = " << bin_map[s] << "\tRef Name = " << s << std::endl;
        // std::cout << i  << "\t" << bin_map[s] << '\n';
    }

    TIndex big_fm_index;
    if (!open(big_fm_index, toCString(options.bigIndexDir), OPEN_RDONLY))
        throw RuntimeError("Error while opening big_fm_index");


    Semaphore thread_limiter(options.threadsCount);
    std::vector<std::future<void>> tasks;

    for (uint32_t binNo = 0; binNo < options.numberOfBins; ++binNo)
    {
        tasks.emplace_back(std::async([=, &thread_limiter, &uniq_counts, &kmer_counts, &big_fm_index] {
        Critical_section _(thread_limiter);

        uint32_t batchSize = 10000;
        Finder<TIndex> finder(big_fm_index);

        CharString fastaFile;
        appendFileName(fastaFile, options.kmersDir, binNo);
        append(fastaFile, comExt);
        SeqFileIn seqFileIn;
        if (!open(seqFileIn, toCString(fastaFile)))
        {
            std::cerr <<"Unable to open contigs File: " << toCString(fastaFile) << std::endl;
            exit(1);
        }
        StringSet<CharString> ids;
        StringSet<Dna5String> seqs;

        while(!atEnd(seqFileIn))
        {
            readRecords(ids, seqs, seqFileIn, batchSize);
            // reverse(seqs);
            uint32_t len = length(seqs);
            kmer_counts[binNo] += len;
            for(uint32_t i = 0; i<len; ++i)
            {
                bool found = false;
                //std::cout << "There should be exactly 100 lines of this" << '\n';
                bool uniq = false;
                while (find(finder, seqs[i]))
                {
                    if (!found)
                    {
                        uniq = true;
                    }
                    found = true;
                    //std::cout << "There should be at least 100 lines of this" << '\n';
                    auto pos = position(finder);
                    uint32_t rID = getValueI1(pos);
                    auto bi = big_bin_map.find(rID);
                    // std::cout << binNo << '\t' << rID << '\t' << bi->second << '\n';
                    // std::cout << "The found rID is " << rID << " The proposed binNo is " << bi->second << '\n';
                    // std::cout << rID << "\t" << bi->second << '\n';

                    if(bi != big_bin_map.end() && binNo != bi->second)
                    {
                        uniq = false;
                        break;
                    }
                }
                if (uniq)
                {
                    ++uniq_counts[binNo];
                }
                if (not found)
                    std::cout << "kmer " << seqs[i] << " in bin " << binNo << " not found." << '\n';
                goBegin(finder);
                clear(finder);    // move Finder to the beginning of the text
            }
            clear(ids);
            clear(seqs);
        }
        close(seqFileIn);
        }));
    }
    for (auto &&task : tasks)
    {
        task.get();
    }

    std::cout << "binNo\tkmer count\tuniq count\n";
    for (uint32_t binNo = 0; binNo < options.numberOfBins; ++binNo)
    {
        std::cout << binNo  << "\t" <<  kmer_counts[binNo] << "\t\t" << uniq_counts[binNo] << std::endl ;
    }
}

template<typename TContigsSize, typename TContigsLen>
inline void get_unique_kmers(Options & options)
{
    if (options.contigsSum < MaxValue<uint32_t>::VALUE)
    {
        get_unique_kmers<TContigsSize, TContigsLen, uint32_t>(options);
    }
    else
    {
        get_unique_kmers<TContigsSize, TContigsLen, uint64_t>(options);
    }
}

template<typename TContigsSize>
inline void get_unique_kmers(Options & options)
{
    if (options.contigsLength < MaxValue<uint32_t>::VALUE)
    {
        get_unique_kmers<TContigsSize, uint32_t>(options);
    }
    else
    {
        get_unique_kmers<TContigsSize, uint64_t>(options);
    }
}

inline void get_unique_kmers(Options & options)
{
    String<uint64_t> limits;

    CharString contigsLimitFile(options.bigIndexDir);
    append(contigsLimitFile, ".txt.size");
    open(limits, toCString(contigsLimitFile), OPEN_RDONLY);

    // std::cout << limits[1] << ", "  << limits[0] << ", "  << limits[2] << "\n";

    options.contigsLength = limits[0];
    options.contigsSize = limits[1];
    options.contigsSum = limits[2];

    if (options.contigsSize < MaxValue<uint8_t>::VALUE)
    {
        get_unique_kmers<uint8_t>(options);
    }
    else
    {
        get_unique_kmers<uint16_t>(options);
    }
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    try
    {
        get_unique_kmers(options);
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
