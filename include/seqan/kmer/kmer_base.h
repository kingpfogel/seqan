// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author:  Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
//          Enrico Seiler <enrico.seiler@fu-berlin.de>
// Adapted for a bachelor thesis by Niklas Neesemann
// ==========================================================================

#ifndef INCLUDE_SEQAN_KMER_KMER_BASE_H_
#define INCLUDE_SEQAN_KMER_KMER_BASE_H_

#include <sdsl/bit_vectors.hpp>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <valarray>
#include <algorithm>
#include <future>
#include <mutex>
// TODO change API
// TODO change types

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================
class Semaphore
{
    std::mutex m;
    std::condition_variable cv;
    int count;

public:
    Semaphore(int n) : count{n} {}
    void notify()
    {
        std::unique_lock<std::mutex> l(m);
        ++count;
        cv.notify_one();
    }
    void wait()
    {
        std::unique_lock<std::mutex> l(m);
        cv.wait(l, [this]{ return count!=0; });
        --count;
    }
};

class Critical_section
{
    Semaphore &s;
public:
    Critical_section(Semaphore &ss) : s{ss} { s.wait(); }
    ~Critical_section() { s.notify(); }
};

template <typename T>
int numDigits(T number)
{
    int digits = 0;
    if (number <= 0) digits = 1; // remove this line if '-' counts as a digit
    while (number) {
        number /= 10;
        digits++;
    }
    return digits;
}

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag k-mer Filter Tags
// --------------------------------------------------------------------------

//!\brief A tag for the IBF.
struct InterleavedBloomFilter_;
typedef Tag<InterleavedBloomFilter_> InterleavedBloomFilter;


//!\brief A tag for direct addressing.
struct DirectAddressing_;
typedef Tag<DirectAddressing_> DirectAddressing;

//!\brief A tag for the uncompressed FilterVector.
struct Uncompressed_;
typedef Tag<Uncompressed_> Uncompressed;

//!\brief A tag for the compressed FilterVector.
struct CompressedSimple_;
typedef Tag<CompressedSimple_> CompressedSimple;

//!\brief A tag for the compressed array FilterVector.
struct CompressedArray_;
typedef Tag<CompressedArray_> CompressedArray;

template<typename TSpec>
class FilterVector;

//template <bool>
//struct KmerOffset {};

//!\brief typedef used for selectHelper for overlapping k-mers
//typedef KmerOffset<true> CompleteCoverage;
//!\brief typedef used for selectHelper for non-overlapping k-mers
//typedef KmerOffset<false> IncompleteCoverage;
// --------------------------------------------------------------------------
// Class KmerFilter
// --------------------------------------------------------------------------

//!\brief The KmerFilter class.
template<typename TValue = Dna, typename TSpec = DirectAddressing, typename TFilterVector = Uncompressed, typename TShapeSpec = SimpleShape, unsigned offset = 1>
class KmerFilter;

// ==========================================================================
// Metafunctions
// ==========================================================================

//!\brief Type definition for variables.
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
struct Value<KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> >
{
    typedef uint32_t noOfBins;
    typedef uint16_t kmerSize;
    typedef uint64_t noOfBits;
    typedef uint64_t noOfBlocks;
    typedef uint16_t binWidth;
    typedef uint32_t blockBitSize;
    typedef uint8_t  intSize;
    typedef uint16_t filterMetadataSize;
    typedef uint8_t  noOfHashFunc;
    typedef uint8_t  shiftValue;
    typedef uint64_t preCalcValues;
    typedef uint64_t seedValue;
};

// --------------------------------------------------------------------------
// Metafunction MetafunctionName
// --------------------------------------------------------------------------

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function functionName()
// --------------------------------------------------------------------------

/*!
 * \brief Adds k-mers from a text to a bin in a given filter.
 * \param me The KmerFilter instance.
 * \param text The text from which the k-mers are to be added.
 * \param binNo The bin to add the k-mers to.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset, typename TString, typename TBin, typename TChunk>
inline void insertKmer(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> & me, TString const & text, TBin && binNo, TChunk && chunkNo)
{
    me.insertKmer(text, binNo, chunkNo);
}

/*!
 * \brief Sets the vectors for given bins to 0.
 * \param me The KmerFilter instance.
 * \param bins A vector containing the bin numbers.
 * \param threads The number of threads to use.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset, typename TInt1, typename TInt2>
inline void clear(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, std::vector<TInt1> & bins, TInt2&& threads)
{
    // me.clear(bins, static_cast<uint64_t>(threads));
    me.clear(bins, threads);
}

/*!
 * \brief Adds all k-mers from a fasta file to a bin of a given KmerFilter.
 * \param me The KmerFilter instance.
 * \param kmer The fasta file to process.
 * \param binNo The bin to add the k-mers to.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset, typename TInt>
inline void insertKmer(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, const char * fastaFile, TInt && binNo, [[maybe_unused]] bool batch=false, [[maybe_unused]] uint8_t batchChunkNo=0)
{
    CharString id;
    String<TValue> seq;
    SeqFileIn seqFileIn;
    for (uint8_t i = 0; i < me.filterVector.noOfChunks; ++i)
    {
        if constexpr (std::is_same_v<TFilterVector, CompressedArray>)
        {
            if (batch)
                i = batchChunkNo;
        }
        if (!open(seqFileIn, fastaFile))
        {
            CharString msg = "Unable to open contigs file: ";
            append(msg, CharString(fastaFile));
            std::cerr << msg << std::endl;
            throw toCString(msg);
        }
        if (!batch)
            me.filterVector.decompress(i);
        while(!atEnd(seqFileIn))
        {
            readRecord(id, seq, seqFileIn);
            if(length(seq) < me.kmerSize)
                continue;
            insertKmer(me, seq, binNo, i);
        }
        if (!batch)
            me.filterVector.compress(i);
        close(seqFileIn); // No rewind() for FormattedFile ?

        if constexpr (std::is_same_v<TFilterVector, CompressedArray>)
        {
            if (batch)
                break;
        }
    }
}

/*!
 * \brief Adds all fasta files from a directory to the respective bins.
 * \param me The KmerFilter instance.
 * \param baseDir The directory containing the fasta files in a "bins" subdirectory.
 * \param threads Number of threads to use.
 *
 * The fasta files are expected to follow the pattern <baseDir>/bins/bin_xxxx.fasta, where xxxx stands for the bin
 * number. All bin numbers must have the same number of digits as the total number of bins. E.g. for 8192 bins, the
 * bins are expected to be named bin_0000.fasta, bin_0001.fasta, ..., bin_8191.fasta; or for 64 bins: bin_00.fasta,
 * bin_01.fasta, ..., bin_63.fasta.
 * Up to <threads> fasta files are added to the filterVector at the same time.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline void insertKmerDir(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, const char * baseDir, uint8_t threads)
{
    Semaphore thread_limiter(threads);
    // std::mutex mtx;
    std::vector<std::future<void>> tasks;

    uint32_t bins = me.noOfBins;
    for (uint8_t c = 0; c < me.filterVector.noOfChunks; ++c)
    {
        me.filterVector.decompress(c);
        for(uint32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/bins/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fasta"));
            tasks.emplace_back(
                std::async(std::launch::async, [=, &thread_limiter, &me] { // &mtx
                    Critical_section _(thread_limiter);
                    insertKmer(me, toCString(file), i, true, c);
                    // mtx.lock();
                    // std::cerr << "IBF Bin " << i << " done." << '\n';
                    // mtx.unlock();
                })
            );
        }

        for (auto &&task : tasks){
            task.get();
        }
        me.filterVector.compress(c);
    }
}

/*!
 * \brief Calculates the k-mer counts of a given text.
 * \param me The KmerFilter instance.
 * \param counts Vector of length binNo to save counts to.
 * \param text A single text to count all contained k-mers for.
 */
template<typename TValue, typename TSpec,  typename TFilterVector, typename TShapeSpec, unsigned offset>
inline void select(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, std::vector<uint32_t> & counts, String<TValue> const & text)
{
    me.select(counts, text);
}

/*!
 * \brief Returns the k-mer counts of a given text.
 * \param me The KmerFilter instance.
 * \param text A single text to count all contained k-mers for.
 * \returns std::vector<uint64_t> of size binNo containing counts.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline std::vector<uint32_t> select(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, String<TValue> const & text)
{
    std::vector<uint32_t> counts(me.noOfBins, 0);
    select(me, counts, text);
    return counts;
}

/*!
 * \brief Checks for which bins the counts of all k-mers in a text exceed a threshold.
 * \param me The KmerFilter instance.
 * \param selected Vector of length binNo to save true/false to.
 * \param text A single text to count all contained k-mers for.
 * \param threshold The minimal number of occurences to return true for the bin.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset, typename TInt>
inline void select(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, std::vector<bool> & selected, String<TValue> const & text, TInt threshold)
{
    // me.select(selected, text, static_cast<uint16_t>(threshold));
    me.select(selected, text, threshold);
}

/*!
 * \brief Returns for which bins the counts of all k-mers in a text exceed a threshold.
 * \param me The KmerFilter instance.
 * \param text A single text to count all contained k-mers for.
 * \param threshold The minimal number of occurences to return true for the bin.
 * \returns std::vector<bool> of size binNo indicating whether the text is in the bin.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset, typename TInt>
inline std::vector<bool> select(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, String<TValue> const & text, TInt && threshold)
{
    std::vector<bool> selected(me.noOfBins, false);
    select(me, selected, text, threshold);
    return selected;
}

/*!
 * \brief Returns the number of bins.
 * \param me The KmerFilter instance.
 * \returns Value<KmerFilter<TValue, TSpec> >::Type Number of bins.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline typename Value<KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> >::noOfBins getNumberOfBins(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me)
{
    return me.noOfBins;
}

/*!
 * \brief Returns the k-mer size.
 * \param me The KmerFilter instance.
 * \returns Value<KmerFilter<TValue, TSpec> >::Type k-mer size.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline typename Value<KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> >::kmerSize getKmerSize(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me)
{
    return me.kmerSize;
}

/*!
 * \brief Reads the metadata.
 * \param me The KmerFilter instance.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline void getMetadata(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me)
{
    //-------------------------------------------------------------------
    //| kmer_size | n_hash_func | n_bins |              bf              |
    //-------------------------------------------------------------------
    me.noOfBits = me.filterVector.noOfBits;

    typename Value<KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> >::noOfBits metadataStart = me.filterVector.noOfBits;
    me.noOfBins = me.filterVector.get_int(metadataStart);
    me.noOfHashFunc = me.filterVector.get_int(metadataStart+64);
    me.kmerSize = me.filterVector.get_int(metadataStart+128);
}

/*!
 * \brief Writes the metadata.
 * \param me The KmerFilter instance.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline void setMetadata(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me)
{
    //-------------------------------------------------------------------
    //| kmer_size | n_hash_func | n_bins |              bf              |
    //-------------------------------------------------------------------

    typename Value<KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> >::noOfBits metadataStart = me.noOfBits;
    // TODO also store TValue (alphabet)
    me.filterVector.set_int(metadataStart, me.noOfBins);
    me.filterVector.set_int(metadataStart + 64, me.noOfHashFunc);
    me.filterVector.set_int(metadataStart+128, me.kmerSize);
}

/*!
 * \brief Returns the of the filter vector in MB.
 * \param me The KmerFilter instance.
 * \returns double filter vector size in MB.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline double size(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me)
{
    return me.filterVector.size_in_mega_bytes();
}

/*!
 * \brief Writes the filter vector to a file.
 * \param me The KmerFilter instance.
 * \param fileName Name of the file to write to.
 * \returns bool Indicates if the operation was successful.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline bool store(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, CharString fileName)
{
    setMetadata(me);
    return me.filterVector.store(fileName);
}

/*!
 * \brief Loads the filter vector from a file.
 * \param me The KmerFilter instance.
 * \param fileName Name of the file to read from.
 * \returns bool Indicates if the operation was successful.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset>
inline bool retrieve(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, CharString fileName)
{
    me.filterVector.retrieve(fileName);
    getMetadata(me);
    me.init();
    return true;
}

/*!
 * \brief Indicates whether a bit is set in an integer.
 * \param me The KmerFilter instance.
 * \param num Integer to check.
 * \param bit bit position to check.
 * \returns bool Indicates if the bit is set.
 */
template<typename TValue, typename TSpec, typename TFilterVector, typename TShapeSpec, unsigned offset, typename TInt1, typename TInt2>
inline bool isBitSet(KmerFilter<TValue, TSpec, TFilterVector, TShapeSpec, offset> &  me, TInt1&& num, TInt2&& bit)
{
    return 1 == ( (num >> bit) & 1);
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_KMER_KMER_BASE_H_

