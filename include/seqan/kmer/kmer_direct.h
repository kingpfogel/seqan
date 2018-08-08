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

#ifndef INCLUDE_SEQAN_KMER_KMER_DIRECT_H_
#define INCLUDE_SEQAN_KMER_KMER_DIRECT_H_

// --------------------------------------------------------------------------
// Class KmerFilter using direct addressing
// --------------------------------------------------------------------------
namespace seqan{

/*!
 * \brief Creates and maintains a k-mer directory using direct addressing.
 * Creates a k-mer directory to store occurrence information for k-mers stemming from different bins.
 * An direct addressing k-mer filter represents a collection of direct addressing directories for the individual bins.
 * A k-mer occurs in a direct addressing directory if the position generated by the hash function return is true.
 * Instead of concatenating the individual directories, we interleave them.
 * This results in blocks, where each block represents a hash value and each position in the block corresponds to a
 * bin.
 *
 * \par example
 *
 * ```cpp
 * #include <seqan/kmer.h>
 * CharString file("sequence.fasta");
 * KmerFilter<Dna, DirectAddressing> da (10, 20);
 * insertKmer(da, toCString(file));
 * ```
 *
 */
template<typename TValue, typename TFilterVector, typename TShapeSpec, unsigned offset>
class KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset>
{
public:
    //!\brief The type of the variables.
    typedef String<TValue>                      TString;
    //!\brief The number of Bins.
    typename Value<KmerFilter>::noOfBins        noOfBins;
    //!\brief The k-mer size.
    typename Value<KmerFilter>::kmerSize        kmerSize;
    //!\brief The size of the bit vector.
    typename Value<KmerFilter>::noOfBits        noOfBits;
    //!\brief The number of possible hash values that can fit into a single block.
    typename Value<KmerFilter>::noOfBlocks      noOfBlocks;
    //!\brief The number of 64 bit blocks needed to represent the number of bins.
    typename Value<KmerFilter>::binWidth        binWidth;
    //!\brief Bits we need to represent noBins bits. Multiple of intSize.
    typename Value<KmerFilter>::blockBitSize    blockBitSize;

    //!\brief The bit vector storing the bloom filters.
    FilterVector<TFilterVector>         filterVector;
    //!\brief The shape to be used in the Filter.
    //KmerShape<TValue, TShapeSpec>         shape;
    //!\brief How many bits we can represent in the biggest unsigned int available.
    static const typename Value<KmerFilter>::intSize    intSize{0x40};
    //!\brief Size in bits of the meta data.
    static const typename Value<KmerFilter>::filterMetadataSize filterMetadataSize{256};
    //!\brief The number of used hash functions. Not used but needed for meta data template functions.
    typename Value<KmerFilter>::blockBitSize    noOfHashFunc{1};
    //!\brief A Shape over our filter alphabet.
    typedef Shape<TValue, TShapeSpec>  TShape;
    //!\brief Offset - either CompleteCoverage or IncompleteCoverage
    //KmerOffset<(offset==1u)> kmerOffset;
    /* rule of six */
    /*\name Constructor, destructor and assignment
     * \{
     */
    //!\brief Default constructor
    KmerFilter():
        noOfBins(0),
        kmerSize(0)
        {}

    /*!
     * \brief Constructs direct addressing directory given parameters.
     * \param n_bins Number of bins. Preferably a multiple of 64.
     * \param kmer_size The Size of the k-mer.
     */
    KmerFilter(typename Value<KmerFilter>::noOfBins n_bins, typename Value<KmerFilter>::kmerSize kmer_size):
        noOfBins(n_bins),
        kmerSize(kmer_size),
        filterVector(n_bins, ipow(ValueSize<TValue>::VALUE, kmerSize) * std::ceil((double)noOfBins / intSize) * intSize)
    {
        init();
    }

    //!\brief Copy constructor
    KmerFilter(KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset> & other)
    {
        *this = other;
    }

    //!\brief Copy assignment
    KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset> & operator=(KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset> & other)
    {
        noOfBins = other.noOfBins;
        kmerSize = other.kmerSize;
        noOfBits = other.noOfBits;
        noOfBlocks = other.noOfBlocks;
        filterVector = other.filterVector;
        binWidth = other.binWidth;
        blockBitSize = other.blockBitSize;
        noOfBlocks = other.noOfBlocks;
        return *this;
    }

    //!\brief Move constrcutor
    KmerFilter(KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset> && other)
    {
        *this = std::move(other);
    }

    //!\brief Move assignment
    KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset> & operator=(KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset> && other)
    {
        noOfBins = std::move(other.noOfBins);
        kmerSize = std::move(other.kmerSize);
        noOfBits = std::move(other.noOfBits);
        noOfBlocks = std::move(other.noOfBlocks);
        filterVector = std::move(other.filterVector);
        binWidth = std::move(other.binWidth);
        blockBitSize = std::move(other.blockBitSize);
        noOfBlocks = std::move(other.noOfBlocks);
        return *this;
    }

    //!\brief Destructor
    // ~KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset>() = default;
    ~KmerFilter<TValue, DirectAddressing, TFilterVector, TShapeSpec, offset>() = default;
    //!\}

    /*!
     * \brief Calculates the power of integer x to integer y.
     * \param base Base (integer).
     * \param exp Exponent (integer).
     * \returns uint64_t base^exp
     */
    uint64_t ipow(uint64_t base, uint64_t exp)
    {
        uint64_t result = 1;
        while (exp)
        {
            if (exp & 1)
                result *= base;
            exp >>= 1;
            base *= base;
        }
        return result;
    }

    /*!
     * \brief Resets the bloom filter to 0 for all given bins.
     * \param bins Vector with the ID of the bins to clear.
     * \param threads Number of threads to use.
     */
    template<typename TInt>
    void clear(std::vector<uint32_t> const & bins, TInt&& threads)
    {
        std::vector<std::future<void>> tasks;
        uint64_t chunkBlocks = filterVector.chunkSize / filterVector.blockBitSize;

        for (uint8_t chunk = 0; chunk < filterVector.noOfChunks; ++chunk)
        {
            tasks.clear();
            filterVector.decompress(chunk);

            // We have so many blocks that we want to distribute to so many threads
            uint64_t batchSize = chunkBlocks / threads;
            if(batchSize * threads < chunkBlocks) ++batchSize;

            for (uint8_t taskNo = 0; taskNo < threads; ++taskNo) // TODO Rather divide by chunks?
            {
                // hashBlock is the number of the block the thread will work on. Each block contains binNo bits that
                // represent the individual bins. Each thread has to work on batchSize blocks. We can get the position in
                // our filterVector by multiplying the hashBlock with noOfBins. Then we just need to add the respective
                // binNo. We have to make sure that the vecPos we generate is not out of bounds, only the case in the last
                // thread if the blocks could not be evenly distributed, and that we do not clear a bin that is assigned to
                // another thread.
                tasks.emplace_back(std::async([=] {
                    for (uint64_t hashBlock=taskNo*batchSize;
                        hashBlock < chunkBlocks && hashBlock < (taskNo +1) * batchSize;
                        ++hashBlock)
                    {
                        uint64_t vecPos = hashBlock * filterVector.blockBitSize;
                        uint8_t chunkNo = vecPos / filterVector.chunkSize;
                        for(uint32_t binNo : bins)
                        {
                            if (chunk == chunkNo)
                                filterVector.unset_pos(vecPos + binNo);
                        }
                    }
                }));
            }
            for (auto &&task : tasks)
            {
                task.get();
            }

            filterVector.compress(chunk);
        }
    }
    /*!
    * \brief resize wrapper for any shape that is not Shape<Dna, SimpleShape> - does nothing
    * \param T can be anything, empty function
    */
    template<typename T>
    void resizeShape(T)
    {
    }
    /*!
    * \brief resize wrapper for shape that is Shape<Dna, SimpleShape>
    * \param shape to be resized
    */
    void resizeShape(Shape<Dna, SimpleShape> & shape)
    {
        resize(shape, kmerSize);
    }

    /*!
     * \brief Counts number of occurences in each bin for a given text.
     * \param counts Vector to be filled with counts.
     * \param text Text to count occurences for.
     */
    void select(std::vector<uint32_t> & counts, TString const & text)
    {
        TShape kmerShape;
        uint16_t possible = length(text) - kmerSize + 1;
        // Supports text lengths up to 65535 + k

        uint16_t x = (length(text) - kmerSize) % offset;
        uint16_t noOfKmerHashes = 1 + (length(text) - kmerSize + offset - x) / offset;
        if (x == 0)
        {
            noOfKmerHashes -= 1;
        }
        std::vector<uint64_t> kmerHashes(noOfKmerHashes, 0);
        resizeShape(kmerShape);
        auto it = begin(text);
        hashInit(kmerShape, it);

        uint32_t c = 0;
        for (uint32_t i = 0; i < possible; ++i)
        {
            uint64_t kmerHash = hashNext(kmerShape, it);
            if(i-c*offset == 0)
            {
                if(c < ((length(text) - kmerSize + offset - x) / offset))
                {
                    kmerHashes[c] = kmerHash;
                }
                ++c;
            }
            if(i == possible-1u && x != 0u)
            {

                kmerHashes[c] = kmerHash;
            }
            ++it;
        }

        for (uint64_t kmerHash : kmerHashes)
        {
            // Move to first bit representing the hash kmerHash for bin 0, the next bit would be for bin 1, and so on
            kmerHash *= blockBitSize;

            uint32_t binNo = 0;
            for (uint32_t batchNo = 0; batchNo < binWidth; ++batchNo)
            {
                binNo = batchNo * intSize;
                // get_int(idx, len) returns the integer value of the binary string of length len starting
                // at position idx, i.e. len+idx-1|_______|idx, Vector is right to left.
                uint64_t tmp = filterVector.get_int(kmerHash, intSize);

                // Behaviour for a bit shift with >= maximal size is undefined, i.e. shifting a 64 bit integer by 64
                // positions is not defined and hence we need a special case for this.
                if (tmp ^ (1ULL<<(intSize-1)))
                {
                    // As long as any bit is set
                    while (tmp > 0)
                    {
                        // sdsl::bits::lo calculates the position of the rightmost 1-bit in
                        // the 64bit integer x if it exists.
                        // For example, for 8 = 1000 it would return 3
                        uint64_t step = sdsl::bits::lo(tmp);
                        // Adjust our bins
                        binNo += step;
                        // Remove up to next 1
                        ++step;
                        tmp >>= step;
                        // Count
                        ++counts[binNo];
                        // ++binNo because step is 0-based, e.g., if we had a hit with the next bit we
                        // would otherwise count it for binNo=+ 0
                        ++binNo;
                    }
                }
                else
                {
                    ++counts[binNo + intSize - 1];
                }
                // We will now start with the next batch if possible, so we need to shift the index.
                kmerHash += intSize;
            }
        }
    }

    /*!
     * \brief Tests for occurence in each bin given a text and count threshold.
     * \param selected Vector to be filled with booleans signalling occurence.
     * \param text Text to count occurences for.
     * \param threshold Minimal count (>=) of containing k-mers to report bin as containing text.
     */
    template<typename TInt>
    inline void select(std::vector<bool> & selected, TString const & text, TInt && threshold)
    {
        std::vector<uint32_t> counts(noOfBins, 0);
        select(counts, text);
        for(uint32_t binNo=0; binNo < noOfBins; ++binNo)
        {
            if(counts[binNo] >= threshold)
                selected[binNo] = true;
        }
    }

    /*!
     * \brief Adds all k-mers from a text to the IBF.
     * \param text Text to process.
     * \param binNo bin ID to insertKmer k-mers in.
     */
    template<typename TBin, typename TChunk>
    inline void insertKmer(TString const & text,TBin && binNo, TChunk && chunkNo)
    {
        TShape kmerShape;
        //shape.resizeShape(kmerSize);
        resizeShape(kmerShape);
        auto it = begin(text);
        hashInit(kmerShape, it);
        for (uint64_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
        {
            uint64_t kmerHash = hashNext(kmerShape, it);

            ++it;
            uint64_t vecIndex = kmerHash * blockBitSize + binNo;
            filterVector.set_pos(vecIndex, chunkNo);
        }
    }

    //! \brief Initialises internal variables.
    inline void init()
    {
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((double)noOfBins / intSize);       
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * intSize;
       
         // How many hash values must we represent
        noOfBlocks = ipow(ValueSize<TValue>::VALUE, kmerSize);
        // Size of the bit vector
        noOfBits = noOfBlocks * blockBitSize;
        // filterVector = FilterVector<TFilterVector>(noOfBins, noOfBits);
    }
};
}

#endif  // INCLUDE_SEQAN_KMER_KMER_DIRECT_H_
