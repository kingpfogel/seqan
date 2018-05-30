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
// ==========================================================================

#ifndef INCLUDE_SEQAN_KMER_KMER_IBF_H_
#define INCLUDE_SEQAN_KMER_KMER_IBF_H_
#include <seqan/index.h>
// --------------------------------------------------------------------------
// Class KmerFilter using an interleaved bloom filter
// --------------------------------------------------------------------------
namespace seqan{

/*!
 * \brief Creates and maintains a k-mer directory using an interleaved bloom filter (IBF).
 * Creates a k-mer directory to store occurrence information for k-mers stemming from different bins.
 * An IBF represents a collection of bloom filters for the individual bins.
 * A k-mer occurs in a bloom filter if all position generated by the hash functions return a 1.
 * A bloom filter may return false positives, but no false negatives.
 * Instead of concatenating the individual bloom filters, we interleave them.
 * This results in blocks, where each block represents a hash value and each position in the block corresponds to a
 * bin.
 *
 * \par example
 *
 * ```cpp
 * #include <seqan/kmer.h>
 * CharString file("sequence.fasta");
 * KmerFilter<Dna, InterleavedBloomFilter> ibf (10, 3, 20, 16777472);
 * insertKmer(ibf, toCString(file));
 * ```
 *
 */
template<typename TValue, typename TFilterVector, typename TSpec2, typename TSelector>
class KmerFilter<TValue, InterleavedBloomFilter, TFilterVector, TSpec2, TSelector>
{
public:
    //!\brief The type of the variables.
    typedef String<TValue>                      TString;
    //!\brief The number of Bins.
    typename Value<KmerFilter>::noOfBins        noOfBins;
    //!\brief The number of hash functions.
    typename Value<KmerFilter>::noOfHashFunc    noOfHashFunc;
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

    //!\brief Randomized values for hash functions.
    std::vector<typename Value<KmerFilter>::preCalcValues>   preCalcValues;
    //!\brief Shift value used for hash function.
    static const typename Value<KmerFilter>::shiftValue      shiftValue{27};
    //!\brief Random seed.
    static const typename Value<KmerFilter>::seedValue       seedValue{0x90b45d39fb6da1fa};
    //!\brief How many bits we can represent in the biggest unsigned int available.
    static const typename Value<KmerFilter>::intSize         intSize{0x40};
    //!\brief The bit vector storing the bloom filters.
    FilterVector<TFilterVector>                              filterVector;
    //!\brief The shape to be used in the Filter.
    //KmerShape<TValue, TSpec2>                                shape;
    //!\brief Size in bits of the meta data.
    static const typename Value<KmerFilter>::filterMetadataSize filterMetadataSize{256};
    //!\brief A ungapped Shape over our filter alphabet.
    typedef Shape<TValue, TSpec2>  TShape;
    typedef TSelector Selector;
    /* rule of six */
    /*\name Constructor, destructor and assignment
     * \{
     */
    //!\brief Default constructor
    KmerFilter():
        noOfBins(0),
        noOfHashFunc(0),
        kmerSize(0),
        noOfBits(0),
        filterVector() {}

    /*!
     * \brief Constructs an IBF given parameters.
     * \param n_bins Number of bins. Preferably a multiple of 64.
     * \param n_hash_func Number of hash functions.
     * \param kmer_size The Size of the k-mer.
     * \param vec_size The size of the bit vector in bits. Preferably a power of two + 256 for metadata.
     */
    KmerFilter(typename Value<KmerFilter>::noOfBins n_bins, typename Value<KmerFilter>::noOfHashFunc n_hash_func, typename Value<KmerFilter>::kmerSize kmer_size, typename Value<KmerFilter>::noOfBits vec_size):
        noOfBins(n_bins),
        noOfHashFunc(n_hash_func),
        kmerSize(kmer_size),
        noOfBits(vec_size),
        filterVector(noOfBins, noOfBits)
    {
        init();
    }

    //!\brief Copy constructor
    KmerFilter(KmerFilter<TValue, InterleavedBloomFilter, TFilterVector, TSpec2, TSelector> & other)
    {
        *this = other;
    }

    //!\brief Copy assignment
    KmerFilter<TValue, InterleavedBloomFilter, TFilterVector, TSpec2, TSelector> & operator=(KmerFilter<TValue, InterleavedBloomFilter, TFilterVector, TSpec2, TSelector> & other)
    {
        noOfBins = other.noOfBins;
        noOfHashFunc = other.noOfHashFunc;
        kmerSize = other.kmerSize;
        noOfBits = other.noOfBits;
        filterVector = other.filterVector;
        init();
        return *this;
    }

    //!\brief Move constrcutor
    KmerFilter(KmerFilter<TValue, InterleavedBloomFilter, TFilterVector, TSpec2, TSelector> && other)
    {
        *this = std::move(other);
    }

    //!\brief Move assignment
    KmerFilter<TValue, InterleavedBloomFilter, TFilterVector, TSpec2, TSelector> & operator=(KmerFilter<TValue, InterleavedBloomFilter, TFilterVector, TSpec2, TSelector> && other)
    {
        noOfBins = std::move(other.noOfBins);
        noOfHashFunc = std::move(other.noOfHashFunc);
        kmerSize = std::move(other.kmerSize);
        noOfBits = std::move(other.noOfBits);
        filterVector = std::move(other.filterVector);
        init();
        return *this;
    }

    //!\brief Destructor
    ~KmerFilter<TValue, InterleavedBloomFilter, TFilterVector, TSpec2, TSelector>() = default;
    //!\}

    /*!
     * \brief Resets the bloom filter to 0 for all given bins.
     * \param bins Vector with the ID of the bins to clear.
     * \param threads Number of threads to use.
     */
    template<typename TInt>
    void clear(std::vector<uint16_t> const & bins, TInt&& threads)
    {
        std::vector<std::future<void>> tasks;
        uint64_t chunkBlocks = filterVector.chunkSize / filterVector.blockBitSize;

        for (uint8_t chunk = 0; chunk < filterVector.noOfChunks; ++chunk)
        {
            tasks.clear();
            filterVector.decompress(chunk);

            // We have so many blocks that we want to distribute to so many threads
            uint64_t batchSize = noOfBlocks / threads;
            if(batchSize * threads < noOfBlocks) ++batchSize;

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
                        uint8_t  chunkNo = vecPos / filterVector.chunkSize;
                        for(uint16_t binNo : bins)
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
    template<typename TShapeSpec>
    void resizeShape(TShapeSpec)
    {
        //std::cout << "resize arbitrary shape" << std::endl;

    }
    void resizeShape(Shape<Dna, SimpleShape> & shape)
    {
       // std::cout << "resize SimpleShape " << kmerSize<< std::endl;
        resize(shape, kmerSize);
    }
    inline std::vector<uint64_t> selectHelper(OverlappingKmers const &, TString const & text, TShape kmerShape)
    {
        uint16_t possible = length(text) - kmerSize + 1; // Supports text lengths up to 65535 + k
        std::vector<uint64_t> kmerHashes(possible, 0);

        //TShape kmerShape;
        //shape.resizeShape(kmerSize);
        resizeShape(kmerShape);
        hashInit(kmerShape, begin(text));
        auto it = begin(text);
        for (uint16_t i = 0; i < possible; ++i)
        {
            kmerHashes[i] = hashNext(kmerShape, it);
            ++it;
        }
        return kmerHashes;
    }
    inline std::vector<uint64_t> selectHelper(NonOverlappingKmers const &, TString const & text, TShape kmerShape)
    {
        uint16_t possible = length(text) - kmerSize + 1; // Supports text lengths up to 65535 + k

        uint16_t x = length(text) % kmerSize;
        uint16_t noOfKmerHashes = ((length(text) - x) / kmerSize)+1;
        if (x == 0)
        {
            noOfKmerHashes -= 1;
        }
        std::vector<uint64_t> kmerHashes(noOfKmerHashes, 0);
        resizeShape(kmerShape);
        hashInit(kmerShape, begin(text));
        auto it = begin(text);
        uint32_t c = 0;
        for (uint32_t i = 0; i < possible; ++i)
        {
            uint64_t kmerHash = hashNext(kmerShape, it);
            if(i-c*kmerSize == 0)
            {
                if(c < ((length(text) - x) / kmerSize))
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
        return kmerHashes;
    }
    /*!
     * \brief Counts number of occurences in each bin for a given text.
     * \param counts Vector to be filled with counts.
     * \param text Text to count occurences for.
     */
    void select(std::vector<uint16_t> & counts, TString const & text) // TODO uint16_t
    {
        TShape kmerShape;
        Selector sel;
        std::vector<uint64_t> kmerHashes = selectHelper(sel, text, kmerShape);

        for (uint64_t kmerHash : kmerHashes)
        {
            std::vector<uint64_t> vecIndices = preCalcValues;
            for(uint8_t i = 0; i < noOfHashFunc ; ++i)
            {
                vecIndices[i] *= kmerHash;
                hashToIndex(vecIndices[i]);
            }

            uint16_t binNo = 0;
            for (uint16_t batchNo = 0; batchNo < binWidth; ++batchNo)
            {
                binNo = batchNo * intSize;
                // get_int(idx, len) returns the integer value of the binary string of length len starting
                // at position idx, i.e. len+idx-1|_______|idx, Vector is right to left.
                uint64_t tmp = filterVector.get_int(vecIndices[0], intSize);

                // A k-mer is in a bin of the IBF iff all hash functions return 1 for the bin.
                for(uint8_t i = 1; i < noOfHashFunc;  ++i)
                {
                    tmp &= filterVector.get_int(vecIndices[i], intSize);
                }

                // Behaviour for a bit shift with >= maximal size is undefined, i.e. shifting a 64 bit integer by 64
                // positions is not defined and hence we need a special case for this.
                if (tmp ^ (1ULL<<(intSize-1)))
                {
                    // As long as any bit is set
                    while (tmp > 0)
                    {
                        // sdsl::bits::lo calculates the position of the rightmost 1-bit in the 64bit integer x if it exists.
                        // For example, for 8 = 1000 it would return 3
                        uint64_t step = sdsl::bits::lo(tmp);
                        // Adjust our bins
                        binNo += step;
                        // Remove up to next 1
                        ++step;
                        tmp >>= step;
                        // Count
                        ++counts[binNo];
                        // ++binNo because step is 0-based, e.g., if we had a hit with the next bit we would otherwise count it for binNo=+ 0
                        ++binNo;
                    }
                }
                else
                {
                    ++counts[binNo + intSize - 1];
                }
                // We will now start with the next batch if possible, so we need to shift all the indices.
                for (uint8_t i = 0; i < noOfHashFunc; ++i)
                {
                    vecIndices[i] += intSize;
                }
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
        std::vector<uint16_t> counts(noOfBins, 0);
        select(counts, text);
        for(uint16_t binNo=0; binNo < noOfBins; ++binNo)
        {
            if(counts[binNo] >= threshold)
                selected[binNo] = true;
        }
    }

    /*!
     * \brief Calculates the first index of a block that corresponds to a hash value.
     * \param hash hash value
     */
    inline void hashToIndex(uint64_t & hash) const
    {
        // We do something
        hash ^= hash >> shiftValue;
        // Bring it back into our vector range (noOfBlocks = possible hash values)
        hash %= noOfBlocks;
        // hash &= (noOfBlocks - 1);
        // Since each block needs blockBitSize bits, we multiply to get the correct location
        hash *= blockBitSize;
    }

    /*!
     * \brief Adds all k-mers from a text to the IBF.
     * \param text Text to process.
     * \param binNo bin ID to insertKmer k-mers in.
     */
    template<typename TBin, typename TChunk>
    inline void insertKmer(TString const & text, TBin && binNo, TChunk && chunkNo)
    {

        TShape kmerShape;
        resizeShape(kmerShape);
        //resize(kmerShape, kmerSize);
        hashInit(kmerShape, begin(text));

        for (uint64_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
        {
            uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);

            for(uint8_t i = 0; i < noOfHashFunc ; ++i)
            {
                uint64_t vecIndex = preCalcValues[i] * kmerHash;
                hashToIndex(vecIndex);
                vecIndex += binNo;
                filterVector.set_pos(vecIndex, chunkNo);
            }
        }
    }

    //! \brief Initialises internal variables.
    inline void init()
    {
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((double)noOfBins / intSize);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * intSize;
        // How many hash values can we represent
        noOfBlocks = (noOfBits - filterMetadataSize) / blockBitSize;

        preCalcValues.resize(noOfHashFunc);
        for(uint8_t i = 0; i < noOfHashFunc ; i++)
            preCalcValues[i] = i ^  (kmerSize * seedValue);
    }
};
}

#endif  // INCLUDE_SEQAN_KMER_KMER_IBF_H_
