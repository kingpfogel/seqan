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

// --------------------------------------------------------------------------
// Class KmerFilter using an interleaved bloom filter
// --------------------------------------------------------------------------
namespace seqan{

// struct Distributor
// {
//     static const uint64_t MAX_MEM = 1ULL<<35; // 4 GiB
//     static const uint64_t MAX_VEC = 1ULL<<32; // 0.5 GiB
//     static const uint64_t MAX_NUM = MAX_MEM / MAX_VEC; // 8
//     static const uint64_t FILTER_METADATA_SIZE = 256;
//     static const uint64_t INT_SIZE = 0x40; // 64
//
//     uint64_t noOfBins;
//     uint64_t noOfBits;
//     uint64_t binWidth;
//     uint64_t blockBitSize;
//     uint64_t noOfBlocks;
//     uint64_t noOfChunks;
//     uint64_t chunkSize;
//
//     std::list<uint64_t> queue;
//
//     void decompress(uint64_t chunk)
//     {
//         if (queue.size() < MAX_NUM)
//         {
//             std::get<1>(filterVector[chunk]) = static_cast<sdsl::bit_vector>(std::get<2>(filterVector[chunk]));
//             std::get<0>(filterVector[chunk]) = false;
//             queue.emplace_back(chunk);
//         }
//         else
//         {
//             compress();
//             decompress(chunk);
//         }
//     }
//
//     void compress()
//     {
//         uint64_t chunk = queue.front();
//         std::get<2>(filterVector[chunk]) = sdsl::sd_vector<>(std::get<1>(filterVector[chunk]));
//         std::get<1>(filterVector[chunk]) = sdsl::bit_vector();
//         std::get<0>(filterVector[chunk]) = true;
//         queue.pop_front();
//     }
//
//     std::vector<std::tuple<bool,sdsl::bit_vector, sdsl::sd_vector<> > > filterVector;
//
//     Distributor() {}
//
//     Distributor(uint64_t bins, uint64_t bits):
//         noOfBins(bins),
//         noOfBits(bits)
//     {
//         // How many blocks of 64 bit do we need to represent our noOfBins
//         binWidth = std::ceil((float)noOfBins / INT_SIZE);
//         // How big is then a block (multiple of 64 bit)
//         blockBitSize = binWidth * INT_SIZE;
//         // How many hash values can we represent
//         noOfBlocks = (noOfBits - FILTER_METADATA_SIZE) / blockBitSize;
//
//         // We need to split at the end of a block.
//         chunkSize = std::ceil((float)MAX_VEC/blockBitSize) * blockBitSize;
//         // This is how many chunks we need.
//         uint64_t noChunks = std::ceil((float) noOfBits/chunkSize);
//         for (uint64_t i = 0; i < noChunks; ++i)
//         {
//             filterVector.push_back(std::make_tuple(true,
//                                                    sdsl::bit_vector(),
//                                                    sdsl::sd_vector<>(sdsl::bit_vector(chunkSize, 0))));
//         }
//     }
//
//     auto get_int(uint64_t idx, uint64_t len)
//     {
//         uint64_t access = idx;
//         uint64_t chunkNo = access / MAX_VEC;
//         uint64_t chunkPos = access - chunkNo * MAX_VEC;
//         if (std::get<0>(filterVector[chunkNo]))
//         {
//             return std::get<2>(filterVector[chunkNo]).get_int(chunkPos, len);
//         }
//         else
//         {
//             return std::get<1>(filterVector[chunkNo]).get_int(chunkPos, len);
//         }
//     }
//
//     uint64_t get_pos(uint64_t vecIndex)
//     {
//         uint64_t access = vecIndex;
//         uint64_t chunkNo = access / MAX_VEC;
//         uint64_t chunkPos = access - chunkNo * MAX_VEC;
//         if (std::get<0>(filterVector[chunkNo]))
//         {
//             return std::get<2>(filterVector[chunkNo])[chunkPos];
//         }
//         else
//         {
//             return std::get<1>(filterVector[chunkNo])[chunkPos];
//         }
//     }
//
//     void set_pos(uint64_t vecIndex)
//     {
//         set_impl(vecIndex, true);
//     }
//
//     void unset_pos(uint64_t vecIndex)
//     {
//         set_impl(vecIndex, false);
//     }
//
//     void set_impl(uint64_t vecIndex, bool value)
//     {
//         uint64_t access = vecIndex;
//         uint64_t chunkNo = access / MAX_VEC;
//         uint64_t chunkPos = access - chunkNo * MAX_VEC;
//
//         if (std::get<0>(filterVector[chunkNo]))
//         {
//             decompress(chunkNo);
//             std::get<1>(filterVector[chunkNo])[chunkPos] = value;
//         }
//         else
//         {
//             std::get<1>(filterVector[chunkNo])[chunkPos] = value;
//         }
//     }
// };
struct Distributor
{
    static const uint64_t MAX_MEM = 1ULL<<35; // 4 GiB
    static const uint64_t MAX_VEC = 1ULL<<32; // 0.5 GiB
    static const uint64_t MAX_NUM = MAX_MEM / MAX_VEC; // 8
    static const uint64_t FILTER_METADATA_SIZE = 256;
    static const uint64_t INT_SIZE = 0x40; // 64
    CharString PREFIX{"test"};

    uint64_t noOfBins;
    uint64_t noOfBits;
    uint64_t binWidth;
    uint64_t blockBitSize;
    uint64_t noOfBlocks;
    uint64_t noOfChunks;
    uint64_t chunkSize;

    std::list<uint64_t> queue;

    inline bool size_check()
    {
        return queue.size() < MAX_NUM;
    }

    void decompress(uint64_t chunk)
    {
        if (size_check())
        {
            sdsl::load_from_file(*std::get<1>(filterVector[chunk]), toCString(PREFIX)+std::to_string(chunk));
            std::get<0>(filterVector[chunk]) = false;
            queue.emplace_back(chunk);
        }
        else
        {
            compress();
            decompress(chunk);
        }
    }

    void compress()
    {
        uint64_t chunk = queue.front();
        std::get<2>(filterVector[chunk]) = std::make_unique<sdsl::sd_vector<> >(*std::get<1>(filterVector[chunk]));
        sdsl::store_to_file(*std::get<1>(filterVector[chunk]), toCString(PREFIX)+std::to_string(chunk));
        std::get<1>(filterVector[chunk]) = std::make_unique<sdsl::bit_vector>(0,0);
        std::get<0>(filterVector[chunk]) = true;
        queue.pop_front();
    }

    std::vector<std::tuple<bool,std::unique_ptr<sdsl::bit_vector>, std::unique_ptr<sdsl::sd_vector<> > > > filterVector;

    Distributor() {}

    Distributor(uint64_t bins, uint64_t bits):
        noOfBins(bins),
        noOfBits(bits)
    {
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((float)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = (noOfBits - FILTER_METADATA_SIZE) / blockBitSize;

        // We need to split at the end of a block.
        chunkSize = std::ceil((float)MAX_VEC/blockBitSize) * blockBitSize;
        // This is how many chunks we need.
        uint64_t noChunks = std::ceil((float) noOfBits/chunkSize);
        for (uint64_t i = 0; i < noChunks; ++i)
        {
            if (size_check())
            {
                filterVector.emplace_back(std::make_tuple(false,
                                                       std::make_unique<sdsl::bit_vector>(chunkSize, 0),
                                                       std::make_unique<sdsl::sd_vector<> >()));
                queue.emplace_back(i);
            }
            else
            {
                compress();
                filterVector.emplace_back(std::make_tuple(false,
                                                       std::make_unique<sdsl::bit_vector>(chunkSize, 0),
                                                       std::make_unique<sdsl::sd_vector<> >()));
                queue.emplace_back(i);
            }
        }
        while (!queue.empty())
        {
            compress();
        }
    }

    auto get_int(uint64_t idx, uint64_t len)
    {
        uint64_t access = idx;
        uint64_t chunkNo = access / MAX_VEC;
        uint64_t chunkPos = access - chunkNo * MAX_VEC;
        if (std::get<0>(filterVector[chunkNo]))
        {
            return std::get<2>(filterVector[chunkNo])->get_int(chunkPos, len);
        }
        else
        {
            return std::get<1>(filterVector[chunkNo])->get_int(chunkPos, len);
        }
    }

    uint64_t get_pos(uint64_t vecIndex)
    {
        uint64_t access = vecIndex;
        uint64_t chunkNo = access / MAX_VEC;
        uint64_t chunkPos = access - chunkNo * MAX_VEC;
        if (std::get<0>(filterVector[chunkNo]))
        {
            return (*std::get<2>(filterVector[chunkNo]))[chunkPos];
        }
        else
        {
            return (*std::get<1>(filterVector[chunkNo]))[chunkPos];
        }
    }

    void set_pos(uint64_t vecIndex)
    {
        set_impl(vecIndex, true);
    }

    void unset_pos(uint64_t vecIndex)
    {
        set_impl(vecIndex, false);
    }

    void set_impl(uint64_t vecIndex, bool value)
    {
        uint64_t access = vecIndex;
        uint64_t chunkNo = access / chunkSize;
        uint64_t chunkPos = access - chunkNo * chunkSize;

        if (std::get<0>(filterVector[chunkNo]))
        {
            decompress(chunkNo);
            (*std::get<1>(filterVector[chunkNo]))[chunkPos] = value;
        }
        else
        {
            (*std::get<1>(filterVector[chunkNo]))[chunkPos] = value;
        }
    }
};

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
 * addFastaFile(ibf, toCString(file));
 * ```
 *
 */
template<typename TValue>
class KmerFilter<TValue, InterleavedBloomFilter>
{
public:
    //!\brief The type of the variables.
    typedef typename Value<KmerFilter>::Type    THValue;
    // limit for splitting vector
    THValue    limit;
    //!\brief The number of Bins.
    THValue    noOfBins;
    //!\brief The number of hash functions.
    THValue    noOfHashFunc;
    //!\brief The k-mer size.
    THValue    kmerSize;
    //!\brief The size of the bit vector.
    THValue    noOfBits;
    //!\brief The number of possible hash values that can fit into a single block.
    THValue    noOfBlocks;
    //!\brief The number of 64 bit blocks needed to represent the number of bins.
    THValue    binWidth;
    //!\brief Bits we need to represent noBins bits. Multiple of intSize.
    THValue    blockBitSize;

    Distributor filterVector;

    //!\brief Randomized values for hash functions.
    std::vector<THValue>   preCalcValues;
    //!\brief Shift value used for hash function.
    static const THValue   shiftValue = 27;
    //!\brief Random seed.
    static const THValue   seedValue = 0x90b45d39fb6da1fa;
    //!\brief How many bits we can represent in the biggest unsigned int available.
    static const THValue   intSize = 0x40;
    //!\brief Maximum size of uncompressed bitvector
    static const uint64_t  SIZE_LIMIT = 1ULL<<32;
    //!\brief A vector holding the bit vectors.
    // std::vector<sdsl::bit_vector> filterVector;
    //!\brief Size in bits of the meta data.
    static const uint32_t               filterMetadataSize{256};
    //!\brief A ungapped Shape over our filter alphabet.
    typedef Shape<TValue, SimpleShape>  TShape;

    /* rule of six */
    /*\name Constructor, destructor and assignment
     * \{
     */
    //!\brief Default constructor
    KmerFilter():
        noOfBins(0),
        noOfHashFunc(0),
        kmerSize(0),
        noOfBits(0)
        {}

    /*!
     * \brief Constructs an IBF given parameters.
     * \param n_bins Number of bins. Preferably a multiple of 64.
     * \param n_hash_func Number of hash functions.
     * \param kmer_size The Size of the k-mer.
     * \param vec_size The size of the bit vector in bits. Preferably a power of two + 256 for metadata.
     */
    KmerFilter(THValue n_bins, THValue n_hash_func, THValue kmer_size, THValue vec_size):
        noOfBins(n_bins),
        noOfHashFunc(n_hash_func),
        kmerSize(kmer_size),
        noOfBits(vec_size)
    {
        init();
    }

    //!\brief Copy constructor
    KmerFilter(KmerFilter<TValue, InterleavedBloomFilter> const & other)
    {
        *this = other;
    }

    //!\brief Copy assignment
    KmerFilter<TValue, InterleavedBloomFilter> & operator=(KmerFilter<TValue, InterleavedBloomFilter> const & other)
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
    KmerFilter(KmerFilter<TValue, InterleavedBloomFilter> && other)
    {
        *this = std::move(other);
    }

    //!\brief Move assignment
    KmerFilter<TValue, InterleavedBloomFilter> & operator=(KmerFilter<TValue, InterleavedBloomFilter> && other)
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
    ~KmerFilter<TValue, InterleavedBloomFilter>() = default;
    //!\}

    /*!
     * \brief Resets the bloom filter to 0 for all given bins.
     * \param bins Vector with the ID of the bins to clear.
     * \param threads Number of threads to use.
     */
    template<typename TInt>
    void clearBins(std::vector<THValue> const & bins, TInt&& threads)
    {
        std::vector<std::future<void>> tasks;

        // We have so many blocks that we want to distribute to so many threads
        uint64_t batchSize = noOfBlocks / threads;
        if(batchSize * threads < noOfBlocks) ++batchSize;

        for (uint32_t taskNo = 0; taskNo < threads; ++taskNo)
        {
            // hashBlock is the number of the block the thread will work on. Each block contains binNo bits that
            // represent the individual bins. Each thread has to work on batchSize blocks. We can get the position in
            // our filterVector by multiplying the hashBlock with noOfBins. Then we just need to add the respective
            // binNo. We have to make sure that the vecPos we generate is not out of bounds, only the case in the last
            // thread if the blocks could not be evenly distributed, and that we do not clear a bin that is assigned to
            // another thread.
            tasks.emplace_back(std::async([=] {
                for (uint64_t hashBlock=taskNo*batchSize;
                    hashBlock < noOfBlocks && hashBlock < (taskNo +1) * batchSize;
                    ++hashBlock)
                {
                    uint64_t vecPos = hashBlock * blockBitSize;
                    for(uint32_t binNo : bins)
                    {
                        filterVector.unset_pos(vecPos + binNo);
                    }
                }
            }));
        }
        for (auto &&task : tasks)
        {
            task.get();
        }
    }

    /*!
     * \brief Counts number of occurences in each bin for a given text.
     * \param counts Vector to be filled with counts.
     * \param text Text to count occurences for.
     */
    template<typename TString>
    void whichBins(std::vector<uint64_t> & counts, TString const & text)
    {
        uint8_t possible = length(text) - kmerSize + 1;

        std::vector<uint64_t> kmerHashes(possible, 0);

        TShape kmerShape;
        resize(kmerShape, kmerSize);
        hashInit(kmerShape, begin(text));
        auto it = begin(text);
        for (uint32_t i = 0; i < possible; ++i)
        {
            kmerHashes[i] = hashNext(kmerShape, it);
            ++it;
        }

        for (uint64_t kmerHash : kmerHashes)
        {
            std::vector<uint64_t> vecIndices = preCalcValues;
            for(uint8_t i = 0; i < noOfHashFunc ; ++i)
            {
                vecIndices[i] *= kmerHash;
                hashToIndex(vecIndices[i]);
            }

            uint64_t binNo = 0;
            for (uint64_t batchNo = 0; batchNo < binWidth; ++batchNo)
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
                for (uint64_t i = 0; i < noOfHashFunc; ++i)
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
    template<typename TString, typename TInt>
    inline void whichBins(std::vector<bool> & selected, TString const & text, TInt && threshold)
    {
        std::vector<uint64_t> counts(noOfBins, 0);
        whichBins(counts, text);
        for(uint32_t binNo=0; binNo < noOfBins; ++binNo)
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
        // Since each block needs blockBitSize bits, we multiply to get the correct location
        hash *= blockBitSize;
    }

    /*!
     * \brief Adds all k-mers from a text to the IBF.
     * \param text Text to process.
     * \param binNo bin ID to insert k-mers in.
     */
    template<typename TString, typename TInt>
    inline void addKmer(TString const & text, TInt && binNo)
    {

        TShape kmerShape;
        resize(kmerShape, kmerSize);
        hashInit(kmerShape, begin(text));

        for (uint64_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
        {
            uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);

            for(uint8_t i = 0; i < noOfHashFunc ; i++)
            {
                uint64_t vecIndex = preCalcValues[i] * kmerHash;
                hashToIndex(vecIndex);
                vecIndex += binNo;
                filterVector.set_pos(vecIndex);
            }
        }
    }

    //! \brief Initialises internal variables.
    inline void init()
    {
        // // How many blocks of 64 bit do we need to represent our noOfBins
        // binWidth = std::ceil((float)noOfBins / intSize);
        // // How big is then a block (multiple of 64 bit)
        // blockBitSize = binWidth * intSize;
        // // How many hash values can we represent
        // noOfBlocks = (noOfBits - filterMetadataSize) / blockBitSize;
        //
        // // We need to split at the end of a block.
        // limit = std::ceil((float)SIZE_LIMIT/blockBitSize) * blockBitSize;
        // // This is how many chunks we need.
        // uint64_t noChunks = std::ceil((float) noOfBits/limit);
        // for (uint64_t i = 0; i < noChunks; ++i)
        // {
        //     filterVector.push_back(sdsl::bit_vector(limit, 0));
        // }

        filterVector = Distributor(noOfBins, noOfBits);
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((float)noOfBins / intSize);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * intSize;
        // How many hash values can we represent
        noOfBlocks = (noOfBits - filterMetadataSize) / blockBitSize;

        preCalcValues.resize(noOfHashFunc);
        for(uint64_t i = 0; i < noOfHashFunc ; i++)
            preCalcValues[i] = i ^  (kmerSize * seedValue);
    }
};
}

#endif  // INCLUDE_SEQAN_KMER_KMER_IBF_H_
