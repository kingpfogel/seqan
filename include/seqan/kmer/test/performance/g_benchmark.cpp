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
#include <random>
#include <benchmark/benchmark.h>
#include <seqan/kmer.h>
#include <atomic>
#include <seqan/index.h>
using namespace seqan;

// CharString baseDir{"/srv/public/enricoseiler/benchmark/"}; // lncrna
CharString baseDir{"/home/pfogel/devel/testenv_exp/"};
uint64_t e{2};
uint64_t readNo{1000000};
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
template <typename TAlphabet, typename TFilter, typename TSpec2>
static void insertKmer_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter, TSpec2> ibf (bins, hash, k, (1ULL<<bits));

    for (auto _ : state)
    {
        double elapsed_seconds{0.0};
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/bins/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fasta"));

            auto start = std::chrono::high_resolution_clock::now();
            insertKmer(ibf, toCString(file), i);
            auto end   = std::chrono::high_resolution_clock::now();
            elapsed_seconds += (std::chrono::duration_cast<std::chrono::duration<double> >(end - start)).count();
            std::cerr << "IBF Bin " << i << " done." << '\n';
        }
        state.SetIterationTime(elapsed_seconds);
    }

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bits)));

    append(storage, CharString("_ibf_A1.filter"));
    store(ibf, storage);

    state.counters["Size"] = ibf.filterVector.size_in_mega_bytes();
}


template <typename TAlphabet, typename TFilter, typename TSpec2>
static void select_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    //uint16_t t = std::floor(100/k)-e;
    uint16_t t = 100-k+1 - k*e;
    KmerFilter<TAlphabet, InterleavedBloomFilter/*NoOverlaps*/, TFilter, TSpec2> ibf(bins, hash, k, (1ULL<<bits));

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bits)));

    append(storage, CharString("_ibf.filter"));
    retrieve(ibf, storage);

    double verifications{0};
    double tp{0};
    double p{0};
    double fp{0};
    double fn{0};
    double readNo{0};

    for (auto _ : state)
    {
        double elapsed_seconds{0.0};
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fastq"));

            CharString id;
            String<TAlphabet> seq;
            SeqFileIn seqFileIn;
            uint64_t c{0};
            if (!open(seqFileIn, toCString(file)))
            {
                CharString msg = "Unable to open contigs file: ";
                append(msg, CharString(file));
                throw toCString(msg);
            }
            while(!atEnd(seqFileIn))
            {
                readRecord(id, seq, seqFileIn);
                auto start = std::chrono::high_resolution_clock::now();
                auto res = select(ibf, seq, t);
                auto end   = std::chrono::high_resolution_clock::now();
                ++readNo;
                elapsed_seconds += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
                if (res[i])
                    ++tp;
                else
                    ++fn;
                c = count(res.begin(), res.end(), true);
                verifications += c;
                if (c > 1)
                {
                    if (res[i])
                        fp += c - 1;
                    else
                        fp += c;
                }
                p += c;
            }
        }
        state.SetIterationTime(elapsed_seconds);
        state.counters["5_TP"] = tp;
        state.counters["6_FN"] = fn;
        state.counters["7_FP"] = fp;
        state.counters["8_P"] = p;
        state.counters["99_readNo"] = readNo;
        state.counters["9_verifications"] = verifications;
        state.counters["0_Verifications"] = verifications/readNo;
        state.counters["1_Sensitivity"] = tp/readNo;
        state.counters["2_Precision"] = tp/p;
        state.counters["3_FNR"] = fn/readNo;
        state.counters["4_FDR"] = fp/p;
    }
}

template <typename TAlphabet, typename TFilter, typename TSpec2>
static void insertKmer_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressingNoOverlaps, TFilter, TSpec2> da (bins, k);

    for (auto _ : state)
    {
        double elapsed_seconds{0.0};
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/bins/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fasta"));

            auto start = std::chrono::high_resolution_clock::now();
            insertKmer(da, toCString(file), i);
            auto end   = std::chrono::high_resolution_clock::now();
            elapsed_seconds += (std::chrono::duration_cast<std::chrono::duration<double> >(end - start)).count();
            std::cerr << "DA Bin " << i << " done." << '\n';
        }
        state.SetIterationTime(elapsed_seconds);
    }

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));

    append(storage, CharString("_da.filter"));
    store(da, storage);

    state.counters["Size"] = da.filterVector.size_in_mega_bytes();
}

template <typename TAlphabet, typename TFilter>
static void select_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing, TFilter> da (bins, k);

    uint64_t readNo{1};

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));

    append(storage, CharString("_da.filter"));
    retrieve(da, storage);
    da.filterVector.compress(0);

    uint64_t verifications{0};
    uint64_t tp{0};

    for (auto _ : state)
    {
        double elapsed_seconds{0.0};
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fastq"));

            CharString id;
            String<TAlphabet> seq;
            SeqFileIn seqFileIn;
            if (!open(seqFileIn, toCString(file)))
            {
                CharString msg = "Unable to open contigs file: ";
                append(msg, CharString(file));
                throw toCString(msg);
            }
            while(!atEnd(seqFileIn))
            {
                readRecord(id, seq, seqFileIn);
                auto start = std::chrono::high_resolution_clock::now();
                auto res = select(da, seq, 100-k+1 - k*e);
                auto end   = std::chrono::high_resolution_clock::now();
                elapsed_seconds += (std::chrono::duration_cast<std::chrono::duration<double> >(end - start)).count();

                if (res[i])
                {
                    ++tp;
                }
                verifications += count(res.begin(), res.end(), true);
            }
        }
        state.SetIterationTime(elapsed_seconds);
        state.counters["Verifications"] = verifications/readNo;
        state.counters["Sensitivity"] = tp/verifications;
    }
}

static void IBFArguments(benchmark::internal::Benchmark* b)
{
    //b->Args({64, 18, 31, 3});
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo <= 64) || binNo==128 || binNo==512 || binNo==1024 || binNo==2048 || binNo==4096 || binNo==8192)
            continue;
        for (int32_t k = 19; k < 20; ++k)
        {
            // 35 = 4GiB, 36 = 8GiB, 37 = 16GiB
            for (int32_t bits = 30; bits <= 31; ++bits )
            {
                for (int32_t hashNo = 3; hashNo < 4; ++hashNo)
                {

                   b->Args({binNo, k, bits, hashNo});
                }
            }
        }
    }
}

//[[maybe_unused]]
/*static void DAArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 13; k <= 15; ++k)
        {
            b->Args({binNo, k});
        }a
    }
}*/

//BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed, Simple)->Apply(IBFArguments);
BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed, A1)->Apply(IBFArguments);
//BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, Simple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedSimple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedArray)->Apply(IBFWhichArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_DA, Dna, CompressedSimple)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_DA, Dna, CompressedArray)->Apply(DAWhichArguments)->UseManualTime();

BENCHMARK_MAIN();
