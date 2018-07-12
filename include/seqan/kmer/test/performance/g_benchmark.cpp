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

//CharString baseDir{"/srv/public/enricoseiler/share/small_benchmark/"};
CharString baseDir{"/storage/mi/kingpfogel/small_benchmark/"};
//uint64_t e{2};
HammingDistance hd;
ThreshHeuristic th;

//11101011001100101111
typedef GappedShape<
    HardwiredShape< 1, 1, 2, 2, 1, 3, 1, 3, 2, 1, 1, 1 >
> ShapePaper1;

//11101001101001110111
typedef GappedShape<
    HardwiredShape< 1, 1, 2, 3, 1, 2, 3, 1, 1, 2, 1, 1 >
> ShapePaper2;

//11110101011001101111
typedef GappedShape<
    HardwiredShape< 1, 1, 1, 2, 2, 2, 1, 3, 1, 2, 1, 1, 1 >
> ShapePaper3;

//1110110110010101111
typedef GappedShape<
    HardwiredShape< 1, 1, 2, 1, 2, 1, 3, 2, 2, 1, 1, 1 >
> ShapePaper4;
//111010110100110111 weight 12
typedef GappedShape<
    HardwiredShape< 1,1,2,2,1,2,3,1,2,1,1 >
> ShapePaper5;

//110111010110100111
typedef GappedShape<
    HardwiredShape< 1,2,1,1,2,2,1,2,3,1,1 >
> SpEEDw12_s18;
//1110100100010010001010111
typedef GappedShape<
    HardwiredShape< 1,1,2,3,4,3,4,2,2,1,1 >
> SpEEDw12_s25;

//1101101001110101111 
typedef GappedShape<
    HardwiredShape< 1,2,1,2,3,1,1,2,2,1,1,1 >
> SpEEDw13_s19;

//11110101000100010010011011 
typedef GappedShape<
    HardwiredShape< 1,1,1,2,2,4,4,3,3,1,2,1 >
> SpEEDw13_s26;

//111110101110011010110111 
typedef GappedShape<
    HardwiredShape< 1,1,1,1,2,2,1,1,3,1,2,2,1,2,1,1 >
> SpEEDw17_s24;

//11101101001100001010100001111011 
typedef GappedShape<
    HardwiredShape< 1,1,2,1,2,3,1,5,2,2,5,1,1,1,2,1 >
> SpEEDw17_s32;

//1110110110101011100111111 
typedef GappedShape<
    HardwiredShape< 1,1,2,1,2,1,2,2,2,1,1,3,1,1,1,1,1 >
> SpEEDw18_s25;

//111010110011000011010000101110111 
typedef GappedShape<
    HardwiredShape< 1,1,2,2,1,3,1,5,1,2,5,2,1,1,2,1,1>
> SpEEDw18_s33;

//11110110110101110111001111 
typedef GappedShape<
    HardwiredShape< 1,1,1,2,1,2,1,2,1,1,2,1,1,3,1,1,1>
> SpEEDw19_s26;

//1111010100110010000011010110110111 
typedef GappedShape<
    HardwiredShape< 1,1,1,2,2,3,1,3,6,1,2,2,1,2,1,2,1,1 >
> SpEEDw19_s34;


template <typename TAlphabet, typename TFilter, typename TShapeSpec, unsigned offset>
static void insertKmer_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);

    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter, TShapeSpec, offset> ibf(bins, hash, k, (1ULL << bits));

    for (auto _ : state)
    {
        insertKmerDir(ibf, toCString(baseDir), 8);
    }
    CharString shape;
    if constexpr (std::is_same_v<TShapeSpec, SimpleShape>) {
        shape = CharString("_SimpleShape");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper1>) {
        shape = CharString("ShapePaper1");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper2>) {
        shape = CharString("ShapePaper2");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper3>) {
        shape = CharString("ShapePaper3");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper4>) {
        shape = CharString("ShapePaper4");
    }

    CharString storage("/srv/public/neesemann/ibf/");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    if constexpr(std::is_same_v<TShapeSpec, SimpleShape>){
        append(storage, CharString(std::to_string(k)));
        append(storage, CharString("_"));
    }
    append(storage, CharString(std::to_string(bits)));
    append(storage, shape);
    append(storage, CharString("_Uncompressed_ibf.filter"));
    store(ibf, storage);
    
    state.counters["Size"] = ibf.filterVector.size_in_mega_bytes();
}


template <typename TAlphabet, typename TFilter, typename TShapeSpec, unsigned offset>
static void select_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    auto e = state.range(4);

    int16_t t;
    CharString shape;
    if constexpr (std::is_same_v<TShapeSpec, SimpleShape>) {
        shape = CharString("_SimpleShape");
        t = std::ceil((100-k+1)/offset)-std::ceil(k/offset)*e;
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper1>) {
        shape = CharString("ShapePaper1");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper2>) {
        shape = CharString("ShapePaper2");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;        
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper3>) {
        shape = CharString("ShapePaper3");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper4>) {
        shape = CharString("ShapePaper4");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else {
        t=0;
        shape=CharString("");
    }
    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter, TShapeSpec, offset> ibf(bins, hash, k, (1ULL<<bits));

    CharString storage("/srv/public/neesemann/ibf/");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    if constexpr(std::is_same_v<TShapeSpec, SimpleShape>){
        append(storage, CharString(std::to_string(k)));
        append(storage, CharString("_"));
    }
    append(storage, CharString(std::to_string(bits)));
    append(storage, shape);
    append(storage, CharString("_Uncompressed_ibf.filter"));
   // std::cerr << storage << std::endl;
    retrieve(ibf, storage);
if(t>0){
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
            append(file, CharString{"/reads_e"});
            append(file, CharString{std::to_string(e)});
            append(file, CharString{"/bin_"});
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
        state.counters["10_Threshold"] = t;
        state.counters["4_TP"] = tp;
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
}


template <typename TAlphabet, typename TFilter, typename TShapeSpec, unsigned offset>
static void insertKmer_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing, TFilter, TShapeSpec, offset> da (bins, k);

    for (auto _ : state)
    {
        insertKmerDir(da, toCString(baseDir), 8);

    }

    CharString storage("/srv/public/neesemann/direct_filter/");
    CharString shape;
    if constexpr (std::is_same_v<TShapeSpec, SimpleShape>) {
        shape = CharString("SimpleShape");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper1>) {
        shape = CharString("ShapePaper1");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper2>) {
        shape = CharString("ShapePaper2");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper3>) {
        shape = CharString("ShapePaper3");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper4>) {
        shape = CharString("ShapePaper4");
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper5>) {
        shape = CharString("ShapePaper5");
    }
    else {
        shape=CharString("xxx");
    }
    append(storage, shape);
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString("da.filter"));
    store(da, storage);

    state.counters["Size"] = da.filterVector.size_in_mega_bytes();
}

template <typename TAlphabet, typename TFilter, typename TShapeSpec, unsigned offset>
static void select_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto e = state.range(2);
//    int16_t t = ((100-k+offset - ((100-k) % offset))/offset) - (1+std::ceil((k-offset)/offset))*2;
    KmerFilter<TAlphabet, DirectAddressing, TFilter, TShapeSpec, offset> da (bins, k);

    //uint64_t readNo{1};

    CharString storage("/srv/public/neesemann/direct_filter/");
    int16_t t;
    CharString shape;
    if constexpr (std::is_same_v<TShapeSpec, SimpleShape>) {
        shape = CharString("SimpleShape");
        t = std::ceil((100-k+1)/offset)-std::ceil(k/offset)*e;
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper1>) {
        shape = CharString("ShapePaper1");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper2>) {
        shape = CharString("ShapePaper2");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;        
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper3>) {
        shape = CharString("ShapePaper3");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper4>) {
        shape = CharString("ShapePaper4");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else if constexpr(std::is_same_v<TShapeSpec, ShapePaper5>) {
        shape = CharString("ShapePaper5");
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
    }
    else {
        typedef Shape<Dna, TShapeSpec> TShapeSpec2;
        TShapeSpec2 sia2;
        t = qgramThreshold(sia2, 100, e, hd, th);
        shape=CharString("xxx");
    }
    append(storage, shape);
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString("da.filter"));

    auto fullTime = std::chrono::high_resolution_clock::now();

    double loadingTime{0.0};
    auto start = std::chrono::high_resolution_clock::now();
    retrieve(da, storage);
    auto end   = std::chrono::high_resolution_clock::now();
    loadingTime = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();


    //double compressionTime{0.0};
    start = std::chrono::high_resolution_clock::now();
//    da.filterVector.compress(0); // Loading automatically compresses
    end   = std::chrono::high_resolution_clock::now();
//    compressionTime = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();

    std::atomic<uint64_t> verifications{0};
    std::atomic<uint64_t> tp{0};
    std::atomic<uint64_t> p{0};
    std::atomic<uint64_t> fp{0};
    std::atomic<uint64_t> fn{0};
    std::atomic<uint64_t> readNo{0};
if(t>0){
    for (auto _ : state)
    {

        double selectTime{0.0};
        double ioTime{0.0};
        Semaphore thread_limiter(8);
        std::mutex mtx;
        std::mutex mtx2;
        std::vector<std::future<void>> tasks;

        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads_e"});
            append(file, CharString{std::to_string(e)});
            append(file, CharString{"/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fastq"));

            tasks.emplace_back(
                std::async(std::launch::async, [&, file, i] {
                    Critical_section _(thread_limiter);
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
                        auto start = std::chrono::high_resolution_clock::now();
                        readRecord(id, seq, seqFileIn);
                        auto end   = std::chrono::high_resolution_clock::now();
                        mtx2.lock();
                        ioTime += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
                        mtx2.unlock();

                        start = std::chrono::high_resolution_clock::now();
                        auto res = select(da, seq, t);
                        end   = std::chrono::high_resolution_clock::now();
                        ++readNo;
                        mtx.lock();
                        selectTime += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
                        mtx.unlock();

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
                })
            );
        }

        for (auto &&task : tasks){
            task.get();
        }

        auto fullTime2   = std::chrono::high_resolution_clock::now();
        state.counters["10_Threshold"] = t;
        state.counters["5_TP"] = tp.load();
        state.counters["6_FN"] = fn.load();
        state.counters["7_FP"] = fp.load();
        state.counters["8_P"] = p.load();
        state.counters["99_readNo"] = readNo.load();
        state.counters["9_verifications"] = verifications.load();
        state.counters["0_Verifications"] = static_cast<double>(verifications.load())/readNo.load();
        state.counters["1_Sensitivity"] = static_cast<double>(tp.load())/readNo.load();
        state.counters["2_Precision"] = static_cast<double>(tp.load())/p.load();
        state.counters["3_FNR"] = static_cast<double>(fn.load())/readNo.load();
        state.counters["4_FDR"] = static_cast<double>(fp.load())/p.load();
  //      state.counters["compressionTime"] = compressionTime;
        state.counters["loadingTime"] = loadingTime;
        state.counters["ioTime"] = ioTime;
        state.counters["selectTime"] = selectTime;
        state.counters["vectorSize"] = da.filterVector.size_in_mega_bytes();
        state.counters["fullTime"] = std::chrono::duration_cast<std::chrono::duration<double> >(fullTime2 - fullTime).count();
        break;
    }
}//ifend
}

[[maybe_unused]]
static void IBFSelectArguments(benchmark::internal::Benchmark* b)
{
    //b->Args({64, 19, 30, 3, 10});
    //b->Args({64, 19, 30, 3, 20});
    //b->Args({64, 19, 30, 3, 30});
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 17; k < 20; ++k)
        {
            // 35 = 4GiB, 36 = 8GiB, 37 = 16GiB
            for (int32_t bits = 26; bits <= 32; ++bits )
            {
                for (int32_t hashNo = 3; hashNo < 4; ++hashNo)
                {
                   for (int16_t error: {2,3,4,5,10})
                   {
                        b->Args({binNo, k, bits, hashNo, error}); 
                    }
                }
            }
        }
    }  
}

[[maybe_unused]]
static void IBFInsertArguments(benchmark::internal::Benchmark* b)
{
//    b->Args({64, 20, 30, 3});
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 17; k < 20; ++k)
        {
            // 35 = 4GiB, 36 = 8GiB, 37 = 16GiB
            for (int32_t bits = 26; bits <= 32; ++bits )
            {
                for (int32_t hashNo = 3; hashNo < 4; ++hashNo)
                {
                   b->Args({binNo, k, bits, hashNo});
                }
            }
        }
    }  
}

[[maybe_unused]]
static void DASelectArguments(benchmark::internal::Benchmark* b)
{

    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 12; k <= 14; ++k)
        {
            for (int16_t error: {2,3,4,5,10})
            {   
                b->Args({binNo, k, error});
            }
        }
    }
}

[[maybe_unused]]
static void DAInsertArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 12; k <= 14; ++k)
        {
            b->Args({binNo, k});
        }
    }
}


//BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed, ShapeIlieA1, 1u)->Apply(IBFArguments);
//BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed, SimpleShape)->Apply(IBFArguments);



//SimpleShapes#insertKmer
//BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed, SimpleShape, 1u)->Apply(IBFInsertArguments);
BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed, SimpleShape, 1u)->Apply(DAInsertArguments);
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, SimpleShape, 1u)->Apply(DASelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, SimpleShape, 12u)->Apply(DASelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, SimpleShape, 14u)->Apply(DASelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, SimpleShape, 13u)->Apply(DASelectArguments)->UseManualTime();
//SimpleShapes#select
//BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, SimpleShape, 1u)->Apply(IBFSelectArguments)->UseManualTime();
/*BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, SimpleShape, 19u)->Apply(IBFSelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, SimpleShape, 17u)->Apply(IBFSelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, SimpleShape, 18u)->Apply(IBFSelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, SimpleShape, 2u)->Apply(IBFSelectArguments)->UseManualTime();
*/

//GappedShapes 

/*BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed, ShapePaper1, 1u)->Apply(IBFInsertArguments);
BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed, ShapePaper2, 1u)->Apply(IBFInsertArguments);
BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed, ShapePaper3, 1u)->Apply(IBFInsertArguments);
BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed, ShapePaper4, 1u)->Apply(IBFInsertArguments);
*/
/*BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, ShapePaper1, 1u)->Apply(IBFSelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, ShapePaper2, 1u)->Apply(IBFSelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, ShapePaper3, 1u)->Apply(IBFSelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, ShapePaper4, 1u)->Apply(IBFSelectArguments)->UseManualTime();
*/
/*BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed, SpEEDw12_s18, 1u)->Apply(DAInsertArguments);
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, SpEEDw12_s18, 1u)->Apply(DASelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed, SpEEDw12_s25, 1u)->Apply(DAInsertArguments);
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, SpEEDw12_s25, 1u)->Apply(DASelectArguments)->UseManualTime();
BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed, SpEEDw13_s19, 1u)->Apply(DAInsertArguments);
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, SpEEDw13_s19, 1u)->Apply(DASelectArguments)->UseManualTime();*/
BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed, SpEEDw13_s26, 1u)->Apply(DAInsertArguments);
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, SpEEDw13_s26, 1u)->Apply(DASelectArguments)->UseManualTime();

//BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed, ShapePaper5, 1u)->Apply(DAInsertArguments);
//BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed, ShapePaper5, 1u)->Apply(DASelectArguments)->UseManualTime();
//BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, ShapePaper1, 1u)->Args({64, 20, 30, 3, 10})->Args({64, 20, 30, 2, 20})->Args({64, 20, 30, 2, 30})->UseManualTime();
//BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed, ShapePaper2, 1u)->Args({64, 20, 30, 3, 10})->Args({64, 20, 30, 2, 20})->Args({64, 20, 30, 2, 30})->UseManualTime();


BENCHMARK_MAIN();

