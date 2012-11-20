// ==========================================================================
//                                mini_bowtie
// ==========================================================================

#include <iostream>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/index_fm.h>
#include <seqan/store.h>

using namespace seqan;

struct ForwardTag {};
struct ReverseTag {};

template <typename TStore, typename TIter, typename TPatternIt>
void addMatchToStore(TStore & fragStore, TPatternIt const & patternIt, TIter const & localIt, ForwardTag)
{
    typedef FragmentStore<>::TAlignedReadStore TAlignedReadStore;
    typedef Value<TAlignedReadStore>::Type TAlignedRead;

    for (unsigned num = 0; num < countOccurrences(localIt); ++num)
    {
        unsigned pos = getOccurrences(localIt)[num].i2;
        TAlignedRead match(length(fragStore.alignedReadStore), position(patternIt), getOccurrences(localIt)[num].i1 ,
            pos,  pos + length(value(patternIt)));
        appendValue(fragStore.alignedReadStore, match);
    }
}

template <typename TStore, typename TIter, typename TPatternIt>
void addMatchToStore(TStore & fragStore, TPatternIt const & patternIt, TIter const & localIt, ReverseTag)
{
    typedef FragmentStore<>::TAlignedReadStore TAlignedReadStore;
    typedef Value<TAlignedReadStore>::Type TAlignedRead;

    for (unsigned num = 0; num < countOccurrences(localIt); ++num)
    {
        unsigned contigLength = length(fragStore.contigStore[getOccurrences(localIt)[num].i1].seq);
        unsigned pos = contigLength - getOccurrences(localIt)[num].i2 - length(value(patternIt));
        TAlignedRead match(length(fragStore.alignedReadStore), position(patternIt), getOccurrences(localIt)[num].i1,
            pos, pos + length(value(patternIt)));
        appendValue(fragStore.alignedReadStore, match);
    }
}

template <typename TIter, typename TStringSet, typename TStore, typename DirectionTag>
void search(TIter & it, TStringSet const & pattern, TStore & fragStore, DirectionTag /*tag*/)
{
    typedef typename Iterator<TStringSet const, Standard>::Type TPatternIter;

    for (TPatternIter patternIt = begin(pattern, Standard()); patternIt != end(pattern, Standard()); ++patternIt)
    {
        // exact search on pattern half
        unsigned startApproxSearch = length(value(patternIt)) / 2;
        if (goDown(it, infix(value(patternIt), startApproxSearch + 1, length(value(patternIt)))))
        {
            for (unsigned i = startApproxSearch; ; --i)
            {
                Dna character = getValue(patternIt)[i];
                for (Dna5 c = MinValue<Dna>::VALUE; c < +ValueSize<Dna>::VALUE; ++c)
                {
                    if (c != character)
                    {
                        TIter localIt = it;
                        if (goDown(localIt, c)){
                            if (goDown(localIt, infix(value(patternIt), 0, i)))
                            {
                                addMatchToStore(fragStore, patternIt, localIt, DirectionTag());
                            }
                        }
                    }
                }
                if (!goDown(it, character))
                    break;
                else if (i == 0)
                {
                    if(IsSameType<DirectionTag, ForwardTag>::VALUE)
                        addMatchToStore(fragStore, patternIt, it, DirectionTag());
                    break;
                }
            }
        }
        goRoot(it);
    }
}

int main(int argc, char *argv[]) 
{
    typedef String<Dna5> TString;
    typedef StringSet<String<Dna5> > TStringSet;
    typedef Index<StringSet<TString>, FMIndex<> > TIndex;
    typedef Iterator<TIndex, TopDown<> >::Type TIter;

    // 0) Handle command line arguments.
    if (argc < 3) {
        std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE: minimapper GENOME.fasta READS.fasta OUT.sam" << std::endl;
        return 1;
    }
    std::cerr << "Phase -2" << std::endl;
    // 1) Load contigs and reads.
    FragmentStore<> fragStore;
    if (!loadContigs(fragStore, argv[1])) return 1;
    if (!loadReads(fragStore, argv[2])) return 1;
    
    std::cerr << "Phase -1" << std::endl;

    StringSet<TString> text;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
    {
        if(length(fragStore.contigStore[i].seq) != 0)
            appendValue(text, fragStore.contigStore[i].seq);
    }
       
    std::cerr << "Phase 0" << std::endl;

    TIndex fmIndex(text);
    std::cerr << "Phase 1" << std::endl;
    TIter it(fmIndex);
    std::cerr << "Phase 2" << std::endl;
    search(it, fragStore.readSeqStore, fragStore, ForwardTag());
    std::cerr << "Phase 3" << std::endl;
    clear(fmIndex);
    clear(it);

    std::cerr << "Phase 4" << std::endl;
    reverse(text);
    reverse(fragStore.readSeqStore);

    std::cerr << "Phase 5" << std::endl;
    fmIndex = TIndex(text);
    std::cerr << "Phase 6" << std::endl;
    it = TIter(fmIndex);
    std::cerr << "Phase 7" << std::endl;
    search(it, fragStore.readSeqStore, fragStore, ReverseTag());
    std::cerr << "Phase 8" << std::endl;
    clear(fmIndex);
    clear(it);

    std::cerr << "Phase 9" << std::endl;
    reverse(text);
    reverse(fragStore.readSeqStore);
    std::ofstream samFile(argv[3], std::ios_base::out);
    write(samFile, fragStore, Sam());

    return 0;
}
