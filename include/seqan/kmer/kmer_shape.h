#ifndef KMER_SHAPE
#define KMER_SHAPE

namespace seqan {

template<typename TAlphabet>
struct KmerShape<SimpleShape>
{
    typedef Shape<TAlphabet, SimpleShape> TShape;
    Shape<TAlphabet, SimpleShape> x;
    void resize(TShape shape, uint64_t kmerSize)
    {
        resize(shape, kmerSize);
    }
};

template<typename TAlphabet>
struct KmerShape<TAlphabet, ShapeIlieB3>
{
    Shape<TAlphabet, ShapeIlieB3> x;
    void resize(TShape &, uint64_t &)
    {
    }
};

template<typename TAlphabet>
struct KmerShape<TAlphabet, ShapePatternHunter>
{
    Shape<TAlphabet, ShapePatternHunter> x;
    void resize(TShape &, uint64_t &)
    {
    }
};

template<typename TAlphabet>
struct KmerShape<TAlphabet, ShapeIlieA1>
{
    Shape<TAlphabet, ShapeIlieA1> x;
    void resize(TShape &, uint64_t &)
    {
    }
};

#endif // KMER_SHAPE

