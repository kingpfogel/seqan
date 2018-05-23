#ifndef INCLUDE_SEQAN_KMER_KMER_SHAPE_H_
#define INCLUDE_SEQAN_KMER_KMER_SHAPE_H_

#include <seqan/index.h>
namespace seqan {

    template<typename TAlphabet>
    struct KmerShape<TAlphabet, SimpleShape>
    {
        KmerShape() {}
        ~KmerShape() = default;
        typedef Shape<TAlphabet, SimpleShape> TShape;
        TShape shape;
        void resizeShape(uint64_t kmerSize)
        {
            resize(shape, kmerSize);
        }
    };

    template<typename TAlphabet>
    struct KmerShape<TAlphabet, ShapeIlieB3>
    {
        KmerShape() {}
        ~KmerShape() = default;
        typedef Shape<TAlphabet, ShapeIlieB3> TShape;
        TShape shape;
        void resizeShape(uint64_t kmerSize)
        {
        }
    };


    template<typename TAlphabet>
    struct KmerShape<TAlphabet, ShapePatternHunter>
    {
        KmerShape() {}
        ~KmerShape() = default;
        typedef Shape<TAlphabet, ShapePatternHunter> TShape;
        TShape shape;
        void resizeShape(uint64_t &)
        {
        }
    };

    template<typename TAlphabet>
    struct KmerShape<TAlphabet, ShapeIlieA1>
    {
        KmerShape() {}
        ~KmerShape() = default;
        typedef Shape<TAlphabet, ShapeIlieA1> TShape;
        TShape shape;
        void resizeShape(uint64_t &)
        {
        }
    };
    template<typename TAlphabet>
    struct KmerShape<TAlphabet, ShapeIlieA2>
    {
        KmerShape() {}
        ~KmerShape() = default;
        typedef Shape<TAlphabet, ShapeIlieA2> TShape;
        TShape shape;
        void resizeShape(uint64_t &)
        {
        }
    };
    template<typename TAlphabet>
    struct KmerShape<TAlphabet, ShapeIlieA3>
    {
        KmerShape() {}
        ~KmerShape() = default;
        typedef Shape<TAlphabet, ShapeIlieA3> TShape;
        TShape shape;
        void resizeShape(uint64_t &)
        {
        }
    };
}
#endif  // INCLUDE_SEQAN_KMER_KMER_SHAPE_H_

