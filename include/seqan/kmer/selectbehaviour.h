//default
struct OverlappingKmers_;
typedef Tag<OverlappingKmers_> OverlappingKmers;

//specification: take only every k-th kmer
struct NonOverlappingKmers_;
typedef Tag<NonOverlappingKmers_> NonOverlappingKmers;

std::vector<uint64_t> selectHelper(OverlappingKmers const &, TString const & text)
{
    uint16_t possible = length(text) - kmerSize + 1; // Supports text lengths up to 65535 + k
    std::vector<uint64_t> kmerHashes(possible, 0);

    //TShape kmerShape;
    shape.resizeShape(kmerSize);
    //resize(kmerShape, kmerSize);
    hashInit(shape.shape, begin(text));
    auto it = begin(text);
    for (uint16_t i = 0; i < possible; ++i)
    {
        kmerHashes[i] = hashNext(shape.shape, it);
        ++it;
    }
    return kmerHashes;
}
std::vector<uint64_t> selectHelper(NonOverlappingKmers const &, TString const & text)
{
    uint16_t possible = length(text) - kmerSize + 1; // Supports text lengths up to 65535 + k

    uint16_t x = length(text) % kmerSize;
    uint16_t noOfKmerHashes = ((length(text) - x) / kmerSize)+1;
    if (x == 0)
    {
        noOfKmerHashes -= 1;
    }
    std::vector<uint64_t> kmerHashes(noOfKmerHashes, 0);
    shape.resizeShape(kmerSize);
    hashInit(shape.shape, begin(text));
    auto it = begin(text);
    uint32_t c = 0;
    for (uint32_t i = 0; i < possible; ++i)
    {
        uint64_t kmerHash = hashNext(shape.shape, it);
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

int main()
{
    myFunction(TagA());  // (3)
    myFunction(TagB());  // (4)
    return 0;
}
