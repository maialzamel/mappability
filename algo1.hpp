using namespace seqan;

template <unsigned errors, typename TIndex, typename TContainer>
inline void runAlgoTrivial(TIndex & index, auto const & text, unsigned const length, TContainer & c, unsigned const /*overlap*/, unsigned const threads)
{
    constexpr uint64_t max_val = (1 << 8) - 1;
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length);

    uint64_t textLength = seqan::length(text);

    #pragma omp parallel for schedule(dynamic, 1000000) num_threads(threads)
    for (uint64_t i = 0; i < textLength - length + 1; ++i)
    {
        unsigned hits = 0;
        auto delegate = [&hits](auto const &it, auto const & /*read*/, unsigned const /*errors*/) {
            if (hits + countOccurrences(it) <= max_val)
                hits += countOccurrences(it);
            else
                hits = max_val;
        };

        auto const & needle = infix(text, i, i + length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
        c[i] = hits;
    }
}
