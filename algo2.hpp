using namespace seqan;

template <unsigned errors, typename TIndex>
inline void runAlgo2Prototype(TIndex & index, auto const & text, unsigned const length, sdsl::int_vector<16> & c, unsigned const threads)
{
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length - 1);

    uint64_t textLength = seqan::length(text);

    #pragma omp parallel for schedule(dynamic, 1000000) num_threads(threads)
    for (uint64_t i = 0; i < textLength - length + 1; i += 2)
    {
        unsigned hitsL = 0, hitsR = 0;
        auto delegate = [&hitsL, &hitsR, i, length, textLength, &text](auto /*const &*/it, auto const & /*read*/, unsigned const errors_spent) {
            // cout << "i = " << i << ": " << representative(it, Fwd()) << " (" << countOccurrences(it) << ")" << endl;
            if (errors_spent == errors)
            {
                if (i + length < textLength)
                {
                    auto it2 = it;
                    if (goDown(it2, text[i + length], Rev()))
                    {
                        hitsR = std::min((uint64_t) countOccurrences(it2) + hitsR, (uint64_t) (1 << 16) - 1);
                    }
                }

                if (goDown(it, text[i], Fwd()))
                {
                    hitsL = std::min((uint64_t) countOccurrences(it) + hitsL, (uint64_t) (1 << 16) - 1);
                }
            }
            else
            {
                if (i + length < textLength)
                {
                    auto it2 = it;
                    if (goDown(it2, Rev()))
                    {
                        do {
                            hitsR = std::min((uint64_t) countOccurrences(it2) + hitsR, (uint64_t) (1 << 16) - 1);
                        } while (goRight(it2, Rev()));
                    }
                }

                if (goDown(it, Fwd()))
                {
                    do {
                        hitsL = std::min((uint64_t) countOccurrences(it) + hitsL, (uint64_t) (1 << 16) - 1);
                    } while (goRight(it, Fwd()));
                }
            }
        };

        auto const & needle = infix(text, i + 1, i + length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
        c[i] = hitsL;
        if (i + 1 < c.size())
            c[i + 1] = hitsR;
    }
}

template <typename TIter>
void extendExact(TIter it, unsigned * hits, auto & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
    if (b - a + 1 == length)
    {
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], (uint64_t) (1 << 16) - 1);
        return;
    }
    //if (b + 1 <= bb)
    //{
        auto it2 = it;
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) / 2); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            bool success = true;
            for (uint64_t i = b + 1; i <= b_new && success; ++i)
            {
                success = goDown(it2, text[i], Rev());
            }
            if (success)
                extendExact(it2, hits, text, length, a, b_new, ab, bb);
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) / 2, 0l);
        bool success = true;
        for (int64_t i = a - 1; i >= a_new && success; --i)
        {
            success = goDown(it, text[i], Fwd());
        }
        if (success)
            extendExact(it, hits, text, length, a_new, b, ab, bb);
    }
}

// forward
template <typename TIter>
void extend(TIter it, unsigned * hits, unsigned errors_left, auto & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb // entire interval
);

// TODO: remove text everywhere: auto & text = indexText(index(it));
template <typename TIter>
void approxSearch(TIter it, unsigned * hits, unsigned errors_left, auto & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb, // entire interval
            uint64_t b_new,
            Rev const & /*tag*/
)
{
    if (b == b_new)
    {
        extend(it, hits, errors_left, text, length, a, b, ab, bb);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1]);
                approxSearch(it, hits, errors_left - delta, text, length, a, b + 1, ab, bb, b_new, Rev());
            } while (goRight(it, Rev()));
        }
    }
    else
    {
        for (uint64_t i = b + 1; i <= b_new; ++i)
        {
            if (!goDown(it, text[i], Rev()))
                return;
        }
        extendExact(it, hits, text, length, a, b_new, ab, bb);
        //approxSearch(it, hits, errors_left, text, length, a, b_new, ab, bb, Rev());
    }
}
template <typename TIter>
void approxSearch(TIter it, unsigned * hits, unsigned errors_left, auto & text, unsigned const length,
                  uint64_t a, uint64_t b, // searched interval
                  uint64_t ab, uint64_t bb, // entire interval
                  int64_t a_new,
                  Fwd const & /*tag*/
)
{
    if (a == a_new)
    {
        extend(it, hits, errors_left, text, length, a, b, ab, bb);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1]);
                approxSearch(it, hits, errors_left - delta, text, length, a - 1, b, ab, bb, a_new, Fwd());
            } while (goRight(it, Fwd()));
        }
    }
    else
    {
        for (int64_t i = a - 1; i >= a_new; --i)
        {
            if (!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact(it, hits, text, length, a_new, b, ab, bb);
        //approxSearch(it, hits, errors_left, text, length, a_new, b, ab, bb, Fwd());
    }
}

template <typename TIter>
void extend(TIter it, unsigned * hits, unsigned errors_left, auto & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
    if (errors_left == 0)
    {
        extendExact(it, hits, text, length, a, b, ab, bb);
        return;
    }
    if (b - a + 1 == length)
    {
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], (uint64_t) (1 << 16) - 1);
        return;
    }
    //if (b + 1 <= bb)
    //{
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) / 2); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            approxSearch(it, hits, errors_left, text, length,
                         a, b, // searched interval
                         ab, bb, // entire interval
                         b_new,
                         Rev()
            );
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) / 2, 0l);
        approxSearch(it, hits, errors_left, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd()
        );
    }
}

template <unsigned errors, typename TIndex>
inline void runAlgo2(TIndex & index, auto const & text, unsigned const length, sdsl::int_vector<16> & c, unsigned const overlap, unsigned const threads)
{
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, overlap);

    uint64_t textLength = seqan::length(text);

    const uint64_t max_i = textLength - length + 1;
    const uint64_t step_size = length - overlap + 1;
    #pragma omp parallel for schedule(guided) num_threads(threads)
    //#pragma omp parallel for schedule(dynamic, 1000000)
    // TODO: if we choose a multiple of step_size * |cyclic_rotations of int_vector| as chunks (not just chunksize), we dont need to worry about locking
    for (uint64_t i = 0/*startPos*/; i < max_i; i += step_size)
    {
        unsigned hits[length - overlap + 1] = {};
        auto delegate = [&hits, i, length, textLength, overlap, &text](auto /*const &*/it, auto const & /*read*/, unsigned const errors_spent) {
            uint64_t bb = std::min(textLength - 1, i + length - 1 + length - overlap);
            //cout << "i = " << i << ": " << representative(it, Fwd()) << " (" << countOccurrences(it) << ")" << endl;
            extend(it, hits, errors - errors_spent, text, length,
                i + length - overlap, i + length - 1, // searched interval
                i, bb // entire interval
            );
        };

        auto const & needle = infix(text, i + length - overlap, i + length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
        uint64_t max_pos = std::min(i + length - overlap, textLength - length/* + 1*/);
        for (uint64_t j = i; j <= max_pos; ++j)
            c[j] = hits[j - i];
    }
}
