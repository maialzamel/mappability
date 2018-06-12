#include "common.h"

template <typename TIV, typename TRng>
void randText(TIV & iv, TRng & rng, unsigned const len)
{
    resize(iv, len);
    for (uint64_t i = 0; i < len; ++i)
        iv[i] = Dna(rng() % 4); // 0 is reserved for sentinel character!
}

int main(int argc, char *argv[])
{
    auto now = chrono::system_clock::now();
    auto seed = chrono::duration_cast<chrono::nanoseconds>(now.time_since_epoch()).count();
    cout << "Seed: " << seed << '\n';
    mt19937_64 rng(seed);

    constexpr unsigned errors = 1;

    while (true)
    {
        unsigned textLength = 100000 + (rng() % 13);

        DnaString genome; // = "TGGCTCTTTG";
        randText(genome, rng, textLength);
        Index<DnaString, TIndexConfig> index(genome);
        indexCreate(index, FibreSALF());
        auto & text = indexText(index);

        for (unsigned j = 0; j < 100; ++j)
        {
            unsigned length = (rng() % 25) + 4;

            sdsl::int_vector<16> c1(seqan::length(text) - length + 1);
            runAlgo2Prototype<errors>(index, genome, length, c1);

            for (unsigned overlap = 0; overlap <= length - errors - 2; ++overlap) // because there have to be enough characters for the infix using search schemes
            {
                sdsl::int_vector<16> c2(seqan::length(text) - length + 1);
                runAlgo2<errors>(index, genome, length, c2, length - overlap);

                for (unsigned i = 0; i < c1.size(); ++i)
                {
                    if (c1[i] != c2[i])
                    {
                        cout << "\nOverlap: " << overlap << ", Length: " << length << '\n';
                        for (unsigned j = 0; j < seqan::length(genome); ++j)
                            cout << genome[j] << ' ';
                        cout << endl;
                        for (unsigned j = 0; j < c1.size(); ++j)
                            cout << c1[j] << ' ';
                        cout << endl;
                        for (unsigned j = 0; j < c2.size(); ++j)
                            cout << c2[j] << ' ';
                        cout << endl;
                        exit(1);
                    }
                }
                cout << "." << flush;
            }
        }
    }
}
