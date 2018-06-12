#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include <chrono>

#include <sdsl/int_vector.hpp>

#include "algo1.hpp"
#include "algo2.hpp"

using namespace std;
using namespace seqan;

// reduce space consumption.
// requires genome not to have more than ~ 4 gigabases
// multi-sequence fasta file must contain less than ~64k sequences of at most ~ 4 gigabases in total

namespace seqan {

template <typename TChar, typename TOwner>
struct SAValue<StringSet<String<TChar>, TOwner> >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TOwner>
struct SAValue<String<TChar, TOwner> >
{
    typedef uint32_t Type;
};

};

std::string mytime()
{
    auto r = time(nullptr);
    auto c = ctime(&r);
    std::string buf(c);
    buf.insert(0, "[");
    buf.append("] ");
    buf.erase(remove(buf.begin(), buf.end(), '\n'), buf.end());
    return buf;
}

// index type

// typedef StringSet<String<Dna>, Owner<ConcatDirect<> > > TDnaText;
// typedef StringSet<String<Dna, MMap<> >, Owner<ConcatDirect<> > > TDnaTextMMap;

// typedef Concatenator<TText>::Type TConcat;
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
typedef BidirectionalIndex<FMIndex<void, TMyFastConfig> > TIndexConfig;
// typedef Index<TText, TIndexConfig> TIndex;
// typedef Index<TTextMMap, TIndexConfig> TIndexMMap;
// typedef Index<DnaString, TIndexConfig> TIndex;
// typedef Iter<TIndex, VSTree<TopDown<> > > TIter;
