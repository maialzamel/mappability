// reduce space consumption.
// requires genome not to have more than ~ 4 gigabases
// multi-sequence fasta file must contain less than ~64k sequences of at most ~ 4 gigabases in total

namespace seqan {

template <typename TChar, typename TAlloc, typename TOwner>
struct SAValue<StringSet<String<TChar, TAlloc>, TOwner> >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TAlloc>
struct SAValue<String<TChar, TAlloc> >
{
    typedef uint32_t Type;
};

};

struct Options
{
    unsigned errors;
    unsigned length;
    unsigned overlap;
    unsigned threads;
    bool mmap;
    bool indels;
    bool singleIndex;
    seqan::CharString indexPath;
    seqan::CharString outputPath;
    seqan::CharString alphabet;
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

using TMyFastConfig = seqan::FastFMIndexConfig<void, uint32_t, 2, 1>;
using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<void, TMyFastConfig> >;

template <typename TText>
using TIndex = seqan::Index<TText, TIndexConfig>;

// template <typename T>
// struct is_int_vector
// {
//     static constexpr bool value = false;
// };
//
// template <uint8_t width_t>
// struct is_int_vector<sdsl::int_vector<width_t> >
// {
//     static constexpr bool value = true;
// };
//
// template <typename T>
// constexpr bool is_int_vector<T>::value;
//
// template <uint8_t width_t>
// constexpr bool is_int_vector<sdsl::int_vector<width_t> >::value;
