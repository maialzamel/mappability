#include "common.h"

#include <sdsl/int_vector.hpp>

template <unsigned errors, typename TDistance, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, StringSet<CharString> const & ids,
                CharString const & outputPath, unsigned const length, unsigned const chromosomeId)
{
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length - 1);

    uint64_t textLength = seqan::length(text);
    // TODO: is there an upper bound? are we interested whether a k-mer has 60.000 or 70.000 hits?
    sdsl::int_vector<16> c(textLength - length + 1);
    cout << mytime() << "Vector initialized (size: " << c.size() << ")." << endl;

    // TODO(cpockrandt): think about scheduling strategy
    #pragma omp parallel for schedule(dynamic, 1000000)
    for (uint64_t i = 0; i < textLength - length + 1; i += 2)
    { // TODO: abort criteria in _optimalSearchScheme
        unsigned hitsL = 0, hitsR = 0;
        // TODO(cpockrandt): for more than 2 errors we need to filter duplicates when counting or
        // choose disjunct search schemes. also we need to do this for EditDistance!
        auto delegate = [&hitsL, &hitsR, i, length, textLength, &text](auto /*const &*/it, auto const & /*read*/, unsigned const errors_spent) {
            if (errors_spent == errors)
            {
                if (i + length < textLength)
                {
                    auto it2 = it;
                    if (goDown(it2, text[i + length], Rev()))
                        hitsR = std::min((uint64_t) countOccurrences(it2) + hitsR, (uint64_t) (1 << 16) - 1);
                }

                if (goDown(it, text[i], Fwd()))
                    hitsL = std::min((uint64_t) countOccurrences(it) + hitsL, (uint64_t) (1 << 16) - 1);
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
                    } while(goRight(it, Fwd()));
                }
            }
        };

        auto const & needle = infix(text, i + 1, i + length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        _optimalSearchScheme(delegate, it, needle, scheme, TDistance());
        c[i] = hitsL;
        if (i + 1 < c.size())
            c[i + 1] = hitsR;
    }
    cout << mytime() << "Done.\n";

    std::string output_path = toCString(outputPath);
    output_path += std::to_string(chromosomeId);
    store_to_file(c, output_path);
    cout << mytime() << "Saved to disk: " << output_path << '\n';
}

template <typename TDistance, typename TIndex, typename TText>
inline void run(TIndex /*const*/ & index, TText const & text, StringSet<CharString> const & ids,
                CharString const & outputPath, unsigned const errors, unsigned const length,
                unsigned const chromosomeId)
{
    switch (errors)
    {
        case 0: run<0, TDistance>(index, text, ids, outputPath, length, chromosomeId);
                break;
        case 1: run<1, TDistance>(index, text, ids, outputPath, length, chromosomeId);
                break;
        case 2: run<2, TDistance>(index, text, ids, outputPath, length, chromosomeId);
                break;
        case 3: run<3, TDistance>(index, text, ids, outputPath, length, chromosomeId);
               break;
        case 4: run<4, TDistance>(index, text, ids, outputPath, length, chromosomeId);
               break;
        default: cout << "E = " << errors << " not yet supported.\n";
                 exit(1);
    }
}

template <typename TChar, typename TAllocConfig, typename TDistance>
inline void run(StringSet<CharString>  & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const singleIndex)
{
    typedef String<TChar, TAllocConfig> TString;
    if (singleIndex)
    {
        Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
        open(index, toCString(indexPath), OPEN_RDONLY);

        // TODO(cpockrandt): replace with a ConcatView
        auto & text = indexText(index);
        typename Concatenator<StringSet<TString, Owner<ConcatDirect<> > >>::Type concatText = concat(text);

        cout << mytime() << "Index loaded." << endl;
        run<TDistance>(index, concatText, ids, outputPath, errors, length, 0);
    }
    else
    {
        for (unsigned i = 0; i < seqan::length(ids); ++i)
        {
            std::string _indexPath = toCString(indexPath);
            _indexPath += "." + to_string(i);
            Index<TString, TIndexConfig> index;
            open(index, toCString(_indexPath), OPEN_RDONLY);
            auto & text = indexText(index);
            cout << mytime() << "Index of " << ids[i] << " loaded." << endl;
            run<TDistance>(index, text, ids, outputPath, errors, length, i);
        }
    }
}

template <typename TChar, typename TAllocConfig>
inline void run(StringSet<CharString> & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const indels, bool const singleIndex)
{
    if (indels)
        run<TChar, TAllocConfig, EditDistance>(ids, indexPath, outputPath, errors, length, singleIndex);
    else
        run<TChar, TAllocConfig, HammingDistance>(ids, indexPath, outputPath, errors, length, singleIndex);
}

template <typename TChar>
inline void run(StringSet<CharString> & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const indels, bool const singleIndex,
                bool const mmap)
{
    if (mmap)
        run<TChar, MMap<> >(ids, indexPath, outputPath, errors, length, indels, singleIndex);
    else
        run<TChar, Alloc<> >(ids, indexPath, outputPath, errors, length, indels, singleIndex);
}

int main(int argc, char *argv[])
{
    // Argument parser
    ArgumentParser parser("Mappabilty");
    addDescription(parser,
        "App for calculating the mappability values. Only supports Dna4/Dna5 so far.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("i", "indels", "Turns on indels (EditDistance). "
        "If not selected, only mismatches will be considered."));

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString indexPath, outputPath, _indexPath;
    unsigned errors, length;
    getOptionValue(errors, parser, "errors");
    getOptionValue(length, parser, "length");
    getOptionValue(indexPath, parser, "index");
    getOptionValue(outputPath, parser, "output");
    bool indels = isSet(parser, "indels");
    bool mmap = isSet(parser, "mmap");

    cout << mytime() << "Program started." << endl;

    bool singleIndex;
    CharString alphabet;
    StringSet<CharString> ids;

    _indexPath = indexPath;
    _indexPath += ".singleIndex";
    open(singleIndex, toCString(_indexPath));

    _indexPath = indexPath;
    _indexPath += ".alphabet";
    open(alphabet, toCString(_indexPath));

    _indexPath = indexPath;
    _indexPath += ".ids";
    open(ids, toCString(_indexPath), OPEN_RDONLY);

    if (alphabet == "dna4")
        run<Dna>(ids, indexPath, outputPath, errors, length, indels, singleIndex, mmap);
    else
        run<Dna5>(ids, indexPath, outputPath, errors, length, indels, singleIndex, mmap);
}
