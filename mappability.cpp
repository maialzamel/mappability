//#include <inttypes.h>

//uint64_t startPos = 0;

#include "common.h"

template <typename TDistance, typename TIndex, typename TText>
inline void run(TIndex /*const*/ & index, TText const & text, StringSet<CharString> const & ids,
                CharString const & outputPath, unsigned const errors, unsigned const length,
                unsigned const chromosomeId, unsigned const overlap, unsigned const threads)
{
    sdsl::int_vector<16> c(seqan::length(text) - length + 1);
    // TODO: is there an upper bound? are we interested whether a k-mer has 60.000 or 70.000 hits?
    cout << mytime() << "Vector initialized (size: " << c.size() << ")." << endl;
    switch (errors)
    {
        case 0: runAlgo2<0/*, TDistance*/>(index, text, length, c, length - overlap, threads);
                break;
        case 1: runAlgo2<1/*, TDistance*/>(index, text, length, c, length - overlap, threads);
                break;
        case 2: runAlgo2<2/*, TDistance*/>(index, text, length, c, length - overlap, threads);
                break;
        // case 3: run<3, TDistance>(index, text, ids, outputPath, length, chromosomeId);
        //        break;
        // case 4: run<4, TDistance>(index, text, ids, outputPath, length, chromosomeId);
        //        break;
        default: cout << "E = " << errors << " not yet supported.\n";
                 exit(1);
    }
    cout << mytime() << "Done.\n";

    std::string output_path = toCString(outputPath);
    output_path += std::to_string(chromosomeId);
    store_to_file(c, output_path);
    cout << mytime() << "Saved to disk: " << output_path << '\n';

    for (unsigned i = 0; i < 200; ++i)
        std::cout << c[i] << ' ';
    std::cout << '\n';
}

template <typename TChar, typename TAllocConfig, typename TDistance>
inline void run(StringSet<CharString>  & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const singleIndex, unsigned const overlap, unsigned const threads)
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
        run<TDistance>(index, concatText, ids, outputPath, errors, length, 0, overlap, threads);
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
            run<TDistance>(index, text, ids, outputPath, errors, length, i, overlap, threads);
        }
    }
}

template <typename TChar, typename TAllocConfig>
inline void run(StringSet<CharString> & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const indels, bool const singleIndex,
                unsigned const overlap, unsigned const threads)
{
    if (indels) {
        cout << "E = " << errors << " not yet supported.\n";
        exit(1);
        // run<TChar, TAllocConfig, EditDistance>(ids, indexPath, outputPath, errors, length, singleIndex, overlap, threads);
    }
    else
        run<TChar, TAllocConfig, HammingDistance>(ids, indexPath, outputPath, errors, length, singleIndex, overlap, threads);
}

template <typename TChar>
inline void run(StringSet<CharString> & ids, CharString const & indexPath, CharString const & outputPath,
                unsigned const errors, unsigned const length, bool const indels, bool const singleIndex,
                bool const mmap, unsigned const overlap, unsigned const threads)
{
    if (mmap)
        run<TChar, MMap<> >(ids, indexPath, outputPath, errors, length, indels, singleIndex, overlap, threads);
    else
        run<TChar, Alloc<> >(ids, indexPath, outputPath, errors, length, indels, singleIndex, overlap, threads);
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

    addOption(parser, ArgParseOption("o", "overlap", "Length of overlap region (usually: the bigger, the faster)", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "overlap");

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", omp_get_max_threads());

    //addOption(parser, ArgParseOption("x", "startPos", "startPos (for debugging large genomes", ArgParseArgument::INT64, "INT64"));
    //setDefaultValue(parser, "startPos", 0);

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString indexPath, outputPath, _indexPath;
    unsigned errors, length, overlap, threads;
    getOptionValue(errors, parser, "errors");
    getOptionValue(length, parser, "length");
    getOptionValue(overlap, parser, "overlap");
    getOptionValue(indexPath, parser, "index");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(threads, parser, "threads");
    bool indels = isSet(parser, "indels");
    bool mmap = isSet(parser, "mmap");
    //getOptionValue(startPos, parser, "startPos");

    if (overlap > length - errors - 2)
    {
        std::cout << "ERROR: overlap should be <= K - E - 2\n";
    }

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
        run<Dna>(ids, indexPath, outputPath, errors, length, indels, singleIndex, mmap, overlap, threads);
    else
        run<Dna5>(ids, indexPath, outputPath, errors, length, indels, singleIndex, mmap, overlap, threads);
}
