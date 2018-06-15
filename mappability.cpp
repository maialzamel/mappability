//#include <inttypes.h>

//uint64_t startPos = 0;

#include "common.h"

struct Options
{
    unsigned errors;
    unsigned length;
    unsigned overlap;
    unsigned threads;
    bool mmap;
    bool indels;
    bool singleIndex;
    CharString indexPath;
    CharString outputPath;
    CharString alphabet;
};

template <typename TDistance, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, Options & opt, signed chromosomeId)
{
    // sdsl::int_vector<16> c(seqan::length(text) - opt.length + 1);
    std::vector<uint8_t> c(seqan::length(text) - opt.length + 1, 0);

    // TODO: is there an upper bound? are we interested whether a k-mer has 60.000 or 70.000 hits?
    cout << mytime() << "Vector initialized (size: " << c.size() << ")." << endl;
    switch (opt.errors)
    {
        case 0: runAlgo2<0/*, TDistance*/>(index, text, opt.length, c, opt.length - opt.overlap, opt.threads);
                break;
        case 1: runAlgo2<1/*, TDistance*/>(index, text, opt.length, c, opt.length - opt. overlap, opt.threads);
                break;
        case 2: runAlgo2<2/*, TDistance*/>(index, text, opt.length, c, opt.length - opt.overlap, opt.threads);
                break;
        // case 3: run<3, TDistance>(index, text);
        //        break;
        // case 4: run<4, TDistance>(index, text);
        //        break;
        default: cerr << "E = " << opt.errors << " not yet supported.\n";
                 exit(1);
    }
    cout << mytime() << "Done.\n";

    std::string output_path = toCString(opt.outputPath);
    output_path += "_" + to_string(opt.errors) + "_" + to_string(opt.length) + "_" + to_string(opt.overlap);
    if (chromosomeId >= 0)
        output_path += "-" + to_string(chromosomeId);

    // store_to_file(c, output_path);

    std::ofstream outfile(output_path, std::ios::out | std::ofstream::binary);
    std::copy(c.begin(), c.end(), std::ostream_iterator<uint8_t>(outfile));
    outfile.close();

    cout << mytime() << "Saved to disk: " << output_path << '\n';
}

template <typename TChar, typename TAllocConfig, typename TDistance>
inline void run(Options & opt)
{
    typedef String<TChar, TAllocConfig> TString;
    if (opt.singleIndex)
    {
        Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index;
        open(index, toCString(opt.indexPath), OPEN_RDONLY);

        // TODO(cpockrandt): replace with a ConcatView
        auto & text = indexText(index);
        typename Concatenator<StringSet<TString, Owner<ConcatDirect<> > >>::Type concatText = concat(text);

        cout << mytime() << "Index loaded." << endl;
        run<TDistance>(index, concatText, opt, -1 /*no chromosomeId*/);
    }
    else
    {
        StringSet<CharString> ids;
        CharString _indexPath = opt.indexPath;
        _indexPath += ".ids";
        open(ids, toCString(_indexPath), OPEN_RDONLY);

        for (unsigned i = 0; i < seqan::length(ids); ++i)
        {
            std::string _indexPath = toCString(opt.indexPath);
            _indexPath += "." + to_string(i);
            Index<TString, TIndexConfig> index;
            open(index, toCString(_indexPath), OPEN_RDONLY);
            auto & text = indexText(index);
            cout << mytime() << "Index of " << ids[i] << " loaded." << endl;
            run<TDistance>(index, text, opt, i);
        }
    }
}

template <typename TChar, typename TAllocConfig>
inline void run(Options & opt)
{
    if (opt.indels) {
        cerr << "Indels are not yet supported.\n";
        exit(1);
        // run<TChar, TAllocConfig, EditDistance>(opt);
    }
    else
        run<TChar, TAllocConfig, HammingDistance>(opt);
}

template <typename TChar>
inline void run(Options & opt)
{
    if (opt.mmap)
        run<TChar, MMap<> >(opt);
    else
        run<TChar, Alloc<> >(opt);
}

int main(int argc, char *argv[])
{
    // Argument parser
    ArgumentParser parser("Mappabilty");
    addDescription(parser,
        "App for calculating the mappability values. Only supports Dna4/Dna5 so far.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory (error number, length and overlap will be appended to the output file)", ArgParseArgument::OUTPUT_FILE, "OUT"));
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
    Options opt;
    getOptionValue(opt.errors, parser, "errors");
    getOptionValue(opt.length, parser, "length");
    getOptionValue(opt.overlap, parser, "overlap");
    getOptionValue(opt.threads, parser, "threads");
    getOptionValue(opt.indexPath, parser, "index");
    getOptionValue(opt.outputPath, parser, "output");
    opt.mmap = isSet(parser, "mmap");
    opt.indels = isSet(parser, "indels");

    if (opt.overlap > opt.length - opt.errors - 2)
    {
        cerr << "ERROR: overlap should be <= K - E - 2\n";
        exit(1);
    }

    cout << mytime() << "Program started." << endl;

    CharString _indexPath;
    _indexPath = opt.indexPath;
    _indexPath += ".singleIndex";
    open(opt.singleIndex, toCString(_indexPath));

    _indexPath = opt.indexPath;
    _indexPath += ".alphabet";
    open(opt.alphabet, toCString(_indexPath));

    if (opt.alphabet == "dna4")
    {
        run<Dna>(opt);
    }
    else
    {
        // run<Dna5>(opt);
        cerr << "Dna5 alphabet has not been tested yet. Please do and remove this error message afterwards.\n";
        exit(1);
    }
}
