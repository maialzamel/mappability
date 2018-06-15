#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include "common.h"

using namespace seqan;

template <typename TString, typename TStringSetConfig>
void buildIndex(StringSet<TString, TStringSetConfig> const & chromosomes, CharString const & indexPath)
{
    // TODO(cpockrandt): avoid copying
    StringSet<TString, Owner<ConcatDirect<> > > chromosomesConcat(chromosomes);
    Index<StringSet<TString, Owner<ConcatDirect<> > >, TIndexConfig> index(chromosomesConcat);
    indexCreate(index, FibreSALF());
    save(index, toCString(indexPath));
}

template <typename TString>
void buildIndex(TString & chromosomes, CharString const & indexPath)
{
    Index<TString, TIndexConfig> index(chromosomes);
    indexCreate(index, FibreSALF());
    save(index, toCString(indexPath));
}

template <typename TText>
void run(StringSet<CharString> const & ids, TText const & chromosomes, bool const singleIndex, CharString const & indexPath)
{
    if (singleIndex/* && length(chromosomes) > 1*/)
    {
        // Build and save index
        buildIndex(chromosomes, indexPath);
    }
    else
    {
        // Chromosome by chromosome
        for (unsigned i = 0; i < length(chromosomes); ++i)
        {
            // Build and save index
            std::string _indexPath = toCString(indexPath);
            // if (length(chromosomes) > 1)
            _indexPath += "." + std::to_string(i);
            buildIndex(chromosomes[i], _indexPath);
            std::cout << ids[i] << " indexed.\n" << std::flush;
        }
    }
    CharString _indexPath = indexPath;
    _indexPath += ".singleIndex";
    save(singleIndex/* && length(chromosomes) > 1*/, toCString(_indexPath));
}

int main(int argc, char *argv[])
{
    // Argument Parser
    ArgumentParser parser("Index Creation");
    addDescription(parser, "App for creating an index. Only supports Dna (with and without N's). At most 4 gigabases in total allowed.");

    addOption(parser, ArgParseOption("G", "genome", "Path to the genome", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "genome", "fa fasta fastq");
	setRequired(parser, "genome");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("s", "single-index", "Builds a single index instead of one for each chromosome."));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString indexPath, genomePath;
    getOptionValue(indexPath, parser, "index");
    getOptionValue(genomePath, parser, "genome");
    bool singleIndex = isSet(parser, "single-index");

    // Read fasta input file
    StringSet<CharString> ids;
    StringSet<Dna5String> chromosomes;
    SeqFileIn seqFileIn(toCString(genomePath));
    readRecords(ids, chromosomes, seqFileIn);
    std::cout << "Number of sequences: " << length(chromosomes) << '\n';

    // check whether it can be converted to Dna4
    unsigned i = 0;
    for (; i < length(chromosomes); ++i)
    {
        unsigned j = 0;
        for (; j < length(chromosomes[i]); ++j)
        {
            if (chromosomes[i][j] == 'N')
                break;
        }
        if (j < length(chromosomes[i]) && chromosomes[i][j] == 'N')
            break;
    }
    bool canConvert = i == length(chromosomes);
    std::cout << "Index will be constructed using Dna" << (5 - canConvert) << " alphabet.\n" << std::flush;

    if (canConvert)
    {
        StringSet<DnaString> chromosomes4(chromosomes);
        clear(chromosomes);
        run(ids, chromosomes4, singleIndex, indexPath);
    }
    else
    {
        run(ids, chromosomes, singleIndex, indexPath);
    }
    CharString _indexPath = indexPath;
    _indexPath += ".alphabet";
    CharString alphabet = canConvert ? "dna4" : "dna5";
    save(alphabet, toCString(_indexPath));

    _indexPath = indexPath;
    _indexPath += ".ids";
    save(ids, toCString(_indexPath));
    std::cout << "Index created successfully.\n";

    return 0;
}
