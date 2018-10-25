#include <vector>

#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;

template <typename value_t>
void dump(CharString const & inputPath, CharString const & outputPath)
{
    vector<value_t> v;
    ifstream file(toCString(inputPath), std::ios::binary);
    if (!file.eof() && !file.fail())
    {
        file.seekg(0, std::ios_base::end);
        std::streampos fileSize = file.tellg();
        std::cout << "File size: " << fileSize << '\n';
        std::cout << "Sizeof: " << sizeof(value_t) << '\n';
        v.resize(fileSize / sizeof(value_t));
        file.seekg(0, std::ios_base::beg);
        file.read((char*) & v[0], fileSize);
        file.close();

        cout << "Load successful\n";

        cout << "v.size(): " << v.size() << '\n';

        // TODO: progress bar
        ofstream outfile(toCString(outputPath), std::ios::out | std::ofstream::binary);
        copy(v.begin(), v.end(), (std::ostream_iterator<value_t>(outfile), std::ostream_iterator<int>(outfile, " ")));

        cout << "Done.\n";

        return;
    }

    cout << "Something went wrong ...\n";
}

int main(int argc, char *argv[])
{
    // Argument Parser
    // TODO: allow more than 4 gigabases
    ArgumentParser parser("Mappability Dumper");
    addDescription(parser, "Transforms the output of the mappability program into human readable format. Please be aware that the file size will increase significantly and most likely multiply.");

    addOption(parser, ArgParseOption("I", "input", "Path to the mappability file", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "input", "gmapp8 gmapp16");
	setRequired(parser, "input");

    addOption(parser, ArgParseOption("O", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString inputPath, outputPath;
    getOptionValue(inputPath, parser, "input");
    getOptionValue(outputPath, parser, "output");

    string _inputPath = toCString(inputPath);

    if (_inputPath.substr(_inputPath.find_last_of(".") + 1) == "gmapp8")
        dump<uint8_t>(inputPath, outputPath);
    else
        dump<uint16_t>(inputPath, outputPath);

    // {
    //     vector<uint8_t>  v1 {0, 1, 255, 255};
    //     vector<uint16_t> v2 {0, 1, 255, 256, 16000, 65535};
    //
    //     CharString v1_path = "/dev/shm/v1";
    //     CharString v2_path = "/dev/shm/v2";
    //
    //     ofstream outfile1(toCString(v1_path), ios::out | ios::binary);
    //     outfile1.write((const char*) &v1[0], v1.size() * sizeof(uint8_t));
    //     outfile1.close();
    //
    //     ofstream outfile2(toCString(v2_path), ios::out | ios::binary);
    //     outfile2.write((const char*) &v2[0], v2.size() * sizeof(uint16_t));
    //     outfile2.close();
    // }

    return 0;
}
