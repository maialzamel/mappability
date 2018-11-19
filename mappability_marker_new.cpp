#include <vector>
#include <seqan/arg_parse.h>
#include <fstream>
using namespace std;
using namespace seqan;


void marker(CharString const & inputPath1, CharString const & inputPath2)
{


 vector<pair <string,int>> v1;

 vector<pair <string,int>> v2;
string chr_info, chr_info2;
char singleCharacter;
int pos=0;
 ifstream file2(toCString("./outputdump/GCF_000177135.fasta"), std::ios::binary);
 ifstream file(toCString("./data/GCF_000177135.fasta"), std::ios::binary);

if (!file.eof() && !file.fail())
{

while(file.get(singleCharacter))
{
if (singleCharacter=='>'){
   getline(file,chr_info2);
v1.push_back(make_pair(chr_info2,-1));}
else{
string s;
if (singleCharacter != '\n'){
s.append(1, singleCharacter);
v1.push_back(make_pair(s,pos));
pos++;
}

}

}
file.close();
}//end if 

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
pos=0;

if (!file2.eof() && !file2.fail())
{

while(file2 >>chr_info)
{

v2.push_back(make_pair(chr_info,pos));
pos++;


}
file2.close();
}

int j=0;
for (int i=0; i<v1.size();i++){

cout << v1.at(i).second << " ";
cout << v1.at(i).first <<  " ";
if (!(v1.at(i).second  == -1)&& i <v2.size()){
cout << v2.at(j).second << " ";
cout << v2.at(j++).first ;
}
cout <<endl;
}
cout << v1.size() << " "<< v2.size()<<endl;

}

int main(int argc, char *argv[])
{
    
    ArgumentParser parser("Mappability Marker");
    addDescription(parser, "**");

    addOption(parser, ArgParseOption("I", "input", "Path to the mappability file", ArgParseArgument::INPUT_FILE, "IN"));

	setRequired(parser, "input");

    addOption(parser, ArgParseOption("O", "compared", "Path to compared file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "compared");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
     CharString inputPath, comparedPath;
     getOptionValue(inputPath, parser, "input");
    getOptionValue(comparedPath, parser, "compared");
    string _inputPath = toCString(inputPath);
   string _inputPath2 = toCString(comparedPath);

    marker(inputPath,comparedPath);
    return 0;
}
