#include <vector>
#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;


void marker(CharString const & inputPath1, CharString const & inputPath2)
{

int number;
vector <int> file_size;
 vector<int> v1;

for(int i=0;i<=31;i++){
 ifstream file(toCString("./outputdump/"+std::to_string(i)+".fasta"), std::ios::binary);

if (!file.eof() && !file.fail())
{


while(file >> number)
{
   
v1.push_back(number);



}
file.close();
file_size.push_back(v1.size());
v1.clear();

}

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ifstream original_file(toCString(inputPath1), std::ios::binary);
     ifstream contacted_file(toCString(inputPath2), std::ios::binary);

vector<int> contacted_v;
  vector<int> v;
if (!contacted_file.eof() && !contacted_file.fail())
{


while(contacted_file >> number)
{
 
contacted_v.push_back(number);

}
contacted_file.close();


}



if (!original_file.eof() && !original_file.fail())
{


while(original_file >> number)
{
  
v.push_back(number);

}
original_file.close();

}
  
  
  
   
         int i=0;
        cout << "number of runs\n";

        for (int j=0; j< v.size();j++){

            while (v[i+j]< v.size() && v[i+j]==contacted_v[i+j]){
             
        i++;
            }//end while 



        if (i>100){

            int genome_pos=0;

           for(int k=0; k<=31;k++){
            if (j>=genome_pos && j<= (file_size.at(k)+genome_pos)){
              if (k==0)
               cout <<"  genome "<< k <<" Starting Pos of Marker: " << (j) <<" run's len "<< i <<'\n';
//else
//cout <<"  genome "<< k <<" Starting Pos of Marker: " << (j)-file_size.at(k-1) <<" run's len "<< i <<'\n';

              break;
                 }//end if 
            else{  
  
               genome_pos+=file_size.at(k);


}
                                  }//end loop


       
 
                         }//end if


        j+=i;
        i=0;

      
}
int sum=0;
for (int i=0; i< file_size.size();i++){
sum+=file_size.at(i);
}


cout << v.at(v.size()-1) << " "<< contacted_v.at(contacted_v.size()-1)<<" "<<sum<<endl;

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
