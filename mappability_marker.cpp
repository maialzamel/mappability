#include <vector>
#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;


void marker(CharString const & inputPath1, CharString const & inputPath2)
{

int number;
vector <int> file_size;
 vector<int> v1;
int prev_size=0;

for(int i=1;i<=30;i++){
 ifstream file(toCString("./outputdump/"+std::to_string(i)+".fasta"), std::ios::binary);

if (!file.eof() && !file.fail())
{


while(file >> number)
{
   
v1.push_back(number);



}
file.close();
if (i==1){
file_size.push_back(v1.size());
prev_size=v1.size();

cout << v1.size()<<endl;
}
else{
file_size.push_back(v1.size()+prev_size);
prev_size+=v1.size();
cout <<v1.size()<<endl;

}

v1.clear();

}

}


            

for (int i=0;i< file_size.size();i++){


cout << file_size.at(i)<<endl;
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
  
  cout <<file_size.size()<<endl;
  
   
         int i=0;
        cout << "number of runs\n";

        for (int j=0; j< v.size();j++){

            while (v[i+j]< v.size() && v[i+j]==contacted_v[i+j]){
             
        i++;
            }//end while 



        if (i==100){




          for (int k=0;k<=29;k++){
                 if(k==0)
                if (j >= 0 && j < file_size.at(k) ){
                cout <<j;
                cout << " genome 1"<< endl; 

               cout << " Pos "<<j << endl; 
                                                   }
                if (k==1)

                 if (j >= file_size.at(k-1) && j < file_size.at(k) ){
        
                 cout << " genome "<<k+1 <<" "<< file_size.at(k) << endl; 

                  cout << " Pos "<<j-file_size.at(k-1) << endl; 

                                                                    }
                                      }//end loop
 
              }//end if


        j+=i;
        i=0;

      
}//end big for 



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
