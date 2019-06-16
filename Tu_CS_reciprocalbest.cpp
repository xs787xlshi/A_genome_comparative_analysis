/*NOV 21 2017
  Written by Xiaoli Shi
  Pick best much from  tblastn output .
  Usage: ./a.out 
*/

#include <cstdio>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <valarray>
#include <map>
#include <iomanip>
#include <vector>
#include <sstream>
#include <list>
using namespace std;

typedef map<string,string> StrStrMap;
typedef map<unsigned int, string> IntStrMap;
typedef multimap<unsigned int, string> IntStrMMap;
typedef multimap<string, string> StrStrMMap;

int hit(int exon_s,int exon_e,int read_s,int read_e)
{
  int len;
  len = min(exon_e, read_e) -max(exon_s, read_s) +1;
  return len;
}
    
void read_gff(StrStrMap& tu_dic);
void tu_best(StrStrMMap& tubest_dic, StrStrMap& tu_dic);
void CS_best(IntStrMMap& CSbest_dic, string& CSchr, StrStrMap& tu_dic);
void out_dic(StrStrMMap& tubest_dic, IntStrMMap& CSbest_dic, StrStrMap& tu_dic,string& CSchr);
void Tokenize(const string& str,vector<string>& tokens,const string& delimiters );
double string2double (const string& str);
int string2int (const string& str);

int main ()
{ 
  double total_read;
  vector <string> samfile;

  StrStrMap tu_dic;
  IntStrMMap CSbest_dic;
  StrStrMMap tubest_dic;
  

  int mi;
  string CSchr[] = {"1A","1B", "1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D","Un"}; 

  read_gff(tu_dic);
  /*
  for(StrStrMap::iterator rpos =tu_dic.begin();rpos != tu_dic.end();++ rpos)
    cout<<rpos->first<<"\t"<<rpos->second<<endl;
  exit(0);
  */
  tu_best(tubest_dic, tu_dic);
 
  /*
  for(StrStrMMap::iterator tpos =tubest_dic.begin();tpos != tubest_dic.end();++ tpos)
    cout<<tpos->first<<"\t"<<tpos->second<<endl;
  exit(0);
  */
  for (mi=0;mi<=21;mi++)
  {
      CSbest_dic.clear();
      CS_best(CSbest_dic,CSchr[mi], tu_dic);
      /*
      for(IntStrMMap::iterator bpos = CSbest_dic.begin();bpos != CSbest_dic.end();++ bpos)
	cout<<bpos->first<<"\t"<<bpos->second<<endl;
	exit(0);*/
      out_dic(tubest_dic, CSbest_dic,tu_dic,CSchr[mi]);
  }
  cout<<"end of all subroutine"<<endl;
  tu_dic.clear();
  CSbest_dic.clear();
}


void read_gff(StrStrMap &tu_dic)
{
  vector <string> rm_tokens;
  string infile,line,tuchr, temp_str;
  ifstream ifile;
  int rgk;
  stringstream rgconvert;
  
  for(rgk =1;rgk<=8;rgk++){
    if(rgk < 8)
      {
	rgconvert.str("");
	rgconvert << rgk;
	tuchr = rgconvert.str();
	infile = "/mnt/Data_diskD/diploid_tetraploid_hexaploid_Agenome/protein_comparison/MCScanTuEmmer/Tu_CDS_Emmer_DNA/Tu"+tuchr+"_gene_seq.out";
      }
    else 
      infile = "/mnt/Data_diskD/diploid_tetraploid_hexaploid_Agenome/protein_comparison/MCScanTuEmmer/Tu_CDS_Emmer_DNA/TuUngrouped_gene_seq.out";
    //cout<<infile<<endl;
    ifile.open(infile.c_str());
    if(!ifile)
      {
	cerr << "can't open input tu seq file."<<endl;
	exit(EXIT_FAILURE);
      }
    while(getline(ifile,line))
      {
	if(line.substr(0,1) == ">")
	  {
	    rm_tokens.clear();
	    Tokenize(line,rm_tokens,"\t");
	    //key = geneid; value = chr start end
	    if(string2int(rm_tokens.at(2)) > string2int(rm_tokens.at(3))){
	      temp_str = rm_tokens.at(2);
	      rm_tokens[2] = rm_tokens.at(3);
	      rm_tokens[3] = temp_str;
	    }
	    tu_dic.insert(make_pair(rm_tokens.at(0).substr(1,rm_tokens.at(0).size()-1),rm_tokens.at(1)+"\t"+rm_tokens.at(2)+"\t"+rm_tokens.at(3)));
	  }
      }
    ifile.clear();
    ifile.close();
  }
 /*cout<<"test unique_read content."<<endl;
  StrStrMap::iterator test_pos;
  for(test_pos = tu_dic.begin();test_pos != tu_dic.end();++test_pos){
    cout<<test_pos->first<<" "<<test_pos->second<<endl;
  }
  exit(0);*/
}

void tu_best(StrStrMMap& tubest_dic, StrStrMap& tu_dic)
{
  float left, right;
  int marker = 0;
  string line,infile,tuid,temp_str, temp_strand;
  vector <string> tokens1, tokens2;

  ifstream ifile;
  int tbi,tbk;
  
  string CSchr[] = {"1A","1B", "1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D","Un"}; 
  string tuchr;
  stringstream tbconvert;

  for (tbi=0;tbi<=21;tbi++)
    {
      for(tbk =1;tbk<=8;tbk++)
	{
	  if(tbk < 8)
	    {
	      tbconvert.str("");
	      tbconvert << tbk;
	      tuchr = tbconvert.str();
	      infile = "./Tu_CDS_Cs_DNA_alignment/CS_chr"+CSchr[tbi]+"_Tu"+tuchr+"_gene.blastn.outfmt6";
	    }
	  else 
	    infile = "./Tu_CDS_Cs_DNA_alignment/CS_chr"+CSchr[tbi]+"_TuUngrou.blastn.outfmt6";
	  //cout<<infile<<endl;
	  ifile.open(infile.c_str());
	  if(!ifile)
	    {
	      cerr << "can't open input blast1 file."<<endl;
	      exit(EXIT_FAILURE);
	    }
	  while(getline(ifile,line))
	    {
	      tokens1.clear();
	      Tokenize(line,tokens1,"\t");
	      if(string2double(tokens1.at(10)) > 0.00001 || string2double(tokens1.at(2)) < 30){continue;}
	      tuid = tokens1.at(0);
	      //tu_dic: key = geneid; value = chr start end (42,821 genes)
	      if(tu_dic.find(tuid) == tu_dic.end()){continue;}
	      temp_strand = "+";
	      if(string2int(tokens1.at(8)) > string2int(tokens1.at(9)))
		{
		  temp_str = tokens1.at(8);
		  tokens1[8] = tokens1.at(9);
		  tokens1[9] = temp_str;
		  temp_strand = "-";
		}
	      
	      marker = 0;
	      //check if the tubest_dic has contained a match for this tu gene
	      if(tubest_dic.find(tuid) != tubest_dic.end())
		{
		  for (StrStrMMap::iterator it = tubest_dic.lower_bound(tuid); it != tubest_dic.upper_bound(tuid);++it )
		    {
		      tokens2.clear();
		      Tokenize(it->second,tokens2,"\t");
		      //score threshold: [score-score*5%, score+score*5%]
		      left = string2double(tokens2.at(11))-string2double(tokens2.at(11))*.05;
		      right = string2double(tokens2.at(11))+string2double(tokens2.at(11))*.05;
		      if(string2double(tokens1.at(11)) > right)
			{//this is better than the deposited old matches 
			  marker = 1;
			  break;
			}
		      else if(string2double(tokens1.at(11)) < left)
			{//this is worse than the deposited old matches
			  marker = -1;
			  break;
			}
		    }
		}
	      if(marker == 1)
		{
		  tubest_dic.erase(tuid);
		}
	      if(marker > -1)
		{//put new match to tubest_dic: key = tuid; value = alignment line with query start < query end.
		  temp_str =  tokens1.at(0)+"\t"+tokens1.at(1)+"\t"+tokens1.at(2)+"\t"+tokens1.at(3)+"\t"+tokens1.at(4)+"\t"+tokens1.at(5)+"\t"+tokens1.at(6)+"\t"+tokens1.at(7)+"\t"+tokens1.at(8)+"\t"+tokens1.at(9)+"\t"+tokens1.at(10)+"\t"+tokens1.at(11)+"\t"+temp_strand;
		  tubest_dic.insert(make_pair(tuid, temp_str));
		}
	    }
	  ifile.clear();
	  ifile.close();  
	}
    }
}

void CS_best(IntStrMMap& CSbest_dic,string& CSchr,StrStrMap& tu_dic)
{
  int line_count =0,marker = 0,over_len,over_tu,ref_len1,ref_len2,rbk;
  string line,infile,temp_str,tuchr, temp_strand;
  vector <string> tokens1, tokens2, tu_tokens1, tu_tokens2;
  ostringstream convert;
   
  ifstream ifile;
  for(rbk =1;rbk<=8;rbk++)
    {
      if(rbk < 8){
	convert.str("");
	convert << rbk;
	tuchr = convert.str();
	infile = "./Tu_CDS_Cs_DNA_alignment/CS_chr"+CSchr+"_Tu"+tuchr+"_gene.blastn.outfmt6";
      }
      else{
	infile = "./Tu_CDS_Cs_DNA_alignment/CS_chr"+CSchr+"_TuUngrou.blastn.outfmt6";
      }
      ifile.open(infile.c_str());
      if(!ifile)
	{
	  cerr << "can't open input blast2 file."<<endl;
	  exit(EXIT_FAILURE);
	}
      cout<<infile<<endl;
      while(getline(ifile,line))
	{
	  tokens1.clear();
	  Tokenize(line,tokens1,"\t");
	  if(string2double(tokens1.at(10)) > 0.00001 || string2double(tokens1.at(2)) < 30){continue;}
	  temp_strand = "+";
	  if(string2int(tokens1.at(8)) > string2int(tokens1.at(9))){
	    temp_str = tokens1.at(8);
	    tokens1[8] = tokens1.at(9);
	    tokens1[9] = temp_str;
	    temp_strand = "-";
	  }
	  //check if the CS dna segment is overlaped with exist pairs
	  line_count++;
	  //cout << line_count << endl;
	  //cout << line << endl;
	  //CSbest_dic: key= start;value= alignment1+"\n"+alignment2+....
	  for (IntStrMMap::iterator it = CSbest_dic.begin(); it != CSbest_dic.end();)
	    {
	      if(it->first > string2int(tokens1.at(9)))
		break;
	      tokens2.clear();
	      Tokenize(it->second,tokens2,"\t");
	      if(string2int(tokens2.at(9)) < string2int(tokens1.at(8))){
		++it;
		continue;
	      }
	      over_len = hit(string2int(tokens1.at(8)),string2int(tokens1.at(9)),string2int(tokens2.at(8)),string2int(tokens2.at(9)));
	      //80% of the longer segment
	      ref_len1 = 0.8*max(string2int(tokens1.at(9))-string2int(tokens1.at(8)),string2int(tokens2.at(9))-string2int(tokens2.at(8)));
	      //95% of the shorter segment
	      ref_len2 = 0.95*min(string2int(tokens1.at(9))-string2int(tokens1.at(8)),string2int(tokens2.at(9))-string2int(tokens2.at(8)));
	      //cout<<over_len<<"==\t=="<<ref_len1<<"==\t=="<<ref_len2<<endl;
	      marker = 0;
	      if(over_len > ref_len1 || over_len > ref_len2)//overlap of the two segments must longer than 80% of the longest segment or longer than 95% of the shortest segment. 
		{
		  //Some tu genes have overlap, in this case both alignment should be kept
		  tu_tokens1.clear();
		  Tokenize(tu_dic[tokens1.at(0)], tu_tokens1, "\t");
		  tu_tokens2.clear();
		  Tokenize(tu_dic[tokens2.at(0)], tu_tokens2, "\t");
		  
		  if (tu_tokens1.at(0) == tu_tokens2.at(0) && hit(string2int(tu_tokens1.at(1)),string2int(tu_tokens1.at(2)),string2int(tu_tokens2.at(1)),string2int(tu_tokens2.at(2))) > min(ref_len1, ref_len2)){
		    ++it;
		    continue;
		  }
		  if(string2double(tokens1.at(11)) > string2double(tokens2.at(11)))
		    { //new is better than old, delete old line
		      CSbest_dic.erase(it++);
		      marker = 1;
		      continue;
		    } 
		  else if(string2double(tokens1.at(11)) == string2double(tokens2.at(11)))
		    {//if scores are equal, check identity
		      if(string2double(tokens1.at(2)) > string2double(tokens2.at(2)))
			{//new is better than old
			  CSbest_dic.erase(it++);
			  marker = 1;
			  continue;
			}
		      else if(string2double(tokens1.at(2)) == string2double(tokens2.at(2)))
			marker = 2;
		      else
			marker = -1;
		    } 
		  else
		    marker = -1;
		}
	      if(marker == -1){
		break;
	      }//find an exist element better than the new one, go out from the CSbest_dic loop
	      ++it;
	    }//end of for loop of map
	  //marker == 0: no overlap; marker == 1: better than the exist overlp item; marker == 2:equvalent with the exist overlap item; marker == -1: find a overlap better than the new one
	  if(marker != -1){
	    temp_str = tokens1.at(0)+"\t"+tokens1.at(1)+"\t"+tokens1.at(2)+"\t"+tokens1.at(3)+"\t"+tokens1.at(4)+"\t"+tokens1.at(5)+"\t"+tokens1.at(6)+"\t"+tokens1.at(7)+"\t"+tokens1.at(8)+"\t"+tokens1.at(9)+"\t"+tokens1.at(10)+"\t"+tokens1.at(11)+"\t"+temp_strand;
	    CSbest_dic.insert(make_pair(string2int(tokens1.at(8)),temp_str));   
	  }
	  marker = 0;
	}
      ifile.clear();
      ifile.close();
    }
}


void out_dic (StrStrMMap& tubest_dic, IntStrMMap& CSbest_dic,StrStrMap& tu_dic,string& CSchr)
{ 
  vector <string> tokens;
  StrStrMMap::iterator bpos;

  string outname = "Tu_CS"+CSchr+"_reciprocalbest.out";
  ofstream ofile(outname.c_str());
  if(!ofile)
  {
      cerr << "can't open output file"<<endl;
      exit(EXIT_FAILURE);
  }
  
  for(bpos = tubest_dic.begin();bpos != tubest_dic.end();++ bpos)
  {
      tokens.clear();
      Tokenize(bpos->second,tokens,"\t");
      //if alignments not occurred on CSchr, go to next
      if (tokens.at(1) != "chr"+CSchr)
	continue;
      if (CSbest_dic.find(string2int(tokens.at(8))) == CSbest_dic.end()){
	//cout<<bpos->second<<endl;
	//cout<<tokens.at(8)<<endl;
	continue;
      }
      for (IntStrMMap::iterator pos = CSbest_dic.lower_bound(string2int(tokens.at(8))); pos != CSbest_dic.upper_bound(string2int(tokens.at(8)));++pos )
	{
	  if(pos->second == bpos->second)
	    ofile<<tokens.at(0)+"\t"+tu_dic[tokens.at(0)]+"\t"+tokens.at(1)+"\t"+tokens.at(8)+"\t"+tokens.at(9)+"\t"+tokens.at(2)+"\t"+tokens.at(3)+"\t"+tokens.at(10)+"\t"+tokens.at(11)+"\t"+tokens.at(12)<<endl;
	}
  }
  ofile.clear();
  ofile.close();
}


int string2int (const string& str)
{
  int num1=0,num2=0,strlen;
  int i;
  strlen = str.length();

  for(i=0;i<strlen;i++)
    {
      num2 = str[i]-'0';
      num1 = num1*10+num2;
    }
  return num1;
}

void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters )
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

double string2double (const string& str)
{
  double num1=0,num2=0;
  int i,j,strlen,counter=-1;
  strlen = str.length();

  for(i=0;i<strlen;i++)
    {
      //cout<<"str[i]="<<str[i]<<endl;
      if(str[i] == 'e')
	{
	  num1 = num1+num2;
	  num2 = 0;
	  
	  for(j=i+2;j<strlen;j++)
	    {
	      num2 = (num2*10)+(str[j]-'0');
	    }
	  //cout<<num1<<"*"<<num2<<endl;
	  if(str[i+1] == '+')
	    num1 = num1*pow10(num2);	  
	  else if(str[i+1] == '-')
	    num1 = num1*pow10(num2*(-1));
	  else
	    {
	      cout<<"something wrong with the value:"<< str<<endl;
	      exit(EXIT_FAILURE);
	    }
	  return num1;
	}
      else if (str[i] == '.')
	{
	  counter++;
	}
      else if (counter > -1)
	{
	  counter++;
	  num2 = num2+(str[i]-'0')*pow10(counter*(-1));
	}
      else
	{
	  num1 = (num1*10)+(str[i]-'0');
	}
    }//end of i cycle.
  num1 = num1 + num2;
  return num1;
}

int max(int a, int b)
{
  if(a > b)
    return a;
  else
    return b;
}


int min(int a, int b)
{
  if(a > b)
    return b;
  else
    return a;
}


