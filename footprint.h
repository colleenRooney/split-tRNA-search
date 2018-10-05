//------------------------------------------------------------------
// footprint.h
//------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <regex.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>
#include <cctype>
#include "matrix.h"

using namespace std;

class footprint {

struct pattern_type {

  //general paramaters:
  unsigned char type;        //1: weight matrix
                             //2: IUPAC consensus
                             //3: regular expression
                             //4: spacer

  int no;                    //references the index of the subpattern

  unsigned char orientation; //0: restricted to pattern orientation
                             //1: allow both orientations
                             
  string name;               //name of subpattern                             
                             
  //weight matrix parameters:                           
  unsigned char model;       //scoring model
  float msens;               //weight matrix sensitivity
  float csens;               //core sensitivity
  float mth;                 //matrix threshold
  float cth;                 //core threshold
  unsigned char core;        //core size (0: no core)
  
  //IUPAC consensus paramaters:
  unsigned char mism;        //number of mismatches

};

struct spacer_type {
  //spacer parameters:
  int min_space;   //minimal space
  int max_space;   //maximal space
};

struct result_type {
string  id;
  int   no;                         //index of an individual subpattern (PWM, IUPAC, RegEx)
  int   pattern_no;                 //number of the subpattern
  string name;                      //name of the subpattern
  int   type;                       //type of subpattern (1: PWM, 2: IUPAC, 3: RegEx)
  long  start_pos;                  //absolute start position of (sub)pattern
  long  end_pos;                    //absolute end position of (sub)pattern
  short orientation;                //(sub)pattern orientation (0/1: +/- strand
  string seq;                       //(sub)pattern sequence
  float score;                      //weight matrix score
  float core_score;                 //core score
};

struct series_type{
  long  start_pos;
  long  end_pos;
  short orientation;
  string seq;
  float score;                      //overall score
  float se_score;                   //stacking energy
  vector<int> result;
  int last;
};

private:

  struct result_order{ //define compare operator for result_type
    bool operator()(const result_type &a, const result_type &b)  {
	  return a.start_pos < b.start_pos; //order by start position (ascending)
    }
  };

  struct series_order{ //define compare operator for result_type
    bool operator()(const series_type &a, const series_type &b)  {
	  return a.start_pos < b.start_pos; //order by start position (ascending)
    }
  };

  struct series_order_score{ //define compare operator (score) for result_type
    bool operator()(const series_type &a, const series_type &b)  {
	  return a.score > b.score; //order by score (descending)
    }
  };
  
  string patternXML;         //XML raw data from pattern file
  int patternXML_pos;        //cursor during parsing the XML file

  string cons2regex(string cons, int mm);
  string strval(int a);

public:

  string sequence_file, pattern_file, output_file;

  vector<pattern_type> pattern;  //pattern rules vector
  vector<spacer_type>  spacer;   //spacer vector
  vector<result_type>  result;   //result vector of positive matches
  vector<series_type>  series;   //defines patterns as series of subpatterns
  vector<matrix>       mat;      //vector of weight matrices

  int pattern_size;              //number of subpatterns
  float cum_threshold;           //cumulative threshold
  float se_threshold;            //stacking energy threshold
  
  fasta genome;                  //DNA sequence

  string cons[10];               //IUPAC consensi

  bool verbose_mode;             //true: show massages
  bool redundancy_mode;          //true: remove redundant (palindromic) matches
  bool nop_mode;                 //true: non-occurance penalty in PWM
  bool log_mode;                 //true: write log files

  footprint();                                    //empty constructor

  footprint(string sf, string pf,
            bool nop, bool mode, 
	    bool redundancy,
	    long truncate,
	    bool logm);                           //constructor

 ~footprint();                                    //destructor

  string upper(string s);                         //convert string to uppercase string

  time_t showtime();                              //show timer

  void verbose_cout(string cout_str);             //print messages in verbose mode

  string getnexttag (int pos);                    //return next tag

  string getattrib (string tag, string att_name); //return tag attribute

  string getXMLdata (int pos);                    //return XML data between two tags

  void openpattern(string patternfile, int seqno);//open searchpattern

  void showinfo();                                //show pattern information

  string showCSVmatch(int listnr, string sep);    //return match as CSV result line

  string showXMLmatch(int listnr);                //return match as XML tag

  void execute(string resultfile,
               unsigned char file_format);        //execute pattern search

  void saveresults(string resultfile,
                   int seq_nr,
                   unsigned char format);         //save results of hits

  long analyze(int seq_nr);                       //analyse subpatterns
  
  void remove_lower_scores();
  
  void remove_redundancies();                     //remove redundant palindromic matches

  long msearch(int mnr, int gnr, int pnr,
               float Smin, float Cmin, int core); //weight matrix search

  long csearch(int mnr, int gnr, int pnr, int mm);                  //consensus matrix search:

  long regexsearch(int mnr, int gnr, int pnr);              //regular expression search

protected:

};
