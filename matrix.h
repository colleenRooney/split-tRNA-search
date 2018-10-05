//------------------------------------------------------------------
// matrix.h
//------------------------------------------------------------------

// edited 2018 Colleen Rooney
#include <stdlib.h>
#include <math.h>
// end edit
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdbool.h>
#include <vector>
#include <algorithm>
#include "fasta.h"

using namespace std;

class matrix {


  struct mweight {               //individual matrix weights:
    int number;                  //absolute base number
    float freq;                  //frequency
    float number_bias;           //bias corrected base number
    float freq_bias;             //bias corrected frequency
    float score;                 //scores
  };

private:
  int buflen;                    //file reading buffer length
  int maxsites; 				 //maximum number of sites


public:



  //weight matrix content
  fasta site;                    //matrix sites (fasta class)
  mweight value[50][20];         //matrix weight 2D array (cols x rows)
  int colsum[50];                //matrix column sum of sequence characters
  float colsum_bias[50];         //bias corrected matrix column sum of sequence characters
  float Rs[50];                  //information content at a certain matrix position
  float Rs_bias[50];             //same as Rs, with consideration of the sequence bias
  float Ci[50];                  //Ci(l) vector (MatInspector)
  float Rsequence;               //matrix information content
  float Rsequence_bias;
  float maxscore;                //maximal reachable score of a weight matrix
  float minscore;                //minimal reachable score of a weight matrix
  vector<float> site_score;      //scores of individual sites comprising the weight matrix
  vector<float> core_score;      //core_scores of individual sites comprising the weight matrix
  vector<int>   core;            //most conserved positions
  float Pb[20];                  //frequency of sequence characters

  //pattern parameters
  unsigned char msize;           //weight matrix size (columns)
  string mname;                  //weight matrix name
  unsigned char model;           //weight matrix model
  float  msens;                  //overall weight matrix sensitivity
  float  csens;                  //core sensitivity
  unsigned char core_size;       //number of positions defining the core
  float  mth;                    //weight matrix score threshold
  float  cth;                    //core score threshold
  string consensus;              //IUPAC consensus string
  string regex;                  //regular expression
  int mism;                      //number of mismatches
  int nop;                       //non-occurrence panelty

  unsigned char  charval[90];    //assigns ascii value of sequence letter to weight matrix array
  unsigned char rcharval[90];    //assigns reverese sequence letter to weight matrix array


  matrix();
  ~matrix();

  void addseqcomp(long sc[],int num);       //enter sequence character composition
  void init(string i_name,
            unsigned char i_model,
			float i_msens, float i_csens,
			unsigned char i_core_Size,
			float i_score, int i_nop);      //calculate matrix values from sites and parameters
  float getRsequence();                     //returns the information content of the site (Rsequence)
  void getCore(unsigned char core_number);
  float getScore(const char *seq,
                 int offset, int dir);      //returns the Score of a sequence
  float getCoreScore(const char *seq,
                     int offset, int dir);
  void getSiteScores();
  int getMismatch(const char *seq, int offset, int dir, int mm_max);  //returns the number of mismatches to a IUPAC consensus
  void addsite(string seq, string ID);      //add site to matrix
  float calehnb (int n);                    //calculates small sample size correction
  void cons2matrix(string cons);            //calculate matrix content from consensus
  void show();              //print matrix contents and features
  float getthreshold(unsigned char type,
                     float sens);           //get threshold (core-)score from sensitifity
  void hist(string histfile,
            int Smin, int Smax, float step);//create histogram file

  float ilog(float n, float a);    //logarithm of base 2 where log2(0)=log2(1/(a+2))
  float elog(float n);           //natural logarithm where ln(0)=0;
  float round(float n, int d);



protected:

};

