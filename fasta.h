//------------------------------------------------------------------
// fasta.h
//------------------------------------------------------------------

// edited 2018 Colleen Rooney
#include <stdlib.h>
// end edit
#include <fstream>
#include <iostream>

using namespace std;

class fasta {

private:

public:

  struct seqtype {  //record of sequences
    string ID;      //sequence ID
    string seq;     //sequence
    long c[20];    //sequence composition (A;C;G;T) or amino acids
  };
 
  seqtype sequence[3000];            //array of FASTA sequences
  short type;                       //1: DNA; 2: protein
  int number;                       //number of sequences
  string file;                      //fasta path&filename
  int buflen;                       //file reading buffer length
  //float SE_array[8000000];

    
  fasta();                          //default constructor DNA
  fasta(short seq_type);            //1: DNA; 2: protein; 3: RNA
  fasta(short seq_type,
        string file);               //sequence type and FASTA filename

  ~fasta();                         //destructor

  void open(string filename);       //open and read fasta file
  
  void add (string id, string seq); //add new sequence
  
  void set (int SeqNo, 
            string id, string seq); //enter or modify sequence
  
  int getnumber();                  //get number of sequences   
  
  string getID (int SeqNo);         // return FASTA ID
  
  string getseq (int SeqNo);        // return sequence  
   
  string subseq(int SeqNo,
                int SeqPos,
                int SeqLength);     //return subsequence

  string revcomp(int SeqNo,
                 int SeqPos,
                 int SeqLength);    //return reverse complement subsequence
  
  string uppercase(int SeqNo);      //return uppercase sequence               
                   
  long getSC(int SeqNr, int ChrNr);    //return sequence composition of a certain sequence

  float stacking_energy(int SeqNo,
                        long start_pos,
						long end_pos,
						short ori,
						int up,
						int down,
						int window);   //return medium stacking energy of a subsequence

protected:

};
