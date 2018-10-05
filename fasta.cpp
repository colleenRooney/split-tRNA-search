//------------------------------------------------------------
//fasta.cpp
//------------------------------------------------------------

#include "fasta.h"

fasta::fasta(){
//------------------------------------------------------------
// Constructor:
//------------------------------------------------------------
  type=1;      //DNA as default
  buflen=1000; //set file reading buffer
  number=0;    //set number of sequences=0 
}

fasta::fasta(short seq_type) {
//------------------------------------------------------------
// Constructor:
// Set type: 1: DNA; 2: protein
//------------------------------------------------------------
  type=seq_type;
  buflen=1000; //set file reading buffer
  number=0;    //set number of sequences=0 
}


fasta::fasta(short seq_type, string file) {
//------------------------------------------------------------
// Constructor:
// Set type: 1: DNA; 2: protein
// Open FASTA file
//------------------------------------------------------------
  type=seq_type;
  buflen=1000; //set file reading buffer
  fasta::open(file);
}

fasta::~fasta(){
//------------------------------------------------------------
// Destructor:
//------------------------------------------------------------
}

void fasta::open(string filename) {
//------------------------------------------------------------
// Open FASTA file:
// filename: filename of FASTA file
//------------------------------------------------------------
char buf[buflen]; //file reading buffer
int counter=-1;   //index of sequence array
int pos, old_pos; //position
const char *seq;  //pointer to sequence
string BC_str;    //base composition string in CSV format
string sep=",";   //separator for predefined base composition
string element;   //part of predefined base composition

//character sum A, C, G, T:
long ChrSumA=0; long ChrSumC=0; long ChrSumG=0; long ChrSumT=0;

ifstream InFile(filename.data()); //open file
  while (InFile.getline(buf,buflen)) {
    if (buf[0]=='>') {
      counter++;                //increase array index
      sequence[counter].ID=buf; //set ID
    } //if
    else {
      sequence[counter].seq+=buf;
      //sequence[counter].seq.append(buf);
    } //if
  } //while
InFile.close(); //close file


//calculate sequence composition of all sequences:
for (unsigned int n=0; n<=counter; n++) {

  //cut > from FASTA ID:
  sequence[n].ID=sequence[n].ID.substr(1,(sequence[n].ID.size()-1));

  pos=sequence[n].ID.find("$BC={",0); //position of base composition id
  //cout << "position: " << pos << endl;
  if (pos>=0) { //predefined base composition
    old_pos=sequence[n].ID.find("}",pos);
	BC_str=sequence[n].ID.substr(pos+5,(old_pos-pos-5));
	//cout << BC_str << endl;
    old_pos=0;
    for (unsigned int i=0; i<4; i++) {
      pos=BC_str.find(sep,old_pos);
      if (pos>=0) {
        element=BC_str.substr(old_pos,(pos-old_pos));
	    old_pos=pos+1;
      } else element=BC_str.substr(old_pos,(BC_str.size()-old_pos));
      //cout << element << endl;
	  sequence[n].c[i]=atoi(element.data());
    } //for i
  }
  else { //calculate base composition
    seq=sequence[n].seq.data(); //pointer to sequence
    for (unsigned int i=0; i<sequence[n].seq.size(); i++) {
      if (seq[i]==65) ChrSumA++;
      if (seq[i]==67) ChrSumC++;
      if (seq[i]==71) ChrSumG++;
      if (seq[i]==84) ChrSumT++;
    } //for i

    sequence[n].c[0]=ChrSumA;
    sequence[n].c[1]=ChrSumC;
    sequence[n].c[2]=ChrSumG;
    sequence[n].c[3]=ChrSumT;
  } //else if

} //for n

  /*
  charval[65]=0; //A
  charval[67]=1; //C
  charval[71]=2; //G
  charval[84]=3; //T
  */

file=filename;     //set class-variable: filename
number=counter+1;  //set class-variable: number of sequences

}

void fasta::add (string id, string seq) {
//------------------------------------------------------------
// Add new sequence
// id: FASTA ID
// seq: sequence
//------------------------------------------------------------
  sequence[number].ID=id;
  sequence[number].seq=seq;
  for (unsigned int i=0; i<4; i++) {
    sequence[number].c[i]=-1;
  } //for i
  number++;
}

void fasta::set (int SeqNo, string id, string seq) {
//------------------------------------------------------------
// Enter or modify sequence
// SeqNo: Sequence Number
// id: FASTA ID
// seq: squence
//------------------------------------------------------------
  sequence[SeqNo].ID=id;
  sequence[SeqNo].seq=seq;
  number=SeqNo+1;
}


int fasta::getnumber() {
//------------------------------------------------------------
// Get number of sequences   
//------------------------------------------------------------
  return number;
}


string fasta::getID (int SeqNo) {
//------------------------------------------------------------
// Return FASTA ID
// SeqNo: sequence number of FASTA array
//------------------------------------------------------------
  return sequence[SeqNo].ID;
}  

string fasta::getseq (int SeqNo) {
//------------------------------------------------------------
// Return sequence  
// SeqNo: sequence number of FASTA array
//------------------------------------------------------------
  return sequence[SeqNo].seq;
}   

string fasta::subseq(int SeqNo, int SeqPos, int SeqLength) {
//------------------------------------------------------------
// Return subsequence:
// SeqNo:     sequence number of FASTA array
// SeqPos:    start position
// SeqLEngth  sequence length
//------------------------------------------------------------ 
  string temp;
  temp=sequence[SeqNo].seq.substr(SeqPos,SeqLength);
  return temp.data();
}


string fasta::revcomp(int SeqNo, int SeqPos, int SeqLength) {
//------------------------------------------------------------
// Returns reverse complement DNA subsequence:
// SeqNo:     sequence number of FASTA array
// SeqPos:    start position
// SeqLEngth  sequence length
//------------------------------------------------------------ 
  string temp;
  string rc;
  string rc_out;
  //if no SeqLength then take whole sequence length:
  if (SeqLength==0) SeqLength=sequence[SeqNo].seq.size();
  temp=sequence[SeqNo].seq.substr(SeqPos,SeqLength);
  rc=temp.data();
  
  //make DNA complement
  for (unsigned int i=0; i<rc.size(); i++) {
    switch(rc.at(i)) {
      case 65: rc_out="T"+rc_out; //A->T
      break;
      case 67: rc_out="G"+rc_out; //C->G
      break;
      case 71: rc_out="C"+rc_out; //G->C
      break;
      case 84: rc_out="A"+rc_out; //T->A
      break;
      case 97: rc_out="t"+rc_out; //a->t
      break;
      case 99: rc_out="g"+rc_out; //c->g
      break;
      case 103: rc_out="c"+rc_out; //g->c
      break;
      case 116: rc_out="a"+rc_out; //t->a
      break;
     }; //switch
  }//for i
  return rc_out;
}

string fasta::uppercase(int SeqNo) {
//------------------------------------------------------------------
// Return uppercase string
// SeqNo: sequence number
//------------------------------------------------------------------
  char ch;
  string tempStr;

  for (unsigned int i=0; i<sequence[SeqNo].seq.size(); i++) {
    ch=sequence[SeqNo].seq.at(i);
    tempStr+=char(toupper(ch));
  }
  return tempStr;
}


long fasta::getSC(int SeqNo, int ChrNo) {
//---------------------------------------------------------------
// Return sequence composition of a certain sequence and base
// SeqNo:  sequence number
// ChrNo: sequence character DNA (A:0,C:1,G:2,T:3)
//---------------------------------------------------------------
return sequence[SeqNo].c[ChrNo];
}


float fasta::stacking_energy(int SeqNo, long start_pos, long end_pos, short ori, int up, int down, int window) {
//---------------------------------------------------------------
// Return the medium stacking energy of a DNA stretch
// start_pos:  sequence start_position
// end_pos: sequence end_position
// stacking energy in Kcal/mol bp
//---------------------------------------------------------------

string subseq;
string duplet;
long p1,p2,m;
float SE=0.0;
float SE_max=-9999.9;
char c1,c2;
//string *seq

//seq=sequence[SeqNo].seq.data();

m=start_pos+(long)((end_pos-start_pos)/2);
if (!ori) {
  p1=m-up;
  p2=m+down;
} else {
  p1=m-down;
  p2=m+up;
}
if (p1<0) p1=0;
if (p2>(sequence[SeqNo].seq.size()-1)) p2=sequence[SeqNo].seq.size()-1;

//cout << "p1:" << p1 << "p2:" << p2; cout.flush();

//subseq=sequence[SeqNo].seq.substr(p1,p2-p1+1);
if ((p2-p1+1)<window) {
  return -10.0;
} else {
for (unsigned int w=p1; w<=(p2-window+1); w++) {
  SE=0.0;
  for (unsigned int i=w; i<(w+window-1); i++) {

	c1=sequence[0].seq.at(i);
	c2=sequence[0].seq.at(i+1);

	if  ((c1==71)&&(c2==67))  SE+=-14.59; //GC
	if  ((c1==67)&&(c2==71))  SE+= -9.69; //CG
	if  ((c1==65)&&(c2==84))  SE+= -6.57; //AT
	if  ((c1==84)&&(c2==65))  SE+= -3.82; //TA
	if (((c1==65)&&(c2==67))||
	    ((c1==71)&&(c2==84))) SE+=-10.51; //AC|GT
	if (((c1==84)&&(c2==67))||
	    ((c1==71)&&(c2==65))) SE+= -9.81; //TC|GA
	if (((c1==71)&&(c2==71))||
	    ((c1==67)&&(c2==67))) SE+= -8.26; //GG|CC
	if (((c1==84)&&(c2==71))||
	    ((c1==67)&&(c2==65))) SE+= -6.57; //TG|CA
	if (((c1==65)&&(c2==71))||
	    ((c1==67)&&(c2==84))) SE+= -6.78; //AG|CT
	if (((c1==65)&&(c2==65))||
	    ((c1==84)&&(c2==84))) SE+= -5.37; //AA|TT

  }
  if (SE>SE_max) SE_max=SE;
}
return SE_max/(float)(window-1);
}

/*

$stacking_par= array("GC"=>-14.59,
                     "AC"=>-10.51, "GT"=>-10.51,
                     "TC"=> -9.81, "GA"=> -9.81,
                     "CG"=> -9.61,
                     "GG"=> -8.26, "CC"=> -8.26,
                     "AT"=> -6.57,
                     "TG"=> -6.57, "CA"=> -6.57,
                     "AG"=> -6.78, "CT"=> -6.78,
                     "AA"=> -5.37, "TT"=> -5.37,
                     "TA"=>  3.82);
*/

}
