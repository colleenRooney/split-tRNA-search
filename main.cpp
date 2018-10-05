//---------------------------------------------------------------------
// Virtual Footprint
// Version: 0.70 alpha
// Description: Search tool of complex DNA patterns
// Author(s): Richard Muench
// Requirements: C++ (gcc version 3.3.1 was used)
// Email: r.muench@tu-bs.de
// Address: Technical University Braunschweig, Institute of Microbiology
//          Spielmannstr. 7, 38106 Brauschweig, GERMANY
// Web: http://www.prodoric.com (.net, .de)
// last change: 2004-07-10
// example call: ./vfp -p pattern.xml -s sequence.seq -o results.xml -v
//---------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include "footprint.h"

using namespace std;

string vfp_version="0.71";
string vfp_release="alpha release 2004-08-26";

bool file_exist(string filename) {
//---------------------------------------------------------------------------
// Check if file exists:
//---------------------------------------------------------------------------
  ifstream test(filename.data());
  if (!test) {
    return false;
  }
  else {
    return true;
    test.close();
  }
}

void showhead(unsigned char type) {
//---------------------------------------------------------------------------
// Print program head with command line options:
// type 0: Virtual Footprint (default)
//      1: tDNA search
//---------------------------------------------------------------------------
switch (type) {
case 0:
  cout <<"Virtual Footprint " << vfp_version << " (" << vfp_release << ")" << endl;
  cout <<"This program is part of the PRODORIC package http://www.prodoric.com (.net .de)" << endl;
  cout <<"(C) Richard Muench, 2001-2004, Technical University Braunschweig, GERMANY" << endl;
  cout <<"email: r.muench@tu-bs.de" << endl << endl;
  cout <<"Usage: vfp [options] [-p <pattern file>] [-s <sequence file>] [-o <output file>]" << endl;
  cout <<"-x               XML formatted output file (default)" << endl;
  cout <<"-c               CSV formatted output file" << endl;
  cout <<"-n               set non-occurance penalty in PWM" << endl;
  cout <<"-t <size>        truncated sequence <size in bp>" << endl;
  cout <<"-r               remove redundant palindromic matches" << endl;
  cout <<"-l               write log file(s)" << endl;
  cout <<"-v               verbose mode" << endl;
  cout <<"-show            show pattern information" << endl;
  cout <<"-?               this help" << endl;
  cout <<"--help           this help" << endl;
  cout << endl;
  break;
} //switch
}

int main(int argc, char *argv[]) {
//---------------------------------------------------------------------------
// Main program: Virtual Footprint
//---------------------------------------------------------------------------
string arg;                    //command line arguments
string pattern_file,
       sequence_file,
       output_file;            //in- and output files

long truncate=0;               //truncated sequence size (0: not truncated)
unsigned char file_format=0;   //resultfile format (0:XML, 1:CSV)
bool nop_mode=false;           //non-occurance penalty in PWM
bool verbose_mode=false;       //verbose mode (show messages)
bool log_mode=false;
bool redundancy=false;         //remove redundant (palindromic) matches

unsigned char action=1;        //type of action
//0: show program head
//1: do standard search
//2: show matrix information

//command line arguments handling:
if (argc<2) {
  showhead(0); //if no arguments print program head
  action=0;
} //if

else {

  for (unsigned int i=1; i<argc; i++) {
    arg=argv[i];
    if (arg=="-p") {
      pattern_file=argv[i+1];  //assign pattern filename
    }
    if (arg=="-s") {
      sequence_file=argv[i+1]; //assign sequence filename
    }
    if (arg=="-o") {
      output_file=argv[i+1];   //assign output filename
    }
    if (arg=="-x") {
      file_format=0;           //assign output fileformat
    }
    if (arg=="-c") {
      file_format=1;           //assign output fileformat
    }
    if (arg=="-v") {
      verbose_mode=true;       //set verbose mode
    }
    if (arg=="-l") {
      log_mode=true;           //write log files
    }
    if (arg=="-r") {
      redundancy=true;         //remove redundant palindromic matches
    }
    if (arg=="-n") {
      nop_mode=true;           //set non-occurance penalty in PWM
    }    
    if (arg=="-t") {
      truncate=atoi(argv[i+1]);//assign output filename
    }
    if (arg=="-show") {
      action=2;                //alter action number (show info)
    }
    if (arg=="-?") {
        showhead(0);           //print program head as help
        action=0;
    }
    if (arg=="--help") {
      showhead(0);             //print program head as help
      action=0;
    }
  } //for

} // else if

//check if pattern and sequence files are existing:
if ((file_exist(pattern_file))&&(file_exist(sequence_file))) {
  footprint vfp(sequence_file, pattern_file, nop_mode, verbose_mode, redundancy, truncate, log_mode); //create footprint object

  //start action:
  switch (action) {
    case 1: cout <<"Virtual Footprint " << vfp_version << " (" << vfp_release << ")" << endl;
         cout.flush();
         vfp.execute(output_file, file_format);
    break;

    case 2: vfp.showinfo();
    break;
  } //switch

} //if
  else {
    if (action) {
      showhead(0);
      cout << "File open error in pattern or sequence file" << endl;
      action=0;
    }

} //else if

}
