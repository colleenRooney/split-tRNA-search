//---------------------------------------------------------------------------
// footprint.cpp
//---------------------------------------------------------------------------

#include "footprint.h"

footprint::footprint() {

}

footprint::footprint(string sf, string pf, bool nop, bool mode, bool redundancy, long truncate, bool logm){
//---------------------------------------------------------------------------
// Constructor:
// sf: sequence file
// pf: pattern file
// nop: non-occurance panalty
// mode: true: verbose mode (show messages), false: silent mode
// redundancy: remove redundant (palindromic) matches
// truncate: size of truncated sequence (0: do not truncate)
//---------------------------------------------------------------------------

//(un)set non-occurance penalty mode
nop_mode=nop;

//(un)set verbose-mode (show messages):
verbose_mode=mode;

//(un)set log-mode (write log files):
log_mode=logm;

//(un)set remove-redundancy-mode (remove redundant palindromic matches)
redundancy_mode=redundancy;

//get sequence and pattern files:
sequence_file=sf; pattern_file=pf;
genome.open(sequence_file.data()); //open sequence file

//truncate (first) sequence (for test mode)
if ((genome.sequence[0].seq.size()>truncate)&&(truncate)) {
  genome.sequence[0].seq=genome.sequence[0].seq.substr(0,truncate);
}
openpattern(pattern_file.data(),0);  //open pattern file and apply on first sequence of sequence file

}

footprint::~footprint(){
}

string footprint::strval(int a) {
  char s[80];
  sprintf(s,"%i",a);
  return string(s);
}

string footprint::upper(string s) {
//---------------------------------------------------------------------------
// Convert a string to uppercase letters
// s: string
//---------------------------------------------------------------------------
for (int i=0; i<s.length(); i++) {
  s[i]=toupper(s[i]);
}
return s;
}


time_t footprint::showtime() {
//---------------------------------------------------------------------------
// Show timer
//---------------------------------------------------------------------------
time_t now = time(0);
cout << ctime(&now);
return now;
}

void footprint::verbose_cout(string cout_str) {
//---------------------------------------------------------------------------
// Print messages immediately in verbose mode
// cout_str: output string
//---------------------------------------------------------------------------
if (verbose_mode) {
  cout << cout_str;
  cout.flush(); //show immediately
} //if
}

string footprint::cons2regex(string cons, int mm) {
//---------------------------------------------------------------------------
// Convert consensus to regular expression
// cons=IUPAC consensus, mm=number of mismatches
//---------------------------------------------------------------------------
string re;
string temp;
string consstr;
if (mm) {
  for (unsigned int j=0; j<cons.size(); j++) {
    temp=cons;
    temp.replace(j,1,"N");
    //cout<<temp<<endl;
    if (consstr.size()) consstr+="|";    
    consstr+=temp;
  } //for j
} //if 
else consstr=cons;

  for (unsigned int i=0; i<consstr.size(); i++) {
    switch(consstr.at(i)) {
      case 65: re+="A"; //A->A
      break;
      case 66: re+="[CGT]"; //B->[CGT]
      break;
      case 67: re+="C"; //C->C
      break;
      case 68: re+="[AGT]"; //D->[AGT]
      break;
      case 71: re+="G"; //G->G
      break;
      case 72: re+="[ACT]"; //H->[ACT]
      break;
      case 75: re+="[GT]"; //K->[GT]
      break;
      case 77: re+="[AC]"; //M->[AC]
      break;
      case 78: re+="[AGCT]"; //N->[AGCT]
      break;
      case 82: re+="[AG]"; //R->[AG]
      break;
      case 83: re+="[GC]"; //S->[GC]
      break;
      case 84: re+="T"; //T->T
      break;
      case 86: re+="[ACG]"; //V->[ACG]
      break;
      case 87: re+="[AT]"; //W->[AT]
      break;
      case 89: re+="[CT]"; //Y->[CT]
      break;
      case 124: re+="|"; //|->|
      break;

    }; //switch
  } //for i 

return re;  
}

void footprint::openpattern(string patternfile, int seqno) {
//------------------------------------------------------------
// Open pattern file:
// patternfile: filename of XML pattern file
// seqno: sequence number of sequence file
//------------------------------------------------------------

  char buf[1000]; //file reading buffer
  patternXML="";  //delete XML raw data string

  //(1) open and read XML pattern file:
  ifstream InFile(patternfile.data());
    while (!InFile.eof()) {
      InFile.getline(buf,1000); //read stream line by line
      patternXML+=buf; //concatenate XML raw data string
    } //while
  InFile.close();

  //(2) parse XML pattern file:
  matrix tmp_matrix;        //template weight matrix
  pattern_type tmp_pattern; //template pattern
  spacer_type  tmp_spacer;  //template spacer
  string tag;        //contains one XML tag
  patternXML_pos=0;  //set XML data cursor to position 0 (start position)
  int counter;       //universl count variable for (pattern und matrix vectors)

  //pattern features from pattern file:
  string p_name;                                   //name of subpattern
  float  p_msens, p_csens, p_score;                //weigth matrix and core sensitivity, score
  unsigned char p_model, p_core_size, p_ori, p_mm; //weight matrix model, core size and subpattern orientation, mismatch number
  string p_string;                                 //IUPAC string
  int s_min, s_max;                                //minimal and maximal spacer distance

  patternXML=upper(patternXML); //convert to uppercase

  //parse XML string:
  while (patternXML_pos>=0) {
    tag=getnexttag(patternXML_pos);
    
    // 0) pattern features
    //se_threshold=-100.0; //default value
    if (tag.substr(0,7)=="PATTERN") {
      cum_threshold=atof(getattrib(tag.data(), "THRESHOLD").data());   //get cumulative threshold
      se_threshold=atof(getattrib(tag.data(), "SE_THRESHOLD").data()); //get stacking energy threshold
      if (se_threshold==0) se_threshold=-100.0; //default value
    }

    // I) pattern type: Weight Matrix:
    if (tag.substr(0,6)=="MATRIX") {
      mat.push_back(tmp_matrix);
	  counter=mat.size()-1;

	  //get pattern parameters from attributes:
      p_name      =     getattrib(tag.data(), "NAME" );               //weight matrix name
 	  p_model     =atoi(getattrib(tag.data(), "MODEL").data());       //weight matrix model
	  if (!p_model) p_model=5; //if no model assigened set default model
	  
      p_msens     =atof(getattrib(tag.data(), "MSENS").data());       //weight matrix (overall) sensitivity
	  p_csens     =atof(getattrib(tag.data(), "CSENS").data());       //core sensitivity
	  p_core_size =atoi(getattrib(tag.data(), "CORE" ).data());       //core size
	  p_score     =atof(getattrib(tag.data(), "THRESHOLD").data());       //score (optionally)
      p_ori       =atoi(getattrib(tag.data(), "ORIENTATION").data()); //subpattern orientation

	  //get sites comprising the weight matrix:
      while ((patternXML_pos>=0)&&(tag.substr(0,7)!="/MATRIX")) {
        tag=getnexttag(patternXML_pos);
        if (tag.substr(0,8)=="SEQUENCE") {
          mat[counter].addsite(getXMLdata(patternXML_pos),
		                       getattrib(tag.data(),"ID")); //add sites
        } //if
      } //while

	  //get sequence composition :

	  //long temp_array[4]={1,1,1,1};

      long temp_array[4]={genome.getSC(seqno,0),
                          genome.getSC(seqno,1),
                          genome.getSC(seqno,2),
                          genome.getSC(seqno,3)};

	  /*
	  for (unsigned int z=0; z<4; z++) {
	    cout << "seq " << z << " " << temp_array[z] << endl;
	  }
	  */

	  //position weight matrix initialization:
	  //add sequence composition to weight matrix model:
      mat[counter].addseqcomp(temp_array,4);
	  //calculate weight values and scores
	  mat[counter].init(p_name, p_model, p_msens, p_csens, p_core_size, p_score, nop_mode);

	  //add new subpattern into pattern
	  pattern.push_back(tmp_pattern); //add new subpattern into pattern
      counter=pattern.size()-1;
	  pattern[counter].type=1;
	  pattern[counter].no=mat.size()-1;    //references matrix index
	  pattern[counter].orientation=p_ori;  // orientation of subpattern

    } //if Matrix

    //(2) Consensus:
    if (tag.substr(0,9)=="CONSENSUS") {

      //get pattern parameters from attributes:
	  p_name  =getattrib(tag.data(), "NAME" );                     //consensus name
	  p_string=getattrib(tag.data(), "STRING" );                   //consensus name
      p_mm  =atoi(getattrib(tag.data(), "MISMATCH" ).data());      //mismatch number
	  p_ori =atoi(getattrib(tag.data(), "ORIENTATION" ).data());   //subpattern orientation

	  //consensus initialization:
      mat.push_back(tmp_matrix);
	  counter=mat.size()-1;
	  mat[counter].mname=p_name;
      mat[counter].consensus=p_string;
	  mat[counter].cons2matrix(p_string);
	  mat[counter].mism=p_mm;

	  //add new subpattern into pattern
	  pattern.push_back(tmp_pattern);
      counter=pattern.size()-1;
	  pattern[counter].type=2;
      pattern[counter].no=mat.size()-1;    //references matrix index
	  pattern[counter].orientation=p_ori;  // orientation of subpattern
	} // if Consensus

    //(3) Regular Expression:
    if (tag.substr(0,5)=="REGEX") {

	  //get pattern parameters from attributes:
	  p_name  =getattrib(tag.data(), "NAME" );                     //consensus name
	  p_string=getattrib(tag.data(), "REGEX" );                   //consensus name
	  p_ori =atoi(getattrib(tag.data(), "ORIENTATION" ).data());   //subpattern orientation

	  //regular expression initialization:
	  mat.push_back(tmp_matrix);
	  counter=mat.size()-1;
	  mat[counter].mname=p_name;
      mat[counter].regex=p_string;

      //add new subpattern into pattern
	  pattern.push_back(tmp_pattern);
      counter=pattern.size()-1;
	  pattern[counter].type=3;
	  pattern[counter].no=mat.size()-1;    //references matrix index
	  pattern[counter].orientation=p_ori;  // orientation of subpattern

    }  // if regular expression

    //(4) Spacer
    if (tag.substr(0,6)=="SPACER") {

	  //get pattern parameters from attributes:
      s_min=atoi(getattrib(tag.data(),"MIN").data()); //minimal spacer distance
 	  s_max=atoi(getattrib(tag.data(),"MAX").data()); //maximal spacer distance

	  //add new subpattern into pattern
	  spacer.push_back(tmp_spacer);    //create new subpattern
      counter=spacer.size()-1;          //get number of subpatterns
	  spacer[counter].min_space=s_min;  //set minimal
	  spacer[counter].max_space=s_max;  //and maximal spacer distance

    }  // if Spacer

  } //while

}

string footprint::getnexttag (int pos) {
//------------------------------------------------------------
// Return XML tag next to position pos
// note: <> are excluded, attributes are included
//------------------------------------------------------------
  int tag_pos1=patternXML.find("<", pos);
  int tag_pos2=patternXML.find(">", pos);
  patternXML_pos=tag_pos2+1;
  if (patternXML_pos==0) {
    patternXML_pos=-1;
    return "";
  } else return patternXML.substr(tag_pos1+1,tag_pos2-tag_pos1-1);
}

string footprint::getattrib (string tag, string att_name) {
//------------------------------------------------------------
// Return attribute
// XML tag: tag
// attribute name: att_name
//------------------------------------------------------------
  string temp_str;
  int att_pos =tag.find(att_name, 0);
  int att_pos1=tag.find("\"", att_pos);
  int att_pos2=tag.find("\"", att_pos1+1);
  if (att_pos) temp_str=tag.substr(att_pos1+1,att_pos2-att_pos1-1); else temp_str="";
  return temp_str;
}

string footprint::getXMLdata (int pos) {
//------------------------------------------------------------
// Return XML data between two tags
// pos: start of XML data
//------------------------------------------------------------
  int tag_pos=patternXML.find("<", pos);
  return patternXML.substr(pos,tag_pos-pos);
}

void footprint::showinfo() {
//------------------------------------------------------------
// Show information about pattern features:
//------------------------------------------------------------
  cout << "Virtual Footprint Pattern Information" << endl;
  cout << "-------------------------------------" << endl;
  cout << "Patternfile: " << pattern_file << endl;
  cout << "Sequencefile: " << sequence_file << endl;

  for (int m=0; m<mat.size(); m++) {
    mat[m].show();
  }

}

long footprint::analyze(int seq_nr) {
//------------------------------------------------------------
// Analyse subpatterns with spacers and filter out subpatterns
// that don't fulfill the criteria of the given pattern
// seq_nr: corresponding (genome-) sequence
// returns the number of pattern hits
//------------------------------------------------------------
series_type tmp_series;       //empty series to expand result result vector
vector<series_type> p_series; //result series vector on plus strand
vector<series_type> t_series; //temporary result series vector
vector<series_type> n_series; //result series vector on minus strand
vector<series_type> m_series; //merged series vector of plus and minus strand
int dist_space;
int rnum, tnum, ssize, s1, s2;
bool pattern_hit, stop_test;
int q=0;
vector<series_type>::iterator it;
int series_pos;


if (pattern.size()>1) { //pattern consists of multiple subpatterns?

//------------------------------------------------------------
//(1) Analyze + strand:
//------------------------------------------------------------

//write matches of first subpattern (indexes result) into result series (p_series):
//cout << endl;
for (unsigned int i=0; i<result.size(); i++) {
    
  if (result[i].pattern_no==0) {
    //if subpattern orientation restricted to one strand
	//only write matches in + orientation into series
    //if subpattern orientation is allowed on both strands
	//take matches in + and - orientation
    if (((!pattern[0].orientation) && (!result[i].orientation)) ||
	    (pattern[0].orientation)){
//	cout << result[i].pattern_no << " ! " << result[i].start_pos << " ! " << result[i].end_pos << " ! " << result[i].seq << endl;  
      p_series.push_back(tmp_series);
	  series_pos=p_series.size()-1;
	  p_series[series_pos].result.push_back(i);             //result index
	  p_series[series_pos].start_pos=result[i].start_pos;   //
	  p_series[series_pos].orientation=0;
	  
	} //if/elseif
  } //if (first subpattern)
} //for i

for (unsigned int p=1; p<pattern.size(); p++) { //pth pattern
  s1=spacer[p-1].min_space; //get min size of pth spacer
  s2=spacer[p-1].max_space; //get max size of pth spacer
  ssize=p_series.size();    //get size of pth result series

  for (unsigned int s=0; s<ssize; s++) { //expand result series
    pattern_hit=false;
    stop_test=false;
    rnum=p_series[s].result[p-1]; //last result number
    tnum=rnum+1;                  //next result number to test
    while ((tnum<result.size())&&(!stop_test)) {
      //check maximal subpattern distance
      if ((result[tnum].start_pos - result[rnum].end_pos) <= s2) {
        //check minimal subpattern distance 
        if ((result[tnum].start_pos - result[rnum].end_pos)>=s1) {
	  //check subpattern number  
	  if (result[tnum].pattern_no==p) { 
	    //check subpattern orientation   
	    if (((!result[tnum].orientation) && (!pattern[p].orientation)) ||
	        (pattern[p].orientation)) {
              if (!pattern_hit) { //expand series (first subpattern combination)
	        p_series[s].result.push_back(tnum); //expand result series
	        t_series.push_back(p_series[s]);    //store expanded result series tomporary
	        pattern_hit=true;                   //store that current result series was expanded
	      } //if
	      else { //add new series (if more subpattern combinations)
	        p_series.push_back(p_series[s]);                //copy current result as new series at the end of the vector
                p_series[p_series.size()-1].result[p]=tnum;     //replace last series
	        t_series.push_back(p_series[p_series.size()-1]);//store new result series tomporary
	      } //elseif
	    } //if subpattern orientation
	  } //if subpattern number   
	} //if minimal subpattern distance
	tnum++; 
      } //if maximal subpattern distance
      else stop_test=true; //stop pattern combination test with this series
      
    } //while

  } //for s (series number loop)

  //delete all series entries, with no combinstion of subpatterns:
  //note: erase is too slow -> replace p_series with t_series
  p_series=t_series;
  t_series.clear();  //clear temporary series vector

} //for p (pattern number loop)

//------------------------------------------------------------
//(2) Analyze - strand:
//------------------------------------------------------------

//write matches of first subpattern into result series:
for (int i=result.size()-1; i>=0; i--) {
  if (result[i].pattern_no==0) { //0=first subpattern
    //if subpattern orientation restricted to one strand
	//only write matches in - orientation into series
    //if subpattern orientation is allowed on both strands
	//take matches in + and - orientation
    if (((!pattern[0].orientation) && (result[i].orientation)) ||
	    (pattern[0].orientation)){
      n_series.push_back(tmp_series);
	  series_pos=n_series.size()-1;
	  n_series[series_pos].result.push_back(i);       //result index
      n_series[series_pos].end_pos=result[i].end_pos;
	  n_series[series_pos].orientation=1;

	} //if/elseif
  } //if (first subpattern)
} //for i


for (unsigned int p=1; p<pattern.size(); p++) { //pth pattern
  s1=spacer[p-1].min_space;  //get min size of pth spacer
  s2=spacer[p-1].max_space;  //get max size of pth spacer
  ssize=n_series.size();     //get size of pth result series

  for (unsigned int s=0; s<ssize; s++) { //sth  series result
    pattern_hit=false;
    stop_test=false;
    rnum=n_series[s].result[p-1]; //last result number
    tnum=rnum-1;                  //next result number to test
    //cout << s << "(" << tnum << ")" << "|"; cout.flush();
    while ((tnum>=0)&&(!stop_test)) { //test pattern combinations 
      //check maximal subpattern distance
      if ((result[rnum].start_pos - result[tnum].end_pos) <= s2) {
        //check minimal subpattern distance 
        if ((result[rnum].start_pos - result[tnum].end_pos)>=s1) {
	  //check subpattern number  
	  if (result[tnum].pattern_no==p) {
       	    //check subpattern orientation
	    if (((result[tnum].orientation) && (!pattern[p].orientation)) ||
	        (pattern[p].orientation)) {
              if (!pattern_hit) { //expand series (first subpattern combination)
                n_series[s].result.push_back(tnum); //expand existing result series
	        t_series.push_back(n_series[s]);    //store expanded result series tomporary
	        pattern_hit=true;                   //store that current result series was expanded
	      } //if	    
	      else { //add new series (if more subpattern combinations)
	        n_series.push_back(n_series[s]);                //copy current result as new series at the end of the vector
                n_series[n_series.size()-1].result[p]=tnum;     //replace last series
                t_series.push_back(n_series[n_series.size()-1]);//store new result series tomporary
              } //elseif
	    } //if subpattern orientation
	  } //if subpattern number
	} //if minimal subpattern distance
	tnum--;
      } //if maximal subpattern distance
      else stop_test=true; //stop pattern combination test with this series
      
    } //while
  } //for s (series number loop)

  //delete all series entries, with no combinstion of subpatterns:
  //note: erase is too slow -> replace n_series with t_series
  n_series=t_series;
  t_series.clear();  //clear temporary series vector

} //for p (pattern number loop)


//supplement start position, end position, overall-score and overall-sequence of series:
for (unsigned int i=0; i<p_series.size(); i++) {
  p_series[i].end_pos=result[p_series[i].result[p_series[i].result.size()-1]].end_pos;
  p_series[i].score=0;
  for (unsigned int j=0; j<p_series[i].result.size(); j++) {
    p_series[i].score+=result[p_series[i].result[j]].score;
  }
  p_series[i].seq=genome.subseq(seq_nr, p_series[i].start_pos-1, p_series[i].end_pos-p_series[i].start_pos+1); 
  if (p_series[i].score>=cum_threshold) series.push_back(p_series[i]);
  //series.push_back(p_series[i]);
}
for (unsigned int i=0; i<n_series.size(); i++) {
  n_series[i].start_pos=result[n_series[i].result[n_series[i].result.size()-1]].start_pos;
  n_series[i].score=0;
  for (unsigned int j=0; j<n_series[i].result.size(); j++) {
    n_series[i].score+=result[n_series[i].result[j]].score;
  }
  n_series[i].seq=genome.revcomp(seq_nr, n_series[i].start_pos-1, n_series[i].end_pos-n_series[i].start_pos+1); 
  if (n_series[i].score>=cum_threshold) series.push_back(n_series[i]);
  //series.push_back(n_series[i]);
}




//merge and order dataset:
sort (series.begin(), series.end(), series_order());

} //if pattern.size()

else { //single pattern -> write as series

for (unsigned int i=0; i<result.size(); i++) {
  series.push_back(tmp_series);
  series_pos=series.size()-1;
  series[series_pos].result.push_back(i);               //result index
  series[series_pos].start_pos=result[i].start_pos;     //
  series[series_pos].end_pos=result[i].end_pos;         //
  series[series_pos].orientation=result[i].orientation; //
  series[series_pos].seq=result[i].seq; //
  series[series_pos].score=result[i].score;             //
  
} //for

} //else

return series.size(); //return number of pattern hits

}

void footprint::remove_lower_scores() {
//------------------------------------------------------------
// Keep only highest scoring composite matches (series)
// with same start position:
// Note: series must be sorted by start position
//------------------------------------------------------------
vector<series_type> tmp_series; //temporary series vector
int old_i=-1;

for (unsigned int i=0; i<series.size(); i++) {
  if (old_i==-1) {
    old_i=i;
  }  
  else if (series[i].start_pos==series[old_i].start_pos) { //same start positions
    //keep match with highest score
    if (series[i].score>=series[old_i].score) {
      old_i=i; 
    }
  } // if
  else {
    tmp_series.push_back(series[old_i]);
    old_i=-1;
  }
} //for  i
 
series=tmp_series;
cout << "Result: " << series.size();
}
  
 

void footprint::remove_redundancies() {
//------------------------------------------------------------
// Remove redundant palindromic (subpattern-) matches:
// Note: keeps the match with the highest score
//------------------------------------------------------------
vector<result_type> tmp_result; //temporary result vector
int old_i;                      //resultnumber to compare

if (result.size()>1) { //must have more than one result
  for (unsigned int p=0; p<pattern.size(); p++) { //handle all subpattern matches
    old_i=-1; //initialize with "unknown" result number
    for (int i=1; i<result.size(); i++) { //results loop
      if (result[i].pattern_no==p) { //check if match belongs to right subpattern
        if (old_i==-1) {
	  old_i=i; //initialize with result number
	}
	else {
          if (result[i].start_pos==result[old_i].start_pos) { //same start positions
	    //keep match with highest score
            if (result[i].score>=result[old_i].score) { 
	      tmp_result.push_back(result[i]);
	    } else tmp_result.push_back(result[old_i]);
	    old_i=-1; //(re-)initialize with "unknown" result number
	  }
	  else { //different matches
	    tmp_result.push_back(result[old_i]); 
	    old_i=i;
	  }
        }
      } //if
    }  //for i
    if (old_i!=-1) tmp_result.push_back(result[old_i]); //last result
  } //for p  
  result=tmp_result;  //replace result vector with temporary result vector
  tmp_result.clear(); //clear temporary result vector
  sort (result.begin(), result.end(), result_order()); //sort again
} //if resultsize>1

}

void footprint::execute(string resultfile, unsigned char file_format) {
//------------------------------------------------------------
// Execute pattern search:
// resultfile    : filename and path of resultfile
// file_format   : 0: XML, 1: CSV
// scoring method: 0 : individual information score
//                 1 : logo score
//                 2 : MatInspector score
//                 3 : MatrixSearch score
//                 4 : TargetExplorer score
//                 5 : biased logo score
//                 6 : to do: Berg and von Hippel score
//------------------------------------------------------------
  int mat_counter=0;
  long hit_nr;                             //number of overall hits
  long pat_nr;                             //number of hits that fulfill the pattern criteria
  time_t T1, T2;                           //time points (start and end of virtual footprint analysis
  vector<series_type> tmp_series;          //temporary series vector

  if (verbose_mode) {
      cout << endl << "start pattern search ... "; T1=showtime(); 
  }

  for  (unsigned int s=0; s<genome.getnumber(); s++) {  //sequences (later for more sequences)
  //for  (unsigned int s=0; s<1; s++) {                   //sequences (currently restriced to one sequence)

    if (verbose_mode) {
	  cout << endl << "Sequence: " << genome.getID(s) << endl;
    }

	if (s) { //reload pattern if new sequence
	  openpattern(pattern_file, s);
	}

    for (unsigned int p=0; p<pattern.size(); p++) {     //sub-patterns (weight matrices, IUPAC, regular expressions)

	  switch (pattern[p].type) {
	    case 1: //weight matrix
		  if (verbose_mode) {
		    cout << "starting weight matrix search with " << mat[mat_counter].mname << "..." << endl; cout.flush();
		    cout << "score sensitivity: " << mat[mat_counter].msens
		         << " set threshold score: " << mat[mat_counter].mth << endl;
		    if (mat[mat_counter].core_size) { 
		    cout << "core  sensitivity: " << mat[mat_counter].csens
		         << " set threshold score: " << mat[mat_counter].cth << endl;
		    }
          }
          hit_nr=msearch(mat_counter, s, p, mat[mat_counter].mth, mat[mat_counter].cth, mat[mat_counter].core_size);
		  mat_counter++;
		  if (verbose_mode) cout<< "number of hits: " <<hit_nr<<endl;

		break;

		case 2: //IUPAC code
		  if (verbose_mode) {
		    cout << "starting consensus search with " << mat[mat_counter].mname << "..." << endl;
            cout << "number of allowed mismatches: "<< mat[mat_counter].mism << endl;
          }
          hit_nr=csearch(mat_counter, s, p, mat[mat_counter].mism);
		  mat_counter++;
		  if (verbose_mode) cout<< "number of hits: " <<hit_nr<<endl;

		break;

		case 3: //regular expression
		  if (verbose_mode) {
		    cout << "starting regular expression search with " << mat[mat_counter].mname << "..." << endl;
          }
          hit_nr=regexsearch(mat_counter, s, p);
		  mat_counter++;
		  if (verbose_mode) cout<< "number of hits: " <<hit_nr<<endl;
		break;

		case 4: //spacer
		break;
	  } //switch
    } //for p (subpatterns)

    verbose_cout("sorting results...");
    sort (result.begin(), result.end(), result_order());
    verbose_cout("done\n");

    if (redundancy_mode) { 
      verbose_cout("remove redundancies...");
      remove_redundancies();
      
	verbose_cout("done\n");
}	

	verbose_cout("analyzing results...");
    pat_nr=analyze(s);
	verbose_cout("done\n");
    verbose_cout("number of pattern hits: ");
	verbose_cout(strval(pat_nr)+"\n");
	
	
	
	verbose_cout("restrict to high-scoring composites...");
	//remove_lower_scores();
	verbose_cout("done\n");
    

	verbose_cout("calculating stacking energy...");
	cout << "threshold: " << se_threshold << endl;
	for (int i=0; i<series.size(); i++) {
	//cout << i << "!"; cout.flush();
      series[i].se_score=genome.stacking_energy(s, series[i].start_pos-1, series[i].end_pos-1, series[i].orientation, 80, 120, 31);
          if (se_threshold<=series[i].se_score) {
	    tmp_series.push_back(series[i]);
	  }
	}
	series=tmp_series;
	
	//========================================================================
	sort (tmp_series.begin(), tmp_series.end(), series_order_score());

	ofstream output;   //file open stream
	//open result file (if more than 1 sequence is examinded, open for append):
	output.open("/var/lib/nobody/tmp/size.log");
	for (unsigned int i=0; i<series.size(); i++) {
  	output << tmp_series[i].seq.size() << "\t" << tmp_series[i].score << "\t" << tmp_series[i].se_score << endl;
	}  
	output.close(); //close file
	//========================================================================
	
	
	
	verbose_cout("done\n");

	verbose_cout("saving results...");
	saveresults(resultfile, s, file_format); //save result in file
	verbose_cout("done\n");

	series.clear(); //clear series
	tmp_series.clear(); //clear series
	result.clear(); //clear results
	mat.clear();    //cleat weight matrices
	spacer.clear(); //clear spacers
	pattern.clear();//clear pattern
	mat_counter=0;

  } //for s (sequences)

  if (verbose_mode) {
    cout << endl << "finished pattern search "; T2=showtime();
	cout <<"length of time: " << difftime(T2, T1) << " seconds" << endl;
	cout <<"thank you for using Virtual Footprint!" << endl;
  }

}

string footprint::showCSVmatch(int listnr, string sep) {
//------------------------------------------------------------
// Return result line as CSV (comma separated values)
// listnr: position in resultlist
// sep:  list separator
//------------------------------------------------------------
string ori;            //orientation
stringstream CSV_line; //concatenates one line as CSV

  //replace 0/1 orientation with +/-:
  if (result[listnr].orientation==0) ori="+"; else ori="-";

  //concatenate by use of a stringstream:
  CSV_line << result[listnr].pattern_no  << sep
           << result[listnr].no          << sep
           << result[listnr].name        << sep
           << result[listnr].type        << sep
           << result[listnr].start_pos   << sep
           << result[listnr].end_pos     << sep
           << ori                        << sep
		   << result[listnr].seq         << sep
           << result[listnr].score       << sep
           << result[listnr].core_score;

return CSV_line.str();
}

string footprint::showXMLmatch(int listnr) {
//------------------------------------------------------------
// Return result line as XML (eXtensible Markup Language)
// listnr: position in resultlist
// sep:  list separator
//------------------------------------------------------------
string ori;            //orientation
stringstream XML_line; //concatenates one line as XML

  //replace 0/1 orientation with +/-:
  if (result[listnr].orientation==0) ori="+"; else ori="-";

  //concatenate match as XML by use of a stringstream:
  XML_line << "<Site "
   	       << "name=\""          << result[listnr].name        << "\" "
     	   << "type=\""          << result[listnr].type        << "\" "
           << "startPosition=\"" << result[listnr].start_pos   << "\" "
           << "endPosition=\""   << result[listnr].end_pos     << "\" "
           << "orientation=\""   << ori                   << "\" "
           << "score=\""         << result[listnr].score       << "\" "
           << "core_score=\""    << result[listnr].core_score  << "\" "
           << ">"
           << result[listnr].seq
           << "</Site>";

return XML_line.str();
}


void footprint::saveresults(string resultfile, int seq_nr, unsigned char format) {
//---------------------------------------------------------------------------
// Save results of pattern search:
// resultfile: name of resultfile
// seq_nr: number of examined sequence (genome FASTA file)
// format: file format of resultfile 0: XML (default), 1: CSV
//---------------------------------------------------------------------------
  string ori;        //orientation
  ofstream output;   //file open stream

  //open result file (if more than 1 sequence is examinded, open for append):
  if (!seq_nr) output.open(resultfile.data());
      else output.open(resultfile.data(), ios::app);

  switch (format) {

    case 0: //XML format

      if (!seq_nr) { //print header and VfpResult Tag at the top of the file (if first sequence)
	    output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
	    output << "<VfpResult>" << endl;
	  }

      output << "<Sequence "
	         << "ID=\"" <<genome.getID(seq_nr) << "\" matches=\"" << series.size() << "\">" << endl;

      for (unsigned int i=0; i<series.size(); i++) {

        //replace 0/1 orientation with +/-:
        if (series[i].orientation==0) ori="+"; else ori="-";
        output << "<Match "
               << "number=\"" << i << "\" "
		       << "startPosition=\"" << series[i].start_pos   << "\" "
               << "endPosition=\""   << series[i].end_pos     << "\" "
               << "orientation=\""   << ori                   << "\" "
	       << "sequence=\""   << series[i].seq         << "\" "
			   << "score=\""         << series[i].score       << "\" "
			   << "se_score=\""      << series[i].se_score    << "\">"
			   << endl;

        for (unsigned int j=0; j<series[i].result.size(); j++) {
          output << showXMLmatch(series[i].result[j]) << endl;
        } //for j

        output << "</Match>" << endl;
      }//for i

      output << "</Sequence>" << endl;

	  //if last sequence write close tag:
	  if (seq_nr==(genome.getnumber()-1)) output <<"</VfpResult>" << endl;

    break;

	case 1: //CSV format
    /*
    for (unsigned int i=0; i<result.size(); i++) {
       output << showCSV(i,",")<<endl;
    } //for
    */
    int pattern_size=pattern.size();
    for (unsigned int i=0; i<series.size(); i++) {
	  if (pattern_size>1) { //print header for composite patterns
        //replace 0/1 orientation with +/-:
        if (series[i].orientation==0) ori="+"; else ori="-";
	    output << "pattern match," << series[i].start_pos << "," << series[i].end_pos << "," << ori << endl;
	  }
      for (unsigned int j=0; j<series[i].result.size(); j++) {
        output << showCSVmatch(series[i].result[j],",")<<endl;
      } //for j

    } //for i

  } //switch

  output.close(); //close file
}

long footprint::msearch(int mnr, int gnr, int pnr, float Smin, float Cmin, int core) {
//---------------------------------------------------------------------------
// Start weight matrix search:
// mnr: matrix number
// gnr: genome number
// Smin/Smax: minimal and maximal score
//---------------------------------------------------------------------------
  result_type hit;
  fasta tmp;

  float C, S;        //core score and weight matrix score
  long site_nr=0;    //counter for number of hits
  const char *seq;

  seq=genome.sequence[gnr].seq.data();

  for (unsigned int i=0; i<(genome.getseq(gnr).size()-mat[mnr].msize+1); i++) {

    if (core) {
      //plus direction:
	C=mat[mnr].getCoreScore(seq,i,0);
	//if (C>9.8) cout << C << "/";
	if (C>=Cmin) {
        S=mat[mnr].getScore(seq,i,0);
	
        if (S>=Smin) {
		  hit.no=mnr;
		  hit.pattern_no=pnr;
		  hit.name=mat[mnr].mname;
		  hit.type=1;
          hit.start_pos=i+1;
          hit.end_pos=i+mat[mnr].msize;
          hit.orientation=0;
          hit.seq=genome.subseq(gnr,i,mat[mnr].msize);
          hit.score=S;
		  hit.core_score=C;
          result.push_back(hit);
          site_nr++;
        } //if (R>=Smin)
      } //if (C>=Cmin)

      //reverse complement direction:
      C=mat[mnr].getCoreScore(seq,i,1);
	  if (C>=Cmin) {
        S=mat[mnr].getScore(seq,i,1);
        if (S>=Smin) {
		  hit.no=mnr;
		  hit.pattern_no=pnr;
		  hit.name=mat[mnr].mname;
		  hit.type=1;
          hit.start_pos=i+1;
          hit.end_pos=i+mat[mnr].msize;
          hit.orientation=1;
          tmp.set(0,"",genome.subseq(gnr,i,mat[mnr].msize));
          hit.seq=tmp.revcomp(0,0,0);
          hit.score=S;
		  hit.core_score=C;
          result.push_back(hit);
          site_nr++;
        } //if (R>=Smin)
      } //if (C>=Cmin)
    } //if core

    else { // no core
      //plus direction:
      S=mat[mnr].getScore(seq,i,0);
      if (S>=Smin) {
	    hit.no=mnr;
		hit.pattern_no=pnr;
	    hit.name=mat[mnr].mname;
		hit.type=1;
        hit.start_pos=i+1;
        hit.end_pos=i+mat[mnr].msize;
        hit.orientation=0;
        hit.seq=genome.subseq(gnr,i,mat[mnr].msize);
        hit.score=S;
		hit.core_score=C;
        result.push_back(hit);
        site_nr++;
      } //if (R>=Smin)

      //reverse complement direction:
      S=mat[mnr].getScore(seq,i,1);
      if (S>=Smin) {
	    hit.no=mnr;
		hit.pattern_no=pnr;
	    hit.name=mat[mnr].mname;
		hit.type=1;
        hit.start_pos=i+1;
        hit.end_pos=i+mat[mnr].msize;
        hit.orientation=1;
        tmp.set(0,"",genome.subseq(gnr,i,mat[mnr].msize));
        hit.seq=tmp.revcomp(0,0,0);
        hit.score=S;
		hit.core_score=C;
        result.push_back(hit);
        site_nr++;
      } //if (R>=Smin)
    } //else (no core)
  } //for i

  return site_nr; // return number of hits
}


long footprint::csearch(int mnr, int gnr, int pnr, int mm) {
//---------------------------------------------------------------------------
// Start consensus matrix search:
// mnr: consensus number
// mm: number of allowed mismatches
//---------------------------------------------------------------------------

  result_type hit;
  fasta tmp;
  float S;           //mismatch score

  long site_nr=0;    //counter for number of hits
  const char *seq;

  seq=genome.sequence[gnr].seq.data();

  for (unsigned int i=0; i<(genome.getseq(gnr).size()-mat[mnr].msize+1); i++) {

    //plus direction:
    S=mat[mnr].getMismatch(seq,i,0,mm);
    if (S<=mm) {
	  hit.no=mnr;
	  hit.pattern_no=pnr;
	  hit.name=mat[mnr].mname;
	  hit.type=2;
      hit.start_pos=i+1;
      hit.end_pos=i+mat[mnr].msize;
      hit.orientation=0;
      hit.seq=genome.subseq(gnr,i,mat[mnr].msize);
	  //cout << i << ":" << hit.seq <<endl;
      hit.score=S;
	  hit.core_score=0;
      result.push_back(hit);
      site_nr++;
    } //if (R>=Smin)

      //reverse complement direction:
    S=mat[mnr].getMismatch(seq,i,0,mm);
    if (S<=mm) {
	  hit.no=mnr;
	  hit.pattern_no=pnr;
	  hit.name=mat[mnr].mname;
	  hit.type=2;
      hit.start_pos=i+1;
      hit.end_pos=i+mat[mnr].msize;
      hit.orientation=1;
      tmp.set(0,"",genome.subseq(gnr,i,mat[mnr].msize));
      hit.seq=tmp.revcomp(0,0,0);
	  //cout << i << ":" << hit.seq <<endl;
      hit.score=S;
	  hit.core_score=0;
      result.push_back(hit);
      site_nr++;
    } //if (R>=Smin)

  } //for i

return site_nr; // return number of hits

}


long footprint::regexsearch(int mnr, int gnr, int pnr) {
//---------------------------------------------------------------------------
// Start regular expression search:
//---------------------------------------------------------------------------
regex_t re;                   //regular expression object
regmatch_t regex_hit;         //regular expression matching object

int regex_nr=0;               //counter for number of hits
int start=0;                  //start at position 0
long match_pos, match_length; //match position and langth

string rc_sequence;           //reverse complement sequence
const char *seq;              //pointer to sequence

string restr=mat[mnr].regex;                     //get regular expression
seq=genome.sequence[gnr].seq.data();             //get pointer to sequence
long seq_length=genome.sequence[gnr].seq.size(); //sequence length
result_type hit;

regcomp(&re, restr.data(), REG_EXTENDED);        //initialize regular expression object



//plus direction
while (regexec(&re, seq+start, 1, &regex_hit, 0) == 0) {    // while matches found
  // substring found between hit.rm_so and hit.rm_eo
  match_pos=start+regex_hit.rm_so;        //match position
  match_length=regex_hit.rm_eo-regex_hit.rm_so; //match length
  start+=regex_hit.rm_eo;                 //new start position
  //cout <<"Position: "<<match_pos+1<<"-"<<match_pos+match_length<<"; "<<genome.subseq(gnr,match_pos,match_length)<<endl;
  hit.no=mnr;
  hit.pattern_no=pnr;
  hit.name=mat[mnr].mname;
  hit.type=3;
  hit.start_pos=match_pos+1;
  hit.end_pos=match_pos+match_length;
  hit.orientation=0;
  hit.seq=genome.subseq(gnr,match_pos,match_length);
  hit.score=0;
  hit.core_score=0;
  result.push_back(hit);
  regex_nr++; //increase number of hits
} //while

//get reverse complement
rc_sequence=genome.sequence[gnr].seq;                  //get sequence
reverse(rc_sequence.begin(), rc_sequence.end());       //reverse sequence

replace(rc_sequence.begin(), rc_sequence.end(),65,90); //A->Z
replace(rc_sequence.begin(), rc_sequence.end(),84,65); //T->A
replace(rc_sequence.begin(), rc_sequence.end(),90,84); //Z->T
replace(rc_sequence.begin(), rc_sequence.end(),67,90); //C->Z
replace(rc_sequence.begin(), rc_sequence.end(),71,67); //G->C
replace(rc_sequence.begin(), rc_sequence.end(),90,71); //Z->G

seq=rc_sequence.data(); //get pointer to sequence
start=0;                //start at first position

//minus direction
while (regexec(&re, seq+start, 1, &regex_hit, 0) == 0) {    // while matches found
  // substring found between pm.rm_so and pm.rm_eo
  match_pos=start+regex_hit.rm_so;        //match position
  match_length=regex_hit.rm_eo-regex_hit.rm_so; //match length
  start+=regex_hit.rm_eo;                 //new start position
  //cout <<"Position: "<<seq_length-match_pos-match_length+1<<"-"<<seq_length-match_pos<<"; "<<rc_sequence.substr(match_pos,match_length)<<endl;
  hit.no=mnr;
  hit.pattern_no=pnr;
  hit.name=mat[mnr].mname;
  hit.type=3;
  hit.start_pos=seq_length-match_pos-match_length+1;
  hit.end_pos=seq_length-match_pos;
  hit.orientation=1;
  hit.seq=rc_sequence.substr(match_pos,match_length);
  hit.score=0;
  hit.core_score=0;
  result.push_back(hit);
  regex_nr++; //increase number of hits
} //while

regfree(&re); //deallocate memory

return regex_nr; //return number of hits

}

