//---------------------------------------------------------------------------
// matrix.cpp
//---------------------------------------------------------------------------

#include "matrix.h"

matrix::matrix(){
  buflen=100000; //set file reading buffer
  // default genomic base frequency:
  Pb[0]=0.25; Pb[1]=0.25; Pb[2]=0.25; Pb[3]=0.25;
  maxsites=1000; //maximum number of sites
  
  //assigns ascii value of sequence letter to weight matrix array
  for (unsigned int i=0; i<90; i++) {
    charval[i]=4;
   rcharval[i]=4;
  }
  charval[65]=0; //A
 rcharval[65]=3;
  charval[67]=1; //C
 rcharval[67]=2;
  charval[71]=2; //G
 rcharval[71]=1;
  charval[84]=3; //T
 rcharval[84]=0;
}

matrix::~matrix(){
}


float matrix::ilog(float n, float a){
//---------------------------------------------------------------------------
// Logarithm of base 2 of n; define log2(0)=log2(1/(a+2))  
//---------------------------------------------------------------------------
if (n) return log(n)/log((float)2);
  else return log(1/((float)a+2))/log((float)2);
}

float matrix::elog(float n){           
//---------------------------------------------------------------------------
// Natural logarithm (base e) of n; define ln(0)=0
//---------------------------------------------------------------------------
if (n) return log(n);
  else return 0;
}


float matrix::round(float n, int d) { 
//---------------------------------------------------------------------------
// Round float number
// n: float
// d: number of digits
//---------------------------------------------------------------------------
  n *= powf(10, (float)d);
  if (n >= 0) 
      n=floorf(n + 0.5);
  else 
      n=ceilf(n - 0.5); 
  n /= powf(10, (float)d); 
  
  return n; 
}

void matrix::addseqcomp(long sc[],int num) {
//---------------------------------------------------------------------------
// Add sequence composition (optionally)
// if not, the composition is 0.25 for DNA {A,C,G,T}
//---------------------------------------------------------------------------
  long ChrSum=0; //character sum

  //(1) calculate character sum
  for (unsigned int i=0; i<num; i++) {
    ChrSum+=sc[i];
  } //for

  //(2) calculate character frequency
  for (unsigned int i=0; i<num; i++) {
    Pb[i]=(float)sc[i]/(float)ChrSum;
  } //for
}

void matrix::init(string i_name, unsigned char i_model, float i_msens,
                  float i_csens, unsigned char i_core_size, float i_score, int i_nop) {
//---------------------------------------------------------------------------
// Calculates weight matrix values (numbers, sums, frequencies and scores)
// overall information (Rsequence)
//---------------------------------------------------------------------------
//to do: replace i->b, s->l to get consistent s(b,l) (score of base b at position l)
//       replace maxcolscore

float maxcolscore=-999999.9;  //maximal score for a certain column
float mincolscore=+999999.9;  //minimal score for a certain column
unsigned char mat_row;        //weight matrix row number
string tmp_site;

//set global matrix variables:
mname=i_name;                 //weight matrix name
msize=site.getseq(0).size();  //weight matrix size
model=i_model;                //define weight matrix scoring model
msens=i_msens;                //weight matrix (overall) sensitivity
csens=i_csens;                //core sensitivity
core_size=i_core_size;        //core size

Rsequence=0.0;
Rsequence_bias=0.0;
/*
for (unsigned int i=0; i<site.getnumber(); i++) {
  cout << i << "   " << site.uppercase(i) << endl;
}//for
*/

for (unsigned int l=0; l<msize; l++) { //go by column

  //Reset all weight matrix values to zero
  for (unsigned int i=0; i<4; i++) {
    value[l][i].number=0;
	value[l][i].number_bias=0;
    colsum[l]=0;
	colsum_bias[l]=0;
  }

  //(1) calculate alignment numbers and sums:
  for (unsigned int i=0; i<site.getnumber(); i++) {
    mat_row=charval[site.uppercase(i).at(l)];
    value[l][mat_row].number++;
    if (mat_row<4) colsum[l]++; //increase column sum if valid letter
  } //for i

  //(2) calculate frequencies,
  //    overall information content (Rsequence),
  //    information content at position l (Rsequence(l))
  Rs[l]=0.0;
  for (unsigned int i=0; i<4; i++) {
    //frequencies:
    value[l][i].freq=(float)value[l][i].number/(float)colsum[l];
    //Rsequence(l):
    Rs[l]+=-ilog(value[l][i].freq, colsum[l])*value[l][i].freq;
  } //for i

  Rs[l]=2-Rs[l];    //Rsequence(l)

  Rsequence+=Rs[l]; //Rsequence (sum up Rsequence(l))

} //for l (numbers, frequencies, Rsequence, Rsequence(l)




/*
  for (unsigned int i=0; i<site.getseq(0).size(); i++) {
    cout << setw(6) << i+1;
    for (unsigned int n=0; n<4; n++) {
      cout << setw(6)  << round(value[i][n].number,2);
      //cout << setw(6)  << round(value[i][n].number_bias,2);
	  //cout << setw(6)  << round(value[i][n].freq_bias,2);
    } // for n
   cout << setw(6) << round(Rs[i],1);
	//cout << setw(6) << round(Rs[i],1);
    cout << endl;
  } //for i


cout << "Rsequence :" << Rsequence << endl;
*/




//Calculate matrix scores dependent of the underlying scoring model

float score_number=0.01; //for Patser/Target Explorer
maxscore=0; minscore=0;

switch (model) {

  case 0: //individual information scoring model
    for (unsigned int l=0; l<msize; l++) { //go by column
	  maxcolscore=-999999.9; mincolscore=+999999.9; //reset max/min
      for (unsigned int i=0; i<4; i++) {
	    //individual information weights (Riw(b,l)),
        value[l][i].score=2+ilog(value[l][i].freq, value[l][i].number);
        if (value[l][i].score>maxcolscore) maxcolscore=value[l][i].score;
		if (value[l][i].score<mincolscore) mincolscore=value[l][i].score;
      } //for i
	  maxscore+=maxcolscore; //sum up maximal score
	  minscore+=mincolscore; //sum up minimal score
    } //for l
  break;

  case 1: //sequence logo scoring model
    for (unsigned int l=0; l<msize; l++) { //go by column
      maxcolscore=-999999.9; mincolscore=+999999.9; //reset max/min
	  for (unsigned int i=0; i<4; i++) {
        //sequence logo weights (equivalent to the heigth of a certain letter in a sequence logo)
        value[l][i].score=value[l][i].freq*Rs[l];
        if (value[l][i].score>maxcolscore) maxcolscore=value[l][i].score;
		if (value[l][i].score<mincolscore) mincolscore=value[l][i].score;
      } //for i
	  maxscore+=maxcolscore; //sum up maximal score
	  minscore+=mincolscore; //sum up minimal score
    } //for l
  break;

  case 2: //MatInspector/Match scoring model
    for (unsigned int l=0; l<msize; l++) { //go by column
      Ci[l]=0.0; //calculate ci vector at position s
      for (unsigned int i=0; i<4; i++) {
        Ci[l]+=value[l][i].freq*elog(value[l][i].freq);
      } //for i
      Ci[l]=(100/elog(4))*(Ci[l]+elog(4)); //normalize ci vector

	  maxcolscore=-999999.9; mincolscore=+999999.9; //reset max/min
      for (unsigned int i=0; i<4; i++) {
	    //MatInspector score weight (equivalent to a normalized sequence logo score)
        value[l][i].score=value[l][i].number*Ci[l];
        if (value[l][i].score>maxcolscore) maxcolscore=value[l][i].score;
		if (value[l][i].score<mincolscore) mincolscore=value[l][i].score;
      } //for i
	  maxscore+=maxcolscore; //sum up maximal score
	  minscore+=mincolscore; //sum up minimal score
    } //for l
  break;

  case 3: //Matrix Search scoring model

    for (unsigned int l=0; l<msize; l++) { //go by column
	  maxcolscore=-999999.9; mincolscore=+999999.9; //reset max/min
      for (unsigned int i=0; i<4; i++) {
        value[l][i].score=(value[l][i].number+score_number)/((colsum[l]+score_number)*Pb[i]);
        if (value[l][i].score>maxcolscore) maxcolscore=value[l][i].score;
        if (value[l][i].score<mincolscore) mincolscore=value[l][i].score;
      } //for i
	  maxscore+=maxcolscore; //sum up maximal score
	  minscore+=mincolscore; //sum up minimal score
    } //for l
  break;

  case 4: //Patser/Target Explorer score
    for (unsigned int l=0; l<msize; l++) { //go by column
	  maxcolscore=-999999.9; mincolscore=+999999.9; //reset max/min
      for (unsigned int i=0; i<4; i++) {
        value[l][i].score=elog(((value[l][i].number+Pb[i])/(colsum[l]+1))/Pb[i]);
        if (value[l][i].score>maxcolscore) maxcolscore=value[l][i].score;
		if (value[l][i].score<mincolscore) mincolscore=value[l][i].score;
      } //for i
	  maxscore+=maxcolscore; //sum up maximal score
	  minscore+=mincolscore; //sum up minimal score
    } //for l
  break;

  case 5: //Sequence Logo scoring model with bias correction

    for (unsigned int l=0; l<msize; l++) { //go by column
      for (unsigned int i=0; i<4; i++) {
        value[l][i].number_bias=value[l][i].freq*((float)0.25/Pb[i])*(float)colsum[l];
        colsum_bias[l]+=value[l][i].number_bias;
      } //for i
	  Rs_bias[l]=0.0;
      for (unsigned int i=0; i<4; i++) {
        value[l][i].freq_bias=(float)value[l][i].number_bias/(float)colsum_bias[l];
        //Rsequence_bias(l)
        Rs_bias[l]+=-ilog(value[l][i].freq_bias, colsum_bias[l])*value[l][i].freq_bias;
        //value[l][i].score=value[l][i].freq_bias*Rs_bias[l]; //?????????????????
      } //for i

      Rs_bias[l]=2-Rs_bias[l];
      Rsequence_bias+=Rs_bias[l];

      maxcolscore=-999999.9; mincolscore=+999999.9; //reset max/min
      for (unsigned int i=0; i<4; i++) {
        value[l][i].score=value[l][i].freq_bias*Rs_bias[l];
		if ((!value[l][i].freq_bias)&&(i_nop)) {
		  value[l][i].score=float((1-colsum_bias[l])/(colsum_bias[l]+2))*Rs_bias[l];
		}
        if (value[l][i].score>maxcolscore) maxcolscore=value[l][i].score;
		if (value[l][i].score<mincolscore) mincolscore=value[l][i].score;
      }// for i
	  maxscore+=maxcolscore; //sum up maximal score
	  //cout << "maxi " << maxcolscore << " sum " << maxscore << endl;
	  minscore+=mincolscore; //sum up minimal score
    } //for l
  break;

} //switch (model)

getCore(core_size);                      //get weight matrix core
getSiteScores();                         //calculate individual (core-)scores of sites
if (!i_score) mth=getthreshold(0,msens); //calculate matrix threshold score
  else mth=i_score;                      //take score as threshold (if defined)
cth=getthreshold(1,csens);               //calculate core threshold score
nop=i_nop;                               //non-occurrence panelty

//show();

}

void matrix::getSiteScores() {
//---------------------------------------------------------------------------
// Get (core-)scores of individual sites comprising the weight matrix
//---------------------------------------------------------------------------

for (unsigned int i=0; i<site.getnumber(); i++) {
  site_score.push_back(getScore(site.getseq(i).data(),0,0));    // site-scores
  core_score.push_back(getCoreScore(site.getseq(i).data(),0,0)); // core-scores
  //cout << "Core Score: " << core_score[i] << endl;
  //cout << "Site Score: " << site_score[i] << endl;
} //for

}

void matrix::getCore(unsigned char core_number) {
//---------------------------------------------------------------------------
// Get the most conserved matrix positions:
// core_number: size of core
//---------------------------------------------------------------------------
  vector<float> co; //temporary vector
  float co_max;
  int pos=0;

  core_size=core_number; //save core_number in class variable

  //copy information content in temporary vector
  for (unsigned int i=0; i<msize; i++) {
    co.push_back(Rs[i]);
  } //for


  for (unsigned int i=0; i<core_number; i++) {

    //search for the highest information value:
    co_max=-99999.9; //very low value as default

    for (unsigned int j=0; j<msize; j++) {
      if (co[j]>co_max) {
        co_max=co[j];
        pos=j;
      } //if
    } //for j

    co[pos]=-99999.9;  //"remove" this value from the information vector
    core.push_back(pos);
    //cout << "core: " << pos << endl;
  } //for


}

float matrix::getthreshold(unsigned char type, float sens) {
//---------------------------------------------------------------------------
// Get threshold score from sensitifity
// type: 0: overall score, 1: core-score
// sens: sensitifity value should be between 0 and 1
//---------------------------------------------------------------------------
  vector<float> th;       //vector with threshold scores
  vector<float> ss;       //vector with sensitivities   
  float site_sens;        //sensitivity of PWM sites
  int snr1=-1;            //number of sites with flanking sensitivities
  int snr2=-1;            //number of sites with flanking sensitivities
  float m, t;             //parameters of a straight line function y=mx+t
  float tolerance=0.001;  //tolerance to avoid mistakes by rounding

  //copy (core-)scores of all weight matrix sequences into a temporary vector:
  switch(type) {
    case 0: //overall scores
      for (unsigned int i=0; i<site.getnumber(); i++) {
        th.push_back(site_score[i]); //attach scores into vector
      } //for i
    break;

    case 1: //core scores
      for (unsigned int i=0; i<site.getnumber(); i++) {
        th.push_back(core_score[i]); //attach core-scores into vector
      } //for i
    break;

  } //switch

  sort(th.begin(), th.end());  //sort scores using STL
  
  //find flanking sensitifities:
  for (unsigned int i=0; i<th.size(); i++) {
    site_sens=(float)(th.size()-i)/(float)th.size();
    ss.push_back(site_sens); //attach sensitivity into vector
    if ((sens>=site_sens)&&(snr1==-1)) {
      snr1=i;
    }   
  //if (!type) cout << ss[i] << " ! " << th[i] << endl;
  }
  if (snr1==0) snr1=1;                    //site with lowest scoree           
    else if (snr1==-1) snr1=th.size()-1;  //no site found (score too high)
  snr2=snr1-1;                            //take site with next lower score as other flanking site 
  
  //if (!type) cout << snr1 << " ? " << snr2 << endl;
  
  //calculate parameters of a straight line function y=mx+t     
  m=(th[snr2]-th[snr1])/(ss[snr2]-ss[snr1]);
  t=(ss[snr2]*th[snr1]-ss[snr1]*th[snr2])/(ss[snr2]-ss[snr1]);
  //if (!type) cout << m << " x " << t << endl;

  return (m*sens+t-tolerance); //return score (mimus a small tolerance to avoid mistakes by rounding)

}


void matrix::show() {
//---------------------------------------------------------------------------
// Print matrix contents and features
//---------------------------------------------------------------------------

string ch;
cout << "Pattern Information" << endl;
cout << "name:  " << mname << endl;
cout << "size:  " << (int)msize << endl;
cout << "model: " << (int)model << endl;
cout << "weight matrix sensitivity: " << msens << endl;
cout << "weight matrix score threshold: " << mth << endl;
cout << "core sensitivity: "; if (core_size) cout << csens <<endl; else cout <<"-"<< endl;
cout << "core size: " << (int)core_size << endl;
cout << "core positions: ";
if (core_size) {
  for (unsigned int i=0; i<core.size(); i++) {
    cout << (core[i]+1)<< ", ";
  }
} else cout << "-"; cout << endl;
cout << "core score threshold: " << cth << endl;
cout << "Rsequence: " << Rsequence << endl;
cout << "Rsequence (bias correction): " << Rsequence_bias << endl;
cout << "Maximal Score: " << round(maxscore,2) << endl;
cout << "Minimal Score: " << round(minscore,2) << endl;
cout << endl;
cout << "Position Weight Matrix:" << endl;
cout << setw(6) << "pos(l)" << setw(6) << "A" << setw(6) << "C" << setw(6) << "G" << setw(6) << "T" << setw(6) << "Rs(l)" << endl;

  for (unsigned int i=0; i<site.getseq(0).size(); i++) {
    cout << setw(6) << i+1;
    for (unsigned int n=0; n<4; n++) {
      cout << setw(6)  << round(value[i][n].score,2);
      //cout << setw(6)  << round(value[i][n].number_bias,2);
	  //cout << setw(6)  << round(value[i][n].freq_bias,2);
    } // for n
    if (model==5) cout << setw(6) << round(Rs_bias[i],1); else cout << setw(6) << round(Rs[i],1);
	//cout << setw(6) << round(Rs[i],1);
    cout << endl;
  } //for i
cout << endl;
cout << "Binding Sites:" << endl;
cout << "Sequence/Score/Similarity/Core-Score:" << endl;
for (unsigned int i=0; i<site.getnumber(); i++) {
  cout << site.uppercase(i) << setw(6) << round(site_score[i],2)<< setw(6)
       << round(site_score[i]/maxscore*100,2)<<"%"<< setw(6)
       << round(core_score[i],2) << endl;
}//for
cout << "Number of sites: " << site.getnumber() << endl;

cout << endl;
cout << "Sequence Information:" <<endl;
cout << "Genomic composition: A:" << round(Pb[0],2) <<
                            " C:" << round(Pb[1],2) <<
                            " G:" << round(Pb[2],2) <<
                            " T:" << round(Pb[3],2) << endl;

}

float matrix::getRsequence() {
//---------------------------------------------------------------------------
// Returns the information content of the whole site model (Rsequence)
//---------------------------------------------------------------------------
  return Rsequence;
}

float matrix::getScore(const char *seq, int offset, int dir) {
//---------------------------------------------------------------------------
// Returns the score of a sequence to the weigth matrix
// seq: pointer to start of genome sequence
// offset: distance from start of genome sequence
// score_type: type of weight matrix model
// dir: direction (+ or - strand)
//---------------------------------------------------------------------------

float Ri=0.0;
signed char mat_row;  //matrix row number

switch (dir) {
  case 0: //(+ direction)
    for (unsigned int i=0; i<msize; i++) {
      mat_row=charval[seq[i+offset]];
      Ri+=value[i][mat_row].score;         //build sum of score values
    }
  break;
  case 1: //(- direction=reverse complement)
    for (unsigned int i=0; i<msize; i++) {
      mat_row=3-charval[seq[i+offset]];
      if (mat_row>=0) Ri+=value[msize-i-1][mat_row].score; //build sum of score values (reverse)
    }
  break;
}
return Ri;
}

float matrix::getCoreScore(const char *seq, int offset, int dir) {
//---------------------------------------------------------------------------
// Returns the score of a sequence to the weigth matrix
// seq: sequence stretch
// score_type: type of weight matrix model
// dir: direction (+ or - strand)
//---------------------------------------------------------------------------

float Ri=0.0;
signed char mat_row;  //matrix row number
int pos;


switch (dir) {
  case 0: //(+ direction)
    for (unsigned int i=0; i<core_size; i++) {
      pos=core[i];
      mat_row=charval[seq[pos+offset]];
      Ri+=value[pos][mat_row].score;         //build sum of score values
	  //cout << "pos:" <<pos << "mat_row:"<<(int)mat_row<<"score:"<<value[pos][mat_row].score <<"|";
    }
  break;

  case 1: //(- direction=reverse complement)
    for (unsigned int i=0; i<core_size; i++) {
      pos=core[i];
      mat_row=rcharval[seq[msize-pos+offset-1]];
	  //cout << (int)mat_row << " - ";
      Ri+=value[pos][mat_row].score; //build sum of score values (reverse)
	  //cout << (int)pos << " - " << (int)mat_row << " - " << value[pos][mat_row].score <<endl;
    }
  break;

}

return Ri;
}


void matrix::addsite(string seq, string ID) {
//---------------------------------------------------------------------------
// Add site to matrix
//---------------------------------------------------------------------------
site.add(ID, seq);

//site[sitesum].seq=uppercase(seq);

}

float matrix::calehnb (int n) {
//---------------------------------------------------------------------------
// Calculates small sample size correction
// Algorithm taken with some modifications from rseq.c
// (Schneider T.D. et. al 1986, J. Mol. Biol. 188: 415-431)
//---------------------------------------------------------------------------
double hg, ehnb, varhnb;
double logfact[maxsites];
double mplog2p[maxsites];

//get logs of genome probabilities:
double logpa=log(Pb[0]);
double logpc=log(Pb[1]);
double logpg=log(Pb[2]);
double logpt=log(Pb[3]);

//find genomic entropy
hg=-(Pb[0]*logpa + Pb[1]*logpc + Pb[2]*logpg + Pb[3]*logpt)/log(2.0);
ehnb=0.0; //start error entropy at zero

//make table of log of n factorial up to n and entropies for nb/n
logfact[0]=0.0;
mplog2p[0]=0.0;

double logi;
double logn=log((double)n);
double nlog2=n*log(2.0);
for (unsigned int i=1; i<=n; i++) {
  logi=log((double)i);
  logfact[i]=logfact[i-1]+logi;
  mplog2p[i]=i*(logn-logi)/nlog2;
} //for i

long na, nc=0, ng=0, nt=0; //number of bases in a site

na=n; //begin by looking at the combination with all a: na=n

float pnb; // multinomial probability of a combination of na, nc, ng, nt
float hnb;       //entropy for a combination of na..nt
float pnbhnb;    //pnb*hnb, an intermediate result;  
float sshnb=0.0; //sum of squares of hnb
float total=0.0;
long counter=0;
bool done=false;
do {
  pnb=exp(logfact[n]-logfact[na]-logfact[nc]-logfact[ng]-logfact[nt]+
          na*logpa+nc*logpc+ng*logpg+nt*logpt); 
  hnb=mplog2p[na]+mplog2p[nc]+mplog2p[ng]+mplog2p[nt];
  pnbhnb=pnb*hnb;    //
  ehnb+=pnbhnb;      
  sshnb+=pnbhnb*hnb; //sum of squares of hnb
  counter++;
  total+=pnb;

  if (nt>0) {
    if (ng>0) {
      ng--;
      nt++;
    } else if (nc>0) {
      nc--;
      ng=nt+1;
      nt=0;
      } else if (na>0) {
        na--;
        nc=nt+1;
        nt=0;
        } else done=true;
      }else {
       if (ng>0) {
         ng--;
         nt++;
       } else if (nc>0) {
         nc--;
         ng++;
       } else {
         na--;
         nc++;
       }
     }          
} while (!done); //do

varhnb=sshnb-ehnb*ehnb;
return (float)ehnb;
}

void matrix::cons2matrix(string cons) {
//---------------------------------------------------------------------------
// Calculates a consensus matrix from a IUPAC consensus
// cons=IUPAC string
//---------------------------------------------------------------------------
  msize=cons.size(); //define consensus matrix size
  cout << "cons: " << cons << endl;
  for (unsigned int i=0; i<cons.size(); i++) {

    //set all numbers to zero as default:  
    for (unsigned int j=0; j<4; j++) {
      value[i][j].number=0;
    }
  
    switch(cons.at(i)) {
      case 65: value[i][0].number=1; //A->A
      break;
      case 66: value[i][1].number=1; //B->[CGT]
               value[i][2].number=1;
               value[i][3].number=1;  
      break;
      case 67: value[i][1].number=1; //C->C
      break;
      case 68: value[i][0].number=1; //D->[AGT]
               value[i][2].number=1;
               value[i][3].number=1;
      break;
      case 71: value[i][2].number=1; //G->G
      break;
      case 72: value[i][0].number=1; //H->[ACT]
               value[i][1].number=1;
               value[i][3].number=1;
      break;
      case 75: value[i][2].number=1; //K->[GT]
               value[i][3].number=1;
      break;
      case 77: value[i][0].number=1; //M->[AC]
               value[i][1].number=1; 
      break;
      case 78: value[i][0].number=1; //N->[ACGT]
               value[i][1].number=1;
               value[i][2].number=1;
               value[i][3].number=1;
      break;
      case 82: value[i][0].number=1; //R->[AG]
               value[i][2].number=1;
      break;
      case 83: value[i][1].number=1; //S->[CG]
               value[i][2].number=1;
      break;
      case 84: value[i][3].number=1; //T->T
      break;
      case 86: value[i][0].number=1; //V->[ACG]
               value[i][1].number=1;
               value[i][2].number=1; 
      break;
      case 87: value[i][0].number=1; //W->[AT]
               value[i][3].number=1;
      break;
      case 89: value[i][1].number=1; //Y->[CT]
               value[i][3].number=1;
      break;

    }; //switch

  } //for i

}

int matrix::getMismatch(const char *seq, int offset, int dir, int mm_max) {
//---------------------------------------------------------------------------
// Returns the number of mismatches to a IUPAC consensus
// seq: pointer to start of genome sequence
// offset: distance from start of genome sequence
// score_type: type of weight matrix model
// dir: direction (+ or - strand)
// mm_max: maximum mismatches
//---------------------------------------------------------------------------
int cm=0;       //consensus match
int mmcount=0;  //number of mismatches

switch (dir) {
  case 0: //(+ direction)
    for (unsigned int i=0; i<msize; i++) {
	  cm=0; //reset consensus match
      switch (seq[i+offset]) {
        case 65: cm+=value[i][0].number; //match for A
        break;
        case 67: cm+=value[i][1].number; //match for C
        break;
        case 71: cm+=value[i][2].number; //match for G
        break;
        case 84: cm+=value[i][3].number; //match for T
        break;
      }
      if (!cm) mmcount++;
	  if (mmcount>mm_max) i=msize; //abort loop if maximum mismaches
    }
  break;
  case 1: //(- direction=reverse complement)
    for (unsigned int i=0; i<msize; i++) {
	  cm=0; //reset consensus match
      switch (seq[i+offset]) {
        case 65: cm+=value[msize-i-1][3].number; //match for A(rev T)
        break;
        case 67: cm+=value[msize-i-1][2].number; //match for C(rev G)
        break;
        case 71: cm+=value[msize-i-1][1].number; //match for G(rev C)
        break;
        case 84: cm+=value[msize-i-1][0].number; //match for T(rev A)
        break;
      }
	  if (!cm) mmcount++;
	  if (mmcount>mm_max) i=msize; //abort loop if maximum mismaches
    }
  break;
}
return mmcount;
}


void matrix::hist(string histfile, int Smin, int Smax, float step) {
//---------------------------------------------------------------------------
// Write histogram file of individual scores
// histfile:   file path/name
// score_type: score type
// Smin:       minimum histogram score
// Smax:       maximum histogram score
// step:       histogram steps
//---------------------------------------------------------------------------
/*
  int histnr=int(Smax/step);     // number of histogram entries
  int hist[histnr];              // histogram data array
  int h=0;                       // histogram position

  //Clear histogram data:
  for (unsigned int i=0; i<=histnr; i++) {
    hist[i]=0;
  } //for i

  //calculate histogram:
  for (float s=Smin; s<=(float)Smax; s+=step) {
    for (unsigned int i=0; i<site.getnumber(); i++) {
        if ((sitescore[i]<(s+step/2))&&(sitescore[i]>=(s-step/2))) {
        hist[h]++; //increase histogram entry
      }
    } //for i
    h++; //next histogram position
  } //for s


  //write histogram data to file:
  h=0;
  ofstream output(histfile.data()); //open histogram datafile
  for (float s=Smin; s<=(float)Smax; s+=step) {
    output << s << "," << hist[h]<<endl; //write comma separated data
    h++; //next histogram position
  }
  output.close(); //close histogram datafile
*/
}


