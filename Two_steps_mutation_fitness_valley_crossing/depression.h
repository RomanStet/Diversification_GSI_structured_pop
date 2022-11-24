#ifndef DEPRESSION_H
#define DEPRESSION_H

#include <vector>
#include <iostream>
#include "MersenneTwister.h"
using namespace std;


// global variables

#define fileRead "parameters.csv"     // names of input

// definition of structure "chr" representing a chromosome:
// "S" is the S-allele

struct chr
{
   int S_pol; //S-locus pollen part
   int S_pis; //S-locus pistil part
   int neutral_all; // allele at the neutral locus, idenpendent from the S-locus
};

struct new_SC_haplo
{
	int female_spe;
	int male_spe;
	int presence;
	int mutation_type;
};

struct Nall // allows to keep track of of the identity and the frequency of neutral alleles in the population
{
	int all; // identity of the allele
	int freq; // frequency of the allele
};

struct result
{
	double frSC; //frequency of SC
    double delta; // inbreeding depression
	int nb_haplo; // final number of S-haplotypes
};

// Prototypes of functions

void openFileE();
void closeFileE();

bool readFile(int &Numr, int &Nr,int &pr, int &max_allelesr, int &alleles_initr, double &dpr, double &dsr, double &ar,
                     int &itr, double &stepr, double &Uscr, int &mUscr,
                     int &NbGenr, double &deltar, int &stepgr, double &Unr, int &position_stream, int &state);


result recursion(int numsimulv, int iv, int Nv, double av,
                double Uscv, int mUscv, int NbGenv, double deltav, int stepgv, double Unv,
                double dpv , double dsv , int pv, int nbr_allelesv, int nbr_alleles_initv);


double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
void cntl_c_handler(int bidon);


#endif
