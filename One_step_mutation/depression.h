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
   int S; //S-locus
};

struct result
{
	int nb_haplo;
};

// Prototypes of functions

void openFileE();
void closeFileE();

bool readFile(int &Numr, int &Nr, int &pr, double &dpr, double &dsr, int &itr, double &stepir, double &Uscr, int &NbGenr, double &deltar, int &stepscr,
                 int &position_stream, int &state);


result recursion(int numsimulv, int iv, int Nv,
                double Uscv, int NbGenv, double deltav, int stepv,
                double dpv , double dsv , int pv);


double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
void cntl_c_handler(int bidon);


#endif
