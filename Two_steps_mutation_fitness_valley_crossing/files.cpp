// Functions to open input and output files,
// read parameter values and write them in output file

#include "depression.h"
#include <iostream>
#include <fstream>
#include <sstream>
#ifdef __unix__
#include <unistd.h>
#else
#include <Windows.h>
#endif
using namespace std;

extern FILE * fileE;
FILE * filestop;

//Opens input file:

void openFileE()
{
 int counter_attempts=0;

    while (counter_attempts<10) {

        filestop = fopen("stop.txt","r");

        if (filestop==NULL) { // if the "stop" file does not exist, it is created and the "result" file is opened
            ofstream foutstop("stop.txt");
            foutstop.close();

            fileE = fopen(fileRead,"r+");
            if (!fileE) {
                cout << "The file " << fileRead << " doesn't exist!" << endl;
                remove("stop.txt");
            }

            counter_attempts = 10; // the "parameter" file either does not exist or is already open, we may stop

        } else { // if the "stop" file exists
            // we wait a bit to retry
            cout << "Sleep " << counter_attempts << endl;
            #ifdef __unix__
            usleep(10000);
            #else
            Sleep(10000);
            #endif

            counter_attempts++;
        }
    }
}

void closeFileE()
{
    fclose(fileE);

    #ifdef __unix__
    usleep(2000);
    #else
    Sleep(2000);
    #endif

    remove("stop.txt");

}


// Reads parameter values from input file,
// returns 1 if it reaches the end of input file, else returns 0:

bool readFile(int &Numr, int &Nr,int &pr, int &max_allelesr, int &alleles_initr, double &dpr, double &dgr, double &ar,
                     int &itr, double &stepr, double &Uscr, int &mUscr,
                     int &NbGenr, double &deltar, int &stepgr, double &Unr, int &position_stream, int &stater)
{
	int x;
	bool term;
	do {
        x = fgetc(fileE); //read caracters one by one
	}
	while ( !( (x == '*') || (x == ':') || (x == EOF) ) ); // until you reach a * or the end of the file
		// Lines with parameter sets must begin with *

	if (x == EOF) {
		cout << "\nEnd of input file\n";
		term = true;
	} else if  ((x == '*') || (x == ':')) {

        fseek(fileE,-1,SEEK_CUR);
        position_stream = ftell(fileE);
        putc('!',fileE);
        fseek(fileE,0,SEEK_CUR);

        fscanf(fileE,",%d,",&Numr); // simulation number
        fscanf(fileE,"%d,",&Nr); // population size
        fscanf(fileE,"%d,",&pr); // number of demes
        fscanf(fileE,"%d,",&max_allelesr); // maximum number of SI haplotypes
        fscanf(fileE,"%d,",&alleles_initr); // initial number of SI haplotypes
        fscanf(fileE,"%lf,",&dpr); // pollen dispersion
        fscanf(fileE,"%lf,",&dgr); // seed dispersion
        fscanf(fileE,"%lf,",&ar); // selfing rat
        fscanf(fileE,"%d,",&itr); // number of iterations
        fscanf(fileE,"%lf,",&stepr); // step to search for delta
        fscanf(fileE,"%lf,",&Uscr); // mutation rate of the SC alleles
		fscanf(fileE,"%d,",&mUscr); // to allow for multiple mutations generating SC alleles
        fscanf(fileE,"%d,",&NbGenr); // number of generations
        fscanf(fileE,"%lf,",&deltar); // inbreeding depression
        fscanf(fileE,"%d,",&stepgr);  // number of generations between the generation where results are written
		fscanf(fileE,"%lf",&Unr); // mutation rate at the neutral locus

        if (x == '*') {
            stater = 0;
        }
        if (x == ':') {
            stater = 1;
        }
        term = false;
    }
	return term;
}
