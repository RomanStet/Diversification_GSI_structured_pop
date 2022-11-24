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
 int attempts_counter=0;

    while (attempts_counter<10) {

        filestop = fopen("stop.txt","r");

        if (filestop==NULL) { // if the stop file does not exist it is created
            ofstream foutblocage("stop.txt");
            foutblocage.close();

            fileE = fopen(fileRead,"r+");
            if (!fileE) {
                cout << "The file " << fileRead << " doesn't exist!" << endl;
                remove("stop.txt");
            }

            attempts_counter = 10;  // the "parameter" file either does not exist or is already open, we may stop

        } else { // if the "stop" file exists
            // we wait a bit to retry
            cout << "Sleep " << attempts_counter << endl;
            #ifdef __unix__
            usleep(10000);
            #else
            Sleep(10000);
            #endif

            attempts_counter++;
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

bool readFile(int &Numr, int &Nr,int &pr, double &dpr, double &dsr,
                     int &itr, double &stepir, double &Uscr,
                     int &NbGenr, double &deltar, int &stepscr, int &position_stream, int &stater)
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
        fscanf(fileE,"%lf,",&dpr); // pollen dispersal
        fscanf(fileE,"%lf,",&dsr); // seed dispersal
        fscanf(fileE,"%d,",&itr); // number of iterations
        fscanf(fileE,"%lf,",&stepir); // step to search for delta
        fscanf(fileE,"%lf,",&Uscr); // mutation rate for S-haplotypes
        fscanf(fileE,"%d,",&NbGenr); // number of generations
        fscanf(fileE,"%lf,",&deltar); // inbreeding depression
        fscanf(fileE,"%d,",&stepscr);  // number of generations at which variables are measured

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
