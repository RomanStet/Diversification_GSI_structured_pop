#include "depression.h"
#include "MersenneTwister.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#ifdef __unix__
#include <sys/stat.h>
#else
#include <direct.h>
#endif
using namespace std;

FILE * fileP;

// Random number generator:

MTRand rnd; // rnd tout court pour rétablir une graine aléatoire rnd(X) avec x un int pour une graine fixe

// Pointers on input and output files:

FILE * fileE;

int main()
{

    // Parameters:
    int Nt, NbGen, stepsc, it, p, position_stream, state;
	double stepi, Usc, deltai, delta, dp, ds;
	result res;

	// Opens input and output files:

	bool end;
	end = false;
	bool test;
	bool stop;

	int i;
    int no;
    int numsimul;

	do	{
        //reads parameter values from input file:
        openFileE();
 		end = readFile(numsimul ,Nt, p, dp, ds, it, stepi, Usc, NbGen, deltai, stepsc, position_stream, state);
        closeFileE();
		// end = true if end of input file

        // actual population size
        int size_pop = Nt-Nt%p;

		delta = deltai;
        stop = false;

		if (!end==true) { // if the end of the file in not reached

            if (state==0) { // if the simulation has never been launched

                // Writes parameter values in output file:
                char namefolder[50];
                sprintf(namefolder,"%d_folder",numsimul);

                #ifdef __unix__
                mkdir(namefolder,07777);
                #else
                mkdir(namefolder);
                #endif

                char nameFileStart[256];
                stringstream nameStart;
                nameStart << numsimul << "_folder/" << "start.txt";
                nameStart >> nameFileStart;

                if(fopen(nameFileStart,"r")==NULL) {
					
                    ofstream foutstart(nameFileStart);
                    foutstart << "Start" << endl;
                    foutstart.close();

                    cout << "Processing " << numsimul << endl;

                    // creating name of output file (which indicates the parameter values):
                    char nameFileParamGeneral[256];
                    stringstream nameFPG;
                    nameFPG << numsimul << "_folder/" << numsimul << "_paramgeneral.txt";
                    nameFPG >> nameFileParamGeneral;
                    ofstream foutparamgene(nameFileParamGeneral);
                    // writing the parameter file of the simulation
					foutparamgene  << "no" << " " << "iterations" << " " << "delta" << " " << "nbr_haplotype_end" << endl;
                    foutparamgene.close();

                    no = 0; // initialisation of a counter of the number of iterations of the simulation

                    int itp=it;
                    for (i = 0; i < itp; i++) {
                        no++;
                        // Simulation:
							cout << no << ", delta: " << delta << endl;
							res = recursion(numsimul, no, size_pop,
											Usc, NbGen, delta, stepsc,
											dp, ds, p);

							foutparamgene.open(nameFileParamGeneral,std::ofstream::app);
							foutparamgene << numsimul << " " << no << " " << delta << " " << res.nb_haplo << endl;
							foutparamgene.close();

							delta += stepi;
                    }

                    openFileE();
                    fseek(fileE,position_stream,SEEK_CUR);
                    putc('#',fileE);
                    closeFileE();

                    char nameFileEnd[256];
                    stringstream nameEnd;
                    nameEnd << numsimul << "_folder/" << "end.txt";
                    nameEnd >> nameFileEnd;
                    ofstream foutend(nameFileEnd);
                    foutend << "Finished" << endl;
                    foutend.close();

                    remove(nameFileStart);

                    cout << endl;
                } else { // if there already is a start file
                    cout << "Work already processing on  " << numsimul << endl;
                }
            } else if (state==1) { // the simulation has already been launched before and it is resumed
                cout << "Restart processing " << numsimul << endl;

                no = 0; // initialisation of a counter of the number of iterations of the simulation

                char nameFileParamGeneral[256];
                stringstream nameFPG;
                nameFPG << numsimul << "_folder/" << numsimul << "_paramgeneral.txt";
                nameFPG >> nameFileParamGeneral;


					cout << "delta scanned" << endl;
					int itp=it;

					fileP = fopen(nameFileParamGeneral,"r");

					char testchar[20];
					for (int k=0;k<7;k++) {
						fscanf(fileP,"%s",testchar);
					}
					int y=1;
					int idepart=0;

						while (y==1) {
							cout << "\n";
							y = fscanf(fileP,"%d",&numsimul); // simulation number
							if (y!=-1) {
								fscanf(fileP,"%d",&no); // iteration number
								fscanf(fileP,"%lf",&delta); // value of delta
								fscanf(fileP,"%d",&res.nb_haplo); // final number of haplotypes
								fscanf(fileP,"%d",&test); // test
								fscanf(fileP,"%d",&stop); // stop
								cout << numsimul <<  " " << no <<  " " << delta <<  " " << no <<  " " << res.nb_haplo << "\n";
								if (delta!=floor(delta)) {
									idepart++;
								}
							}
						}

						fclose(fileP);

						delta += stepi;

					for (i = idepart; i < itp; i++) {
						 no++;
						 // Simulation:

							cout << no << ", delta: " << delta << endl;
							res = recursion(numsimul, no, size_pop,
											Usc, NbGen, delta, stepsc,
											dp, ds, p);


							fileP = fopen(nameFileParamGeneral,"a");
							fprintf(fileP, "%i %i %lf %i \n",numsimul,no,delta,res.nb_haplo);
							fclose(fileP);

							delta += stepi;

					}
			}

				openFileE();
				fseek(fileE,position_stream,SEEK_CUR);
				putc('#',fileE);
				closeFileE();

				char nameFileEnd[256];
				stringstream nameEnd;
				nameEnd << numsimul << "_folder/" << "end.txt";
				nameEnd >> nameFileEnd;
				ofstream foutend(nameFileEnd);
				foutend << "Finished" << endl;
				foutend.close();

				char nameFileStart[256];
				stringstream nameStart;
				nameStart << numsimul << "_folder/" << "start.txt";
				nameStart >> nameFileStart;
				remove(nameFileStart);

				cout << endl;
        }

    // end = true if end of input file
	} while (!end);
	return 0 ;
}
