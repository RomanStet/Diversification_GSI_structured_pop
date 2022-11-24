#include "depression.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <cmath>
#include <csignal>
#include <algorithm>
using namespace std;

extern MTRand rnd;
extern FILE * fileE;


/* function recursion: iterates the life cycle.

   During NbGenv generations: specificity mutations

   Other parameters are:
    numsimul: simulation number
    iv: number of iterations for the value of delta (not used and set to 1)
    Nv: number of diploid individuals
	Uscv: mutation rate from mutation rate from Si to Sj
	stepv: number of generations between measures
	dpv: pollen dispersal
	dsv: seed dispersal
	pv: number of demes
*/

result recursion(int numsimulv, int iv, int Nv,
                double Uscv, int NbGenv, double deltav, int stepv,
                double dpv, double dsv, int pv)
{
    //i,j,k,gen:counters, chrmut: mutated chr, chrm= maternal chr, chrd= paternal chr, ind = individual
	unsigned int un_i, un_j, un_k;
	int i, gen, chrm, chrd, ind, nb, ind_chr, Si, Sj, mom, dad, modif_pos, nbr_haplotypes;

	double rd, D, D_patch, D_patch_var, div_temp_patch;
	result Res;

	bool stop_signal = false;

	vector <int> S_alleles; S_alleles.clear();  //keeps track of the number of S-alleles in the population
	vector <int> S_freq_alleles ; S_freq_alleles.clear(); //keeps track of the number of S-alleles with a certain frequency (>0.01) in the population
    vector <int> S_alleles_patch; S_alleles_patch.clear(); //keeps track of the number of S-alleles in a deme
	vector <double> freq; // vector of the frequency of each haplotype in the population
	vector <double> freq_loc; // vector of the frequency of each haplotype in a deme
	freq.resize(1000000); // nbr of possible S-haplotypes in the population
	freq_loc.resize(1000000); // nbr of possible S-haplotypes in a deme
	vector <double> expected_hetero; // expected heterozygosity in each deme to compute FST at the S-locus
	expected_hetero.resize(pv);

	//deme size
	int size_patch = Nv/pv; // Nv has been previously resized so that it gives an int when divided by pv
	int two_size_patch = 2*size_patch; // Nv has been previously resized so that it gives an int when divided by pv
	int twoN = 2*Nv; // Nv has been previously resized so that it gives an int when divided by pv

	int compt_patch =0; // deme counter

	// to compute the diversity in the demes
	double nS_patch_bar; // for the mean number of SI per deme
	double nS_patch_var; // for the variance in the number of SI per deme
	double tot_expected_hetero; // expected heterozygosity at the S-locus in the population
	double sum_expected_hetero; // sum of expected heterozygosities at the S-locus in all demes
	double FST; // FST for S-haplotypes
	double p_disp_bar; // to compute the actual pollen dispersal rate
    double s_disp_bar; // to compute the actual seed dispersal rate

	// population: table of 2Nn chromosomes:
	chr * pop = new chr [twoN]; // vector of the population
   	chr * temp = new chr [twoN]; // temporary vector of the population

	// creating name of output file (which indicates the parameter values):
    char nameFileParam[256];
	stringstream nameFP;
    nameFP << numsimulv << "_folder/" << numsimulv << "_param_" << iv <<".txt";
    nameFP >> nameFileParam;
	ofstream foutparam(nameFileParam);
    // writing the parameters of the simulation
	foutparam  << "no"<<" " << "i" <<" "<< "N"<< " " <<
                     "Usc" << " " << "nb_gen"<< " " << "delta"<< " " <<
                     "step"<< " " << "dp"<< " " << "ds"<< " " << "p"<< endl;
    foutparam << numsimulv << " " << iv << " " << Nv << " " <<
                      Uscv << " " << NbGenv << " " << deltav << " " <<
                     stepv <<  " " << dpv << " " << dsv << " " << pv << endl;
    foutparam.close();


    // creating the result file
	char nameFileResult[256];
	stringstream nameFResult;
	nameFResult << numsimulv << "_folder/" << numsimulv << "_result_" << iv <<".txt";
	nameFResult >> nameFileResult;
	ofstream foutresult(nameFileResult);
	//writing the headers of the result file
    foutresult <<"generation"
        << " " << "haplotype_number" << " " << "haplotype_diversity" << " " << "neutral_diversity"
        << " " << "FST" << " " << "FST_neutre" << " " << "local_haplotype_number" << " " << "local_haplotype_diversity" << " " << "local_neutral_diversity" << " " << "local_haplotype_number_variance" << " " << "local_haplotype_diversity_variance"<< " " << "local_neutral_diversity_variance"
        << " " << "disp_pol" << " " << "disp_graine" <<  endl;
    foutresult.close();

    // creating a file allowing to take a picture of the population
	char nameFilePicture[256];
	stringstream namePicture;
	namePicture << numsimulv << "_folder/" << numsimulv << "_picture_" << iv <<".txt";
	namePicture >> nameFilePicture;
	ofstream foutPicture(nameFilePicture);
	// headers
	foutPicture <<"phase" << " " << "patch" << " " << "S_1" << " " << "S_2" << endl;
    foutPicture.close();

    //////////// Population initialization
	//// All individuals are heterozygous
    ///////////

    for (i = 0; i < Nv; i++) { // for each individual in the population
        nb = 2*i; // postion of the first chromosome of the individual
        pop[nb].S = rnd.randInt(1000000); // drawing a haplotype among all possible haplotypes
        do {
            pop[nb + 1].S = rnd.randInt(1000000);
        } while (pop[nb + 1].S == pop[nb].S); //different S-haplotype at the second chromosome
    }

    foutPicture.open(nameFilePicture,std::ofstream::app);
    for (i = 0; i < twoN; i++)
	{
        foutPicture << "S" << " " << ((i-(i%two_size_patch))/two_size_patch)+1 << " " << pop[i].S << " " << pop[i+1].S << " " << endl;
		i++;
	}
    foutPicture.close();

    if (!stop_signal==true) //if no stop signal


    for ( gen = 0; gen<NbGenv; gen++)
	{ // for all generations
		
		//////////////
		////mutations
		//////////////
        for (i = 0; i < twoN; i++)
		{ // for each chromosome

            // mutation at the S-locus
            rd = rnd.randExc();
            if (rd < Uscv) { // if there is a mutation at the S-locus
                        pop[i].S = rnd.randInt(1000000);
                    }
        }


		//////////////
		//// Creating next generation
		//////////////

        p_disp_bar = 0;
        s_disp_bar = 0;

        // sampling the next generation:
        for (ind = 0; ind < Nv; ind++)
		{ // for each individual
            ind_chr = 2*ind; // to work with chromosomes
            // sampling the mother:

					if (rnd.randExc() >= dsv) { // if no seed dispersal
						chrm = (ind_chr-ind_chr%two_size_patch)+rnd.randInt(two_size_patch-1); // a chromosome is drawn in the same deme
						mom = chrm/2; // the individual carying this chromosome
					} else { // if seed dispersal, the mother is drawn from a different deme
						do{ // sampling the mother
							chrm = rnd.randInt(twoN-1); // a chromosome is drawn in the population
							mom = chrm/2; // the individual carying this chromosome
						} while (ind_chr-ind_chr%two_size_patch==chrm-chrm%two_size_patch); // while the individual is from the same deme
					}

				// to check if the mother can be fertilized
				if (chrm % 2 == 0) { // if the chromosome is the first of this mother
					modif_pos = 1; // complementary chromosome
				} else { // if the chromosome is the second of this mother
					modif_pos = -1; // complementary chromosome
				}
				Si = pop[chrm].S; // records the haplotype of one chromosome of the mother
				Sj = pop[chrm + modif_pos].S; // records the haplotype of the other chromosome of the mother

			// dispersal
            if (ind_chr-ind_chr%two_size_patch!=chrm-chrm%two_size_patch) {
                s_disp_bar++;
            }

            // sampling the father:
                do {
                    if (rnd.randExc() >= dpv) { // if no pollen dispersal
                        chrd = (chrm-chrm%two_size_patch)+rnd.randInt(two_size_patch-1); // drawing the parternal haplotype in the local deme
                        dad = chrd/2; // individual carrying this haplotype
                    } else {// if pollen dispersal
                        do {
                            chrd = rnd.randInt(twoN-1); // drawing the parternal haplotype in the population
                            dad = chrd/2; // individual carrying this haplotype
                        } while (chrm-chrm%two_size_patch==chrd-chrd%two_size_patch); // while no pollen from a different deme
                    }
                } while ( (mom == dad) || (pop[chrd].S == Si) || (pop[chrd].S == Sj) ); //while the haplotype of the father is not compatible with the haplotypes of the mother

                if (chrm-chrm%two_size_patch!=chrd-chrd%two_size_patch) { // if pollen dispersal
                    p_disp_bar++;
                }

				// drawing the maternal chromosome and assigning it to its position in the next generation
				if (rnd.rand() < 0.5)
				{
					temp[ind_chr] = pop[chrm]; // if it is one chromosome
				}
				else
				{
					temp[ind_chr] = pop[chrm+modif_pos]; // if it is the other chromosome
				}
				temp[ind_chr+1] = pop[chrd]; // assigning the paternal chromosome to its position in the next generation
		}

		// dispersal rate
        s_disp_bar /= Nv;
        p_disp_bar /= Nv;

        for (i = 0; i < twoN; i++) // for each chromosome
            pop[i] = temp[i]; // replacing the chromosomes by the new ones

        //////////////
		//// Various measures and writing of the results
		//////////////
		//// Computing of frequency of S-alleles
		// Initializing
		S_alleles.clear();
		S_freq_alleles.clear();
        nS_patch_bar=0;
		nS_patch_var=0;
        D_patch=0;
        D_patch_var=0;
		compt_patch=0;
		tot_expected_hetero=0;
		sum_expected_hetero=0;
		FST=0;

        for (un_i = 0; un_i < freq.size(); un_i++)
            freq[un_i] = 0;
		for (un_i = 0; un_i < expected_hetero.size(); un_i++)
			expected_hetero[un_i] = 0;

		for (i = 0; i < twoN; i++)
		{ // for each chromosome
            Si = pop[i].S; // S-haplotype of the chromosome

            // S-haplotype in the local deme
            // Initialization if it is the first individual of the deme
            if ( (i%two_size_patch)==0) {
                S_alleles_patch.clear();
                for (un_i = 0; un_i < freq_loc.size(); un_i++)
                    freq_loc[un_i] = 0;
            }
                for (un_j = 0; un_j < S_alleles_patch.size(); un_j++) // for each SI haplotype in the local deme
                    if (Si == S_alleles_patch[un_j]) // if it is already in the vector
                        break; // stop searching
                    if (un_j == S_alleles_patch.size()) // if it is not yet in the vector
                    S_alleles_patch.push_back(Si); // it is added
                    freq_loc[Si-1]++;
                for (un_j = 0; un_j < S_alleles.size(); un_j++) // for each SI haplotype in the population
                    if (Si == S_alleles[un_j]) // if it is already in the vector
                        break; // stop searching
                    if (un_j == S_alleles.size()) // if it is not yet in the vector
                    S_alleles.push_back(Si); // it is added
                    freq[Si-1]++; // to compute its frequency


            if ( (i%two_size_patch)==(two_size_patch-1) ) { // if last chromosome of the deme
                nS_patch_bar += S_alleles_patch.size(); // the number of SI is recorded
                nS_patch_var += S_alleles_patch.size()*S_alleles_patch.size();

                div_temp_patch = 0;

                if ((two_size_patch)>0)
				{
                    for (un_i = 0; un_i < freq_loc.size(); un_i++) // for each SI haplotype
					{
                        div_temp_patch += pow((double(freq_loc[un_i]) / double(two_size_patch)),2);// expected homozygosity in the deme at the S-locus
                    }
                    expected_hetero[compt_patch] = 1- div_temp_patch; // expected heterozygosity in the deme at the S-locus
					div_temp_patch = 1/div_temp_patch; // effective number of alleles (Nei's diversity)

                }
                D_patch += div_temp_patch; /// sum of diversities
                // variance in diversity
                D_patch_var += div_temp_patch*div_temp_patch; // to compute the variance
            }
        }

		// Adding SI haplotypes with frequency > 0.01 in the vector
		for (un_k  = 0 ; un_k<freq.size() ; un_k++)
			{
				if (freq[un_k]> (0.01*double(twoN)))
					S_freq_alleles.push_back(freq[un_k]);
			}
		nbr_haplotypes = S_freq_alleles.size();

        //// Number of SI haplotypes per deme and diversity per deme
        nS_patch_bar = nS_patch_bar/double(pv); // dividing by the number of demes
        nS_patch_var = (nS_patch_var/double(pv))-(nS_patch_bar*nS_patch_bar); // König-Huygens' formula

        D_patch = D_patch/double(pv); // dividing by the number of demes
        D_patch_var = (D_patch_var/double(pv))-(D_patch*D_patch); // König-Huygens' formula

        D = 0;
        if ((twoN)>0)
		{
            for (un_i = 0; un_i < freq.size(); un_i++) // for each frequency of S-allele
                D += pow((double(freq[un_i]) / double(twoN)),2); // Simpson index
            tot_expected_hetero = 1-D; // expected heterozygosity in the population
			D = 1/D; //Gini-Simpson index, probability that two type are different rather than being the same
        }

		// every "stepv" generations:
        if ((gen % stepv == 0) && (gen > 60000))
		{ //si on est sur une génération de mesure, cf précédement

			// Computing FST

				// Computing expected heterozygosity in each deme

				for (un_j = 0; un_j < expected_hetero.size() ; un_j++)
					sum_expected_hetero += expected_hetero[un_j];
				FST = 1 - ((sum_expected_hetero/double(expected_hetero.size()))/tot_expected_hetero);

			// writing some variables:
            /* foutresult <<"generation"
            << " " << "haplotype_number" << " " << "haplotype_diversity"
            << " " << "FST" << " " << "local_haplotype_number" << " " << "local_haplotype_diversity" << " " << "local_haplotype_number_variance" << " " << "local_haplotype_diversity_variance"
            << " " << "disp_pol" << " " << "disp_seed" <<  endl;
            foutresult.close();*/

            foutresult.open(nameFileResult,std::ofstream::app);
			foutresult << gen
                << " " << S_alleles.size() << " " << D // global alleles
                << " " << FST << " " << nS_patch_bar << " " << D_patch << " " << nS_patch_var << " " << D_patch_var // local alleles
                << " " << p_disp_bar << " " << s_disp_bar << endl; // dispersal
            foutresult.close();
		}
	}

	delete [] pop;
	delete [] temp;

	// The number of S-haplotypes is put in the ParamGeneral file
	Res.nb_haplo = nbr_haplotypes;

    return Res;
}
