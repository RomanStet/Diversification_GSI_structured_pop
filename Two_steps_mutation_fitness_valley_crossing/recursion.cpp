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

   During NbGenv generations: specificity mutations + deleterious mutations

   Other parameters are:
    numsimul: simulation number
    iv: number of iterations for the value of delta
    Nv: number of diploid individuals
	av: self-pollination rate
	deltav: inbreeding depression for SC individuals
	Usiv: mutation rate from Si to Sj
	Uscv: mutation rate from SI to SC and from SC to Si
	Une: mutation rate of enutral allele
	mUscv: allow one mutation (0) or multiple mutations (1) to segregate at the same time
	stepgv: number of generations between the generation where results are written
	dpv: pollen dispersal
	dsv: seed dispersal
	pv: number of demes
*/

result recursion(int numsimulv, int iv, int Nv, double av,
                double Uscv, int mUscv, int NbGenv, double deltav, int stepgv, double Unv,
                double dpv, double dsv, int pv, int nbr_allelesv, int nbr_alleles_initv)
{
    //i,j,k,gen:counters, chrmut: mutated chr, chrm= maternal chr, chrd= paternal chr, ind = individual
	unsigned int un_i, un_j, un_k;
	int i, j, k, gen, chrm, chrd, ind, nb, ind_chr, Si, Sii, Sj, Ni, mom, dad, SCch, SCch_patch, mut_state, modif_pos, modif_pos_dad, patch, nbr_haplotypes, male_spe, female_spe
	, male_spe_j, new_male_spe_i, new_female_spe_i, female_spe_j, num_deme, haplo_already_present, SI_anc_freq, SC_haplo_freq, iter, numb_col, mutation_type;
	SCch=0 ; SCch_patch=0;
	double w, wbar, varw, wmax, rd, fSCbar, slf, D, D_neutral, D_patch, D_patch_n, D_patch_var, D_patch_var_n, div_temp_patch, div_temp_n_patch, expected_sc, SC_haplo_freq_double, SI_anc_freq_double;
	result Res;

	bool stop_signal = false;
	bool pollen_avail = false;

	vector <int> S_alleles; S_alleles.clear(); //keeps track of the number of S-alleles in the population
	vector <int> neutral_alleles; neutral_alleles.clear(); //keeps track of the number of alleles at the neutral in the population
	vector <int> S_freq_alleles ; S_freq_alleles.clear(); //keeps track of the number of S-alleles with a certain frequency (>0.01) in the population
    vector <int> S_alleles_patch; S_alleles_patch.clear(); //keeps track of the number of S-alleles in a deme
	vector <int> neutral_alleles_patch; neutral_alleles_patch.clear(); //keeps track of the number of alleles at the neutral locus in a deme
    vector <int> SC_alleles; SC_alleles.clear(); //keeps track of the number of SC alleles in the population
    vector <int> SC_alleles_patch; SC_alleles_patch.clear(); //keeps track of the number of SC alleles in a deme
	vector <double> fSC; fSC.clear(); //vector of the frequency of SC alleles
	vector <double> dep; dep.clear(); //vector of the level of indreeding depression
	vector <double> freq; // vector of the frequency of each haplotype in the population
	vector <double> freq_loc; // vector of the frequency of each haplotype in a deme
	vector <double> freq_n; // vector of the frequency of each allele at the neutral locus in the population
	vector <double> freq_n_loc; // vector of the frequency of each allele at the neutral locus in a deme
	freq.resize(nbr_allelesv); //number of haplotypes in the population
	freq_loc.resize(nbr_allelesv); //number of haplotypes in a deme
	vector<Nall> Neut_glo; // vector containing the identity and frequency of alleles at the neutral locus in the population
	vector<Nall> Neut_loc; // vector containing the identity and frequency of alleles at the neutral locus in a deme
	vector<new_SC_haplo> created_SC_haplo; // vector that keeps track of new SC haplotypes from which they are removed when lost from the population
	new_SC_haplo haplo_temp; // temporary records the informations of a new SC haplotype
	Nall allTemp; // temporary records the identity and frequency of a new neutral allele

	// count the number of alleles at the neutral locus
	vector <double> expected_hetero; // expected heterozygosity in each deme to compute FST at the S-locus
	expected_hetero.resize(pv);
	vector <double> expected_hetero_n; // expected heterozygosity in each deme to compute FST at the neutral locus
	expected_hetero_n.resize(pv);

	// measuring the fitness valley crossing process

	vector < vector< int> > tab_caract; // table of vectors recording the haplotype, the generation and the deme for each event
	// creating pv columns begining by 1
	numb_col = 6;
	for (i=0; i< numb_col+1; i++)
	{
		tab_caract.push_back(vector<int>(1,1));
	}

	vector < vector< double> > tab_meca; // table of vectors recording the frequency of SC haplotypes and ancestral SI haplotypes
	// creating as many columns begining by one as the number of demes
	for (i=0; i< pv*2; i++)
	{
		tab_meca.push_back(vector<double>(1,1));
	}

	//deme size
	int size_patch = Nv/pv; // Nv has been previously resized so that it gives an int when divided by pv
	int two_size_patch = 2*size_patch; // Nv has been previously resized so that it gives an int when divided by pv
	int twoN = 2*Nv; // Nv has been previously resized so that it gives an int when divided by pv

    // counters for the results
	int cmptSelf = 0; //selfing counter
	int compt_patch =0; //deme counter

	// to compute the diversity in the demes
	double nS_patch_bar; // for the mean number of SI per deme
	double nS_patch_var; // for the variance in the number of SI per deme
    double nSC_patch_bar; // for the mean number of SC per deme
	double nSC_patch_var; // for the variance in the number of SC per deme
    double SC_patch_bar; // for the mean frequency of SC per deme
	double SC_patch_var; // for the variance in the frequency of SC per deme
	double tot_expected_hetero; // expected heterozygosity at the S-locus in the population
	double sum_expected_hetero; // sum of expected heterozygosities at the S-locus in all demes
	double tot_expected_hetero_n; // expected heterozygosity at the neutral locus in the population
	double sum_expected_hetero_n; // sum of expected heterozygosities at the neutral locus in all demes
	double FST; // FST for S-haplotypes
	double FST_neutral; // FST for the neutral locus
	double p_disp_bar; // to compute the actual pollen dispersal rate
    double s_disp_bar; //to compute the actual seed dispersal rate

	// population: table of 2Nn chromosomes:
	chr * pop = new chr [twoN]; // vector of the population
   	chr * temp = new chr [twoN]; // temporary vector of the population

	double * Wij = new double [Nv]; // table of fitness values (for all individuals)
	double * selfing = new double [Nv]; // table of selfing values for all individuals (= 1 if produced by selfing and = 0 si produced by outcrossing)
	for (i = 0; i < Nv; i++) //initializing the table
            selfing[i] = 0;
	double * selfing_temp = new double [Nv]; // temporary table of selfing values
	for (i = 0; i < Nv; i++) //initializing the table
            selfing_temp[i] = 0;
	// for each S-allele, sum of fitnesses of compatible chromosomes
	vector<double> Pfer; //potential of ferilization
	Pfer.resize(nbr_allelesv+1); // for each haplotype (+1 because position 0 is the sum of fitnesses for each haplotype)
    double tab_Pfer[pv][nbr_allelesv+1]; //potential of ferilization per deme
    double wbar_p[pv]; //table of the mean fitness per deme

	// creating name of output file (which indicates the parameter values):
    char nameFileParam[256];
	stringstream nameFP;
    nameFP << numsimulv << "_folder/" << numsimulv << "_param_" << iv <<".txt";
    nameFP >> nameFileParam;
	ofstream foutparam(nameFileParam);
    // writing the parameters of the simulation
	foutparam  << "no"<<" " << "i" <<" "<< "N"<< " " << "a"<< " " <<
                     "nmaxSI"<< " " << "SI_init" << " " << "Usc"<< " " << "mUsc"<< " "<< "nb_gen"<< " " << "delta"<< " " <<
                     "stepg"<< " " << "Un" << " " << "dp"<< " " << "ds"<< " " << "p"<< endl;
    foutparam << numsimulv << " " << iv << " " << Nv << " " << av <<  " " <<
                     nbr_allelesv << " " << nbr_alleles_initv << " " << Uscv << " " << mUscv << " " << NbGenv << " " << deltav << " " <<
                     stepgv << " " << Unv << " " << dpv << " " << dsv << " " << pv << endl;
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
        << " " << "FST" << " " << "FST_neutral" << " " << "local_haplotype_number" << " " << "local_haplotype_diversity" << " " << "local_neutral_diversity" << " " << "local_haplotype_number_variance" << " " << "local_haplotype_diversity_variance"<< " " << "local_neutral_diversity_variance"
        << " " << "SC_frequency" << " " << "SC_expected" << " " << "selfing_rate" << " " << "SC_mean" << " " << "SC_variance" << " " << "SC_haplotype_number" << " " << "SC_haplotype_diversity" << " " << "SC_haplotype_variance"
        << " " << "fitness_mean" << " " << "fitness_variance"
        << " " << "disp_pol" << " " << "disp_graine" <<  endl;
    foutresult.close();

	// creating the file containing the events from the fitness valley crossing process
	char nameFileMeca[256];
	stringstream nameMeca;
	nameMeca << numsimulv << "_folder/" << numsimulv << "_tabmeca_" << iv <<".txt";
	nameMeca >> nameFileMeca;
	ofstream foutMeca(nameFileMeca);
	//headers
	foutMeca << "event" << " " << "mutationnal_path" << " " << "female_spe" << " " << "male_spe" << " " << "generation" << " " << "deme" ;
    for (i=0; i<(2*pv); i++)
		{
			if (i<pv)
            {
                foutMeca << " " << "SC_haplo_freq_deme_" << i+1 ;
            }
            else
            {
                foutMeca << " " << "SI_anc_freq_deme_" << (i+1)-pv ;
            }
		}
	foutMeca << endl;
	foutMeca.close();

    // creating a file allowing to take a picture of the population
	char nameFilePicture[256];
	stringstream namePicture;
	namePicture << numsimulv << "_folder/" << numsimulv << "_picture_" << iv <<".txt";
	namePicture >> nameFilePicture;
	ofstream foutPicture(nameFilePicture);
	// headers
	foutPicture <<"phase" << " " << "deme" << " " << "all_neutral_1" << " " << "all_neutral_2" << " " << "S_pol_1" << " " << "S_pis_1" << " " << "S_pol_2" << " " << "S_pis_2" << " " << "fitness" << " " << "selfing" << endl;
    foutPicture.close();

    //////////// Population initialization
	//// All individuals are heterozygous
    ///////////

    for (i = 0; i < Nv; i++) { // for each individual in the population
        nb = 2*i;// postion of the first chromosome of the individual
        pop[nb].S_pol = rnd.randInt(nbr_alleles_initv-1)+1; // drawing a haplotype among all possible haplotypes
        pop[nb].S_pis = pop[nb].S_pol; // same at pollen and pistil to have only S-haplotypes
        do {
                pop[nb + 1].S_pol = rnd.randInt(nbr_alleles_initv-1)+1;
        } while (pop[nb + 1].S_pol == pop[nb].S_pol); //different S-haplotype at the second chromosome
        pop[nb + 1].S_pis = pop[nb + 1].S_pol; // same at pollen and pistil to have only S-haplotypes

		pop[nb].neutral_all = 0; // no variation at the neutral locus at the begining
		pop[nb + 1].neutral_all = 0; // no variation at the neutral locus at the begining

    }

    foutPicture.open(nameFilePicture,std::ofstream::app);
    for (i = 0; i < twoN; i++)
	{
        foutPicture << "S" << " " << ((i-(i%two_size_patch))/two_size_patch)+1 << " " << pop[i].neutral_all << " " << pop[i+1].neutral_all << " " << pop[i].S_pol << " " << pop[i].S_pis << " " << pop[i+1].S_pol << " " << pop[i+1].S_pis << " " << 1 << " " << selfing[i/2] << endl;
		i++;
	}
    foutPicture.close();

    if (!stop_signal==true) //if no stop signal

    for ( gen = 0; gen<NbGenv; gen++)
	{ // for all genrations

		//////////////
		////mutations
		//////////////
		mut_state=0;
        for (i = 0; i < twoN; i++)
		{ // for each chromosome
		male_spe = pop[i].S_pol; // recording male and female specificity
		female_spe = pop[i].S_pis;
            // mutation at the S-locus
            rd = rnd.randExc();
            if (rd < Uscv) { //if there is a mutation at the S-locus
                    if ( mUscv==1 || ( (SCch==0) & (mut_state==0) ) ) { //if mutiple mutations are allowed or if there is no mutation
                        if ( rnd.rand() < 0.5 ) {
                            pop[i].S_pis = rnd.randInt(nbr_allelesv-1)+1; // mutation on the pistil part
                            mut_state = 1;
							if ((pop[i].S_pis != female_spe) && (pop[i].S_pis != male_spe) ) // if the mutation produces a different specificity that is also different from the male specificity
							{
								num_deme = ((i-(i%two_size_patch))/two_size_patch) ; // deme of the SC haplotype
								haplo_already_present = 0;
								for (j = 0; j < twoN; j++)  // for each chromosome in the population
								{
									if ((pop[i].S_pis == pop[j].S_pis) && (male_spe == pop[j].S_pol)) // if the same SC haplotype is found
									{
										haplo_already_present ++;
									}
								}

								if (haplo_already_present <= 1)	// if the SC haplotype is not in the population
								{
									new_female_spe_i = pop[i].S_pis;; // the new male specificity is recorded
									compt_patch = 1;
									tab_caract[0].push_back(0); // the event is recorded as a first mutation (coded with a 0)
									tab_caract[1].push_back(1) ; // we record that it is a female-part mutation (coded with a 1)
									tab_caract[2].push_back(new_female_spe_i) ; // the female specificity is recorded
									tab_caract[3].push_back(male_spe) ; // the male specificity is recorded
									tab_caract[4].push_back(gen) ; // the generation is recorded
									tab_caract[5].push_back(num_deme+1) ; // the deme is recorded
									SC_haplo_freq = 0;
									SI_anc_freq = 0;
									haplo_temp.presence = 1 ; // this new haplotype is present
									haplo_temp.female_spe =  new_female_spe_i ; // female specificity of this haplotype
									haplo_temp.male_spe = male_spe ; // male specificity of this haplotype
									haplo_temp.mutation_type = 1 ; // it is a female-part mutation (coded by a 1)
									created_SC_haplo.push_back(haplo_temp) ; // added to the vector

									for (j = 0; j < twoN; j++) // to compute the frequency of this SC haplotype and the frequency of the ancestral SI haplotype in each deme
									{
										male_spe_j = pop[j].S_pol;
										female_spe_j = pop[j].S_pis;

										if ((new_female_spe_i == female_spe_j) && (male_spe == male_spe_j)) // if it is the same SC haplotype
										{
											SC_haplo_freq ++;
										}

										if ((male_spe_j == male_spe) && (female_spe_j == male_spe)) // if it is the ancestral SI haplotype
										{
											SI_anc_freq ++;
										}

										if ((j%two_size_patch)==(two_size_patch-1)) // if it is the last chromosme of the deme
										{
											SC_haplo_freq_double = SC_haplo_freq/double(two_size_patch) ; // frequency of the SC haplotype in the deme
											SI_anc_freq_double = SI_anc_freq/double(two_size_patch) ; // frequency of the ancestral SI haplotype in the deme
											tab_meca[compt_patch-1].push_back(SC_haplo_freq_double) ; // frequency of the SC haplotype in the table
											tab_meca[(compt_patch+pv)-1].push_back(SI_anc_freq_double) ; // frequency of the ancestral SI haplotype in the table
											SC_haplo_freq = 0;
											SI_anc_freq = 0;
											compt_patch ++;
										}
									}
								}
							}
							if ((female_spe != pop[i].S_pis) && (male_spe == pop[i].S_pis)) // if compensatory mutation: different male and female specificites and the new female specificity is the same as the male specificity
							{
								num_deme = ((i-(i%two_size_patch))/two_size_patch) ; // number of the deme
								compt_patch = 1;
								tab_caract[0].push_back(2); // the event is recorded as a compensatory mutation (coded with a 2)
								tab_caract[1].push_back(2) ; // we record that it is a compensatory mutation on a male-part mutation (coded with a 2)
								tab_caract[2].push_back(female_spe) ; // the female specificity is recorded
								tab_caract[3].push_back(male_spe) ; // the male specificity is recorded
								tab_caract[4].push_back(gen) ; // the generation is recorded
								tab_caract[5].push_back(num_deme+1) ; // the deme is recorded
								SC_haplo_freq = 0;
								SI_anc_freq = 0;

								for (j = 0; j < twoN; j++) // to compute the frequency of this SC haplotype and the frequency of the ancestral SI haplotype in each deme
								{
									male_spe_j = pop[j].S_pol;
									female_spe_j = pop[j].S_pis;

									if ((female_spe == female_spe_j) && (pop[i].S_pol == male_spe_j)) // if it is the same SC haplotype
									{
										SC_haplo_freq ++;
									}

									if ((female_spe_j == female_spe) && (male_spe_j == female_spe)) // if it is the ancestral SI haplotype
									{
										SI_anc_freq ++;
									}

									if ((j%two_size_patch)==(two_size_patch-1)) // if it is the last chromosme of the deme
									{
										SC_haplo_freq_double = SC_haplo_freq/double(two_size_patch) ; // frequency of the SC haplotype in the deme
										SI_anc_freq_double = SI_anc_freq/double(two_size_patch) ; // frequency of the ancestral SI haplotype in the deme
										tab_meca[compt_patch-1].push_back(SC_haplo_freq_double) ; // frequency of the SC haplotype in the table
										tab_meca[(compt_patch+pv)-1].push_back(SI_anc_freq_double) ; // frequency of the ancestral SI haplotype in the table
										SC_haplo_freq = 0;
										SI_anc_freq = 0;
										compt_patch ++;
									}
								}
							}
                        } else { // mutation on the pollen part

							pop[i].S_pol = rnd.randInt(nbr_allelesv-1)+1;
                            mut_state = 1;

							if ((pop[i].S_pol != male_spe) && (pop[i].S_pol != female_spe) ) // if the mutation produces a different specificity that is also different from the female specificity
							{
								num_deme = ((i-(i%two_size_patch))/two_size_patch) ; // number of the deme
								haplo_already_present = 0;
								for (j = 0; j < twoN; j++) // for each chromosome in the population
								{
									if ((pop[i].S_pol == pop[j].S_pol) && (female_spe == pop[j].S_pis)) // if the same SC haplotype is found
									{
										haplo_already_present ++;
									}
								}

								if (haplo_already_present <= 1)	// if the SC haplotype is not in the population
								{
									new_male_spe_i = pop[i].S_pol;  // the new male specificity is recorded
									compt_patch = 1;
									tab_caract[0].push_back(0); // the event is recorded as a first mutation (coded with a 0)
									tab_caract[1].push_back(2) ; // we record that it is a male-part mutation (coded with a 2)
									tab_caract[2].push_back(female_spe) ; // the female specificity is recorded
									tab_caract[3].push_back(new_male_spe_i) ; // the male specificity is recorded
									tab_caract[4].push_back(gen) ; // the generation is recorded
									tab_caract[5].push_back(num_deme+1) ; // the deme is recorded
									SC_haplo_freq = 0;
									SI_anc_freq = 0;
									haplo_temp.presence = 1 ;// this new haplotype is present
									haplo_temp.female_spe =  female_spe ; // female specificity of this haplotype
									haplo_temp.male_spe = new_male_spe_i ; // male specificity of this haplotype
									haplo_temp.mutation_type = 2 ; // it is a male-part mutation (coded by a 2)
									created_SC_haplo.push_back(haplo_temp) ; // added to the vector

									for (j = 0; j < twoN; j++) // to compute the frequency of this SC haplotype and the frequency of the ancestral SI haplotype in each deme
									{
										male_spe_j = pop[j].S_pol;
										female_spe_j = pop[j].S_pis;

										if ((new_male_spe_i == male_spe_j) && (female_spe == female_spe_j)) // if it is the same SC haplotype
										{
											SC_haplo_freq ++;
										}

										if ((female_spe_j == female_spe) && (male_spe_j == female_spe)) // if it is the ancestral SI haplotype
										{
											SI_anc_freq ++;
										}

										if ((j%two_size_patch)==(two_size_patch-1)) // if it is the last chromosme of the deme
										{
											SC_haplo_freq_double = SC_haplo_freq/double(two_size_patch) ; // frequency of the SC haplotype in the deme
											SI_anc_freq_double = SI_anc_freq/double(two_size_patch) ; // frequency of the ancestral SI haplotype in the deme
											tab_meca[compt_patch-1].push_back(SC_haplo_freq_double) ; // frequency of the SC haplotype in the table
											tab_meca[(compt_patch+pv)-1].push_back(SI_anc_freq_double) ; // frequency of the ancestral SI haplotype in the table
											SC_haplo_freq = 0;
											SI_anc_freq = 0;
											compt_patch ++;
										}
									}
								}
							}

							if ((male_spe != pop[i].S_pol) && (female_spe == pop[i].S_pol)) // if compensatory mutation: different male and female specificites and the new male specificity is the same as the female specificity
							{
								num_deme = ((i-(i%two_size_patch))/two_size_patch) ; // number of the deme
								compt_patch = 1;
								tab_caract[0].push_back(2); // the event is recorded as a compensatory mutation (coded with a 2)
								tab_caract[1].push_back(1) ;  // we record that it is a compensatory mutation on a female-part mutation (coded with a 1)
								tab_caract[2].push_back(female_spe) ; // the female specificity is recorded
								tab_caract[3].push_back(male_spe) ; // the male specificity is recorded
								tab_caract[4].push_back(gen) ; // the generation is recorded
								tab_caract[5].push_back(num_deme+1) ; // the deme is recorded
								SC_haplo_freq = 0;
								SI_anc_freq = 0;

								for (j = 0; j < twoN; j++) // to compute the frequency of this SC haplotype and the frequency of the ancestral SI haplotype in each deme
								{
									male_spe_j = pop[j].S_pol;
									female_spe_j = pop[j].S_pis;

									if ((female_spe == female_spe_j) && (male_spe == male_spe_j)) // if it is the same SC haplotype
									{
										SC_haplo_freq ++;
									}

									if ((female_spe_j == male_spe) && (male_spe_j == male_spe)) // if it is the ancestral SI haplotype
									{
										SI_anc_freq ++;
									}

									if ((j%two_size_patch)==(two_size_patch-1)) // if it is the last chromosme of the deme
									{
										SC_haplo_freq_double = SC_haplo_freq/double(two_size_patch) ; // frequency of the SC haplotype in the deme
										SI_anc_freq_double = SI_anc_freq/double(two_size_patch) ; // frequency of the ancestral SI haplotype in the deme
										tab_meca[compt_patch-1].push_back(SC_haplo_freq_double) ; // frequency of the SC haplotype in the table
										tab_meca[(compt_patch+pv)-1].push_back(SI_anc_freq_double) ; // frequency of the ancestral SI haplotype in the table
										SC_haplo_freq = 0;
										SI_anc_freq = 0;
										compt_patch ++;
									}
								}
							}
                        }
                    }
                }
			rd = rnd.randExc();
			if (rd < Unv) // if mutation at the neutral locus
			{
				pop[i].neutral_all = rnd.randInt(1000000);
			}
        }

		//////////////
		//// Creating next generation
		//////////////

        /* Computing fitness */
		//fills table Wij with fitnesses of all individuals,
		// table wmax with maximal fitness,
		// table Pfer with sums of fitnesses of compatible chromosomes:

			wbar = 0;
			wmax = 0;
			varw = 0;

			for (k=0; k < pv ; k++)
				wbar_p[k]=0;

			for (un_k = 0; un_k < Pfer.size(); un_k++)
				Pfer[un_k] = 0;

			for (int l=0; l < pv; l++) {
				for (un_k = 0; un_k < Pfer.size();un_k++) {
					tab_Pfer[l][un_k]=0;
				}
			}

			for (j = 0; j < Nv; j++)
			{ // for each individual
				nb = 2 * j; //position at first chromosome
				w = 1.0-selfing[j]*deltav; // fitness
				Wij[j] = w; // put in the table
				varw += w * w; // to compute the variance in fitness

				if (wmax < w) // if the fitness is higher than the maximum fitness
					wmax = w; // it becomes the maximum fitness

				// S-alleles:
				for (k = 0; k < 2; k++) { // for each chromosome of this individual
					Si = pop[nb + k].S_pol; // male specificity
					Pfer[Si] += w; // fitness of this haplotype
					Pfer[0] += w; // sum of fitnesses
					tab_Pfer[(j-j%size_patch)/size_patch][Si] += w; // put in the table
					tab_Pfer[(j-j%size_patch)/size_patch][0] += w; // put in the table
				}
			}
			varw /= Nv; // to compute the variance in fitness

			varw /= Nv; // calcul préliminaire de variance de la fitness
			// mean fitness
			wbar += Pfer[0]; // sum of mean fitnesses
			for (int l=0; l < pv; l++) {
				for (un_k = 1; un_k < Pfer.size();un_k++) {
					wbar_p[l] += tab_Pfer[l][un_k];
				}
			}
			wbar /= twoN; // mean
			varw -= wbar * wbar; // variance


        cmptSelf = 0;
        expected_sc = 0;
        p_disp_bar = 0;
        s_disp_bar = 0;

        // sampling the next generation:
        for (ind = 0; ind < Nv; ind++)
		{ // for each individual
            ind_chr = 2*ind; // to work with chromosomes
            do
			{ // sampling the mother according to her fitness
					if (rnd.randExc() >= dsv) { // if no seed dispersal
						chrm = (ind_chr-ind_chr%two_size_patch)+rnd.randInt(two_size_patch-1); // a chromosome is drawn in the same deme
						mom = chrm/2; // the individual carying this chromosome
					} else { // if seed dispersal, the mother is drawn from a different deme
						do{  // sampling the mother according to her fitness
							chrm = rnd.randInt(twoN-1); // a chromosome is drawn in the population
							mom = chrm/2; // the individual carying this chromosome
						} while ((ind_chr-ind_chr%two_size_patch==chrm-chrm%two_size_patch) || (rnd.rand() >= (Wij[mom] / wmax) ));  // while the relative fitness of the individual is not sufficient or individual from the same deme
						s_disp_bar++;
					}

				// to check if the mother can be fertilized
				if (chrm % 2 == 0) { // if the chromosome is the first if this mother
					modif_pos = 1; // complementary chromosome
				} else { // if the chromosome is the second if this mother
					modif_pos = -1; // complementary chromosome
				}
				Si = pop[chrm].S_pis; // records the female specificity of one chromosome
				Sj = pop[chrm + modif_pos].S_pis; // records the female specificity of the other chromosome
				if (Si!=Sj) { // if both female specificities are different
					pollen_avail = (Pfer[0]-Pfer[Si]-Pfer[Sj])>(Pfer[0]/(10*twoN)); // computes the proportion of compatible pollen in the pollen pool
				} else { // if both female specificities are the same
					pollen_avail = (Pfer[0]-Pfer[Si])>(Pfer[0]/(10*twoN)); // if the proportion of compatible pollen is too low we can neglect it
				}
            } while ( (rnd.rand() >= (Wij[mom] / wmax) ) || (!pollen_avail) ); // while the relative fitness of the individual is not sufficient or if there is no compatible pollen

            // mother's deme
            patch = (mom-mom%size_patch)/size_patch;

            //self-incompatibility
            if ( ( (pop[chrm].S_pol!=Si) && (pop[chrm].S_pol!=Sj) ) || ( (pop[chrm+modif_pos].S_pol!=Si) & (pop[chrm+modif_pos].S_pol!=Sj) ) ) {
                // if the male specificity is different from both female specificity = self-compatible
                //compute the actual selfing rate according to compatible pollen in the pollen pool
                if ( ( (pop[chrm].S_pol!=Si) && (pop[chrm].S_pol!=Sj) ) && ( (pop[chrm+modif_pos].S_pol!=Si) & (pop[chrm+modif_pos].S_pol!=Sj) ) )
				{ // if both female specificities are compatible
                    slf = av*Wij[mom]*2/( av*Wij[mom]*2+ (1-av) * ((1-dpv)*((tab_Pfer[patch][0]-tab_Pfer[patch][Si]-tab_Pfer[patch][Sj]-Wij[mom]*2)/(size_patch-1))+dpv*((Pfer[0]-Pfer[Si]-Pfer[Sj]-Wij[mom]*2)/(Nv-1)) ) );
                }
				else
				{ // if one female specificity is compatible
                    slf = av*Wij[mom]/( av*Wij[mom]+ (1-av) * ((1-dpv)*((tab_Pfer[patch][0]-tab_Pfer[patch][Si]-tab_Pfer[patch][Sj]-Wij[mom])/(size_patch-1))+dpv*((Pfer[0]-Pfer[Si]-Pfer[Sj]-Wij[mom])/(Nv-1)) ) );
                }
            } else { // none of female specificites are compatible
                slf = 0; // no selfing
            }
            expected_sc += slf;

            // sampling the father:
            //selfing:
            if (rnd.randExc() < slf)
			{ // if selfing
                cmptSelf++;
				selfing[ind] = 1.0;
                if ( ( (pop[chrm].S_pol!=Si) && (pop[chrm].S_pol!=Sj) ) && ( (pop[chrm+modif_pos].S_pol!=Si) && (pop[chrm+modif_pos].S_pol!=Sj) ) )
				{ // if both haplotypes are SC
                    // drawing the first chromosome
                    if (rnd.rand() < 0.5)
                        temp[ind_chr] = pop[chrm]; // if it is the maternal chromosome
                    else
                        temp[ind_chr] = pop[chrm+modif_pos]; // if it is the other chromosome

					// drawing the second chromosome
					if (rnd.rand() < 0.5)
                        temp[ind_chr+1] = pop[chrm]; // if it is the maternal chromosome
                    else
                        temp[ind_chr+1] = pop[chrm+modif_pos]; // if it is the other chromosome
                }
				else
				{ // one haplotype is SC
                    if ( (pop[chrm].S_pol!=Si) && (pop[chrm].S_pol!=Sj) ) // if it is the maternal chromosome which is compatible
					{
						temp[ind_chr] = pop[chrm]; // it is chosen as the first chromosome of the offspring
						// drawing the second chromosome of the offspring
						if (rnd.rand() < 0.5)
                        temp[ind_chr+1] = pop[chrm]; // if it is the maternal chromosome
						else
                        temp[ind_chr+1] = pop[chrm+modif_pos]; // if it is the other chromosome
					}
					else // if it is the other chromosome which is compatible
					{
                        temp[ind_chr] = pop[chrm+modif_pos]; // it is chosen as the first chromosome of the offspring
						// drawing the second chromosome of the offspring
						if (rnd.rand() < 0.5)
                        temp[ind_chr+1] = pop[chrm]; // if it is the maternal chromosome
						else
                        temp[ind_chr+1] = pop[chrm+modif_pos]; // if it is the other chromosome
                    }
                }
            }
			else
			{ //outcrossing
                do { // searching for a paternal haplotype with specificites compatible with those of the marternal haplotype

                    if (rnd.randExc() >= dpv) { // if no pollen dispersal
                        chrd = (chrm-chrm%two_size_patch)+rnd.randInt(two_size_patch-1); // drawing the parternal haplotype in the local deme
                        dad = chrd/2; // individual carrying this haplotype
                    } else {// if pollen dispersal
                        do {
                            chrd = rnd.randInt(twoN-1); // drawing the parternal haplotype in the population
                            dad = chrd/2; // individual carrying this haplotype
                        } while (chrm-chrm%two_size_patch==chrd-chrd%two_size_patch); // while no pollen from a different deme
                    }

                } while ( (mom == dad) || (rnd.rand() >= Wij[dad] / wmax) || (pop[chrd].S_pol == Si && pop[chrd].S_pol * Si != 0) || (pop[chrd].S_pol == Sj && pop[chrd].S_pol * Sj != 0) );
                // while the mother and father are the same individual or fitness conditions not fulfilled or incompatibility if no SC haplotype

                if (chrm-chrm%two_size_patch!=chrd-chrd%two_size_patch) { // if pollen dispersal
                    p_disp_bar++;
                }
				selfing[ind] = 0.0;

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

			// idependant transmission of the neutral locus
			if (rnd.rand() < 0.5)
				{
					temp[ind_chr].neutral_all = pop[chrm].neutral_all; // if it is one chromosome
				}
				else
				{
					temp[ind_chr].neutral_all = pop[chrm+modif_pos].neutral_all; // if it is the other chromosome
				}

			if (chrd % 2 == 0) { // if the paternal chromosome is its first chromosome
					modif_pos_dad = 1; // the other chromosome
				} else { // if the paternal chromosome is its second chromosome
					modif_pos_dad = -1; // the other chromosome
				}

			if (rnd.rand() < 0.5)
				{
					temp[ind_chr+1].neutral_all = pop[chrd].neutral_all; // if it is one chromosome
				}
				else
				{
					temp[ind_chr+1].neutral_all = pop[chrd+modif_pos_dad].neutral_all; // if it is the other chromosome
				}
        }

		// dispersal rate
        s_disp_bar /= Nv;
        p_disp_bar /= Nv;

        for (i = 0; i < twoN; i++) // for each chromosome
            pop[i] = temp[i]; // replacing the chromosomes by the new ones


		//////////////
		//// Various measures and writing of the results
		//////////////
		//// Computing of frequency of S-alleles and stoping the simulation if the diversity is not sufficient
		// Initializing
		S_alleles.clear();
		Neut_glo.clear();
		S_freq_alleles.clear();
		SC_alleles.clear();
		SCch = 0;
        nS_patch_bar=0;
		nS_patch_var=0;
        nSC_patch_bar=0;
		nSC_patch_var=0;
        SC_patch_bar=0;
        SC_patch_var=0;
        D_patch=0;
		D_patch_n=0;
        D_patch_var=0;
		D_patch_var_n=0;
		compt_patch=0;
		tot_expected_hetero=0;
		sum_expected_hetero=0;
		tot_expected_hetero_n=0;
		sum_expected_hetero_n=0;
		FST=0;
		FST_neutral=0;

        for (un_i = 0; un_i < freq.size(); un_i++)
            freq[un_i] = 0;
		for (un_i = 0; un_i < expected_hetero.size(); un_i++)
			expected_hetero[un_i] = 0;
		for (un_i = 0; un_i < expected_hetero_n.size(); un_i++)
			expected_hetero_n[un_i] = 0;
		for (un_i = 0; un_i < created_SC_haplo.size(); un_i++)
			created_SC_haplo[un_i].presence = 0;

		for (i = 0; i < twoN; i++)
		{ // for each chromosome
            Si = pop[i].S_pol; // male specificity
            Sii = pop[i].S_pis; // female specificity
			Ni = pop[i].neutral_all; // allele at neutral locus

            // S-locus in the local deme
            // Initialization if it is the first individual of the deme
            if ( (i%two_size_patch)==0) {
                S_alleles_patch.clear();
                SC_alleles_patch.clear();
				Neut_loc.clear();
                SCch_patch=0;
                for (un_i = 0; un_i < freq_loc.size(); un_i++)
                    freq_loc[un_i] = 0;
            }

            if (Si!=Sii) { // if it is a SC haplotype

				for (un_j=0; un_j<created_SC_haplo.size(); un_j++) // for each created SC haplotye
				{
					if ((Si == created_SC_haplo[un_j].male_spe) && (Sii == created_SC_haplo[un_j].female_spe)) // if it is already created
					{
						created_SC_haplo[un_j].presence = 1 ;
					}
				}

				SCch_patch++;
                for (un_j = 0; un_j < SC_alleles_patch.size(); un_j++) // for each SC haplotype in the local deme
                    if (Si == SC_alleles_patch[un_j]) // if it is already in the vector
                        break; // stop searching
                    if (un_j == SC_alleles_patch.size()) // if it is not yet in the vector
                    SC_alleles_patch.push_back(Si); // it is added

                SCch++;
                for (un_j = 0; un_j < SC_alleles.size(); un_j++) // for each SC haplotype in the population
                    if (Si == SC_alleles[un_j]) // if it is already in the vector
                        break; // stop searching
                if (un_j == SC_alleles.size()) // if it is not yet in the vector
                    SC_alleles.push_back(Si); // it is added

            } else { // if it is a SI haplotype
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
            }

			// Alleles at neutral locus in the local deme

			for (un_j = 0; un_j < Neut_loc.size(); un_j++) // for each neutral locus in the local deme
			{
                if (Ni == Neut_loc[un_j].all) // if it is already in the vector
				{
					Neut_loc[un_j].freq++; // we add 1 to its count
					break; // stop searching
				}
			}
            if (un_j == Neut_loc.size()) // if it is not yet in the vector
			{
				allTemp.all = Ni; // its identity is recorded
                allTemp.freq = 1; // 1 is added to its frequency
                Neut_loc.push_back(allTemp); // added to the vector
			}

            // Alleles at neutral locus in the population
			for (un_j = 0; un_j < Neut_glo.size(); un_j++) // for each neutral locus in the population
			{
                if (Ni == Neut_glo[un_j].all) // if it is already in the vector
				{
					Neut_glo[un_j].freq++; // we add 1 to its count
					break; // stop searching
				}
			}
            if (un_j == Neut_glo.size()) // if it is not yet in the vector
			{
				allTemp.all = Ni; // its identity is recorded
                allTemp.freq = 1; // 1 is added to its frequency
                Neut_glo.push_back(allTemp); // added to the vector
			}

            if ( (i%two_size_patch)==(two_size_patch-1) ) { // if last chromosome of the deme
                nS_patch_bar += S_alleles_patch.size(); // number of SI is recorded
                nS_patch_var += S_alleles_patch.size()*S_alleles_patch.size();

                nSC_patch_bar += SC_alleles_patch.size(); // number of SC is recorded
                nSC_patch_var += SC_alleles_patch.size()*SC_alleles_patch.size();

                SC_patch_bar += SCch_patch/double(two_size_patch); // proportion of SC
                SC_patch_var += (SCch_patch/double(two_size_patch))*(SCch_patch/double(two_size_patch));

                div_temp_patch = 0;
				div_temp_n_patch = 0;

                if ((two_size_patch-SCch_patch)>0)
				{
                    for (un_i = 0; un_i < freq_loc.size(); un_i++) // for each SI haplotype
					{
                        div_temp_patch += pow((double(freq_loc[un_i]) / double(two_size_patch-SCch_patch)),2); // expected homozygosity in the deme at the S-locus
                    }
                    expected_hetero[compt_patch] = 1- div_temp_patch; // expected heterozygosity in the deme at the S-locus
					div_temp_patch = 1/div_temp_patch; // effective number of alleles (Nei's diversity)
                }
                D_patch += div_temp_patch; // sum of diversities
                // variance in diversity
                D_patch_var += div_temp_patch*div_temp_patch; // to compute the variance

				for (un_i = 0; un_i < Neut_loc.size(); un_i++) // for each neutral allele in the local deme
				{
					div_temp_n_patch += pow((double(Neut_loc[un_i].freq) / double(two_size_patch)),2); // expected homozygosity in the deme at the neutral locus
				}

                expected_hetero_n[compt_patch] = 1- div_temp_n_patch; // expected heterozygosity in the deme at the neutral locus
				div_temp_n_patch = 1/div_temp_n_patch; // effective number of alleles (Nei's diversity)

				D_patch_n += div_temp_n_patch; // sum of diversities
				D_patch_var_n += div_temp_n_patch*div_temp_n_patch;
				compt_patch ++;
            }
        }

		// Checking if any loss of SC haplotype

		for (un_i = 0; un_i < created_SC_haplo.size(); un_i++) // for all new SC haplotypes
		{
			if (created_SC_haplo[un_i].presence == 0) // if it is not found in the population
			{
				compt_patch = 1;
				tab_caract[0].push_back(1) ; // the event is recorded as a loss of SC haplotype (coded with a 1)
				tab_caract[1].push_back(0) ; // no associated ancestral SI haplotype (coded with a 0)
				tab_caract[2].push_back(created_SC_haplo[un_i].female_spe) ; // the female specificity is recorded
				tab_caract[3].push_back(created_SC_haplo[un_i].male_spe) ; // the male specificity is recorded
				tab_caract[4].push_back(gen) ; // the generation is recorded
				tab_caract[5].push_back(0) ; // non associated deme (coded with a 0)
				SC_haplo_freq = 0;
				SI_anc_freq = 0;

				for (j = 0; j < twoN; j++) // to compute the frequency of this SC haplotype and the frequency of the ancestral SI haplotype in each deme
				{
					male_spe_j = pop[j].S_pol;
					female_spe_j = pop[j].S_pis;
					female_spe = created_SC_haplo[un_i].female_spe;
					male_spe = created_SC_haplo[un_i].male_spe;
					mutation_type = created_SC_haplo[un_i].mutation_type ;

					if (mutation_type == 1) // si c'est un haplotype mutant SC pistil
					{
						if ((male_spe_j == male_spe) & (female_spe_j == male_spe)) // if it is the same SC haplotype
						{
							SI_anc_freq ++; // on ajoute 1 au compteur de l'haplotype SI ancestral
						}
					}
					else if (mutation_type == 2)
					{
						if ((male_spe_j == female_spe) & (female_spe_j == female_spe)) // if it is the ancestral SI haplotype
						{
							SI_anc_freq ++;
						}
					}

					if ((j%two_size_patch)==(two_size_patch-1)) // if it is the last chromosme of the deme
					{
						SI_anc_freq_double = SI_anc_freq/double(two_size_patch) ; // frequency of the ancestral SI haplotype in the deme
						tab_meca[compt_patch-1].push_back(0) ; // the frequency of the SC haplotype in null
						tab_meca[(compt_patch+pv)-1].push_back(SI_anc_freq_double) ;  // frequency of the ancestral SI haplotype in the table
						SI_anc_freq = 0;
						compt_patch ++;
					}
				}
				created_SC_haplo.erase(created_SC_haplo.begin() + un_i); // the haplotype is removed from the vector
			}
		}

		// Adding SI haplotypes with frequency > 0.01 in the vector
		for (un_k  = 0 ; un_k<freq.size() ; un_k++)
			{
				if (freq[un_k]> (0.01*double(twoN)))
					S_freq_alleles.push_back(freq[un_k]);
			}
		nbr_haplotypes = S_freq_alleles.size();

        //// Number of SI, SC haplotypes, proportion of SC per deme and diversity per deme
        nS_patch_bar = nS_patch_bar/double(pv); // dividing by the number of demes
        nS_patch_var = (nS_patch_var/double(pv))-(nS_patch_bar*nS_patch_bar); // König-Huygens' formula

        nSC_patch_bar = nSC_patch_bar/double(pv); // dividing by the number of demes
        nSC_patch_var = (nSC_patch_var/double(pv))-(nSC_patch_bar*nSC_patch_bar); // König-Huygens' formula

        SC_patch_bar = SC_patch_bar/double(pv); // dividing by the number of demes
        SC_patch_var = (SC_patch_var/double(pv))-(SC_patch_bar*SC_patch_bar); // König-Huygens' formula

        D_patch = D_patch/double(pv); // dividing by the number of demes
        D_patch_var = (D_patch_var/double(pv))-(D_patch*D_patch); // König-Huygens' formula

		D_patch_n = D_patch_n/double(pv); // dividing by the number of demes
		D_patch_var_n = (D_patch_var_n/double(pv))-(D_patch_n*D_patch_n); // König-Huygens' formula

        D = 0;
		D_neutral = 0;
        if ((twoN-SCch)>0)
		{
            for (un_i = 0; un_i < freq.size(); un_i++) // for each frquency of neutral allele
                D += pow((double(freq[un_i]) / double(twoN-SCch)),2); // Simpson index
            tot_expected_hetero = 1-D; // expected heterozygosity in the population
			D = 1/D; //Gini-Simpson index, probability that two type are different rather than being the same
        }

		for (un_i = 0; un_i < Neut_glo.size(); un_i++)
		{
			D_neutral += pow((double(Neut_glo[un_i].freq) /double(twoN)),2);
		}
		tot_expected_hetero_n = 1 - D_neutral;
		D_neutral = 1/D_neutral;

		// measuring inbreeding depression every "stepg" generations:
        if (gen % stepgv == 0)
		{
			// Computing FST

				// Computing expected heterozygosity in each deme

				for (un_j = 0; un_j < expected_hetero.size() ; un_j++)
					sum_expected_hetero += expected_hetero[un_j];
				FST = 1 - ((sum_expected_hetero/double(expected_hetero.size()))/tot_expected_hetero);

				for (un_j = 0; un_j < expected_hetero_n.size() ; un_j++)
					sum_expected_hetero_n += expected_hetero_n[un_j];
				FST_neutral = 1 - ((sum_expected_hetero_n/double(expected_hetero_n.size()))/tot_expected_hetero_n);

			// writing some variables:
            /* foutresult <<"generation"
            << " " << "haplotype_number" << " " << "haplotype_diversity" << " " << "neutral_diversity"
            << " " << "FST" << " " << "FST_neutral" << " " << "local_haplotype_number" << " " << "local_haplotype_diversity" << " " << "local_neutral_diversity" << " " << "local_haplotype_number_variance" << " " << "local_haplotype_diversity_variance" << " " << "local_neutral_diversity_variance"
            << " " << "SC_frequency" << " " << "SC_expected" << " " << "selfing_rate" << " " << "SC_mean" << " " << "SC_variance" << " " << "SC_haplotype_number" << " " << "SC_haplotype_diversity" << " " << "SC_haplotype_variance"
            << " " << "fitness_mean" << " " << "fitness_variance"
            << " " << "pollen_dispersal" << " " << "seed_dispersal" <<  endl;
            foutresult.close();*/

            foutresult.open(nameFileResult,std::ofstream::app);
			foutresult << gen
                << " " << S_alleles.size() << " " << D << " " << D_neutral // global alleles
                << " " << FST << " " << FST_neutral << " " << nS_patch_bar << " " << D_patch << " " << D_patch_n << " " << nS_patch_var << " " << D_patch_var << " " << D_patch_var_n // local alleles
                << " " << double (SCch) / twoN << " " << expected_sc/Nv << " " << double (cmptSelf) / Nv << " " << SC_patch_bar << " " << SC_patch_var << " " << nSC_patch_bar << " " << "NA" << " " << nSC_patch_var   // SC haplotypes
                << " " << wbar << " " << varw // fitness
                << " " << p_disp_bar << " " << s_disp_bar << endl; // dispersal
            foutresult.close();
            if (gen >= NbGenv - 100 * stepgv) // if in the last generations
                fSC.push_back(double (SCch) / twoN); //recording the frequency of SC haplotypes
		}

		if (gen < 500001)
		{
			// Writing the table recording the events concerning the crossing of fitness valleys

			foutMeca.open(nameFileMeca,std::ofstream::app);
			iter = tab_caract[0].size()-1; // the number of lines to write in the table excluding the first one

			for (j=1; j < iter+1; j++)
			{
				foutMeca << tab_caract [0] [j] << " " << tab_caract [1] [j] << " " << tab_caract [2] [j] << " " << tab_caract [3] [j] << " " << tab_caract [4] [j] << " " << tab_caract [5] [j] << " ";
				for (i=0; i<(2*pv); i++)
				{
					foutMeca << tab_meca[i] [j] << " " ;
				}
				foutMeca << endl;
			}
			foutMeca.close();

			// Initializing tables

			for (j=1; j < iter+1; j++)
			{
				tab_caract[0].pop_back() ;
				tab_caract[1].pop_back() ;
				tab_caract[2].pop_back() ;
				tab_caract[3].pop_back() ;
				tab_caract[4].pop_back() ;
				tab_caract[5].pop_back() ;

				for (i=0; i<(2*pv); i++)
				{
					tab_meca[i].pop_back() ;

				}
			}
		}
	}

	delete [] pop;
	delete [] temp;
	delete [] Wij;
	fSCbar = 0;

    // Computing the mean frequency of SC haplotypes
	for (un_i = 0; un_i < fSC.size(); un_i++)
    fSCbar += fSC[un_i];
    fSCbar /= fSC.size();
	Res.frSC = fSCbar;
	// The number of S-haplotypes is put in the ParamGeneral file
	Res.nb_haplo = nbr_haplotypes;

    return Res;
}
