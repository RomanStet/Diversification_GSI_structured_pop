# Diversification_GSI_structured_pop

## Simulation programs Stetsenko et al., 2023

The programs simulate *N* diploid individuals with gametophytic self-incompatibility. Each individual has 2 copies of a chromosome each carrying an *S*-allele (for the program "One_step_mutation") or an *S*-haplotype (for the program "Two_steps_mutation" and "Two_steps_mutation_fitness_valley_crossing"). The population is divided into *p* equally-sized demes. At each generation, pollen disperses between demes with probability *d<sub>p</sub>*.

Each program produces a file "results.txt" which stores parameter values and the time length of the simulation, and another file (result_\'85txt) where each line gives the mean map length, mean fitness, mean number of deleterious alleles per chromosome and number of fixed deleterious alleles, measured at different generations.

The program uses the Mersenne Twister random number generator, and the file MersenneTwister.h (available from http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/MersenneTwister.h) must be included in the program before compilation. Note that depending on your compiler, you may have to replace istream and ostream with std::istream and std:ostream in MersenneTwister.h.

The file depression.h provides function prototypes and global variables, in particular the structure "chr" that represents a chromosome with its *S*-haplotype. 
The main function (in "main.cpp") calls functions (defined in "files.cpp") that read parameters from an input file ("parameters.txt"), write parameters in the file "results.txt" and then calls the function "recursion" that performs the simulation. 
The file "parameters.txt" provides parameter values in the order indicated at the beginning of the file; the line corresponding to the parameter set must begin with *.
Each execution of the program will scan through the "parameters.txt" file from top to bottom and will perform the simulation *i* times (*i* being the parameter defining the number of iterations) for the first line begining with a "*".
While runining the "*" is replaced by a "!". Once it is done, the "!" is replaced by a "# "and the program scans the "parameters.txt" file again to find the first line begining with a "*". Multiple executions can thus be lauched to run in parallele on the same "parameters.txt" file. 

Using gcc, the program can be compiled using the command
`g++ -O3 *.cpp`


### One_step_mutation

Performs the simulation with the 1-step mutation model as in Schierup et al., 1998. It was used to draw figures ...

### Two_steps_mutation

Performs the simulation with the 2-step mutation model as in Gervais et al., 2011. It was used to draw figures ...

### Two_steps_mutation_fitness_valley_crossing

Same as "Two_steps_mutation" but records 
