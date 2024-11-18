// This is a program to calculate the 3D HP model using wl algorithm with some modification
// note : yang diubah adalah jumlah monomer N, urutan asam amino, rentang energi

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <limits.h>

//declaration of functions defined in this program
int get_binactive(double dos[], int range);
void norm_dos(double dos[], int range, double normdos[]);
int check_flatness(long long int his[],double threshold, int range);
void dopivotmove(int simdir[], int tmpsimdir[], int index, int pivottype, int flag[]);

const int N = 46; // number of monomers

int main() {
    
    int amino_sequence[] = {2,2,0,0,2,2,0,0,2,1,2,1,0,1,0,0,1,0,2,2,2,2,1,2,0,0,2,2,2,2,2,0,0,0,0,2,2,2,2,0,2,2,1,2,2,1};; // amino acid sequence of the protein. H = 0 and P =1
    // int dim = 3; // dimension of the protein lattice    

    // number of internal HH bond in the protein chain
    int HHbond = 0;
    for (int i = 0; i < N - 1; i++) {
        if(amino_sequence[i] == 0 && amino_sequence[i+1] == 0) HHbond++;
    }

    // number of internal H0 bond in the protein chain
    int H0bond = 0;
    for (int i = 0; i < N - 1; i++) {
        if((amino_sequence[i] == 0 && amino_sequence[i+1] == 2) || (amino_sequence[i] == 2 && amino_sequence[i+1] == 0)) H0bond++;
    }
    
    // define array to store the direction of the arrow of monomers 
    int simdir[N-1];
    int tmpsimdir[N-1];

    // initialize the simdir array
    for (int i = 0; i < N-1; i++) {  
        simdir[i] = 100;
    }

    /// define array to store the location of monomer in site and the monomer that occupy it 
    int site[N];
    int tmpsite[N];
    int occufield[2*10000*N];

    // initialize the occupancy field
    for (int i = 0; i < 2*10000*N; i++) {
        occufield[i] = -1;
    }

    // initialize site array and place the monomer in occufield
    for (int i = 0; i < N; i++) {
        if (i < N-1) simdir[i] = 100;
        site[i] = 10000*N + i*100;
        tmpsite[i] = site[i];
        occufield[site[i]] = amino_sequence[i]; 
    }

    // flag array for various purposes
    int flag[3] = {0, 1, 1}; // flag[0] for the flagoverlap ;  flag[1] for the flag2  ; flag[2] for the "+" or "-" direction in pivot move

    // variables to calculate the energy of the configuration
    int eHH = 4, nHH = 0, xHH = 0, yHH = 0;
    int eH0 = 2, nH0 = 0, xH0 = 0, yH0 = 0;
    int etheta = 1, ntheta = 0, xtheta = 0, ytheta = 0;
    int totalEnergy = 0;
    int currenergy = 0; 
    int prevenergy = 0;


    // histogram array to record the energy
    int energy_interval[] =  {0, 500, 1}; //the range of energy from 0 to 20 with interval 1
    int range = (energy_interval[1] - energy_interval[0]) * energy_interval[2];
    long long int his[range];           // array for histogram
    double dos[range];                    // array for the estimated dos for each energy level
    double normdos[range];                // normalized dos

    // initialize the value of histogram and dos
    for (int t = 0; t < range; t++) {
        his[t] = 0;
        dos[t] = 0.0;
        normdos[t] = 0.0;
    }
    
    // variables to count the numbers overlap and valid configurations
    long long int countoverlap = 0;
    long long int countvalid = 0;
    long long int moveaccepted = 0;
    long long int moverejected = 0;

    // parameter for wang-landau algorithm
    double modf = 1.0;          //modification factor
    double modfmin = 1.0e-8;     //modification factor minimum
    int checkinterval = 1000000;      //the iterval of mc step done before doing flatness check
    double flatcriteria = 0.8;          //the histogram flatness criteria
    int flatflag = 0;
    double threshold;     // minimum threshold for each histogram to have to satisfy flatness criteria
    int histentries =0;         // number of entries of all histogram bins
    int binactive = 0;          // number of active bins

    
    int picksite;       
    int pivottype;

    double rand2;       // random number variable
    int rngtype = 1;
    unsigned long int rngseed = 30;
    int wlstep = 0;             // i-th wang-landau step
    int flatcount = 0;          // number of sufficient flatness achieved




    gsl_rng* rng;

    switch (rngtype) {
    case  1 : {rng = gsl_rng_alloc(gsl_rng_ranlxd1); break;}
    case  2 : {rng = gsl_rng_alloc(gsl_rng_mt19937); break;}
    default : {rng = gsl_rng_alloc(gsl_rng_ranlxd2);}
    }

    if (rngseed == 0) rngseed = time(NULL) % ULONG_MAX;
    gsl_rng_set(rng, rngseed);

    // WL algorithm start
    while(modf >= modfmin ) {
        wlstep++;
        for (int hnum = 0; hnum < checkinterval; hnum++) {
            histentries++;

            //do pivot move
            // randomly pick the component direction of fold as starting point then do pivot move
            picksite = ((N-1)*gsl_rng_uniform(rng));
            pivottype = ((47)*gsl_rng_uniform(rng)) + 1;

            if (picksite >= N-picksite-2) 
                flag[2] = 1;
            else 
                flag[2] = -1;

            // printf("picksite = %d ", picksite);
            // printf("pivottype = %d \n", pivottype);

            // printf("flag[2] = %d \n", flag[2]);
            if (flag[2] == -1) {
                // save the state of the arrow
                for (int i = picksite; i >= 0; i--) {
                    tmpsimdir[i] = simdir[i];
                }
                
                // save the state of the site
                for (int i = picksite; i >= 0; i--) {
                    tmpsite[i] = site[i];
                }

                // calculate ntheta before
                for (int i = picksite; i >= 0; i--) {
                    if (simdir[i] != simdir[i+1]) xtheta++ ;
                }

                // printf("simdir = ");
                // for (int i = 0; i < N-1; i++ ) {
                //     printf(" %d ", simdir[i]);
                // }
                // printf("\n");
                
                // do pivot move
                dopivotmove(simdir, tmpsimdir, picksite, pivottype, flag);

                // printf("simdir = ");
                // for (int i = 0; i < N-1; i++ ) {
                //     printf(" %d ", simdir[i]);
                // }
                // printf("\n");

                // overlap check part-2
                if (flag[1]) {

                    if (picksite < N/4) {
                        for (int i = picksite; i >= 0; i--) {
                            
                            if(occufield[site[i]] == 0) {
                                // komponen H-H
                                if (occufield[site[i] + 1] == 0) xHH++;
                                if (occufield[site[i] - 1] == 0) xHH++;
                                if (occufield[site[i] + 100] == 0) xHH++;
                                if (occufield[site[i] - 100] == 0) xHH++;
                                if (occufield[site[i] + 10000] == 0) xHH++;
                                if (occufield[site[i] - 10000] == 0) xHH++;
                                // komponen H-0
                                if (occufield[site[i] + 1] == 2) xH0++;
                                if (occufield[site[i] - 1] == 2) xH0++;
                                if (occufield[site[i] + 100] == 2) xH0++;
                                if (occufield[site[i] - 100] == 2) xH0++;
                                if (occufield[site[i] + 10000] == 2) xH0++;
                                if (occufield[site[i] - 10000] == 2) xH0++;
                            }
                            if(occufield[site[i]] == 2) {
                                // komponen 0-H
                                if (occufield[site[i] + 1] == 0) xH0++;
                                if (occufield[site[i] - 1] == 0) xH0++;
                                if (occufield[site[i] + 100] == 0) xH0++;
                                if (occufield[site[i] - 100] == 0) xH0++;
                                if (occufield[site[i] + 10000] == 0) xH0++;
                                if (occufield[site[i] - 10000] == 0) xH0++;
                            }    
                        }
                    }

                    // update the site and occufield
                    for (int i = picksite; i >= 0; i--) {
                        // calculate the site elements (note that the site change is +1 of array index of simdir)
                        site[i] = site[i+1] - simdir[i];
                        
                        // update the occupancy field according to new site. do not forget to erase the previous occufield
                        occufield[tmpsite[i]] = -1;
                    }

                    // update the occupancy field according to new site if the site is not occupied
                    for (int i = picksite; i >= 0; i--) {
                        // if the site occupied then it is overlapped
                        if (occufield[site[i]] != -1) {
                            flag[0] = 1;

                            // erase the current occufield that result from the current site
                            for (int in = picksite; in > i; in--) {
                                occufield[site[in]] = -1;
                            }

                            // revert back occufield and site to previous state
                            for (int in = picksite; in >= 0; in--) {
                                occufield[tmpsite[in]] = amino_sequence[in];
                                site[in] = tmpsite[in];
                            }

                            // revert to previous arrow
                            for (int in = picksite; in >= 0; in--) {
                                simdir[in] = tmpsimdir[in];
                            }
                            break;
                        }
                        occufield[site[i]] = amino_sequence[i];
                    }
                }

                // printf("simdir = ");
                // for (int i = 0; i < N-1; i++ ) {
                //     printf(" %d ", simdir[i]);
                // }
                // printf("\n");

                // printf("site = ");
                // for (int i = 0; i < N; i++ ) {
                //     printf(" %d ", site[i]);
                // }
                // printf("\n");

                // printf("occufield = ");
                // for (int i = 0; i < N; i++ ) {
                //     printf(" %d ", occufield[site[i]]);
                // }
                // printf("\n");

                // printf("flagoverlap = %d  ", flag[0]);

                if (flag[0]) {
                    countoverlap++;
                    // printf("overlap configuration \n\n");
                }

                if (flag[0] == 0) { 
                    // calculate energy        
                    prevenergy = currenergy;

                    if (picksite < N/4) {
                        for (int i = picksite; i >= 0; i--) {
                            if(occufield[site[i]] == 0) {
                                // komponen H-H
                                if (occufield[site[i] + 1] == 0) yHH++;
                                if (occufield[site[i] - 1] == 0) yHH++;
                                if (occufield[site[i] + 100] == 0) yHH++;
                                if (occufield[site[i] - 100] == 0) yHH++;
                                if (occufield[site[i] + 10000] == 0) yHH++;
                                if (occufield[site[i] - 10000] == 0) yHH++;
                                // komponen H-0
                                if (occufield[site[i] + 1] == 2) yH0++;
                                if (occufield[site[i] - 1] == 2) yH0++;
                                if (occufield[site[i] + 100] == 2) yH0++;
                                if (occufield[site[i] - 100] == 2) yH0++;
                                if (occufield[site[i] + 10000] == 2) yH0++;
                                if (occufield[site[i] - 10000] == 2) yH0++;
                            }
                            
                            if(occufield[site[i]] == 2) {
                            // komponen 0-H
                                if (occufield[site[i] + 1] == 0) yH0++;
                                if (occufield[site[i] - 1] == 0) yH0++;
                                if (occufield[site[i] + 100] == 0) yH0++;
                                if (occufield[site[i] - 100] == 0) yH0++;
                                if (occufield[site[i] + 10000] == 0) yH0++;
                                if (occufield[site[i] - 10000] == 0) yH0++;
                            }        
                        }

                        // calculate ntheta after
                        for (int i = picksite; i >= 0; i--) {
                            if (simdir[i] != simdir[i+1]) ytheta++ ;
                        }

                        nHH = (yHH - xHH);
                        nH0 = (yH0 - xH0);
                        ntheta = (ytheta - xtheta);
                        totalEnergy = prevenergy + eHH*nHH + eH0*nH0 + etheta*ntheta;
                        // printf("type2 \n");
                    }

                    else {
                        for (int i = N-1; i >= 0; i--) {
                            
                            if (occufield[site[i]] == 0) {
                                // komponen H-H
                                if (occufield[site[i] + 1] == 0) nHH++;
                                if (occufield[site[i] + 100] == 0) nHH++;
                                if (occufield[site[i] + 10000] == 0) nHH++;
                                // komponen H-0
                                if (occufield[site[i] + 1] == 2) nH0++;
                                if (occufield[site[i] + 100] == 2) nH0++;
                                if (occufield[site[i] + 10000] == 2) nH0++;
                            }

                            // komponen 0-H
                            if (occufield[site[i]] == 2) {
                                if (occufield[site[i] + 1] == 0) nH0++;
                                if (occufield[site[i] + 100] == 0) nH0++;
                                if (occufield[site[i] + 10000] == 0) nH0++;
                            }
                        }

                        // calculate ntheta
                        for (int i = N-2; i >= 1; i--) {
                            if (simdir[i] != simdir[i-1]) ntheta++ ;
                        }

                        nHH = nHH-HHbond;
                        nH0 = nH0-H0bond;
                        totalEnergy = eHH*nHH + eH0*nH0 + etheta*ntheta;
                    }

                    
                    // printf("nHH = %d \n", nHH);
                    // printf("nH0 = %d \n", nH0);
                    // printf("ntheta = %d \n", ntheta);

                    countvalid++;
                    currenergy = totalEnergy;
                    // printf("Total energy = %d \n", totalEnergy);
                    // printf("\n");


                    // *WL acceptance criterion*
                    rand2 = gsl_rng_uniform(rng);
                    
                    if (rand2 < exp(dos[prevenergy]-dos[currenergy])) {
                    moveaccepted++;
                    }
                    else {
                        currenergy = prevenergy;
                        moverejected++;
                        // printf("check \n");
                        // revert to previous arrow
                        for (int in = picksite; in >= 0; in--) {
                            simdir[in] = tmpsimdir[in];
                        }

                        // when overlapped return the occufield and simdir to previous arrays
                        for (int in = picksite; in >= 0; in--) {
                            occufield[site[in]] = -1;
                        }

                        for (int in = picksite; in >= 0; in--) {
                            occufield[tmpsite[in]] = amino_sequence[in];
                        }

                        for (int in = picksite; in >= 0; in--) {
                            site[in] = tmpsite[in];
                        }

                        // printf("site = ");
                        // for (int i = 0; i < N; i++ ) {
                        //     printf(" %d ", site[i]);
                        // }
                        // printf("\n");
                    }
                            
                }
            }

            else if (flag[2] == 1) {
                // save the state of the arrow
                for (int i = picksite; i < N-1; i++) {
                    tmpsimdir[i] = simdir[i];
                }
                
                // save the state of the site
                for (int i = picksite + 1; i < N; i++) {
                    tmpsite[i] = site[i];
                }

                // calculate ntheta before
                for (int i = picksite; i < N-1; i++) {
                    if (simdir[i] != simdir[i-1]) xtheta++ ;
                }

                // printf("simdir = ");
                // for (int i = 0; i < N-1; i++ ) {
                //     printf(" %d ", simdir[i]);
                // }
                // printf("\n");
                
                // do pivot move
                dopivotmove(simdir, tmpsimdir, picksite, pivottype, flag);

                // printf("simdir = ");
                // for (int i = 0; i < N-1; i++ ) {
                //     printf(" %d ", simdir[i]);
                // }
                // printf("\n");

                // overlap check part-2
                if (flag[1]) {

                    if (picksite > 3*N/4) {
                        for (int i = picksite + 1; i < N; i++) {
                            if(occufield[site[i]] == 0) {
                                // komponen H-H
                                if (occufield[site[i] + 1] == 0) xHH++;
                                if (occufield[site[i] - 1] == 0) xHH++;
                                if (occufield[site[i] + 100] == 0) xHH++;
                                if (occufield[site[i] - 100] == 0) xHH++;
                                if (occufield[site[i] + 10000] == 0) xHH++;
                                if (occufield[site[i] - 10000] == 0) xHH++;
                                // komponen H-0
                                if (occufield[site[i] + 1] == 2) xH0++;
                                if (occufield[site[i] - 1] == 2) xH0++;
                                if (occufield[site[i] + 100] == 2) xH0++;
                                if (occufield[site[i] - 100] == 2) xH0++;
                                if (occufield[site[i] + 10000] == 2) xH0++;
                                if (occufield[site[i] - 10000] == 2) xH0++;
                            }
                            if(occufield[site[i]] == 2) {
                                // komponen 0-H
                                if (occufield[site[i] + 1] == 0) xH0++;
                                if (occufield[site[i] - 1] == 0) xH0++;
                                if (occufield[site[i] + 100] == 0) xH0++;
                                if (occufield[site[i] - 100] == 0) xH0++;
                                if (occufield[site[i] + 10000] == 0) xH0++;
                                if (occufield[site[i] - 10000] == 0) xH0++;
                            }    
                        }
                    }

                    // update the site and occufield
                    for (int i = picksite + 1; i < N; i++) {
                        // calculate the site elements (note that the site change is +1 of array index of simdir)
                        site[i] = site[i-1] + simdir[i-1];
                        
                        // update the occupancy field according to new site. do not forget to erase the previous occufield
                        occufield[tmpsite[i]] = -1;
                    }

                    // update the occupancy field according to new site if the site is not occupied
                    for (int i = picksite + 1; i < N; i++) {
                        // if the site occupied then it is overlapped
                        if (occufield[site[i]] != -1) {
                            flag[0] = 1;

                            // erase the current occufield that result from the current site
                            for (int in = picksite + 1; in < i; in++) {
                                occufield[site[in]] = -1;
                            }

                            // revert back occufield and site to previous state
                            for (int in = picksite + 1; in < N; in++) {
                                occufield[tmpsite[in]] = amino_sequence[in];
                                site[in] = tmpsite[in];
                            }

                            // revert to previous arrow
                            for (int in = picksite; in < N - 1; in++) {
                                simdir[in] = tmpsimdir[in];
                            }
                            break;
                        }
                        occufield[site[i]] = amino_sequence[i];
                    }
                }

                // printf("simdir = ");
                // for (int i = 0; i < N-1; i++ ) {
                //     printf(" %d ", simdir[i]);
                // }
                // printf("\n");

                // printf("site = ");
                // for (int i = 0; i < N; i++ ) {
                //     printf(" %d ", site[i]);
                // }
                // printf("\n");

                // printf("occufield = ");
                // for (int i = 0; i < N; i++ ) {
                //     printf(" %d ", occufield[site[i]]);
                // }
                // printf("\n");

                // printf("flagoverlap = %d  ", flag[0]);

                if (flag[0]) {
                    countoverlap++;
                    // printf("overlap configuration \n\n");
                }

                if (flag[0] == 0) { 
                    // calculate energy        
                    prevenergy = currenergy;

                    if (picksite > 3*N/4) {
                        for (int i = picksite + 1; i < N; i++) {
                            if(occufield[site[i]] == 0) {
                                // komponen H-H
                                if (occufield[site[i] + 1] == 0) yHH++;
                                if (occufield[site[i] - 1] == 0) yHH++;
                                if (occufield[site[i] + 100] == 0) yHH++;
                                if (occufield[site[i] - 100] == 0) yHH++;
                                if (occufield[site[i] + 10000] == 0) yHH++;
                                if (occufield[site[i] - 10000] == 0) yHH++;
                                // komponen H-0
                                if (occufield[site[i] + 1] == 2) yH0++;
                                if (occufield[site[i] - 1] == 2) yH0++;
                                if (occufield[site[i] + 100] == 2) yH0++;
                                if (occufield[site[i] - 100] == 2) yH0++;
                                if (occufield[site[i] + 10000] == 2) yH0++;
                                if (occufield[site[i] - 10000] == 2) yH0++;
                            }
                            
                            if(occufield[site[i]] == 2) {
                            // komponen 0-H
                                if (occufield[site[i] + 1] == 0) yH0++;
                                if (occufield[site[i] - 1] == 0) yH0++;
                                if (occufield[site[i] + 100] == 0) yH0++;
                                if (occufield[site[i] - 100] == 0) yH0++;
                                if (occufield[site[i] + 10000] == 0) yH0++;
                                if (occufield[site[i] - 10000] == 0) yH0++;
                            }    
                        }

                        // calculate ntheta after
                        for (int i = picksite; i < N-1; i++) {
                            if (simdir[i] != simdir[i-1]) ytheta++ ;
                        }

                        // printf("xHH[0] = %d  xHH[1] = %d \n", xHH, yHH);
                        nHH = (yHH - xHH);
                        nH0 = (yH0 - xH0);
                        ntheta = (ytheta - xtheta);
                        totalEnergy = prevenergy + eHH*nHH + eH0*nH0 + etheta*ntheta;
                        
                        // printf("type2 \n");
                    }

                    else {
                        for (int i = 0; i < N; i++) {
                            if (occufield[site[i]] == 0) {
                                // komponen H-H
                                if (occufield[site[i] + 1] == 0) nHH++;
                                if (occufield[site[i] + 100] == 0) nHH++;
                                if (occufield[site[i] + 10000] == 0) nHH++;
                                // komponen H-0
                                if (occufield[site[i] + 1] == 2) nH0++;
                                if (occufield[site[i] + 100] == 2) nH0++;
                                if (occufield[site[i] + 10000] == 2) nH0++;
                            }

                            // komponen 0-H
                            if (occufield[site[i]] == 2) {
                                if (occufield[site[i] + 1] == 0) nH0++;
                                if (occufield[site[i] + 100] == 0) nH0++;
                                if (occufield[site[i] + 10000] == 0) nH0++;
                            }
                        }

                        // calculate ntheta
                        for (int i = 0; i < N-2; i++) {
                            if (simdir[i] != simdir[i+1]) ntheta++ ;
                        }
                        nHH = nHH-HHbond;
                        nH0 = nH0-H0bond;

                        totalEnergy = eHH*nHH + eH0*nH0 + etheta*ntheta;
                        // printf("type1 \n");
                    }

                    countvalid++;
                    currenergy = totalEnergy;
                    // printf("Total energy = %d \n", totalEnergy);
                    // printf("\n");


                    // *WL acceptance criterion*
                    rand2 = gsl_rng_uniform(rng);
                    
                    if (rand2 < exp(dos[prevenergy]-dos[currenergy])) {
                    moveaccepted++;
                    }
                    else {
                        currenergy = prevenergy;
                        moverejected++;
                        // printf("check \n");
                        // revert to previous arrow
                        for (int in = picksite; in < N-1; in++) {
                            simdir[in] = tmpsimdir[in];
                        }

                        // when overlapped return the occufield and simdir to previous arrays
                        for (int in = picksite + 1; in < N; in++) {
                            occufield[site[in]] = -1;
                        }

                        for (int in = picksite + 1; in < N; in++) {
                            occufield[tmpsite[in]] = amino_sequence[in];
                        }

                        for (int in = picksite + 1; in < N; in++) {
                            site[in] = tmpsite[in];
                        }

                        // printf("site = ");
                        // for (int i = 0; i < N; i++ ) {
                        //     printf(" %d ", site[i]);
                        // }
                        // printf("\n");
                        // printf("sus \n");
                    }
                
                    
                }
            }
            
            // update the dos and histogram
            his[currenergy]++;
            dos[currenergy] += modf;  

            // printf("Accepted energy = %d \n", currenergy);
            

            // check histogram flatness
            if(hnum == checkinterval - 1) {
                binactive = get_binactive(dos, range);
                //check flatness
                threshold = flatcriteria * ((double)(histentries) / (double)(binactive));
                // printf("threshold = %f \n", threshold);
                flatflag = check_flatness(his, threshold, range);
                

                if (flatflag == 0) hnum = 0; 
                else {

                    modf = modf / 2;    
                    printf("\n");
                    printf("wlstep = %d \n", wlstep);
                    printf("modf = %e \n", modf);
                    // print the dos and histogram value
                    norm_dos(dos, range, normdos);
                    for (int i = 0; i < range; i++) {
                        if (dos[i] != 0.0) {
                            printf("%d      %.8e      \n", i, normdos[i]);
                        }
                    }

                    flatcount++;
                    binactive = 0;
                    histentries = 0;
                    flatflag = 0;
                    if (wlstep != 27 )
                    for (int i = 0; i < range; i++) {
                        his[i] = 0;
                    }
                }
            }        

            // set parameters to zero
            nHH = 0; xHH = 0; yHH = 0;
            nH0 = 0; xH0 = 0; yH0 = 0;
            ntheta = 0, xtheta = 0; ytheta=0;
            totalEnergy = 0;
            flag[0] = 0;
            flag[1] = 1;
            
        }
    }

    // calculate final normalized dos
    norm_dos(dos, range, normdos);

    gsl_rng_free(rng);

    printf("\n\n");
    printf("Final DOS =  \n");

    //print dos ke file output
    FILE *f;
    f = fopen("30SF1.txt", "w+");
    for (int i = 0; i < range; i++) {
        if (dos[i] != 0.0) {
            printf("%d      %.8e      \n", i, normdos[i]);
            fprintf(f, "%d      %.8e      \n", i, normdos[i]); 
        }
    }

    printf("\n\n");
    printf("flatcount = %d \n", flatcount);
    fprintf(f,"accepted moves = %lli  \n", moveaccepted);
    fprintf(f,"rejected moves = %lli  \n", moverejected);
    fprintf(f,"valid configuration = %lli  \n", countvalid);
    fprintf(f,"overlap configuration = %lli  \n", countoverlap);
    fprintf(f,"total configuration = %lli  \n", countvalid+countoverlap);
    
    printf("accepted moves = %lli  \n", moveaccepted);
    printf("rejected moves = %lli  \n", moverejected);
    printf("valid configuration = %lli  \n", countvalid);
    printf("overlap configuration = %lli  \n", countoverlap);
    printf("total configuration = %lli  \n", countvalid+countoverlap);
    
    fclose(f);

}

int get_binactive(double dos[], int range) {
    int binactive = 0;
    for (int i = 0; i < range; i++) {
        if (dos[i] != 0) binactive++;
    }
    return binactive;
}

int check_flatness(long long int his[], double threshold, int range) {
    int flag = 0; 
    for (int i = 0; i < range; i++) {
        if ((his[i] != 0) && (his[i] < threshold)) {
            flag = 0;
            break;
        }
        flag = 1; 
    }
    return flag;
}

void norm_dos(double dos[], int range, double normdos[]) {
    double dosmax = 0.0;
    double dum = 0;
    
    for (int i = 0; i < range; i++) {
        if (dos[i] != 0) {
            if (dos[i] > dosmax) dosmax = dos[i];
        }
    }

    for (int i = 0; i < range; i++) {
        dum += exp(dos[i] - dosmax);
    }

    for (int i = 0; i < range; i++) {
        normdos[i] =  exp(dos[i] - dosmax) / dum;
    }
}

void dopivotmove(int simdir[], int tmpsimdir[], int index, int pivottype, int flag[]) {
    int i = index;
    int step, limit;

    if (flag[2] == 1) {
        step = 1;
        limit = N-2;
    }
    else if (flag[2] == -1) {
        step = -1;
        limit = 0;
    }
    
    // printf("step = %d  limit = %d  \n", step, limit);

    while(1) {
    
        switch (pivottype) {
        
        case  0 : {
            simdir[i] *= 1;
            break;
        }

        case  1 : {
            simdir[i] *= -1; 
            break;
        }
        
        case  2 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -1;    
            break;
        }

        case  3 : {
            if(abs(simdir[i]) == 100) simdir[i] *= -1;
            break;
        }
        
        case  4 : {
            if(abs(simdir[i]) == 10000) simdir[i] *= -1;
            break;
        }

        case  5 : {
            if(abs(simdir[i]) != 1) simdir[i] *= -1;
            break;
        }

        case  6 : {
            if(abs(simdir[i]) != 100) simdir[i] *= -1;
            break;
        }

        case  7 : {
            if(abs(simdir[i]) != 10000) simdir[i] *= -1; 
            break;
        }

        case  8 : {
            if(abs(simdir[i]) == 100) simdir[i] *= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 100;
            break;
        }

        case  9 : {
            if(abs(simdir[i]) == 100) simdir[i] *= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -100;
            break;
        }

        case  10 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 10000;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -10000;
            break;
        }

        case  11 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -10000;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 10000;
            break;
        }

        case  12 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -100;
            else if(abs(simdir[i]) == 100) simdir[i] /= 100;
            break;
        }
        
        case  13 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 100;
            else if(abs(simdir[i]) == 100) simdir[i] /= -100;
            break;
        }

        case  14 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -1;
            else if(abs(simdir[i]) == 100) simdir[i] *= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -100;
            break;
        }

        case  15 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -1;
            else if(abs(simdir[i]) == 100) simdir[i] *= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 100;
            break;
        }

        case  16 : {
            if(abs(simdir[i]) == 100) simdir[i] *= -1;
            else if(abs(simdir[i]) == 1) simdir[i] *= -10000;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 10000;
            break;
        }

        case  17 : {
            if(abs(simdir[i]) == 100) simdir[i] *= -1;
            else if(abs(simdir[i]) == 1) simdir[i] *= 10000;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -10000;
            break;
        }

        case  18 : {
            if(abs(simdir[i]) == 10000) simdir[i] *= -1;
            else if(abs(simdir[i]) == 1) simdir[i] *= 100;
            else if(abs(simdir[i]) == 100) simdir[i] /= -100;
            break;
        }
        
        case  19 : {
            if(abs(simdir[i]) == 10000) simdir[i] *= -1;
            else if(abs(simdir[i]) == 1) simdir[i] *= -100;
            else if(abs(simdir[i]) == 100) simdir[i] /= 100;
            break;
        }
        
        case  20 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 100;
            break;
        }

        case  21 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 100;
            else if(abs(simdir[i]) == 100) simdir[i] *= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 10000;
            break;
        }

        case  22 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 100;
            break;
        }

        case  23 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -100;
            else if(abs(simdir[i]) == 100) simdir[i] *= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -10000;
            break;
        }

        case  24 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -100;
            break;
        }

        case  25 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 100;
            else if(abs(simdir[i]) == 100) simdir[i] *= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -10000;
            break;
        }

        case  26 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -100;
            break;
        }

        case  27 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -100;
            else if(abs(simdir[i]) == 100) simdir[i] *= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 10000;
            break;
        }

        case  28 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -100;
            break;
        }

        case  29 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -100;
            else if(abs(simdir[i]) == 100) simdir[i] *= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -10000;
            break;
        }

        case  30 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -100;
            break;
        }

        case  31 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 100;
            else if(abs(simdir[i]) == 100) simdir[i] *= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 10000;
            break;
        }

        case  32 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 100;
            break;
        }

        case  33 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -100;
            else if(abs(simdir[i]) == 100) simdir[i] *= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 10000;
            break;
        }

        case  34 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 100;
            break;
        }

        case  35 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 100;
            else if(abs(simdir[i]) == 100) simdir[i] *= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -10000;
            break;
        }

        case  36 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -1;
            else if(abs(simdir[i]) == 100) simdir[i] *= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 100;
            break;
        }

        case  37 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 10000;
            else if(abs(simdir[i]) == 100) simdir[i] *= -1;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 10000;
            break;
        }
        
        case  38 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 100;
            else if(abs(simdir[i]) == 100) simdir[i] /= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] *= -1;
            break;
        }

        case  39 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -1;
            else if(abs(simdir[i]) == 100) simdir[i] *= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -100;
            break;
        }

        case  40 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -10000;
            else if(abs(simdir[i]) == 100) simdir[i] *= -1;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -10000;
            break;
        }
        
        case  41 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -100;
            else if(abs(simdir[i]) == 100) simdir[i] /= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] *= -1;
            break;
        }

        case  42 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 1;
            else if(abs(simdir[i]) == 100) simdir[i] *= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -100;
            break;
        }

        case  43 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= 1;
            else if(abs(simdir[i]) == 10000) simdir[i] /= -10000;
            break;
        }
        
        case  44 : {
            if(abs(simdir[i]) == 1) simdir[i] *= -100;
            else if(abs(simdir[i]) == 100) simdir[i] /= -100;
            else if(abs(simdir[i]) == 10000) simdir[i] *= 1;
            break;
        }

        case  45 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 1;
            else if(abs(simdir[i]) == 100) simdir[i] *= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 100;
            break;
        }

        case  46 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 10000;
            else if(abs(simdir[i]) == 100) simdir[i] /= 1;
            else if(abs(simdir[i]) == 10000) simdir[i] /= 10000;
            break;
        }
        
        case  47 : {
            if(abs(simdir[i]) == 1) simdir[i] *= 100;
            else if(abs(simdir[i]) == 100) simdir[i] /= 100;
            else if(abs(simdir[i]) == 10000) simdir[i] *= 1;
            break;
        }
        
        }
        // printf("simdir = ");
        // for (int i = 0; i < N-1; i++) {  
        //     printf("  %d", simdir[i]);
        // }
        // printf("\n");
        // overlap check part 1
        if(simdir[i] == -simdir[i-step]) {
            flag[0] = 1;
            flag[1] = 0;
            
            // revert to previous arrow
            int in = index;
            while (1) {
                simdir[in] = tmpsimdir[in];
                if (in == i) break;
                in += step;
            }
            break;
        }

        

        if (i == limit) break;
        i += step;        
    }
    // printf("flagoverlap = %d  ", flag[0]);
    // printf("flag2 = %d \n", flag[1]);
}
