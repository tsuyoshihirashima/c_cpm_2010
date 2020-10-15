// This file is for the simulation of a simple case of cellular Potts model. The code is written in C.
// Be sure to check if the coding is correct by yourself.
// commented by Tsuyoshi Hirashima (tsuyoshi.hirashima @ gmail.com)

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

// Parameters for simulation condition //
#define XMAX (101) 
#define YMAX (101)
#define MCS_MAX (100) // Maximum value of monte-carlo steps
#define SIGMA_MAX (10) // Maximum value of sigma (Cell index), Note that this includes sigma=0(medium).

// Parameters for Cellular Dynamics //
#define J_EE (1)
#define J_EM (2)
#define LAM_AREA (1.0)
#define AREA (100)
#define BETA (1.0)

#define PI (3.14159265358979)

// Function //
int GetRandom(int min, int max){
return min + (int)( rand()*(max-min+1.0)/(1.0+RAND_MAX) );
}

// Declaration of Global Variables //
int sigma[YMAX][XMAX], vol[SIGMA_MAX];	

main(){
	int i, j, t, mcs;
	int rnd_x, rnd_y, c, c_neighbor, nb_x, nb_y;
	double e_adhesion, new_e_adhesion, e_vol, e_all, prob;
	
	// source of randomness; other methods are better, e.g., Mersenne Twister
	srand((unsigned int)time(NULL));  
	
	
	// ----------- Setting up an initial condition -----------
	for(i=0; i<XMAX; i++){
		for(j=0; j<YMAX; j++){
			sigma[j][i]=0;
		}
	}
	for(t=1; t<=SIGMA_MAX; t++){
		i=GetRandom(XMAX/4, XMAX/4*3);
		j=GetRandom(YMAX/4, YMAX/4*3);
		sigma[j][i]=t;
	}
	
	for(i=0; i<SIGMA_MAX; i++)
	vol[i]=0;
	
	for(i=0; i<XMAX; i++){
		for(j=0; j<YMAX; j++){
			vol[sigma[j][i]]++;
		}
	}
	
	
	// -- Cellular Dynamics -- //
	for(mcs=1; mcs<=MCS_MAX; mcs++){
		for(t=1; t<=16*XMAX*YMAX; t++){
			// ------------ Random choosing one site (1st chosen site)----------------	
			rnd_x= GetRandom(1, XMAX-2);
			rnd_y= GetRandom(1, YMAX-2);
			c = sigma[rnd_y][rnd_x]; 
			
			// ------------ Random choosing one site adjacent to the first chosen site (2nd chosen site)----------------	
			nb_x=GetRandom(1,3)-2; 
			nb_y=GetRandom(1,3)-2; 
			c_neighbor = sigma[rnd_y+nb_y][rnd_x+nb_x];
			
			if(c!= c_neighbor){ // You don't have to consider a case, in which the labeled value of 1st chosen site and that of 2nd chosen site is same.
				
				// --------------  Calculation of Interfacial Energy --------------  //
				e_adhesion = 0; new_e_adhesion = 0;
				for(i=-1; i<=1; i++){
					for(j=-1; j<=1; j++){
						if( !(i==0 && j==0) ){
							if(c != sigma[rnd_y+j][rnd_x+i]){
								if(c>0 && sigma[rnd_y+j][rnd_x+i]>0)
								e_adhesion += J_EE;
								else if( (c>0 && sigma[rnd_y+j][rnd_x+i]==0) || (c==0 && sigma[rnd_y+j][rnd_x+i]>0) )
								e_adhesion += J_EM;
							}
							if(c_neighbor != sigma[rnd_y+j][rnd_x+i]){
								if(c_neighbor>0 && sigma[rnd_y+j][rnd_x+i]>0)
								new_e_adhesion += J_EE;
								else if( (c_neighbor>0 && sigma[rnd_y+j][rnd_x+i]==0) || (c_neighbor==0 && sigma[rnd_y+j][rnd_x+i]>0) )
								new_e_adhesion += J_EM;
							}
						}
					}
				}
				
				// --------------  Calculation of Area Constraints --------------  //
				e_vol = 0;
				if(c>0)
				e_vol += LAM_AREA*(2*AREA - 2*vol[c] +1);
				if(c_neighbor>0)
				e_vol += LAM_AREA*(2*vol[c_neighbor] - 2*AREA +1);
				
				// --------------  Summation of all energies difference -----------------  //
				e_all = new_e_adhesion - e_adhesion + e_vol;
			
				if(e_all >= 0)
						prob = exp(-e_all*BETA);
				else if(e_all < 0)
						prob = 1.0;
				
				// ---------------------------------------------------------   UPDATE  ---------------------------------------------------------  
				if(prob >= (double) rand()/RAND_MAX ){ 
					sigma[rnd_y][rnd_x] = c_neighbor;
					
					// Update of volume //
					if(c>0)
					vol[c]--;
					if(c_neighbor>0)
					vol[c_neighbor]++;
				}
			}
		}
	}
}
