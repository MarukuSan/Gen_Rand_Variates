/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <stdio.h>
#include <stdlib.h> //--- Pour utiliser "malloc"
#include <math.h> //--- Pour utiliser "log"


/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
#define TAILLE 6

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

//--- Question 2 
double uniform(double a, double b) {
    return (genrand_real2()*(b-a))+a; //--- Génère un nombre pseudo-aléatoire entre a et b
}

//--- Question 3a
int* simul_classes(int repetition) {
    double random;
    int n = 3;

    int* tab = malloc(n*sizeof(int));
    for (int i=0; i<n; i++) {
        tab[i] = 0;
    }


    printf("Drawing : %d \n\n", repetition); //--- Affiche le nombre de drawing

    for (int i=0; i<repetition; i++) {
        random = genrand_real1(); //--- Génère un nombre dans [0-1]

        if (random<=0.5) {
            printf("Classe A : %f \n", random);
            tab[0]++;
        }
        else {
            if (random<=0.65) {
                printf("Classe B : %f \n", random);
                tab[1]++;
            }
            else {
                printf("Classe C : %f \n", random);
                tab[2]++;
            }
        }
    }

    return tab; //--- Repartition des classes
}

//--- Question 3b
void simul_classes_2(int taille, int tab[]) {
    float* proba = malloc(taille*sizeof(float));
    float* proba_cumul = malloc(taille*sizeof(float));

    int somme = 0; //--- Somme des drawing

    //--- Calcule de la somme des drawing
    for (int i=0; i<taille; i++) {
        somme = somme + tab[i];
    }

    //--- Probabilite d'appartenir à chaque classe
    char classe = 'A';
    printf("Tableau des probabilites de chaque classe : \n");
    for (int i=0; i<taille; i++) {
        proba[i] = (float)tab[i]/(float)somme;
        printf(" Classe %c : %f |", classe, proba[i]);
        classe++;
    }
    printf("\n");


    //--- Probabilite cumulee
    proba_cumul[0] = proba[0];
    for (int i=1; i<taille; i++) {
        proba_cumul[i] = proba_cumul[i-1] + proba[i];
    }

    //--- Affichage du tableau
    classe = 'A';
    printf("Tableau des probabilites cumulees de chaque classe : \n");
    for (int i=0; i<taille; i++) {
        printf(" Classe %c : %f |", classe, proba_cumul[i]);
        classe++;
    }
    printf("\n");
}

//--- Question 4a
double negExp(double mean) {
    double a = log(1.0 - genrand_real1());
    return -mean * a;
}



int main(int argc, char const *argv[])
{

    //--- Question 2 (Test)
    /*

    for (int i=0; i<=10; i++) {
        printf("%f \n", uniform(-89.2, 56.7));
    }
    
    */
    //--- Fin Question 2 

    //--- Question 3a (Test)
    /*

    int d1 = 1000;
    simul_classes(d1);

    //int d2 = 10000;
    //simul_classes(d2);
    
    */
    //--- Fin Question 3a

    //--- Question 3b (Test)

    //--- Exemple 1 :
    /*
    int d3 = 1000;
    simul_classes_2(3, simul_classes(d3));
    */

    //--- Exemple 2 (Dans le cours): 
    /*
    int tab[TAILLE] = {100, 400, 600, 400, 100, 200};
    simul_classes_2(TAILLE, tab);
    */
    //--- Fin Question 3b

    //--- Question 4b
    /*
    float mean = 11.0;
    printf("negExp(%f) = %f \n", mean, negExp(mean));

    double somme = 0;
    int repetition = 10000;
    
    for (int i=0; i<repetition; i++) {
        somme = somme + negExp(11);
    }

    printf("Average after %d drawing = %f \n", repetition, somme/(double)repetition);
    */
    //--- Fin Question 4b

    //--- Question 4c
    
    /*
    int tab2[22] = {0};
    int d3 = 1000;
    
    double a; //--- Stocker negExp()
    int b; //--- Stocker la partie entière de a
    
    
    for (int i=0; i<d3; i++) {
    	a = negExp(11);
    	b = (int)a;
    	
    	if (b<21) {
    		//--- Si le nombre est dans [0-1[, il se placera dans i=0.
    	        //--- Et ainsi de suite.
    	    
    	    tab2[b]++;
    	}
    	else {
    	    tab2[21]++;
    	}
    }
    
    int somme = 0;
    for (int i=0; i<22; i++) {
        printf(" %d |", tab2[i]);
        somme = somme + tab2[i];
    }
    printf("\n");
    */
    
    //--- Fin Question 4c
    
    //--- Question 5a


    return 0;
}
