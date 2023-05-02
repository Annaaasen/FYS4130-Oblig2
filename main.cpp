#include<iostream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>

#include<iomanip>
#include<fstream>
#include "writing_functions.h"

using namespace std;
#define PI 3.141592653589793238462643


int power(int base, int exp) {
    int result = 1;
    for (int i = 0; i < exp; i++) {
        result *= base;
    }
    return result;
}

const int q = 3;                                // q spin states
const int L = 16;                               // linear system size 

const double MAX_TEMP = 0.5;
const double TEMP_STEP = 0.25;

const int NDIMS = 1;                            // Number of dimensions 
const int N = power(L, NDIMS);                  // Total number of spins 

const int NCLUSTERS = 1;
const int NESTEPS = 100000;
const int NMSTEPS = 100000;
const int NBINS = 10;                           // # measurement bins (# experiments)

vector<int> S(N);                               // spin array 
vector<int> M(q);                               // # spins in different states 
vector <complex<double> > W(q);                 // returns mj-value given a certain state (s)

// lattice handling: 
enum dirs{RIGHT,LEFT, UP, DOWN};
// int indx(int x){return x;}
int indx(int x, int y){ return y * L + x; }     // Make an index on every site 
int xpos(int i){ return i%L;}
int ypos(int i){ return i/L;}

int Nbr(int i, int dir) {
    // i is the index of the point for which we want to find the neighbouring indeces
    // dir is in which direction we want to find the neighbouring index (is an integer which refers to the 'list' dirs)
    int x=xpos(i);
    int y=ypos(i);
    switch (dir)
    {
        case RIGHT: return indx((x+1)%L, y);
        case LEFT: return indx((x-1+L)%L, y);
        case UP: return indx(x, ((y+1)%L));
        case DOWN: return indx(x,(y-1+L)%L);
    }
    return 0;
}

void FlipandBuildFrom(int s, double pconnect) {
    // In this case s refers to the site (I think), a random site is often sent in 
    int oldstate(S[s]), newstate((S[s]+1)%q); //asigning the value of S[s] to oldstate, the newstate will be one state 'larger'

    S[s] = newstate;
    M[oldstate]--; M[newstate]++;                       // Update spin counts

    for(int dir = 0; dir < NDIMS*2; dir++) {            // Go through neighbours
        int j = Nbr(s, dir);                            // How can j become something?? Nbr returns 0, but only if dir is out of reach of dirs (eg dir=4)
        if(S[j] == oldstate)                            // If the neighbour has the same state as before the flip (shouldn't this also happen for the third possible state all that matters is that they are in different states, right?)
            if( rand()/(RAND_MAX+1.) < pconnect)
                FlipandBuildFrom(j, pconnect); 
    }
}


int main() 
{
    double T = 0.0;
    while(T<MAX_TEMP){
        T += TEMP_STEP;

        double pconnect = 1 - exp(-1.0/T);                        // connection probability 

        //Initialize order parameter weights
        for(int s = 0; s < q; s++) 
            W[s] = complex<double>(cos(2*PI*s/q), sin(2*PI*s/q)); //Gives the mj-value
        for(int i = 0; i < N; i++) S[i] = 0;
        for(int s = 1; s < q; s++) M[s] = 0;
        M[0] = N;                                                 //All particles start in state 0
        srand((unsigned) time(0));

        // equilibriate 
        for(int t = 0; t < NESTEPS; t++) {
            for(int c = 0; c < NCLUSTERS; c++) {
                FlipandBuildFrom(rand()%N, pconnect);                     //rand()%N Make sure a number larger than N-1 is never sent in 
            }
        }

        // measure 
        for(int n = 0; n < NBINS; n++) {
            complex<double> m(0., 0.);                          // Initialize (start with no magnetization)
            double m1 = 0., m2 = 0., m4 = 0.;

            complex<double> m0c = 0.;                           //Initialize the corr_func factors
            complex<double> mr[N],  m0cmr[N];                   // For all r's
            for(int j=0; j<N; j++){                         
                    m0cmr[j] = 0;
                    mr[j] = 0;
            }


            for(int t = 0; t < NMSTEPS; t++) {
                for(int c = 0; c < NCLUSTERS; c++)              //If NCLUSTERS>1, there are more than one cluster, and this will be start several flipandbuildfrom
                    FlipandBuildFrom(rand()%N, pconnect);                 //Determines all of the states in the grid after the Wolff cluser algo

                complex<double> tm(0., 0.);                     //'counts' the magnetizations given the specific grid 
                for(int s = 0; s < q; s++)
                    tm += W[s]*double(M[s]);

                tm /= N;

                double tm1 = abs(tm);
                double tm2 = tm1*tm1;

                m += tm; m1 += tm1; m2 += tm2; m4 += tm2*tm2;

                //factors to calculate correlation:
                m0c += conj(W[S[0]]);  
                for(int j=0; j<N; j++){                         // For all the r's
                    m0cmr[j] += conj(W[S[0]]) * W[S[j]];
                    mr[j] += W[S[j]];
                }

            }

            //Take the average: 
            m /= NMSTEPS; m1 /= NMSTEPS; m2 /= NMSTEPS; m4 /= NMSTEPS; 
            m0c/= NMSTEPS; 
            for(int j=0; j<N; j++){
                m0cmr[j] /= NMSTEPS; mr[j] /= NMSTEPS;
            }

            //calculating correlation function:
            complex<double> corr_func[N];
            for(int j=0; j<N; j++){
                corr_func[j] = m0cmr[j] - m0c * mr[j];
            }


            //Writing to file: 
            string prefix = makePrefix(q, L, T);

            if (NDIMS == 1 && (T == 0.25 || T == 0.5)) {
                if (writeCorrelation(N, prefix, corr_func))
                    return ERR_FILE;
            }
            

            //Printing:
            cout << m << " " << m1 << " " << m2 << " " << m4 << " "<< m0c << endl;
            // cout << m0cmr[5] << endl;
        }
    }

    return 0;
}



// int writeMoment(string m_filename, double T, complex<double> m, double m1, double m2, double m4) {
//     ofstream moment_file(m_filename, ios_base::app);

//     if (!moment_file.is_open())
//         return ERR_FILE;

//     moment_file << T << "," << m.real() << "," << m1 << "," << m2 << "," << m4 << endl ;
//     moment_file.close();
//     cout << "moment written to file: " + m_filename << endl;
//     return 0;
// }