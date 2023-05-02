#include<iostream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>

using namespace std;
#define PI 3.141592653589793238462643

const int q = 3;                                // q spin states
const int L = 16;                               // linear system size 
const double T = 0.5;

const int N = L;                                // Total number of spins 
const double pconnect = 1 - exp(-1.0/T);        // connection probability 

const int NCLUSTERS = 1;
const int NESTEPS = 100000;
const int NMSTEPS = 100000;
const int NBINS = 10;                           // # measurement bins (?)

vector<int> S(N);                               // spin array 
vector<int> M(q);                               // # spins in different states 
vector <complex<double> > W(q);                 // order parameter weights (?)

// lattice handling: 
enum dirs{RIGHT,LEFT};
int indx(int x){return x;}
// int indx(int x, int y){ return y * L + x; }                    // Make an index on every site 
int xpos(int i){ return i%L;}
// int ypos(int i){ return ;}

int Nbr(int i, int dir) {
    int x=xpos(i);
    switch (dir)
    {
        case RIGHT: 
            return indx((x+1)%L);
        case LEFT: 
            return indx((x-1+L)%L);
    }
    return 0;
}


void FlipandBuildFrom(int s) {
    int oldstate(S[s]), newstate((S[s]+1)%q);

    S[s] = newstate;
    M[oldstate]--; M[newstate]++;               // Update spin counts

    for(int dir = 0; dir < 2; dir++) {          // Go through neighbours
        int j = Nbr(s, dir);
        if(S[j] == oldstate)
            if( rand()/(RAND_MAX+1.) < pconnect)
                FlipandBuildFrom(j); 
    }
}


int main() 
{
    //Initialize order parameter weights
    for(int s = 0; s < q; s++) 
        W[s] = complex<double>(cos(2*PI*s/q), sin(2*PI*s/q)); //Gives the mj-value
    for(int i = 0; i < N; i++) S[i] = 0;
    for(int s = 1; s < q; s++) M[s] = 0;
    M[0] = N;
    srand((unsigned) time(0));

    // equilibriate 
    for(int t = 0; t < NESTEPS; t++) {
        for(int c = 0; c < NCLUSTERS; c++) {
            FlipandBuildFrom(rand()%N);
        }
    }

    // measure 
    for(int n = 0; n < NBINS; n++) {
        complex<double> m(0., 0.);
        double m1 = 0., m2 = 0., m4 = 0.;

        for(int t = 0; t < NMSTEPS; t++) {
            for(int c = 0; c < NCLUSTERS; c++)
                FlipandBuildFrom(rand()%N);

            complex<double> tm(0., 0.);
            for(int s = 0; s < q; s++)
                tm += W[s]*double(M[s]);

            tm /= N;

            double tm1 = abs(tm);
            double tm2 = tm1*tm1;

            m += tm; m1 += tm1; m2 += tm2; m4 += tm2*tm2;
        }

        m /= NMSTEPS; m1 /= NMSTEPS; m2 /= NMSTEPS; m4 /= NMSTEPS;
        cout << m << " " << m1 << " " << m2 << " " << m4 << endl;
    }

    return 0;
}