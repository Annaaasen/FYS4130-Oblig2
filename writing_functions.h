#include<iostream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>

#include<iomanip>
#include<fstream>

using namespace std;
// ERROR CODE
#define ERR_FILE -2

string makePrefix(int q, int N, double T);

int writeCorrelation(string prefix, complex<double> *corr_func);


string makePrefix(int q, int L, double T) {
    string filename;
    ostringstream T_2d;
    T_2d << std::fixed << std::setprecision(2) << T;
    std::string T_str = T_2d.str();
    filename = to_string(q) + "," + to_string(L) + "," + T_str;
    return filename;
}

int writeCorrelation(int N, string prefix, complex<double> *corr_func) {
    string corr_filename = "correlation," + prefix + ".csv";
    ofstream corr_file(corr_filename);
    
    if (!corr_file.is_open())
        return ERR_FILE;
    for (int j=0; j<N; j++) {
        //        cout << j << "," << Cr[j].real() << endl;
        corr_file << j << "," << corr_func[j].real() << endl;
    }
    corr_file.close();
    cout << "Correlation written to file: " + corr_filename << endl;
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