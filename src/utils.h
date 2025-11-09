#ifndef UTILS
#define UTILS

#include <fstream>
#include <vector>
#include <string>


inline void WriteCSV(const std::string& filename,
                     const std::vector<double>& x,
                     const std::vector<double>& p,
                     const std::vector<double>& v,
                     const std::vector<double>& r)
{
    std::ofstream fout(filename);
    fout << "x,p,v,r\n";
    for (size_t i = 0; i < x.size(); ++i)
        fout << x[i] << "," << p[i] << "," << v[i] << "," << r[i] << "\n";
    fout.close();
}

inline double max( double a, double b ) { return ( a < b ) ? b : a; }

inline double min( double a, double b ) { return ( a < b ) ? a : b; }

#endif