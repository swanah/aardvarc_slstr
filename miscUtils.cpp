/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <cstdlib>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include "miscUtils.hpp"

void replaceStringInPlace(std::string& subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    if (search.empty()) return;
    while((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}

void trimLine(std::string& s){
    if (s.empty()) return;
    int i=0;
    int j,n;
    j=n=s.size();
    while (i<n && s[i]<=0x20) i++;
    while (j>=0 && s[j]<=0x20) j--;
    if (j<n) s.erase(j+1, n-j-1);
    if (i>0) s.erase(0, i);
}

void splitLine(const std::string& s, std::vector<std::string>& words){
    int i=0, j;
    int n=s.size();
    words.clear();
    do {
        j=i;
        while (i<n && s[i]>0x20) i++;
        if (i>j){
            if (i<=n) words.push_back(s.substr(j, i-j));
            else return;
        }
        while (i<n && s[i]<=0x20) i++;
    } while (i<n);

}

double timer() {
    struct timeval tv;
    struct timezone tz;
    double walltime;

    gettimeofday(&tv, &tz);
    walltime = tv.tv_sec + tv.tv_usec * 1.0e-6;

    return (walltime);
}

void setBit(unsigned short *flag, unsigned short bitvalue, char test) {
    if (test) *flag |= bitvalue;
    else *flag &= ~bitvalue;
}

void setBit(short *flag, short bitvalue, char test) {
    if (test) *flag |= bitvalue;
    else *flag &= ~bitvalue;
}

double intAng(double a, double b, double w){
    double d2r = acos(-1.)/180.;
    a *= d2r;
    b *= d2r;
    double cs = (1-w)*cos(a) + w*cos(b);
    double sn = (1-w)*sin(a) + w*sin(b);
    return atan2(sn,cs)/d2r;
}

bool fileExists(const std::string& s){
    std::ifstream f(s.c_str());
    return f.good();
}

bool isDir(const std::string& s){
    struct stat sStat;
    stat(s.c_str(), &sStat);
    return S_ISDIR(sStat.st_mode);
}

bool createDirsUnix(const std::string& s){
    std::stringstream ss;
    ss << "mkdir -vp " << s;
    int err = system(ss.str().c_str());
    return err==0;
}
