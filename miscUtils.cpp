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
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
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

/**
 * angular interpolation
 * @param a
 * @param b
 * @param w
 * @return 
 */
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

std::string getCurrentTimeStr() {
    boost::posix_time::ptime pt1 = boost::posix_time::second_clock::universal_time();
    boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%Y%m%dT%H%M%SZ");
    //facet->format("%Y%m%dT%H%M%SZ");
    
    std::stringstream stream;
    stream.imbue(std::locale(std::locale::classic(), facet));
    stream << pt1;
    
    return stream.str();
}

std::string convS3DateToCci(const std::string& s3DateTime) {
    boost::posix_time::time_input_facet* inFacet = new boost::posix_time::time_input_facet("%Y-%m-%dT%H:%M:%sZ");
    boost::posix_time::time_facet* outFacet = new boost::posix_time::time_facet("%Y%m%dT%H%M%SZ");
    
    const std::locale loc = std::locale(std::locale::classic(), inFacet);
    std::istringstream is(s3DateTime);
    is.imbue(loc);
    boost::posix_time::ptime pt1;
    is >> pt1;
    
    std::stringstream stream;
    stream.imbue(std::locale(std::locale::classic(), outFacet));
    stream << pt1;
    
    return stream.str();
}

double parseTimeToSec1970(const std::string& s3DateTime){
    const boost::posix_time::ptime time_epoch(boost::gregorian::date(1970, 1, 1));
    const boost::posix_time::time_input_facet* inFacet = new boost::posix_time::time_input_facet("%Y-%m-%dT%H:%M:%sZ");
    const std::locale loc = std::locale(std::locale::classic(), inFacet);
    std::istringstream is(s3DateTime);
    is.imbue(loc);
    boost::posix_time::ptime pt1;
    is >> pt1;
    

    boost::posix_time::time_duration dur = pt1 - time_epoch;
    long ms = dur.total_microseconds();
    //std::cout << "microseconds: " << ms << "\n";
    //boost::format output("%.6f");
    //output % (ms/1000000.0);
    //std::cout << output << std::endl;
    
    return (double)(1e-6 * ms);
}


/**
 * get current system time using high res clock and
 * convert to %Y-%m-%dT%H:%M:%S.microsec
 * @return 
 *
std::string getCurrentTimeStr(const char* format, const bool& appendMirco) {
    std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds ms = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch());
    size_t frac_sec = ms.count() % 1000000;
    
    std::time_t t = std::chrono::high_resolution_clock::to_time_t(now);
    tm* tm_time = std::gmtime(&t);
    char ftm[1000];
    std::strftime(ftm, 1000, format, tm_time);
    if (appendMirco) {
        snprintf(ftm, 1000, "%s.%dZ", ftm, frac_sec);
    }
    std::string tStr(ftm);
    return tStr;
}

/**
 * convert date string from S3 format yyyy-mm-ddThh:mm:ss.ssssssZ 
 * to cci format yyyymmddThhmmssZ
 * @param s3DateTime
 * @return 
 *
std::string convS3DateToCci(const std::string& s3DateTime) {
    std::regex regexp("^([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2}).([0-9]*)Z*$");
    std::smatch m;
    std::regex_match(s3DateTime, m, regexp);
    
    const size_t bufSiz = 100;
    char stmp[bufSiz];
    snprintf(stmp, bufSiz, "%s%s%sT%s%s%sZ", m[1].str().c_str(), m[2].str().c_str(), m[3].str().c_str(), 
                                             m[4].str().c_str(), m[5].str().c_str(), m[6].str().c_str());
    return std::string(stmp);
}

double parseTimeToSec1970(const std::string& s3DateTime){

    std::regex regexp("^([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2}).([0-9]*)Z*$");
    std::smatch m;
    std::regex_match(s3DateTime, m, regexp);
    
    std::tm tm = {0};
    tm.tm_year = std::strtol(m[1].str().c_str(), NULL, 10) - 1900;
    tm.tm_mon  = std::strtol(m[2].str().c_str(), NULL, 10) - 1;
    tm.tm_mday = std::strtol(m[3].str().c_str(), NULL, 10);
    tm.tm_hour = std::strtol(m[4].str().c_str(), NULL, 10);
    tm.tm_min  = std::strtol(m[5].str().c_str(), NULL, 10);
    tm.tm_sec  = std::strtol(m[6].str().c_str(), NULL, 10);
    tm.tm_isdst = -1;
    std::time_t tt = timegm(&tm);
    
    std::chrono::system_clock::time_point tp = std::chrono::system_clock::from_time_t(tt);
    std::chrono::microseconds ms1 = std::chrono::duration_cast<std::chrono::microseconds>(tp.time_since_epoch());
    
    int microSec = std::strtol(m[7].str().c_str(), NULL, 10);
    tp += std::chrono::microseconds(microSec);
    std::chrono::microseconds ms = std::chrono::duration_cast<std::chrono::microseconds>(tp.time_since_epoch());
    
    return 1e-6 * ms.count();
}
/**/

