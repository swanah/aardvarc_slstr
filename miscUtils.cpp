/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <sys/time.h>
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

double timer() {
    struct timeval tv;
    struct timezone tz;
    double walltime;

    gettimeofday(&tv, &tz);
    walltime = tv.tv_sec + tv.tv_usec * 1.0e-6;

    return (walltime);
}

