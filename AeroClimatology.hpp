/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AeroClimatology.hpp
 * Author: akheckel
 *
 * Created on 29. März 2017, 12:14
 */

#ifndef AEROCLIMATOLOGY_HPP
#define AEROCLIMATOLOGY_HPP

#include "InputParameter.hpp"


class AeroClimatology {
public:
    float *fineOfTotal, *weakOfFine, *dustOfCoarse, *aod;
    float lat_start, lon_start;
    float dlat, dlon;
    int nx, ny;

    AeroClimatology(const std::string climName, const int month);
    ~AeroClimatology();
private:
    AeroClimatology();
    AeroClimatology(const AeroClimatology& orig);           // disable copy contructor
    AeroClimatology& operator=(const AeroClimatology& rhs){ throw std::logic_error("assigning AeroClimatology not implemented"); } // disable copy assignment 

};

#endif /* AEROCLIMATOLOGY_HPP */

