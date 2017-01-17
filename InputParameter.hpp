/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputParameter.h
 * Author: akheckel
 *
 * Created on 23. November 2016, 10:39
 */

#ifndef INPUTPARAMETER_HPP
#define INPUTPARAMETER_HPP

#include <iostream>
#include "S3MetaData.hpp"

using std::string;

/**
 * input parameter class to hold parameters and 
 * parse cmd line and read/write par file
 */
class InputParameter {
public:
    InputParameter();
    ~InputParameter();               // std. destructor (empty dtor inlined below))
    
    string parFileName;
    string slstrProductDir;
    string aodOutDir;
    string atmLutFileName;
    string ocnLutFileName;
    string climFileName;
    
    float cloudFraction;
    int winSize, offset, skip;
    float szaLimit;
    float binValidThrs;
    
private:
    InputParameter(const InputParameter& orig){}                // disable copy constructor
    InputParameter& operator=(const InputParameter& rightVal){} // disable copy assignment
};

inline InputParameter::~InputParameter(){}

#endif /* INPUTPARAMETER_HPP */

