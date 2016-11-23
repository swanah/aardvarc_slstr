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

#ifndef INPUTPARAMETER_H
#define INPUTPARAMETER_H

#include <iostream>
#include "S3MetaData.h"

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
    S3MetaData s3MetaData;
    string aodOutDir, aodOutName;
    string atmLutFileName;
    string ocnLutFileName;
    string climFileName;
    
    float cloudFraction;
    int winsize, offset, skip;
    
private:
    InputParameter(const InputParameter& orig){}                // disable copy constructor
    InputParameter& operator=(const InputParameter& rightVal){} // disable copy assignment

    void modifyAodOutNames();
    void replaceStringInPlace(std::string& subject, 
                          const std::string& search, const std::string& replace);
};

inline InputParameter::~InputParameter(){}

#endif /* INPUTPARAMETER_H */

