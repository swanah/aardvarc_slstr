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

//using std::string;

/**
 * input parameter class to hold parameters and 
 * parse cmd line and read/write par file
 */
class InputParameter {
public:
    InputParameter(int argc, char** argv);
    ~InputParameter();               // std. destructor (empty dtor inlined below))
    
    std::string parFileName;
    std::string slstrProductDir;
    std::string aodOutDir, aodOutName;
    std::string cciOutDir, cciOutName;
    std::string atmLutFileName;
    std::string ocnLutFileName;
    std::string climFileName;
    
    bool applySimpleRecalib;
    float simpleRecalFactors[N_SLSTR_VIEWS][N_SLSTR_BANDS];
    
    bool applyCloudFilter;
    bool doCciOutput;
    bool useSCloudS3SU;
    std::string sCloudS3SUpath, sCloudS3SUname;
    int winSize; //size of bin window
    int offset;  //offset center bin = winSize/2
    int skip;    //unused
    float szaLimit;
    float binValidThrs;
    bool doGeoSubset;
    float latLim[2];
    float lonLim[2];
    S3MetaData s3MD;
    
private:
    InputParameter();
    InputParameter(const InputParameter& orig){}                // disable copy constructor
    InputParameter& operator=(const InputParameter& rightVal){} // disable copy assignment
    
    void parsePlaceholders(std::string& s);
    void ensurePathExists(std::string& path);
    void createAodName();
    void createCldName();
    void createCciAodName();
};

inline InputParameter::~InputParameter(){}

#endif /* INPUTPARAMETER_HPP */

