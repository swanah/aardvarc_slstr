/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   S3NcdfData.h
 * Author: akheckel
 *
 * Created on 25. November 2016, 10:02
 */

#ifndef S3NCDFDATA_HPP
#define S3NCDFDATA_HPP

#include <netcdf>
#include "defs.hpp"
#include "S3MetaData.hpp"
#include "Images.hpp"


class S3NcdfData {
public:
    std::string dataDir;
    S3BasicImage<short> s3RadImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS];
    S3BasicImage<double> s3LatImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3LonImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3TpgLatImg;
    S3BasicImage<double> s3TpgLonImg;
    S3BasicImage<double> s3SzaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3SaaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3VzaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3VaaImgs[N_SLSTR_VIEWS];
    
    S3NcdfData();
    ~S3NcdfData();
    
    void readNcdf(S3MetaData& s3md);
    
private:
    S3NcdfData(const S3NcdfData& orig){ throw std::logic_error("S3NcdfData shoudlnt be copied"); }           // disable copy construtor
    S3NcdfData& operator=(const S3NcdfData& rhs){ throw std::logic_error("S3NcdfData shoudlnt be copied"); } // disable copy assignment 

    void readImg(const std::string& ncdfName, const std::string varName, S3BasicImage<short>* s3Img);
    void readImg(const std::string& ncdfName, const std::string varName, S3BasicImage<double>* s3Img);
    void split(const std::string& s, char c, std::vector<std::string>& v);
    void readImageProp(S3MetaData* s3md, ImageProperties* imgProp);
    void readImageProp(const netCDF::NcFile& ncF, ImageProperties* imgProp);
    void getVarAttSafely(S3BasicImage<short>* s3Img, netCDF::NcVar& imgVar);
    void getVarAttSafely(S3BasicImage<double>* s3Img, netCDF::NcVar& imgVar);
    bool hasAtt(const netCDF::NcVar& var, const std::string& attName);};

#endif /* S3NCDFDATA_H */

