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
#include "InputParameter.hpp"


class S3NcdfData {
public:
    std::string s3DataDir;
    std::string aodDataDir, aodOutName;
    S3MetaData s3MetaData;

    S3BasicImage<signed char>   s3LandImg;
    S3BasicImage<signed char>   s3CloudImg;
    S3BasicImage<signed char>   s3ValidImg;
    S3BasicImage<short>  s3RadImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS];
    S3BasicImage<double> s3LatImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3LonImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3TpgLatImg;
    S3BasicImage<double> s3TpgLonImg;
    S3BasicImage<double> s3SzaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3SaaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3VzaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3VaaImgs[N_SLSTR_VIEWS];
    S3BasicImage<short>  flags;

    double s3Irrad[N_SLSTR_VIEWS][N_SLSTR_BANDS];
    S3NcdfData(const InputParameter& inPar);
    ~S3NcdfData();
    
    /*void readNcdf();*/
    void readNcdf(const ImageProperties& outImgProp);
    void convRad2Refl();
    void verifyInput();
    
private:
    
    enum NcdfImageType {Nadir0500, Obliq0500, NadirTpg, ObliqTpg}; 
    
    S3NcdfData();
    S3NcdfData(const S3NcdfData& orig);           // disable copy construtor
    S3NcdfData& operator=(const S3NcdfData& rhs){ throw std::logic_error("S3NcdfData shoudlnt be copied"); } // disable copy assignment 

    void setAodDataDir(const InputParameter& inPar, const S3MetaData& s3MetaData);
    /*void readImg(S3BasicImage<short>* s3Img, const std::string& ncdfName, 
                 const std::string& varName, const NcdfImageType& imgType);
    void readImg(S3BasicImage<double>* s3Img, const std::string& ncdfName, 
                 const std::string& varName, const NcdfImageType& imgType);*/
    void readImgBinned(S3BasicImage<short>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType);
    void readImgBinned(S3BasicImage<ushort>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType);
    void readImgBinned(S3BasicImage<double>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType);
    void getFlagImg(S3BasicImage<ushort>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getBinRadImg(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getBinGeoLocImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    /*void getBinImgShifted(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getBinImgShifted(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getBinImgScaled(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);*/
    void getBinGeomImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);

    /*void getImgScaled(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getImgScaled(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getImgShifted(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getImgShifted(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);*/

    void split(const std::string& s, char c, std::vector<std::string>& v);
    void readImageProp(S3MetaData* s3md, ImageProperties* imgProp);
    void readImageProp(const netCDF::NcFile& ncF, ImageProperties* imgProp);
    void getVarAttSafely(S3BasicImage<short>* s3Img, netCDF::NcVar& imgVar);
    void getVarAttSafely(S3BasicImage<double>* s3Img, netCDF::NcVar& imgVar);
    bool hasAtt(const netCDF::NcVar& var, const std::string& attName);
    void readIrrad(double irrad[][N_SLSTR_BANDS]);
    void createLandMask();
    void createCloudMask();
    void createValidMask();
};

#endif /* S3NCDFDATA_H */

