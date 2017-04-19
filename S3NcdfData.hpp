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
    const InputParameter &pars;
    const std::string& s3DataDir;
    /*const S3MetaData& s3MetaData;
    const float &szaLim;
    const bool &doGeoSubset;
    const float *latLim;
    const float *lonLim;*/

    S3BasicImage<signed char>   s3LandImg;  //L1b nadir resolution image of land mask
    S3BasicImage<signed char>   s3CloudImg; //L1b nadir resolution image of cloud mask
    S3BasicImage<signed char>   s3ValidImg; //L1b nadir resolution image of valid mask
    S3BasicImage<short>  s3RadImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS]; // TOA reflec images [views][bands]
    S3BasicImage<double> s3LatImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3LonImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3TpgLatImg;
    S3BasicImage<double> s3TpgLonImg;
    S3BasicImage<double> s3SzaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3SaaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3VzaImgs[N_SLSTR_VIEWS];
    S3BasicImage<double> s3VaaImgs[N_SLSTR_VIEWS];
    S3BasicImage<short>  flags; // retrieval resolution flag image
    S3BasicImage<short>  s3SdrImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS]; // SDR images [views][bands]
    S3BasicImage<float>  s3AodImgs[N_SLSTR_BANDS];                // AOD images [bands]
    S3BasicImage<float>  s3AerFracImgs[3];                        // fine_of_total, weak_of_fine, dust_of_coarse
    S3BasicImage<float>  s3FminImg;                               // fmin
    S3BasicImage<float>  s3UncImg;                                // uncer
    S3BasicImage<float>  s3RPathImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS]; // rPath images [views][bands]
    S3BasicImage<float>  s3TDownImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS]; // downw. transm images [views][bands]
    S3BasicImage<float>  s3TUpImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS];   // upw. transm images [views][bands]
    S3BasicImage<float>  s3TGasImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS];  // gas transm images [views][bands]
    S3BasicImage<float>  s3SpherAImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS]; // spher Albedo images [views][bands]
    S3BasicImage<float>  s3DifFracImgs[N_SLSTR_VIEWS][N_SLSTR_BANDS]; // diffuse Frac images [views][bands]

    double s3Irrad[N_SLSTR_VIEWS][N_SLSTR_BANDS];
    
    static const std::string SDR_NAMES[N_SLSTR_VIEWS][5];

    static const std::string RPATH_NAMES[N_SLSTR_VIEWS][5];
    static const std::string TDOWN_NAMES[N_SLSTR_VIEWS][5];
    static const std::string TUP_NAMES[N_SLSTR_VIEWS][5];
    static const std::string TGAS_NAMES[N_SLSTR_VIEWS][5];
    static const std::string SPHERA_NAMES[N_SLSTR_VIEWS][5];
    static const std::string DIFFRAC_NAMES[N_SLSTR_VIEWS][5];

    static const std::string AOD_NAMES[5];
    static const std::string AER_FRAC_NAMES[3];
    static const std::string FMIN_NAME;
    static const std::string UNC_NAME;
    
    S3NcdfData(const InputParameter& inPar);
    ~S3NcdfData();
    
    /*void readNcdf();*/
    void readNcdf(const ImageProperties& outImgProp);
    void convRad2Refl();
    void verifyInput();
    void initResultImgs(const ImageProperties& outImgProp);
    bool isValidPixel(const int& idx);
    void getGeoPos(const int& idx, GeoPos* gp);
    void getViewGeom(const int& idx, ViewGeom* vg);
    void getToaReflec(const int& idx, float tarr[][N_SLSTR_VIEWS]);
    void setRetrievalResults(const int& idx, SlstrPixel& pix);
    
private:
    int offCorr[2];
    bool isLandMaskAvailable;
    bool isFlagsImgAvailable;
    enum NcdfImageType {Nadir0500, Obliq0500, NadirTpg, ObliqTpg}; 
    
    S3NcdfData();
    S3NcdfData(const S3NcdfData& orig);           // disable copy construtor
    S3NcdfData& operator=(const S3NcdfData& rhs){ throw std::logic_error("S3NcdfData shoudlnt be copied"); } // disable copy assignment 

    void readImgBinned(S3BasicImage<short>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType);
    void readImgBinned(S3BasicImage<unsigned short>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType);
    void readSCloudS3SU(S3BasicImage<int>* s3Img, const std::string& ncdfName, const std::string varName);
    void readImgBinned(S3BasicImage<double>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType);
    void getFlagImg(S3BasicImage<unsigned short>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getBinRadImg(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getBinGeoLocImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void getBinGeomImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const netCDF::NcVar& imgVar);
    void split(const std::string& s, char c, std::vector<std::string>& v);
    void readImageProp(S3MetaData* s3md, ImageProperties* imgProp);
    void readImageProp(const netCDF::NcFile& ncF, ImageProperties* imgProp);
    void getVarAttSafely(S3BasicImage<short>* s3Img, netCDF::NcVar& imgVar);
    void getVarAttSafely(S3BasicImage<double>* s3Img, netCDF::NcVar& imgVar);
    bool hasAtt(const netCDF::NcVar& var, const std::string& attName);
    void readIrrad(double irrad[][N_SLSTR_BANDS]);
    void createLandMask();
    void createCloudMask();
    void createMyCloudMask();
    void createValidMask();
    void corrTpg(S3BasicImage<double>& tpg);
};

#endif /* S3NCDFDATA_H */

