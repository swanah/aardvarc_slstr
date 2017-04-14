/*
 * File:   S3MetaData.h
 * Author: akheckel
 *
 * Created on 28. Oktober 2016, 21:36
 */

#ifndef S3METADATA_HPP
#define S3METADATA_HPP

#include <iostream>
#include "tinyxml/tinyxml.h"
#include "defs.hpp"


    class S3MetaData {                          // begin declaration of class
    public:                                     // begin public section
        S3MetaData(){}
        ~S3MetaData();                          // destructor (empty dtor inlined below)
        void parseManifest(const std::string& s3ProdDir);
        int getMonth();

        std::string s3FileName;
        std::string productName, productType;

        std::string startTime, stopTime;
        std::string platformFamily, platformNumber, instrument;
        std::string timeliness, baselineColl, creaTime;
        std::string prodUnitType, prodUnitDuration, prodUnitAlongTrackCoord;
        struct{
            ImageProperties nadirImg0500m, obliqImg0500m;
            ImageProperties nadirTpgImg, obliqTpgImg;
            int ncdfXoff, ncdfYoff;
        } slstrPInfo;
        
        static const std::string XFDU_MANIFEST_NAME;
        static const std::string CHANNEL_RAD_NAME[][5];
        static const std::string CHANNEL_QUAL_NAME[][5];
        static const std::string CHANNEL_IRRAD_NAME[][5];
        static const std::string GEODETIC_NAME[];
        static const std::string LAT_NAME[];
        static const std::string LON_NAME[];
        static const std::string GEODETIC_TPG_NAME;
        static const std::string LAT_TPG_NAME;
        static const std::string LON_TPG_NAME;
        static const std::string GEOMETRY_NAME[];
        static const std::string SZA_NAME[];
        static const std::string SAA_NAME[];
        static const std::string VZA_NAME[];
        static const std::string VAA_NAME[];
        static const std::string FLAGS_NAME[];
        static const std::string CONFID_NAME[];
        static const std::string BASIC_CLOUD_NAME[];
        static const std::string TIME_NAME[];
        
    private:
        S3MetaData(const S3MetaData& orig){}     // copy constructor
        S3MetaData& operator=(const S3MetaData& rightVal){}  //copy assignement

        void readAcquiTime(TiXmlElement* acquiTime);
        void readPlatformInfo(TiXmlElement* pf);
        void readGenProdInfo(TiXmlElement* genPInfo);
        void readSlstrProdInfo(TiXmlElement* slstrInfo);
        void readImageInfo(TiXmlNode* node, ImageProperties* imgInfo);
        void assertValidImgProp(const ImageProperties& imgProp);

    };

    inline S3MetaData::~S3MetaData() {
        //std::cout << "meta data class for \"" << s3FileName << "\" destroyed." << std::endl;
    }
    
    
#endif /* S3METADATA_H */

