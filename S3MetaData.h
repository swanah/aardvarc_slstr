/*
 * File:   S3MetaData.h
 * Author: akheckel
 *
 * Created on 28. Oktober 2016, 21:36
 */

#ifndef S3METADATA_H
#define S3METADATA_H

#include <iostream>
#include "tinyxml/tinyxml.h"

    typedef struct {
        int width, height;
        int xOff, yOff;
    } imageInfo;
    
    
    class S3MetaData {                          // begin declaration of class
    public:                                     // begin public section
        S3MetaData(){}
        ~S3MetaData();                          // destructor (empty dtor inlined below)
        void parseManifest(const std::string& s3ProdDir);

        std::string s3FileName;
        std::string productName, productType;

        std::string startTime, stopTime;
        std::string platformFamily, platformNumber, instrument;
        std::string timeliness, baselineColl, creaTime;
        std::string prodUnitType, prodUnitDuration, prodUnitAlongTrackCoord;
        struct{
            int res0500m;
            int resTpg;
            imageInfo nadirImg0500m, obliqImg0500m;
            imageInfo nadirTpgImg, obliqTpgImg;
        } slstrPInfo;
        
        static const std::string XFDU_MANIFEST_NAME;
        static const std::string CHANNEL_RAD_NAME[][5];
        static const std::string CHANNEL_QUAL_NAME[][5];
        static const std::string GEODECTIC_NAME[];
        static const std::string GEOMETRY_NAME[];
        static const std::string FLAGS_NAME[];
        static const std::string TIME_NAME[];
        
    private:
        S3MetaData(const S3MetaData& orig){}     // copy constructor
        S3MetaData& operator=(const S3MetaData& rightVal){}  //copy assignement

        void readAcquiTime(TiXmlElement* acquiTime);
        void readPlatformInfo(TiXmlElement* pf);
        void readGenProdInfo(TiXmlElement* genPInfo);
        void readSlstrProdInfo(TiXmlElement* slstrInfo);
        void readImageInfo(TiXmlNode* node, imageInfo* imgInfo);

    };

    inline S3MetaData::~S3MetaData() {
        //std::cout << "meta data class for \"" << s3FileName << "\" destroyed." << std::endl;
    }


#endif /* S3METADATA_H */

