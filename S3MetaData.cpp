/* 
 * File:   S3MetaData.cpp
 * Author: akheckel
 * 
 * Created on 28. Oktober 2016, 21:36
 */

#include <iostream>
#include <typeinfo>

#include "tinyxml/tinyxml.h"
#include "miscUtils.hpp"
#include "S3MetaData.hpp"

/* const definition */
    const std::string S3MetaData::XFDU_MANIFEST_NAME = "xfdumanifest.xml";
    const std::string S3MetaData::CHANNEL_RAD_NAME[][5] = {{"S1_radiance_an", "S2_radiance_an", "S3_radiance_an", "S5_radiance_an", "S6_radiance_an"},
                                                          {"S1_radiance_ao", "S2_radiance_ao", "S3_radiance_ao", "S5_radiance_ao", "S6_radiance_ao"}};
    const std::string S3MetaData::CHANNEL_QUAL_NAME[][5] = {{"S1_quality_an", "S2_quality_an", "S3_quality_an", "S5_quality_an", "S6_quality_an"},
                                                           {"S1_quality_ao", "S2_quality_ao", "S3_quality_ao", "S5_quality_ao", "S6_quality_ao"}};
    const std::string S3MetaData::CHANNEL_IRRAD_NAME[][5] = {{"S1_solar_irradiance_an", "S2_solar_irradiance_an", "S3_solar_irradiance_an", "S5_solar_irradiance_an", "S6_solar_irradiance_an"},
                                                           {"S1_solar_irradiance_ao", "S2_solar_irradiance_ao", "S3_solar_irradiance_ao", "S5_solar_irradiance_ao", "S6_solar_irradiance_ao"}};
    const std::string S3MetaData::GEODETIC_NAME[] = {"geodetic_an", "geodetic_ao"};
    const std::string S3MetaData::LAT_NAME[] = {"latitude_an", "latitude_ao"};
    const std::string S3MetaData::LON_NAME[] = {"longitude_an", "longitude_ao"};
    const std::string S3MetaData::GEODETIC_TPG_NAME = "geodetic_tx";
    const std::string S3MetaData::LAT_TPG_NAME = "latitude_tx";
    const std::string S3MetaData::LON_TPG_NAME = "longitude_tx";
    const std::string S3MetaData::GEOMETRY_NAME[] = {"geometry_tn", "geometry_to"};
    const std::string S3MetaData::SZA_NAME[] = {"solar_zenith_tn", "solar_zenith_to"};
    const std::string S3MetaData::SAA_NAME[] = {"solar_azimuth_tn", "solar_azimuth_to"};
    const std::string S3MetaData::VZA_NAME[] = {"sat_zenith_tn", "sat_zenith_to"};
    const std::string S3MetaData::VAA_NAME[] = {"sat_azimuth_tn", "sat_azimuth_to"};
    const std::string S3MetaData::FLAGS_NAME[] = {"flags_an", "flags_ao"};
    const std::string S3MetaData::CONFID_NAME[] = {"confidence_an", "confidence_ao"};
    const std::string S3MetaData::BASIC_CLOUD_NAME[] = {"cloud_an", "cloud_ao"};
    const std::string S3MetaData::S3SU_CLOUD_FNAME = "sCloudS3SU";
    const std::string S3MetaData::S3SU_CLOUD_VNAME[] = {"cld_n", "cld_o"};
    const std::string S3MetaData::TIME_NAME[] = {"time_an", "time_ao"};
    
    
/*    
 * no special con-/destructors 
 */  
    
/*
 * function definitions
 */
/**
 * constructor taking sentinel 3 product dir
 * to read xfdu manifest metadata
 * @param s3ProdDir
 */    
void S3MetaData::parseManifest(const std::string& s3ProdDir) {
    
    s3FileName = s3ProdDir;
    
    TiXmlDocument doc(s3ProdDir + "/" + XFDU_MANIFEST_NAME);
    bool ok = doc.LoadFile();
    if (ok) {
        
        /*FILE* f = fopen("E:/Temp/xmlTest.out", "wa");
        dump_to_file(f, &doc);
        fclose(f);*/
        
        TiXmlHandle hDoc(&doc);
        //TiXmlElement* metaDataElement = hDoc.FirstChildElement("xfdu:XFDU").FirstChildElement("metadataSection").Element();

        // get metaDataSection as Element
        TiXmlElement* mdSectionEle = hDoc.FirstChild("xfdu:XFDU").FirstChild("metadataSection").ToElement();
        
        // loop through metadataObjects
        TiXmlElement* elem;
        for (elem = mdSectionEle->FirstChildElement("metadataObject"); elem; elem = elem->NextSiblingElement("metadataObject")){
            if (elem) {
                const std::string* nodeId = elem->Attribute(std::string("ID"));
                //printf("Node: %s  Obj: %s\n", elem->Value(), nodeId->c_str());
                TiXmlElement* subEle;
                if (*nodeId == "acquisitionPeriod"){
                    subEle = elem->FirstChildElement("metadataWrap")->FirstChildElement("xmlData")
                                              ->FirstChildElement("sentinel-safe:acquisitionPeriod");
                    readAcquiTime(subEle);
                }
                else if (*nodeId == "platform"){
                    subEle = elem->FirstChildElement("metadataWrap")->FirstChildElement("xmlData")
                                   ->FirstChildElement("sentinel-safe:platform");
                    readPlatformInfo(subEle);
                }
                else if (*nodeId == "generalProductInformation"){
                    subEle = elem->FirstChildElement("metadataWrap")->FirstChildElement("xmlData")
                                   ->FirstChildElement("sentinel3:generalProductInformation");
                    readGenProdInfo(subEle);
                }
                else if (*nodeId == "slstrProductInformation"){
                    subEle = elem->FirstChildElement("metadataWrap")->FirstChildElement("xmlData")
                                   ->FirstChildElement("slstr:slstrProductInformation");
                    readSlstrProdInfo(subEle);
                }
            }
        }
        
    }
    else {
        throw std::runtime_error("shiiit!");
    }
    
    //TODO: assert i/o files exist and input is correct!!!
    assertValidImgProp(slstrPInfo.nadirImg0500m);
    assertValidImgProp(slstrPInfo.obliqImg0500m);
    assertValidImgProp(slstrPInfo.nadirTpgImg);
    assertValidImgProp(slstrPInfo.obliqTpgImg);
    
}

/**
 * provide month of startTime as Integer
 * @return month of year (startTime)
 */
int S3MetaData::getMonth(){
    // expecting startTime to be of format yyyy-MM-ddThh:mm:ss.sssssZ
    int month = StringToNumber<int>(startTime.substr(5, 2));
    return month;
}

/**
 * read acquitsition time information from sentinel 3 safe manifest
 * @param acquiTime - pointing to TiXmlElement
 */
void S3MetaData::readAcquiTime(TiXmlElement* acquiTime){
    TiXmlElement* e = acquiTime->FirstChildElement("sentinel-safe:startTime");
    startTime = std::string(e->FirstChild()->ToText()->ValueStr());
    e = acquiTime->FirstChildElement("sentinel-safe:stopTime");
    stopTime = std::string(e->FirstChild()->ToText()->ValueStr());
}

/**
 * read platform information from sentinel 3 safe manifest
 * @param platform
 */
void S3MetaData::readPlatformInfo(TiXmlElement* platform){
    TiXmlNode* node = platform->FirstChild();
    for (node; node; node = node->NextSibling()){
        if ( strcmp( node->Value(), "sentinel-safe:familyName" ) == 0 ){
            platformFamily = node->FirstChild()->ToText()->ValueStr();
        }
        else if ( strcmp( node->Value(), "sentinel-safe:number" ) == 0 ){
            platformNumber = node->FirstChild()->ToText()->ValueStr();
        }
        else if ( strcmp( node->Value(), "sentinel-safe:instrument" ) == 0 ){
            instrument = node->FirstChildElement("sentinel-safe:familyName")->Attribute("abbreviation");
        }
    }
}

/**
 * read general product information from sentinel 3 safe manifest
 * @param genPInfo
 */
void S3MetaData::readGenProdInfo(TiXmlElement* genPInfo){
    TiXmlNode* node;
    for (node = genPInfo->FirstChild(); node; node = node->NextSibling()){
        if ( strcmp( node->Value(), "sentinel3:productName" ) == 0 ){
            productName = node->FirstChild()->ToText()->ValueStr();
        }
        else if ( strcmp( node->Value(), "sentinel3:productType" ) == 0 ){
            productType = node->FirstChild()->ToText()->ValueStr();
        }
        else if ( strcmp( node->Value(), "sentinel3:timeliness" ) == 0 ){
            timeliness = node->FirstChild()->ToText()->ValueStr();
        }
        else if ( strcmp( node->Value(), "sentinel3:baselineCollection" ) == 0 ){
            baselineColl = node->FirstChild()->ToText()->ValueStr();
        }
        else if ( strcmp( node->Value(), "sentinel3:creationTime" ) == 0 ){
            creaTime = node->FirstChild()->ToText()->ValueStr();
        }
        else if ( strcmp( node->Value(), "sentinel3:productUnit" ) == 0 ){
            prodUnitType = node->FirstChildElement("sentinel3:type")->FirstChild()->ToText()->ValueStr();
            prodUnitDuration = node->FirstChildElement("sentinel3:duration")
                                   ->FirstChild()->ToText()->ValueStr();
            if (node->FirstChildElement("sentinel3:alongtrackCoordinate")){
                prodUnitAlongTrackCoord = node->FirstChildElement("sentinel3:alongtrackCoordinate")
                                          ->FirstChild()->ToText()->ValueStr();
            }
            else prodUnitAlongTrackCoord = std::string("0000");
        }
        else if ( strcmp( node->Value(), "sentinel-safe:instrument" ) == 0 ){
            instrument = node->FirstChildElement("sentinel-safe:familyName")->Attribute("abbreviation");
        }
    }
}

/**
 * read slstr product information from sentinel 3 safe manifest
 * @param slstrInfo
 */
void S3MetaData::readSlstrProdInfo(TiXmlElement* slstrInfo){
    TiXmlNode* node = slstrInfo->FirstChild();
    for (node; node; node = node->NextSibling()){
        if ( strcmp( node->Value(), "slstr:resolution" ) == 0){
            if ( strcmp( node->ToElement()->Attribute("grid"), "0.5 km stripe A" ) == 0 ) {
                std::string s = node->FirstChild("slstr:spatialResolution")
                                      ->FirstChild()->ToText()->ValueStr();
                slstrPInfo.nadirImg0500m.xRes = StringToNumber<int>(s);
                slstrPInfo.nadirImg0500m.yRes = slstrPInfo.nadirImg0500m.xRes;
                slstrPInfo.obliqImg0500m.xRes = slstrPInfo.nadirImg0500m.xRes;
                slstrPInfo.obliqImg0500m.yRes = slstrPInfo.nadirImg0500m.yRes;
            }
            else if ( strcmp( node->ToElement()->Attribute("grid"), "Tie Points" ) == 0 ) {
                std::string s = node->FirstChild("slstr:spatialResolution")
                                      ->FirstChild()->ToText()->ValueStr();
                slstrPInfo.nadirTpgImg.xRes = StringToNumber<int>(s);
                slstrPInfo.nadirTpgImg.yRes = 1000;
                slstrPInfo.obliqTpgImg.xRes = slstrPInfo.nadirTpgImg.xRes;
                slstrPInfo.obliqTpgImg.yRes = slstrPInfo.nadirTpgImg.yRes;
            }
        }
        else if ( strcmp( node->Value(), "slstr:nadirImageSize" ) == 0){
            if ( strcmp( node->ToElement()->Attribute("grid"), "0.5 km stripe A" ) == 0 ) {
                readImageInfo(node, &slstrPInfo.nadirImg0500m);
            }
            if ( strcmp( node->ToElement()->Attribute("grid"), "Tie Points" ) == 0 ) {
                readImageInfo(node, &slstrPInfo.nadirTpgImg);
            }
        }
        else if ( strcmp( node->Value(), "slstr:obliqueImageSize" ) == 0){
            if ( strcmp( node->ToElement()->Attribute("grid"), "0.5 km stripe A" ) == 0 ) {
                readImageInfo(node, &slstrPInfo.obliqImg0500m);
            }
            if ( strcmp( node->ToElement()->Attribute("grid"), "Tie Points" ) == 0 ) {
                readImageInfo(node, &slstrPInfo.obliqTpgImg);
            }
        }
    }    
}

/**
 * read slstr image information from sentinel 3 safe manifest
 * @param node
 * @param imgInfo
 */
void S3MetaData::readImageInfo(TiXmlNode* node, ImageProperties* imgInfo){
    std::string s = node->FirstChild("sentinel3:startOffset")
            ->FirstChild()->ToText()->ValueStr();
    imgInfo->yOff = StringToNumber<int>(s);
    s = node->FirstChild("sentinel3:trackOffset")
            ->FirstChild()->ToText()->ValueStr();
    imgInfo->xOff = StringToNumber<int>(s);
    s = node->FirstChild("sentinel3:rows")
            ->FirstChild()->ToText()->ValueStr();
    imgInfo->height = StringToNumber<int>(s);
    s = node->FirstChild("sentinel3:columns")
            ->FirstChild()->ToText()->ValueStr();
    imgInfo->width = StringToNumber<int>(s);
    
    imgInfo->nPix = imgInfo->width * imgInfo->height;
}

void S3MetaData::assertValidImgProp(const ImageProperties& imgProp){
    if (  (imgProp.width <= 0) || (imgProp.height <= 0)
            || (imgProp.xRes  <= 0) || (imgProp.yRes   <= 0)
            || (imgProp.xOff  < 0) || (imgProp.yOff   < 0) ){

        std::string msg(typeid(*this).name());
        msg.append(": Image properties could not be read from manifest!\n");
        throw std::domain_error(msg);
    }
}
