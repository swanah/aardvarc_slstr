/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   S3NcdfImg.h
 * Author: akheckel
 *
 * Created on 23. November 2016, 20:43
 */

#ifndef S3NCDFIMG_H
#define S3NCDFIMG_H

#include "S3BasicImage.h"

class S3NcdfImg : public S3BasicImage {
public:
    S3NcdfImg(){}
    S3NcdfImg(const int& newWidth, const int& newHeight) 
        : S3BasicImage(newWidth, newHeight){};
    S3NcdfImg(const int& newWidth, const int& newHeight, const int& xOffset, const int& yOffset) 
        : S3BasicImage(newWidth, newHeight, xOffset, yOffset){};
    ~S3NcdfImg(){}
    
    void readNcdf(std::string ncdfName, std::string varName);
    
private:

};

#endif /* S3NCDFIMG_H */

