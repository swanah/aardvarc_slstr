/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   S3BasicImage.h
 * Author: akheckel
 *
 * Created on 23. November 2016, 19:21
 */

#ifndef S3BASICIMAGE_H
#define S3BASICIMAGE_H

#include <cstddef>


class S3BasicImage {
public:
    S3BasicImage(){ img = NULL; }
    //S3BasicImage(const int& newWidth, const int& newHeight){ S3BasicImage(newWidth, newHeight, 0, 0); }
    S3BasicImage(const int& newWidth, const int& newHeight, const int& xOffset = 0, const int& yOffset = 0);
    S3BasicImage(const S3BasicImage& orig);
    S3BasicImage& operator=(const S3BasicImage& rhs);

    virtual ~S3BasicImage();
    
    int width, height;
    int xOff, yOff;
    double valScale, valOffset;
    short *img;
    short noData;
    
private:

};

#endif /* S3BASICIMAGE_H */

