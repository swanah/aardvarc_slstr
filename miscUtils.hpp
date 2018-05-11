/* 
 * File:   miscUtils.h
 * Author: akheckel
 *
 * Created on 22. November 2016, 18:49
 */

#ifndef MISCUTILS_HPP
#define MISCUTILS_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>



void replaceStringInPlace(std::string& subject, const std::string& search, const std::string& replace);
void trimLine(std::string& s);
void splitLine(const std::string& s, std::vector<std::string>& words);
double timer();
void setBit(unsigned short *flag, unsigned short bitvalue, char test);
void setBit(short *flag, short bitvalue, char test);
double intAng(double a, double b, double w);
bool fileExists(const std::string& s);
bool isDir(const std::string& s);
bool createDirsUnix(const std::string& s);
std::string getCurrentTimeStr();
std::string convS3DateToCci(const std::string& s3DateTime);
double parseTimeToSec1970(const std::string& s3DateTime);
bool isIntegerStr(const std::string& s);
float calcScatAng(float sza, float vza, float raz);


template<typename T>
T StringToNumber(const std::string& numberAsString)
{
   T valor;

   std::stringstream stream(numberAsString);
   stream >> valor;
   if (stream.fail()) {
      std::runtime_error e(numberAsString);
      throw e;
   }
   return valor;
}    


#endif /* MISCUTILS_HPP */

