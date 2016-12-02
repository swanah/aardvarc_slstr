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
#include <functional>



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
