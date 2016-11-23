/* 
 * File:   InputParameter.cpp
 * Author: akheckel
 * 
 * Created on 23. November 2016, 10:39
 */

#include "InputParameter.h"


using std::cout;
using std::endl;

        
InputParameter::InputParameter() {
    slstrProductDir = "e:/sat/S3A_SLSTR/S3A_SL_1_RBT____20160902T072518_20160902T072818_20160903T162446_0179_008_163_2879_LN2_O_NT_002.SEN3";
    s3MetaData.parseManifest(slstrProductDir);
    aodOutDir = "e:/sat/S3A_SLSTR_AOD/%YYYY%/%MM%";
    modifyAodOutNames();
}

/**
 * applies s3MetaData info to output dir and name
 */
void InputParameter::modifyAodOutNames(){
    std::string yyyy = s3MetaData.startTime.substr(0, 4);
    std::string mm = s3MetaData.startTime.substr(5, 2);
    std::string dd = s3MetaData.startTime.substr(8, 2);

    replaceStringInPlace(aodOutDir, "%YYYY%", yyyy);
    replaceStringInPlace(aodOutDir, "%MM%", mm);
    replaceStringInPlace(aodOutDir, "%DD%", dd);
    
    aodOutName = s3MetaData.productName;
    int extIdx = aodOutName.find_last_of(".") + 1;
    int extLen = aodOutName.length() - extIdx;
    aodOutName.replace(extIdx, extLen, "nc");
}

void InputParameter::replaceStringInPlace(std::string& subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    if (search.empty()) return;
    while((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}