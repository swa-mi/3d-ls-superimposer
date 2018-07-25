/*
 * CliqParser.h
 *
 * Created on: 23-Jul-2018
 * Author: Swastik
 */

#include <string>
#include "vector"
#include "bits/stdc++.h"

#ifndef CLIQPARSER_H_
#define CLIQPARSER_H_


void Scan_ATOM_for_coord(const std::vector<std::string> & vec_atomstatements,          // vector with ATOM statements for members of clique
    std::vector<std::pair<std::string, std::vector<double>> > &vecmap_strcoordinates,  //  atomdetails: coordinates pairs
    std::map<std::string, std::vector<double> > &map_str1coordinates                          //the map of --do--
    )

;
void Scan_Cliq_for_Cliqs(const std::vector<std::string> cliquevector,                // the filename of the .cliq file to be read
    std::vector<std::string>& vec_atomstatements,             // map of line number to atom statements in gpdb file, for clique only
    std::map<int, std::string> &map_atomstatements,          //  map of line number to atom statements in gpdb file, for whole file
    std::ifstream& gpdbif
    )
;


#endif /* CLIQPARSER_H_ */
