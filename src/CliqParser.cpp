/*
 * CliqParser.cpp
 *
 * Created on: 22-Jul-2018
 * Author: Swastik
 */

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <string>
#include <float.h>
#include <iomanip>
#include <bits/stdc++.h>
#include <boost/algorithm/string.hpp>
#include "vector_operation.h"

using namespace std;
using namespace boost;

void Scan_ATOM_for_coord(const std::vector<std::string> & vec_atomstatements,          // vector with ATOM statements for members of clique
    std::vector<std::pair<std::string, std::vector<double>> > &vecmap_strcoordinates,  //  atomdetails: coordinates pairs
    std::map<string, vector<double> > &map_strcoordinates                          //the map of --do--
    ){


    vecmap_strcoordinates.clear();
    map_strcoordinates.clear();
/*
Here's a story
Jeff was a cricket
Unfortunately
Jeff was an addict
When the bug spray came
Jeff was high
Jeff couldn't hear them
And Jeff fucking died 
*/

    //vector_operation vo;            // vector operation functions
    try {
        string tempx, tempy, tempz;
        double tmpx, tmpy, tmpz;
        std::vector<double> atom_container;
        //cout << "start for loop" << endl;
        for (auto iter : vec_atomstatements) {
            //loop over vec_atomstatements and use each string (iter) to find relevant data

            tempx = iter.substr(30, 8);
            tempy = iter.substr(38, 8);
            tempz = iter.substr(46, 8);

            tmpx = std::stod(tempx.c_str());
            tmpy = std::stod(tempy.c_str());
            tmpz = std::stod(tempz.c_str());

            string fulldetail = iter.substr(6,20);

            atom_container.push_back(tmpx);
            atom_container.push_back(tmpy);
            atom_container.push_back(tmpz);

            vecmap_strcoordinates.push_back(
                    make_pair(fulldetail,
                            atom_container));
            map_strcoordinates.insert(
                        std::pair<string, std::vector<double>>
                        (fulldetail, atom_container));


            //free memory
            atom_container.clear();
            fulldetail.clear();
        }
        //cout << "end for loop" << endl;
        /*
        for (auto testiter: vecmap_strcoordinates){
            cout << testiter.first << ":";
            vo.print_vector(testiter.second);
            cout << endl;
        }*/


    } // end try block
    catch (std::exception &e) {

        cerr << "\n##########################\n";
        cerr
                << "\nThe exception happened in scanning function for coordinates\n\n"
                << e.what();
        cerr << "\n###########################\n";
    } // end catch block


} // end function


 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Scan_Cliq_for_Cliqs(const std::vector<std::string> cliquevector,                // the filename of the .cliq file to be read
    std::vector<std::string>& vec_atomstatements,             // map of line number to atom statements in gpdb file, for clique only
    std::map<int, std::string> &map_atomstatements,          //  map of line number to atom statements in gpdb file, for whole file
    std::ifstream& gpdbif
    ){


/* 
Takes .cliq file and .gpdb file as input.
Creates a vec-map (int to vector; vector of clique member group indices) of cliques
Creates a vec-map (int to string; string is the line of ATOM statement mapped to int which is the line number,
    starting with zero, in the gpdb file where the clique was found)

For each clique:
    push_back into a 2d_vector, vector containing the /// decide later

*/

    try {
        // mutex to protect file access (shared across threads)
        static std::mutex mutex;

        // lock mutex before accessing file
        std::lock_guard<std::mutex> lock(mutex);


        //vector_operation vo;            // vector operation functions

        //cliquevector looks like [ ./2mta     2021    2020    2019    2010    2005    2009    2004    2003]
        
        int linecount = 0;      // to keep track of line number. .cliq file has line number based clique indexing
        std::string atomstatement;
        //std::string::iterator cliquememberiterator; 
        int cliquemembervalue;

        if (gpdbif.is_open()){ // begin if block for gpdb reading if gpdb operation
            while (std::getline(gpdbif, atomstatement)) {
                if (atomstatement.substr(0,4) == "ATOM") {
                    map_atomstatements.insert(std::pair<int, std::string>(linecount, atomstatement));
                } // end if checking for ATOM statements
                linecount = linecount + 1;
            } // end looping over gpdb file

            for (auto cliquememberiterator=cliquevector.begin()+1; cliquememberiterator != cliquevector.end(); cliquememberiterator++){
                cliquemembervalue = std::stoi(*cliquememberiterator);
                vec_atomstatements.push_back(map_atomstatements[cliquemembervalue]);
            }
        
        } // end if block for gpdb reading if gpdb open
    } // end try block
    catch (std::exception &e) {

        cerr << "\n##########################\n";
        cerr
                << "\nThe exception happened in file scanning function for ATOM statements\n\n"
                << e.what();
        cerr << "\n###########################\n";
    } // end catch block

} // end function block




