/*
 * Main.cpp
 *
 *  Created on: 2017-10-01
 *      Author: Swastik
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>
#include <iomanip>
#include "ctime"

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>

//#include "jacobi_eigenvalue.cpp"
#include "jacobi_eigenvalue.hpp"
//#include "Align_cliques.cpp"
#include "AlignInterface.h"
//#include "AlignmentFunctions.cpp"
#include "vector_operation.h"
//#include "PDBParser.cpp"
#include "CliqParser.h"

using namespace std;
using namespace boost;

inline bool exists_test1 (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}

int main(int argc, char** argv) {


    vector<vector<double> > store_str1repatom; // this will store coordinates for query
    vector<vector<double> > store_str2repatom; // this will store coordinates for target inside the for loop (see below)
    store_str1repatom.clear();
    store_str2repatom.clear();


    vector_operation vo;
    AlignInterface ai;
    

    int clique_size ;
    float sup_cutoff ;
    float penalty ;
    string query_dir;
    string target_dir;


    

    try
    {
        namespace po = boost::program_options;

        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("input-cliqfile,i", po::value<string>(), "path to .cliq file")
        ("query-dir,q", po::value<string>(&query_dir)->required(), "Path to dir where query should be extracted from")
        ("target-dir,t", po::value<string>(&target_dir)->required(), "Path to dir where targets should be extracted from")
        ("size,n", po::value<int>(&clique_size)->default_value(8), "Set size of clique")
        ("rmsd,r", po::value<float>(&sup_cutoff)->default_value(2.0), "Set rmsd cutoff for succesful superimposition/match")
        ("penalty,p", po::value<float>(&penalty)->default_value(10.0), "Set penalty for non-superimposition")
/*
        ("repatom,a", po::value<string>(&repatom)->required(), "representative atoms list, comma separated. e.g. r1,r2,r5,r10")
        ("target_pdbs", po::value<string>(&target_pdbs)->required(), "representative atoms list, comma separated. e.g. r1,r2,r5,r10")
        ("strctre1_file", po::value<string>(&strctre1_file)->required(), "representative atoms list, comma separated. e.g. r1,r2,r5,r10")

*/        ;

        po::variables_map vm; 
        try 
        { 
          po::store(po::parse_command_line(argc, argv, desc),  
                    vm); // can throw 
     
          /** --help option 
           */ 
          if ( vm.count("help")  ) 
          { 
            std::cout << "Superimposition code for cliques" << std::endl 
                      << desc << std::endl; 
            return 0; 
          } 
     
          po::notify(vm); // throws on error, so do after help in case 
                          // there are any problems 
        } 
        catch(po::error& e) 
        { 
          std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
          std::cerr << desc << std::endl; 
          return 1; 
        }

        //start main function actions here


        double rmsd = penalty;
        double rmsd2 = penalty;
        float so=0.0;
        float so2=0.0;

        std::string line;                                // string to store each line while reading stdin or .cliq file
        std::vector<string> vec_atomstatements;     // vector to store atom statements from gpdb file
        std::string query_clique;;                    //first line of .cliq file 
        std::vector<string> targetcliqs;            //rest of the lines of .cliq file

        int linecount = 0;

        // mutex to protect file access (shared across threads)
        static std::mutex mutex;

        // lock mutex before accessing file
        std::lock_guard<std::mutex> lock(mutex);


        if( vm.count("input-cliqfile"))
        {
            ifstream inputcliqfile(vm["input-cliqfile"].as<std::string>());
            if( inputcliqfile.is_open() ){
                while( getline(inputcliqfile, line) ){
                    //cout << line << endl;
                    if (linecount==0){
                        query_clique = query_dir+"/"+line;
                    } // end if block checking for first line or not
                    else{
                        targetcliqs.push_back(target_dir+"/"+line);
                    } // end else block 

                    linecount = linecount+1;
                } // end while block looping over lines in .cliq filestream object
            } // end if block checking if file is open
        } // end if block checking if input-cliqfile was passed as argument
        else { // else read from stdin - actually means if stdin is being redirected
            if ( isatty (0) ){ // there was no stdin redirection
                throw std::invalid_argument("ERROR: STDIN unused and input-cliqfile argument empty.\nSee --help for usage.\nExiting.\n\n");
            }
            else{
                while( getline(cin, line) ){
                    //cout << line << endl;
                    if (line.empty()){
                        throw std::invalid_argument("ERROR: STDIN is empty. Exiting.\n\n");
                    }
                    else{ // stdin isn't empty
                        if (linecount==0){
                            query_clique = query_dir+"/"+line;
                        }
                        else{
                            targetcliqs.push_back(target_dir+"/"+line);
                        }
                    }
                    linecount = linecount+1;
                } // end looping over filestream object
            } // end checking if stdin was redirected
        } // end reading .cliq file

        /* At this point query_clique is the first line of .cliq file
         and targetcliqs is vector<string> with the rest of the lines
        */

        /**if targetcliqs vector is empty, it means that no clique with same composition 
        was found in the database and hence, we should abort, giving this match a penalty for non-matching
        */
        if ( !targetcliqs.empty() ){

            // QUERY CLIQUE PROCESSING

            //Read query_clique ATOM statements into a vector
            std::map<int, std::string> q_map_atomstatements;          //  map of line number to atom statements in gpdb file, for whole file
            std::vector<string> q_vec_atomstatements;
            // q_vec_atomstatements is going to be of the form:
            /*
            ["ATOM      7 r1   GLN H   3      38.621  -27.34   1.484        2.71",
            "ATOM      8 r2   GLN H   3      40.562  -25.93   1.889        2.85",
            ...
            "ATOM     20 r3   ARG H   7      39.363 -31.798   4.613        3.49"]
            */
            std::vector<string> q_cliquevector;

            boost::split(q_cliquevector, query_clique, boost::is_any_of(" \t"), boost::token_compress_on);
            //cliquevector now looks like [ ./2mta     2021    2020    2019    2010    2005    2009    2004    2003]
            string q_gpdbfilename = q_cliquevector.front() + ".gpdb";

            if ( !exists_test1(q_gpdbfilename) ){
                throw std::runtime_error("ERROR: query gpdb file not found with path " + q_gpdbfilename + "\n\n");
            }

            std::ifstream q_gpdbif(q_gpdbfilename, std::ifstream::in);
            if (q_gpdbif.is_open()){
                Scan_Cliq_for_Cliqs(q_cliquevector, q_vec_atomstatements, q_map_atomstatements, q_gpdbif);
                //vo.print_vector(q_vec_atomstatements);
            } // end if block checking for open q_gpdbif filestream object

            //q_rep_atom

            // Read in coordinates and other variables of query vector
            std::vector<pair<string, vector<double>> > q_vecmap_strcoordinates;    //  atomdetails:coordinates pairs
            std::map<string, vector<double> > q_map_strcoordinates;               //   the map of --do--
            Scan_ATOM_for_coord(q_vec_atomstatements,          // vector with ATOM statements for members of clique
                                q_vecmap_strcoordinates,
                                q_map_strcoordinates);     //  atomdetails: coordinates pairs
        
            //check how q_vecmap_strcoordinates looks like:                    
            /*for (auto iter1 : q_vecmap_strcoordinates){
                cout << iter1.first;
                vo.print_vector(iter1.second);
            }*/
            /*q_vecmap_strcoordinates = 
            dict{
                 7 r1   GLN H   3  : 38.621, -27.34, 1.484
                 8 r2   GLN H   3  : 40.562, -25.93, 1.889
                11 r1   ALA H   4  : 35.974, -28.632, 2.453
                12 r8   ALA H   4  : 36.566, -27.78, 4.584
                17 r1   ARG H   7  : 36.603, -32.429, -0.936
                18 r2   ARG H   7  : 38.866, -31.969, -0.025
                19 r2   ARG H   7  : 38.541, -31.49, 1.388
                20 r3   ARG H   7  : 39.363, -31.798, 4.613
            }

            */

            std::vector<pair<string, vector<double>> >::iterator queryiterator;    //to iterate over all entries of vector
            string keys1[clique_size];                       //to store all keys in the pair entries of the vector
            int it;                                         //to store index of iterator
            for (queryiterator=q_vecmap_strcoordinates.begin(); queryiterator!= q_vecmap_strcoordinates.end(); queryiterator++) {
                it=std::distance(q_vecmap_strcoordinates.begin(), queryiterator);
                keys1[it]=queryiterator->first;

            }

           
            std::sort(keys1+1, keys1+clique_size);  //sort to reach lexicographically lowest permutation

            std::vector<std::vector<string>> allpermutvecs;
            std::vector<string> permutedvec;

            do {
                    for (int i=0; i<clique_size; i++){
                        permutedvec.push_back(keys1[i]);
                    }
                    allpermutvecs.push_back(permutedvec);
                    permutedvec.clear();
                } while( std::next_permutation( (keys1+1), keys1+clique_size));
            /**allpermutvecs now has vectors such that
            each vector is a permutation of the atoms in query clique
            however only points except the first point(central point) are permuted, so # of permutations is greatly reduced
            permute only last n-1 atoms/chm-groups
            */

            std::vector< pair<double, std::vector<string> > > allsuperimps;
            //all superimpositions to be stored as a pair of rmsd and vector of strings of atom-details
            /** loop over each string in targetcliqs vector
                for each string in the vector, find coordinates and superimpose on q_vecmap_strcoordinates
                find rmsd and append rmsd,details to a vector
            */

            //Read query_clique ATOM statements into a vector
            std::map<int, std::string> t_map_atomstatements;          //  map of line number to atom statements in gpdb file, for whole file
            std::vector<string> t_vec_atomstatements;
            // t_vec_atomstatements is going to be of the form:
            /*
            ["ATOM      7 r1   GLN H   3      38.621  -27.34   1.484        2.71",
            "ATOM      8 r2   GLN H   3      40.562  -25.93   1.889        2.85",
            ...
            "ATOM     20 r3   ARG H   7      39.363 -31.798   4.613        3.49"]
            */

            std::vector<string> t_cliquevector;
            string t_gpdbfilename;
            // START LOOPING OVER TARGETS
            for (auto target_clique : targetcliqs){
                // target_clique is of type std::string
                boost::split(t_cliquevector, target_clique, boost::is_any_of(" \t"), boost::token_compress_on);
                //cliquevector now looks like [ ./2mta     2021    2020    2019    2010    2005    2009    2004    2003]
                t_gpdbfilename = t_cliquevector.front() + ".gpdb";
                if ( !exists_test1(t_gpdbfilename) ){
                    throw std::runtime_error("ERROR: target gpdb file not found with path " + t_gpdbfilename + "\n\n");
                }

                std::ifstream t_gpdbif(t_gpdbfilename, std::ifstream::in);
                if (t_gpdbif.is_open()){
                    Scan_Cliq_for_Cliqs(t_cliquevector, t_vec_atomstatements, t_map_atomstatements, t_gpdbif);
                    //vo.print_vector(t_vec_atomstatements);
                } // end if block checking for open t_gpdbif filestream object
                //q_rep_atom == t_rep_atom?

                // Read in coordinates and other variables of target vector
                std::vector<pair<string, vector<double>> > t_vecmap_strcoordinates;    //  atomdetails:coordinates pairs
                std::map<string, vector<double> > t_map_strcoordinates;               //the map of --do--
                
                Scan_ATOM_for_coord(t_vec_atomstatements,          // vector with ATOM statements for members of clique
                                    t_vecmap_strcoordinates,
                                    t_map_strcoordinates);     //  atomdetails: coordinates pairs
            
                //check how t_vecmap_strcoordinates looks like:                    
                /*for (auto iter1 : t_vecmap_strcoordinates){
                    cout << iter1.first;
                    vo.print_vector(iter1.second);
                }*/

                //superimpose and find all rmsds for all permutations in all target_pdbs
                string compoarr2;       //to store composition of chm-groups/atom-names for later comparison
                string da_atomname;     //temp atom-name string var

                //place all the atom details (basically columns 12:26)
                //from the keys in vecmap_str1coordinates in
                //keys1 which is an array of size= clique_size
                for (auto targetcoordpair : t_vecmap_strcoordinates){
                    store_str2repatom.push_back(targetcoordpair.second);
                    da_atomname = targetcoordpair.first.substr(6,4);
                    trim(da_atomname);
                    compoarr2 += da_atomname;
                    //cout << targetcoordpair.first << ":" ;
                    //vo.print_vector(targetcoordpair.second);
                }

                string permutarr;       //temp string for the composition of strctre-1 in one of the permutations
                string atom_name;       //temp string for atom_name to be added to permutarr
                std::vector<double> v;  //stores coordinates of strctre-1 in each permutation
                std::vector<string> temporaryvec1;
                string atomids="";
                string tempatomid="";
                for (unsigned int i=0; i < allpermutvecs.size(); i++){

                    for (int j=0; j<clique_size; j++){
                        atom_name=allpermutvecs[i][j];
                        tempatomid = atom_name.substr(0,5);
                        trim(tempatomid);
                        atomids += tempatomid;
                        atom_name=atom_name.substr(6,4);
                        trim(atom_name);
                        permutarr += atom_name;
                        tempatomid.clear();
                    }
                    //cout << "|" << permutarr << "|" << compoarr2 << endl;

                    if (permutarr==compoarr2){
                        //cout << "yay " <<  t_gpdbfilename <<" | " <<permutarr << endl;
                        for (int k=0; k<clique_size; k++){
                            v= q_map_strcoordinates[allpermutvecs[i][k]];
                            store_str1repatom.push_back(v);
                            temporaryvec1.push_back(allpermutvecs[i][k]);
                        }
                        // superimpose str2 and str1, and str1 on str2
                        // take the best of the two superimpositions, since A on B != B on A

                        // since str1 gets transformed after superimposition,
                        // use a copy of str1 and str2 for the second superimposition
                        std::vector<std::vector<double> > store_str1repatom_dup(store_str1repatom);
                        std::vector<std::vector<double> > store_str2repatom_dup(store_str2repatom);

                        //vo.print_2Dvector(store_str1repatom);
                        //vo.print_2Dvector(store_str2repatom);

                        ai.Align_Cliques(store_str1repatom, store_str2repatom, rmsd);
                        ai.find_so(store_str1repatom, store_str2repatom, so, sup_cutoff);

                        ai.Align_Cliques(store_str2repatom_dup, store_str1repatom_dup, rmsd2);
                        ai.find_so(store_str2repatom_dup, store_str1repatom_dup, so2, sup_cutoff);

                        // if either or both of so and so2 are 100, then
                        // add the rmsd:temporaryvec1 pair or rmsd2:temporaryvec1 pair respectively to allsuperimps
                        // later sort for the best rmsd out of those and output that rmsd and temporaryvec1
                        if (so==100){
                            //cout << "|" << permutarr << "|" << endl;
                            //temporaryvec1.push_back(std::to_string(so));
                            temporaryvec1.push_back(t_gpdbfilename);
                            //cout << rmsd <<" | " << rmsd2 << endl;
                            //cout << rmsd << endl;
                            //vo.print_2Dvector(temporaryvec1);
                            //cout << " so 100 " << t_gpdbfilename << endl;
                            allsuperimps.push_back(make_pair(rmsd, temporaryvec1));
                        }
                        if (so2==100){
                            //cout << "|" << permutarr << "|" << endl;
                            //temporaryvec1.push_back(std::to_string(so));
                            temporaryvec1.push_back(t_gpdbfilename);
                            //cout << rmsd <<" | " << rmsd2 << endl;
                            //cout << rmsd << endl;
                            //vo.print_2Dvector(temporaryvec1);
                            //cout << " so 100 " << t_gpdbfilename << endl;
                            allsuperimps.push_back(make_pair(rmsd2, temporaryvec1));
                        }


                        /*else {
                            cout << "so not 100pc " << endl;
                            cout << t_gpdbfilename << endl;

                        }*/
                        temporaryvec1.clear();
                    }
                    store_str1repatom.clear();
                    permutarr.clear();
                    atomids.clear();
                    v.clear();

                } //end for loop for superimposition against all permutations
                t_vecmap_strcoordinates.clear();
                t_map_strcoordinates.clear();
                t_map_atomstatements.clear();
                t_vec_atomstatements.clear();
                store_str2repatom.clear();
                compoarr2.clear();
                target_clique.clear();
                t_cliquevector.clear();
            } // end for loop looping over all strings in targetcliqs vector
        
            //display allsuperimps which is a vector of pairs
            //such that the key in each pair is the rmsd of the superimp and the value
            //is SO+"#"+atomids
            /*std::vector<pair<double, std::vector<string> >>::iterator superimpiter;
            for (superimpiter=allsuperimps.begin(); superimpiter!=allsuperimps.end(); superimpiter++){
                cout << "later" << superimpiter->first << endl;
                vo.print_vector(superimpiter->second);
            }*/

            string best_match;
            if (allsuperimps.empty()){
                rmsd = penalty;
                best_match = "No matches found";
                //cout << "No matches found" << endl;
            }
            else{
                std::sort(allsuperimps.begin(), allsuperimps.end());
                best_match= allsuperimps[0].second[clique_size];
                cout << "Best match is on: \n" << best_match << endl << endl;
                rmsd = (allsuperimps[0]).first;

                for (int j=0; j<clique_size; j++){
                    cout << allsuperimps[0].second[j] << endl;
                }

            }

            // take filenames of both and concat them then add .clique : this is the outfilename
            /*string output_file = strctre1_file.substr(0, strctre1_file.find_last_of("/")) + "/"
                                    + strctre1_file.substr(strctre1_file.find_last_of("/") + 1, strctre1_file.length())
                                    + "-"
                                    + strctre2_file.substr(strctre2_file.find_last_of("/") + 1,strctre2_file.length())
                                    + ".clique";
                                    */
            //display rmsd
            cout << "\nRMSD is: \n" << rmsd << endl;

            //uncomment below block to print out best_match details


            /*Printing the transformed coordinates, un-comment below line to print*/
            //vo.print_2Dvector(store_str1repatom);

            /* uncomment below code block to print output everytime */
            /*clock_t end = clock();

            cout << "\nElapsefrom ‘std::map<std::__cxx11::basic_string<char>, std::vector<double> >’ to ‘std::__cxx11::string {aka std::__cxx11::basic_string<char>}d Time(seconds):" << setprecision(3) << double(
                    end - start) / CLOCKS_PER_SEC
                    << "\nElapsed Time(milliseconds):" << setprecision(3)
                    << ((double) (end - start) / CLOCKS_PER_SEC) * 1000 << endl;*/
        }
        else{ // if targetcliqs is empty
            rmsd = penalty;
            cout << "\nRMSD is: \n" << rmsd << endl;
        }
    } catch (std::exception &e) {

        cout << e.what();
        return EXIT_FAILURE;
    }

    return 0;
}


