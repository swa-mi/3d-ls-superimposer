/*
 * vector_operation.h
 *
 *  Created on: 24-Oct-2015
 *      Author: neelesh
 */

#ifndef VECTOR_OPERATION_H_
#define VECTOR_OPERATION_H_

#include <vector>
#include <iostream>
#include <cstdlib>

using namespace std;

class vector_operation {

public:    
    
/*template <class T1,class T2> void set_vector(T1 &V, T2 x, T2 y, T2 z);

template <class T> void print_vector(T V);

template <class T> void print_2Dvector(T V);

template <class T> T get_minmax_from_vector(T v);*/


/*
 * Set vector T1 with values in vector T2
 */
template <class T1, class T2>
void set_vector(T1 &V, T2 x, T2 y, T2 z) {
    V[0] = x;
    V[1] = y;
    V[2] = z;
}

/*
 * print the single dimension vector T
 */
template <class T>
void print_vector(T V) {
    if (V.size() == 0) {
        cout << "Empty vector";
    }
    for (uint i = 0; i < V.size(); i++) {
        cout << ", " << V[i];
    }
    cout << endl;
}

/*
 * Print the 2D vector
 */
template <class T>
void print_2Dvector(T V) {

    if (V.size() == 0) {
        cout << "Empty vector";
    }
    cout << endl;
    for (uint i = 0; i < V.size(); i++) {
        cout << "Vector " << i << ": ";
        //print_vector(V[i]);
        for(uint j=0;j<V[i].size();j++){
        	cout << V[i][j] << ",";

        }
        cout << endl;

    }
    //cout << endl;
}

/*
 * Get the minimum and maximum elements in the vector
 */
template <class T>T get_minmax_from_vector(T v) {
    double max = -99999999;
    double min = 99999999;

    for (uint i = 0; i < v.size(); i++) {
        if (v[i] > max) {
            max = v[i];
        }
        if (v[i] < min) {
            min = v[i];
        }
    }
    vector<double> minmax;
    minmax.push_back(min);
    minmax.push_back(max);
    return minmax;
}

};

#endif /* VECTOR_OPERATION_H_ */




