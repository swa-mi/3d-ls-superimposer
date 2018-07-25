/*
 * AlignInterface.cpp

 *
 *  Created on: 01-Jul-2015
 *      Author: parichit
 */

#include <vector>
using namespace std;


#ifndef ALIGNINTERFACE_H_
#define ALIGNINTERFACE_H_

class AlignInterface{

public:

void calculatekearsley(vector<vector<double> > &, vector<vector<double> > &, vector<double> &, vector<double> &, vector<vector<double> > &);

void rotationmatrix(vector<double> &, vector<vector<double> > &);

void rotateandtranslatecoords(vector<vector<double> > &, vector<vector<double> > &, vector<double> &, vector<double> &);

void find_deviation(vector<vector<double> > &, vector<vector<double> > &, vector<vector<double> > &);

double findrmsd(vector<vector<double> > &);

void Align_Cliques(vector<vector<double> > &, vector<vector<double> > &, double &);

void find_so(vector<vector<double> > &, vector<vector<double> > &, float &, float);

};

#endif /* ALIGNINTERFACE_H_ */


