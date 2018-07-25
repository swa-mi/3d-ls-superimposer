//============================================================================
// Name        : Click_test-1.cpp
// Author      : parichit
// Version     : just a testing version
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cmath>
#include <limits>
#include <cstdio>
#include <numeric>
#include <iomanip>
#include "AlignInterface.h"
#include "vector_operation.h"
#include "jacobi_eigenvalue.hpp"
using namespace std;


void AlignInterface::Align_Cliques(vector<vector<double> > &store_str1repatom,
		vector<vector<double> > &store_str2repatom, double &rmsd) {

	vector_operation vo;

	if (store_str1repatom.size() != store_str2repatom.size()) {
		cout<< "unequal sizes of query and target matrices, please chesk that you are passing correct repAtoms and no duplicates\n" << endl;
		//vo.print_2Dvector(store_str1repatom);
		//cout << "\ntarget matrix" << endl;
		//vo.print_2Dvector(store_str2repatom);
		exit(1);
	}

	vector<double> querycentroid(3), targetrcentroid(3), translation(3),
			sum(3), diff(3);

	int str1_size = store_str1repatom.size();
	int str2_size = store_str2repatom.size();

	vector<vector<double> > rotation(3, vector<double> (3));
	vector<vector<double> > kmatrix(4, vector<double> (4));
	double egnvec[16];
	double kmat[16];
	double egnval[4];
	int it_num;
	int rot_num;
	vector<double> te;
	vector<vector<double> > dev(str1_size);

	double sum_x = 0, sum_y = 0, sum_z = 0;

	for (uint i = 0; i < store_str1repatom.size(); i++) {

		sum_x = sum_x + store_str1repatom[i][0];
		sum_y = sum_y + store_str1repatom[i][1];
		sum_z = sum_z + store_str1repatom[i][2];

	}

	querycentroid[0] = sum_x / str1_size;
	querycentroid[1] = sum_y / str1_size;
	querycentroid[2] = sum_z / str1_size;

	sum_x = 0;
	sum_y = 0;
	sum_z = 0;

	for (uint j = 0; j < store_str2repatom.size(); j++) {

		sum_x = sum_x + store_str2repatom[j][0];
		sum_y = sum_y + store_str2repatom[j][1];
		sum_z = sum_z + store_str2repatom[j][2];

	}

	targetrcentroid[0] = sum_x / str2_size;
	targetrcentroid[1] = sum_y / str2_size;
	targetrcentroid[2] = sum_z / str2_size;

	sum_x = 0;
	sum_y = 0;
	sum_z = 0;

	//cout << "query matrix is:" << endl;
	//vo.print_2Dvector(store_str1repatom);

	//cout << "Target matrix is:" << endl;
	//vo.print_2Dvector(store_str2repatom);

	//cout << "Query centroid\n" << endl;
	//vo.print_vector(querycentroid);
	//cout << "Target centroid\n" << endl;
	//vo.print_vector(targetrcentroid);

	translation[0] = targetrcentroid[0] - querycentroid[0];
	translation[1] = targetrcentroid[1] - querycentroid[1];
	translation[2] = targetrcentroid[2] - querycentroid[2];

	//cout << "translation vector is:\n" << endl;
	//vo.print_vector(translation);

	sum[0] = querycentroid[0] + targetrcentroid[0];
	sum[1] = querycentroid[1] + targetrcentroid[1];
	sum[2] = querycentroid[2] + targetrcentroid[2];
	diff[0] = querycentroid[0] - targetrcentroid[0];
	diff[1] = querycentroid[1] - targetrcentroid[1];
	diff[2] = querycentroid[2] - targetrcentroid[2];

	try {
		calculatekearsley(store_str1repatom, store_str2repatom, sum, diff,
				kmatrix);
		//cout << "Kearseley matrix is:\n" << endl;
		//vo.print_2Dvector(kmatrix);

	} catch (std::exception& msg) {
		cout << "The exception occured after kearsley calling: " << msg.what()
				<< endl;
	}
	/*
	mat temp(4, 4, fill::zeros);
	for (uint i = 0; i < kmatrix.size(); i++) {

		temp(i, 0) = kmatrix[i][0];
		temp(i, 1) = kmatrix[i][1];
		temp(i, 2) = kmatrix[i][2];
		temp(i, 3) = kmatrix[i][3];
	}
	*/

	kmat[0] = kmatrix[0][0];
	kmat[1] = kmatrix[0][1];
	kmat[2] = kmatrix[0][2];
	kmat[3] = kmatrix[0][3];
	kmat[4] = kmatrix[1][0];
	kmat[5] = kmatrix[1][1];
	kmat[6] = kmatrix[1][2];
	kmat[7] = kmatrix[1][3];
	kmat[8] = kmatrix[2][0];
	kmat[9] = kmatrix[2][1];
	kmat[10] = kmatrix[2][2];
	kmat[11] = kmatrix[2][3];
	kmat[12] = kmatrix[3][0];
	kmat[13] = kmatrix[3][1];
	kmat[14] = kmatrix[3][2];
	kmat[15] = kmatrix[3][3];

	/*
	eig_sym(eignval, eignvec, temp);
	vec temp1 = eignvec.col(0);

	cout << eignval << endl;
	*/

	try {
		jacobi_eigenvalue(4, kmat, 100, egnvec, egnval, it_num, rot_num);
		//cout << "Rotation are" << rot_num << "iterations are" << it_num << endl;

	} catch (exception &e) {
		e.what();
	}

	//cout << egnval[0] << " " << egnval[3] << endl;


	vector<double> t;
	t.push_back(egnvec[0]);t.push_back(egnvec[1]);t.push_back(egnvec[2]);t.push_back(egnvec[3]);

	rotationmatrix(t, rotation);

	//cout << "\n\nThe eigen values are as follows: \n\n" << eignval << endl;
	//cout << "The eigen vectors are as follows: \n\n" << eignvec << endl;
	//cout << "rotation matrix is:\n" << endl;
	//vo.print_2Dvector(rotation);

	/*
	 vo.print_2Dvector(store_str1repatom);
	 cout << "\n" << endl;
	 vo.print_2Dvector(store_str2repatom);
	 cout << "\n" << endl;
	 */

	/*
	 cout << "real coordinates are:\n" << endl;
	 vo.print_2Dvector(store_str2repatom);
	 */

	//cout << store_str1repatom.size() << endl;

	rotateandtranslatecoords(store_str1repatom, rotation, querycentroid,
			translation);

	/*
	 cout << "rotated coordinates are:\n" << endl;
	 vo.print_2Dvector(store_str1repatom);
	 */

	//dev.clear();
	find_deviation(store_str1repatom, store_str2repatom, dev);

	//cout << "Size of the deviation vector is: " << dev.size() << dev[0].size() << endl;

	/*
	 cout << "Deviation matrix is:\n" << endl;
	 vo.print_2Dvector(dev);
	 */

	if (egnval[0] < egnval[3] * numeric_limits<double>::epsilon()) {

		//double rmsd = findrmsd(dev);
		//double rmsd2 = sqrt(rmsd / store_str1repatom.size());
		//cout << "The epsilon value is: \n" << numeric_limits<double>::epsilon()
		//	<< endl;
		//cout << "The RMSD less than epsilon: \n" << rmsd << endl;
		rmsd = 0.0;
		//vo.print_2Dvector(store_str1repatom);

	}

	else {
		rmsd = findrmsd(dev);
		//if(rmsd<1){
                //cout << "\nThe RMSD is: " << rmsd  << endl;
                //"\nTransformed coordinates are: \n"<< endl;
                //vo.print_2Dvector(store_str1repatom);
		//}
	}

}

