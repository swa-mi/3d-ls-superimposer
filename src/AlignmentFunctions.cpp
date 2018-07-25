/*
 * AlignmentFunctions.cpp
 *
 *  Created on: 17-Jun-2015
 *      Author: parichit
 */
#include <iomanip>
#include <math.h>
#include "AlignInterface.h"
#include "vector_operation.h"
using namespace std;


void AlignInterface::calculatekearsley(vector<vector<double> > &querycoords,
		vector<vector<double> > &targetcoords, vector<double> &sum,
		vector<double> &diff, vector<vector<double> > &kmatrix) {

	vector<double> addm(3), diffm(3);
	vector<double> addvec(3), diffvec(3);

	double xadd = 0, yadd = 0, zadd = 0, xdiff = 0, ydiff = 0, zdiff = 0;

	//Starting the creation of kearsley's matrix
	//Actual element calculation start here

	try {

		//cout <<"The addition and difference matrix are as follows:" << addm<< "\n"<< diffm<< end

		for (uint i = 0; i < querycoords.size(); i++) {

			addm[0] = querycoords[i][0] + targetcoords[i][0];
			addm[1] = querycoords[i][1] + targetcoords[i][1];
			addm[2] = querycoords[i][2] + targetcoords[i][2];

			diffm[0] = targetcoords[i][0] - querycoords[i][0];
			diffm[1] = targetcoords[i][1] - querycoords[i][1];
			diffm[2] = targetcoords[i][2] - querycoords[i][2];

			addvec[0] = addm[0] - sum[0];
			addvec[1] = addm[1] - sum[1];
			addvec[2] = addm[2] - sum[2];

			diffvec[0] = diffm[0] + diff[0];
			diffvec[1] = diffm[1] + diff[1];
			diffvec[2] = diffm[2] + diff[2];

			xdiff = diffvec[0];
			ydiff = diffvec[1];
			zdiff = diffvec[2];

			xadd = addvec[0];
			yadd = addvec[1];
			zadd = addvec[2];

			kmatrix[0][0] += xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
			kmatrix[0][1] += yadd * zdiff - ydiff * zadd;
			kmatrix[0][2] += xdiff * zadd - xadd * zdiff;
			kmatrix[0][3] += xadd * ydiff - xdiff * yadd;

			kmatrix[1][1] += yadd * yadd + zadd * zadd + xdiff * xdiff;
			kmatrix[1][2] += xdiff * ydiff - xadd * yadd;
			kmatrix[1][3] += xdiff * zdiff - xadd * zadd;

			kmatrix[2][2] += xadd * xadd + zadd * zadd + ydiff * ydiff;
			kmatrix[2][3] += ydiff * zdiff - yadd * zadd;

			kmatrix[3][3] += xadd * xadd + yadd * yadd + zdiff * zdiff;

			addm.clear();
			diffm.clear();
			addvec.clear();
			diffvec.clear();

		}

		kmatrix[1][0] += kmatrix[0][1];
		kmatrix[2][0] += kmatrix[0][2];
		kmatrix[2][1] += kmatrix[1][2];
		kmatrix[3][0] += kmatrix[0][3];
		kmatrix[3][1] += kmatrix[1][3];
		kmatrix[3][2] += kmatrix[2][3];

	} catch (std::exception& msg) {

		cout << "\n";
		cout << "--------------------------------\n";
		cout
				<< "EXCEPTION --> function: Kearsleymatrix file name: Alignment Functions"
				<< msg.what() << endl;
		cout << msg.what() << "\n";
		cout << "--------------------------------\n";
		cout << "\n";

	}

}

void AlignInterface::rotationmatrix(vector<double> &eignvecarr,
		vector<vector<double> > &rotatate) {

	double q1 = 0.0;
	double q2 = 0.0;
	double q3 = 0.0;
	double q4 = 0.0;

	try {

		/*
		 q1 = eignvecarr(0);
		 q2 = eignvecarr(1);
		 q3 = eignvecarr(2);
		 q4 = eignvecarr(3);
		 */

		q1 = eignvecarr[0];
		q2 = eignvecarr[1];
		q3 = eignvecarr[2];
		q4 = eignvecarr[3];

		//Starting the procedure for initializing the rotation matrix

		rotatate[0][0] = (q1 * q1 + q2 * q2 - q3 * q3 - q4 * q4);
		rotatate[0][1] = 2.0 * (q2 * q3 + q1 * q4);
		rotatate[0][2] = 2.0 * (q2 * q4 - q1 * q3);
		rotatate[1][0] = 2.0 * (q2 * q3 - q1 * q4);
		rotatate[1][1] = (q1 * q1 + q3 * q3 - q2 * q2 - q4 * q4);
		rotatate[1][2] = 2.0 * (q3 * q4 + q1 * q2);
		rotatate[2][0] = 2.0 * (q2 * q4 + q1 * q3);
		rotatate[2][1] = 2.0 * (q3 * q4 - q1 * q2);
		rotatate[2][2] = (q1 * q1 + q4 * q4 - q2 * q2 - q3 * q3);
		//return rotationmatrix;

	} catch (std::exception& e) {

		cout << "\n";
		cout << "--------------------------------\n";
		cout
				<< "EXCEPTION --> function: rotationmatrix file name: Alignment Functions"
				<< e.what() << endl;
		cout << e.what() << "\n";
		cout << "--------------------------------\n";
		cout << "\n";
	}

}

void AlignInterface::rotateandtranslatecoords(
		vector<vector<double> > &coordinatematrix,
		vector<vector<double> > &rotation, vector<double> &centroid,
		vector<double> &translation) {

	vector<double> temp;
	temp.clear();

	try {

		if (coordinatematrix.size() == 0) {

			cout << "\n";
			cout << "--------------------------------\n";
			cout
					<< "Invalid dimensions of the Vector: function name: rotatecoords file name:  AlignmentFunction"
					<< endl;
			cout << "--------------------------------\n";
			cout << "\n";

		} else {

			for (uint i = 0; i < coordinatematrix.size(); i++) {

				//rotate
				temp.push_back(0);
				temp.push_back(0);
				temp.push_back(0);

				coordinatematrix[i][0] = coordinatematrix[i][0] - centroid[0];
				coordinatematrix[i][1] = coordinatematrix[i][1] - centroid[1];
				coordinatematrix[i][2] = coordinatematrix[i][2] - centroid[2];

				temp[0] = ((rotation[0][0] * coordinatematrix[i][0])
						+ (rotation[0][1] * coordinatematrix[i][1])
						+ (rotation[0][2] * coordinatematrix[i][2]));
				temp[1] = ((rotation[1][0] * coordinatematrix[i][0])
						+ (rotation[1][1] * coordinatematrix[i][1])
						+ (rotation[1][2] * coordinatematrix[i][2]));
				temp[2] = ((rotation[2][0] * coordinatematrix[i][0])
						+ (rotation[2][1] * coordinatematrix[i][1])
						+ (rotation[2][2] * coordinatematrix[i][2]));

				//cout << "change:1" << endl;
				//vo.print_vector(temp);

				temp[0] = temp[0] + centroid[0] + translation[0];
				temp[1] = temp[1] + centroid[1] + translation[1];
				temp[2] = temp[2] + centroid[2] + translation[2];

				//cout << "change:2" << endl;
				//vo.print_vector(temp);

				//translate
				coordinatematrix[i][0] = temp[0];
				coordinatematrix[i][1] = temp[1];
				coordinatematrix[i][2] = temp[2];

				/*
				 cout << "changed coordinates" << endl;
				 cout << setw(10) << right << coordinatematrix[i][0] << setw(10)
				 << right << coordinatematrix[i][1] << setw(10) << right
				 << coordinatematrix[i][2] << endl;
				 */

				temp.clear();
			}

		}

	} catch (exception& e) {
		cout << "\n";
		cout << "--------------------------------\n";
		cout
				<< "EXCEPTION --> function: rotate and translatecoords file name: Alignment Functions"
				<< e.what() << endl;
		cout << "--------------------------------\n";
		cout << "\n";
	}

}

void AlignInterface::find_deviation(vector<vector<double> > &str1,
		vector<vector<double> > &str2, vector<vector<double> > &dev) {
	//dev.clear();

	for (uint i = 0; i < str2.size(); i++) {

		dev[i].push_back(str2[i][0] - str1[i][0]);
		dev[i].push_back(str2[i][1] - str1[i][1]);
		dev[i].push_back(str2[i][2] - str1[i][2]);

	}

}

double AlignInterface::findrmsd(vector<vector<double> > &deviationvector) {

	vector_operation vo;
	vector<double> rms;
	double rmsd = 0.0, rms_x = 0.0, rms_y = 0.0, rms_z = 0.0;

	try {

		if (deviationvector.size() == 0) {

			cout << "\n";
			cout << "--------------------------------\n";
			cout
					<< "Invalid dimensions of the Vector: function name: findrmsd file name:  AlignmentFunction"
					<< endl;
			cout << "--------------------------------\n";
			cout << "\n";
		} else {

			for (uint i = 0; i < deviationvector.size(); i++) {

				rms_x += pow(deviationvector[i][0], 2);
				rms_y += pow(deviationvector[i][1], 2);
				rms_z += pow(deviationvector[i][2], 2);
			}

			rmsd = rmsd + rms_x + rms_y + rms_z;
			rmsd = rmsd / deviationvector.size();
			//vo.print_2Dvector(deviationvector);
			//cout << "rmsd without squareroot is:" << rmsd << endl;
			rmsd = sqrt(rmsd);

		}

		return rmsd;
	} catch (exception& e) {
		cout << "\n";
		cout << "--------------------------------\n";
		cout
				<< "Exception --> function name: findrmsd file name: Alignment Functions"
				<< endl;
		cout << "--------------------------------\n";
		cout << "\n";

	}
	return rmsd;
}

void AlignInterface::find_so(vector<vector<double> >&str1,
		vector<vector<double> >&str2, float &so, float sup_cutoff) {
    vector_operation vo;


	if (str1.size() != str2.size()) {
		cout
				<< "Matrix size is not equal in SO calculation, exiting the program."
				<< endl;
		exit(0);
	}

	double x1 = 0, x2 = 0, y1 = 0, y2 = 0, z1 = 0, z2 = 0;
	double distance = 0;
	float counter = 0;

	//vo.print_2Dvector(str1);
	for (uint i = 0; i < str1.size(); i++) {

		x1 = str1[i][0];
		x2 = str2[i][0];
		y1 = str1[i][1];
		y2 = str2[i][1];
		z1 = str1[i][2];
		z2 = str2[i][2];

		x1 = x2 - x1;
		y1 = y2 - y1;
		z1 = z2 - z1;

		distance = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
        //cout << "sup_cutoff is: " << sup_cutoff << " and distance is: " << distance <<  endl;

		if (distance < sup_cutoff) {
			counter++;
		//	cout << "counter is: " << counter << endl;
		}
	}

	so = (counter / str1.size()) * 100;
	//cout << "Counter ,  SO, Size  is: " << counter << ", " << so << ", " << str1.size() << endl;
}

