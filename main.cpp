#include <iostream>
#include <stdio.h>
#include "matrix.h"
#include <fstream>

/*
 *Name: Maliak Green 
 *Course : COSC320
 *Description: Matrix operations program that allows user to input the dimensions for their desired matrix. Program allows user to enter in which operation they would like to perform if the dimensions match the operation conditions. 
 *
 *
 * 
 */


// Read in data from points100.txt file.<



int main(int arg, char **argv){


	srand(time(NULL));
	int operations = 0;

	printf("%s", "***Linear Regression of Oridanry Least Squares***\n\n");
	std::cout << "Based on the file " << argv[1] << " here is the solution.\n\n";
	Matrix <float> A (2,1);
	Matrix <float> tmp(2,1);

	A = tmp.obtainData(argv[1], operations);


	std::cout << "Optimal vector Beta: \n";
	A.printMatrix();

	std::cout << "\n\nNumber of operations to obtain answer: "  << operations <<  std::endl;


	return 0;
}
