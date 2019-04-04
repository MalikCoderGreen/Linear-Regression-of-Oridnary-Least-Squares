#ifndef MATRIX_H
#define MATRIX_H
#include <stdio.h>
#include <string>

template <class T>
class Matrix
{

	private:

		T **m;
		int rows;
		int cols;

	public:

		Matrix();

		Matrix(int ,int);
		
		Matrix(const Matrix &);
	
		Matrix<T> operator+(const Matrix&);
	

		Matrix<T> operator*(const Matrix&); 

		Matrix<T> operator-(const Matrix&);
		
		Matrix<T> operator=(const Matrix&);

	
		void giveRandom();
		void giveSorted();
		void giveBackwards();
		void giveDuplicate();
		void makeIdentity();
		bool checkIdentity();
		void makeDiagonal();
		void upperTriangular();
		void lowerTriangular();
		Matrix<T> Transpose();
		void scalarMultiply(T);
		void setElement(int, int, int);
		Matrix<T> inverse(int&);
		Matrix<T> inversePrivate(int&);
		Matrix<T> padMatrix();
		Matrix<T> Slice(int, int, int, int);
		Matrix<T> Join_Matrices(Matrix<T>, Matrix<T>, Matrix<T>, Matrix<T>);
		bool checkPower(int);
		void printMatrix();
		bool isNonSingular();
		bool isSymmetrical();
		Matrix<T> StrassenAlg(Matrix <T>, Matrix <T>);
		Matrix<T> obtainData(std::string, int&);

		


		~Matrix();
	
};

#include "matrix.cpp"

#endif 

