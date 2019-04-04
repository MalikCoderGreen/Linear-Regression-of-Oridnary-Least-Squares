
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
/*
 *Name: Maliak Green 
 *Course : COSC320
 *Description: Matrix operations program that allows user to input the dimensions for their desired matrix. Program allows user to enter in which operation they would like to perform if the dimensions match the operation conditions. 
 *
 *
 * 
 */

template<class T>
Matrix<T>::Matrix(){




}


// Non-default constructor.
template <class T>
Matrix<T>::Matrix(int  n_rows, int n_cols){

	rows = n_rows; 
	cols = n_cols;

	m = new T*[n_rows];
	for(int i = 0; i < n_rows; i++)
		m[i] = new T[n_cols];

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			m[i][j] = 0;

	//Can access by m[3][1].
	//Or access by *(*(m+2)).

}

template <class T>
// Copy constructor.
Matrix<T>::Matrix(const Matrix &right){

	rows = right.rows;
	cols = right.cols;
	m = new T*[rows];
	for(int j = 0; j < rows; j++){
		m[j] = new T[cols];
		for(int k = 0; k < cols; k++){
			m[j][k] = right.m[j][k];
		}}
}

// Operator overload for matrix addition.
template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& B){

	Matrix C(B.rows, B.cols);
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			C.m[i][j] = m[i][j] + B.m[i][j];	

	return C;
}

// Opeator overload for matrix multiplication.
template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix &B){ 

	Matrix C(rows, B.cols);
	for(int i = 0; i < C.rows; i++)
		for(int j = 0; j < C.cols; j++)
			for(int l = 0; l < cols; l++)
				C.m[i][j] += m[i][l] * B.m[l][j];		

	return C;

}

// Operator overload for matrix subtraction
template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix &B){

	Matrix C(rows, cols);
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			C.m[i][j] = m[i][j] - B.m[i][j];

	return C;


}


template <class T>
Matrix<T> Matrix<T>::operator=(const Matrix &B){

	if(this != &B){

		for(int i = 0; i < rows; i++)
			delete [] m[i];

		delete [] m;

		rows = B.rows;
		cols = B.cols;

		m = new T*[rows];
		for(int i = 0; i < rows; i++){
			m[i] = new T[cols];
			for(int j = 0; j < cols; j++){

				m[i][j] = B.m[i][j];
			}

		}

	}

	return *this;
}



// Function to make the matrix random values for nxm.
template <class T>
void Matrix<T>::giveRandom(){

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			m[i][j] = std::rand() % 9  + 1;		
}


// Function to make the matrix sorted in order ending with (nxm) - 1.
template <class T>
void Matrix<T>::giveSorted(){

	T sortedCounter = 1;

	for(int i = 0; i < rows; i++)
	{	for(int j = 0; j < cols; j++)
		{	m[i][j] = sortedCounter;
			sortedCounter++;
		}
	}
}

// Function to make the matrix backwards starting from (nxm) as the first entry.
template <class T>
void Matrix<T>::giveBackwards(){



	int  backwardsCounter = (rows*cols);
	for(int i = 0; i < rows; i++){
		{
			for(int j = 0; j < cols; j++)
				m[i][j] = backwardsCounter;
			backwardsCounter--;
		}
	}

}

// Function to make matrix full of duplicate values.
template <class T>
void Matrix<T>::giveDuplicate(){

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			m[i][j] = std::rand() % (rows*cols) / 3;



}


// Function to print the matrix A.
template <class T>
void Matrix<T>::printMatrix(){

		
	std::cout.precision(6);
	for(int i = 0; i < rows; i++)
	{	printf("%s", "|");
		for(int j = 0; j < cols; j++)
		{	if(j - (cols - 1))
			std::cout << std::showpoint << m[i][j] << " ";
			else{
				std::cout << std::showpoint<< m[i][j] <<"|";
			}
		}
		printf("%s", "\n");

	}
}





// Function to make Diagonal nxn Matrix.
template <class T>
void Matrix<T>::makeDiagonal(){

	int push = 0;
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){



			// Condition for A sub nn.	
			if(j  == push)
				m[i][j] = std::rand() % 9 + 1;


			// Every other element.
			else
				m[i][j] = 0; 

		}
		push = i + 1;
	}
}



// Function to make identity matrix.
template <class T>
void Matrix<T>::makeIdentity(){


	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			// Condition for A sub nn.	
			if(i == j)
				m[i][j] = 1;


			// Every other element.
			else
				m[i][j] = 0;

		}
	}


}

//Function to check if a given matrix is the indentity.
template<class T>
bool Matrix<T>::checkIdentity(){

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){

			if(i == j){
				if(m[i][j] == 1)
					continue;
				else
					return false;
			}

			else{
				if(m[i][j] == 0)
					continue;
				else
					return false;
			}

		}

	}

	return true;

	}





// Function to make upper-triangular matrix.
template <class T>
void Matrix<T>::upperTriangular(){

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){

			for(int k = 0; k < i; k++)
				m[i][k] = 0;

			m[i][j] = std::rand() % 9 + 1;

		}


	}




}

// Function to make lower-triangular matrix.
template <class T>
void Matrix<T>::lowerTriangular(){

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){

			m[i][j] = std::rand() % 9 + 1;

			for(int k = rows - 1; k > i; k--)
				m[i][k] = 0;		

		}
	}
}

// Function to take transpose of a matrix.
template <class T>
Matrix<T> Matrix<T>::Transpose(){

	Matrix ATranspose(cols, rows);

	for(int i = 0; i < ATranspose.rows; i++)
		for(int j = 0; j < ATranspose.cols; j++)
			ATranspose.m[i][j] = m[j][i];



	return ATranspose;

}

// Function to multiply matrix A by a scalar.
template <class T>
void Matrix<T>::scalarMultiply(T scalar){

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			m[i][j] *= scalar;

}


//Functionm to pad a matrix whose dimensions are not a power of 2.
template <class T>
Matrix<T> Matrix<T>:: padMatrix(){

	int k = pow(2, ceil(log2(rows))) - rows; // Obtain k for new (n+k) x (n+k) matrix.
	Matrix padded (rows+k, rows+k);

	for(int i = 0; i < padded.rows; i++){
		for(int j = 0; j < padded.cols; j++){
			if(i < rows && j < cols)
				padded.m[i][j] = m[i][j];


			else  if(i == j)
				padded.m[i][j] = 1; 
			else 
				padded.m[i][j] = 0;

		}

	}

	return padded;
}

// Function to check the power of a given matrix.
template <class T>
bool Matrix<T>::checkPower(int n){

	bool isPower = true;
	while(n != 1)
	{	
		if(n % 2 == 0){
			n = n / 2;
			continue;

		}
		else{
			isPower = false;
			break;
		}
	}


	return isPower;
}


// Function to check preprocess matrix A for inversePrivate,
// If the matrix's dimensions multiplied do not produce a power of 2, the matrix will be "padded",
template <class T>
Matrix<T> Matrix<T>::inverse(int &ops){



	if(!isSymmetrical()){
	
	Matrix <T> nonSym = checkPower(rows) ? *this : padMatrix();
	ops++;

	Matrix <T> nonSym_I = (nonSym.Transpose() * nonSym).inversePrivate(ops) * nonSym.Transpose();
	ops++;

	return nonSym_I.Slice(0, 0, rows - 1, rows - 1);
	}


	Matrix <T> tmp = checkPower(rows) ? *this : padMatrix();
	ops++;


	Matrix <T> tmp_I = tmp.inversePrivate(ops);
	ops++;

	return tmp_I.Slice(0, 0, rows - 1, rows - 1);
}

// Private function to take the inverse of a Matrix A.
template <class T>
Matrix<T> Matrix<T>::inversePrivate(int &ops){


	// Base case, A is a one element matrix.
	if(rows == 1){
		Matrix <T> A_I(1, 1);

		if(m[0][0] != 0){
			A_I.m[0][0] = 1 / m[0][0];

			return A_I;
		}

	}
	ops++;


	Matrix <T> B = Slice(0, 0, rows/2 - 1, rows/2 - 1); // Top-left.
	ops++;

	Matrix <T> C = Slice((rows/2) , 0, rows - 1, (rows/2) - 1); // Bottom - Left.
	ops++;

	Matrix <T> C_T = Slice(0, rows / 2, rows / 2 - 1, rows - 1); // Top - Right.
	ops++;


	Matrix <T> D = Slice((rows / 2), (rows / 2), rows - 1, rows  - 1); // Bottom - Right.
	ops++;




	// Series of Operations for algorithm,
	Matrix <T> B_I =  B.inversePrivate(ops); 
	ops++;

	Matrix <T> W = C * B_I;
	ops++;
	Matrix <T> W_T = B_I * C_T;
	ops++;

	Matrix <T> X = W * C_T;
	ops++;

	Matrix <T> S = D - X;
	ops++;

	Matrix <T> S_I = S.inversePrivate(ops);
	ops++;

	Matrix <T> V = S_I;
	ops++;

	Matrix <T> Y = S_I * W;
	ops++;

	Matrix <T> Y_T = Y.Transpose();
	ops++;

	Y_T.scalarMultiply(-1);
	ops++;

	Matrix <T> T_M = Y_T;
	ops++;

	Y.scalarMultiply(-1);
	ops++;

	Matrix <T> U = Y;
	ops++;

	Y.scalarMultiply(-1); // Reset back to positive for rest of operations.
	ops++;

	Matrix <T> Z = W_T * Y;
	ops++;

	Matrix <T> R = B_I + Z;
	ops++;



	Matrix <T> A_IN = Join_Matrices(R, T_M, U, V);
	ops++;






	return A_IN;


}

// Function to slice a tmp matrix from  Matrix A.
// Takes in the top corner of the matrix along with the bottom corner seperated into four parameters.
template <class T>
Matrix<T> Matrix<T>::Slice(int tr, int tc, int br, int bc){

	int t_rows = br - tr + 1;
	int t_cols = bc - tc + 1;

	Matrix <T> tmp(t_rows, t_cols);
	for(int i = 0; i < t_rows; i++)
		for(int j = 0; j < t_cols; j++)
			tmp.m[i][j] = m[i + tr][j + tc];

	return tmp;
}



// Function to join 4 n x n matrices into one.
// Takes in four matricies for the top-left, bottom-left, top-right, and bottom right portions of the resultant matrix.
template <class T>
Matrix<T> Matrix<T>::Join_Matrices(Matrix<T> top_L, Matrix<T> top_R, Matrix <T> B_L, Matrix <T> B_R){

	Matrix <T> A( 2 * top_L.rows, 2 * top_L.cols);

	for(int i = 0; i < A.rows / 2; i++)
		for(int j = 0; j < A.cols / 2; j++)
			A.m[i][j] = top_L.m[i][j];

	// Iterators for sub matrices.
	int k = 0 , l = 0;
	for(int i = 0; i < A.rows / 2; i++){
		for(int j = A.rows / 2; j < A.rows; j++){
			A.m[i][j] = top_R.m[k][l];
			l++;
		}
		l = 0;
		k++;
	}

	k = 0;
	l = 0;

	for(int i = A.rows / 2; i < A.rows; i++){
		for(int j = 0; j < A.rows / 2; j++){
			A.m[i][j] = B_L.m[k][l];
			l++;
		}
		l = 0;
		k++;
	}

	k = 0; 
	l = 0;

	for(int i = A.rows / 2; i < A.rows; i++){
		for(int j = A.rows / 2; j < A.rows; j++){
			A.m[i][j] = B_R.m[k][l];
			l++;
		}
		l = 0; 
		k++;
	}




	return A;

}


// Function set element for matrix.
// Takes in the element, and the respective positions for the matrix.
template <class T>
void Matrix<T>::setElement(int element, int pos1, int pos2){

	m[pos1][pos2] = element;
}

template <class T>
Matrix<T>::~Matrix(){

	for(int i = 0; i < rows; i++)
		delete [] m[i];

	delete [] m;
}




// Bonus function to determine if a matrix is non-singular or if isn't.
template <class T>
bool Matrix<T>::isNonSingular(){
	int ops;
	Matrix <T> A_I = inverse(ops);

	Matrix <T> iden = *this * A_I;

	if(iden.checkIdentity())
		return true;

	else
		return false;
}

// Bonus function to determine if an nxn matrix is symmetrical.
template<class T>
bool Matrix<T>::isSymmetrical(){


	Matrix <T> tmp = Transpose();

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			if(m[i][j] != tmp.m[i][j])
				return false;

	return true;
}


// Bonus function to implement Strassen's algortihm to perform matrix multiplication.
template<class T>
Matrix<T> Matrix<T>::StrassenAlg(Matrix<T> A, Matrix <T> B){
	
	if(A.rows == 1)
		return A * B;
	// 1.) Split operand matrices.

	Matrix <T> a11 = A.Slice(0, 0, rows/2 - 1, rows/2 - 1);
	Matrix <T> a12 = A.Slice(0, rows / 2, rows / 2 - 1, rows - 1);
	Matrix <T> a21 = A.Slice((rows/2) , 0, rows - 1, (rows/2) - 1);
	Matrix <T> a22 = A.Slice((rows / 2), (rows / 2), rows - 1, rows  - 1);

	Matrix <T> b11 = B.Slice(0, 0, rows/2 - 1, rows/2 - 1);
	Matrix <T> b12 = B.Slice(0, rows / 2, rows / 2 - 1, rows - 1);
	Matrix <T> b21 = B.Slice((rows/2) , 0, rows - 1, (rows/2) - 1);
	Matrix <T> b22 = B.Slice((rows / 2), (rows / 2), rows - 1, rows  - 1);

	
	// 2.) Create matrices.
	
	Matrix <T> S1 = b12 - b22;
	Matrix <T> S2 = a11 + a12;
	Matrix <T> S3 = a21 + a22;
	Matrix <T> S4 = b21 - b11;
	Matrix <T> S5 = a11 + a22;
	Matrix <T> S6 = b11 + b22;
	Matrix <T> S7 = a12 - a22;
	Matrix <T> S8 = b21 + b22;
	Matrix <T> S9 = a21 - a11;
	Matrix <T> S10 = b11 + b12;

	// 3.) Create recursion.
	
	Matrix <T> P1 = StrassenAlg(a11 , S1);
	Matrix <T> P2 = StrassenAlg(S2 , b22);
	Matrix <T> P3 = StrassenAlg(S3 , b11);
	Matrix <T> P4 = StrassenAlg(a22 , S4);
	Matrix <T> P5 = StrassenAlg(S5 , S6);
	Matrix <T> P6 = StrassenAlg(S7 , S8);
	Matrix <T> P7 = StrassenAlg(S9 , S10);

	// 4.) Result matrix C.
	


	
	Matrix <T> C11 = P5 + P4 - P2 + P6;
	Matrix <T> C12 = P1 + P2;

	Matrix <T> C21 = P3 + P4;
	Matrix <T> C22 = P5 + P1 - P3 - P7;

	Matrix <T> C = Join_Matrices(C11, C12, C21, C22);

	return C;
	
}




// Function to read points data into matrix.
// Takes in a file name as an argument.
template <class T>
Matrix<T> Matrix<T>::obtainData(std::string fileName, int &ops){

	std::ifstream points(fileName, std::ios::in);
	std::stringstream s; // For counting columns.
	std::string line, tmp;

	int fileRows = 0, fileCols = 0;

	s.clear();
	std::getline(points, line);
	s << line;

	// String stream will count amount of columns using first line of matrix.
	while(s >> tmp)	
		fileCols++;


	// Reset ifstream to beginning of file.
	points.clear();
	points.seekg(0);

	if(!points)
		std::cout << "File Not found! \n";


	// Count file rows.
	while(std::getline(points, line)){
		fileRows++;		
	}


	Matrix <T> A(fileRows, fileCols);
	Matrix <T> b(fileRows, 1);

	
	points.clear();
	points.seekg(0);


	int i = 0;

	while(!points.eof() && i < fileRows){
	
		for(int j = 0; j < fileCols; j++){
			points >> A.m[i][j];
		}
		i++;

		if(points.eof())
			break;	
	}

	// Give b vector the last column of A.
	for(int i = 0; i < A.rows; i++){

		b.m[i][0] = A.m[i][A.cols - 1];
		A.m[i][A.cols - 1] = 1;
	}


	return  (A.Transpose() * A).inverse(ops) * A.Transpose() * b;

}
