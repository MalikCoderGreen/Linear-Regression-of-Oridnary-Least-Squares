README 

Maliak Green
COSC 320

Project 1: Linear Regression With Ordinary Least Squares

This program uses my previous Lab - 4 with a matrix.cpp and matrix.h with more methods implemented in order to carry out the project task. The program takes command line arguments for reading in a data file into a matrix i.e "./ols points100.txt". I have included a set-Element function to explicity make a symmetrical matrix and test to see if my inverse function is working properly. Multiplying the symmetrical matrix by its inverse verifies this by producing the identity matrix. All tests for bonus are tested in bonusTest.cpp.

**Bonus Portion**: 

-  (For 5 points) I have implemented The non-singular bonus to determine if a matrix is singular or non-singular. The function takes the inverse of a matrix; if it produces the identity then it will return true. I also had to make a function "checkIdentity()" to verify that the resultant is actually the identity. 

- (For 10 points) I have implemented a function to determine if a matrix is symmetrical. I also have updated the inverse function to take the inverse of any nxn matrix including non-symmetrical matrices.

- (For 20 points) I have implemented Strassen's Algorithm. A function that takes two matrices as input and returns a new matrix.




