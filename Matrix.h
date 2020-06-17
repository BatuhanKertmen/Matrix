#ifndef _AdvancedMath_H
#define _AdvancedMath_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <cstdarg>

class Vector
{
private:
	
	int Dimension;
	float* vector;

public:
	
	Vector();
	Vector(int dimension, ...);
	Vector(int dimension, bool automatic);
	Vector(const Vector& rhs);

	void print() const;


	~Vector();
};


class Matrix
{
private:
	
	int row;
	int column;

	float* matrix;

private:

	float Multiply(const Matrix& rhs, int GivenRow, int GivenCol) const;
	Matrix LowerMatrix(int GivenRow, int GivenColumn) const;
	
	Matrix& EchalonizerDown(int Row, int Col);
	Matrix& EchalonizerUp(int Row, int Col);
	
	// Elementary Row Operations
	Matrix& RowSwap(int Row1, int Row2);
	Matrix& RowMultiplication(float Multiplier, int Row);
	Matrix& RowAddition(int Row1, float Multiplier, int Row2);

public:	

	// Constructors
	Matrix(int Row, int Column, bool Automatic); // empty matrix
	Matrix(int witdh, int length); // user created matrix
	Matrix(const Matrix &mat); // Deep Copy Constructor
	explicit Matrix(int size = 0); //identity matrix constructor

	// Matrix Operations
	Matrix& operator= (const Matrix& rhs);
	Matrix operator+ (const Matrix& rhs) const;
	Matrix operator- (const Matrix& rhs) const;
	Matrix operator* (const Matrix& rhs) const;
	Matrix operator* (const float rhs) const;


	Matrix transpose() const;
	float det() const;
	bool isInvertible() const;
	Matrix inverse() const;
	Matrix ref() const;
	Matrix rref() const;

	void print() const;
	

	~Matrix();
};



#endif 

