#include "AdvancedMath.h"


Vector::Vector(): Dimension(0), vector(nullptr) {}

Vector::Vector(int num, ...) : Dimension(num), vector(new float[num])
{
	va_list arguments;
	va_start(arguments, num);

	for (int i = 0; i < num; i++)
	{
		vector[i] = va_arg(arguments, float);
		vector[i] = 0;
	}
	va_end(arguments);
}

Vector::Vector(int dimension, bool automatic) : Dimension(dimension), vector(new float[dimension])
{
	for (int i = 0; i < dimension; i++)
	{
		vector[i] = 0;
	}
}

Vector::Vector(const Vector& rhs): Dimension(rhs.Dimension), vector(new float[rhs.Dimension])
{
	for (int i = 0; i < Dimension; i++)
	{
		vector[i] = rhs.vector[i];
	}
}

void Vector::print() const
{
	std::cout << "(";
	for (int i = 0; i < Dimension - 1; i++)
	{
		std::cout << std::setw(4) << vector[i] << ", ";
	}
	std::cout << std::setw(4) << vector[Dimension - 1] << ")";
}

Vector::~Vector()
{
	delete[] vector;
}


////////////////////////////////////////////////////////////////////
//															      //
//                BEGINNING OF MATRIX CLASS                       // 
//                                                                //
////////////////////////////////////////////////////////////////////


// creates an empty object but only allocates memory 
Matrix::Matrix(int row, int column, bool automatic) : row(row), column(column), matrix(new float[row * column]) {}

// creates a matrix with user entered inputs
Matrix::Matrix(int width, int length) : row(width), column(length), matrix(nullptr)
{
	if (width < 0 || length < 0)
	{
		std::cout << "Dimensions of a matrix can not be a negative number! \n";
		row = 0;
		column = 0;
		return;
	}
	
	if (width == 0 || length == 0)
	{
		row = 0;
		column = 0;
		return;
	}

	matrix = new float[width * length];
	std::cout << "Please enter the entries of the matrix \n";
	float dummy;

	for (int i = 0; i < (row * column); i++)
	{
		std::cin >> dummy;
		matrix[i] = dummy;
	}
}

// Identity matrix constructor
Matrix::Matrix(int size) : row(size), column(size), matrix(nullptr)
{
	if (size < 0)
	{
		std::cout << "Dimensions of a matrix can not be a negative number! \n";
		row = 0;
		column = 0;
		return;
	}
	if (size == 0)
	{
		return;
	}

	matrix = new float[size * size];
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			if (i == j)
			{
				matrix[j + column * i] = 1;
			}
			else
			{
				matrix[j + column * i] = 0;
			}
		}
	}
}

// Deep copy constructor
Matrix::Matrix(const Matrix& mat) : row(mat.row), column(mat.column)
{
	matrix = new float[row * column];
	for (int i = 0; i < (row * column); i++)
	{
		matrix[i] = mat.matrix[i];
	}

}

Matrix Matrix::operator+(const Matrix& rhs) const
{
	if (column != rhs.column || row != rhs.row)
	{
		std::cout << "Matrices with different dimensions can not be summed! \n";
		return Matrix(0,0);
	}

	Matrix result(*this);
	for (int i = 0; i < (row * column); i++)
	{
		result.matrix[i] += rhs.matrix[i];
	}

	return result;
}

Matrix Matrix::operator-(const Matrix& rhs) const
{
	if (column != rhs.column || row != rhs.row)
	{
		std::cout << "Matrices with different dimensions can not be subtracted! \n";
		return Matrix(0, 0);
	}

	Matrix result(*this);
	for (int i = 0; i < (row * column); i++)
	{
		result.matrix[i] = matrix[i] - rhs.matrix[i];
	}

	return result;
}

Matrix& Matrix::operator=(const Matrix& rhs)
{
	if (matrix == rhs.matrix)
	{
		return *this;
	}
	delete[] matrix;

	row = rhs.row;
	column = rhs.column;
	matrix = new float[rhs.row * rhs.column];
	for (int i = 0; i < (row * column); i++)
	{
		matrix[i] = rhs.matrix[i];
	}
	return *this;
}

Matrix Matrix::operator*(const Matrix& rhs) const
{
	Matrix result(row, rhs.column, true);

	for (int r = 0; r < result.row; r++)
	{
		for (int col = 0; col < result.column; col++)
		{
			result.matrix[col + r * column] = this->Multiply(rhs, r + 1, col + 1);
		}
	}
	return result;
}

Matrix Matrix::operator*(const float rhs) const
{
	Matrix result(row, column, true);

	for (int r = 0; r < result.row; r++)
	{
		for (int col = 0; col < result.column; col++)
		{
			result.matrix[col + r * column] = matrix[col + r * column] * rhs;
		}
	}
	return result;
}

Matrix Matrix::transpose() const
{
	Matrix result(column, row, true);
	
	for (int col = 0; col < column; col++)
	{
		for (int r = 0; r < row; r++)
		{
			result.matrix[r + col * result.column] = matrix[col + r * column];
		}
	}

	return result;
}

void Matrix::print() const
{
	std::cout << "\n";
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			std::cout << std::setw(4) << matrix[j + column * i] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

float Matrix::det() const
{
	if (column != row)
	{
		std::cout << "Determinant can only be used on square matrices! \n";
		return INT_MAX;
	}

	if (column == 1)
	{
		return matrix[0];
	}
	
	float result = 0;
	for (int i = 0; i < column; i++)
	{
		if (((i + 1) % 2) == 1)	result += matrix[i] * this->LowerMatrix(1, i + 1).det();
		else result -= matrix[i] * this->LowerMatrix(1, i + 1).det();
	}
	return result;
}

bool Matrix::isInvertible() const
{
	if (this->det() == 0) return false;
	else return true;
}

Matrix Matrix::inverse() const
{
	if (column != row)
	{
		std::cout << "Only square matrices are invertible! \n";
		return Matrix(0);
	}
	if (!this->isInvertible())
	{
		std::cout << "This matrix is not invertible!\n";
		return Matrix(0);
	}
	
	Matrix result(row, column, true);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			result.matrix[j + i * column] = this->LowerMatrix(i + 1, j + 1).det();
		}
	}

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			if (result.matrix[j + i * column] != 0 && ((j + i) % 2) == 1)
			{
				result.matrix[j + i * column] *= -1;
			}
		}
	}

	result = result.transpose();

	result = result * (1 / this->det());

	return result;
}

Matrix Matrix::ref() const
{
	Matrix result(*this);
	int PivotCounter = 0;
	int Col = 1;

	while (Col <= column && PivotCounter < row)
	{
		bool IsEchonable = false;
		for (int i = PivotCounter; i < row; i++)
		{
			if (result.matrix[(Col - 1) + i * column] != 0)
			{
				IsEchonable = true;
				break;
			}
		}

		if (IsEchonable)
		{
			result = result.EchalonizerDown(PivotCounter + 1, Col);
			PivotCounter++;
		}
		Col++;
	}

	return result;
}

Matrix Matrix::rref() const
{
	Matrix result(*this);
	int PivotCounter = 0;
	int Col = 1;

	while (Col <= column && PivotCounter < row)
	{
		bool IsEchonable = false;
		for (int i = PivotCounter; i < row; i++)
			if (result.matrix[(Col - 1) + i * column] != 0)
				IsEchonable = true;

		if (IsEchonable)
		{
			result = result.EchalonizerDown(PivotCounter + 1, Col);
			result = result.EchalonizerUp(PivotCounter + 1, Col);
			PivotCounter++;
		}
		Col++;
	}

	return result;
}

Matrix::~Matrix()
{
	delete[] matrix;
}

float Matrix::Multiply(const Matrix& rhs, int GivenRow, int GivenCol) const
{
	if (column != rhs.row)
	{
		std::cout << "These matrices can not be multiplied! \n";
		return 0;
	}

	float result = 0;
	for (int i = 0; i < column; i++)
	{
		result += (matrix[i + (column * (GivenRow - 1))] * rhs.matrix[(i * column) + (GivenCol - 1)]);
	}

	return result;
}

Matrix Matrix::LowerMatrix(int GivenRow, int GivenCol) const
{
	Matrix result(row - 1, column - 1, true);
	int RowCounter = 0, ColCounter = 0;

	for (int i = 0; i < row; i++)
	{
		ColCounter = 0;
		if (i == GivenRow - 1)
		{
			continue;
		}
			for (int j = 0; j < column; j++)
			{
				if (j == GivenCol - 1)
				{
					continue;
				}
				result.matrix[ColCounter + RowCounter * result.column] = matrix[j + i * column];

				ColCounter++;
			}
		RowCounter++;
	}

	return result;
}

Matrix& Matrix::RowSwap(int Row1, int Row2)
{
	std::vector<float> SwapRow(column);
	for (int i = 0; i < column; i++)
	{
		SwapRow[i] = matrix[i + column * (Row1 - 1)];
	}

	for (int i = 0; i < column; i++)
	{
		matrix[i + column * (Row1 - 1)] = matrix[i + column * (Row2 - 1)];
		matrix[i + column * (Row2 - 1)] = SwapRow[i];
	}

	return *this;
}

Matrix& Matrix::RowMultiplication(float Multiplier, int Row)
{
	for (int i = 0; i < column; i++)
	{
		if (matrix[i + column * (Row - 1)] != 0)
		{
			matrix[i + column * (Row - 1)] *= Multiplier;
		}
	}

	return *this;
}

Matrix& Matrix::RowAddition(int Row1, float Multiplier, int Row2)
{
	for (int i = 0; i < column; i++)
	{
		matrix[i + column * (Row1 - 1)] += (Multiplier * matrix[i + column * (Row2 - 1)]);
	}

	return *this;
}

Matrix& Matrix::EchalonizerDown(int Row, int Col)
{
	if (matrix[(Col - 1) + column * (Row - 1)] == 0)
	{
		int tempRow = Row;
		while (matrix[(Col - 1) + column * (tempRow - 1)] == 0)
		{
			tempRow++;
		}
		*this = this->RowSwap(Row, tempRow);
	}

	float pivot = matrix[(Col - 1) + column * (Row - 1)];
	*this = this->RowMultiplication(1 / pivot, Row);

	for (int i = (Row - 1) + 1; i < row; i++)
	{
		*this = this->RowAddition(i + 1, (-1) * matrix[(Col - 1) + column * i], Row);
	}

	return *this;
}

Matrix& Matrix::EchalonizerUp(int Row, int Col)
{
	if (matrix[(Col - 1) + column * (Row - 1)] == 0)
	{
		int tempRow = Row;
		while (matrix[(Col - 1) + column * (tempRow - 1)] == 0)
		{
			tempRow++;
		}
		*this = this->RowSwap(Row, tempRow);
	}

	float pivot = matrix[(Col - 1) + column * (Row - 1)];
	*this = this->RowMultiplication(1 / pivot, Row);

	for (int i = (Row - 1) - 1; i >= 0; i--)
	{
		*this = this->RowAddition(i + 1, (-1) * matrix[(Col - 1) + column * i], Row);
	}

	return *this;
}
