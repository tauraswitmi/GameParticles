//---------------------------------------------------------------
// Copyright 2024, Ed Keenan, all rights reserved.
//---------------------------------------------------------------

#ifndef MATRIX_H
#define MATRIX_H

// includes
#include "Enum.h"
#include "Vect4D.h" // Added

// forward declare
class Vect4D;

// class
class Matrix
{
public:

	// local enumerations
	enum class MatrixRow
	{
		MATRIX_ROW_0,
		MATRIX_ROW_1,
		MATRIX_ROW_2,
		MATRIX_ROW_3
	};

	Matrix();	
	Matrix( const Matrix& t );	
	Matrix& operator = (const Matrix& t) = default;
	~Matrix();

	void set( MatrixRow row, Vect4D &t );				//88 was *
	void get( MatrixRow row, Vect4D &vOut ) const;		//88 was *

	void setDefault();
	void setIdentMatrix();
	void setTransMatrix(Vect4D &t );				//88 was *
	void setScaleMatrix(Vect4D &s );				//88 was *
	void setRotZMatrix( float Z_Radians );

	float &operator[]( Index e);
	Matrix operator*( const float s ) const;

	float Determinant() const;
	
	Matrix GetAdjugate() const;

	Matrix& operator*=(const Matrix& t);
	Matrix& operator/=(const float t);
	
	void Inverse( Matrix &out ) const;

	
private:
	
	union
	{
		struct
		{
			Vect4D v0;
			Vect4D v1;
			Vect4D v2;
			Vect4D v3;
		};
	
		struct
		{
			// ROW 0
			float m0;
			float m1;
			float m2;
			float m3;

			// ROW 1
			float m4;
			float m5;
			float m6;
			float m7;

			// ROW 2
			float m8;
			float m9;
			float m10;
			float m11;

			// ROW 3
			float m12;
			float m13;
			float m14;
			float m15;

		};
	};

	friend struct MxM;
	friend struct MxM3;
	friend struct MxM4;
	friend struct MxM5;
};

struct MxM
{
	const Matrix& m1;
	const Matrix& m2;

	MxM(const Matrix& t1, const Matrix& t2)
		:m1(t1), m2(t2) {};

	operator Matrix()
	{
		Matrix tmp;

		// ROW 0 //

		__m128 _x = _mm_set_ps1(m1.v0.x);		// (x0, x0, x0, x0)
		__m128 _y = _mm_set_ps1(m1.v0.y);		// (y0, y0, y0, y0)
		__m128 _z = _mm_set_ps1(m1.v0.z);		// (z0, z0, z0, z0)
		__m128 _w = _mm_set_ps1(m1.v0.w);		// (w0, w0, w0, w0)

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;


		// ROW 1 //

		_x = _mm_set_ps1(m1.v1.x);
		_y = _mm_set_ps1(m1.v1.y);
		_z = _mm_set_ps1(m1.v1.z);
		_w = _mm_set_ps1(m1.v1.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(m1.v2.x);
		_y = _mm_set_ps1(m1.v2.y);
		_z = _mm_set_ps1(m1.v2.z);
		_w = _mm_set_ps1(m1.v2.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(m1.v3.x);
		_y = _mm_set_ps1(m1.v3.y);
		_z = _mm_set_ps1(m1.v3.z);
		_w = _mm_set_ps1(m1.v3.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;


		return tmp;
	}

};

struct MxM3
{
	const Matrix& m1;
	const Matrix& m2;
	const Matrix& m3;

	MxM3(const MxM& t1, const Matrix& t2)
		:m1(t1.m1), m2(t1.m2), m3(t2) {};

	operator Matrix()
	{
		Matrix tmp;

		// (Matrix) * (Matrix) // ===============================================//

		// ROW 0 //

		__m128 _x = _mm_set_ps1(m1.v0.x);	
		__m128 _y = _mm_set_ps1(m1.v0.y);	
		__m128 _z = _mm_set_ps1(m1.v0.z);	
		__m128 _w = _mm_set_ps1(m1.v0.w);	

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;


		// ROW 1 //

		_x = _mm_set_ps1(m1.v1.x);
		_y = _mm_set_ps1(m1.v1.y);
		_z = _mm_set_ps1(m1.v1.z);
		_w = _mm_set_ps1(m1.v1.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(m1.v2.x);
		_y = _mm_set_ps1(m1.v2.y);
		_z = _mm_set_ps1(m1.v2.z);
		_w = _mm_set_ps1(m1.v2.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(m1.v3.x);
		_y = _mm_set_ps1(m1.v3.y);
		_z = _mm_set_ps1(m1.v3.z);
		_w = _mm_set_ps1(m1.v3.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;


		// (MxM) * (Matrix) // ===============================================//

		// ROW 0 //

		_x = _mm_set_ps1(tmp.v0.x);	
		_y = _mm_set_ps1(tmp.v0.y);	
		_z = _mm_set_ps1(tmp.v0.z);	
		_w = _mm_set_ps1(tmp.v0.w);	

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;

		// ROW 1 //

		_x = _mm_set_ps1(tmp.v1.x);
		_y = _mm_set_ps1(tmp.v1.y);
		_z = _mm_set_ps1(tmp.v1.z);
		_w = _mm_set_ps1(tmp.v1.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(tmp.v2.x);
		_y = _mm_set_ps1(tmp.v2.y);
		_z = _mm_set_ps1(tmp.v2.z);
		_w = _mm_set_ps1(tmp.v2.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(tmp.v3.x);
		_y = _mm_set_ps1(tmp.v3.y);
		_z = _mm_set_ps1(tmp.v3.z);
		_w = _mm_set_ps1(tmp.v3.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;


		return tmp;
	}
};

struct MxM4
{
	const Matrix& m1;
	const Matrix& m2;
	const Matrix& m3;
	const Matrix& m4;

	MxM4(const MxM3& t1, const Matrix& t2)
		:m1(t1.m1), m2(t1.m2), m3(t1.m3), m4(t2) {};

	operator Matrix()
	{
		Matrix tmp;

		// (Matrix) * (Matrix) // ===============================================//

		// ROW 0 //

		__m128 _x = _mm_set_ps1(m1.v0.x);
		__m128 _y = _mm_set_ps1(m1.v0.y);
		__m128 _z = _mm_set_ps1(m1.v0.z);
		__m128 _w = _mm_set_ps1(m1.v0.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;


		// ROW 1 //

		_x = _mm_set_ps1(m1.v1.x);
		_y = _mm_set_ps1(m1.v1.y);
		_z = _mm_set_ps1(m1.v1.z);
		_w = _mm_set_ps1(m1.v1.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(m1.v2.x);
		_y = _mm_set_ps1(m1.v2.y);
		_z = _mm_set_ps1(m1.v2.z);
		_w = _mm_set_ps1(m1.v2.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(m1.v3.x);
		_y = _mm_set_ps1(m1.v3.y);
		_z = _mm_set_ps1(m1.v3.z);
		_w = _mm_set_ps1(m1.v3.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;


		// (MxM) * (Matrix) // ===============================================//

		// ROW 0 //

		_x = _mm_set_ps1(tmp.v0.x);
		_y = _mm_set_ps1(tmp.v0.y);
		_z = _mm_set_ps1(tmp.v0.z);
		_w = _mm_set_ps1(tmp.v0.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;

		// ROW 1 //

		_x = _mm_set_ps1(tmp.v1.x);
		_y = _mm_set_ps1(tmp.v1.y);
		_z = _mm_set_ps1(tmp.v1.z);
		_w = _mm_set_ps1(tmp.v1.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(tmp.v2.x);
		_y = _mm_set_ps1(tmp.v2.y);
		_z = _mm_set_ps1(tmp.v2.z);
		_w = _mm_set_ps1(tmp.v2.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(tmp.v3.x);
		_y = _mm_set_ps1(tmp.v3.y);
		_z = _mm_set_ps1(tmp.v3.z);
		_w = _mm_set_ps1(tmp.v3.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;


		// (MxM3) * (Matrix) // ===============================================//

		// ROW 0 //

		_x = _mm_set_ps1(tmp.v0.x);
		_y = _mm_set_ps1(tmp.v0.y);
		_z = _mm_set_ps1(tmp.v0.z);
		_w = _mm_set_ps1(tmp.v0.w);

		_x = _mm_mul_ps(m4.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m4.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m4.v2._m, _z);
		_y = _mm_mul_ps(m4.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;

		// ROW 1 //

		_x = _mm_set_ps1(tmp.v1.x);
		_y = _mm_set_ps1(tmp.v1.y);
		_z = _mm_set_ps1(tmp.v1.z);
		_w = _mm_set_ps1(tmp.v1.w);

		_x = _mm_mul_ps(m4.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m4.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m4.v2._m, _z);
		_y = _mm_mul_ps(m4.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(tmp.v2.x);
		_y = _mm_set_ps1(tmp.v2.y);
		_z = _mm_set_ps1(tmp.v2.z);
		_w = _mm_set_ps1(tmp.v2.w);

		_x = _mm_mul_ps(m4.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m4.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m4.v2._m, _z);
		_y = _mm_mul_ps(m4.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(tmp.v3.x);
		_y = _mm_set_ps1(tmp.v3.y);
		_z = _mm_set_ps1(tmp.v3.z);
		_w = _mm_set_ps1(tmp.v3.w);

		_x = _mm_mul_ps(m4.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m4.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m4.v2._m, _z);
		_y = _mm_mul_ps(m4.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;

		return tmp;
	}
};

struct MxM5
{
	const Matrix& m1;
	const Matrix& m2;
	const Matrix& m3;
	const Matrix& m4;
	const Matrix& m5;

	MxM5(const MxM4& t1, const Matrix& t2)
		:m1(t1.m1), m2(t1.m2), m3(t1.m3), m4(t1.m4), m5(t2) {};

	operator Matrix()
	{
		Matrix tmp;

		// (Matrix) * (Matrix) // ===============================================//

		// ROW 0 //

		__m128 _x = _mm_set_ps1(m1.v0.x);
		__m128 _y = _mm_set_ps1(m1.v0.y);
		__m128 _z = _mm_set_ps1(m1.v0.z);
		__m128 _w = _mm_set_ps1(m1.v0.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;


		// ROW 1 //

		_x = _mm_set_ps1(m1.v1.x);
		_y = _mm_set_ps1(m1.v1.y);
		_z = _mm_set_ps1(m1.v1.z);
		_w = _mm_set_ps1(m1.v1.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(m1.v2.x);
		_y = _mm_set_ps1(m1.v2.y);
		_z = _mm_set_ps1(m1.v2.z);
		_w = _mm_set_ps1(m1.v2.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(m1.v3.x);
		_y = _mm_set_ps1(m1.v3.y);
		_z = _mm_set_ps1(m1.v3.z);
		_w = _mm_set_ps1(m1.v3.w);

		_x = _mm_mul_ps(m2.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m2.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m2.v2._m, _z);
		_y = _mm_mul_ps(m2.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;


		// (MxM) * (Matrix) // ===============================================//

		// ROW 0 //

		_x = _mm_set_ps1(tmp.v0.x);
		_y = _mm_set_ps1(tmp.v0.y);
		_z = _mm_set_ps1(tmp.v0.z);
		_w = _mm_set_ps1(tmp.v0.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;

		// ROW 1 //

		_x = _mm_set_ps1(tmp.v1.x);
		_y = _mm_set_ps1(tmp.v1.y);
		_z = _mm_set_ps1(tmp.v1.z);
		_w = _mm_set_ps1(tmp.v1.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(tmp.v2.x);
		_y = _mm_set_ps1(tmp.v2.y);
		_z = _mm_set_ps1(tmp.v2.z);
		_w = _mm_set_ps1(tmp.v2.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(tmp.v3.x);
		_y = _mm_set_ps1(tmp.v3.y);
		_z = _mm_set_ps1(tmp.v3.z);
		_w = _mm_set_ps1(tmp.v3.w);

		_x = _mm_mul_ps(m3.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m3.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m3.v2._m, _z);
		_y = _mm_mul_ps(m3.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;


		// (MxM3) * (Matrix) // ===============================================//

		// ROW 0 //

		_x = _mm_set_ps1(tmp.v0.x);
		_y = _mm_set_ps1(tmp.v0.y);
		_z = _mm_set_ps1(tmp.v0.z);
		_w = _mm_set_ps1(tmp.v0.w);

		_x = _mm_mul_ps(m4.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m4.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m4.v2._m, _z);
		_y = _mm_mul_ps(m4.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;

		// ROW 1 //

		_x = _mm_set_ps1(tmp.v1.x);
		_y = _mm_set_ps1(tmp.v1.y);
		_z = _mm_set_ps1(tmp.v1.z);
		_w = _mm_set_ps1(tmp.v1.w);

		_x = _mm_mul_ps(m4.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m4.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m4.v2._m, _z);
		_y = _mm_mul_ps(m4.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(tmp.v2.x);
		_y = _mm_set_ps1(tmp.v2.y);
		_z = _mm_set_ps1(tmp.v2.z);
		_w = _mm_set_ps1(tmp.v2.w);

		_x = _mm_mul_ps(m4.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m4.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m4.v2._m, _z);
		_y = _mm_mul_ps(m4.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(tmp.v3.x);
		_y = _mm_set_ps1(tmp.v3.y);
		_z = _mm_set_ps1(tmp.v3.z);
		_w = _mm_set_ps1(tmp.v3.w);

		_x = _mm_mul_ps(m4.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m4.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m4.v2._m, _z);
		_y = _mm_mul_ps(m4.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;


		// (MxM4) * (Matrix) // ===============================================//


		// ROW 0 //

		_x = _mm_set_ps1(tmp.v0.x);
		_y = _mm_set_ps1(tmp.v0.y);
		_z = _mm_set_ps1(tmp.v0.z);
		_w = _mm_set_ps1(tmp.v0.w);

		_x = _mm_mul_ps(m5.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
		_y = _mm_mul_ps(m5.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m5.v2._m, _z);
		_y = _mm_mul_ps(m5.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v0._m = _z;

		// ROW 1 //

		_x = _mm_set_ps1(tmp.v1.x);
		_y = _mm_set_ps1(tmp.v1.y);
		_z = _mm_set_ps1(tmp.v1.z);
		_w = _mm_set_ps1(tmp.v1.w);

		_x = _mm_mul_ps(m5.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
		_y = _mm_mul_ps(m5.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m5.v2._m, _z);
		_y = _mm_mul_ps(m5.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v1._m = _z;


		// ROW 2 //

		_x = _mm_set_ps1(tmp.v2.x);
		_y = _mm_set_ps1(tmp.v2.y);
		_z = _mm_set_ps1(tmp.v2.z);
		_w = _mm_set_ps1(tmp.v2.w);

		_x = _mm_mul_ps(m5.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
		_y = _mm_mul_ps(m5.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m5.v2._m, _z);
		_y = _mm_mul_ps(m5.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v2._m = _z;


		// ROW 3 //

		_x = _mm_set_ps1(tmp.v3.x);
		_y = _mm_set_ps1(tmp.v3.y);
		_z = _mm_set_ps1(tmp.v3.z);
		_w = _mm_set_ps1(tmp.v3.w);

		_x = _mm_mul_ps(m5.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
		_y = _mm_mul_ps(m5.v1._m, _y);
		_x = _mm_add_ps(_x, _y);

		_z = _mm_mul_ps(m5.v2._m, _z);
		_y = _mm_mul_ps(m5.v3._m, _w);
		_y = _mm_add_ps(_z, _y);

		_z = _mm_add_ps(_x, _y);

		tmp.v3._m = _z;

		return tmp;
	}
};


// INLINE FUNCTIONS

inline MxM operator * (const Matrix& m1, const Matrix& m2)
{
	{
		return MxM(m1, m2);
	}
}

inline MxM3 operator * (const MxM& m1, const Matrix& m2)
{
	{
		return MxM3(m1, m2);
	}
}

inline MxM4 operator * (const MxM3& m1, const Matrix& m2)
{
	{
		return MxM4(m1, m2);
	}
}

inline MxM5 operator * (const MxM4& m1, const Matrix& m2)
{
	{
		return MxM5(m1, m2);
	}
}

#endif  

// --- End of File ---
