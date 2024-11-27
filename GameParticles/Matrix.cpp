//---------------------------------------------------------------
// Copyright 2024, Ed Keenan, all rights reserved.
//---------------------------------------------------------------

#include <Math.h>
#include "Vect4d.h"
#include "Matrix.h"

Matrix::Matrix()
{	
	this->v0._m = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
	this->v1._m = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
	this->v2._m = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
	this->v3._m = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
}

Matrix::Matrix(const Matrix& t)
{
	this->v0._m = t.v0._m;
	this->v1._m = t.v1._m;
	this->v2._m = t.v2._m;
	this->v3._m = t.v3._m;
}

Matrix::~Matrix()
{
	
}

void Matrix::setDefault()
{
	this->v0._m = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
	this->v1._m = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
	this->v2._m = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
	this->v3._m = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
}

void Matrix::setIdentMatrix()
{ 
	this->v0._m = _mm_set_ps(0.0, 0.0, 0.0, 1.0);

	this->v1._m = _mm_set_ps(0.0, 0.0, 1.0, 0.0);

	this->v2._m = _mm_set_ps(0.0, 1.0, 0.0, 0.0);
	
	this->v3._m = _mm_set_ps(1.0, 0.0, 0.0, 0.0);
}

void Matrix::setTransMatrix(Vect4D &t)
{ 
	this->v0._m = _mm_set_ps(0.0, 0.0, 0.0, 1.0);

	this->v1._m = _mm_set_ps(0.0, 0.0, 1.0, 0.0);

	this->v2._m = _mm_set_ps(0.0, 1.0, 0.0, 0.0);
	
	this->v3._m = _mm_set_ps(1.0, t.z, t.y, t.x);
}

void Matrix::set( MatrixRow row, Vect4D &t )
{
	switch( row )
	{
	case MatrixRow::MATRIX_ROW_0:
		this->v0._m = t._m;
		break;

	case MatrixRow::MATRIX_ROW_1:
		this->v1._m = t._m;
		break;

	case MatrixRow::MATRIX_ROW_2:
		this->v2._m = t._m;
		break;

	case MatrixRow::MATRIX_ROW_3:
		this->v3._m = t._m;
		break;

	default:
		// should never get here, if we are here something bad has happened
		assert(0);
	}
}

float &Matrix::operator[](Index e)
{
	// get the individual elements
	switch(e)
	{
	case Index::m0:
		return m0;
		break;
	case Index::m1:
		return m1;
		break;
	case Index::m2:
		return m2;
		break;
	case Index::m3:
		return m3;
		break;
	case Index::m4:
		return m4;
		break;
	case Index::m5:
		return m5;
		break;
	case Index::m6:
		return m6;
		break;
	case Index::m7:
		return m7;
		break;
	case Index::m8:
		return m8;
		break;
	case Index::m9:
		return m9;
		break;	
	case Index::m10:
		return m10;
		break;
	case Index::m11:
		return m11;
		break;	
	case Index::m12:
		return m12;
		break;	
	case Index::m13:
		return m13;
		break;	
	case Index::m14:
		return m14;
		break;
	case Index::m15:
		return m15;
		break;
	default:
		assert(0);
		return m0;
		break;
	}
}

void Matrix::get( MatrixRow row, Vect4D &t ) const
{ // get a row of the matrix
	switch( row )
	{
	case MatrixRow::MATRIX_ROW_0:
		t._m = this->v0._m;
		break;

	case MatrixRow::MATRIX_ROW_1:
		t._m = this->v1._m;
		break;

	case MatrixRow::MATRIX_ROW_2:
		t._m = this->v2._m;
		break;

	case MatrixRow::MATRIX_ROW_3:
		t._m = this->v3._m;
		break;

	default:
		assert(0);
	}
}



Matrix& Matrix::operator*=( const Matrix& rhs) // Check if (w remains 1)
{ 
	//Matrix tmp;

	// ROW 0 //

	__m128 _x = _mm_set_ps1(this->v0.x);	// (x0, x0, x0, x0)
	__m128 _y = _mm_set_ps1(this->v0.y);	// (y0, y0, y0, y0)
	__m128 _z = _mm_set_ps1(this->v0.z);	// (z0, z0, z0, z0)
	__m128 _w = _mm_set_ps1(this->v0.w);	// (w0, w0, w0, w0)
	
	_x = _mm_mul_ps(rhs.v0._m, _x);			// (m0, m1, m2, m3) (x0, x0, x0, x0) == (x0m0, x0m1, x0m2, x0m3)
	_y = _mm_mul_ps(rhs.v1._m, _y);
	_x = _mm_add_ps(_x, _y);
	
	_z = _mm_mul_ps(rhs.v2._m, _z);
	_y = _mm_mul_ps(rhs.v3._m, _w);
	_y = _mm_add_ps(_z, _y);
	
	_z = _mm_add_ps(_x, _y);
	
	this->v0._m = _z;
	
	
	// ROW 1 //
	
	_x = _mm_set_ps1(this->v1.x);
	_y = _mm_set_ps1(this->v1.y);
	_z = _mm_set_ps1(this->v1.z);
	_w = _mm_set_ps1(this->v1.w);
	
	_x = _mm_mul_ps(rhs.v0._m, _x);			// (m0, m1, m2, m3) (x1, x1, x1, x1) == (x1m0, x1m1, x1m2, x1m3)
	_y = _mm_mul_ps(rhs.v1._m, _y);
	_x = _mm_add_ps(_x, _y);
	
	_z = _mm_mul_ps(rhs.v2._m, _z);
	_y = _mm_mul_ps(rhs.v3._m, _w);
	_y = _mm_add_ps(_z, _y);
	
	_z = _mm_add_ps(_x, _y);
	
	this->v1._m = _z;
	
	
	// ROW 2 //
	
	_x = _mm_set_ps1(this->v2.x);
	_y = _mm_set_ps1(this->v2.y);
	_z = _mm_set_ps1(this->v2.z);
	_w = _mm_set_ps1(this->v2.w);
	
	_x = _mm_mul_ps(rhs.v0._m, _x);			// (m0, m1, m2, m3) (x2, x2, x2, x2) == (x2m0, x2m1, x2m2, x2m3)
	_y = _mm_mul_ps(rhs.v1._m, _y);
	_x = _mm_add_ps(_x, _y);
	
	_z = _mm_mul_ps(rhs.v2._m, _z);
	_y = _mm_mul_ps(rhs.v3._m, _w);
	_y = _mm_add_ps(_z, _y);
	
	_z = _mm_add_ps(_x, _y);
	
	this->v2._m = _z;
	
	
	// ROW 3 //
	
	_x = _mm_set_ps1(this->v3.x);
	_y = _mm_set_ps1(this->v3.y);
	_z = _mm_set_ps1(this->v3.z);
	_w = _mm_set_ps1(this->v3.w);
	
	_x = _mm_mul_ps(rhs.v0._m, _x);			// (m0, m1, m2, m3) (x3, x3, x3, x3) == (x3m0, x3m1, x3m2, x3m3)
	_y = _mm_mul_ps(rhs.v1._m, _y);
	_x = _mm_add_ps(_x, _y);
	
	_z = _mm_mul_ps(rhs.v2._m, _z);
	_y = _mm_mul_ps(rhs.v3._m, _w);
	_y = _mm_add_ps(_z, _y);
	
	_z = _mm_add_ps(_x, _y);
	
	this->v3._m = _z;
	

	return *this;
}

Matrix& Matrix::operator/=(const float rhs)
{ 
	// divide each element by a value
	// using inverse multiply trick, faster that individual divides
	//float inv_rhs = 1.0f/rhs;

	__m128 div = _mm_set_ps1(1.0f / rhs);
	
	this->v0._m = _mm_mul_ps(this->v0._m, div);
	this->v1._m = _mm_mul_ps(this->v1._m, div);
	this->v2._m = _mm_mul_ps(this->v2._m, div);
	this->v3._m = _mm_mul_ps(this->v3._m, div);


	return *this;
}

float Matrix::Determinant() const
{
	// Can improve with RVO (will look messy though)

	float ta = (m10 * m15) - (m11 * m14);
	// tb = (kq - mo)
	float tb = (m9 * m15) - (m11 * m13);
	// tc = (kp - lo)
	float tc = (m9 * m14) - (m10 * m13);
	// td = (jq - mn)
	float td = (m8 * m15) - (m11 * m12);
	// te = (jo - kn)
	float te = (m8 * m13) - (m9 *  m12);
	// tf = (jp - ln)
	float tf = (m8 * m14) - (m10 * m12);
	
	return ((m0 * ((m5 * ta) - (m6 * tb) + (m7 * tc)))
			- (m1 * ((m4 * ta) - (m6 * td) + (m7 * tf)))
			+ (m2 * ((m4 * tb) - (m5 * td) + (m7 * te)))
			- (m3 * ((m4 * tc) - (m5 * tf) + (m6 * te))));
	

	// v0 (
}

Matrix Matrix::GetAdjugate() const
{
	Matrix tmp;
	
	// load		ABC		(3)		ABC--
	float t1 = (m10*m15) - (m11*m14);
	float t2 = (m9*m15) - (m11*m13);
	float t3 = (m9*m14) - (m10*m13);
	
	// a = f(ta) - g(tb) + h(tc)
	tmp.m0 = (m5*t1) - (m6*t2) + (m7*t3);				// m0	0
	// b = -( b(ta) - c(tb) + d(tc))
	tmp.m1 = -((m1*t1) - (m2*t2) + (m3*t3));			// m1	0
	
	// load		JK		(5)		ABCJK
	float t4 = (m8*m15) - (m11*m12);
	float t5 = (m8*m14) - (m10*m12);
	// e = - ( e(ta) - g(tj) + h(tk))
	tmp.m4 = -( (m4*t1) - (m6*t4) + (m7*t5));			// m4	1
	// f = a(ta) - c(tj) + d(tk)
	tmp.m5 = (m0*t1) - (m2*t4) + (m3*t5);				// m5	1
	
	// unload	AJ		(3)		-BC-K
	// load		P		(4)		PBC-K
	t1 = (m8*m13) - (m9*m12);
	// n = -( e(tc) - f(tk) + g(tp))
	tmp.m12 = -((m4*t3) - (m5*t5) + (m6*t1));			// m12	3
	// o = a(tc) - b(tk) + c(tp)
	tmp.m13 = (m0*t3) - (m1*t5) + (m2*t1);				// m13	3
	
	// unload	KC		(2)		PB---
	// load		J		(3)		PBJ--
	t3 = (m8*m15) - (m11*m12);
	
	// j = e(tb) - f(tj) + h(tp)
	tmp.m8 = (m4*t2) - (m5*t3) + (m7*t1);				// m8	2
	// k = - ( a(tb) - b(tj) + d(tp))
	tmp.m9 = -((m0*t2) - (m1*t3) + (m3*t1));			// m9	2
	
	// unload	BPJ		(0)		-----
	// load		DLM		(3)		DLM--
	t1 = (m6*m15) - (m7*m14);
	t2 = (m4*m15) - (m7*m12);
	t3 = (m4*m14) - (m6*m12);
	
	// g = - ( a(td) - c(tl) + d(tm))
	tmp.m6 = -((m0*t1) - (m2*t2) + (m3*t3));			// m6	1
	
	// load		FR		(5)		DLMFR
	t4 = (m5*m14) - (m6*m13);
	t5 = (m4*m13) - (m5*m12);
	
	// p = - ( a(tf) - b(tm) + c(tr))
	tmp.m14 = -( (m0*t4) - (m1*t3) + (m3*t5));			// m14	3
	
	// unload	M		(4)		DL-FR
	// load		E		(5)		DLEFR
	t3 = (m5*m15) - (m7*m13);
	
	// l = a(te) - b(tl) + d(tr)
	tmp.m10 = (m0*t3) - (m1*t2) + (m3*t5);				// m10	2
	
	// unload	LR		(3)		D-EF-
	// c = b(td) - c(te) + d(tf)
	tmp.m2 = (m1*t1) - (m2*t3) + (m3*t4);				// m2	0
	
	// unload	DEF		(0)		-----
	// load		GHI		(3)		GHI--
	t1 = (m6*m11) - (m7*m10);
	t2 = (m5*m11) - (m7*m9);
	t3 = (m5*m10) - (m6*m9);
	
	// d = -( b(tg) - c(th) + d(ti) )
	tmp.m3 = -( (m1*t1) - (m2*t2) + (m3*t3));			// m3	0
	
	// load		NO		(5)		GHINO
	t4 = (m4*m11) - (m7*m8);
	t5 = (m4*m10) - (m6*m8);
	
	// h = a(tg) - c(tn) + d(to)
	tmp.m7 = (m0*t1) - (m2*t4) + (m3*t5);				// m7	1
	
	// unload	G		(4)		-HINO
	// load		Q		(5)		QHINO
	t1 = (m4*m9) - (m5*m8);
	
	// m = -( a(th) - b(tn) + d(tq))
	tmp.m11 = -((m0*t2) - (m1*t4) + (m3*t1));			// m11	2
	
	// unload	HN		(3)		Q-I-O
	// q = a(ti) - b(to) + c(tq)
	tmp.m15 = (m0*t3) - (m1*t5) + (m2*t1);				// m15	3
	
	return tmp;
}

void Matrix::Inverse( Matrix &out ) const 
{
	Matrix tmp;
	float det = Determinant();
	if(fabs(det) < 0.0001)
	{
		// do nothing, Matrix is not invertable
		//tmp.setDefault();
	}
	else
	{
		tmp = GetAdjugate();
		tmp /= det;
	}

	out = tmp;
}

void Matrix::setScaleMatrix(Vect4D &scale)
{
	//	{	sx		0		0		0	}
	//	{	0		sy		0		0	}
	//	{	0		0		sz		0	}
	//	{	0		0		0		1	}
	
	//Matrix tmp;

	this->v0._m = _mm_set_ps(0, 0, 0, scale.x); 
	
	this->v1._m = _mm_set_ps(0, 0, scale.y, 0);
	
	this->v2._m = _mm_set_ps(0, scale.z, 0, 0);
	
	this->v3._m = _mm_set_ps(1.0, 0, 0, 0);
}

void Matrix::setRotZMatrix(float az)
{
	//	{	cos		-sin	0		0	}
	//	{	sin		cos		0		0	}
	//	{	0		0		1		0	}
	//	{	0		0		0		1	}
	
	this->v0._m = _mm_set_ps(0, 0, -sin(az), cos(az));
	
	this->v1._m = _mm_set_ps(0, 0, cos(az), sin(az));
	
	this->v2._m = _mm_set_ps(0, 1.0, 0, 0);
	
	this->v3._m = _mm_set_ps(1.0, 0, 0, 0);
}



// --- End of File ---
