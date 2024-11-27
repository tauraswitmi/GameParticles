//---------------------------------------------------------------
// Copyright 2024, Ed Keenan, all rights reserved.
//---------------------------------------------------------------

#ifndef Vect4D_H
#define Vect4D_H

// includes
#include "Enum.h"

// Foward Declarations
class Matrix;

// class
class Vect4D
{
public:

	friend class Matrix;

	Vect4D();
	Vect4D(const float tx, const float ty, const float tz, const float tw = 1.0f);
	Vect4D(const Vect4D& r);							
	Vect4D(const __m128& r);
	Vect4D& operator = (const Vect4D& r) = default;				
	~Vect4D();

	void Cross(Vect4D &vin, Vect4D &vout) const;
	float &operator[](Vect e);
	void norm(Vect4D &out);
	void set(float tx, float ty, float tz, float tw = 1.0f);


	Vect4D operator + (const Vect4D& t) const;
	Vect4D operator - (const Vect4D& t) const;
	//Vect4D operator * (const float scale) const;
	

	Vect4D operator *= (const float scale) const; 
	Vect4D& operator += (const Vect4D& r); 


public:

	union
	{
		__m128	_m;
	
		struct
		{
			float x;
			float y;
			float z;
			float w;
		};
	};
	
	friend struct VxF;
	friend struct VxFxF;
};



struct VxF
{
	const Vect4D& v1;
	const float f1;

	VxF(const Vect4D& t1, const float tf1)
		:v1(t1), f1(tf1)
	{
	}

	operator Vect4D()
	{
		__m128 tmp = _mm_set1_ps(f1);
		return Vect4D(_mm_mul_ps(v1._m, tmp));
	}
};

struct  VxFxF
{
	const Vect4D& v1;
	const float f1;
	const float f2;

	VxFxF(const VxF& t1, const float tf1)
		:v1(t1.v1), f1(t1.f1), f2(tf1)
	{
	}

	operator Vect4D()
	{
		__m128 tmp = _mm_set1_ps(f2 * f1);
		return Vect4D(_mm_mul_ps(v1._m, tmp));
	}
};

inline VxF operator * (const Vect4D& a1, const float& af1)
{
	return VxF(a1, af1);
}

inline VxFxF operator * (const VxF& a1, const float& af1)
{
	return VxFxF(a1, af1);
}


#endif


// --- End of File ---
