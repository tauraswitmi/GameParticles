//---------------------------------------------------------------
// Copyright 2024, Ed Keenan, all rights reserved.
//---------------------------------------------------------------

#include <math.h>
#include <assert.h>
#include "Vect4D.h"

Vect4D::Vect4D()
	:x(0.0f), y(0.0f), z(0.0f), w(1.0f)
{
}

Vect4D::Vect4D(const float tx, const float ty, const float tz, const float tw)
	:x(tx), y(ty), z(tz), w(tw)
{
}

Vect4D::Vect4D(const Vect4D& r)
	:x(r.x), y(r.y), z(r.z), w(r.w)
{
}

Vect4D::Vect4D(const __m128& r)
	: _m(r)
{
}


Vect4D::~Vect4D()
{
	// nothing to delete
}

void Vect4D::norm(Vect4D &v)	
{
	float mag = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);


	if(0.0f < mag)
	{
		mag = 1 / mag;

		v.x = this->x * mag;
		v.y = this->y * mag;
		v.z = this->z * mag;
		v.w = 1.0;

	}
}

Vect4D Vect4D:: operator + (const Vect4D &t) const
{
	return Vect4D(_mm_add_ps(this->_m, t._m));
}

float &Vect4D::operator[](Vect e)
{
	switch(e)
	{
	case Vect::X:
		return x;
		break;
	case Vect::Y:
		return y;
		break;
	case Vect::Z:
		return z;
		break;
	case Vect::W:
		return w;
		break;
	default:
		assert(0);
		return x;
	}
}



Vect4D Vect4D::operator - (const Vect4D &t) const
{
	return Vect4D(_mm_sub_ps(this->_m, t._m));
}


// Called once
void Vect4D::Cross(Vect4D &vin, Vect4D &vout) const
{
	//vout._m = _mm_set_ps((y * vin.z - z * vin.y), (z * vin.x - x * vin.z), (x * vin.y - y * vin.x), 1.0f);

	vout.x = (y * vin.z - z * vin.y);
	vout.y = (z * vin.x - x * vin.z);
	vout.z = (x * vin.y - y * vin.x);
	vout.w = 1.0f;
}

void Vect4D::set(float tx, float ty, float tz, float tw)
{
	this->_m = _mm_set_ps(tw, tz, ty, tx);
}


//============================================================//

Vect4D Vect4D::operator*= (const float r) const
{
	return Vect4D(_mm_mul_ps(this->_m, _mm_set1_ps(r)));;
}


Vect4D& Vect4D::operator+=(const Vect4D& r)
{
	//this->x += r.x;
	//this->y += r.y;
	//this->z += r.z;

	//float last = this->w;
	//
	//this->_m = _mm_add_ps(this->_m, r._m);
	//this->w = 1.0f;

	this->_m = _mm_add_ps(this->_m, r._m);

	return *this;
}




// --- End of File ---
