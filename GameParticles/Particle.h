//---------------------------------------------------------------
// Copyright 2024, Ed Keenan, all rights reserved.
//---------------------------------------------------------------

#ifndef PARTICLE_H
#define PARTICLE_H

// include
#include "Vect4D.h"
#include "Matrix.h"
#include <new>

class Matrix;

class Particle
{
public:
	friend class ParticleEmitter;
	
	Particle();	
	Particle(const Particle& r) = default;					//888
	Particle& operator= (const Particle& r) = default;		//888
	~Particle();


	void Update(const float& time_elapsed);
	void CopyDataOnly( const Particle &p );

	#undef new
	void* operator new(size_t i)
	{
		return _mm_malloc(i, 16);
	}

	void operator delete(void* p)
	{
		return _mm_free(p);
	}

private:
	Particle *next;			// 8 8
	Particle *prev;			
	
	//float pad;
	//float pad2;
	//float pad3;
	
	// Padding makes it a factor of 16 (might be wrong because it could be doing allignment instead

	Matrix tmp;

	Vect4D	prev_Row0;		// 64 16	//	72	24
	Vect4D	prev_Row1;		
	Vect4D  prev_Row2;		
	Vect4D  prev_Row3;		

	Vect4D  diff_Row0;		// 64 16	//	136	40
	Vect4D  diff_Row1;		
	Vect4D  diff_Row2;		
	Vect4D  diff_Row3;		
	
	Vect4D	curr_Row0;		// 64 16	//	200 56
	Vect4D	curr_Row1;		
	Vect4D  curr_Row2;		
	Vect4D  curr_Row3;	

	Vect4D	position;		// 60 24	//	260	80
	Vect4D	velocity;
	Vect4D	scale;

	float	life;
	float	rotation;
	float	rotation_velocity;

	

	// Add float padd here
};

#endif 

// --- End of File ---
