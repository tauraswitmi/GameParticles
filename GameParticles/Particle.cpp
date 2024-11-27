//---------------------------------------------------------------
// Copyright 2024, Ed Keenan, all rights reserved.
//---------------------------------------------------------------

#include "Particle.h"
#include "Matrix.h"

Particle::Particle()
{
	// construtor
	this->life = 0.0f;
	this->position.set( 0.0f, 0.0f,  -10.0f );
	this->velocity.set( -3.0f, 0.0f,  0.0f );
	this->scale.set( 1.0f, 1.0f, 1.0f );
	this->rotation = 0.0f;
	this->rotation_velocity = -0.25f;
	this->next = nullptr;
	this->prev = nullptr;
}


Particle::~Particle()
{
	// nothing to do
}

void Particle::CopyDataOnly(const Particle &p )
{
	// copy the data only
	// this way of copying data is more efficient that element by element
	this->position = p.position;
	this->velocity = p.velocity;
	this->scale    = p.scale;
	this->rotation = p.rotation;
	this->rotation_velocity = p.rotation_velocity;
	this->life     = p.life;
}

void Particle::Update(const float& time_elapsed)
{
	// Rotate the matrices


	this->prev_Row0 = this->curr_Row0;
	this->prev_Row1 = this->curr_Row1;
	this->prev_Row2 = this->curr_Row2;
	this->prev_Row3 = this->curr_Row3;


	tmp.set(Matrix::MatrixRow::MATRIX_ROW_0, this->diff_Row0); //88
	tmp.set(Matrix::MatrixRow::MATRIX_ROW_1, this->diff_Row1); //88
	tmp.set(Matrix::MatrixRow::MATRIX_ROW_2, this->diff_Row2); //88
	tmp.set(Matrix::MatrixRow::MATRIX_ROW_3, this->diff_Row3); //88

	float MatrixScale = -3.0f*tmp.Determinant();

	// serious math below - magic secret sauce
	this->life += time_elapsed;
	this->position += (this->velocity * time_elapsed);		// Implemented +=  88888


	Vect4D z_axis(0.0f, 0.0f, 3.0f);			
	Vect4D v(2,2,0);
	position.Cross( z_axis, v);	// Might Change so that it's not a function call 888


	v.norm(v);
	position += (v * 0.05f * life); // += IMPLEMENT HERE 88888888

	if( MatrixScale > 1.0 )
	{
		MatrixScale = 1.0f/MatrixScale;
	};



	// Changes the rotation of the particle
	rotation += MatrixScale + rotation_velocity * time_elapsed * 2;
}

// --- End of File ---
