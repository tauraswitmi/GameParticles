//---------------------------------------------------------------
// Copyright 2024, Ed Keenan, all rights reserved.
//---------------------------------------------------------------

#ifndef PARTICLEEMITTER_H
#define PARTICLEEMITTER_H

#include "Matrix.h"
#include "Vect4D.h"
#include "Particle.h"

//#include <list>

class ParticleEmitter
{
public:
	ParticleEmitter();
	ParticleEmitter(const ParticleEmitter& r) = default;
	ParticleEmitter& operator= (const ParticleEmitter& r) = default;
	~ParticleEmitter();
	
	void SpawnParticle();
	void update();
	void draw();

	void addParticleToList(Particle *p );
	void removeParticleFromList( Particle *p );

	void Execute(Vect4D& pos, Vect4D& vel, Vect4D& sc);

private:

	Particle* pNewParticle;
	Particle* headParticle;

	Vect4D	start_position;
	Vect4D	start_velocity;

	float	spawn_frequency;
	float	last_spawn;
	float	last_loop;	
	float	max_life;
	int		max_particles;
	int		last_active_particle;
	int bufferCount;
	float	scale_variance;
	
	Vect4D	vel_variance;
	Vect4D	pos_variance;
	
	Matrix cameraMatrix;
	Matrix transMatrix;
	Matrix tmp;
	Matrix transCamera;
	Matrix transParticle;
	Matrix rotParticle;
	Matrix pivotParticle;
	Matrix scaleMatrix;

};

#endif 

// --- End of File ---
