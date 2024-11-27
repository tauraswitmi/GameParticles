//---------------------------------------------------------------
// Copyright 2024, Ed Keenan, all rights reserved.
//---------------------------------------------------------------

#include "OpenGLDevice.h"
#include "ParticleEmitter.h"
#include "Settings.h"

PerformanceTimer globalTimer;

ParticleEmitter::ParticleEmitter()
:	start_position( 0.0f, 2.0f, 2.0f ),
	start_velocity( -4.0f, 4.0f, 0.0f), 
	spawn_frequency(0.00001f),		
	last_spawn(globalTimer.GetGlobalTime()),		
	last_loop(globalTimer.GetGlobalTime()),
	max_life( MAX_LIFE ),
	max_particles( NUM_PARTICLES ),
	last_active_particle(-1),
	bufferCount(0),
	headParticle(nullptr),
	vel_variance(15.0f, 0.70f, -1.0f),
	pos_variance(1.50f, 0.50f, 10.0f),
	scale_variance(3.0f)
{
	// nothing to do
	bufferCount = 0;
}

ParticleEmitter::~ParticleEmitter()
{
	Particle *pTmp = this->headParticle;

	//Particle* pDeleteMe;
	while (pTmp != nullptr)
	{
		Particle *pDeleteMe = pTmp;
		//pDeleteMe = pTmp;
		pTmp = pTmp->next;
		delete pDeleteMe;
	}
}

//999
void ParticleEmitter::SpawnParticle()
{
	// eate another particle if there are ones free
	if( last_active_particle < max_particles-1 )
	{
		// create a new particle
		pNewParticle = new Particle();

		// 999 
		// initialize the particle
		pNewParticle->life     = 0.0f;
		pNewParticle->position = start_position;
		pNewParticle->velocity = start_velocity;
		pNewParticle->scale    = Vect4D(-1.0, -1.0, -1.0, 1.0);

		// apply the variance
		this->Execute(pNewParticle->position, pNewParticle->velocity, pNewParticle->scale);

		// increment count
		last_active_particle++;

		// add to list
		this->addParticleToList( pNewParticle );

	}
}




// OLD UPDATE
void ParticleEmitter::update()
{

	// get current time
	float current_time = globalTimer.GetGlobalTime();

	// spawn particles
	float time_elapsed = current_time - this->last_spawn;
	
	// update
	while( spawn_frequency < time_elapsed )
	{
		// spawn a particle
		this->SpawnParticle();
		// adjust time
		time_elapsed -= spawn_frequency;
		// last time
		this->last_spawn = current_time;
	}
	
	// total elapsed
	time_elapsed = current_time - last_loop;

	Particle *p = this->headParticle;
	// walk the particles

	while( p != nullptr )
	{
		// call every particle and update its position 
		p->Update(time_elapsed);

		// if life is greater that the max_life 
		// and there is some left on the list
		// remove node
		if((last_active_particle > 0) && (p->life > max_life))
		{
			// particle to remove
			Particle *s = p;

			// need to squirrel it away.
			p=p->next;

			// remove last node
			this->removeParticleFromList(s);

			// update the number of particles
			last_active_particle--;
		}
		else
		{
			// increment to next point
			p = p->next;
		}
	}

	last_loop = current_time;
}
	   











//999 I wonder if changing the little bit from the unopt code would make a diff?
void ParticleEmitter::addParticleToList(Particle *p )
{
	assert(p);
	if( this->headParticle == nullptr )
	{ // first on list (empty)
		this->headParticle = p;
		p->next = 0;
		p->prev= 0;
	}
	else 
	{ // add to front of list
		headParticle->prev = p;
		p->next = headParticle;
		p->prev = 0;
		headParticle = p;
	}

}


// NO CHANGES NEEDED //
void ParticleEmitter::removeParticleFromList( Particle *p )
{
	// make sure we are not screwed with a null pointer
	assert(p);

	if( p->prev == nullptr && p->next == nullptr  )
	{ // only one on the list
		this->headParticle = nullptr;
	}
	else if( p->prev == nullptr && p->next != nullptr  )
	{ // first on the list
		p->next->prev = nullptr;
		this->headParticle = p->next;
	}
	else if( p->prev!= nullptr && p->next == nullptr )
	{ // last on the list 
		p->prev->next = nullptr;
	}
	else//( p->next != nullptr  && p->prev !=nullptr )
	{ // middle of the list
		p->prev->next = p->next;
		p->next->prev = p->prev;
	}
	
	// bye bye
	delete p;
}


void ParticleEmitter::draw()
{
	// initialize the camera matrix
	//Matrix cameraMatrix;
	cameraMatrix.setIdentMatrix();

	// setup the translation matrix
	//Matrix transMatrix;
	Vect4D trans(0.0f, 5.0f, 40.0f);
	transMatrix.setTransMatrix(trans); //888

	//Matrix tmp;
	// multiply them together
	tmp = cameraMatrix * transMatrix;

	// get the inverse matrix
	Matrix inverseCameraMatrix;
	tmp.Inverse(inverseCameraMatrix);


	// iterate throught the list of particles
	Particle* temp = this->headParticle;

	while (temp != 0)
	{

		// get the position from this matrix
		Vect4D camPosVect;
		inverseCameraMatrix.get(Matrix::MatrixRow::MATRIX_ROW_3, camPosVect); //88

		// camera position
		//Matrix transCamera;
		transCamera.setTransMatrix(camPosVect); //88

		// particle position
		//Matrix transParticle;
		transParticle.setTransMatrix(temp->position); //88

		// rotation matrix
		//Matrix rotParticle;
		rotParticle.setRotZMatrix(temp->rotation);

		// scale Matrix
		scaleMatrix.setScaleMatrix(temp->scale); //55

		// total transformation of particle
		tmp = scaleMatrix * transCamera * transParticle * rotParticle * scaleMatrix; ///99999 PROXIES!!!!!

		// ------------------------------------------------
		//  Set the Transform Matrix and Draws Triangle
		//  Note: 
		//       this method is using doubles... 
		//       there is a float version (hint)
		// ------------------------------------------------
		OpenGLDevice::SetTransformMatrixFloat((const float*)&tmp);

		// squirrel away matrix for next update
		tmp.get(Matrix::MatrixRow::MATRIX_ROW_0, temp->curr_Row0);		  //88
		tmp.get(Matrix::MatrixRow::MATRIX_ROW_1, temp->curr_Row1);		  //88
		tmp.get(Matrix::MatrixRow::MATRIX_ROW_2, temp->curr_Row2);		  //88
		tmp.get(Matrix::MatrixRow::MATRIX_ROW_3, temp->curr_Row3);		  //88

		// difference vector
		temp->diff_Row0 = temp->curr_Row0 - temp->prev_Row0;
		temp->diff_Row1 = temp->curr_Row1 - temp->prev_Row1;
		temp->diff_Row2 = temp->curr_Row2 - temp->prev_Row2;
		temp->diff_Row3 = temp->curr_Row3 - temp->prev_Row3;


		// Get the next particle in the list!
		temp = temp->next;
	}

}

void ParticleEmitter::Execute(Vect4D& pos, Vect4D& vel, Vect4D& sc)
{
	// Ses it's ugly - I didn't write this so don't bitch at me
	// Sometimes code like this is inside real commerical code ( so know you now how it feels )
	
	// x - variance
	float var = static_cast<float>(rand() % 1000) * 0.001f;
	float sign = static_cast<float>(rand() % 2);
	float *t_pos = reinterpret_cast<float*>(&pos);
	float *t_var = &pos_variance[Vect::X];
	if(sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	// y - variance
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if(sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;
	
	// z - variance
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if(sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;
	
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	
	// x  - add velocity
	t_pos = &vel[Vect::X];
	t_var = &vel_variance[Vect::X];
	if(sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;
	
	// y - add velocity
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if(sign == 0)
	{
		var *= -2.0;
	}
	*t_pos += *t_var * var;
	
	// z - add velocity
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if(sign == 0)
	{
		var *= -2.0;
	}
	*t_pos += *t_var * var;
	
	// correct the sign
	var = 2.0f * static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	
	if(sign == 0)
	{
		var *= -2.0;
	}
	sc = sc * var;	// *= 8888
}

// --- End of File ---
