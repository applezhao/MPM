#pragma once
#ifndef FF_PARTICLE_HH
#define FF_PARTICLE_HH
#include <vector>
#include <cinder/Vector.h>
#include "FF_material.h"
using namespace std;


class particle_body
{
public:
	material* mat;
	MaterialModel* material_model;
	int material_id;
	int component_id;//comid
	cinder::Vec3f gravity;//Gravp
	int begin;//the beginning of particle of body  par_begin
	int end;//the end of particle of body par_end
	particle_body()
	{
		gravity.set(0,0,0);
	}

	void set_material(const int& material_id, material* mat, MaterialModel* material_model, const int& component_id)
	{
		this->material_id=material_id;
		this->mat=mat;
		this->material_model=material_model;
		this->component_id=component_id;
	}
};

class particle_data
{
public:
	int num_particle;
	vector<particle*> particle_list;
	int num_component;
	int num_body;
	vector<particle_body*> body_list;

	//method
	bool MUSL;
	bool USL;
	bool USF;
	bool GIMP;
	bool contact;

	bool gravity;

	//time
	int current_step;
	float time_interval;
	float current_time;
	float end_time;
	float time_interval_factor;//time step size factor (<= 1.0)

	float energy_internal;
	float energy_kinetic;
	cinder::Vec3f momentum;
	cinder::Vec3f momentum_body1;
	cinder::Vec3f momentum_body2;

	void init_particles(int& num_particle)
	{
		//this->num_particle=num_particle;
		//particle_list.resize(this->num_particle,new particle());
	}

	void init_bodies(int& num_body)
	{
		//this->num_body=num_body;
		//body_list.resize(this->num_body, new particle_body());
	}

	//calculate momentum, kinematic energy  and internal energy
	void calcEnergy()
	{
		energy_internal=0;
		energy_kinetic=0;
		momentum.set(0,0,0);

		
		for(int i=0;i<num_particle;i++)
		{
			energy_kinetic+=particle_list[i]->velocity.dot(particle_list[i]->velocity)*particle_list[i]->particle_mass*0.5f;
			energy_internal+=particle_list[i]->internel_energy;
			momentum+=particle_list[i]->velocity*particle_list[i]->particle_mass;
		}
	}

	void init_method(const bool& MUSL=false,const bool& USL=false,const bool& USF=false,const bool& GIMP=false,const bool& contact=false,const bool& gravity=false)
	{
		this->MUSL=MUSL;
		this->USL=USL;
		this->USF=USF;
		this->GIMP=GIMP;
		this->contact=contact;
		this->gravity=gravity;
	}

	void init_time_param(const float& time_interval=0,const float& end_time=0,const float& time_interval_factor=0)
	{
		this->current_step=0;
		this->time_interval=time_interval;
		this->current_time=0;
		this->end_time=end_time;
		this->time_interval_factor=time_interval_factor;
	}
};
#endif