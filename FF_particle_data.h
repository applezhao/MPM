#pragma once
#ifndef FF_PARTICLE_DATA_HH
#define FF_PARTICLE_DATA_HH
#include <vector>
#include <cinder/Vector.h>
#include "FF_material.h"
#include "FF_material_model.h"
using namespace std;


class particle_body
{
public:
	material* mat;
	MaterialModel* material_model;
	int material_id;//mat
	int component_id;//comid
	cinder::Vec3f gravity;//Gravp
	int begin;//the beginning of particle of body  par_begin
	int end;//the end of particle of body par_end
	particle_body()
	{
		gravity.set(0,0,0);
	}

	void set_material(const int& material_id, material* mat, MaterialModel* material_model, const int& component_id);
};

class particle_data
{
public:
	int num_particle;//nb_particle
	vector<particle*> particle_list;
	int num_component;//nb_component
	int num_body;//nb_body
	vector<particle_body*> body_list;

	//method
	bool MUSL;
	bool USL;
	bool USF;
	bool GIMP;
	bool contact;
	bool gravity;

	//time
	int current_step;//istep
	float time_interval;//DT
	float current_time;//CurrentTime
	float end_time;//EndTime
	float time_interval_factor;//DTScale time step size factor (<= 1.0)

	float energy_internal;//EngInternal
	float energy_kinetic;//EngKinetic
	cinder::Vec3f momentum;//Momentum
	cinder::Vec3f momentum_body1;//Mombody1
	cinder::Vec3f momentum_body2;//Mombody2

	//calculate momentum, kinematic energy  and internal energy
	void calcEnergy();

	void init_method(const bool& MUSL=false,const bool& USL=false,const bool& USF=false,const bool& GIMP=false,const bool& contact=false,const bool& gravity=false);
	void init_time_param(const float& time_interval=0,const float& end_time=0,const float& time_interval_factor=0);
};
#endif