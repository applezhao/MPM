#include <vector>
#include <cinder/Vector.h>
#include "FF_particle_data.h"
using namespace std;


void particle_body::set_material(const int& material_id, material* mat, MaterialModel* material_model, const int& component_id)
{
	this->material_id=material_id;
	this->mat=mat;
	this->material_model=material_model;
	this->component_id=component_id;
}

//calculate momentum, kinematic energy  and internal energy
void particle_data::calcEnergy()
{
	energy_internal=0;
	energy_kinetic=0;
	momentum.set(0,0,0);

		
	for(int i=0;i<num_particle;i++)
	{
		energy_kinetic+=particle_list[i]->velocity.dot(particle_list[i]->velocity)*particle_list[i]->particle_mass*0.5f;//0.5*mv^2
		energy_internal+=particle_list[i]->internel_energy;
		momentum+=particle_list[i]->velocity*particle_list[i]->particle_mass;//mv
	}
}

void particle_data::init_method(const bool& MUSL,const bool& USL,const bool& USF,const bool& GIMP,const bool& contact,const bool& gravity)
{
	this->MUSL=MUSL;
	this->USL=USL;
	this->USF=USF;
	this->GIMP=GIMP;
	this->contact=contact;
	this->gravity=gravity;
}

void particle_data::init_time_param(const float& time_interval,const float& end_time,const float& time_interval_factor)
{
	this->current_step=0;
	this->time_interval=time_interval;
	this->current_time=0;
	this->end_time=end_time;
	this->time_interval_factor=time_interval_factor;
}
