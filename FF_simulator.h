#pragma once
#ifndef FF_SIMULATOR_HH
#define FF_SIMULATOR_HH
#include <vector>
#include <cinder/Vector.h>
#include <boost/multi_array.hpp>
#include "FF_grid.h"
using namespace std;
class simulator
{
public:
	grid_data* p_grid_data;
	particle_data* p_particle_data;
	materialData* p_material_data;

	simulator()
	{}

	void initialize();


	void init_grid_momentum();

	void update_grid_momentum();

	void calc_grid_momentum_MUSL();

	void update_particle_position_velocity();

	void update_particle_stress();

	void Lagr_NodContact();

	//for each frame or time intervel
	void update();
};
#endif