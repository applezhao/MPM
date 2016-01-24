#pragma once
#ifndef FF_GRID_HH
#define FF_GRID_HH
#include <vector>
#include <cinder/Vector.h>
#include <boost/multi_array.hpp>
#include "FF_particle.h"
using namespace std;


enum e_contact_type
{
	e_ct_no=0,
	e_ct_langrange=1,
	e_ct_penalty=2
};
enum e_grid_edge_condition
{
	e_gec_free=0,
	e_gec_fix=1,
	e_gec_symmetry=2
};

class grid_node
{
public:
	cinder::Vec3f position;//xg
	cinder::Vec3i is_fix;//fix_x y z edge condition
	grid_node()
	{
		is_fix.set(0,0,0);
	}
};

class grid_node_property
{
public:
	float mass;//mg
	cinder::Vec3f momentum;//pxg
	cinder::Vec3f force;// fxg internal/external force on grid node
};

class contact_grid_node_property
{
public:
	cinder::Vec3f normal;//ndir the normal direction of contact grid node
	cinder::Vec3f tangent;//sdir the tangential unit vetors of contact grid node
};

class grid_data
{
public:
	vector<grid_node*> grid_node_list;//node_list
	//vector<grid_node_property*> grid_node_property_list;//grid list
	boost::multi_array<grid_node_property*,2> grid_node_property_list;//grid_list
	//vector<contact_grid_node_property*> contact_grid_node_property_list;//cp_list
	boost::multi_array<contact_grid_node_property*,2> contact_grid_node_property_list;//cp_list

	float frictional_coefficient;//fricfa
	int flag_normal_compute;//normbody
	e_contact_type contact_type;
	cinder::Vec3f contact_total_force;//tot_cont_for;

	//grid
	cinder::Vec3f min_SPAN;
	cinder::Vec3f max_SPAN;

	float grid_intervel;//DCell
	float mass_cutoff;//cutoff value of grid mass
	vector<e_grid_edge_condition> fixS;

	int num_gridnode;
	cinder::Vec3i num_gridnode_axis;//NGX NGY NGZ
	int num_gridnode_xy;//ngxy

	//cell
	int num_cell;
	cinder::Vec3i num_cell_axis;
	int num_cell_xy;
	boost::multi_array<int,2> cell_node_list;

	int num_influence_node;
	vector<int> influenced_node_list;
	vector<cinder::Vec3f> distance_particle_grid;//rpg

	//for computing
	float iJacobi, iJacobi4;
	float igrid_intervel;//idcell
	vector<float> value_shape_func;//SHP
	vector<cinder::Vec3f> d_shape_func_axis;//DNDX DNDY DNDZ
	vector<cinder::Vec3i> const_coord;//snx sny snz

	grid_data();

	grid_data(const cinder::Vec3f& min_SPAN,const cinder::Vec3f& max_SPAN,const float& grid_intervel,const vector<e_grid_edge_condition>& fixS, 
		const bool& isGimp, const int& num_component,const float& mass_cutoff, const float& frictional_coefficient=0, const int& flag_normal_compute=0, 
		const e_contact_type& contact_type=e_ct_no);

	int in_which_cell(cinder::Vec3f& particle);

	void init_grid_and_cell_list(const int& num_component);
	void calc_shape_function(const grid_node* node, const particle* p, int compute_flag);

	void calc_shape_function_GIMP(const particle* p);

	
	void integrate_grid_momentum(const int& current_step, const float& time_interval);

	

	void apply_boundary_condition();



	void init_grid_list();

	void init_cell_list();

	//FOR GIMP ONLY
	void get_influenced_node_list(const particle* p, const int& cell_index);

	void update_influenced_node_list();


	void calc_sh_dn(const float& x, float& out_sh, float& out_dn);

	
	inline float sign(float x, float y)
	{
		if(y>0)
			return abs(x);
		else
			return -abs(x);
	}
};
#endif