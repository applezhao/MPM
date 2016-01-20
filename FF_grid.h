#pragma once
#ifndef FF_GRID_HH
#define FF_GRID_HH
#include <vector>
#include <cinder/Vector.h>
#include <boost/multi_array.hpp>
#include "FF_particle.h"
using namespace std;
class grid_node
{
public:
	cinder::Vec3f position;//xg
	cinder::Vec3i is_fix;//edge condition fix_y y z
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
class grid_data
{
public:
	vector<grid_node*> grid_node_list;//node list
	//vector<grid_node_property*> grid_node_property_list;//grid list
	boost::multi_array<grid_node_property*,2> grid_node_property_list;//grid list
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
	cinder::Vec3i num_gridnode_axis;
	int num_gridnode_xy;

	//cell
	int num_cell;
	cinder::Vec3i num_cell_axis;
	int num_cell_xy;
	boost::multi_array<int,2> cell_node_list;

	int num_influence_node;
	vector<int> influenced_node_list;
	vector<cinder::Vec3f> distance_particle_grid;//rpg

	//for computing
	float iJacobi, iJacobi4, igrid_intervel;
	vector<float> value_shape_func;
	vector<cinder::Vec3f> d_shape_func_axis;
	vector<cinder::Vec3i> const_coord;

	grid_data()
	{
		frictional_coefficient=0;
		flag_normal_compute=0;
		contact_type=e_ct_no;
		min_SPAN.set(0,0,0);
		max_SPAN.set(0,0,0);
		grid_intervel=0;
		mass_cutoff=0;
		fixS.resize(6,e_gec_free);
		
		num_gridnode=0;
		num_cell=0;
		num_cell_axis.set(0,0,0);
		num_cell_xy=0;

		num_influence_node=8;
		influenced_node_list.resize(27,0);
		distance_particle_grid.resize(27);

		const_coord.resize(8);
		const_coord[0].set(-1,-1,-1);
		const_coord[1].set(1,-1,-1);
		const_coord[2].set(1,1,-1);
		const_coord[3].set(-1,1,-1);
		const_coord[4].set(-1,-1,1);
		const_coord[5].set(1,-1,1);
		const_coord[6].set(1,1,1);
		const_coord[7].set(-1,1,1);
	}

	grid_data(const cinder::Vec3f& min_SPAN,const cinder::Vec3f& max_SPAN,const float& grid_intervel,const vector<e_grid_edge_condition>& fixS, const bool& isGimp, const int& num_component,const float& mass_cutoff, const float& frictional_coefficient=0, const int& flag_normal_compute=0, const e_contact_type& contact_type=e_ct_no)
	{
		this->frictional_coefficient=frictional_coefficient;
		this->flag_normal_compute=flag_normal_compute;
		this->contact_type=contact_type;
		this->mass_cutoff=mass_cutoff;

		num_influence_node=8;
		influenced_node_list.resize(27,0);
		distance_particle_grid.resize(27);

		const_coord.resize(8);
		const_coord[0].set(-1,-1,-1);
		const_coord[1].set(1,-1,-1);
		const_coord[2].set(1,1,-1);
		const_coord[3].set(-1,1,-1);
		const_coord[4].set(-1,-1,1);
		const_coord[5].set(1,-1,1);
		const_coord[6].set(1,1,1);
		const_coord[7].set(-1,1,1);

		this->grid_intervel=grid_intervel;
		this->iJacobi=2.0f/this->grid_intervel;
		this->iJacobi4=iJacobi*0.125f;
		this->igrid_intervel=1.0f/this->grid_intervel;

		if(isGimp)
		{
			value_shape_func.resize(27,0);
			d_shape_func_axis.resize(27);
		}
		else
		{
			value_shape_func.resize(8,0);
			d_shape_func_axis.resize(8);
		}
		
		this->fixS.resize(6,e_gec_free);
		#pragma omp parallel for
		for(int i=0;i<6;i++)
			this->fixS[i]=fixS[i];

		//cell
		this->num_cell_axis.set(((max_SPAN-min_SPAN)/this->grid_intervel+cinder::Vec3f(0.5f,0.5f,0.5f)));
		this->min_SPAN.set(min_SPAN);
		this->max_SPAN.set(min_SPAN+num_cell_axis*grid_intervel);
		this->num_cell_xy=num_cell_axis.x*num_cell_axis.y;
		this->num_cell=num_cell_xy*num_cell_axis.z;

		//grid
		this->num_gridnode_axis=num_cell_axis+cinder::Vec3i(1,1,1);
		this->num_gridnode_xy=num_gridnode_axis.x*num_gridnode_axis.y;
		this->num_gridnode=num_gridnode_axis.z*num_gridnode_xy;
	}

	int in_which_cell(cinder::Vec3f& particle)
	{
		if(particle.x<min_SPAN.x||particle.y<min_SPAN.y||particle.z<min_SPAN.z||
			particle.x>max_SPAN.x||particle.y>max_SPAN.y||particle.z>max_SPAN.z)
			return -1;

		int x=int((particle.x-min_SPAN.x)/grid_intervel);
		int y=int((particle.y-min_SPAN.y)/grid_intervel);
		int z=int((particle.z-min_SPAN.z)/grid_intervel);

		int ret=num_cell_xy*z+num_cell_axis.x*z+x;
		if(ret<0 || ret>=num_cell)
			return -1;
		return ret;
	}	

	void init_grid_and_cell_list(const int& num_component)
	{
		//init grid
		init_grid_list();

		//init cell
		init_cell_list();
		
		//init grid property list
		this->grid_node_property_list.resize(boost::extents[num_component][this->num_gridnode]);
		for(int i=0;i<grid_node_property_list.shape()[0];i++)
		{
			for(int j=0;j<grid_node_property_list.shape()[1];j++)
			{
				this->grid_node_property_list[i][j]=new grid_node_property();
			}
		}

		//init contact grid list
		this->contact_grid_node_property_list.resize(boost::extents[num_component][this->num_gridnode]);
		for(int i=0;i<contact_grid_node_property_list.shape()[0];i++)
		{
			for(int j=0;j<contact_grid_node_property_list.shape()[1];j++)
			{
				this->contact_grid_node_property_list[i][j]=new contact_grid_node_property();
			}
		}
	}

	void calc_shape_function(const grid_node* node, const particle* p, int compute_flag)
	{
		cinder::Vec3f nature_coord=(p->position_t-node->position)*iJacobi-cinder::Vec3f(1,1,1);//-1-1
		vector<cinder::Vec3f> s_nature_coord;
		#pragma omp parallel for
		for(int i=0;i<8;i++)
		{
			s_nature_coord[i]=nature_coord.operator*(const_coord[i])+cinder::Vec3f(1,1,1);
		}

		if(compute_flag==0||compute_flag==2)
		{
			#pragma omp parallel for
			for(int i=0;i<8;i++)
				value_shape_func[i]=s_nature_coord[i].x*s_nature_coord[i].y*s_nature_coord[i].z*0.125f;
		}
		if(compute_flag==1||compute_flag==2)
		{
			#pragma omp parallel for
			for(int i=0;i<8;i++)
			{
				d_shape_func_axis[i].x=const_coord[i].x*s_nature_coord[i].y*s_nature_coord[i].z*iJacobi4;
				d_shape_func_axis[i].y=const_coord[i].y*s_nature_coord[i].x*s_nature_coord[i].z*iJacobi4;
				d_shape_func_axis[i].z=const_coord[i].z*s_nature_coord[i].x*s_nature_coord[i].y*iJacobi4;
			}
		}
	}

	void calc_shape_function_GIMP(const particle* p, const int cell_index)
	{
		get_influenced_node_list(p, cell_index);

		#pragma omp parallel for
		for(int i=0;i<num_influence_node;i++)
		{
			cinder::Vec3f sh, dn;
			if(influenced_node_list[i]<num_gridnode&&influenced_node_list[i]>=0)
			{
				calc_sh_dn(distance_particle_grid[i].x,sh.x,dn.x);
				calc_sh_dn(distance_particle_grid[i].y,sh.y,dn.y);
				calc_sh_dn(distance_particle_grid[i].z,sh.z,dn.z);
			}
			value_shape_func[i]=sh.x*sh.y*sh.z;
			d_shape_func_axis[i].x=dn.x*sh.y*sh.z*igrid_intervel;
			d_shape_func_axis[i].y=dn.y*sh.x*sh.z*igrid_intervel;
			d_shape_func_axis[i].z=dn.z*sh.x*sh.y*igrid_intervel;
		}
	}

	
	void integrate_grid_momentum(const int& current_step, const float& time_interval)
	{
		for(int i=0;i<grid_node_property_list.shape()[0];i++)//component
		{
			for(int j=0;j<grid_node_property_list.shape()[1];j++)//grid node
			{
				grid_node_property* gnp=grid_node_property_list[i][j];
				if(current_step==1)
					gnp->momentum+=gnp->force*time_interval*0.5f;
				else
					gnp->momentum+=gnp->force*time_interval;
			}
		}
		apply_boundary_condition();
	}

	

	void apply_boundary_condition()
	{
		for(int i=0;i<grid_node_property_list.shape()[0];i++)//component
		{
			for(int j=0;j<grid_node_property_list.shape()[1];j++)//grid node
			{
				grid_node_property* gnp=grid_node_property_list[i][j];
				grid_node* gn=grid_node_list[j];
				if(gn->is_fix.x)
					gnp->momentum.x=0;
				if(gn->is_fix.y)
					gnp->momentum.y=0;
				if(gn->is_fix.z)
					gnp->momentum.z=0;
			}
		}
	}



	void init_grid_list()
	{
		grid_node_list.resize(num_gridnode, new grid_node());
		#pragma omp parallel for
		for(int z=0;z<num_gridnode_axis.z;z++)
		{
			for(int y=0;y<num_gridnode_axis.y;y++)
			{
				for(int x=0;x<num_gridnode_axis.x;x++)
				{
					int grid_index=num_gridnode_xy*z+num_gridnode_axis.x*y+x;
					grid_node_list[grid_index]->position.set(min_SPAN+cinder::Vec3f(x,y,z)*grid_intervel);

					if( (x==0&&fixS[0]==e_gec_fix) || (x==num_gridnode_axis.x-1&&fixS[1]==e_gec_fix) ||
						(y==0&&fixS[2]==e_gec_fix) || (y==num_gridnode_axis.y-1&&fixS[3]==e_gec_fix) ||
						(z==0&&fixS[4]==e_gec_fix) || (z==num_gridnode_axis.z-1&&fixS[5]==e_gec_fix) )
						grid_node_list[grid_index]->is_fix.set(1,1,1);

					if( (x==0&&fixS[0]==e_gec_symmetry) || (x==num_gridnode_axis.x-1&&fixS[1]==e_gec_symmetry))
						grid_node_list[grid_index]->is_fix.x=1;

					if( (y==0&&fixS[2]==e_gec_symmetry) || (y==num_gridnode_axis.y-1&&fixS[3]==e_gec_symmetry))
						grid_node_list[grid_index]->is_fix.y=1;

					if( (z==0&&fixS[4]==e_gec_symmetry) || (z==num_gridnode_axis.z-1&&fixS[5]==e_gec_symmetry))
						grid_node_list[grid_index]->is_fix.z=1;
				}
			}
		}
	}

	void init_cell_list()
	{
		cell_node_list.resize(boost::extents[num_cell][8]);
		#pragma omp parallel for
		for(int z=0;z<num_cell_axis.z;z++)
		{
			for(int y=0;y<num_cell_axis.y;y++)
			{
				for(int x=0;x<num_cell_axis.x;x++)
				{
					int cell_index=num_cell_xy*z+num_cell_axis.x*y+x;
					int grid1=z*num_gridnode_xy+y*num_gridnode_axis.x+x;
					int gird5=grid1+num_gridnode_xy;

					cell_node_list[cell_index][0]=grid1;
					cell_node_list[cell_index][1]=grid1+1;
					cell_node_list[cell_index][2]=grid1+1+num_gridnode_axis.x;
					cell_node_list[cell_index][3]=grid1+num_gridnode_axis.x;

					cell_node_list[cell_index][4]=gird5;
					cell_node_list[cell_index][5]=gird5+1;
					cell_node_list[cell_index][6]=gird5+1+num_gridnode_axis.x;
					cell_node_list[cell_index][7]=gird5+num_gridnode_axis.x;
				}
			}
		}
	}

	//FOR GIMP ONLY
	void get_influenced_node_list(const particle* p, const int& cell_index)
	{
		#pragma omp parallel for
		for(int i=0;i<8;i++)
		{
			influenced_node_list[i]=cell_node_list[cell_index][i];
			distance_particle_grid[i]=(p->position_t-grid_node_list[influenced_node_list[i]]->position)*igrid_intervel;
		}

		update_influenced_node_list();

		//update distance
		#pragma omp parallel for
		for(int i=8;i<num_influence_node;i++)
		{
			if(influenced_node_list[i]<num_gridnode&&influenced_node_list[i]>=0)
				distance_particle_grid[i]=p->position_t-grid_node_list[influenced_node_list[i]]->position;
			else
				distance_particle_grid[i].set(100,100,100);
			distance_particle_grid[i]*=igrid_intervel;
		}
	}

	void update_influenced_node_list()
	{
		if(distance_particle_grid[0].z<0.25 ||distance_particle_grid[0].z>0.75)
		{
			if(distance_particle_grid[0].z<0.25)
			{
				influenced_node_list[8]=influenced_node_list[0]-num_gridnode_xy;
				influenced_node_list[9]=influenced_node_list[1]-num_gridnode_xy;
				influenced_node_list[10]=influenced_node_list[2]-num_gridnode_xy;
				influenced_node_list[11]=influenced_node_list[3]-num_gridnode_xy;
			}
			else//(distance_particle_grid[0].z>0.75)
			{
				influenced_node_list[8]=influenced_node_list[4]+num_gridnode_xy;
				influenced_node_list[9]=influenced_node_list[5]+num_gridnode_xy;
				influenced_node_list[10]=influenced_node_list[6]+num_gridnode_xy;
				influenced_node_list[11]=influenced_node_list[7]+num_gridnode_xy;
			}

			if(distance_particle_grid[0].y<0.25)
			{
				influenced_node_list[12]=influenced_node_list[8]-num_gridnode_axis.x;
				influenced_node_list[13]=influenced_node_list[9]-num_gridnode_axis.x;
				influenced_node_list[14]=influenced_node_list[0]-num_gridnode_axis.x;
				influenced_node_list[15]=influenced_node_list[1]-num_gridnode_axis.x;
				influenced_node_list[16]=influenced_node_list[4]-num_gridnode_axis.x;
				influenced_node_list[17]=influenced_node_list[5]-num_gridnode_axis.x;

				if(distance_particle_grid[0].x<0.25)
				{
					influenced_node_list[18]=influenced_node_list[0]-1;
					influenced_node_list[19]=influenced_node_list[3]-1;
					influenced_node_list[20]=influenced_node_list[4]-1;
					influenced_node_list[21]=influenced_node_list[7]-1;
					influenced_node_list[22]=influenced_node_list[8]-1;
					influenced_node_list[23]=influenced_node_list[11]-1;
					influenced_node_list[24]=influenced_node_list[12]-1;
					influenced_node_list[25]=influenced_node_list[14]-1;
					influenced_node_list[26]=influenced_node_list[16]-1;
					num_influence_node=27;
				}
				else if(distance_particle_grid[0].x>0.75)
				{
					influenced_node_list[18]=influenced_node_list[1]+1;
					influenced_node_list[19]=influenced_node_list[2]+1;
					influenced_node_list[20]=influenced_node_list[5]+1;
					influenced_node_list[21]=influenced_node_list[6]+1;
					influenced_node_list[22]=influenced_node_list[9]+1;
					influenced_node_list[23]=influenced_node_list[10]+1;
					influenced_node_list[24]=influenced_node_list[13]+1;
					influenced_node_list[25]=influenced_node_list[15]+1;
					influenced_node_list[26]=influenced_node_list[17]+1;
					num_influence_node=27;
				}
				else //(x middle)
					num_influence_node=18;
			}
			else if(distance_particle_grid[0].y>0.75)
			{
				influenced_node_list[12]=influenced_node_list[10]+num_gridnode_axis.x;
				influenced_node_list[13]=influenced_node_list[11]+num_gridnode_axis.x;
				influenced_node_list[14]=influenced_node_list[2]+num_gridnode_axis.x;
				influenced_node_list[15]=influenced_node_list[3]+num_gridnode_axis.x;
				influenced_node_list[16]=influenced_node_list[6]+num_gridnode_axis.x;
				influenced_node_list[17]=influenced_node_list[7]+num_gridnode_axis.x;
				if(distance_particle_grid[0].x<0.25)
				{
					influenced_node_list[18]=influenced_node_list[0]-1;
					influenced_node_list[19]=influenced_node_list[3]-1;
					influenced_node_list[20]=influenced_node_list[4]-1;
					influenced_node_list[21]=influenced_node_list[7]-1;
					influenced_node_list[22]=influenced_node_list[8]-1;
					influenced_node_list[23]=influenced_node_list[11]-1;
					influenced_node_list[24]=influenced_node_list[13]-1;
					influenced_node_list[25]=influenced_node_list[15]-1;
					influenced_node_list[26]=influenced_node_list[17]-1;
					num_influence_node=27;
				}
				else if(distance_particle_grid[0].x>0.75)
				{
					influenced_node_list[18]=influenced_node_list[1]+1;
					influenced_node_list[19]=influenced_node_list[2]+1;
					influenced_node_list[20]=influenced_node_list[5]+1;
					influenced_node_list[21]=influenced_node_list[6]+1;
					influenced_node_list[22]=influenced_node_list[9]+1;
					influenced_node_list[23]=influenced_node_list[10]+1;
					influenced_node_list[24]=influenced_node_list[12]+1;
					influenced_node_list[25]=influenced_node_list[14]+1;
					influenced_node_list[26]=influenced_node_list[15]+1;
					num_influence_node=27;
				}
				else //(x middle)
					num_influence_node=18;
			}
			else //(y middle)
			{
				if(distance_particle_grid[0].x<0.25)
				{
					influenced_node_list[12]=influenced_node_list[8]-1;
					influenced_node_list[13]=influenced_node_list[11]-1;
					influenced_node_list[14]=influenced_node_list[0]-1;
					influenced_node_list[15]=influenced_node_list[3]-1;
					influenced_node_list[16]=influenced_node_list[4]-1;
					influenced_node_list[17]=influenced_node_list[7]-1;
					num_influence_node=18;
				}
				else if(distance_particle_grid[0].x>0.75)
				{
					influenced_node_list[12]=influenced_node_list[9]+1;
					influenced_node_list[13]=influenced_node_list[10]+1;
					influenced_node_list[14]=influenced_node_list[1]+1;
					influenced_node_list[15]=influenced_node_list[2]+1;
					influenced_node_list[16]=influenced_node_list[5]+1;
					influenced_node_list[17]=influenced_node_list[6]+1;
					num_influence_node=18;
				}
				else //(x middle)
					num_influence_node=12;
			}
		}
		else //(z middle)
		{
			if(distance_particle_grid[0].y<0.25)
			{
				influenced_node_list[8]=influenced_node_list[0]-num_gridnode_axis.x;
				influenced_node_list[9]=influenced_node_list[1]-num_gridnode_axis.x;
				influenced_node_list[10]=influenced_node_list[4]-num_gridnode_axis.x;
				influenced_node_list[11]=influenced_node_list[5]-num_gridnode_axis.x;

				if(distance_particle_grid[0].x<0.25)
				{
					influenced_node_list[12]=influenced_node_list[0]-1;
					influenced_node_list[13]=influenced_node_list[3]-1;
					influenced_node_list[14]=influenced_node_list[4]-1;
					influenced_node_list[15]=influenced_node_list[7]-1;
					influenced_node_list[16]=influenced_node_list[8]-1;
					influenced_node_list[17]=influenced_node_list[10]-1;
					num_influence_node=18;
				}
				else if(distance_particle_grid[0].x>0.75)
				{
					influenced_node_list[12]=influenced_node_list[1]+1;
					influenced_node_list[13]=influenced_node_list[2]+1;
					influenced_node_list[14]=influenced_node_list[5]+1;
					influenced_node_list[15]=influenced_node_list[6]+1;
					influenced_node_list[16]=influenced_node_list[9]+1;
					influenced_node_list[17]=influenced_node_list[11]+1;
					num_influence_node=18;
				}
				else //(x middle)
					num_influence_node=12;
			}
			else if(distance_particle_grid[0].y>0.75)
			{
				influenced_node_list[8]=influenced_node_list[2]+num_gridnode_axis.x;
				influenced_node_list[9]=influenced_node_list[3]+num_gridnode_axis.x;
				influenced_node_list[10]=influenced_node_list[6]+num_gridnode_axis.x;
				influenced_node_list[11]=influenced_node_list[7]+num_gridnode_axis.x;

				if(distance_particle_grid[0].x<0.25)
				{
					influenced_node_list[12]=influenced_node_list[0]-1;
					influenced_node_list[13]=influenced_node_list[3]-1;
					influenced_node_list[14]=influenced_node_list[4]-1;
					influenced_node_list[15]=influenced_node_list[7]-1;
					influenced_node_list[16]=influenced_node_list[9]-1;
					influenced_node_list[17]=influenced_node_list[11]-1;
					num_influence_node=18;
				}
				else if(distance_particle_grid[0].x>0.75)
				{
					influenced_node_list[12]=influenced_node_list[1]+1;
					influenced_node_list[13]=influenced_node_list[2]+1;
					influenced_node_list[14]=influenced_node_list[5]+1;
					influenced_node_list[15]=influenced_node_list[6]+1;
					influenced_node_list[16]=influenced_node_list[8]+1;
					influenced_node_list[17]=influenced_node_list[10]+1;
					num_influence_node=18;
				}
				else //(x middle)
					num_influence_node=12;
			}
			else //(y middle)
			{
				if(distance_particle_grid[0].x<0.25)
				{
					influenced_node_list[8]=influenced_node_list[0]-1;
					influenced_node_list[9]=influenced_node_list[3]-1;
					influenced_node_list[10]=influenced_node_list[4]-1;
					influenced_node_list[11]=influenced_node_list[7]-1;
					num_influence_node=12;
				}
				else if(distance_particle_grid[0].x>0.75)
				{
					influenced_node_list[8]=influenced_node_list[1]+1;
					influenced_node_list[9]=influenced_node_list[2]+1;
					influenced_node_list[10]=influenced_node_list[5]+1;
					influenced_node_list[11]=influenced_node_list[6]+1;
					num_influence_node=12;
				}
				else //(x middle)
					num_influence_node=8;
			}
		}//END IF
	}

	inline float sign(float x, float y)
	{
		if(y>0)
			return abs(x);
		else
			return -abs(x);
	}

	void calc_sh_dn(const float& x, float& out_sh, float& out_dn)
	{
		float absx=abs(x);
		if(absx<0.25)
		{
			out_sh=(7-16*absx*absx)*0.125f;
			out_dn=-4*x;
		}
		else if(absx<0.75)
		{
			out_sh=1-x;
			out_dn=-sign(1,x);
		}
		else if(absx<1.25)
		{
			out_sh=powf((5-4*x)*0.25,2);
			out_dn=2*x-sign(1,x)*2.5;
		}
		else
		{
			out_sh=0;
			out_dn=0;
		}
	}
};
#endif