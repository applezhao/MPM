#include "FF_simulator.h"
void simulator::initialize_explosion()
{
	int num_particle=125000;
	float endtime=0.015f;
	float outtime=endtime;

	cinder::Vec3f grid_span_min(0,0,0);
	cinder::Vec3f grid_span_max(3,3,3);
	float grid_space=0.05f;

	float time_interval_factor=0.5;
	float time_interval=1;
	outtime=0.001;
	float report_time=0.0001;

	vector<e_grid_edge_condition> fixS(6);
	fixS[0]=(e_grid_edge_condition)2;
	fixS[1]=(e_grid_edge_condition)0;
	fixS[2]=(e_grid_edge_condition)2;
	fixS[3]=(e_grid_edge_condition)2;
	fixS[4]=(e_grid_edge_condition)2;
	fixS[5]=(e_grid_edge_condition)2;

	int num_material=1;
	p_material_data=new materialData();
	p_material_data->num_material=num_material;
	p_material_data->material_list.resize(p_material_data->num_material);
	p_material_data->moterial_model_list.resize(p_material_data->num_material);
	for(int i=0;i<num_material;i++)
	{
		p_material_data->material_list[i]=new material();
		int materialtype=8;
		int eostype=3;
		switch(materialtype)
		{
		case 8:
			p_material_data->material_list[i]->set_high_explosive_burn(1.63*0.001,6930);
			p_material_data->moterial_model_list[i]=new HighExplosiveBurnModel();
		}

		switch(eostype)
		{
		case 3:
			p_material_data->material_list[i]->set_JWL_eos(3.712*100000, 3.21*1000, 4.15, 0.95, 0.3, 6993);
		}
	}

	float cutoff=0;
	for(int i=0;i<num_material;i++)
	{
		cutoff+=p_material_data->material_list[i]->density;
	}
	cutoff=cutoff*pow(grid_space,3)*0.00001f/num_material;

	bool jaum=1;
	int num_detonation_point=1;
	vector<cinder::Vec3f> detonation_points;
	detonation_points.push_back(cinder::Vec3f(0,0,0));
	//detonation_points.push_back(cinder::Vec3f(0,0.5,0.5));
	cinder::Vec2f bulk_viscosity(1.5f,0.06f);
	for(int i=0;i<num_material;i++)
	{
		p_material_data->material_list[i]->set_data_for_all_materials(jaum,num_detonation_point, bulk_viscosity.x, bulk_viscosity.y, detonation_points);
		p_material_data->set_data(jaum,num_detonation_point, bulk_viscosity.x, bulk_viscosity.y, detonation_points);
	}

	//get time intervel
	for(int i=0;i<num_material;i++)
	{
		material* m=p_material_data->material_list[i];
		if(m->matType!=e_undeformable)
		{
			float sound_speed=sqrt( m->young*(1-m->poisson)/(1+m->poisson)/(1-2*m->poisson)/m->density );
			sound_speed=max(sound_speed,m->detonationVelocity);
			sound_speed=max(sound_speed,m->wavespeed);
			time_interval=min(time_interval,grid_space/sound_speed);
		}
	}
	time_interval*=time_interval_factor;
	cout<<"time intervel:"<<time_interval<<endl;

		
	//particle&body
	p_particle_data=new particle_data();
	int body_counter=0;
	int particle_counter=0;

	int componentid=0;
	int materialid=0;

	particle_body* pb=new particle_body();
	pb->set_material(materialid,p_material_data->material_list[materialid], p_material_data->moterial_model_list[materialid], componentid); 
	pb->begin=particle_counter;
		

	float density_block=1.63*0.001;
	float distance_between_particles=0.025;
	float particle_mass=pow(density_block*distance_between_particles,3);
	cinder::Vec3f block_min(0,0,0);
	block_min+=distance_between_particles*0.5;
	cinder::Vec3i num_particle_axis(50,50,50);

	//p_particle_data->particle_list.resize(num_particle);
	for(int x=0;x<num_particle_axis.x;x++)
	{
		for(int y=0;y<num_particle_axis.y;y++)
		{
			for(int z=0;z<num_particle_axis.z;z++)
			{
				particle* p=new particle();
				p->particle_id=particle_counter;
				p->particle_mass=particle_mass;
				p->position_t.set(x*distance_between_particles+block_min.x,
					y*distance_between_particles+block_min.y,
					z*distance_between_particles+block_min.z);
				p_particle_data->particle_list.push_back(p);
				particle_counter++;

			}
		}
	}
	pb->end=particle_counter;
	p_particle_data->body_list.push_back(pb);
	body_counter++;



	bool MUSL=1;
	bool GIMP=true;
	bool contact=false;
	bool USF=false;
	bool USL=false;
	bool gravity=false;
	p_particle_data->init_method(MUSL, USL, USF, GIMP, contact, gravity);

	int num_component=1;
	int num_body=1;
	p_particle_data->num_body=num_body;
	p_particle_data->num_component=num_component;
	p_particle_data->num_particle=num_particle;

	cout<<"body:"<<p_particle_data->num_body<<" "<<p_particle_data->body_list.size()<<endl;
	cout<<"particle:"<<p_particle_data->num_particle<<" "<<p_particle_data->particle_list.size()<<endl;
	cout<<"comp:"<<p_particle_data->num_component<<endl;

	//init particles
	for(int b=0;b<num_body;b++)
	{
		for(int i=p_particle_data->body_list[b]->begin;i<p_particle_data->body_list[b]->end;i++)
		{
			particle* p=p_particle_data->particle_list[i];
			material* m=p_particle_data->body_list[b]->mat;
			p->volume=p->particle_mass/m->density;
			p->yield_stress=m->yield0;
			p->internel_energy=m->cEOS[0]*p->volume;
			p->lighting_time=10;
			for(int j=0;j<num_detonation_point;j++)
			{
				p->lighting_time=min(p->lighting_time,(p->position_t-detonation_points[j]).length()/m->detonationVelocity);
			}

		}
	}

	p_particle_data->init_time_param(time_interval, endtime, time_interval_factor);

	//grid
	p_grid_data=new grid_data(grid_span_min,grid_span_max,grid_space,fixS,GIMP,num_component,cutoff);
	p_grid_data->init_grid_and_cell_list(num_component);
}


void simulator::init_grid_momentum()
{
	//reset grid prpperty
	for(int i=0;i<p_grid_data->grid_node_property_list.shape()[0];i++)
	{
		for(int j=0;j<p_grid_data->grid_node_property_list.shape()[1];j++)
		{
			p_grid_data->grid_node_property_list[i][j]->mass=0;
			p_grid_data->grid_node_property_list[i][j]->momentum.set(0,0,0);
		}
	}
	int comID=0;
	for(int i=0;i<p_particle_data->num_body;i++)
	{
		if(p_particle_data->contact)
			comID=p_particle_data->body_list[i]->component_id;
		for(int p=p_particle_data->body_list[i]->begin;p<p_particle_data->body_list[i]->end;p++)
		{
			particle* pp=p_particle_data->particle_list[p];
			pp->icell=p_grid_data->in_which_cell(pp->position_t);
			if(pp->icell<0)
				continue;

			for(int ss=0;ss<8;ss++)
				p_grid_data->influenced_node_list[i]=p_grid_data->cell_node_list[pp->icell][ss];

			if(p_particle_data->GIMP)
			{
				p_grid_data->get_influenced_node_list(pp,pp->icell);
				p_grid_data->calc_shape_function_GIMP(pp,pp->icell);
			}
			else
				p_grid_data->calc_shape_function(p_grid_data->grid_node_list[p_grid_data->influenced_node_list[0]],pp,0);

			for(int m=0;m<p_grid_data->num_influence_node;m++)
			{
				if(p_grid_data->influenced_node_list[i]>=0&&p_grid_data->influenced_node_list[i]<p_grid_data->num_gridnode)
				{
					grid_node_property* gnp=p_grid_data->grid_node_property_list[comID][p_grid_data->influenced_node_list[i]];
					gnp->mass+=p_grid_data->value_shape_func[m]*pp->particle_mass;
					gnp->momentum+=pp->velocity*p_grid_data->value_shape_func[m]*pp->particle_mass;
				}
			}
		}
	}
}

void simulator::update_grid_momentum()
{
	//reset grid prpperty
	for(int i=0;i<p_grid_data->grid_node_property_list.shape()[0];i++)
	{
		for(int j=0;j<p_grid_data->grid_node_property_list.shape()[1];j++)
		{
			p_grid_data->grid_node_property_list[i][j]->force.set(0,0,0);
		}
	}

	if(p_particle_data->contact)
	{
		for(int i=0;i<p_grid_data->contact_grid_node_property_list.shape()[0];i++)
		{
			for(int j=0;j<p_grid_data->contact_grid_node_property_list.shape()[1];j++)
			{
				p_grid_data->contact_grid_node_property_list[i][j]->normal.set(0,0,0);
			}
		}
	}
			

	int comID=0;
	for(int i=0;i<p_particle_data->num_body;i++)
	{
		if(p_particle_data->contact)
			comID=p_particle_data->body_list[i]->component_id;
		for(int p=p_particle_data->body_list[i]->begin;p<p_particle_data->body_list[i]->end;p++)
		{
			particle* pp=p_particle_data->particle_list[p];
			pp->icell=p_grid_data->in_which_cell(pp->position_t);
			if(pp->icell<0)
				continue;

			cinder::Vec3f stress_axis=pp->deviatoric_stress_asix+cinder::Vec3f(pp->mean_stress,pp->mean_stress,pp->mean_stress);
			cinder::Vec3f stress_plane=pp->deviatoric_stress_plane;
			cinder::Vec3f external_force=pp->force_load;

			if(p_particle_data->gravity)
				external_force+=pp->particle_mass*p_particle_data->body_list[i]->gravity;

			for(int ss=0;ss<8;ss++)
				p_grid_data->influenced_node_list[i]=p_grid_data->cell_node_list[pp->icell][ss];

			if(p_particle_data->GIMP)
			{
				p_grid_data->get_influenced_node_list(pp,pp->icell);
				p_grid_data->calc_shape_function_GIMP(pp,pp->icell);
			}
			else
				p_grid_data->calc_shape_function(p_grid_data->grid_node_list[p_grid_data->influenced_node_list[0]],pp,2);

			//loop
			for(int m=0;m<p_grid_data->num_influence_node;m++)
			{
				if(p_grid_data->influenced_node_list[i]>=0&&p_grid_data->influenced_node_list[i]<p_grid_data->num_gridnode)
				{
					grid_node_property* gnp=p_grid_data->grid_node_property_list[comID][p_grid_data->influenced_node_list[i]];
					cinder::Vec3f internal_force;
					internal_force.x=-pp->volume*cinder::Vec3f(stress_axis.x,stress_plane.x,stress_plane.z).dot(p_grid_data->d_shape_func_axis[m]);
					internal_force.y=-pp->volume*cinder::Vec3f(stress_plane.x,stress_axis.y,stress_plane.y).dot(p_grid_data->d_shape_func_axis[m]);
					internal_force.z=-pp->volume*cinder::Vec3f(stress_plane.z,stress_plane.y,stress_axis.z).dot(p_grid_data->d_shape_func_axis[m]);

					cinder::Vec3f external_force_a=external_force*p_grid_data->value_shape_func[m];

					gnp->force+=internal_force+external_force_a;
					if(p_particle_data->contact)
					{
						p_grid_data->contact_grid_node_property_list[comID][p_grid_data->influenced_node_list[i]]->normal+=p_grid_data->d_shape_func_axis[m]*pp->particle_mass;
					}
				}
			}
		}
	}
}

void simulator::calc_grid_momentum_MUSL()
{
	//reset grid prpperty
	for(int i=0;i<p_grid_data->grid_node_property_list.shape()[0];i++)
	{
		for(int j=0;j<p_grid_data->grid_node_property_list.shape()[1];j++)
		{
			p_grid_data->grid_node_property_list[i][j]->momentum.set(0,0,0);
		}
	}
	int comID=0;
	for(int i=0;i<p_particle_data->num_body;i++)
	{
		if(p_particle_data->contact)
			comID=p_particle_data->body_list[i]->component_id;
		for(int p=p_particle_data->body_list[i]->begin;p<p_particle_data->body_list[i]->end;p++)
		{
			particle* pp=p_particle_data->particle_list[p];
			pp->icell=p_grid_data->in_which_cell(pp->position_t);
			if(pp->icell<0)
				continue;

			for(int ss=0;ss<8;ss++)
				p_grid_data->influenced_node_list[i]=p_grid_data->cell_node_list[pp->icell][ss];

			if(p_particle_data->GIMP)
			{
				p_grid_data->get_influenced_node_list(pp,pp->icell);
				p_grid_data->calc_shape_function_GIMP(pp,pp->icell);
			}
			else
				p_grid_data->calc_shape_function(p_grid_data->grid_node_list[p_grid_data->influenced_node_list[0]],pp,0);

			for(int m=0;m<p_grid_data->num_influence_node;m++)
			{
				if(p_grid_data->influenced_node_list[i]>=0&&p_grid_data->influenced_node_list[i]<p_grid_data->num_gridnode)
				{
					grid_node_property* gnp=p_grid_data->grid_node_property_list[comID][p_grid_data->influenced_node_list[i]];
					gnp->momentum+=pp->velocity*p_grid_data->value_shape_func[m]*pp->particle_mass;
				}
			}
		}
	}

	p_grid_data->apply_boundary_condition();
}

void simulator::update_particle_position_velocity()
{
	int comID=0;
	for(int i=0;i<p_particle_data->num_body;i++)
	{
		if(p_particle_data->contact)
			comID=p_particle_data->body_list[i]->component_id;
		for(int p=p_particle_data->body_list[i]->begin;p<p_particle_data->body_list[i]->end;p++)
		{
			particle* pp=p_particle_data->particle_list[p];
			pp->icell=p_grid_data->in_which_cell(pp->position_t);
			if(pp->icell<0)
				continue;
			cinder::Vec3f position_x=pp->position_t;
			cinder::Vec3f velocity_x(0,0,0);
			cinder::Vec3f acceleration_x(0,0,0);
			for(int ss=0;ss<8;ss++)
				p_grid_data->influenced_node_list[ss]=p_grid_data->cell_node_list[pp->icell][ss];

			if(p_particle_data->GIMP)
			{
				p_grid_data->get_influenced_node_list(pp,pp->icell);
				p_grid_data->calc_shape_function_GIMP(pp,pp->icell);
			}
			else
				p_grid_data->calc_shape_function(p_grid_data->grid_node_list[p_grid_data->influenced_node_list[0]],pp,2);

			for(int m=0;m<p_grid_data->num_influence_node;m++)
			{
				if(p_grid_data->influenced_node_list[i]>=0&&p_grid_data->influenced_node_list[i]<p_grid_data->num_gridnode)
				{
					grid_node_property* gnp=p_grid_data->grid_node_property_list[comID][p_grid_data->influenced_node_list[i]];
					if(gnp->mass>p_grid_data->mass_cutoff)
					{
						velocity_x+=gnp->momentum/gnp->mass*p_grid_data->value_shape_func[m];
						acceleration_x+=gnp->force/gnp->mass*p_grid_data->value_shape_func[m];
					}
				}
			}

			//update particle
			pp->position_t_1=position_x+velocity_x*p_particle_data->time_interval;
			if(p_particle_data->current_step==1)
				pp->velocity+=acceleration_x*p_particle_data->time_interval*0.5f;
			else
				pp->velocity+=acceleration_x*p_particle_data->time_interval;
			if(p_particle_data->USF)
				pp->position_t=pp->position_t_1;
		}
	}
}

void simulator::update_particle_stress()
{
	int comID=0;
	for(int i=0;i<p_particle_data->num_body;i++)
	{
		if(p_particle_data->contact)
			comID=p_particle_data->body_list[i]->component_id;
		for(int p=p_particle_data->body_list[i]->begin;p<p_particle_data->body_list[i]->end;p++)
		{
			particle* pp=p_particle_data->particle_list[p];
			pp->icell=p_grid_data->in_which_cell(pp->position_t);
			if(pp->icell<0)
				continue;
			vector<float> d_strain(6,0.f);
			cinder::Vec3f d_vorticity(0,0,0);
			for(int ss=0;ss<8;ss++)
				p_grid_data->influenced_node_list[ss]=p_grid_data->cell_node_list[pp->icell][ss];

			if(p_particle_data->GIMP)
			{
				p_grid_data->get_influenced_node_list(pp,pp->icell);
				p_grid_data->calc_shape_function_GIMP(pp,pp->icell);
			}
			else
				p_grid_data->calc_shape_function(p_grid_data->grid_node_list[p_grid_data->influenced_node_list[0]],pp,2);

			for(int m=0;m<p_grid_data->num_influence_node;m++)
			{
				if(p_grid_data->influenced_node_list[i]>=0&&p_grid_data->influenced_node_list[i]<p_grid_data->num_gridnode)
				{
					grid_node_property* gnp=p_grid_data->grid_node_property_list[comID][p_grid_data->influenced_node_list[i]];
					if(gnp->mass>p_grid_data->mass_cutoff)
					{
						cinder::Vec3f grid_velocity=gnp->momentum/gnp->mass;
						d_strain[0]+=grid_velocity.x*p_grid_data->d_shape_func_axis[m].x;
						d_strain[1]+=grid_velocity.y*p_grid_data->d_shape_func_axis[m].y;
						d_strain[2]+=grid_velocity.z*p_grid_data->d_shape_func_axis[m].z;
						d_strain[3]+=grid_velocity.z*p_grid_data->d_shape_func_axis[m].y+grid_velocity.y*p_grid_data->d_shape_func_axis[m].z;
						d_strain[4]+=grid_velocity.y*p_grid_data->d_shape_func_axis[m].x+grid_velocity.x*p_grid_data->d_shape_func_axis[m].y;
						d_strain[5]+=grid_velocity.y*p_grid_data->d_shape_func_axis[m].x+grid_velocity.x*p_grid_data->d_shape_func_axis[m].y;

						d_vorticity.x+=grid_velocity.z*p_grid_data->d_shape_func_axis[m].y-grid_velocity.y*p_grid_data->d_shape_func_axis[m].z;
						d_vorticity.y+=grid_velocity.x*p_grid_data->d_shape_func_axis[m].z-grid_velocity.z*p_grid_data->d_shape_func_axis[m].x;
						d_vorticity.z+=grid_velocity.y*p_grid_data->d_shape_func_axis[m].x-grid_velocity.x*p_grid_data->d_shape_func_axis[m].y;
					}
				}
			}

			d_vorticity=d_vorticity*p_particle_data->time_interval*0.5f;
			for(int kk=0;kk<6;kk++)
				d_strain[kk]*=p_particle_data->time_interval;
			p_particle_data->body_list[i]->material_model->model_compute(pp,p_particle_data->body_list[i]->mat,d_strain,d_vorticity, p_particle_data->time_interval, p_particle_data->current_time, p_grid_data->grid_intervel);
			if(!p_particle_data->USF)
				pp->position_t=pp->position_t_1;
		}
	}
}

void simulator::Lagr_NodContact()
{
	grid_node_property* gnp0;
	grid_node_property* gnp1;
	contact_grid_node_property* cgnp0;
	contact_grid_node_property* cgnp1;
	for(int i=0;i<p_grid_data->num_gridnode;i++)
	{
		cinder::Vec3f normal;
		cinder::Vec3f contact_force;
		gnp0=p_grid_data->grid_node_property_list[0][i];
		gnp1=p_grid_data->grid_node_property_list[1][i];
		cgnp0=p_grid_data->contact_grid_node_property_list[0][i];
		cgnp1=p_grid_data->contact_grid_node_property_list[1][i];
		if(p_grid_data->flag_normal_compute==0)
			normal=cgnp0->normal-cgnp1->normal;
		else if(p_grid_data->flag_normal_compute==1)
			normal=cgnp0->normal;
		else if(p_grid_data->flag_normal_compute==2)
			normal=-cgnp1->normal;
		normal.normalize();

		cgnp0->normal=normal;
		cgnp1->normal=-cgnp0->normal;

		float contact_criteria=0;
		if(gnp0->mass>p_grid_data->mass_cutoff&&gnp1->mass>p_grid_data->mass_cutoff)
			contact_criteria=normal.dot(gnp0->momentum*gnp1->mass-gnp1->momentum*gnp0->mass);
		if(contact_criteria>0)
		{
			float normal_force=contact_criteria/(gnp0->mass+gnp1->mass)/p_particle_data->time_interval;
			if(p_grid_data->frictional_coefficient>0)
			{
				contact_force=(gnp0->momentum*gnp1->mass-gnp1->momentum*gnp0->mass)/(gnp0->mass+gnp1->mass)/p_particle_data->time_interval;
				cinder::Vec3f tangent_stick_force=contact_force-normal_force*cgnp0->normal;
				float tangent_slip_force=p_grid_data->frictional_coefficient*fabs(normal_force);
				if(tangent_slip_force<tangent_stick_force.length())
					contact_force=normal_force*cgnp0->normal+tangent_slip_force*tangent_stick_force.normalized();
			}
			else
				contact_force=normal_force*cgnp0->normal;
		}

		gnp0->force+=-contact_force;
		gnp1->force+=contact_force;
		gnp0->momentum+=-contact_force*p_particle_data->time_interval;
		gnp1->momentum+=contact_force*p_particle_data->time_interval;
	}
		
}

//for each frame or time intervel
void simulator::update()
{
	test();
	p_particle_data->current_step++;
	p_particle_data->current_time+=p_particle_data->time_interval;
	p_particle_data->energy_internal=0;
	test();
	//cout<<"step0"<<endl;
	//p_particle_data->particle_list[0]->debug();
	//p_grid_data->init_grid_momentum(p_particle_data->num_body,p_particle_data);
	init_grid_momentum();
	//cout<<"step1"<<endl;
	//p_particle_data->particle_list[0]->debug();
	if(p_particle_data->USF)
	{
		p_grid_data->apply_boundary_condition();
		//p_grid_data->update_particle_stress(p_particle_data->num_body,p_particle_data);
		update_particle_stress();
	}

	//p_grid_data->update_grid_momentum(p_particle_data->num_body,p_particle_data);
	update_grid_momentum();
	p_grid_data->integrate_grid_momentum(p_particle_data->current_step,p_particle_data->time_interval);
	//cout<<"step2"<<endl;
	//p_particle_data->particle_list[0]->debug();

	if(p_grid_data->contact_type==e_ct_langrange)
		//p_grid_data->Lagr_NodContact(p_particle_data);
			Lagr_NodContact();
	//p_grid_data->update_particle_position_velocity(p_particle_data->num_body,p_particle_data);
	update_particle_position_velocity();
	//cout<<"step3"<<endl;
	//p_particle_data->particle_list[0]->debug();
	if(p_particle_data->MUSL)
	{
		//p_grid_data->calc_grid_momentum_MUSL(p_particle_data->num_body,p_particle_data);
		calc_grid_momentum_MUSL();
		p_grid_data->apply_boundary_condition();
	}
	//cout<<"step4"<<endl;
	//p_particle_data->particle_list[0]->debug();
	if(!p_particle_data->USF)
	{
		//p_grid_data->update_particle_stress(p_particle_data->num_body,p_particle_data);
		update_particle_stress();
	}
	//cout<<"step5"<<endl;
	//p_particle_data->particle_list[0]->debug();
	p_particle_data->calcEnergy();
}