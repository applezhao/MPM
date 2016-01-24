#pragma once
#include "FF_material_model.h"
class HighExplosiveBurnModel : public MaterialModel //8
{
public:
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel)
	{
		//cout<<"HighExplosiveBurnModel"<<endl;
		update_stress_pre(p,mat,dstrain);
		/*if(p->particle_id==0)
		{
			cout<<"step5.1"<<endl;
			p->debug();
			cout<<mean_stress<<" "<<energy_internal<<" "<<sound_speed<<endl;
			cout<<current_time<<" "<<p->lighting_time<<endl;
		}*/
		//depends on model type
		if(current_time > p->lighting_time)
			burn_fraction=1;
		else
			burn_fraction=0;
		
		update_pressure_energy_EOS(mat, failure, time_intervel, grid_intervel);
		
		update_stress_post(p);
		
	}
};
