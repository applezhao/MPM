#pragma once
#include "FF_material_model.h"
class JohnsonCookWithFractureModel : public MaterialModel//7
{
public:
	void compute/*m3dm8*/()
	{
		
	}
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel)
	{
		cout<<"JohnsonCookWithFractureModel"<<endl;
		update_stress_pre(p,mat,dstrain);

		//depends on model type
		update_rotate_stress(vorticity);
		if(!p->failure)
		{
			compute();
			update_pressure_energy_EOS(mat, failure, time_intervel, grid_intervel);
		}
		else
		{
			null_model_compute();
			update_pressure_energy_EOS(mat, failure, time_intervel, grid_intervel);
		}
		if(effective_plastic_strain>mat->epf)
		{
			effective_plastic_strain=mat->epf+0.0000001;
			failure=true;
		}
		if(failure&&mean_stress<0)
		{
			mean_stress=0;
			p->volume=volume_n;
		}
		p->failure=failure;

		update_stress_post(p);
	}
};
