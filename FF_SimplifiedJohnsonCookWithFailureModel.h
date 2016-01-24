#pragma once
#include "FF_material_model.h"
class SimplifiedJohnsonCookWithFailureModel : public MaterialModel//6
{
public:
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel)
	{
		cout<<"SimplifiedJohnsonCookWithFailureModel"<<endl;
		update_stress_pre(p,mat,dstrain);

		//depends on model type
		update_rotate_stress(vorticity);
		if(!p->failure)
			simplified_johnson_cook_model_compute(mat,time_intervel);
		else
			null_model_compute();
		update_pressure_energy_EOS(mat, failure, time_intervel, grid_intervel);
		if(effective_plastic_strain>mat->epf)
		{
			effective_plastic_strain=mat->epf+0.0000001;
			p->failure=true;
		}

		update_stress_post(p);
	}
};
