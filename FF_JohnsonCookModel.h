#pragma once
#include "FF_material_model.h"
class JohnsonCookModel : public MaterialModel
{
public:
	void compute()//4
	{
		
	}
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel)
	{
		cout<<"JohnsonCookModel"<<endl;
		update_stress_pre(p,mat,dstrain);

		//depends on model type
		update_rotate_stress(vorticity);
		compute();
		update_internal_energy();

		//check 0 divide
		p->celsius_temperature+=equivalent_stress*deffective_plastic_strain/density_current/mat->specifiedHeatCapacity;

		update_stress_post(p);
	}
};
