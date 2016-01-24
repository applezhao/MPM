#pragma once
#include "FF_material_model.h"
class SimplifiedJohnsonCookModel : public MaterialModel//5
{
public:
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel)
	{
		cout<<"SimplifiedJohnsonCookModel"<<endl;
		update_stress_pre(p,mat,dstrain);

		//depends on model type
		update_rotate_stress(vorticity);
		simplified_johnson_cook_model_compute(mat,time_intervel);
		update_pressure_energy_EOS(mat, failure, time_intervel, grid_intervel);

		update_stress_post(p);
	}
};
