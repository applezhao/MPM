#pragma once
#include "FF_material_model.h"
class IsotropicHardenPlasticModel : public MaterialModel
{
public:
	void compute()//3
	{
		
	}
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel)
	{
		cout<<"IsotropicHardenPlasticModel"<<endl;
		update_stress_pre(p,mat,dstrain);

		//depends on model type
		update_rotate_stress(vorticity);
		compute();
		update_internal_energy();

		update_stress_post(p);
	}
};
