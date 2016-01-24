#pragma once
#include "FF_material_model.h"
class ElasticModel : public MaterialModel
{
public:
	void compute()//1
	{
		denergy_internal=stress_compoent[0]*dstrain[0]+
			stress_compoent[1]*dstrain[1]+
			stress_compoent[2]*dstrain[2]+
			stress_compoent[3]*dstrain[3]+
			stress_compoent[4]*dstrain[4]+
			stress_compoent[5]*dstrain[5];

		update_devitoric_stress_by_elastic_relation();
		update_mean_stress_by_elastic_relation();
		equivalent_stress=get_equivalent_stress();

		denergy_internal+=(deviatoric_stress[0]+mean_stress)*dstrain[0]+
			(deviatoric_stress[1]+mean_stress)*dstrain[1]+
			(deviatoric_stress[2]+mean_stress)*dstrain[2]+
			deviatoric_stress[3]*dstrain[3]+
			deviatoric_stress[4]*dstrain[4]+
			deviatoric_stress[5]*dstrain[5];
	}

	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity)
	{
		cout<<"ElasticModel"<<endl;
		update_stress_pre(p,mat,dstrain);
		//depends on model type
		update_rotate_stress(vorticity);
		compute();
		update_internal_energy();

		update_stress_post(p);
	}
};
