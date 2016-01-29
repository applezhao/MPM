#pragma once
#ifndef FF_MATERIAL_MODEL_HH
#define FF_MATERIAL_MODEL_HH
#include <vector>
#include <cinder/Vector.h>
#include "FF_material.h"
using namespace std;
class MaterialModel
{
public:
	int material_id;//mid
	float young_t;//young_ change

	float density_current;//den_

	float volume_0;//vol0_
	float volume_n;//vold
	float volume_current;//vol_
	float dvolume;//dvol volume increment*0.5

	vector<float> dstrain;//dinc strain increment
	float mean_stress;//sm
	float dmean_stress;// dsm ms increment
	vector<float> deviatoric_stress; //sd
	vector<float> stress_compoent; //sig
	vector<float> deviatoric_stress_n;//sold ds at step n
	float bulk_viscosity_force;//bqf

	float equivalent_stress; //seqv
	float effective_plastic_strain; //epeff_
	float deffective_plastic_strain;//depeff
	float yield_stress;//sig_y_
	float ratio;//for hardening calculation

	float G2, K3, plastic_modulus;

	float specifiedHeatCapacity;//specheat_
	float temperature;//tmprt

	float energy_internal;//iener
	float denergy_internal;//ieninc
	float energy_internal_per_initial_volume;//specen

	float mu;//mu=density_current/density_init-1;
	float relative_volume;//rv=volume_current/volume_0;
	float burn_fraction;//bfac
	float sound_speed;//cp

	bool failure;

	MaterialModel()
	{
		dstrain.resize(6);
		deviatoric_stress.resize(6);
		stress_compoent.resize(6);
		deviatoric_stress_n.resize(6);
		energy_internal=denergy_internal=energy_internal_per_initial_volume=burn_fraction=sound_speed=0;
	}
	virtual void model_compute(particle* p,material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel=0, const float& current_time=0, const float& grid_intervel=0)
	{cout<<"virtualcompute"<<endl;}

	void update_stress_pre(particle* p, material* mat,const vector<float>& dstrain);

	void update_stress_post(particle* p);
	void update_devitoric_stress_by_elastic_relation/*elastic_devi*/();

	void update_mean_stress_by_elastic_relation/*elastic_p*/();

	float get_equivalent_stress/*EquivalentStress*/();

	void update_rotate_stress/*sigrot*/(cinder::Vec3f& in_vorticity);

	void update_internal_energy/*lieupd*/();

	void update_internal_energy_EOS/*hieupd*/();

	void update_pressure_energy_EOS/*seleos*/(material* mat, bool& failure, const float& time_intervel=0, const float& grid_intervel=0);

	void update_pressure_energy_EOS1/*eos1*/(material* mat, bool& failure);

	void update_pressure_energy_EOS2/*eos1*/(material* mat, bool& failure, const float& time_intervel, const float& grid_intervel);

	void update_pressure_energy_EOS3/*eos1*/(material* mat, bool& failure, const float& time_intervel, const float& grid_intervel);

	void calc_bulk_viscosity/*bulkq*/(material* mat, const float& time_intervel, const float& grid_intervel);

	void null_model_compute()/*m3md7*/;

	void simplified_johnson_cook_model_compute(material* mat,const float& time_intervel)/*m3md5*/;
};

class materialData
{
public:
	int num_material;//nb_mat
	vector<material*> material_list;
	vector<MaterialModel*> moterial_model_list;
	bool isJaumannRate;//Jaum default true
	int num_DetonationPoint;//nDeto
	vector<cinder::Vec3f> position_DetonationPoint;//Deto
	float bq1;//artificial bulk viscosity coefficients
	float bq2;//artificial bulk viscosity coefficients

	materialData()
	{
		isJaumannRate=true;
		num_DetonationPoint=0;
		num_material=0;
		bq1=bq2=0;
	}

	void set_data(const bool& isJaumannRate, int& num_DetonationPoint, float& bq1, float& bq2, vector<cinder::Vec3f>& position_DetonationPoint)
	{
		this->isJaumannRate=isJaumannRate;
		this->num_DetonationPoint=num_DetonationPoint;
		this->bq1=bq1;
		this->bq2=bq2;

		this->position_DetonationPoint.resize(num_DetonationPoint);
		for(int i=0;i<num_DetonationPoint;i++)
			this->position_DetonationPoint[i]=position_DetonationPoint[i];

	}
};
#endif
