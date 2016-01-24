#pragma once
#ifndef FF_MATERIAL_HH
#define FF_MATERIAL_HH
#include <vector>
#include <cinder/Vector.h>
#include "FF_particle.h"
using namespace std;
enum e_material
{
	e_elastic=1,
	e_perfect_plastic=2,
	e_isotropic_harden_plastic=3,
	e_johnson_cook=4,
	e_simplified_johnson_cook=5,
	e_simplified_johnson_cook_with_failure=6,
	e_johnson_cook_with_fracture=7,
	e_high_explosive_burn=8,
	e_null=9,
	e_drucker_prager=10,
	e_undeformable=12
};

enum e_EOS
{
	e_linear_polynomial_eos=1,
	e_gruneisen_eos=2,
	e_JWL_eos=3
};
class material
{
public:
	e_material matType;
	e_EOS EOSType;
	float density;
	float young;
	float poisson;
	float particleMass;//Mp
	float yield0;
	float tangentialModulus;//TangMod
	float roomTemperature;//roomt
	float meltTemperature;//melt
	float specifiedHeatCapacity;//specHeat
	
	//for johnson-cook material
	float b_jc;
	float n_jc;
	float c_jc;
	float m_jc;

	//strain rate normalization factor used in johnson-cook model
	float epso;
	float prd;//tensile pressure to begin damage
	float epf;//tensile pressure and equivalent plastic strain at failure
	float detonationVelocity;//D
	vector<float> cEOS;//constants in equation of state
	float wavespeed;//wavespd

	//for Drucker-prager soil material
	float q_fai;
	float k_fai;
	float q_psi;
	float ten_f;

	//data for all materials. they are stored in material for easy access as they are const value in the while scene
	bool isJaumannRate;//Jaum default true
	int num_DetonationPoint;//nDeto
	vector<cinder::Vec3f> position_DetonationPoint;//Deto
	float bq1;//artificial bulk viscosity coefficients
	float bq2;//artificial bulk viscosity coefficients

	//set date for all materials
	void set_data_for_all_materials(const bool& isJaumannRate, int& num_DetonationPoint, float& bq1, float& bq2, vector<cinder::Vec3f>& position_DetonationPoint);

	material();
	//for all kinds of eos
	void set_linear_polynomial_eos(const float& c0, const float& c1, const float& c2, const float& c3, const float& c4, const float& c5, const float& c6, const float& e0);

	void set_gruneisen_eos(const float& c0, const float& lambda, const float& gamma0, const float& e0);

	void set_JWL_eos(const float& a, const float& b, const float& r1, const float& r2, const float& w, const float& e0);

	//set all kinds of material
	void set_elastic(const float& density, const float& young, const float&poisson);
	
	void set_perfect_plastic(const float& density, const float& young, const float&poisson, const float& yield0);

	void set_isotropic_harden_plastic(const float& density, const float& young, const float&poisson, const float& yield0, const float& tangentialModulus);

	void set_johnson_cook(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc, const float& m_jc, 
		const float& roomTemperature, const float& meltTemperature, const float& specifiedHeatCapacity, const float& epso);

	void set_simplified_johnson_cook(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc,  
		const float& epso);

	void set_simplified_johnson_cook_with_failure(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc,  
		const float& epso, const float& epf);

	void set_johnson_cook_with_failure(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc, const float& m_jc, 
		const float& roomTemperature, const float& meltTemperature, const float& specifiedHeatCapacity, const float& epso, const float& epf);

	void set_high_explosive_burn(const float& density, const float& detonationVelocity);

	void set_null(const float& density);

	void set_drucker_prager(const float& density, const float& young, const float&poisson, const float& q_fai, const float& k_fai, const float& q_psi, const float&  ten_f);
};


#endif