#include "FF_material.h"

//set date for all materials
void material::set_data_for_all_materials(const bool& isJaumannRate, int& num_DetonationPoint, float& bq1, float& bq2, vector<cinder::Vec3f>& position_DetonationPoint)
{
	this->isJaumannRate=isJaumannRate;
	this->num_DetonationPoint=num_DetonationPoint;
	this->bq1=bq1;
	this->bq2=bq2;

	this->position_DetonationPoint.resize(num_DetonationPoint);
	for(int i=0;i<num_DetonationPoint;i++)
		this->position_DetonationPoint[i]=position_DetonationPoint[i];

}

material::material()
{
	EOSType=e_linear_polynomial_eos;
	cEOS.resize(10);
	isJaumannRate=true;
	num_DetonationPoint=0;
	bq1=0;
	bq2=0;
}
//for all kinds of eos
void material::set_linear_polynomial_eos(const float& c0, const float& c1, const float& c2, const float& c3, const float& c4, const float& c5, const float& c6, const float& e0)
{
	EOSType=e_linear_polynomial_eos;
	cEOS[0]=c0;
	cEOS[1]=c1;
	cEOS[2]=c2;
	cEOS[3]=c3;
	cEOS[4]=c4;
	cEOS[5]=c5;
	cEOS[6]=c6;
	cEOS[9]=e0;
}

void material::set_gruneisen_eos(const float& c0, const float& lambda, const float& gamma0, const float& e0)
{
	EOSType=e_gruneisen_eos;
	cEOS[0]=density*c0*c0;
	cEOS[1]=cEOS[0]*(2*lambda-1);
	cEOS[2]=cEOS[0]*(lambda-1)*(3*lambda-1);
	cEOS[3]=density*gamma0;
	cEOS[4]=gamma0;
	cEOS[9]=e0;
}

void material::set_JWL_eos(const float& a, const float& b, const float& r1, const float& r2, const float& w, const float& e0)
{
	EOSType=e_JWL_eos;
	cEOS[0]=a;
	cEOS[1]=b;
	cEOS[2]=r1;
	cEOS[3]=r2;
	cEOS[4]=w;
	cEOS[9]=e0;
}

//set all kinds of material
void material::set_elastic(const float& density, const float& young, const float&poisson)
{
	matType=e_elastic;
	this->density=density;
	this->young=young;
	this->poisson=poisson;
}
	
void material::set_perfect_plastic(const float& density, const float& young, const float&poisson, const float& yield0)
{
	matType=e_perfect_plastic;
	this->density=density;
	this->young=young;
	this->poisson=poisson;
	this->yield0=yield0;
}

void material::set_isotropic_harden_plastic(const float& density, const float& young, const float&poisson, const float& yield0, const float& tangentialModulus)
{
	matType=e_isotropic_harden_plastic;
	this->density=density;
	this->young=young;
	this->poisson=poisson;
	this->yield0=yield0;
	this->tangentialModulus=tangentialModulus;
}

void material::set_johnson_cook(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc, const float& m_jc, 
	const float& roomTemperature, const float& meltTemperature, const float& specifiedHeatCapacity, const float& epso)
{
	matType=e_johnson_cook;
	this->density=density;
	this->young=young;
	this->poisson=poisson;
	this->yield0=yield0;
	this->b_jc=b_jc;
	this->n_jc=n_jc;
	this->c_jc=c_jc;
	this->m_jc=m_jc;
	this->roomTemperature=roomTemperature;
	this->meltTemperature=meltTemperature;
	this->specifiedHeatCapacity=specifiedHeatCapacity;
	this->epso=epso;
}

void material::set_simplified_johnson_cook(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc,  
	const float& epso)
{
	matType=e_simplified_johnson_cook;
	this->density=density;
	this->young=young;
	this->poisson=poisson;
	this->yield0=yield0;
	this->b_jc=b_jc;
	this->n_jc=n_jc;
	this->c_jc=c_jc;
	this->epso=epso;
}

void material::set_simplified_johnson_cook_with_failure(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc,  
	const float& epso, const float& epf)
{
	matType=e_simplified_johnson_cook_with_failure;
	this->density=density;
	this->young=young;
	this->poisson=poisson;
	this->yield0=yield0;
	this->b_jc=b_jc;
	this->n_jc=n_jc;
	this->c_jc=c_jc;
	this->epso=epso;
	this->epf=epf;
}

void material::set_johnson_cook_with_failure(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc, const float& m_jc, 
	const float& roomTemperature, const float& meltTemperature, const float& specifiedHeatCapacity, const float& epso, const float& epf)
{
	matType=e_johnson_cook_with_fracture;
	this->density=density;
	this->young=young;
	this->poisson=poisson;
	this->yield0=yield0;
	this->b_jc=b_jc;
	this->n_jc=n_jc;
	this->c_jc=c_jc;
	this->m_jc=m_jc;
	this->roomTemperature=roomTemperature;
	this->meltTemperature=meltTemperature;
	this->specifiedHeatCapacity=specifiedHeatCapacity;
	this->epso=epso;
	this->epf=epf;
}

void material::set_high_explosive_burn(const float& density, const float& detonationVelocity)
{
	matType=e_high_explosive_burn;
	this->density=density;
	this->detonationVelocity=detonationVelocity;
}

void material::set_null(const float& density)
{
	matType=e_null;
	this->density=density;
}

void material::set_drucker_prager(const float& density, const float& young, const float&poisson, const float& q_fai, const float& k_fai, const float& q_psi, const float&  ten_f)
{
	matType=e_drucker_prager;
	this->density=density;
	this->young=young;
	this->poisson=poisson;
	this->q_fai=q_fai;
	this->k_fai=k_fai;
	this->q_psi=q_psi;
	this->ten_f=ten_f;
}