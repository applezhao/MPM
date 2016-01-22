#pragma once
#ifndef FF_MATERIAL_HH
#define FF_MATERIAL_HH
#include <vector>
#include <cinder/Vector.h>
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
	float tangentialModulus;
	float roomTemperature;
	float meltTemperature;
	float specifiedHeatCapacity;
	
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
	float wavespeed;

	//for Drucker-prager soil material
	float q_fai;
	float k_fai;
	float q_psi;
	float ten_f;

	//data for all materials. they are stored in material for easy access as they are const value in the while scene
	bool isJaumannRate;
	int num_DetonationPoint;
	vector<cinder::Vec3f> position_DetonationPoint;
	float bq1;//artificial bulk viscosity coefficients
	float bq2;//artificial bulk viscosity coefficients

	//set date for all materials
	void set_data_for_all_materials(const bool& isJaumannRate, int& num_DetonationPoint, float& bq1, float& bq2, vector<cinder::Vec3f>& position_DetonationPoint)
	{
		this->isJaumannRate=isJaumannRate;
		this->num_DetonationPoint=num_DetonationPoint;
		this->bq1=bq1;
		this->bq2=bq2;

		this->position_DetonationPoint.resize(num_DetonationPoint);
		for(int i=0;i<num_DetonationPoint;i++)
			this->position_DetonationPoint[i]=position_DetonationPoint[i];

	}

	material()
	{
		EOSType=e_linear_polynomial_eos;
		cEOS.resize(10);
		isJaumannRate=true;
		num_DetonationPoint=0;
		bq1=0;
		bq2=0;
	}
	//for all kinds of eos
	void set_linear_polynomial_eos(const float& c0, const float& c1, const float& c2, const float& c3, const float& c4, const float& c5, const float& c6, const float& e0)
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

	void set_gruneisen_eos(const float& c0, const float& lambda, const float& gamma0, const float& e0)
	{
		EOSType=e_gruneisen_eos;
		cEOS[0]=density*c0*c0;
		cEOS[1]=cEOS[0]*(2*lambda-1);
		cEOS[2]=cEOS[0]*(lambda-1)*(3*lambda-1);
		cEOS[3]=density*gamma0;
		cEOS[4]=gamma0;
		cEOS[9]=e0;
	}

	void set_JWL_eos(const float& a, const float& b, const float& r1, const float& r2, const float& w, const float& e0)
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
	void set_elastic(const float& density, const float& young, const float&poisson)
	{
		matType=e_elastic;
		this->density=density;
		this->young=young;
		this->poisson=poisson;
	}
	
	void set_perfect_plastic(const float& density, const float& young, const float&poisson, const float& yield0)
	{
		matType=e_perfect_plastic;
		this->density=density;
		this->young=young;
		this->poisson=poisson;
		this->yield0=yield0;
	}

	void set_isotropic_harden_plastic(const float& density, const float& young, const float&poisson, const float& yield0, const float& tangentialModulus)
	{
		matType=e_isotropic_harden_plastic;
		this->density=density;
		this->young=young;
		this->poisson=poisson;
		this->yield0=yield0;
		this->tangentialModulus=tangentialModulus;
	}

	void set_johnson_cook(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc, const float& m_jc, 
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

	void set_simplified_johnson_cook(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc,  
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

	void set_simplified_johnson_cook_with_failure(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc,  
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

	void set_johnson_cook_with_failure(const float& density, const float& young, const float&poisson, const float& yield0, const float& b_jc, const float& n_jc, const float&  c_jc, const float& m_jc, 
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

	void set_high_explosive_burn(const float& density, const float& detonationVelocity)
	{
		matType=e_high_explosive_burn;
		this->density=density;
		this->detonationVelocity=detonationVelocity;
	}

	void set_null(const float& density)
	{
		matType=e_null;
		this->density=density;
	}

	void set_drucker_prager(const float& density, const float& young, const float&poisson, const float& q_fai, const float& k_fai, const float& q_psi, const float&  ten_f)
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
};

class particle
{
public:
	cinder::Vec3f position_t_1; //xx
	cinder::Vec3f position_t;//xp
	cinder::Vec3f velocity;//vxp
	cinder::Vec3f force_load;//fxp

	float volume;//vol
	float yield_stress;//sig_y
	float mean_stress;//sm
	float mises_stress;//seqv
	cinder::Vec3f deviatoric_stress_asix;//sdxx yy zz
	cinder::Vec3f deviatoric_stress_plane;//sdxy yz xz
	float effective_plastic_strain;//epeff
	float celsius_temperature;//celsius_t
	
	bool skipthis;//for plot less
	bool failure;
	
	int icell;//cell number
	float damage;//dmg
	float lighting_time;//lighting time for explosive simulation lt
	float internel_energy;//ie
	float particle_mass;//mass
	float sound_speed;//cp

	//debug
	int particle_id;

	particle()
	{
		mean_stress=mises_stress=effective_plastic_strain=damage=internel_energy=0;
		celsius_temperature=293.0;
		skipthis=false;
		failure=false;
		icell=0;
		velocity.set(0,0,0);
		force_load.set(0,0,0);
		deviatoric_stress_asix.set(0,0,0);
		deviatoric_stress_plane.set(0,0,0);
	}

	void debug()
	{
		cout<<velocity.x<<" "<<velocity.y<<" "<<velocity.z<<endl;
		cout<<particle_mass<<" "<<volume<<" "<<celsius_temperature<<" "<<icell<<endl;
		cout<<mean_stress<<" "<<yield_stress<<" "<<mises_stress<<" "<<effective_plastic_strain<<" "<<internel_energy<<endl;
	}
};



class MaterialModel
{
public:

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
	float effective_plastic_strain; //epeff
	float deffective_plastic_strain;//depeff
	float yield_stress;//sig_y
	float ratio;//for hardening calculation

	float G2, K3, plastic_modulus;

	float specifiedHeatCapacity;
	float temperature;

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
	}
	virtual void model_compute(particle* p,material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel=0, const float& current_time=0, const float& grid_intervel=0)
	{cout<<"virtualcompute"<<endl;}

	void update_stress_pre(particle* p, material* mat, vector<float>& dstrain)
	{
		for(int i=0;i<6;i++)
			this->dstrain[i]=dstrain[i];

		volume_n=p->volume;
		p->volume=volume_n*(1+dstrain[0]+dstrain[1]+dstrain[2]);
		volume_current=p->volume;
		dvolume=(volume_current-volume_n)*0.5f;
		BOOST_ASSERT_MSG(volume_current>=0, "volume is negative");

		effective_plastic_strain=p->effective_plastic_strain;
		yield_stress=p->yield_stress;
		equivalent_stress=p->mises_stress;
		//lighting_time=p->lighting_time;
		temperature=p->celsius_temperature;
		energy_internal=p->internel_energy;
		//mass=p->particle_mass;
		failure=p->failure;
		//density_0=body->mat->density;

		volume_0=p->particle_mass/mat->density;
		density_current=p->particle_mass/volume_current;
		mu=density_current/mat->density-1;
		relative_volume=volume_current/volume_0;
		energy_internal_per_initial_volume=energy_internal/volume_0;

		float young=mat->young*(1-p->damage);
		G2=young/(1+mat->poisson);
		K3=young/(1-2*mat->poisson);
		plastic_modulus=young*mat->tangentialModulus/(young-mat->tangentialModulus);

		mean_stress=p->mean_stress;

		deviatoric_stress[0]=p->deviatoric_stress_asix.x;
		deviatoric_stress[1]=p->deviatoric_stress_asix.y;
		deviatoric_stress[2]=p->deviatoric_stress_asix.z;
		deviatoric_stress[3]=p->deviatoric_stress_plane.x;
		deviatoric_stress[4]=p->deviatoric_stress_plane.y;
		deviatoric_stress[5]=p->deviatoric_stress_plane.z;

		//cauchy stress
		stress_compoent[0]=deviatoric_stress[0]+mean_stress;
		stress_compoent[1]=deviatoric_stress[1]+mean_stress;
		stress_compoent[2]=deviatoric_stress[2]+mean_stress;
		stress_compoent[3]=deviatoric_stress[3];
		stress_compoent[4]=deviatoric_stress[4];
		stress_compoent[5]=deviatoric_stress[5];

		for(int i=0;i<6;i++)
			deviatoric_stress_n[i]=deviatoric_stress[i];
	}

	void update_stress_post(particle* p)
	{
		p->mean_stress=mean_stress;
		p->mises_stress=equivalent_stress;
		p->deviatoric_stress_asix.set(deviatoric_stress[0],deviatoric_stress[1],deviatoric_stress[2]);
		p->deviatoric_stress_plane.set(deviatoric_stress[3],deviatoric_stress[4],deviatoric_stress[5]);
		p->yield_stress=yield_stress;
		p->effective_plastic_strain=effective_plastic_strain;
		p->internel_energy=energy_internal;
		p->sound_speed=sound_speed;
	}
	void update_devitoric_stress_by_elastic_relation/*elastic_devi*/()
	{
		float temp=(dstrain[0]+dstrain[1]+dstrain[2])/3;
		deviatoric_stress[0]+=G2*(dstrain[0]-temp);
		deviatoric_stress[1]+=G2*(dstrain[1]-temp);
		deviatoric_stress[2]+=G2*(dstrain[2]-temp);
		deviatoric_stress[3]+=G2*0.5f*dstrain[3];
		deviatoric_stress[4]+=G2*0.5f*dstrain[4];
		deviatoric_stress[5]+=G2*0.5f*dstrain[5];
	}

	void update_mean_stress_by_elastic_relation/*elastic_p*/()
	{
		float temp=(dstrain[0]+dstrain[1]+dstrain[2])/3;
		dmean_stress=temp*K3;
		mean_stress+=dmean_stress;
	}

	float get_equivalent_stress/*EquivalentStress*/()
	{
		float second_stree_invariant=0.5*(deviatoric_stress[0]*deviatoric_stress[0]+deviatoric_stress[1]*deviatoric_stress[1]+deviatoric_stress[2]*deviatoric_stress[2])+
			deviatoric_stress[3]*deviatoric_stress[3]+deviatoric_stress[4]*deviatoric_stress[4]+deviatoric_stress[5]*deviatoric_stress[5];
		return sqrtf(second_stree_invariant*3.0f);
	}

	void update_rotate_stress/*sigrot*/(cinder::Vec3f& in_vorticity)
	{
		vector<float> rot(6);
		vector<float> q(3);
		q[0]=2*stress_compoent[5]*in_vorticity.z;
		q[1]=2*stress_compoent[4]*in_vorticity.y;
		q[2]=2*stress_compoent[3]*in_vorticity.x;

		rot[0]=-q[0]+q[1];
		rot[1]=q[0]-q[2];
		rot[2]=-q[1]+q[2];
		rot[3]=in_vorticity.x*(stress_compoent[1]-stress_compoent[2])+in_vorticity.z*stress_compoent[4]-in_vorticity.y*stress_compoent[5];
		rot[4]=in_vorticity.y*(stress_compoent[2]-stress_compoent[0])+in_vorticity.x*stress_compoent[5]-in_vorticity.z*stress_compoent[3];
		rot[5]=in_vorticity.z*(stress_compoent[0]-stress_compoent[1])+in_vorticity.y*stress_compoent[3]-in_vorticity.x*stress_compoent[4];

		//update
		for(int i=0;i<6;i++)
			stress_compoent[i]+=rot[i];

		mean_stress=(stress_compoent[0]+stress_compoent[1]+stress_compoent[2])/3.0f;

		deviatoric_stress[0]=stress_compoent[0]-mean_stress;
		deviatoric_stress[1]=stress_compoent[1]-mean_stress;
		deviatoric_stress[2]=stress_compoent[2]-mean_stress;
		deviatoric_stress[3]=stress_compoent[3];
		deviatoric_stress[4]=stress_compoent[4];
		deviatoric_stress[5]=stress_compoent[5];

	}

	void update_internal_energy/*lieupd*/()
	{
		energy_internal+=0.25*denergy_internal*(volume_current+volume_n);
	}

	void update_internal_energy_EOS/*hieupd*/()
	{
		denergy_internal=0;
		for(int i=0;i<6;i++)
			denergy_internal+=dstrain[i]*(deviatoric_stress_n[i]+deviatoric_stress[i]);

		energy_internal+=0.25*denergy_internal*(volume_current+volume_n)+mean_stress*dvolume;
		energy_internal_per_initial_volume=energy_internal/volume_0;
	}

	void update_pressure_energy_EOS/*seleos*/(material* mat, bool& failure, const float& time_intervel=0, const float& grid_intervel=0)
	{
		switch((int)mat->EOSType)
		{
		case 1:
			update_pressure_energy_EOS1(mat,failure);
			break;
		case 2:
			update_pressure_energy_EOS2(mat,failure,time_intervel, grid_intervel);
			break;
		case 3:
			update_pressure_energy_EOS3(mat,failure);
			break;
		}
	}

	void update_pressure_energy_EOS1/*eos1*/(material* mat, bool& failure)
	{
		float A=mat->cEOS[0]+mu*(mat->cEOS[1]+mu*(mat->cEOS[2]+mu*mat->cEOS[3]));
		float B=mat->cEOS[4]+mu*(mat->cEOS[5]+mu*mat->cEOS[6]);
		update_internal_energy_EOS();
		float p_new=(A+B*energy_internal_per_initial_volume)/(1+B*dvolume/volume_0);
		if(p_new<0&&failure)
			p_new=0;
		energy_internal-=dvolume*p_new;
		mean_stress=-p_new;
	}

	void update_pressure_energy_EOS2/*eos1*/(material* mat, bool& failure, const float& time_intervel, const float& grid_intervel)
	{
		float A,B;
		if(mu>0)
		{
			A=(mat->cEOS[0]*mu+mat->cEOS[1]*mu*mu+mat->cEOS[2]*mu*mu*mu)*(1-mat->cEOS[3]*mu*0.5/density_current);
			B=mat->cEOS[4];
		}
		else
		{
			A=mat->cEOS[0]*mu;
			B=0;
		}
		//bulkq();
		calc_bulk_viscosity(mat, time_intervel, grid_intervel);
		update_internal_energy_EOS();
		float p_new=(A+B*energy_internal_per_initial_volume)/(1+B*dvolume/volume_0);
		if(p_new<0&&failure)
			p_new=0;
		energy_internal-=dvolume*p_new;
		mean_stress=-p_new;
	}

	void update_pressure_energy_EOS3/*eos1*/(material* mat, bool& failure)
	{
		//cout<<"step5.2.1"<<endl;
		//cout<<mean_stress<<" "<<energy_internal<<" "<<sound_speed<<endl;
		float er1v=exp(-mat->cEOS[2]*relative_volume);
		float er2v=exp(-mat->cEOS[3]*relative_volume);
		float wdr1v=mat->cEOS[0]-mat->cEOS[0]*mat->cEOS[4]/(mat->cEOS[2]*relative_volume);
		float wdr2v=mat->cEOS[1]-mat->cEOS[1]*mat->cEOS[4]/(mat->cEOS[3]*relative_volume);

		float A=wdr1v*er1v+wdr2v*er2v+bulk_viscosity_force;
		float B=mat->cEOS[4]/relative_volume;
		update_internal_energy_EOS();
		//cout<<"step5.2.2"<<endl;
		//cout<<mean_stress<<" "<<energy_internal<<" "<<sound_speed<<endl;
		A*=burn_fraction;
		B*=burn_fraction;
		float p_new=(A+B*energy_internal_per_initial_volume)/(1+B*dvolume/volume_0);
		if(p_new<0&&failure)
			p_new=0;
		//cout<<dvolume<<" "<<p_new<<" "<<volume_0<<" "<<energy_internal_per_initial_volume<<endl;
		//cout<<burn_fraction<<endl;
		//cout<<A<<" "<<B<<endl;
		energy_internal-=dvolume*p_new;
		mean_stress=-p_new;
		//cout<<"step5.2.3"<<endl;
		//cout<<mean_stress<<" "<<energy_internal<<" "<<sound_speed<<endl;

	}

	void calc_bulk_viscosity/*bulkq*/(material* mat, const float& time_intervel, const float& grid_intervel)
	{
		float bulk_strain_rate=(dstrain[0]+dstrain[1]+dstrain[2])/time_intervel;
		bulk_strain_rate=min(bulk_strain_rate,0.0f);
		bulk_viscosity_force=density_current*grid_intervel*grid_intervel*mat->bq1*bulk_strain_rate*bulk_strain_rate-mat->bq2*density_current*grid_intervel*sound_speed*bulk_strain_rate;
	}

	void null_model_compute()/*m3md7*/
	{
		for(int i=0;i<6;i++)
			deviatoric_stress[i]=0;
	}

	void simplified_johnson_cook_model_compute(material* mat,const float& time_intervel)/*m3md5*/
	{
		plastic_modulus=mat->b_jc*mat->n_jc*pow(effective_plastic_strain+0.0001, mat->n_jc-1);
		update_devitoric_stress_by_elastic_relation();
		equivalent_stress=get_equivalent_stress();

		if(equivalent_stress>yield_stress)
		{
			deffective_plastic_strain=(equivalent_stress-yield_stress)/(1.5f*G2+plastic_modulus);
			effective_plastic_strain+=deffective_plastic_strain;
			yield_stress=( mat->yield0+mat->b_jc*pow(effective_plastic_strain, mat->n_jc) )*
				(1+mat->c_jc*log(deffective_plastic_strain/mat->epso/time_intervel));
			ratio=yield_stress/equivalent_stress;

			for(int i=0;i<6;i++)
				deviatoric_stress[i]*=ratio;

			equivalent_stress*=ratio;
		}
	}
};

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

class PerfectPlasticModel : public MaterialModel
{
public:
	void compute()//2
	{
		
	}
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel, const float& bq1, const float& bq2)
	{
		cout<<"PerfectPlasticModel"<<endl;
		update_stress_pre(p,mat,dstrain);

		//depends on model type
		update_rotate_stress(vorticity);
		compute();
		update_internal_energy();

		update_stress_post(p);
	}
};

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
		p->celsius_temperature+=equivalent_stress*deffective_plastic_strain/density_current/mat->specifiedHeatCapacity;

		update_stress_post(p);
	}
};

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

class NullModel : public MaterialModel //9
{
public:
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel)
	{
		cout<<"NullModel"<<endl;
		update_stress_pre(p,mat,dstrain);

		//depends on model type
		null_model_compute();
		update_pressure_energy_EOS(mat, failure, time_intervel, grid_intervel);

		update_stress_post(p);
	}
};

class DruckerPragerModel : public MaterialModel //10
{
public:
	void compute/*m3dm9*/()
	{
		
	}
	virtual void model_compute(particle* p, material* mat, vector<float>& dstrain, cinder::Vec3f& vorticity, const float& time_intervel, const float& current_time, const float& grid_intervel)
	{
		cout<<"DruckerPragerModel"<<endl;
		update_stress_pre(p,mat,dstrain);

		//depends on model type
		update_rotate_stress(vorticity);
		compute();
		update_internal_energy();

		update_stress_post(p);
	}
};

class materialData
{
public:
	int num_material;
	vector<material*> material_list;
	vector<MaterialModel*> moterial_model_list;
	bool isJaumannRate;
	int num_DetonationPoint;
	vector<cinder::Vec3f> position_DetonationPoint;
	float bq1;//artificial bulk viscosity coefficients
	float bq2;//artificial bulk viscosity coefficients

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