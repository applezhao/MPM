#include "FF_material_model.h"
#include <boost/assert.hpp>


void MaterialModel::update_stress_pre(particle* p, material* mat, const vector<float>& dstrain)
{
	for(int i=0;i<6;i++)
		this->dstrain[i]=dstrain[i];

	volume_n=p->volume;
	p->volume=volume_n*(1+dstrain[0]+dstrain[1]+dstrain[2]);
	volume_current=p->volume;
	dvolume=(volume_current-volume_n)*0.5f;
	//BOOST_ASSERT_MSG(volume_current>=0, "volume is negative");
	if(volume_current<0)
		cout<<p->particle_id<<" "<<"volume is negative"<<endl;

	effective_plastic_strain=p->effective_plastic_strain;
	yield_stress=p->yield_stress;
	equivalent_stress=p->mises_stress;
	//lighting_time=p->lighting_time;
	temperature=p->celsius_temperature;
	energy_internal=p->internel_energy;
	//mass=p->particle_mass;
	failure=p->failure;
	//density_0=body->mat->density;
	young_t=mat->young;
	volume_0=p->particle_mass/mat->density;
	density_current=p->particle_mass/volume_current;
	mu=density_current/mat->density-1;
	relative_volume=volume_current/volume_0;
	energy_internal_per_initial_volume=energy_internal/volume_0;

	young_t=young_t*(1-p->damage);
	G2=young_t/(1+mat->poisson);
	K3=young_t/(1-2*mat->poisson);
	plastic_modulus=young_t*mat->tangentialModulus/(young_t-mat->tangentialModulus);

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

void MaterialModel::update_stress_post(particle* p)
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
void MaterialModel::update_devitoric_stress_by_elastic_relation/*elastic_devi*/()//elastic model
{
	float temp=(dstrain[0]+dstrain[1]+dstrain[2])/3;
	deviatoric_stress[0]+=G2*(dstrain[0]-temp);
	deviatoric_stress[1]+=G2*(dstrain[1]-temp);
	deviatoric_stress[2]+=G2*(dstrain[2]-temp);
	deviatoric_stress[3]+=G2*0.5f*dstrain[3];
	deviatoric_stress[4]+=G2*0.5f*dstrain[4];
	deviatoric_stress[5]+=G2*0.5f*dstrain[5];
}

void MaterialModel::update_mean_stress_by_elastic_relation/*elastic_p*/()//elastic model
{
	float temp=(dstrain[0]+dstrain[1]+dstrain[2])/3;
	dmean_stress=temp*K3;
	mean_stress+=dmean_stress;
}

float MaterialModel::get_equivalent_stress/*EquivalentStress*/()
{
	float second_stree_invariant=0.5*(deviatoric_stress[0]*deviatoric_stress[0]+deviatoric_stress[1]*deviatoric_stress[1]+deviatoric_stress[2]*deviatoric_stress[2])+
		deviatoric_stress[3]*deviatoric_stress[3]+deviatoric_stress[4]*deviatoric_stress[4]+deviatoric_stress[5]*deviatoric_stress[5];
	return sqrtf(second_stree_invariant*3.0f);
}

void MaterialModel::update_rotate_stress/*sigrot*/(cinder::Vec3f& in_vorticity)
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

void MaterialModel::update_internal_energy/*lieupd*/()
{
	energy_internal+=0.25*denergy_internal*(volume_current+volume_n);
}

void MaterialModel::update_internal_energy_EOS/*hieupd*/()
{
	denergy_internal=0;
	for(int i=0;i<6;i++)
		denergy_internal+=dstrain[i]*(deviatoric_stress_n[i]+deviatoric_stress[i]);

	energy_internal+=0.25*denergy_internal*(volume_current+volume_n)+mean_stress*dvolume;
	energy_internal_per_initial_volume=energy_internal/volume_0;
}

void MaterialModel::update_pressure_energy_EOS/*seleos*/(material* mat, bool& failure, const float& time_intervel, const float& grid_intervel)
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
		update_pressure_energy_EOS3(mat,failure,time_intervel, grid_intervel);
		break;
	}
}

void MaterialModel::update_pressure_energy_EOS1/*eos1*/(material* mat, bool& failure)
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

void MaterialModel::update_pressure_energy_EOS2/*eos2*/(material* mat, bool& failure, const float& time_intervel, const float& grid_intervel)
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

void MaterialModel::update_pressure_energy_EOS3/*eos3*/(material* mat, bool& failure, const float& time_intervel, const float& grid_intervel)
{
	//cout<<"step5.2.1"<<endl;
	//cout<<mean_stress<<" "<<energy_internal<<" "<<sound_speed<<endl;
	float er1v=exp(-mat->cEOS[2]*relative_volume);
	float er2v=exp(-mat->cEOS[3]*relative_volume);
	float wdr1v=mat->cEOS[0]-mat->cEOS[0]*mat->cEOS[4]/(mat->cEOS[2]*relative_volume);
	float wdr2v=mat->cEOS[1]-mat->cEOS[1]*mat->cEOS[4]/(mat->cEOS[3]*relative_volume);

	//calc_bulk_viscosity(mat, time_intervel, grid_intervel);
	float A=wdr1v*er1v+wdr2v*er2v;//+bulk_viscosity_force;
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

void MaterialModel::calc_bulk_viscosity/*bulkq*/(material* mat, const float& time_intervel, const float& grid_intervel)
{
	float bulk_strain_rate=(dstrain[0]+dstrain[1]+dstrain[2])/time_intervel;
	bulk_strain_rate=min(bulk_strain_rate,0.0f);
	bulk_viscosity_force=density_current*grid_intervel*grid_intervel*mat->bq1*bulk_strain_rate*bulk_strain_rate-mat->bq2*density_current*grid_intervel*sound_speed*bulk_strain_rate;
}

void MaterialModel::null_model_compute()/*m3md7*/
{
	for(int i=0;i<6;i++)
		deviatoric_stress[i]=0;
}

void MaterialModel::simplified_johnson_cook_model_compute(material* mat,const float& time_intervel)/*m3md5*/
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