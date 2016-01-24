#pragma once
#ifndef FF_PARTICLE_HH
#define FF_PARTICLE_HH
#include <vector>
#include <cinder/Vector.h>
using namespace std;
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
	cinder::Vec3f deviatoric_stress_plane;//sdyz xz xy
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

	void debug(int particle_id=0)
	{
		if(this->particle_id==particle_id)
		{
			cout<<"volume:"<<volume<<" "<<"mass:"<<particle_mass<<" "<<"density:"<<particle_mass/volume<<endl;
			cout<<"fail:"<<failure<<" skip:"<<skipthis<<" icell:"<<icell<<endl;
			cout<<"LT:"<<lighting_time<<" damage:"<<damage<<" celsius_t"<<celsius_temperature<<endl;
			cout<<"yield:"<<yield_stress<<" mean:"<<mean_stress<<" mises:"<<mises_stress<<endl;
			cout<<"sound:"<<sound_speed<<" epeff:"<<effective_plastic_strain<<" inter:"<<internel_energy<<endl;
			cout<<"velocity:"<<velocity.x<<" "<<velocity.y<<" "<<velocity.z<<endl;
			cout<<"position:"<<position_t.x<<" "<<position_t.y<<" "<<position_t.z<<endl;
			cout<<"load:"<<force_load.x<<" "<<force_load.y<<" "<<force_load.z<<endl;
			cout<<"sd:"<<deviatoric_stress_asix.x<<" "<<deviatoric_stress_asix.y<<" "<<deviatoric_stress_asix.z<<" "
				<<deviatoric_stress_plane.x<<" "<<deviatoric_stress_plane.y<<" "<<deviatoric_stress_plane.z<<endl;
		}
	}
};
#endif