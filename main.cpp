
#include "FF_simulator.h"
#include "FF_material.h"
#include "FF_particle.h"
#include "FF_grid.h"
#include <vector>

using namespace std;

void main()
{
	simulator s;
	s.initialize();
	int i=0;
	while(i<10)
	{
		s.update();
		i++;
	


		vector<particle*>& particles = s.p_particle_data->particle_list;
		int nParticles = particles.size();
		for (int i = 0; i < 5; i++) {
			particle* p = particles[i];
			cout<<i<<" "<<p->position_t.x<<" "<<p->position_t.y<<" "<<p->position_t.z<<endl;
		}
	}
}