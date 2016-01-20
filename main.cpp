#include <vector>
#include "FF_simulator.h"

using namespace std;

void main()
{
	simulator* s=new simulator();
	s->initialize();
	s->update();

	vector<particle*>& particles = s->p_particle_data->particle_list;
	int nParticles = particles.size();
	for (int i = 0; i < nParticles; i++) {
		particle* p = particles[i];
		cout<<i<<" "<<p->position_t.x<<" "<<p->position_t.y<<" "<<p->position_t.z<<endl;
	}
}