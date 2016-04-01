#include <iostream>
using namespace std;

class Particle {
	static int eventN, charge;
  public:
    Particle (int,int);
};

Particle::Particle(int eventNumber, int charge0) {
  eventN = eventNumber;
  charge = charge0;
}

int main () {
	Particle part = Particle(1,2);
	cout << "rectb area: " << part.eventN << ", " << part.charge << endl;
	return 0;
}
