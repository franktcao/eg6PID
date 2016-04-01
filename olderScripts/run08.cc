#include "EG5_test.C"

void run()
{
	EG5_test t;
	//cout << " Class Initialized " << endl;
	t.Begin();
	//t.printNumEvents();
	t.Loop();
	t.drawHist();
	t.writeFiles();
//	return nullptr;
}
