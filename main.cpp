#include "Qsystem.cpp"
#include <iostream>



using namespace std;

int main()
{
	int nw =100;					//starting number of walkers
	int ncrit = 500;				//critical number of walkers
	int ns = 200;					//number of steps
	int nb = 1000;					//number of blocks
	int nbSkip = 0;					//number of blocks to skip
	double tau = 0.001;				//DMC timestep

	bool VMC = true;
	bool DMC = true;

/*
	Qsystem * s = new Qsystem(nw,ncrit,ns,nb,nbSkip,tau);
		if(VMC)
		{
			s->VMC();
			s->SaveConfig("configVMC.dat");
		}
		if(DMC)
		{
			s->InitConfig("configVMC.dat"); 
			s->DMC();
			s->SaveConfig("configDMC.dat");
		}
	
	delete s;
	

*/


string filename;
stringstream ss;
Qsystem * sustav = new Qsystem(nw,ncrit,ns,nb,nbSkip,tau); 

int nc=sustav->GetNc();
double np=sustav->GetNp()/(double)nc;

	ss <<"../../out/"<<nc<<"c"<<np<<"p/energy-scattering_length.dat";
	filename = ss.str();
	ss.str("");
	cout<<filename<<endl;
	ofstream fout(filename.c_str(), ofstream::out | ofstream::app);
	fout << showpoint << scientific << setprecision(8);	

	ss <<"../../out/"<<nc<<"c"<<np<<"p/configVMC.dat";
	filename=ss.str();

sustav->VMC();
sustav->SaveConfig(filename.c_str());
	
for(double a = -0.1; a >= -4; a-=.1)
{
	sustav->InitConfig(filename.c_str());
	sustav->SetScatteringLength(a);
	sustav->DMC(); 

	fout<<a<<setw(30)<<sustav->GetE()<<setw(30)<<sustav->GetSigmaE()<<endl;
}


	/*
	ofstream fout("../DMCdata_fermi-fermi/2+1/walkers/walkers.dat", ofstream::out | ofstream::app);
	

	

	for(ncrit=400; ncrit<=5000; ncrit+=200)
	{
		Qsystem * sustav = new Qsystem(nw,ncrit,ns,nb,nbSkip,tau); 
		sustav->InitConfig("configVMC.dat");
		sustav->DMC(); 
		fout<<ncrit<<setw(30)<<sustav->GetE()<<setw(30)<<sustav->GetSigmaE()<<endl;
		delete sustav;
	}
	*/
	
/*
	
	ofstream fout("../DMCdata_fermi-fermi/2+1/timestep/timestep.dat", ofstream::out | ofstream::app);
	fout << showpoint << scientific << setprecision(8);	
	
	//nb=150;
	//Qsystem * sustav = new Qsystem(np,nw,ncrit,ns,nb,nbSkip,nd,min,max,isoStep,tau,folder);
	//sustav->VMC();
	//sustav->SaveConfig("configVMC.dat");
	//delete sustav;

	
	for(tau = 0.00025; tau < .0025; tau+=.00025)
	{
		Qsystem * sustav = new Qsystem(nw,ncrit,ns,nb,nbSkip,tau); 
		sustav->InitConfig("configVMC.dat");
		sustav->DMC(); 
		fout<<tau<<setw(30)<<sustav->GetE()<<setw(30)<<sustav->GetSigmaE()<<endl;
		delete sustav;
	}
	
*/
	return 0;
}
