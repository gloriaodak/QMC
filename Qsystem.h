#ifndef Qsystem_hpp

#include <iostream>  
#include <new>      
#include <cmath>  
#include<stdlib.h>
#include <vector>  
#include "walker.h"

using namespace std; 

class Qsystem{
public:
	Qsystem(int, int, int, int, int, double);
	~Qsystem();


	void VMC();
	
	void DMC();
	void Diffusion();
	double Drift();
	int Branch(int *,double, double);
	

	void Metropolis(int,double*);

	void AdjustStep(double);

	void SetScatteringLength(double);
	
	void SetEt(double et) {Et=et;};
	double GetEt(){return Et;};
	
	int runNo;
	string FileName(string);

	
	double var(double d) {return sqrt(2.0*d*tau);};
	
	void SaveConfig(string);
	void InitConfig(string);
	
	void CopyWalker(int,int);
	void BuryDeadWalkers();
	
	double GetE() {return AvgE;};
	double GetSigmaE() {return sigmaE;};	

	int GetNp() {return np;};
	int GetNc() {return nc;};


private:
	int np; // number of particles
	int nw; // number of walkers
	int ncrit; //critical number of walkers for DMC branching
	int ns; // number of steps 
	int nb; // number of blocks
	int nc; //number of components
	int nbSkip; // number of blocks for equilibration
	int nd; // number of dimensions

	
	int nwmin, nwmax;
	double reduce, amplify;
	int nSons;
	
	int nwNew, nwDead;
	
	double w;
	
	double Et; //referent energy
	
	double AvgE, sigmaE;

	double min, max; 

	double * maxStep;	// max step in Metropolis algorithm  
	double isoStep;		// for isotropic systems
	
	double tau;	//timestep for DMC

	string folder; // for saving data

	walker * walkers;
	walker tmp;
	//vector<walker> walkers;
};

#endif
