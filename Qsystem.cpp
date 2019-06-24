#define Qsystem_cpp

#include "Qsystem.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>


using namespace std;

bool verbal = true;	//cout on/off
bool file = true;	//

double slen;




//--------------------------------------------------------
// Constructor
//--------------------------------------------------------
Qsystem::Qsystem(int i_nw, int i_ncrit, int i_ns, int i_nb, int i_nbSkip, double i_tau): 
nw(i_nw),
ncrit(i_ncrit),
ns(i_ns),
nb(i_nb),
nbSkip(i_nbSkip),
tau(i_tau)
{
	walkers = new walker [(int)(ncrit*2)];							
	 np=walkers[0].GetNp();
	 nc=walkers[0].GetNc();
	

	nd=1;
	isoStep=0.5;

	maxStep = new double [nd];	// max step in Metropolis algorithm  

	walker tmpp;
	tmp=tmpp;
	tmpp.CleanUp();
	for(int k = 0; k < nd; k++)
		maxStep[k]=isoStep;
	for(int iw = 0; iw< nw; iw++)
	{
		walker wlk;
		walkers[iw] = wlk;
		wlk.CleanUp();
	}
	

	ifstream fin("runNo");
	fin>>runNo;
	fin.close();
	ofstream fout("runNo");
	fout<<runNo+1;
	fout.close();


		
	stringstream ss;

	ss <<"../../out/"<<nc<<"c"<<(double)np/nc<<"p/";
	folder = ss.str();
	cout<<folder<<endl;

}//-------------------end constructor---------

//--------------------destructor--------------
Qsystem::~Qsystem()
{
	delete [] maxStep;
	
	for(int iw=0; iw<nw; iw++)
		walkers[iw].CleanUp();
	
	delete [] walkers;
	
}
//--------------------end destructor----------


//Creation of output files 
//--------------------------------------------------------

string Qsystem::FileName(string name)
{

	string extension(".dat");
	string dat = folder;
	dat+=name;

	
	stringstream ss;

	name = "_";
	dat += name;
	ss <<walkers[0].GetA(0)<<name<<ncrit<<name<<tau;
	dat += ss.str();
	ss.str("");
	
	
	dat += extension;

	return dat;
}//-------------------------------------------------------



//--------------------------------------------------------------------------------------------
//-------------------VMC LOOP-----------------------------------------------------------------
//--------------------------------------------------------------------------------------------
void Qsystem::VMC()
{

	double SwE, SsE, SbE, SbE2;	
	double accept = 0.0, acc;

	double **Psi = new double *[nc];
	for(int ic=0; ic<nc; ic++)
	{
		Psi[ic]=new double[200];
		for(int i=0; i<200; i++)
			Psi[ic][i]=0;
	}
	double psi;
	int ib,is,iw,k,i;



	ofstream fout(FileName("E").c_str());
	ofstream pout(FileName("psiVMC").c_str());
	pout << showpoint << scientific << setprecision(8);
	fout << showpoint << scientific << setprecision(8);
	cout << showpoint << scientific << setprecision(2);
	

	
	SbE = 0.0; 
	SbE2 = 0.0;
	

	
	for(ib = 0; ib < nb; ib++)
	{

		SsE=0.0;
		for(is = 0; is < ns; is++)
		{ 
			SwE=0.0;
			for(iw = 0; iw < nw; iw++)
			{

				Metropolis(iw, &accept);
				SwE+=walkers[iw].GetEloc();
			}//walkers

			if(ib>nbSkip) //accumulation 
			{
				SsE+=SwE/nw;
			}
		}//steps
		acc=accept/((ib+1)*nw*ns*np);
		AdjustStep(acc);
		if(ib>nbSkip)  
		{
			SbE+=SsE/ns;
			SbE2+=SsE*SsE/ns/ns;
			fout 
				<< setw(7)  << ib-nbSkip
				<< setw(16) << SsE/ns
				<< setw(16) << SbE/(ib-nbSkip) << endl;

			
			for(iw=0; iw<nw; iw++)
			{
				int count=0;
				for(int ic=0; ic<nc; ic++)
				{
					
					int npc=walkers[iw].GetComponent(ic);
					for(int ip=0; ip<npc; ip++)
					{
						psi = walkers[iw].GetPosition(ip+count,0);
						if(psi<5 &&psi>-5)
						{
							i = int((psi+5)/10.0*200);
							Psi[ic][i]+=1.0/nw/npc;
						
						}
					}
					count+=npc;
				}
					
			}
				
			
		}

		if(verbal)
			cout << setw(6) << ib-nbSkip << ". blok:  " << int(round(acc*100.))
			<< "% prihvacenih,  Eb = "<< setw(10) << SsE/ns << endl;
	}//blocks
	AvgE = SbE/(nb-nbSkip);
	sigmaE = sqrt((SbE2/(nb-nbSkip)-AvgE*AvgE)/(nb-nbSkip-1.));
	accept /= (nw*ns*nb*np);

	fout << "#konacni maksimalni pomak: " ;
	for(k=0;k<nd;k++)
		fout<< setw(10) << maxStep[0];
		
	fout << endl << endl;
	fout << "#E = " << AvgE << " +- " << sigmaE << endl <<
		"#gamma = "<<g<< endl;
	
	fout.close();
	
	for(int ic=0; ic<nc; ic++)
	{
		for(i=0;i<200;i++)
			pout<<i<<setw(30)<<Psi[ic][i]/(nb-nbSkip)<<endl;
		cout<<endl;
	}

	pout.close();
	
	
}
//--------------------------------------------------------------------------------------------
//--------------------VMC LOOP END------------------------------------------------------------
//--------------------------------------------------------------------------------------------



//Metropolis Algorithm
//--------------------------------------------------------
void Qsystem::Metropolis(int iw, double *accept)
{
	walker wlk;
	double p, P;
	double w;
	double dx, xt;




	for(int ip = 0; ip < np; ip++)
	{
		wlk=walkers[iw];
		wlk.Update();

		P=wlk.GetP();

		for(int k = 0; k < nd; k++)	//coordinates
		{
			dx = 2. * (ran1(&idum) - 0.5) * maxStep[k];
			xt = wlk.GetPosition(ip, k) + dx;
			wlk.SetPosition(ip, k, xt);
		}//end coordinates
	
		wlk.Update();
		p=wlk.GetP();
			
		w=p/P;
		
		if(w>=1)
		{
			for(int k=0;k<nd;k++)
				walkers[iw].SetPosition(ip,k,wlk.GetPosition(ip,k));
			*accept+=1; 
		}
		else if(ran1(&idum)<=w)
		{
			for(int k=0;k<nd;k++)
				walkers[iw].SetPosition(ip,k,wlk.GetPosition(ip,k));
			*accept+=1;
		}

	}
	walkers[iw].Update();	
	wlk.CleanUp();
}//Metropolis end


//--------------------------------------------------------------------------------------------
//-------------------DMC LOOP-----------------------------------------------------------------
//--------------------------------------------------------------------------------------------
void Qsystem::DMC()
{
	double SwE, SsE, SbE, SbE2;	

	int ib,is,iw,i;
	
	int nws; //total number of walkers*steps in current block
	
	nwmin=(int)(ncrit*0.9);
	nwmax=(int)(ncrit*1.1);
	reduce = 0.5*((double)(nwmin+nwmax)/(double)nwmax);
	amplify = 0.5*((double)(nwmin+nwmax)/(double)nwmin);
	
	double **Psi = new double *[nc];
	for(int ic=0; ic<nc; ic++)
	{
		Psi[ic]=new double[200];
		for(int i=0; i<200; i++)
			Psi[ic][i]=0;
	}
	double psi;
	
	double E,Eold;
	
	
	ofstream fout(FileName("Edmc").c_str());
	ofstream pout(FileName("psiDMC").c_str());
	fout << showpoint << scientific << setprecision(8);
	pout << showpoint << scientific << setprecision(8);
	
	cout << showpoint << scientific << setprecision(2);
	
	fout	<<"#################################################"<<endl
	<<"#"<<setw(5)<<"nw = "<<setw(5)<<nw<<endl
	<<"#"<<setw(5)<<"mcs = "<<setw(5)<<ns*(nb-nbSkip)<<endl
	<<"#"<<setw(5)<<"ns = "<<setw(5)<<ns<<endl
	<<"#"<<setw(5)<<"nb = "<<setw(5)<<nb<<endl
	<<"#"<<setw(5)<<"nbSkip = "<<setw(5)<<endl
	<<"#"<<setw(5)<<"tau = "<<setw(5)<<tau<<endl
	<<"############################################"<<endl
	<<"#	ib"<<setw(15)<<"<Eb>"<<setw(15)<<"<E>"<<setw(15)<<"nw"<<endl;	
	
	
	
	Et=AvgE;
	cout<<Et<<endl;
	
	SbE = 0.0; 
	SbE2 = 0.0;

	nbSkip=-1; //change if VMC is not used as a preparation step
	
	for(ib = 0; ib < nb; ib++)
	{ 
		SsE=0.0; 
		nws=0;
		for(is = 0; is < ns; is++)
		{	
			nwNew = nw;
			nwDead = 0;
			SwE=0.0;
			for(iw = 0; iw < nw; iw++)
			{
				if(walkers[iw].IsAlive())
				{  
					Eold=walkers[iw].GetEloc();
					tmp = walkers[iw];
					Diffusion();
					Drift(); 
					E=tmp.GetEloc();
					nSons = Branch(&iw,Eold,E); 

					SwE+=nSons*E;
					nwNew += (nSons-1);	
					if(nSons==0)
						nwDead++;

				}//if walker alive

			}//walkers
			
			BuryDeadWalkers(); 
			nw=nwNew;
			
			if(ib>nbSkip) //accumulation 
			{
				nws+=nw;
				SsE+=SwE;
			}
			else
				Et=SwE/nw;
			nw+=nwDead;


		}//steps

		if(ib>nbSkip)  
		{
			SbE+=SsE/nws;	
			SbE2+=SsE*SsE/nws/nws;
			
			AvgE = SbE/(ib-nbSkip);
			
			Et=AvgE;
			
			if(file)
			{
				fout 
					<< setw(7)  << ib-nbSkip
					<< setw(16) << SsE/nws
					<< setw(16) << SbE/(ib-nbSkip) 
					<< setw(16) << nw
					<< endl;
			}
			
			for(iw=0; iw<nw; iw++)
			{
				int count=0;
				for(int ic=0; ic<nc; ic++)
				{
					
					int npc=walkers[iw].GetComponent(ic);
					for(int ip=0; ip<npc; ip++)
					{
						psi = walkers[iw].GetPosition(ip+count,0);
						if(psi<5 &&psi>-5)
						{
							i = int((psi+5)/10.0*200);
							Psi[ic][i]+=1.0/nw/npc;
						
						}
					}
					count+=npc;
				}
					
			}
			
		}

		if(verbal)
			cout << setw(6) << ib-nbSkip << ". blok:  Eb = "<< setw(10) << SsE/nws << setw(30)<<nw<<endl;
	}//blocks
	AvgE = SbE/(nb-1-nbSkip);
	sigmaE = sqrt((SbE2/(nb-nbSkip-1)-AvgE*AvgE)/(nb-1-nbSkip));
	
	if(file)
	{
		fout << endl << endl;
		fout << "#E = " << AvgE << " +- " << sigmaE << endl << endl;
		fout.close();

		for(int ic=0; ic<nc; ic++)
		{
			for(i=0;i<200;i++)
				pout<<i<<setw(30)<<Psi[ic][i]/(nb-nbSkip)<<endl;
			cout<<endl;
		}
	
		pout.close();
	}
	
	
	
}
//--------------------------------------------------------------------------------------------
//--------------------DMC LOOP END------------------------------------------------------------
//--------------------------------------------------------------------------------------------

void Qsystem::Diffusion()
{
	double x,dx;
	double D;
	for(int ip=0; ip<np; ip++)
	{
		D=tmp.GetD(ip); 
		for(int k = 0; k < nd; k++)
		{
			dx=var(D)*gasdev(&idum);
			x=tmp.GetPosition(ip,k);
			tmp.SetPosition(ip,k,x+dx);
		}
	}
}

double Qsystem::Drift()
{
	double **F1, **F2, x, **xold, dx;
	double Eloc;
	double * D = new double [np];
	
	F1=new double * [np];
	F2=new double * [np];
	xold=new double * [np];

	
	for(int ip=0; ip<np; ip++)
	{
		F1[ip]=new double[nd];
		F2[ip]=new double[nd];
		xold[ip]=new double[nd];
		D[ip]=tmp.GetD(ip);

		for(int k=0; k<nd; k++)
		{
			F1[ip][k]=tmp.GetF(ip,k);	//Calculating the drift force
			x=tmp.GetPosition(ip,k);
			dx=D[ip]*0.5*tau*F1[ip][k];
			xold[ip][k]=x;
			tmp.SetPosition(ip,k,x+dx);	//temporary drift move
		}
	}

	tmp.Update();
	for(int ip=0; ip<np; ip++)
		for(int k=0; k<nd; k++)
			F2[ip][k]=tmp.GetF(ip,k);	//Calculating the drift force
	
	for(int ip=0; ip<np; ip++)
	{
		for(int k=0; k<nd; k++)
		{
			dx=D[ip]*0.25*tau*(F1[ip][k]+F2[ip][k]);
			x=xold[ip][k];				//return to the previous position
			tmp.SetPosition(ip,k,x+dx);		//srednji driftni pomak
		}
	}
	tmp.Update();
	for(int ip=0; ip<np; ip++)
		for(int k=0; k<nd; k++)
			F1[ip][k]=tmp.GetF(ip,k);	//Calculating the drift force
	Eloc=tmp.GetEloc();				//and the local energy
	
	for(int ip=0; ip<np; ip++)
	{
		for(int k=0; k<nd; k++)
		{
			dx=D[ip]*tau*F1[ip][k];
			x=xold[ip][k];
			tmp.SetPosition(ip,k,x+dx);		//final drift move
		}
	}
	for(int ip=0;ip<np;ip++)
	{
		delete [] F1[ip];
		delete [] F2[ip];
		delete [] xold[ip];
	}
	delete [] F1;
	delete [] F2;
	delete [] D;
	delete [] xold;
	return Eloc;
}

int Qsystem::Branch(int * iw, double Eold, double E)
{
	w=exp(tau * (Et -0.5*(E+Eold)));

	double r = ran1(&idum);
	nSons = (int)(w+r); 
	
	if(nSons>0)
	{	
		walkers[*iw]=tmp;
		if(nw>nwmax)
			nSons = (int)(nSons*reduce+ran1(&idum));
		else if(nw<nwmin)
			nSons = (int)(nSons*amplify+ran1(&idum));
		
		if (nSons>1)
		{
			CopyWalker(*iw,nSons); 
			return nSons;
		} 
		if(nSons==1)
			return nSons;

	} 
	nSons=0;
	walkers[*iw].Die();
	return nSons;
}

//--------------------------------------------------------
void Qsystem::AdjustStep(double acc) 
{
	if(acc>0.5) for(int k=0; k<nd; k++) {maxStep[k] *= 1.05;}
	if(acc<0.5) for(int k=0; k<nd; k++) {maxStep[k] *= 0.95;}
}
//---------------------------------------------------------
void Qsystem::SaveConfig(string name)
{
	ofstream out(name.c_str());
	out<<AvgE<<setw(30)<<sigmaE<<endl<<endl;
	for(int iw=0; iw<nw; iw++)
	{
		for(int ip=0; ip<np; ip++)
		{
			for(int k=0; k<nd; k++)
				out<<walkers[iw].GetPosition(ip,k)<<setw(30);
			out<<endl;
		}
		out<<endl;
	}
	out.close();
}
//--------------------------------------------------------
void Qsystem::InitConfig(string name)
{
	ifstream in(name.c_str());
	double x;
	in>>AvgE>>sigmaE;
	for(int iw=0; iw<nw; iw++)
		for(int ip=0; ip<np; ip++)
			for(int k=0; k<nd; k++)
			{
				in>>x;
				walkers[iw].SetPosition(ip,k,x);
			}	
}
//--------------------------------------------------------
void Qsystem::CopyWalker(int iw, int n)
{
	walker copy;
	copy = walkers[iw];
	
	int N=nwNew+nwDead; 
	

	
	for(int i=0;i<n-1;i++)
	{
		walkers[N+i]=copy; 
	}
	copy.CleanUp();
}
//--------------------------------------------------------
void Qsystem::BuryDeadWalkers()
{	
	walker copy;
	int N=nwNew+nwDead;				
	for(int iw=0; iw<N && nwDead>0; iw++)
	{
		if(!walkers[iw].IsAlive())
		{
			while(!walkers[N-1].IsAlive())
			{
				nwDead--;
				N--;
			}
			if(iw<N-1)
			{	
				copy = walkers[N-1]; 
				walkers[iw] = copy;
				N--;
				nwDead--;
			}
		}
	}	
	
	copy.CleanUp();
}



void Qsystem::SetScatteringLength(double A)
{
	for(int iw=0; iw<nw; iw++)
		for(int ip=0; ip<np; ip++)
			walkers[iw].SetA(A,ip);
}
