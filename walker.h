#ifndef walker_h
 
#include <new>  
#include <fstream>
  
#include "ran1.c"
#include "gasdev.c"
#define EPSILON 0.000001

using namespace std; 

long idum = -236;

bool test=false;
//------------variational parameters------------------------
double g = 0.5;
//----------------------------------------------------------


class walker{
public:
	walker();
	//~walker(){CleanUp();};

	int InitData();
	
	void SetPosition(int ip, int k, double xk) { x[ip][k]=xk;}; 
	double GetPosition(int ip, int k) {return x[ip][k];};
	
	void UpdateEloc();
	void UpdatePsi();
	void UpdateR(int);
	
	double GetEloc() {return Eloc;};
	double GetPsi() {return Psi;};
	
	
	void UpdateP() {UpdatePsi(); P=Psi*Psi;};
	double GetP() {return P;};
	
	void UpdateF();
	double GetF(int ip, int k) {return F[ip][k];};
	
	double GetD(int ip) {return D[ip];};
	double GetHO(int ip) {return ho[ip];};
	double GetA(int ip) {return a[ip];};

	void SetA(double A, int ip) {a[ip]=A;};
	
	int GetSpin(int ip) {return spin[ip];};


	
	int GetNd() {return nd;};
	int GetNp() {return np;};
	int GetNc() {return nc;};

	int GetComponent(int ic) {return component[ic];};

	bool IsAlive() {return alive;};
	void Die() {alive=false;};
	
	void CleanUp();
	
	
	void Update() { UpdateP(); UpdateF(); UpdateEloc();}
	
	//overloads
	void operator = (walker wlk) 
	{
		CleanUp(); 
		nd = (int)wlk.GetNd();
		np = (int)wlk.GetNp();
		nc = (int)wlk.GetNc();
		D=new double [np];
		x=new double *[np];
		F=new double *[np];
		ho=new double [np];
		a=new double [np];
		spin=new int[np];
		component=new int[nc];
		
		for(int ic=0; ic<nc; ic++)
			component[ic]=wlk.GetComponent(ic);
		
		for(int ip=0;ip<np;ip++)
		{
			x[ip]=new double[nd];
			F[ip]=new double[nd];
			D[ip]=(double)wlk.GetD(ip);
			ho[ip]=wlk.GetHO(ip);
			a[ip]=wlk.GetA(ip);
			spin[ip]=wlk.GetSpin(ip);

			for(int k=0;k<nd;k++)
			{
				x[ip][k]=(double)wlk.GetPosition(ip,k);
				F[ip][k]=(double)wlk.GetF(ip,k);	
			}

		} 
		Eloc=(double)wlk.GetEloc();
		alive = (double)wlk.IsAlive();
	}
	


private:
	int nd; //number of space dimensions
	
	int np; //number of partices

	int nc; //number of components
	

	double xMin, xMax;
	
	double ** x; //np nd-dimensional position vectors
	
	int * component;	

	double Psi, Eloc;
	
	double P; //probability = Psi*Psi
	
	double ** F; //drift force
	
	double * ho; //=m*omega**2 / 2
	
	double * a; //scattering length
	
	int * spin;
	
	double * D; // =hbar**2 / 2m
	
	bool alive;

};
//---------------- default constructor -------------
walker::walker(): 
nd(1),
xMin(-1),
xMax(1),
alive(true)
{
	InitData();
	x = new double * [np];
	F = new double * [np];

	
	for(int ip = 0; ip < np; ip++)
	{
		x[ip] = new double [nd];
		F[ip] = new double [nd];
		for(int k = 0; k < nd; k++)
			x[ip][k] = ran1(&idum) * (xMax - xMin) + xMin;

	}
	
	Update();
}//end const----------------------------------------


//clean up unnecessary pointers---------------------
void walker::CleanUp()
{ 
	for(int ip=0; ip<np; ip++)
	{
		delete [] x[ip];
		delete [] F[ip];
	}
	delete [] D;
	delete [] x;
	delete [] F;
	delete [] ho;
	delete [] spin;
	delete [] a;
	delete [] component;


}
//--------------------------------------------------



//--------------------------------------------------------
//Wavefunction
//--------------------------------------------------------

void walker::UpdatePsi()				
{
	Psi=1.0;
	int n=0; 
	int m=0;

	for(int ic=0; ic<nc; ic++)
	{
		for(int ip=0; ip<component[ic]; ip++)
		{
			for(int jp=0; jp<ip; jp++)
				Psi *= fabs( x[ip+n][0] - x[jp+n][0] );


			m=0;
			for(int jc=0; jc<ic; jc++)
			{
				for(int jp=0; jp<component[jc];jp++)
				{
					Psi *= fabs( fabs( x[ip+n][0] - x[jp+m][0] ) - a[ip+n] );
				}
				m+=component[jc];
			}
			Psi *= exp( -g *  x[ip+n][0] * x[ip+n][0]);
		}
		n+=component[ic];
	}
	
}


//--------------------------------------------------------
//Local Energy
//--------------------------------------------------------
void walker::UpdateEloc()
{
	double e,r;
	Eloc=0;
	int n=0;
	int m=0;

	for(int ic=0; ic<nc; ic++)
	{
		for(int ip=0; ip<component[ic]; ip++)
		{	
			e=2*g;
			for(int jp=0; jp<ip; jp++)
			{
				r=x[ip+n][0]-x[jp+n][0];
				e+=1/r/r;
			}
			for(int jp=ip+1; jp<component[ic]; jp++)
			{
				r=x[ip+n][0]-x[jp+n][0];
				e+=1/r/r;
			}
			m=0;
			for(int jc=0; jc<ic; jc++)
			{
				for(int jp=0; jp<component[jc]; jp++)
				{
					r=fabs(x[ip+n][0]-x[jp+m][0]);
					e+=1/(r-a[ip+n])/(r-a[ip+n]);
				}
				m+=component[jc];
			}
			m+=component[ic];
			for(int jc=ic+1; jc<nc; jc++)
			{
				for(int jp=0; jp<component[jc]; jp++)
				{
					r=fabs(x[ip+n][0]-x[jp+m][0]);
					e+=1/(r-a[ip+n])/(r-a[ip+n]);
				}
				m+=component[jc];
			}
			
			e-=F[ip+n][0]*F[ip+n][0]/4.0;
			e*=D[ip+n];
			Eloc+=e;
			Eloc+=ho[ip+n]*x[ip+n][0]*x[ip+n][0];
		}
		n+=component[ic];
	}
}

//--------------------------------------------------------
//Drift Force
//--------------------------------------------------------
void walker::UpdateF() 
{
	int n=0;
	int m=0;

	for(int ic=0; ic<nc; ic++)
	{
		for(int ip=0; ip<component[ic]; ip++)
		{
			F[ip+n][0] = -2*g*x[ip+n][0];
			for(int jp=0; jp<ip; jp++)
				F[ip+n][0] += 1/( x[ip+n][0] - x[jp+n][0] );
			for(int jp=ip+1; jp<component[ic]; jp++)
				F[ip+n][0] += 1/( x[ip+n][0] - x[jp+n][0] );
			
			m=0;
			for(int jc=0; jc<ic; jc++)
			{
				for(int jp=0; jp<component[jc]; jp++)
				{
					double r = x[ip+n][0]-x[jp+m][0];
					//F[ip+n][0] += r / (fabs(r-a[ip+n])*fabs(r));
					F[ip+n][0] += (r-a[ip+n]*r/fabs(r))/(fabs(r)-a[ip+n])/(fabs(r)-a[ip+n]);
				}
				m+=component[jc];
			}
			m+=component[ic];
			for(int jc=ic+1; jc<nc; jc++)
			{
				for(int jp=0; jp<component[jc]; jp++)
				{
					double r = x[ip+n][0]-x[jp+m][0];
					//F[ip+n][0] += r / (fabs(r-a[ip+n])*fabs(r));
										F[ip+n][0] += (r-a[ip+n]*r/fabs(r))/(fabs(r)-a[ip+n])/(fabs(r)-a[ip+n]);

				}
				m+=component[jc];
			}	
			F[ip+n][0]*=2;
		}
		n+=component[ic];
	}
}
//--------------------------------------------------------
int walker::InitData()
{
	ifstream fin("IN");
	
	if (!fin.is_open()) 
		return 1;	

	fin>>np;
	fin>>nc;

	component=new int[nc];
	ho=new double[np];
	spin=new int[np];
	a=new double[np];
	D = new double [np];



	int n,ni,ic=0;
	double mass,w,s,slen;
	
	n=0;
	
	while(n<np)
	{
		fin.ignore(300, '\n');
		fin>>ni>>mass>>w>>s>>slen;
		component[ic]=ni;
		for(int i=0; i<ni; i++)
		{
			D[n+i] = 1.0 / (2.0 * mass);
			ho[n+i]=0.5*mass*w*w;
			spin[n+i]=s;
			a[n+i]=slen;
		}	
		n+=ni;
		ic++;
		//cout<<n<<endl;
	}



	
	fin.close();
	return 0;
}


#endif
