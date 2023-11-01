
/******************************************************************************
ricker_wavelet -- Compute the   time response of a source function as
    a Ricker wavelet with peak frequency "fpeak" Hz.    
*******************************************************************************
Input: 
    int nt      number samples in output wavelet
    float dt    time step
    float fpeak peak frequency of the Ricker wavelet
    
*******************************************************************************
Output: 50
    float wavelet   array[nt] of computed wavelet

*******************************************************************************
Author: Tong Fei, 1993, Colorado School of Mines.
******************************************************************************/
/*void ricker1_wavelet (int nt, double dt, double fpeak, double *wavelet)
{
    int it;
    double   t1, t0; 

    t0=1.0/fpeak;

    for (it=0; it<nt; it++) {
        t1=it*dt;
        wavelet[it] = exp(-PI*PI*fpeak*fpeak*(t1-t0)*(t1-t0))*
            (1.0-2.*PI*PI*fpeak*fpeak*(t1-t0)*(t1-t0));
    }
}*/


/*void ricker1_wavelet (int nt, double dt, double fpeak, double *wavelet)
{
    int it;
    double   t1, t0; 

    //t0=1.0/fpeak;

    for (it=0; it<nt; it++) {
        t1=it*dt;
        wavelet[it] = exp(-4*PI*PI*fpeak*fpeak*t1*t1*0.25*0.25)*
            cos(2*PI*fpeak*t1);
    }
}*/

double * ricker(double Freq,double dt,int *Npoint)
{
	int i; /* they are the dummy counter*/

        double Bpar,t,u,*Amp;
	int Np1,N;
	
	if(Freq==0.0)Freq=30.0;
	if(dt==0.0)dt=0.004;
	Bpar=sqrt(6.0)/(PI*Freq);
	N=ceil(1.35*Bpar/dt);
        //printf("N=%d\n",N);
	Np1=N;
	*Npoint=2*N+1;
	 
	Amp=alloc1double(*Npoint);
	
	Amp[Np1]=1.0;
  
	for(i=1;i<=N;i++) {
		t=dt*(double)i;
		u=2.0*sqrt(6.0)*t/Bpar;
		Amp[Np1+i]=Amp[Np1-i]=0.5*(2.0-u*u)*exp(-u*u/4.0);
	}

	return Amp;

}

