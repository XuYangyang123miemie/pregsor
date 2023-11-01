
#include "stdlib.h"

/*compute the limit of integration*/

double gaus2_CW(x,z,fxr,fzr,h,vi,vj,kk,n,js,ss,f)
int vi,vj,n,js[];
//double x,z,fxr,fzr,kk,h;
double x,z,fxr,fzr,h;
dcomplex kk; //kk:double->complex 20201112 13:05
double (*f)();
void (*ss)();
{
	int m,j,k,q,l,*is;
	double y[2],p,s,*x1,*a,*b;
	static double t[5]={-0.9061798459,-0.5384693101,0.0,0.5384693101,0.90617984459};
	static double c[5]={0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851};
	is=malloc(2*(n+1)*sizeof(int));
	x1=malloc(n*sizeof(double));
	a=malloc(2*(n+1)*sizeof(double));
	b=malloc((n+1)*sizeof(double));
	m=1;
	l=1;
	a[n]=1.0;
	a[2*n+1]=1.0;
	while(l==1)
	{
		for(j=m;j<=n;j++)
		{
			(*ss)(fxr,fzr,h,vi,vj,j-1,n,y);
			a[j-1]=0.5*(y[1]-y[0])/js[j-1];
			b[j-1]=a[j-1]+y[0];
			x1[j-1]=a[j-1]*t[0]+b[j-1];
			a[n+j]=0.0;
			is[j-1]=1;
			is[n+j]=1;
			//printf("j=%d\n",j);
			//printf("x1[%d]=%e\n",j-1,x1[j-1]);
		}
		j=n;
		q=1;
		while(q==1)
		{
			k=is[j-1];	
			//printf("k=%d\n",k);
			if(j==n) p=(*f)(n,x,z,kk,x1);
			else p=1.0;
			a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];
			//printf("a[%d]=%e ",n+j,a[n+j]);
			is[j-1]=is[j-1]+1;
			if(is[j-1]>5)
			{
		        if(is[n+j]>=js[j-1])
		        {
			         j=j-1;
					 q=1;
					 if(j==0)
					 {
						 s=a[n+1]*a[0];
						 //printf("555555555555\n");
						 //printf("s=%e\n",s);
						 free(is);
						 free(x1);
						 free(a);
						 free(b);
						 return(s);
					 }
		        }else
				{
					is[n+j]=is[n+j]+1;
					b[j-1]=b[j-1]+a[j-1]*2.0;
					is[j-1]=1;
					k=is[j-1];
					x1[j-1]=a[j-1]*t[k-1]+b[j-1];
					if(j==n) q=1;
					else q=0;
					//printf("33333333 ");
					//printf("x1[%d]=%e\n",j-1,x1[j-1]);
					//printf("next subinterval\n");
				}
			}else
			{
				k=is[j-1];
				x1[j-1]=a[j-1]*t[k-1]+b[j-1];
				if(j==n) q=1;
				else q=0;
				//printf("222222 ");
				//printf("x1[%d]=%e\n",j-1,x1[j-1]);
				//printf("go on next node\n");
				
			}
		}
		m=j+1;
	}
	
    return(s);	
}

double gaus4_CW(fxr,fzr,h,vi,vj,vi1,vj1,kk,n,js,ss,f)
int vi,vj,vi1,vj1,n,js[];
//double fxr,fzr,kk,h;
double fxr,fzr,h;
dcomplex kk;
double (*f)();
void (*ss)();
{
	int m,j,k,q,l,*is;
	double y[2],p,s,*x1,*a,*b;
	static double t[5]={-0.9061798459,-0.5384693101,0.0,0.5384693101,0.90617984459};
	static double c[5]={0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851};
	is=malloc(2*(n+1)*sizeof(int));
	x1=malloc(n*sizeof(double));
	a=malloc(2*(n+1)*sizeof(double));
	b=malloc((n+1)*sizeof(double));
	m=1;
	l=1;
	a[n]=1.0;
	a[2*n+1]=1.0;
	while(l==1)
	{
		for(j=m;j<=n;j++)
		{
			(*ss)(fxr,fzr,h,vi,vj,vi1,vj1,j-1,n,y);
			a[j-1]=0.5*(y[1]-y[0])/js[j-1];
			b[j-1]=a[j-1]+y[0];
			x1[j-1]=a[j-1]*t[0]+b[j-1];
			a[n+j]=0.0;
			is[j-1]=1;
			is[n+j]=1;
			//printf("j=%d\n",j);
			//printf("x1[%d]=%e\n",j-1,x1[j-1]);
		}
		j=n;
		q=1;
		while(q==1)
		{
			k=is[j-1];	
			//printf("k=%d\n",k);
			if(j==n) p=(*f)(n,kk,x1);
			else p=1.0;
			a[n+j]=a[n+j+1]*a[j]*p*c[k-1]+a[n+j];
			//printf("a[%d]=%e ",n+j,a[n+j]);
			is[j-1]=is[j-1]+1;
			if(is[j-1]>5)
			{
		        if(is[n+j]>=js[j-1])
		        {
			         j=j-1;
					 q=1;
					 if(j==0)
					 {
						 s=a[n+1]*a[0];
						 //printf("555555555555\n");
						 //printf("s=%e\n",s);
						 free(is);
						 free(x1);
						 free(a);
						 free(b);
						 return(s);
					 }
		        }else
				{
					is[n+j]=is[n+j]+1;
					b[j-1]=b[j-1]+a[j-1]*2.0;
					is[j-1]=1;
					k=is[j-1];
					x1[j-1]=a[j-1]*t[k-1]+b[j-1];
					if(j==n) q=1;
					else q=0;
					//printf("33333333 ");
					//printf("x1[%d]=%e\n",j-1,x1[j-1]);
					//printf("next subinterval\n");
				}
			}else
			{
				k=is[j-1];
				x1[j-1]=a[j-1]*t[k-1]+b[j-1];
				if(j==n) q=1;
				else q=0;
				//printf("222222 ");
				//printf("x1[%d]=%e\n",j-1,x1[j-1]);
				//printf("go on next node\n");
				
			}
		}
		m=j+1;
	}
	
	return(s);
}

/* integrand function*/


void ss2_CW(double fxr,double fzr,double h,int vi,int vj,int j,int n, double lim[])
{
    n=n;
	switch(j){
		case 0:{lim[0]=fxr+vi*h;lim[1]=lim[0]+h; break;}
		case 1:{lim[0]=fzr+vj*h;lim[1]=lim[0]+h; break;}
		
		default:{}
	}
	return;
	
}

void ss4_CW(double fxr,double fzr,double h,int vi,int vj,int vi1,int vj1,int j,int n, double lim[])
{
    n=n;
	switch(j){
		case 0:{lim[0]=fxr+vi*h;lim[1]=lim[0]+h; break;}
		case 1:{lim[0]=fzr+vj*h;lim[1]=lim[0]+h; break;}
		case 2:{lim[0]=fxr+vi1*h;lim[1]=lim[0]+h;; break;}
		case 3:{lim[0]=fzr+vj1*h;lim[1]=lim[0]+h; break;}
		default:{}
	}
	return;
	
}




/*Gauss numerical integration*/


//double fre2(int n,double x,double z,double kk,double var[])
double fre2_CW(int n,double x,double z,dcomplex kk,double var[])
{
	double r;
	double func;
	n=n;
        r=sqrt((var[0]-x)*(var[0]-x)+(var[1]-z)*(var[1]-z));
	//r=kk*r;
	if(r==0.0)
	{
           func=2*PI*0.00001*0.00001*GreenCWr(dcrmul(kk,0.00001));//20201207 15:00 (kk.i*0.00001)--> (-kk.i*0.00001)
        }else
        {
           func=GreenCWr(dcrmul(kk,r));//20201207 15:00 (kk.i*r)-->(-kk.i*r)
        }
	
	return (func);
}

//double fim2(int n,double x,double z,double kk,double var[])
double fim2_CW(int n,double x,double z,dcomplex kk,double var[])
{
	double r;
	double func;
	n=n;
        r=sqrt((var[0]-x)*(var[0]-x)+(var[1]-z)*(var[1]-z));
	//r=kk*r;
	if(r==0.0)
	{
           func=2*PI*0.00001*0.00001*GreenCWi(dcrmul(kk,0.00001));
        }else
        {
           func=GreenCWi(dcrmul(kk,r));
        }
	
	return (func);
}


//double fre4(int n,double kk,double var[])
double fre4_CW(int n,dcomplex kk,double var[])
{
	double r;
	double func;
	n=n;
    
        r=sqrt((var[0]-var[2])*(var[0]-var[2])+(var[1]-var[3])*(var[1]-var[3]));
	//r=kk*r;
        if(r==0.0)
	{
           func=2*PI*0.00001*0.00001*GreenCWr(dcrmul(kk,0.00001));
        }else
        {
           func=GreenCWr(dcrmul(kk,r));
        }
	
	return (func);
}

//double fim4(int n,double kk,double var[])
double fim4_CW(int n,dcomplex kk,double var[])
{
	double r;
	double func;
	n=n;
        r=sqrt((var[0]-var[2])*(var[0]-var[2])+(var[1]-var[3])*(var[1]-var[3]));
	//r=kk*r;
	if(r==0.0)
	{
           func=2*PI*0.00001*0.00001*GreenCWi(dcrmul(kk,0.00001));
        }else
        {
           func=GreenCWi(dcrmul(kk,r));
        }
	
	return (func);
}



