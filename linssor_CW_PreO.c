/*SSOR method LU=S -> (diag-G)U=S*/

void linssor_CW_PreO(complex **g,complex **diag,complex **M,complex **s,int nx,int nz,complex **u0,int itmax)
{
	int i,j;
	int iter;
	double norm,norm0,norm_S,err;
	complex alpha0,alpha;	
	complex **mutips,**r0,**Lr0,**mutipri,**u,**diagc,**r,**Lr,**mutipu0,**mutipr0,**mutipr;
	complex **u00,**r00,**rr;
	complex **utemp,**mutiputemp,**mutipr00;
	
	mutips=alloc2complex(nz,nx);
    r0=alloc2complex(nz,nx);
    Lr0=alloc2complex(nz,nx);
    mutipri=alloc2complex(nz,nx);
	
	u00=alloc2complex(nz,nx);
	r00=alloc2complex(nz,nx);
    /************************ compute alpha0 **************************************/
                             
    /*FFT---Cf=ifft(fft(C1)*fft(f))*/

    /* r0,r0=sre-L*sre */ 
    //BTmatrixmutip(g,u0,nx,nz,mutips);//G*u0,FFT1
	for(i=0;i<nx;i++){
		for(j=0;j<nz;j++){
			u00[i][j]=cmul(u0[i][j],diag[i][j]);			
		}				 					 										 		
	}
	
	BTTB_y(g,u00,nx,nz,mutips);//FFT2
				 
	for(i=0;i<nx;i++){
		for(j=0;j<nz;j++){
			r0[i][j]=cadd(csub(s[i][j],u0[i][j]),mutips[i][j]);
			r00[i][j]=cmul(r0[i][j],diag[i][j]);
		}				 					 										 		
	}
    /* alpha_r0,alpha_r0=(r0,Lr0)/||Lr0||^2*/
	//BTmatrixmutip(g,r0,nx,nz,mutipri);//G*r0,FFT1
	
	BTTB_y(g,r00,nx,nz,mutipri);	//FFT2
	
	for(i=0;i<nx;i++){
		for(j=0;j<nz;j++){
			Lr0[i][j]=csub(r0[i][j],mutipri[i][j]);			
		}				 					 										 		
	}
								 
	norm0=0;
    alpha0.r=0;
    alpha0.i=0;
                
    for(i=0;i<nx;i++){
		for(j=0;j<nz;j++){
			alpha0=cadd(alpha0,cmul(r0[i][j],conjg(Lr0[i][j])));			
            norm0 +=Lr0[i][j].r*Lr0[i][j].r+Lr0[i][j].i*Lr0[i][j].i;
		}
	}
				  					 
	alpha0.r=alpha0.r/norm0;
	alpha0.i=alpha0.i/norm0;
				
	free2complex(mutips);
	free2complex(Lr0); 			                                        		                        			  
    free2complex(mutipri);
    
    /*LU=Uinc,L=diag-A,Un=(I-alpha*L)Un-1+alpha*Uinc*/
    u=alloc2complex(nz,nx);
    diagc=alloc2complex(nz,nx);
    r=alloc2complex(nz,nx);
    Lr=alloc2complex(nz,nx);
    mutipu0=alloc2complex(nz,nx);
    mutipr0=alloc2complex(nz,nx);
    mutipr=alloc2complex(nz,nx);
	
	rr=alloc2complex(nz,nx);
	
	//utemp=alloc2complex(nz,nx);	
	//mutiputemp=alloc2complex(nz,nx);
    //mutipr00=alloc2complex(nz,nx);
    /************************ iter start **************************************/
	iter=0;
	while(iter<=itmax)
    {
        ++(iter);
		//printf("alpha0.r=%lf,alpha0.i=%lf ",alpha0.r,alpha0.i);	

        //BCCB_y(M,r0,nx,nz,mutipr00,-1);	   //M:BCCB_TC matrix
		for(i=0;i<nx;i++){  
			for(j=0;j<nz;j++){
				diagc[i][j]=cmul(M[i][j],alpha0);//M: diag matrix				
				u00[i][j]=cmul(u0[i][j],diag[i][j]);
				//r00[i][j]=cmul(r0[i][j],cmul(diagc[i][j],diag[i][j]));//M: diag matrix
				r00[i][j]=cmul(diag[i][j],cmul(M[i][j],r0[i][j]));//M: diag matrix
				//r00[i][j]=cmul(cmul(alpha0,diag[i][j]),mutipr00[i][j]);//M:BCCB_TC matrix
			}  				    
		}
				            
		//BTmatrixmutip(g,u0,nx,nz,mutipu0);//G*U0,FFT1
		//BTmatrixmutip(g,r0,nx,nz,mutipr0);//G*r0,FFT1
		BTTB_y(g,u00,nx,nz,mutipu0);
		BTTB_y(g,r00,nx,nz,mutipr0);

		#ifdef _OPENMP
        #pragma omp parallel for default(none) num_threads (2)	\
            private(i,j)	\
            shared(nx,nz,u,alpha0,s,diagc,u0,mutipu0,r,r0,mutipr0)
        #endif	
		
        /*M: diag matrix*/
		for(i=0;i<nx;i++){
			for(j=0;j<nz;j++){				
				u[i][j]=csub(u0[i][j],cmul(diagc[i][j],csub(u0[i][j],cadd(s[i][j],mutipu0[i][j]))));
				r[i][j]=csub(r0[i][j],csub(cmul(diagc[i][j],r0[i][j]),cmul(alpha0,mutipr0[i][j])));	
				//u[i][j]=cadd(u0[i][j],cmul(diagc[i][j],r0[i][j]));
				//r[i][j]=csub(r0[i][j],csub(cmul(diagc[i][j],r0[i][j]),cmul(alpha0,mutipr0[i][j])));									
			}
		}
		
		/*M:BCCB_TC matrix*/
		/*
		for(i=0;i<nx;i++){
			for(j=0;j<nz;j++){								
				utemp[i][j]=cmul(alpha0,csub(u0[i][j],cadd(s[i][j],mutipu0[i][j])));								
			}
		}
		
		BCCB_y(M,utemp,nx,nz,mutiputemp,-1);
		
		for(i=0;i<nx;i++){
			for(j=0;j<nz;j++){								
				u[i][j]=csub(u0[i][j],mutiputemp[i][j]);
				r[i][j]=csub(r0[i][j],csub(cmul(alpha0,mutipr00[i][j]),mutipr0[i][j]));				
			}
		}
		*/		
		
		//BTmatrixmutip(g,r,nx,nz,mutipr);//G*r
		for(i=0;i<nx;i++){  
			for(j=0;j<nz;j++){			
				rr[i][j]=cmul(r[i][j],diag[i][j]);
			}  				    
		}
		BTTB_y(g,rr,nx,nz,mutipr);

		#ifdef _OPENMP
        #pragma omp parallel for default(none) num_threads (2)	\
                private(i,j)	\
                shared(nx,nz,Lr,diag,r,mutipr)
        #endif								
		for(i=0;i<nx;i++){
			for(j=0;j<nz;j++){
				Lr[i][j]=csub(r[i][j],mutipr[i][j]);
			}
		}
							
		/*compute alpha_n*/
							
		err=0.0;
		norm=0.0;
        norm_S=0.0;
		alpha.r=0.0;
		alpha.i=0.0;	
													
		/*#ifdef _OPENMP
        #pragma omp parallel for default(none) num_threads (2)	\
                private(i,j)	\
                shared(nx,nz,Lr,r,s) \
                reduction(+:alpha.r,alpha.i,norm,norm_S,err )
        #endif*/
        for(i=0;i<nx;i++){
			for(j=0;j<nz;j++){
				alpha=cadd(alpha,cmul(r[i][j],conjg(Lr[i][j])));
				//alpha.r +=r[i][j].r*Lr[i][j].r+r[i][j].i*Lr[i][j].i;  
                //alpha.i +=r[i][j].i*Lr[i][j].r-r[i][j].r*Lr[i][j].i;
                norm +=Lr[i][j].r*Lr[i][j].r+Lr[i][j].i*Lr[i][j].i;
                norm_S +=s[i][j].r*s[i][j].r+s[i][j].i*s[i][j].i;
                err +=r[i][j].r*r[i][j].r+r[i][j].i*r[i][j].i;
				//err +=(u[i][j].r-u0[i][j].r)*(u[i][j].r-u0[i][j].r)
						//	+(u[i][j].i-u0[i][j].i)*(u[i][j].i-u0[i][j].i);
			}
		} 
    
		alpha0.r=alpha.r/norm;
		alpha0.i=alpha.i/norm;							                    
        err=sqrt(err/norm_S);
		//err=sqrt(err);
                            
        #ifdef _OPENMP
        #pragma omp parallel for default(none)  num_threads (2) \
                private(i,j) \
                shared(nx,nz,u,r,u0,r0)
        #endif
		for(i=0;i<nx;i++){
			for(j=0;j<nz;j++){
				u0[i][j].r=u[i][j].r;
				u0[i][j].i=u[i][j].i;
				r0[i][j].r=r[i][j].r;
				r0[i][j].i=r[i][j].i;
			}
		}  
		
		if(iter%100==0)
			printf("iter=%d, err=%lf\n",iter,err);
		if(err<EPS) break;
		
	}
	
	free2complex(r0);
    free2complex(u);
    free2complex(diagc);
    free2complex(r);
    free2complex(Lr);
    free2complex(mutipr0);
    free2complex(mutipr);
    free2complex(mutipu0);
    
	free2complex(u00);
	free2complex(r00);
	free2complex(rr);
	
	//free2complex(utemp);
	//free2complex(mutiputemp);
	//free2complex(mutipr00);
}
