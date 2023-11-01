void matmul_0(double *a,int r1,int c1,double *b,double *c)
{
     int i,j;
    
	 for(i=0;i<r1;i++)
	 {
		 for(j=0;j<c1;j++)
		 {
			c[i*c1+j] =a[i]*b[j];
		 }
	 }

} 


void matmul_1(double **a,int r1,int c1,double *b,double *c)
{
     int i,j;
     double temp;
    
	 for(i=0;i<r1;i++)
	 {
         temp=0;
		 for(j=0;j<c1;j++)
		 {
			temp+=a[i][j]*b[j];
		 }
         c[i] =temp;
	 }
}

void matmul_2(double **a,int r1,int c1,double **b,int c2,double **c)
{
     int i,j,k;
     double *temp;
     temp=alloc1double(c2);
	 
     for(k=0;k<c2;k++)
     {
         temp[k]=0; 
		
     }	 
	 
	 for(i=0;i<r1;i++)
	 {
                 for(k=0;k<c2;k++)
                 {
                    temp[k]=0;
                 }
		 for(j=0;j<c1;j++)
		 {
			for(k=0;k<c2;k++)
			{
			    temp[k] +=a[i][j]*b[j][k];
			}				
		 }
		 
		 for(k=0;k<c2;k++)
	         {
			c[i][k]=temp[k];			
		 }
		    
	 }
     free1double(temp);
}

complex *mattransp(complex *p,int m,int n)
 {   
	int i,j;

	
   complex* q = (complex *) malloc ( sizeof(complex) * m*n);
     
  
	
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			*(q+i+j*m)=*(p++);
		}
	    	
	}
	
    return q;
 }
 
 
/*Block Toeplitz matrix G(r,r') multiply vector y.___FFT1*/
void BTmatrixmultip(complex **G,complex **y,int Nx,int Nz,complex **multip)
{
	int ix,iz,i,N,Nfft;
	complex **GT=NULL,**yT=NULL,*GC=NULL,*yC=NULL;
	complex *GCfft=NULL,*yCfft=NULL,*multipfft=NULL;
	
	/*The matrix G inserts into a 2*Nx*Nz order Toeplitz matrix GT*/
	   GT=alloc2complex(2*Nz,Nx);
	   yT=alloc2complex(2*Nz,Nx);
	   
	   for(ix=0;ix<Nx;ix++){
	   	   for(iz=0;iz<(2*Nz);iz++){
	   	  	   if(ix==(Nx-1)){
	   	  	  	   if(iz<Nz){
	   	  	  	    	GT[ix][iz]=G[ix][iz];
	   	  	  	    	yT[ix][iz]=y[ix][iz];
					}else if(iz==Nz){
						GT[ix][iz].r=0;
						GT[ix][iz].i=0;
						yT[ix][iz].r=0;
						yT[ix][iz].i=0;
					}else{
						GT[ix][iz].r=0;
						GT[ix][iz].i=0;
						yT[ix][iz].r=0;
						yT[ix][iz].i=0;
					}
				}else{
					if(iz<Nz){
	   	  	  	    	GT[ix][iz]=G[ix][iz];
	   	  	  	    	yT[ix][iz]=y[ix][iz];
					}else if(iz==Nz){
						GT[ix][iz].r=0;
						GT[ix][iz].i=0;
						yT[ix][iz].r=0;
						yT[ix][iz].i=0;
					}else{
						GT[ix][iz]=G[ix+1][2*Nz-iz];
						yT[ix][iz].r=0;
						yT[ix][iz].i=0;
					}
				}
				
		   }
	    }
	    /*printf("\n");
	    printf("GT\n");
	    for(ix=0;ix<Nx;ix++){
	   	   for(iz=0;iz<(2*Nz);iz++){
	           printf("%4.2lf+%4.2lf*i ",GT[ix][iz].r,GT[ix][iz].i);
	        }
	    }
	    printf("\n");
	    printf("yT\n");
	    for(ix=0;ix<Nx;ix++){
	   	   for(iz=0;iz<(2*Nz);iz++){
	           printf("%4.2lf+%4.2lf*i ",yT[ix][iz].r,yT[ix][iz].i);
	        }
	    }*/
	/*The Toeplitz matrix GT inserts into a (2*Nz)*(2*Nx) order Circulant matrix GC*/
	   GC=alloc1complex(4*Nz*Nx);
       yC=alloc1complex(4*Nz*Nx);
       
       for(ix=0;ix<(2*Nx);ix++){
       	   for(iz=0;iz<(2*Nz);iz++){
       	   	   if((ix*2*Nz+iz)<(2*Nz*Nx)){
       	   	   	    GC[ix*2*Nz+iz]=GT[ix][iz];
					yC[ix*2*Nz+iz]=yT[ix][iz];
				}else if((ix*2*Nz+iz)==(2*Nz*Nx)){
					GC[ix*2*Nz+iz].r=0;
					GC[ix*2*Nz+iz].i=0;
					yC[ix*2*Nz+iz].r=0;
					yC[ix*2*Nz+iz].i=0;
				}else{
					GC[ix*2*Nz+iz]=GC[4*Nx*Nz-(ix*2*Nz+iz)];
					yC[ix*2*Nz+iz].r=0;
					yC[ix*2*Nz+iz].i=0;
				}
			}			  
	   }
	   /*printf("\n");
	   printf("GC\n");
	    for(ix=0;ix<2*Nx;ix++){
	   	   for(iz=0;iz<(2*Nz);iz++){
	           printf("%4.2lf+%4.2lf*i ",GC[ix*2*Nx+iz].r,GC[ix*2*Nx+iz].i);
	        }
	    }
	    printf("\n");
	    printf("yC\n");
	    for(ix=0;ix<2*Nx;ix++){
	   	   for(iz=0;iz<(2*Nz);iz++){
	           printf("%4.2lf+%4.2lf*i ",yC[ix*2*Nx+iz].r,yC[ix*2*Nx+iz].i);
	        }
	    }*/
	    free2complex(GT);
	    free2complex(yT);
	
	/*FFT---Cf=ifft(fft(C1)*fft(f))*/
	    N=4*Nz*Nx;
		//printf("N=%d\n",N);
	    Nfft=npfa(N);
        if (Nfft > PFA_MAX)   printf("Padded Nfft=%d--too big\n", Nfft);
        //printf("Nfft=%d \n",Nfft);
        
        GCfft=alloc1complex(Nfft);
        yCfft=alloc1complex(Nfft);
        multipfft=alloc1complex(Nfft);
        //multipC=alloc1complex(4*Nz1*Nx1);
        
        for(i=0;i<Nfft;i++){
            if(i<(4*Nz*Nx)){
                GCfft[i]=GC[i];
                yCfft[i]=yC[i];
		   }else{
			    GCfft[i].r=0;
		    	GCfft[i].i=0;
                yCfft[i].r=0;
                yCfft[i].i=0;
		   }
	    }
	    
	    free1complex(GC);
	    free1complex(yC);
	    
	    pfacc(1,Nfft,GCfft); 
    	pfacc(1,Nfft,yCfft);
    	
    	for(i=0;i<Nfft;i++){
            multipfft[i].r=(GCfft[i].r*yCfft[i].r-GCfft[i].i*yCfft[i].i)/sqrt(Nfft);
            multipfft[i].i=(GCfft[i].r*yCfft[i].i+GCfft[i].i*yCfft[i].r)/sqrt(Nfft);
	    }
	    
	    free1complex(GCfft);
	    free1complex(yCfft);
	    
	    pfacc(-1,Nfft,multipfft);
	    
        for(ix=0;ix<Nx;ix++){
        	for(iz=0;iz<(2*Nz);iz++){
        		if(iz<Nz){
        			multip[ix][iz].r=multipfft[ix*2*Nz+iz].r/sqrt(Nfft);
        		    multip[ix][iz].i=multipfft[ix*2*Nz+iz].i/sqrt(Nfft);
				}else
        		   continue;
			}
		}
		
		free1complex(multipfft);
	        		    	   
 } 



/*a Nx*Nz order BTTB matrix G multiply vector y.__FFT2*/
void BTTB_y(complex **GT,complex **y,int Nx,int Nz,complex **multip_GT) 
{
	int ix,iz;
	complex **GC=NULL,**yC=NULL,**multip_GC=NULL;
	
	/*a Nx*Nz order BTTB matrix G inserts into a 2Nx*2Nz order BCCB matrix*/
	GC=alloc2complex(2*Nz,2*Nx);
	yC=alloc2complex(2*Nz,2*Nx);
	multip_GC=alloc2complex(2*Nz,2*Nx);
	
	for(ix=0;ix<(2*Nx);ix++){
		if(ix<Nx){
			for(iz=0;iz<(2*Nz);iz++){
				if(iz<Nz){
				   	GC[ix][iz]=GT[ix][iz];
				   	yC[ix][iz]=y[ix][iz];
				}else if(iz==Nz){
					GC[ix][iz]=cmplx(0.0,0.0);
				   	yC[ix][iz]=cmplx(0.0,0.0);
				}else if(iz>Nz){
					GC[ix][iz]=GC[ix][2*Nz-iz];
				   	yC[ix][iz]=cmplx(0.0,0.0);
				}			
		    }
		}else if(ix==Nx){
			for(iz=0;iz<(2*Nz);iz++){
				GC[ix][iz]=cmplx(0.0,0.0);
				yC[ix][iz]=cmplx(0.0,0.0);
			}
		}else if(ix>Nx){
			for(iz=0;iz<(2*Nz);iz++){
				GC[ix][iz]=GC[2*Nx-ix][iz];//2022.11.10 modify
				yC[ix][iz]=cmplx(0.0,0.0);
			}
		}		
	}//BCCB matrix construction finish
	
	/*BCCB matrix GC multiply vector yC*/
	BCCB_y(GC,yC,2*Nx,2*Nz,multip_GC,1);
	
	for(ix=0;ix<Nx;ix++){
		for(iz=0;iz<Nz;iz++){
			multip_GT[ix][iz]=multip_GC[ix][iz];
			
		}
	}	
	free2complex(GC);
	free2complex(yC);
	free2complex(multip_GC);	
}

/*make a Strang BCCB precondition matrix from a BTTB matrix*/
void BCCB_S(complex **GT,int Nx,int Nz,complex **Pre_S)
{
	int ix,iz;
	complex **pre_s=NULL;
	
	pre_s=alloc2complex(Nz,Nx);
		
	for(ix=0;ix<Nx;ix++){
		for(iz=0;iz<Nz;iz++){
			if(iz<(Nz/2+1)){
				pre_s[ix][iz]=GT[ix][iz];
			}else{
				pre_s[ix][iz]=GT[ix][Nz-iz];
			}
		}
	}
	
	for(ix=0;ix<Nx;ix++){
		if(ix<(Nx/2+1)){
			for(iz=0;iz<Nz;iz++){
				Pre_S[ix][iz]=pre_s[ix][iz];
			}
		}else{
			for(iz=0;iz<Nz;iz++){
				Pre_S[ix][iz]=pre_s[Nx-ix][iz];
			}
		}		
	}
	
	free2complex(pre_s);
}

/*make a T.Chen BCCB precondition matrix from a BTTB matrix*/
void BCCB_TC(complex **GT,int Nx,int Nz,complex **Pre_C)
{
	int ix,iz;
	complex **pre_c=NULL;
	
	pre_c=alloc2complex(Nz,Nx);
	
	for(ix=0;ix<Nx;ix++){
		for(iz=0;iz<Nz;iz++){
			if(iz==0){
				pre_c[ix][iz]=GT[ix][iz];
			}else{
				pre_c[ix][iz]=crmul(cadd(crmul(GT[ix][iz],(Nz-iz)),crmul(GT[ix][Nz-iz],iz)),1.0/Nz);	
			} 
		}
	}
	
	for(ix=0;ix<Nx;ix++){
		if(ix==0){
			for(iz=0;iz<Nz;iz++){
				Pre_C[ix][iz]=pre_c[ix][iz];	
			} 
		}else{
			for(iz=0;iz<Nz;iz++){
				Pre_C[ix][iz]=crmul(cadd(crmul(pre_c[ix][iz],(Nx-ix)),crmul(pre_c[Nx-ix][iz],ix)),1.0/Nx);
			}
		}
    }
    
    free2complex(pre_c);
}

/*BCCB matrix GC multiply vector y,or the inverse of preconditon BCCB matrix multiply vector.___FFT2*/
/* GC*y ,or M(-1)*b */
void BCCB_y(complex **GC,complex **yC,int Nx,int Nz,complex **multip_GC,int isign)
{
	int ix,iz,Nxfft,Nzfft;
	complex **GCfft=NULL,**yCfft=NULL,**multipfft=NULL;
	
	//Nzfft=npfa(Nz);
	//Nxfft=npfa(Nx);
	Nzfft=Nz;//guarantee the result right 2023.05.11
	Nxfft=Nx;
	if ((Nzfft > PFA_MAX) || (Nxfft > PFA_MAX))   
	    printf("Padded Nzfft=%d or Nxfft=%d--too big\n", Nzfft,Nxfft);
	
	GCfft=alloc2complex(Nzfft,Nxfft);
    yCfft=alloc2complex(Nzfft,Nxfft);
    multipfft=alloc2complex(Nzfft,Nxfft);
	
	for(ix=0;ix<Nxfft;ix++){
		for(iz=0;iz<Nzfft;iz++){
			if((ix<Nx) && (iz<Nz)){
				GCfft[ix][iz]=GC[ix][iz];
				yCfft[ix][iz]=yC[ix][iz];
			}else{
				GCfft[ix][iz]=cmplx(0.0,0.0);
				yCfft[ix][iz]=cmplx(0.0,0.0);
			}
			
			multipfft[ix][iz]=cmplx(0.0,0.0);
		}
	}
    /*
    for(ix=0;ix<Nx;ix++){
    	for(iz=0;iz<Nz;iz++){   		
    		GCfft[ix][iz]=GC[ix][iz];
    		yCfft[ix][iz]=yC[ix][iz];						  		
		}
	}*/
	
	pfa2cc(-1,1,Nzfft,Nxfft,GCfft[0]); 
	pfa2cc(-1,2,Nzfft,Nxfft,GCfft[0]);
	pfa2cc(-1,1,Nzfft,Nxfft,yCfft[0]);
	pfa2cc(-1,2,Nzfft,Nxfft,yCfft[0]);
	/*
	for(ix=0;ix<Nxfft;ix++){
			for(iz=0;iz<Nzfft;iz++){
				printf("GCfft[%d][%d].r=%lf\n",ix,iz,GCfft[ix][iz].r);			
			}
		}
	*/
	/**********************old codes*******************/
	/*if(isign>0){
		for(ix=0;ix<Nxfft;ix++){
			for(iz=0;iz<Nzfft;iz++){
				multipfft[ix][iz]=crmul(cmul(GCfft[ix][iz],yCfft[ix][iz]),1/sqrt(Nxfft*Nzfft));//modified 20220621
			}
		}
	}else if(isign<0){
		for(ix=0;ix<Nxfft;ix++){
			for(iz=0;iz<Nzfft;iz++){
				multipfft[ix][iz]=crmul(cdiv(yCfft[ix][iz],GCfft[ix][iz]),1/sqrt(Nxfft*Nzfft));//modified 20220621
			}
		}
	}
	
	pfa2cc(1,2,Nzfft,Nxfft,multipfft[0]);
	pfa2cc(1,1,Nzfft,Nxfft,multipfft[0]);
	
	for(ix=0;ix<Nx;ix++){
		for(iz=0;iz<Nz;iz++){
			multip_GC[ix][iz]=crmul(multipfft[ix][iz],1/sqrt(Nxfft*Nzfft));//modified 20220621
		}
	}*/
	
	/*************new codes 20220621***************/	
	if(isign>0){
		for(ix=0;ix<Nxfft;ix++){
			for(iz=0;iz<Nzfft;iz++){
				multipfft[ix][iz]=cmul(GCfft[ix][iz],yCfft[ix][iz]);//modified 20220621				
			}
		}
	}else if(isign<0){
		for(ix=0;ix<Nxfft;ix++){
			for(iz=0;iz<Nzfft;iz++){
				multipfft[ix][iz]=cdiv(yCfft[ix][iz],GCfft[ix][iz]);//modified 20220621
			}
		}
	}
	
	pfa2cc(1,1,Nzfft,Nxfft,multipfft[0]);
	pfa2cc(1,2,Nzfft,Nxfft,multipfft[0]);
	
	for(ix=0;ix<Nx;ix++){
		for(iz=0;iz<Nz;iz++){
			multip_GC[ix][iz]=crmul(multipfft[ix][iz],1.0/(Nxfft*Nzfft));//modified 20220621
		}
	}
	
	free2complex(GCfft);
	free2complex(yCfft);
	free2complex(multipfft);
	
}


/*void arymul1(int a[4][5], int b[5][3], int c[4][3])
{
	int i, j, k;
	int temp[3] = {0};
	for(i = 0; i < 4; i++){
		for(k = 0; k < 3; k ++)
			temp[k] = 0;
		for(j = 0; j < 5; j++){//every element in current cow
			for(k = 0; k < 3; k++){
				temp[k] += a[i][j] * b[j][k];
			}
		}
		for(k = 0; k < 3; k++){
			c[i][k] = temp[k];
			printf("%d/t", c[i][k]);
		}
		printf("%d/n");
	}
}*/
