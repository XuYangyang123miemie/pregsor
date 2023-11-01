#include"main.h"
#include"alloc.c"
#include"ricker1_wavelet.c"
#include"pfafft.c"
#include"green.c"
#include"greenCW.c"
#include"gauss.c"
#include"gaussCW.c"
#include"diag.c"
#include"interp.c"

//#include"matrix.c"
//#include"asolve.c"
//#include"linbcg.c"
//#include"linbcg_openmp.c"

#include"matrixBT.c"
#include"linssor_CW_PrePL_0.c"
#include"linssor_CW_PreO.c"
#include"linssor.c"
//#include"linbcg_fft2_CW.c"
//#include"linssor_preOCW.c"
//#include"linbcg_pc_fft2.c"
//#include"linbcgstab_pc_fft2.c"

#include"mpi.h"
#define FROM_MASTER 1
#define FROM_SLAVE 2

#ifdef _OPENMP
#include <omp.h>
#endif

MPI_Status status;

//void omp_set_num_threads (int t);


int main(int argc,char *argv[])
{
    /* parameters declaration:supplement */
    int nt,ntw,nxs,nxo,nx,nz,nx1,nz1,ntfft,nw;
    int it,ixo,ix,iz,ix1,iz1,iw,i,j;//ixs,m,n;
    int nbx,nbz;
    int itol,iter,iteru,itmax;
	
    int process_id,process_num,slave_num,dest,source,rows_aver,rows_w,remainder,offset;//MPI para
    MPI_Datatype COMPLEX;
    int blocklens[2];
    MPI_Aint indices[2];
    MPI_Datatype old_types[2]; 
    
    double fpeak,ft,dt,fxs,dxs,zs,zg,fxo,dxo,fw,dw,fx,fz,xr,zr,rs,v0,vmin,vmax,ornum,h;
    //double ex,ez,exs;
	double xs,xg;
    double k0,w;
    double tempv;    
	double a,b;
    //double norm,norm0,norm_S,
	double err,uerr,uerrmax;
    double *wavelet,*rt,*umt,*trace,**umt2;
    float **v;
    double **or;
    double **window;
    
    complex **s,**g,**gpre,**sg,**u0,**M,**oor,**diag,**green,**mutiptest;   //20201112 add **gpre
	complex **s0,**g0,**r0,**d0,**u1;
    complex dr,kk,gtemp,greentemp,*ct,**um,**umtemp,*p;
	dcomplex epsilon,kk0,pkk0,hank;//20201112 add epsilon,kk0
    //complex test;
    
    FILE *outfp,*outfpsg0,*outfpsg1,*outfpsg2,*outfpsg3,*outfpsg4,*outfpsg5,*outfpsg6;
	FILE *outfpfrewf1,*outfpfrewf2,*fpara,*vfp;

    //clock_t start,end;
    double start,end,runtime;

    //int t;
    //t=omp_get_num_procs ();

    static int js2[2]={4,4};  //update js[2]={nxr,nzr}
    ///static int jsfft[4]={2,2,1,1};
    
    /******************************** input and output files***************************************/
    
    	if((fpara=fopen("parafile","r+"))==NULL)
    	{
       	 	fprintf(stderr,"can not open the parameter file\n");
        	return 1;
    	}
    
            /*********************** get information for seismogram traces ***********************/
    
     		fscanf(fpara,"nt=%d,dt=%lf,ft=%lf,fpeak=%lf\n",&nt,&dt,&ft,&fpeak);
     		       		
    		fscanf(fpara,"nxs=%d,dxs=%lf,fxs=%lf,zs=%lf\n",&nxs,&dxs,&fxs,&zs);
     		
     		//exs=fxs+(nxs-1);
	
     		fscanf(fpara,"nxo=%d,dxo=%lf,fxo=%lf,zg=%lf\n",&nxo,&dxo,&fxo,&zg);    			
	
    		fscanf(fpara,"v0=%lf\n",&v0);     		
	
     		fscanf(fpara,"nx=%d,fx=%lf\n",&nx,&fx);     		
     
     		fscanf(fpara,"nz=%d,fz=%lf\n",&nz,&fz);     		
     
     		nx1=nx-1;
     		nz1=nz-1;
     		//ex=fx+(nx-1)*h;
	 		//ez=fz+(nz-1)*h;

     		fscanf(fpara,"h=%lf,nbx=%d,nbz=%d\n",&h,&nbx,&nbz); 

     		fscanf(fpara,"itmax=%d,itol=%d\n",&itmax,&itol);     		

     		fscanf(fpara,"a=%lf,b=%lf\n",&a,&b);
     		
     		
		fclose(fpara);
		
		/*matrix multiply Nfft*/
		int NX,NZ,NXFFT,NZFFT;
		NX=2*(nx-1);
		NZ=2*(nz-1);
		NXFFT=npfa(NX);
		NZFFT=npfa(NZ);		 
		
		if((vfp=fopen("vel","r+"))==NULL)
    	{
        	fprintf(stderr,"can not open the velocity file\n");
        	return 1;
    	}
			
	/*********************** check source and receiver coordinates *************************/
 			/*for (ixs=0; ixs<nxs; ++ixs) 
			{
				xs = fxs+ixs*dxs;
				for (ixo=0; ixo<nxo; ++ixo) 
				{
					xg = fxo+ixo*dxo;
					if (fx>xs || ex<xs || fx>xg || ex<xg) 
			        printf("shot or receiver lie outside of specified (x,z) grid\n");
				}
			} */

    /*********************** Set up pfa fft **************************************************/
			  ntfft = npfaro(nt, LOOKFAC * nt);
    		      if (ntfft >= SU_NFLTS || ntfft >= PFA_MAX)   printf("Padded nt=%d--too big", ntfft);
			  nw = ntfft/2 + 1;
			  //dw = 1/(ntfft*dt);
			  dw = 1.0/(ntfft*dt);
    		  //fw=0.5;			
    /************************************ make wavelet ************************************/
	  	   
    		/* Allocate fft arrays */
			  rt = alloc1double(ntfft);
			  ct = alloc1complex(nw);

    		/*new20190925:make wavelet*/    
    		  wavelet=ricker(fpeak,dt,&ntw);
    		  memset((void *)rt,0,ntfft*FSIZE);

    		    for(it=0;it<ntw;it++)
    		    {
    		        rt[(it-ntw/2+ntfft)%ntfft]=wavelet[it];          
    		    }
	    
    		/**FFT**/
    		  pfarc(-1, ntfft, rt, ct);  
        
			  free1double(wavelet);
    		  free1double(rt);		

    /***********************************velocity************************************/
  
			/*alloc space*/
    		  v=alloc2float(nz,nx);        
			  or=alloc2double(nz1,nx1);  
        
    		/*read velocity file*/
    		  fread(v[0],sizeof(double),nz*nx,vfp);
    		  
    /*********************Select the max and min vel*********************/
            
			vmin=vmax=v[0][0];
			for(j=0;j<nz;j++)
			{
				for(i=0;i<nx;i++)
				{
					vmin=((vmin<v[i][j])?vmin:v[i][j]);
					vmax=((vmax>v[i][j])?vmax:v[i][j]);
				}
			}
                      
    		/****or****/
    		for(ix=0;ix<nx1;ix++)
    		{           
			    for(iz=0;iz<nz1;iz++)
			    {
    		        //tempv=0.0;
    		        //tempv=0.25*(v[ix][iz]+v[ix][iz+1]+v[ix+1][iz]+v[ix+1][iz+1]);
    		        //or[ix][iz]=v0*v0/(tempv*tempv)-1;
					or[ix][iz]=v0*v0/(v[ix][iz]*v[ix][iz])-1.0;
					//if(or[ix][iz]==0) or[ix][iz]=-0.001;//20221015 16:35 add             
    		    }   	   
    		}
    		
    		ornum=(vmin*vmin)/(vmax*vmax)-1.0;//20201103 09:24
		        //printf("or[200][10]=%f,or[300][10]=%f\n",or[200][10],or[300][10]);      
		    free2float(v); 
	
	/*******************************Hanning window**************************************/
	
		    window=alloc2double(nz1,nx1);
		    for(ix=0;ix<nx1;ix++){
	  	 		for(iz=0;iz<nz1;iz++){
					if(iz<nz1-nbx){
						if((ix>(nbx-1))&&(ix<(nx-nbx-1)))
							window[ix][iz]=1;
						else if(ix<nbx)
							window[ix][iz]=(1+cos(PI*(nbx-ix)/nbx))/2;
						else if(ix>(nx-nbx-2))
							window[ix][iz]=(1+cos(PI*(ix-nx+2+nbx)/nbx))/2;	  	 	
					}else{
						if((ix>(nbx-1))&&(ix<(nx-nbx-1)))
							window[ix][iz]=(1+cos(PI*(iz-nz+2+nbz)/nbz))/2;
						else if(ix<nbx)
							window[ix][iz]=((1+cos(PI*(nbx-ix)/nbx))/2)*((1+cos(PI*(iz-nz+2+nbz)/nbz))/2);
						else if(ix>(nx-nbx-2))
							window[ix][iz]=((1+cos(PI*(ix-nx+2+nbx)/nbx))/2)*((1+cos(PI*(iz-nz+2+nbz)/nbz))/2);						
					}
				}
			} 
    
    /**************************************MPI*******************************************/
	
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&process_id);
    MPI_Comm_size(MPI_COMM_WORLD,&process_num);
	slave_num=process_num-1;
   
    
    if(process_id==0){    			
     	printf("nt=%d,dt=%lf,ft=%lf,fpeak=%lf\n",nt,dt,ft,fpeak);
		printf("nxs=%d,dxs=%lf,fxs=%lf,zs=%lf\n",nxs,dxs,fxs,zs);
		printf("nxo=%d,dxo=%lf,fxo=%lf,zg=%lf\n",nxo,dxo,fxo,zg);
		printf("v0=%lf\n",v0);
		printf("nx=%d,fx=%lf\n",nx,fx);
		printf("nz=%d,fz=%lf\n",nz,fz);
		printf("h=%lf,nbx=%d,nbz=%d\n",h,nbx,nbz);
		printf("itmax=%d,itol=%d\n",itmax,itol);
		printf("a=%lf,b=%lf\n",a,b);
		printf("nt=%d,dt=%f,fpeak=%f,ntfft=%d,nw=%d\n",nt,dt,fpeak,ntfft,nw);
		printf("vmin=%f,vmax=%f\n",vmin,vmax);//20201103 09:27	
		printf("NX=%d,NZ=%d,NXFFT=%d,NZFFT=%d\n", NX,NZ,NXFFT,NZFFT);
	}
	
	/*One value of each type*/
	blocklens[0]=1;
	blocklens[1]=1;
	/*The base type*/
	old_types[0]=MPI_DOUBLE;
	old_types[1]=MPI_DOUBLE;
	/*location*/
	//MPI_Address(&complex.r,&indices[0]);
	//MPI_Address(&complex.i,&indices[1]);
	/*new COMPLEX:offset*/
	indices[1]=8;
	//indices[1]=indices[1]-indices[0];
	indices[0]=0;
	
	/*generate new MPI_Datatype COMPLEX*/
	MPI_Type_struct(2,blocklens,indices,old_types,&COMPLEX);
	
	/*commit*/
	MPI_Type_commit(&COMPLEX);
	
	xs=fxs;

	MPI_Barrier(MPI_COMM_WORLD);
	
	
	/*******************************primary process*********************************/
	
	    if(process_id==0){
	    	
	    	outfp=fopen("NySOR_T_SCpreo.txt","w+");  //record
      		fseek(outfp,(off_t) 0,SEEK_CUR);

      		//outfpsg1=fopen("NySOR_Wr_SCpreo.txt","w+"); //reciever
      		//fseek(outfpsg1,(off_t) 0,SEEK_CUR);

      		//outfpsg2=fopen("NySOR_Wi_SCpreo.txt","w+"); //reciever
      		//fseek(outfpsg2,(off_t) 0,SEEK_CUR);
      		
      		/*alloc space*/ 
			um=alloc2complex(nxo,nw);
			
			/*initialize*/ 
			for(iw=0;iw<nw;iw++)
			{
				for(ixo=0;ixo<nxo;ixo++)
				{
					um[iw][ixo].r=0;
					um[iw][ixo].i=0;
				}
			}
      		
           /*time start*/
            start=MPI_Wtime();
      		//start=clock();
	        								
			rows_aver=(nw/2)/slave_num;//compute 160 frequency samples
			remainder=(nw/2)%slave_num;
			offset=0;
			
			for(dest=1;dest<=slave_num;dest++)
			{
				rows_w=(dest<=remainder)?rows_aver+1:rows_aver;
				if(dest==1 || dest==slave_num) printf("sending %d w-rows to process %d\n",rows_w,dest);
				MPI_Send(&offset, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
				MPI_Send(&rows_w, 1, MPI_INT, dest, FROM_MASTER, MPI_COMM_WORLD);
				offset += rows_w;
			}
			
			for(source=1; source<=slave_num; source++)
			{
				MPI_Recv(&offset, 1, MPI_INT, source, FROM_SLAVE, MPI_COMM_WORLD, &status); 
				MPI_Recv(&rows_w, 1, MPI_INT, source, FROM_SLAVE, MPI_COMM_WORLD, &status); 
				MPI_Recv(&um[offset][0], rows_w*nxo, COMPLEX, source, FROM_SLAVE, MPI_COMM_WORLD, &status); 
				//MPI_Barrier(MPI_COMM_WORLD);
			}
					   
			printf("Freqency domain finish\n"); 
	        /*time end*/
           //end=clock();            
            end=MPI_Wtime();         
            runtime=(double)(end-start);
           	printf("%f seconds\n",runtime);
           	
			/*
           	for(iw=0;iw<nw;iw++){
           		for(ixo=0;ixo<nxo;ixo++){
           			fprintf(outfpsg1,"%16.8f",um[iw][ixo].r);     
                    fprintf(outfpsg2,"%16.8f",um[iw][ixo].i);
				}
				fprintf(outfpsg1,"\n");   
                fprintf(outfpsg2,"\n");
			}*/
           	
           	umt2=alloc2double(nxo,ntfft);
           	pfa2cr(-1,2,nxo,ntfft,um[0],umt2[0]); //2022.11.10 old:pfa2cr(-1,2,....)
           	
           	for (it=0; it<nt; ++it){
           		for(ixo=0;ixo<nxo;ixo++){
           			fprintf(outfp,"%16.8f",umt2[it][ixo]/ntfft);
				}
				fprintf(outfp,"\n");	
			}
			
			free2double(umt2);       	           
            free2complex(um);
	        //free1complex(p);
	        //fclose(outfpsg1);
            //fclose(outfpsg2);
	        fclose(outfp);

		}
		
		/*****************************subprocess*******************************************/
		
		else if(process_id > 0){
		
			MPI_Recv(&offset, 1, MPI_INT, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
		    MPI_Recv(&rows_w, 1, MPI_INT, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
		     			
	        umtemp=alloc2complex(nxo,rows_w);
			
		    /* loop over frequency samples */
			
	        for(iw=0;iw<rows_w;iw++)
		    {
                fw=35.0;
                w=fw+(iw+offset)*dw;	
				
                if(process_id==1) printf("iw=%d,w=%lf in process %d\n",iw,w,process_id);               

                 k0=2*PI*w/v0; 				 				
				 
				 /*****new*****/
				 epsilon=dcsqrt(dcmplx(1.0,a*fabs(ornum)));
				 
				 kk0=dcrmul(epsilon,k0);//20201112 complex wavenumber kk0
				 pkk0=dcmul(kk0,kk0);
				 
                /* allocate space */				
                s=alloc2complex(nz1,nx1);                				
	            g=alloc2complex(nz1,nx1);
				s0=alloc2complex(nz1,nx1);                				
	            g0=alloc2complex(nz1,nx1);
				//gpre=alloc2complex(nz1,nx1);//20201112 add
				oor=alloc2complex(nz1,nx1);
				diag=alloc2complex(nz1,nx1);
				//green=alloc2complex(nz1*nx1,nz1*nx1);
				//mutiptest=alloc2complex(nz1,nx1);
				
	            /* loop over r */ 	
			                				
                for(ix=0;ix<nx1;ix++)
		        {
			        for(iz=0;iz<nz1;iz++)
		    	    {
			           /* compute Si, G(w,r,rs)=i*H0(k0*r) */
                        xr=fx+0.5*h+ix*h;
                        zr=fz+0.5*h+iz*h;
                        rs=sqrt((xs-xr)*(xs-xr)+(zs-zr)*(zs-zr));//20201208 09:08
						
						s0[ix][iz]=crmul(cmplx(Greenr(k0*rs),Greeni(k0*rs)),window[ix][iz]);
						//s[ix][iz]=crmul(cmul(ct[iw+offset],cmplx(GreenCWr(dcrmul(kk0,rs)),GreenCWi(dcrmul(kk0,rs)))),window[ix][iz]);
						s[ix][iz]=crmul(cmplx(GreenCWr(dcrmul(kk0,rs)),GreenCWi(dcrmul(kk0,rs))),window[ix][iz]);
						//s[ix][iz]=cmplx(GreenCWr(dcrmul(kk0,rs)),GreenCWi(dcrmul(kk0,rs)));
						//hank=hankelH1(0,dcrmul(kk0,rs));
						//s[ix][iz].r=(-0.25)*hank.i;
						//s[ix][iz].i=0.25*hank.r;
						//s[ix][iz]=crmul(s[ix][iz],window[ix][iz]);
												
						/*****new******/
						oor[ix][iz]=cmplx(or[ix][iz],(-a)*fabs(ornum));
						
				        diag[ix][iz]=oor[ix][iz];//cinv(oor[ix][iz]);							
						
				        /* compute the first column of the integral operator G, L=diag-G*/
				        /* loop over r' */
				        if((ix==0)&&(iz==0)){			        	
			               for(ix1=0;ix1<nx1;ix1++)
			                {                        
			                   for(iz1=0;iz1<nz1;iz1++)
			                   {                               					 					            
				                 //20201208 09:08
								 g0[ix1][iz1].r=k0*k0*gaus2(xr,zr,fx,fz,h,ix1,iz1,k0,2,js2,ss2,fre2);
				                 g0[ix1][iz1].i=k0*k0*gaus2(xr,zr,fx,fz,h,ix1,iz1,k0,2,js2,ss2,fim2);	
								 
				                 g[ix1][iz1].r=k0*k0*gaus2_CW(xr,zr,fx,fz,h,ix1,iz1,kk0,2,js2,ss2_CW,fre2_CW);
				                 g[ix1][iz1].i=k0*k0*gaus2_CW(xr,zr,fx,fz,h,ix1,iz1,kk0,2,js2,ss2_CW,fim2_CW);								 
								 //gtemp=cmplx(gaus2_CW(xr,zr,fx,fz,h,ix1,iz1,kk0,2,js2,ss2_CW,fre2_CW)
								 //				,gaus2_CW(xr,zr,fx,fz,h,ix1,iz1,kk0,2,js2,ss2_CW,fim2_CW));							 								 
                                 //g[ix1][iz1]=cmul(cmplx(pkk0.r,pkk0.i),gtemp);													               						                       
			                    }							
			                }//r'
			            }
						/*Green output test*/
						/*
						for(ix1=0;ix1<nx1;ix1++)
			            {                        
			                for(iz1=0;iz1<nz1;iz1++)
			                {
								greentemp=cmplx(gaus2_CW(xr,zr,fx,fz,h,ix1,iz1,kk0,2,js2,ss2_CW,fre2_CW)
												,gaus2_CW(xr,zr,fx,fz,h,ix1,iz1,kk0,2,js2,ss2_CW,fim2_CW));	
								green[ix*nz1+iz][ix1*nz1+iz1]=cmul(cmplx(pkk0.r,pkk0.i),greentemp);	
							}
						}
						*/
		            }                   
		        }//r
		        
				
              /************************ PBCG **************************************/
                u0=alloc2complex(nz1,nx1);
                M=alloc2complex(nz1,nx1);
                
                /*precondition matrix M is a Strang/T.Chen BCCB matrix */ 
				//BCCB_S(g,nx1,nz1,M);  //20201103 8:36
				//BCCB_TC(g,nx1,nz1,M);
                
                /*#ifdef _OPENMP
                #pragma omp parallel for default(none) num_threads (8)	\
                     private(i,j)	\
                     shared(nx1,nz1,u0,s,or,M,diag,g)
                #endif*/
                 for(i=0;i<nx1;i++){
				    for(j=0;j<nz1;j++){
									    						 
						/*****new*****/
						 if(a==0.0 || w<5) 
						 	 M[i][j]=cinv(csub(cmplx(1.0,0.0),cmul(g[0][0],oor[i][j])));//
							 //M[i][j]=cmplx(1.0,or[i][j]/(1.0*fabs(ornum)));
						 else 
						 	 M[i][j]=cmplx(1.0,or[i][j]/(b*fabs(ornum)));

						//M[i][j]=cinv(csub(cmplx(1.0,0.0),cmul(g[0][0],oor[i][j])));
						u0[i][j]=cmul(M[i][j],s[i][j]);
						//u0[i][j]=s[i][j];
						/*precondition BCCB_TC*/
						//if(i==0&&j==0) M[i][j]=csub(cmplx(1.0,0.0),crmul(M[i][j],abs(ornum)));
						//else M[i][j]=cneg(crmul(M[i][j],abs(ornum)));
						/*mutiptest code*/
						//M[i][j]=s[i][j];//cmul(s[i][j],oor[i][j]);
					}				 					 										 		
				}
				/*testcode*/
				//BTTB_y(g,s,nx1,nz1,mutiptest);
				/*testcode output*/
				
				if(w==fw){
					outfpsg0=fopen("matrixM1_1.txt","w+"); //oor
					fseek(outfpsg0,(off_t) 0,SEEK_CUR);
					
					//outfpsg1=fopen("matrixtest.txt","w+"); //mutiptest
					//fseek(outfpsg1,(off_t) 0,SEEK_CUR);
					
					outfpsg1=fopen("matrixpre.txt","w+"); //mutiptest
					fseek(outfpsg1,(off_t) 0,SEEK_CUR);

					outfpsg2=fopen("matrixN1.txt","w+"); //g
					fseek(outfpsg2,(off_t) 0,SEEK_CUR);
					
					outfpsg3=fopen("matrixS.txt","w+"); 
					fseek(outfpsg3,(off_t) 0,SEEK_CUR);
					
					//outfpsg4=fopen("matrixG.txt","w+"); 
					//fseek(outfpsg4,(off_t) 0,SEEK_CUR);
					
					for(i=0;i<nx1;i++){
						for(j=0;j<nz1;j++){
							
							fprintf(outfpsg0,"%16.12f+i*%16.12f",oor[i][j].r,oor[i][j].i); 
							//fprintf(outfpsg1,"%16.12f+i*%16.12f",mutiptest[i][j].r,mutiptest[i][j].i);
							fprintf(outfpsg1,"%16.12f+i*%16.12f",M[i][j].r,M[i][j].i);							
							fprintf(outfpsg2,"%16.12f+i*%16.12f",g[i][j].r,g[i][j].i);
							fprintf(outfpsg3,"%16.12f+i*%16.12f",s[i][j].r,s[i][j].i);
							
							/*
							fprintf(outfpsg0,"%16.12f+i*%16.12f",or[i][j],0.0); 
							//fprintf(outfpsg1,"%16.12f+i*%16.12f",mutiptest[i][j].r,mutiptest[i][j].i);
							fprintf(outfpsg1,"%16.12f+i*%16.12f",M[i][j].r,M[i][j].i);							
							fprintf(outfpsg2,"%16.12f+i*%16.12f",g0[i][j].r,g0[i][j].i);
							fprintf(outfpsg3,"%16.12f+i*%16.12f",s0[i][j].r,s0[i][j].i);
							*/
						}
						fprintf(outfpsg0,"\n"); 
						fprintf(outfpsg1,"\n");   
						fprintf(outfpsg2,"\n");
						fprintf(outfpsg3,"\n");
					}
					
					/*
					for(ix=0;ix<nx1;ix++)
		        	{
						for(iz=0;iz<nz1;iz++)
						{	
                            				
			               for(ix1=0;ix1<nx1;ix1++)
			               {                        
			                  for(iz1=0;iz1<nz1;iz1++)
			                  {
							  	
								fprintf(outfpsg4,"%16.12f+i*%16.12f",green[ix*nz1+iz][ix1*nz1+iz1].r,green[ix*nz1+iz][ix1*nz1+iz1].i);                             																			                						                       
			                  }							
			               }//r'
			               fprintf(outfpsg4,"\n"); 			                                            
		                }                   
		        	}//r
					*/
				
		        	fclose(outfpsg0);
					fclose(outfpsg1);
					fclose(outfpsg2);	
					fclose(outfpsg3);
					//fclose(outfpsg4);
				}
				
				printf("here3\n");
				
				/*SSOR*/
				linssor_CW_PreO(g,diag,M,s,nx1,nz1,u0,itmax);	
				//linssor_CW_PreP(g,diag,s,nx1,nz1,u0,itmax);	
				/*
				for(i=0;i<nx1;i++){
					for(j=0;j<nz1;j++){
						u0[i][j]=u0[i][j];//cmul(u0[i][j],cinv(oor[i][j]));
					}
				}*/
				
				printf("here4\n");
				
				/*Osnabrugge preconditioner*/
				//linssor_preOCW(g,diag,M,s,nx1,nz1,u0,itmax);
				/* precondition matrix M is the main diagonal elements of the coefficient matrix */ 
				//linbcg_CW(diag,g,M,nx1,nz1,s,u0,itol,EPS,itmax,&iter,&err);								 				                            
				//linbcg_psc(diag,g,M,nx1,nz1,s,u0,itol,EPS,itmax,&iter,&err);
				//linbcgstab_psc(diag,g,M,nx1,nz1,s,u0,itol,EPS,itmax,&iter,&err);  
					   
				/* free space */
                free2complex(s);
                free2complex(g);
				//free2complex(gpre);//20201112 add
				free2complex(oor);
                free2complex(M);
				//free2complex(green);				
                //free2complex(mutiptest);
				
			 /************** iteration improvement***********************/
			 /*
			 complex **u00=NULL,**mutipu0=NULL;
			 complex uerrtemp;
			 
			 r0=alloc2complex(nz1,nx1);
			 d0=alloc2complex(nz1,nx1);
			 u00=alloc2complex(nz1,nx1);
			 mutipu0=alloc2complex(nz1,nx1);
			 u1=alloc2complex(nz1,nx1);
			 
				for(i=0;i<nx1;i++){
					for(j=0;j<nz1;j++){
						u00[i][j]=crmul(u0[i][j],or[i][j]);			
					}				 					 										 		
				}
	
				BTTB_y(g0,u00,nx1,nz1,mutipu0);//FFT2
				 
				for(i=0;i<nx1;i++){
					for(j=0;j<nz1;j++){
						r0[i][j]=cadd(csub(s0[i][j],u0[i][j]),mutipu0[i][j]);
						d0[i][j]=r0[i][j];
					}				 					 										 		
				}
				
				iteru=0;
								
				while(iteru<10){
					
					++(iteru);
					//************Ld=r0
					linssor(g0,or,r0,nx1,nz1,d0,50);
					
					for(i=0;i<nx1;i++){
						for(j=0;j<nz1;j++){
							u1[i][j]=cadd(u0[i][j],d0[i][j]);
							u00[i][j]=crmul(u1[i][j],or[i][j]);
						}				 					 										 		
					}
					
					BTTB_y(g0,u00,nx1,nz1,mutipu0);//FFT2
				 
					for(i=0;i<nx1;i++){
						for(j=0;j<nz1;j++){
							r0[i][j]=cadd(csub(s0[i][j],u1[i][j]),mutipu0[i][j]);
							d0[i][j]=r0[i][j];
						}				 					 										 		
					}
					
					uerr=0.0;
					uerrmax=0.0;
					uerrtemp=cmplx(0.0,0.0);
					
					
					for(i=0;i<nx1;i++){
						for(j=0;j<nz1;j++){
							uerrtemp=csub(u1[i][j],u0[i][j]);
							uerr =rcabs(uerrtemp);
							if(uerr>uerrmax) uerrmax=uerr;
						}				 					 										 		
					}
					
					//printf("iteru=%d, uerr=%lf\n",iteru,uerr);
					if(uerr<EPS2) break;
					
					for(i=0;i<nx1;i++){
						for(j=0;j<nz1;j++){
							u0[i][j]=u1[i][j];			
						}				 					 										 		
					}										
				}
				
				//********* free space 
				
                free2complex(s0);
                free2complex(g0);	
				
				free2complex(r0);
                free2complex(d0);
				free2complex(u1);
				free2complex(u00);
				free2complex(mutipu0);
				
				printf("iter impro finish\n");
				
			 /************************ Freqency domain wavefeid output **************************************/		     
                    if(w==fw){                                                     //(iw+offset)==(nw/4)
                    	outfpfrewf1=fopen("NySOR_frewf_re_SCpreo.txt","w+"); //frequency wavefeild real
       					fseek(outfpfrewf1,(off_t) 0,SEEK_CUR);

       					outfpfrewf2=fopen("NySOR_frewf_im_SCpreo.txt","w+"); //frequency wavefeild image
       					fseek(outfpfrewf2,(off_t) 0,SEEK_CUR);
       					
                        for(i=0;i<nx1;i++){
                            for(j=0;j<nz1;j++){
                                 fprintf(outfpfrewf1,"%16.8f",u0[i][j].r);
                                 fprintf(outfpfrewf2,"%16.8f",u0[i][j].i);								                               
                            }
                            fprintf(outfpfrewf1,"\n"); 
                            fprintf(outfpfrewf2,"\n"); 
                        }
                        
                       fclose(outfpfrewf1);
            		   fclose(outfpfrewf2);
                    }    
					 					 
			    /* loop over offsets */
		        for(ixo=0; ixo<nxo; ixo++)
	    	    {
		            xg=fxo+ixo*dxo;
                   
                    /* allocate space */

	                sg=alloc2complex(nz1,nx1);	                
	            
	            	/* compute Gmj(w,r,rg)=h*h*Gmj */				      
			        /*#ifdef _OPENMP			        
                    #pragma omp parallel for default(none)  num_threads (8) \
                    private(ix,iz) \
                    shared(nx1,nz1,xg,zg,fx,fz,h,k0,js2,sg)
                    #endif	*/	        
                    for(ix=0;ix<nx1;ix++)
		            {
			            for(iz=0;iz<nz1;iz++)
		    	        {
			        	     //sg[ix][iz].r=k0*k0*gaus2_CW(xg,zg,fx,fz,h,ix,iz,kk0,2,js2,ss2_CW,fre2_CW);
				             //sg[ix][iz].i=k0*k0*gaus2_CW(xg,zg,fx,fz,h,ix,iz,kk0,2,js2,ss2_CW,fim2_CW);
				             sg[ix][iz].r=window[ix][iz]*k0*k0*gaus2(xg,zg,fx,fz,h,ix,iz,k0,2,js2,ss2,fre2);
				             sg[ix][iz].i=window[ix][iz]*k0*k0*gaus2(xg,zg,fx,fz,h,ix,iz,k0,2,js2,ss2,fim2);							 	
		                }
					}
		         /* compute D(r) */
			
		            dr.r=0.0;
                    dr.i=0.0;
			        /*#ifdef _OPENMP			        
                    #pragma omp parallel for default(none)  num_threads (t) \
                    private(i,j) \
                    shared(nx1,nz1,sg,u0) \
                    reduction(+:dr.r,dr.i)
                    #endif	*/	       
		            for(i=0;i<nx1;i++)
		            {
                        for(j=0;j<nz1;j++)
                        {
			                 dr=cadd(dr,cmul(crmul(sg[i][j],or[i][j]),u0[i][j]));
							// dr=cadd(dr,cmul(cmul(sg[i][j],diag[i][j]),u0[i][j]));
                        }			
		            }

                    //umtemp[iw][ixo].r=dr.r;   
                    //umtemp[iw][ixo].i=dr.i;
					umtemp[iw][ixo]=dr;
					//umtemp[iw][ixo]=cmul(ct[iw+offset],dr);
										
                    free2complex(sg);
          
	            }//offset over
				
					if(w==fw){
						outfpsg5=fopen("NySOR_Wr_SCpreo.txt","w+"); //reciever
						fseek(outfpsg5,(off_t) 0,SEEK_CUR);

						outfpsg6=fopen("NySOR_Wi_SCpreo.txt","w+"); //reciever
						fseek(outfpsg6,(off_t) 0,SEEK_CUR);
						
						
							for(ixo=0;ixo<nxo;ixo++){
								fprintf(outfpsg5,"%16.8f",umtemp[iw][ixo].r);     
								fprintf(outfpsg6,"%16.8f",umtemp[iw][ixo].i);
							}
	
						fclose(outfpsg5);
						fclose(outfpsg6);
					}
                    
                printf("Singel Gsg finish\n");
                /* free space */
                free2complex(u0);
                free2complex(diag);
				
                
                //}
	        }//iw over
	       
	        MPI_Send(&offset, 1, MPI_INT, 0, FROM_SLAVE, MPI_COMM_WORLD); 
			MPI_Send(&rows_w, 1, MPI_INT, 0, FROM_SLAVE, MPI_COMM_WORLD); 
			MPI_Send(&umtemp[0][0], rows_w*nxo,COMPLEX, 0, FROM_SLAVE, MPI_COMM_WORLD); 
			
			free2complex(umtemp);
			//MPI_Barrier(MPI_COMM_WORLD);
		}
		
		free1complex(ct);
        free2double(or);
	    free2double(window);
	    fclose(vfp);
			
	   	/*clean up the MPI_Type COMPLEX*/
		MPI_Type_free(&COMPLEX);  
       
		MPI_Finalize(); 		
		
	return 0;
}
