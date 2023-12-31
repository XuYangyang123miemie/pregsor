
double Greenr(double r)
{
     double gr;
     gr=(-0.25)*Hanker_first_zero_imag( r );
	 return gr;
}
double Greeni(double r)
{
     double gi;
     gi=(0.25)*Hanker_first_zero_real( r );
     return gi;
}

// ****** 贝塞尔 Bessel 函数 ******

// ****** 第 一 类 ******

// ****** 第 零 阶 ****** 计算 此函数是个子函数 需要通过主函数调用


double Bessel_first_zero(double x)
{
	//-3<x<3时，多项式系数
	double a1[6] = {-2.2499997, +1.2656208, -0.3163866, +0.0444479, -0.0039444, +0.0002100};

	//x>3时，f0系数
	double a2[7] = {+0.79788456, -0.00000077, -0.00552740, -0.00009512, +0.00137237, -0.00072805, +0.00014476};
	
	//x>3时，s0系数
	double a3[7] = {-0.78539816, -0.04166397, -0.00003954, +0.00262573, -0.00054125, -0.00029333, +0.00013558};
	
	//x的某次幂的计算结果
	double b1[6], b2[7];

	double j0, x1, x2, f0, s0;

	
	if ( fabs(x) <= 3 )

	{
		x1 = x / 3 ;
		b1[0] = x1 * x1; //2次幂
		b1[1] = b1[0] * b1[0]; //4次幂
		b1[2] = b1[1] * b1[0]; //6次幂
		b1[3] = b1[1] * b1[1]; //8次幂
		b1[4] = b1[2] * b1[1]; //10次幂;'
		b1[5] = b1[2] * b1[2]; //12次幂
		j0 = 1 + a1[0]*b1[0] + a1[1]*b1[1] + a1[2]*b1[2] + a1[3]* b1[3] + a1[4]*b1[4] + a1[5]*b1[5];
	}

	else

	{
		x2 = 3 / x ;
		b2[0] = x2 * x2 ;        //2次幂
		b2[1] = b2[0] * x2 ;     //3次幂
		b2[2] = b2[0] * b2[0] ;  //4次幂
		b2[3] = b2[2] * x2 ;     //5次幂
		b2[4] = b2[1] * b2[1];   //6次幂
		f0 = a2[0] + a2[1]*x2 + a2[2]*b2[0] + a2[3]*b2[1] + a2[4]*b2[2] + a2[5]*b2[3] + a2[6]*b2[4];
		s0 = x + a3[0] + a3[1]*x2 + a3[2]*b2[0] + a3[3]*b2[1] + a3[4]*b2[2] + a3[5]*b2[3] + a3[6]*b2[4];
		j0 = pow(x,(-0.5)) * f0 * cos(s0) ;
	}
	
	return  j0;
}





// ****** 贝塞尔 Bessel 函数 ******

// ****** 第 一 类 ******

// ****** 第 一 阶 ****** 计算 此函数是个子函数 需要通过主函数调用


double Bessel_first_first(double x)
{
	//-3<x<3时，多项式系数
	double a1[6] = {-0.56249985, +0.21093573, -0.03954289, +0.00443319, -0.00031761, +0.00001109};

	//x>3时，f0系数
	double a2[7] = {+0.79788456, +0.00000156, +0.01659667, +0.00017105, -0.00249511, +0.00113653, -0.00020033};
	
	//x>3时，s0系数
	double a3[7] = {-2.35619449, +0.12499612, +0.00005650, -0.00637879, +0.00074348, -0.00079824, -0.00029166};
	
	//x的某次幂的计算结果
	double b1[6], b2[7];

	double j1, j11, x1, x2, f1, s1;

	
	if ( fabs(x) <= 3 )

	{
		x1 = x / 3 ;
		b1[0] = x1 * x1; //2次幂
		b1[1] = b1[0] * b1[0]; //4次幂
		b1[2] = b1[1] * b1[0]; //6次幂
		b1[3] = b1[1] * b1[1]; //8次幂
		b1[4] = b1[2] * b1[1]; //10次幂
		b1[5] = b1[2] * b1[2]; //12次幂
		j11 = 0.5 + a1[0]*b1[0] + a1[1]*b1[1] + a1[2]*b1[2] + a1[3]* b1[3] + a1[4]*b1[4] + a1[5]*b1[5];
		j1 = j11 * x ;
	}

	else

	{
		x2 = 3 / x ;
		b2[0] = x2 * x2 ;        //2次幂
		b2[1] = b2[0] * x2 ;     //3次幂
		b2[2] = b2[0] * b2[0] ;  //4次幂
		b2[3] = b2[2] * x2 ;     //5次幂
		b2[4] = b2[1] * b2[1];   //6次幂
		f1 = a2[0] + a2[1]*x2 + a2[2]*b2[0] + a2[3]*b2[1] + a2[4]*b2[2] + a2[5]*b2[3] + a2[6]*b2[4];
		s1 = x + a3[0] + a3[1]*x2 + a3[2]*b2[0] + a3[3]*b2[1] + a3[4]*b2[2] + a3[5]*b2[3] + a3[6]*b2[4];
		j1 = pow(x, (-0.5)) * f1 * cos( s1 ) ;
	}
	
	return  j1;
}





// ****** 贝塞尔 Bessel 函数 ******

// ****** 第 二 类 ******

// ****** 第 零 阶 ****** 计算 此函数是个子函数 需要通过主函数调用


double Bessel_second_zero(double x)
{
	//-3<x<3时，多项式系数
	double a1[7] = {+0.36746691, +0.60559366, -0.74350384, +0.25300117, -0.04261214, +0.00427916, -0.00024846};

	//x>3时，f0系数
	double a2[7] = {+0.79788456, -0.00000077, -0.00552740, -0.00009512, +0.00137237, -0.00072805, +0.00014476};
	
	//x>3时，s0系数
	double a3[7] = {-0.78539816, -0.04166397, -0.00003954, +0.00262573, -0.00054125, -0.00029333, +0.00013558};
	
	//x的某次幂的计算结果
	double b1[6], b2[7];

	double y0, y01, x1, x2, f0, s0;

	
	if ( fabs(x) <= 3 )

	{
		x1 = x / 3 ;
		b1[0] = x1 * x1; //2次幂
		b1[1] = b1[0] * b1[0]; //4次幂
		b1[2] = b1[1] * b1[0]; //6次幂
		b1[3] = b1[1] * b1[1]; //8次幂
		b1[4] = b1[2] * b1[1]; //10次幂
		b1[5] = b1[2] * b1[2]; //12次幂
		y01 = ( 2 / 3.1415926535897932384626 ) * log( 0.5 * x ) * Bessel_first_zero( x );
		y0 = y01 + a1[0] + a1[1]*b1[0] + a1[2]*b1[1] + a1[3]*b1[2] + a1[4]* b1[3] + a1[5]*b1[4] + a1[6]*b1[5];
	}

	else

	{
		x2 = 3 / x ;
		b2[0] = x2 * x2 ;        //2次幂
		b2[1] = b2[0] * x2 ;     //3次幂
		b2[2] = b2[0] * b2[0] ;  //4次幂
		b2[3] = b2[2] * x2 ;     //5次幂
		b2[4] = b2[1] * b2[1];   //6次幂
		f0 = a2[0] + a2[1]*x2 + a2[2]*b2[0] + a2[3]*b2[1] + a2[4]*b2[2] + a2[5]*b2[3] + a2[6]*b2[4];
		s0 = x + a3[0] + a3[1]*x2 + a3[2]*b2[0] + a3[3]*b2[1] + a3[4]*b2[2] + a3[5]*b2[3] + a3[6]*b2[4];
		y0 = pow(x, (-0.5)) * f0 * sin( s0 ) ;
	}
	
	return  y0;
}





// ****** 贝塞尔 Bessel 函数 ******

// ****** 第 二 类 ****** 

// ****** 第 一 阶 ****** 计算 此函数是个子函数 需要通过主函数调用


double Bessel_second_first(double x)
{
	//-3<x<3时，多项式系数
	double a1[7] = {-0.6366198, +0.2212091, +2.1682709, -1.3164827, +0.3123951, -0.0400976, +0.0027873};

	//x>3时，f0系数
	double a2[7] = {+0.79788456, -0.00000077, -0.00552740, -0.00009512, +0.00137237, -0.00072805, +0.00014476};
	
	//x>3时，s0系数
	double a3[7] = {-0.78539816, -0.04166397, -0.00003954, +0.00262573, -0.00054125, -0.00029333, +0.00013558};
	
	//x的某次幂的计算结果
	double b1[6], b2[7];

	double y1, y11, y12, x1, x2, f1, s1;

	
	if ( fabs(x) <= 3 )

	{
		x1 = x / 3 ;
		b1[0] = x1 * x1; //2次幂
		b1[1] = b1[0] * b1[0]; //4次幂
		b1[2] = b1[1] * b1[0]; //6次幂
		b1[3] = b1[1] * b1[1]; //8次幂
		b1[4] = b1[2] * b1[1]; //10次幂
		b1[5] = b1[2] * b1[2]; //12次幂
		y12 = ( 2 / 3.1415926535897932384626 ) * log( 0.5 * x ) * Bessel_first_first( x );
		y11 = y12 + a1[0] + a1[1]*b1[0] + a1[2]*b1[1] + a1[3]*b1[2] + a1[4]* b1[3] + a1[5]*b1[4] + a1[6]*b1[5];
		y1 = y11 / x;
	}

	else

	{
		x2 = 3 / x ;
		b2[0] = x2 * x2 ;        //2次幂
		b2[1] = b2[0] * x2 ;     //3次幂
		b2[2] = b2[0] * b2[0] ;  //4次幂
		b2[3] = b2[2] * x2 ;     //5次幂
		b2[4] = b2[1] * b2[1];   //6次幂
		f1 = a2[0] + a2[1]*x2 + a2[2]*b2[0] + a2[3]*b2[1] + a2[4]*b2[2] + a2[5]*b2[3] + a2[6]*b2[4];
		s1 = x + a3[0] + a3[1]*x2 + a3[2]*b2[0] + a3[3]*b2[1] + a3[4]*b2[2] + a3[5]*b2[3] + a3[6]*b2[4];
		y1 = pow(x, (-0.5)) * f1 * sin( s1 ) ;
	}
	
	return  y1;
}





// ****** 汉克尔 Hanker 函数******

// ****** 第 一 类 ******

// ****** 第 零 阶 ******

// ****** 实 部 ******


double Hanker_first_zero_real(double x)
{
	double Bessel_first_zero(double x);
	
	double Hr;
	Hr = Bessel_first_zero( x );
	return Hr;
}



                                                                     

// ****** 汉克尔 Hanker 函数******

// ****** 第 一 类 ******

// ****** 第 零 阶 ******
                                                                                                                                                                                                                        
// ****** 虚 部 ******


double Hanker_first_zero_imag(double x)
{
	double Bessel_second_zero(double x);
	
	double Hi;
	Hi = Bessel_second_zero( x );
	return Hi;
}
                            
