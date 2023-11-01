
double GreenCWr(dcomplex r)
{
     double gr;
	 dcomplex h;
	 h=hankelH1(0,r);
     gr=(-0.25)*h.i;
	 return gr;
}
double GreenCWi(dcomplex r)
{
     double gi;
	 dcomplex h;
	 h=hankelH1(0,r);
     gi=(0.25)*h.r;
     return gi;
}

