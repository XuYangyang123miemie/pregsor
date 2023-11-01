/* Prime Factor FFTs (double version)*/
int npfa_d (int nmin);
int npfao_d (int nmin, int nmax);
int npfar_d (int nmin);
int npfaro_d (int nmin, int nmax);
void pfacc_d (int isign, int n, dcomplex z[]);
void pfacr_d (int isign, int n, dcomplex cz[], double rz[]);
void pfarc_d (int isign, int n, double rz[], dcomplex cz[]);
void pfamcc_d (int isign, int n, int nt, int k, int kt, dcomplex z[]);
void pfa2cc_d (int isign, int idim, int n1, int n2, dcomplex z[]);
void pfa2cr_d (int isign, int idim, int n1, int n2, dcomplex cz[],
double rz[]);
void pfa2rc_d (int isign, int idim, int n1, int n2, double rz[],
dcomplex cz[]);