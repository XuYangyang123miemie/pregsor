#ifndef PFAFFT_H
#define PFAFFT_H

/* Prime Factor FFTs */
int npfa (int nmin);
int npfao (int nmin, int nmax);
int npfar (int nmin);
int npfaro (int nmin, int nmax);
void pfacc (int isign, int n, complex z[]);
void pfarc (int isign, int n, double rz[], complex cz[]);
void pfacr (int isign, int n, complex cz[], double rz[]);
void pfa2cc (int isign, int idim, int n1, int n2, complex z[]);
void pfa2rc (int isign, int idim, int n1, int n2, double rz[], complex cz[]);
void pfa2cr (int isign, int idim, int n1, int n2, complex cz[], double rz[]);
void pfamcc (int isign, int n, int nt, int k, int kt, complex z[]);

#endif
