#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

#define PARAMETERFILE   "parameters.dat"
#define PLOTFILE        "plot_"
#define SAMPLESIZE      "sample_sizes.dat"
#define INFILE          "infile.dat"
#define PREFIX          "Pair_"
#define SUFFIX          ".dat"
#define ST              0

struct DATA {
	int L;                                                                         // Number of loci
	int ns;                                                                        // Number of sampled subpopoulations
	int *k;                                                                        // Number of alleles (per locus)
	int **T;                                                                       // Sum of allele counts (per locus and per subpopulation)
	int ***n;                                                                      // Allele counts (per locus, per subpopulation, and per allele)
};

struct STAT {
	double **F[3];
	double *multilocus_F[3];  
	double **maf;
	double **H;
	int **K;
};

void GetData() {
	
	void ReadParameterFileHeader(char datafile_name[32],int *dominance,double *MAF,double *A,double *B);
	int ReadData(char datafile_name[32],struct DATA *D);
	int Estimate(struct DATA D,struct STAT *S,double MAF,double A,double B,int dominance);
	void AllocateMemory(struct DATA D,struct STAT *S);
	void ReleaseMemoryData(struct DATA *D);
	void ReleaseMemoryStatistics(struct DATA D,struct STAT *S);
	void WriteOutputs(struct DATA D,struct STAT S,double MAF);	
	
	struct DATA D;
	struct STAT S;
	double MAF,A,B;
	int dominance;
	char datafile_name[32];
	
	ReadParameterFileHeader(datafile_name,&dominance,&MAF,&A,&B);
	if (ReadData(datafile_name,&D) == 0) {
		AllocateMemory(D,&S);
		if (Estimate(D,&S,MAF,A,B,dominance) == 0) {
			WriteOutputs(D,S,MAF);
		}
		ReleaseMemoryStatistics(D,&S);
	}
	ReleaseMemoryData(&D);
}

void ReadParameterFileHeader(char datafile_name[32],
							 int *dominance,
							 double *lMAF,
							 double *A,
							 double *B)
{
	FILE *parameters;	
	int X;
	
	parameters = fopen(PARAMETERFILE,"r");
	fscanf(parameters,"%s",datafile_name);
	while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	fscanf(parameters,"%d",dominance);
	while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	fscanf(parameters,"%lf",lMAF);
	while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	if (*dominance) {
		fscanf(parameters,"%lf",A);
		while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
		fscanf(parameters,"%lf",B);
	}
	fclose(parameters);
}	


int ReadData(char datafile_name[32],struct DATA *D)

{
	FILE *infile;
	int X;
	int ByAlleles;
	int i,j,k,l,m;
	int dummy,error;
	
	infile = fopen(datafile_name,"r");
	if ((error = fscanf(infile,"%d",&ByAlleles)) == 0 || error == EOF) {           // 1 is by allele 0 is by population
		Rprintf("STOPPED: problem reading file...\n");
		error = 1;
		goto end;
	}
	if (ByAlleles > 1) {
		Rprintf("STOPPED: problem reading file...\n");
		error = 1;
		goto end;
	}
	while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
	if ((error = fscanf(infile,"%d",&D -> ns)) == 0 || error == EOF) {             // Read the number of sampled demes 
		Rprintf("STOPPED: problem reading file...\n");
		error = 1;
		goto end;
	}
	while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
	if ((error = fscanf(infile,"%d",&D -> L)) == 0 || error == EOF) {              // Read the number of loci 
		Rprintf("STOPPED: problem reading file...\n");
		error = 1;
		goto end;
	}
	while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
	D -> k = (int *) malloc(D -> L * sizeof(int));                                 // Allocate memory for the number of alleles, the allele coounts and the sum of allele counts
	D -> T = (int **) malloc(D -> L * sizeof(int *));
	D -> n = (int ***) malloc(D -> L * sizeof(int **));
	for(i = 0;i < D -> L; ++i) {
		D -> T[i]=(int *) malloc(D -> ns * sizeof(int));
		D -> n[i]=(int **) malloc(D -> ns * sizeof(int *));
	}
	for(i = 0; i < D -> L; ++i) {                                                  // Loop over all loci in the dataset
		if ((error = fscanf(infile,"%d",&D -> k[i])) == 0 || error == EOF) {       // Read the number of alleles at the ith locus 
			Rprintf("STOPPED: problem reading file...\n");
			for(l = 0; l < D -> L; ++l) {
				for(m = 0; m < D -> ns; ++m) {
					D -> n[l][m] = NULL;
				}				
			}
			error = 1;
			goto end;
		}
		while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
		for(j = 0; j < D -> ns; ++j) {                                             // Allocate memory for the allele counts
			D -> n[i][j] = (int *) malloc(D -> k[i] * sizeof(int));
		}
		if (ByAlleles) {
			for(j = 0; j < D -> k[i]; ++j) {                                       // Loop over all alleles for the ith locus
				for(k = 0; k < D -> ns;++k) {                                      // Loop over all populations for the ith locus and the jth allele
					if ((error = fscanf(infile,"%d",&D -> n[i][k][j])) == 0 || error == EOF) { // Read allele count for the ith locus, jth population and kth allele
						Rprintf("STOPPED: problem reading file...\n");
						for(l = 0; l < D -> L; ++l) {
							for(m = 0; m < D -> ns; ++m) {
								D -> n[l][m] = NULL;
							}				
						}						
						error = 1;
						goto end;
					}
				}                                                                  // End of the loop over populations
			}                                                                      // End of the loop over alleles
		}
		else {
			for(j = 0; j < D -> ns; ++j) {                                         // Loop over all populations for the ith locus
				for(k = 0; k < D -> k[i]; ++k) {                                   // Loop over all alleles for the ith locus and the jth population
					if ((error = fscanf(infile,"%d",&D -> n[i][j][k])) == 0 || error == EOF) { // Read allele count for the ith locus, jth population and kth allele
						Rprintf("STOPPED: problem reading file...\n");
						for(l = 0; l < D -> L; ++l) {
							for(m = 0; m < D -> ns; ++m) {
								D -> n[l][m] = NULL;
							}				
						}
						error = 1;
						goto end;
					}
				}                                                                  // End of the loop over alleles
			}                                                                      // End of the loop over populations
		}
		for(j = 0; j < D -> ns; ++j) {                                             // Loop over all populations for the ith locus
			D -> T[i][j] = 0;                                                      // Initialize the allele counts
			for(k = 0; k < D -> k[i]; ++k) {                                       // Loop over all alleles for the ith locus and the jth populations
				D -> T[i][j] += D -> n[i][j][k];                                   // Calculate the sum of allele counts, per locus and per population
			}                                                                      // End of the loop over alleles
		}                                                                          // End of the loop over populations
	}                                                                              // End of the loop over loci	
	if ((error = fscanf(infile,"%d",&dummy)) != 0 && error != EOF) {
		Rprintf("STOPPED: problem reading file...\n");
		error = 1;
		goto end;
	}
	error = 0;                                                                     // Normal ending
end :
	fclose(infile);
	return(error);
}

int Estimate(struct DATA D,
             struct STAT *S,
			 double MAF,
			 double A,
			 double B,
			 int dominance)

{
	int i,j,k,l,u;
	int cpt,error;
	double sum_num1,sum_num2,sum_den;
	double num,num_1,num_2,den,n_c;
	double yy,q2,q3;
	double var = 0.0;
	double **p;
	
	double nsum,nsum2,nbar,xx;
	double sum_num_WC,sum_den_WC;
	
	cpt = 0;
	for (i = 0; i < (D.ns - 1); i++) {
		for (j = i + 1; j < D.ns; j++) {
			sum_num1 = sum_num2 = sum_den = 0;
			sum_num_WC = sum_den_WC = 0;
			for (l = 0; l < D.L; l++) {
				S -> K[cpt][l] = 0;
				p = (double **) malloc(3 * sizeof (double *));
				for (u = 0; u < 3; ++u) {
					p[u] = (double *) malloc(D.k[l] * sizeof (double));
				}
				S -> maf[cpt][l] = -2147483647L;
				for (u = 0; u < D.k[l]; ++u) {
					if ((D.n[l][i][u] + D.n[l][j][u]) != 0) S -> K[cpt][l] += 1;
					p[2][u] = (double) (D.n[l][i][u] + D.n[l][j][u]) / (D.T[l][i] + D.T[l][j]);
					if (p[2][u] > S -> maf[cpt][l]) {
						S -> maf[cpt][l] = p[2][u];
					}
				}
				if ((dominance) & (S -> K[cpt][l] > 2)) {
					Rprintf("STOPPED: one AFLP marker has more than 2 alleles...\n");
					for (k = 0; k < 3; ++k) {
						free(p[k]);
					}
					free(p);
					error = 1;
					goto end;
				}
				if (S -> maf[cpt][l] <= MAF) {
					if (dominance) {
						p[0][0] = exp(lgammafn(D.n[l][i][0] + A + 0.5) + lgammafn(D.T[l][i] + A + B) - lgammafn(D.n[l][i][0] + A) - lgammafn(D.T[l][i] + A + B + 0.5)); // Zhivotovsky 1999 eq. (10)
						p[1][0] = exp(lgammafn(D.n[l][j][0] + A + 0.5) + lgammafn(D.T[l][j] + A + B) - lgammafn(D.n[l][j][0] + A) - lgammafn(D.T[l][j] + A + B + 0.5));
						p[2][0] = exp(lgammafn(D.n[l][i][0] + D.n[l][j][0] + A + 0.5) + lgammafn(D.T[l][i] + D.T[l][j] + A + B) - lgammafn(D.n[l][i][0] + D.n[l][j][0] + A) - lgammafn(D.T[l][i] + D.T[l][j] + A + B + 0.5));
						for (k = 0; k < 3; ++k) {
							p[k][1] = 1.0 - p[k][0];
						}
						var = exp(lgammafn(D.n[l][i][0] + D.n[l][j][0] + A + 1.0) + lgammafn(D.T[l][i] + D.T[l][j] + A + B) - lgammafn(D.n[l][i][0] + D.n[l][j][0] + A) - lgammafn(D.T[l][i] + D.T[l][j] + A + B + 1.0)) - p[2][0] * p[2][0];
					} else {
						for (u = 0; u < D.k[l]; ++u) {
							p[0][u] = (double) D.n[l][i][u] / D.T[l][i];
							p[1][u] = (double) D.n[l][j][u] / D.T[l][j];
						}
					}
					num = num_1 = num_2 = den = 0;
					nsum = (double) D.T[l][i] + D.T[l][j];
					nsum2 = pow(D.T[l][i],2) + pow(D.T[l][j],2);
					n_c = nsum - nsum2 / nsum;
					for (u = 0; u < D.k[l]; ++u) {
						num_1 += p[0][u] * (1 - p[0][u]);
						num_2 += p[1][u] * (1 - p[1][u]);
						num += (double) pow(D.T[l][i],2) / (D.T[l][i] - 1) * p[0][u] * (1 - p[0][u]);
						num += (double) pow(D.T[l][j],2) / (D.T[l][j] - 1) * p[1][u] * (1 - p[1][u]);
						den += (double) D.T[l][i] * pow((p[0][u] - p[2][u]),2) + (D.T[l][i] - pow(D.T[l][i],2) / (D.T[l][i] + D.T[l][j])) * p[0][u] * (1 - p[0][u]);
						den += (double) D.T[l][j] * pow((p[1][u] - p[2][u]),2) + (D.T[l][j] - pow(D.T[l][j],2) / (D.T[l][i] + D.T[l][j])) * p[1][u] * (1 - p[1][u]);
					}
					num_1 *= n_c * D.T[l][i] / (D.T[l][i] - 1);
					num_2 *= n_c * D.T[l][j] / (D.T[l][j] - 1);
					num *= n_c;
					if (den > 0) S -> F[1][cpt][l] = 1 - (num_1 / den); else S -> F[1][cpt][l] = 9.0;
					if (den > 0) S -> F[2][cpt][l] = 1 - (num_2 / den); else S -> F[2][cpt][l] = 9.0;
					sum_num1 += num_1;
					sum_num2 += num_2;
					sum_den += den;
					if (dominance) {
						S -> H[cpt][l] = 2 * p[2][0] * (1 - p[2][0]) + 2 * var;
					} else {
						S -> H[cpt][l] = 0;
						for (u = 0; u < D.k[l]; ++u) {
							S -> H[cpt][l] += p[2][u] * (1 - p[2][u]);
						}
						S -> H[cpt][l] *= (double) (D.T[l][i] + D.T[l][j]) / (D.T[l][i] + D.T[l][j] - 1);
					}
					nbar = nsum / 2; // FROM fdist2's function thetacal !!!! IDEM WEIR 1996 (5.8)-(5.9)
					xx = yy = 0.0;
					for (u = 0; u < D.k[l]; ++u) {
						xx += D.T[l][i] * pow(p[0][u],2);
						xx += D.T[l][j] * pow(p[1][u],2);
						yy += pow((p[0][u] * D.T[l][i] + p[1][u] * D.T[l][j]),2);
					}
					q2 = (xx - 2) / (2 * (nbar - 1.0));
					q3 = 1.0 / (2 * nbar * n_c) * (yy - nbar * (n_c - 1.0) / (nbar - 1.0) * xx) + (nbar - n_c) / (n_c * (nbar - 1.0)) * (1.0 - xx);
					sum_num_WC += (1 - q2);
					sum_den_WC += (1 - q3);
					if ((1 - q3) > 0) S -> F[ST][cpt][l] = 1.0 - (1 - q2) / (1 - q3); else S -> F[ST][cpt][l] = 9.0;
				}				
				for (k = 0; k < 3; ++k) {
					free(p[k]);
				}
				free(p);
			}
			S -> multilocus_F[1][cpt] = 1.0 - sum_num1 / sum_den;
			S -> multilocus_F[2][cpt] = 1.0 - sum_num2 / sum_den;
			S -> multilocus_F[ST][cpt] = 1.0 - sum_num_WC / sum_den_WC;
			cpt += 1;
		}
	}
	error = 0;
end:
	return(error);
}

void WriteOutputs(struct DATA D,
                  struct STAT S,
                  double MAF)

{
	int i,j,l;
	int cpt;
	char OutFileName[32];
	char PlotFileName[32];
	char x[3];
	FILE *outfile,*plotfile,*samplesize;
	
	samplesize = fopen(SAMPLESIZE,"w");
	fprintf(samplesize,"%6d\n",D.L);
	fprintf(samplesize,"%6d\n",D.ns);
	for (l = 0; l < D.L; l++) {
		fprintf(samplesize,"%6d",(l + 1));
		for (i = 0; i < D.ns; i++) {
			fprintf(samplesize,"%6d",D.T[l][i]);
		}
		fprintf(samplesize,"\n");
	}
	fclose(samplesize);
	outfile = fopen(INFILE,"w");
	cpt = 0;
	fprintf(outfile,"%s\t %3s %8s %15s %15s %15s\n","Filename","i","j","F_1","F_2","F_ST(WC)");
	for (i = 0; i < (D.ns - 1); i++) {
		for (j = i + 1; j < D.ns; j++) {
			strcpy(PlotFileName,PLOTFILE);
			sprintf(x,"%d",i + 1);
			strcat(PlotFileName,x);
			strcat(PlotFileName,"_");
			sprintf(x,"%d",j + 1);
			strcat(PlotFileName,x);
			strcat(PlotFileName,SUFFIX);
			plotfile = fopen(PlotFileName,"w");
			for (l = 0; l < D.L; l++) {
				if (S.maf[cpt][l] <= MAF) {
					fprintf(plotfile,"%15.9f %15.9f %15.9f %15.9f %15d %15d\n",S.F[1][cpt][l],S.F[2][cpt][l],S.F[ST][cpt][l],S.H[cpt][l],S.K[cpt][l],(l + 1));
				}
			}
			strcpy(OutFileName,PREFIX);
			sprintf(x,"%d",i + 1);
			strcat(OutFileName,x);
			strcat(OutFileName,"_");
			sprintf(x,"%d",j + 1);
			strcat(OutFileName,x);
			fprintf(outfile,OutFileName,"\t");
			fprintf(outfile,"%8d %8d %15.9f %15.9f %15.9f \n",(i + 1),(j + 1),S.multilocus_F[1][cpt],S.multilocus_F[2][cpt],S.multilocus_F[ST][cpt]);
			++cpt;
			fclose(plotfile);
		}
	}
	fclose(outfile);
}

void AllocateMemory(struct DATA D,
                    struct STAT *S)

{
	int i,j;
	int cpt;
	
	cpt = (int) (D.ns * (D.ns - 1) / 2);
	for(j = 0; j < 3; ++j) {
		S -> multilocus_F[j] = (double *) malloc (cpt * sizeof(double));
	}
	for(j = 0; j < 3; ++j) {
		S -> F[j] = (double **) malloc (cpt * sizeof(double *));
	}
	S -> H = (double **) malloc (cpt * sizeof(double *));
	S -> maf = (double **) malloc (cpt * sizeof(double *));
	S -> K = (int **) malloc (cpt * sizeof(int *));
	for(i = 0; i < cpt; ++i) {
		for(j = 0; j < 3; ++j) {
			S -> F[j][i] = (double *) malloc (D.L * sizeof(double));
		}
		S -> H[i] = (double *) malloc (D.L * sizeof(double));
		S -> maf[i] = (double *) malloc (D.L * sizeof(double));
		S -> K[i] = (int *) malloc (D.L * sizeof(int));
	}
}

void ReleaseMemoryData(struct DATA *D)

{
	int i,j;
	
	if ((D -> k) != NULL) {free(D -> k);}
	if ((D -> T) != NULL) {	
		for(i = 0; i < D -> L; ++i) {
			free(D -> T[i]);
		}
		free(D -> T);
	}
	if ((D -> n) != NULL) {	
		for(i = 0; i < D -> L; ++i) {
			for(j = 0; j < D -> ns; ++j) {
				if ((D -> n[i][j]) != NULL) {free(D -> n[i][j]);}
			}
		}
		for(i = 0; i < D -> L; ++i) {
			free(D -> n[i]);
		}
		free(D -> n);
	}	
}

void ReleaseMemoryStatistics(struct DATA D,
							 struct STAT *S)

{
	int i,j;
	int cpt;
	
	cpt = (int) (D.ns * (D.ns - 1) / 2);
	for(j = 0; j < 3; ++j) {
		free(S -> multilocus_F[j]);
	}
	for(i = 0; i < cpt; ++i) {
		for(j = 0; j < 3; ++j) {
			free(S -> F[j][i]);
		}
		free(S -> H[i]);
		free(S -> maf[i]);
		free(S -> K[i]);
	}
	for(j = 0; j < 3; ++j) {
		free(S -> F[j]);
	}
	free(S -> K);
	free(S -> H);
	free(S -> maf);
}
