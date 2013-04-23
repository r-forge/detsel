#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

#define PARAMETERFILE   "parameters.dat"
#define INFILE          "infile.dat"
#define OUT             "out.dat"
#define SUFFIX          ".dat"
#define ST              0
#define TEST            10000

typedef struct {
	double F[3];
	double num_1;
	double num_2;
	double den;
	double num_FST;
	double den_FST;
	double H;
	double maf;
	int K[3];
} Stats;

typedef struct {
	double F1;
	double F2;
	int n1;
	int n2;
} Data;

typedef struct {
	double *Ne;
	double *N0;
	double N1;
	double N2;
	double mubar;
	double mu;
	double *t0;
	double *t;
} Parameters;

struct Node {
	double time;
	int allele;
	struct Node *descendant_1;
	struct Node *descendant_2;
	struct Node *ancestor;
};

struct Tree {
	int n;
	int nbr_ancestors;
	double begin;
	double end;
	double N;
	double m;
	struct Node *tree;
	struct Node **list;
};

int k_iam;
struct Tree deme[4];

void SimulDiv(int *p1,int *p2,int *n1,int *n2) {
	void Make_Samples(Parameters P,Data D,int set,int dominance,int MUTMOD);
	void NewStatistics(Stats *S,Data D,double MAF,double A,double B,int dominance);	
	void ReadParameterFile(char filename[32],int *dominance,double *MAF,double *A,double *B,int *totiter,Parameters *P,int *MUTMOD,int *nsets);

	char datafile_name[32];
	char filename[32];
	int i,j,iter,totiter,set,nsets,polymorphic,polymorphicsets;
	int *skip = NULL;
	FILE *outfile,*infile,*out;
	char y[8];			// BUG 17-09-2012 (was y[3])
    int X;
	Stats S;
	Data D;
	Parameters P;
	double dummy;
	double sum_num_1,sum_num_2,sum_den,sum_num_FST,sum_den_FST;
	double MAF,A,B;
	int dominance,MUTMOD;
	int pop1,pop2;
	
	pop1 = *p1;
	pop2 = *p2;
	D.n1 = *n1;
	D.n2 = *n2;
	GetRNGstate();
	if ((out = fopen(OUT,"r")) == NULL) {
		out = fopen(OUT,"w");
		fprintf(out,"%s\t %15s %15s %15s\n","Filename","F_1","F_2","F_ST(WC)");
	} else {
		fclose(out);
		out = fopen(OUT,"a");		
	}
	ReadParameterFile(filename,&dominance,&MAF,&A,&B,&totiter,&P,&MUTMOD,&nsets);
	infile = fopen(INFILE, "r");
	while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
	i = j = 0;
	while (fscanf(infile,"%s %d %d %lf %lf %lf",datafile_name,&i,&j,&D.F1,&D.F2,&dummy) != EOF) {
		if ((i == pop1) && (j == pop2)) break;
	}	
	if ((i == 0) && (j == 0)) {
		Rprintf("STOPPED: populations %d and %d do not exist in the dataset; cannot complete the simulations.\n",pop1,pop2);
		goto end;
	}
	strcat(datafile_name,"_");
	sprintf(y,"%d",D.n1);
	strcat(datafile_name,y);
	strcat(datafile_name,"_");
	sprintf(y,"%d",D.n2);
	strcat(datafile_name,y);
	strcat(datafile_name,SUFFIX);
	sum_num_1 = sum_num_2 = sum_den = sum_num_FST = sum_den_FST = 0;
	skip = (int *) malloc (nsets * sizeof (int));
	if (D.F1 > 0 && D.F2 > 0) {
		polymorphicsets = 0;
		for (set = 0; set < nsets; ++set) {                                           // Test for polymorphism
			P.N1 = 1 / (1 - pow((1 - D.F1),(1 / P.t[set])));
			P.N2 = 1 / (1 - pow((1 - D.F2),(1 / P.t[set])));
			polymorphic = 0;
			for (i = 0; i < TEST; ++i) {
				P.mu = rgamma(2,P.mubar / 2);                                          // Draw individual mutation rates per locus (see Excoffier et al. 2005)
				Make_Samples(P,D,set,dominance,MUTMOD);
				NewStatistics(&S,D,MAF,A,B,dominance);					
				for (j = 0; j < 4; ++j) {
					free(deme[j].tree);
					free(deme[j].list);
				}
				if (S.maf <= MAF) ++polymorphic;
			}
			if ((double) polymorphic / TEST < 0.10) {
				Rprintf("WARNING: cannot generate polymorphic data in file '%s' with parameter set %d; skipping that set.\n",datafile_name,(set + 1));
				skip[set] = 1;
			} else {
				skip[set] = 0;
				++polymorphicsets;
			}
		}
		if (!polymorphicsets) {
			Rprintf("STOPPED: cannot generate polymorphic data in file '%s' with any set of parameters.\n",datafile_name);
			goto end;
		}
		iter = (int) totiter / polymorphicsets;
		outfile = fopen(datafile_name, "w");
		for (set = 0; set < nsets; ++set) {
			if (!skip[set]) {
				P.N1 = 1 / (1 - pow((1 - D.F1),(1 / P.t[set])));
				P.N2 = 1 / (1 - pow((1 - D.F2),(1 / P.t[set])));
				if (set == (nsets - 1)) iter += totiter - polymorphicsets * iter;
				for (i = 0; i < iter; ++i) {
					do {                                                                       // the gamma distribution is parametrized by scale rather than rate.
						P.mu = rgamma(2,P.mubar / 2);                                          // Draw individual mutation rates per locus (see Excoffier et al. 2005)
						Make_Samples(P,D,set,dominance,MUTMOD);
						NewStatistics(&S,D,MAF,A,B,dominance);					
						for (j = 0; j < 4; ++j) {
							free(deme[j].tree);
							free(deme[j].list);
						}
					}
					while (S.maf > MAF);
					sum_num_1 += S.num_1;
					sum_num_2 += S.num_2;
					sum_den += S.den;
					sum_num_FST += S.num_FST;
					sum_den_FST += S.den_FST;
					fprintf(outfile,"%15.9f %15.9f %15.9f %15.9f %15d\n",S.F[1],S.F[2],S.F[ST],S.H,S.K[2]);
				}
			}
		}
		fclose(outfile);
		fprintf(out,datafile_name,"\t");
		fprintf(out,"%15.9f %15.9f %15.9f\n",(1 - sum_num_1 / sum_den),(1 - sum_num_2 / sum_den),(1 - sum_num_FST / sum_den_FST));			
	}
	else {
		Rprintf("STOPPED: multilocus estimates of differentiation in populations %d and %d are negative; cannot complete the simulations.\n",pop1,pop2);
	}
	end :
	fclose(infile);
	fclose(out);
	free(skip);
	free(P.t);
	free(P.N0);
	free(P.t0);
	free(P.Ne);
	PutRNGstate();
}


void ReadParameterFile(char filename[32],
                       int *dominance,
					   double *MAF,
					   double *A,
					   double *B,
					   int *totiter,
					   Parameters *P,
					   int *MUTMOD,
					   int *nsets)

{
	FILE *parameters;	
	int X;
	int n;
	
	parameters = fopen(PARAMETERFILE,"r");
	fscanf(parameters,"%s",filename);
	while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	fscanf(parameters,"%d",dominance);
	while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	fscanf(parameters,"%lf",MAF);
	while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	if (*dominance) {
		fscanf(parameters,"%lf",A);
		while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
		fscanf(parameters,"%lf",B);
		while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	}
	fscanf(parameters,"%d",totiter);
	while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	fscanf(parameters,"%lf",&(P -> mubar));
	while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	if (!(*dominance)) {
		fscanf(parameters,"%d",MUTMOD);
		while(!((X = getc(parameters)) == '\n' || X == '\f' || X == '\r'));
	} else {
		*MUTMOD = 2;
	}
	n = 0;
	P -> t = (double *) malloc (9999 * sizeof (double));
	P -> N0 = (double *) malloc (9999 * sizeof (double));
	P -> t0 = (double *) malloc (9999 * sizeof (double));
	P -> Ne = (double *) malloc (9999 * sizeof (double));
	while (fscanf(parameters,"%lf %lf %lf %lf",&(P -> t[n]),&(P -> N0[n]),&(P -> t0[n]),&(P -> Ne[n])) != EOF) {
		n++;
	}
	*nsets = n;
	P -> t = (double *) realloc (P -> t,n * sizeof (double));
	P -> N0 = (double *) realloc (P -> N0,n * sizeof (double));
	P -> t0 = (double *) realloc (P -> t0,n * sizeof (double));
	P -> Ne = (double *) realloc (P -> Ne,n * sizeof (double));
	fclose(parameters);
}

void Make_Samples(Parameters P,
				  Data D,
				  int set,
				  int dominance,
				  int MUTMOD)
{
	void Make_Tree(int t);
	void Add_Mutation(int MUTMOD,double mu, struct Node *node);
	
	int i,cpt,max_lineages,max_nodes[2];
	
	deme[0].begin = 0.0;
	deme[1].begin = 0.0;
	deme[0].end = deme[1].end = deme[2].begin = P.t[set];
	deme[2].end = deme[3].begin = P.t[set] + P.t0[set];
	deme[3].end = 1.7976931348623157e+308;                                         // Maximum number for a double (this is to ensure that all the lineages in the ancestral population coalesce)
	deme[2].N = P.N0[set];
	deme[0].N = P.N1;
	deme[1].N = P.N2;
	if (dominance) {
		deme[0].n = 2 * D.n1;                                                          // This is twice the number of individuals
		deme[1].n = 2 * D.n2;
		
	} else { 
		deme[0].n = D.n1;                                                          // This is twice the number of individuals
		deme[1].n = D.n2;
	}
	deme[3].N = P.Ne[set];
	max_lineages = max_nodes[0] = max_nodes[1] = 2 * (deme[0].n + deme[1].n) - 1;  // This is the maximum number of nodes in the complete genealogy
	deme[0].list = (struct Node **) malloc (max_lineages * sizeof (struct Node *));// These are the pointers to the lineages left in the deme
	deme[1].list = (struct Node **) malloc (max_lineages * sizeof (struct Node *));
	deme[0].tree =(struct Node *) malloc(max_nodes[0] * sizeof(struct Node));      // This is the genealogy
	deme[1].tree =(struct Node *) malloc(max_nodes[1] * sizeof(struct Node));
	for (i = 0; i < max_nodes[0]; ++i) {                                           // Initialization of terminal nodes in population 0
		(deme[0].tree[i]).ancestor = NULL;
		(deme[0].tree[i]).descendant_1 = NULL;
		(deme[0].tree[i]).descendant_2 = NULL;
	}
	for (i = 0; i < max_nodes[1]; ++i) {                                           // Initialization of terminal nodes in population 1
		(deme[1].tree[i]).ancestor = NULL;
		(deme[1].tree[i]).descendant_1 = NULL;
		(deme[1].tree[i]).descendant_2 = NULL;
	}
	for (i = 0; i < 2; ++i) {
		Make_Tree(i);
	}
	deme[2].n = deme[0].nbr_ancestors + deme[1].nbr_ancestors;                     // Initialization of tree[2] parameters (sample_size)
	deme[2].list = (struct Node **) malloc (deme[2].n * sizeof(struct Node *));
	deme[2].tree =(struct Node *) malloc((2 * deme[2].n - 1) * sizeof(struct Node));
	for (i = 0; i < 2 * deme[2].n - 1; ++i) {                                      // Initialization of all nodes in population 1
		(deme[2].tree[i]).ancestor = NULL;
		(deme[2].tree[i]).descendant_1 = NULL;
		(deme[2].tree[i]).descendant_2 = NULL;
	}
	for (i = 0; i < deme[0].nbr_ancestors; ++i) {                                  // Initialization of BOTH free and intermediate nodes in population 2
		deme[2].tree[i].descendant_1 = deme[0].list[i];
		deme[2].tree[i].descendant_2 = NULL;
		deme[0].list[i] -> ancestor = deme[2].tree + i;
	}
	for (i = 0; i < deme[1].nbr_ancestors; ++i) {                                  // Initialization of BOTH free and intermediate nodes in population 2
		deme[2].tree[i + deme[0].nbr_ancestors].descendant_1 = NULL;
		deme[2].tree[i + deme[0].nbr_ancestors].descendant_2 = deme[1].list[i];
		deme[1].list[i] -> ancestor = deme[2].tree + deme[0].nbr_ancestors + i;
	}
	Make_Tree(2);                                                                  // Simulation of coalescent process in population  2 (bottlenecked population)
	deme[3].n = deme[2].nbr_ancestors;                                             // Initialization of tree[0] parameters (sample_size value, and allocation of memory)
	deme[3].list = (struct Node **) malloc (deme[3].n * sizeof(struct Node *));
	deme[3].tree =(struct Node *) malloc((2 * deme[3].n - 1) * sizeof(struct Node)); // CHECK !!!!!!!!  WAS: deme[3].tree =(struct Node *) malloc(2 * deme[2].n * sizeof(struct Node));
	for (i = 0; i < 2 * deme[3].n - 1; ++i) {                                      // Initialization of all nodes in population 3
		(deme[3].tree[i]).ancestor = NULL;
		(deme[3].tree[i]).descendant_1 = NULL;
		(deme[3].tree[i]).descendant_2 = NULL;
	}
	cpt = (int) deme[3].n / 2;
	for (i = 0; i < cpt; ++i) {                                                    // Initialization of BOTH free and intermediate nodes in population 2
		deme[3].tree[i].descendant_1 = deme[2].list[i];
		deme[3].tree[i].descendant_2 = NULL;
		deme[2].list[i] -> ancestor = deme[3].tree + i;
	}
	for (i = cpt; i < deme[3].n; ++i) {                                            // Initialization of BOTH free and intermediate nodes in population 2
		deme[3].tree[i].descendant_1 = NULL;
		deme[3].tree[i].descendant_2 = deme[2].list[i];
		deme[2].list[i] -> ancestor = deme[3].tree + i;
	}
	Make_Tree(3);                                                                  // Simulation of coalescent process in ancestral population at mutation-drift equilibrium
	if (MUTMOD > 1) {
		deme[3].tree[2 * deme[3].n - 2].allele = (int) (MUTMOD * unif_rand());
	}
	else {
		deme[3].tree[2 * deme[3].n - 2].allele = k_iam = 0;
	}
	Add_Mutation(MUTMOD,P.mu,&deme[3].tree[2 * deme[3].n - 2]);	
}

void Make_Tree(int t)                                                            // Realization of the coalescent process for a sample of N alleles.
{                                                                                // The process lasts until time = final_time.
	int i,number1,number2;
	double time = deme[t].begin;
	
	for (i = 0; i < deme[t].n; ++i) {
		deme[t].tree[i].time = time;
		deme[t].list[i] = deme[t].tree + i;
	}
	while (i > 1 && time < deme[t].end) {
		time += deme[t].N * (-2.0 * log(1 - unif_rand()) / (((double)i) * (i - 1)));
		if (time < deme[t].end) {
			number1 = (int) (unif_rand() * i);
			do {
				number2 = (int) (unif_rand() * i);
			}
			while (number2 == number1);
			deme[t].tree[2 * deme[t].n - i].time = time;
			deme[t].list[number1] -> ancestor = deme[t].tree + 2 * deme[t].n - i;
			deme[t].list[number2] -> ancestor = deme[t].tree + 2 * deme[t].n - i;
			deme[t].tree[2 * deme[t].n - i].descendant_1 = deme[t].list[number1];
			deme[t].tree[2 * deme[t].n - i].descendant_2 = deme[t].list[number2];
			deme[t].list[number1] = deme[t].tree + 2 * deme[t].n - i;
			deme[t].list[number2] = deme[t].list[--i];			
		}
	}
	deme[t].nbr_ancestors = i;
}

void Add_Mutation(int MUTMOD,double mu, struct Node *node)                                  // Add mutation, reccursively, troughout the total tree
{
	int Mutation(int MUTMOD,int allele);
	
	int i;
	int nbr_mut,tmp;
	double time;
	
	if (node -> descendant_1 != NULL) {
		(node -> descendant_1) -> allele = node -> allele;
		time = node -> time - ((node -> descendant_1) -> time);
		nbr_mut = (int) rpois(time * mu);                                            // Calculate the number of mutations along the branch
		for (i = 0; i < nbr_mut; ++i) {
			do {                                                                       // Give a new allelic state, [...]
				tmp = Mutation(MUTMOD,(node -> descendant_1) -> allele);                        // bug fixed 17-01-2007 : was 'tmp = Mutation(node -> allele)'
			}
			while (tmp == (node -> descendant_1) -> allele);                           // [...] different from the current state
			(node -> descendant_1) -> allele = tmp;                                     // Copy the mutation into the current node
		}
		Add_Mutation(MUTMOD,mu, node -> descendant_1);                                      // Play it again, starting from the current node
	}
	if (node -> descendant_2 != NULL) {
		(node -> descendant_2) -> allele = node -> allele;
		time = node -> time - ((node -> descendant_2) -> time);
		nbr_mut = (int) rpois(time * mu);                                            // Calculate the number of mutations along the branch
		for (i = 0; i < nbr_mut; ++i) {
			do {                                                                       // Give a new allelic state, [...]
				tmp = Mutation(MUTMOD,(node -> descendant_2) -> allele);                        // bug fixed 17-01-2007 : was 'tmp = Mutation(node -> allele)'
			}
			while (tmp == (node -> descendant_2) -> allele);                           // [...] different from the current state
			(node -> descendant_2) -> allele = tmp;                                     // Copy the mutation into the current node
		}
		Add_Mutation(MUTMOD,mu, node -> descendant_2);                                      // Play it again, starting from the current node
	}
}

int Mutation(int MUTMOD,
			 int allele)                                                         // Generation of a new allele, through mutation
{
	int new_allele;  
	
	if (!MUTMOD) {
		++k_iam;
		return(k_iam);
	}
	if (MUTMOD == 1) {
		new_allele = allele + 2 * (int) (2 * unif_rand()) - 1;
		return(new_allele);
	}
	if (MUTMOD == 2) {
		new_allele = 1 - allele;
		return(new_allele);
	}
	else {
		new_allele = (int) (MUTMOD * unif_rand());
		return(new_allele);
	}
}

void NewStatistics(Stats *S,
                   Data D,
				   double MAF,
				   double A,
				   double B,
				   int dominance)
{ 
	int **tmp,**n;                                                                 // tmp is a dummy array to list all possible allelic states
	double **p;                                                                     // n is the [2 x cpt] array that gives allele counts in both pops
	int i,j,k,u;
	long min,max;
	double n_c,nsum,nsum2,nbar,xx,yy,q2,q3;
	double var = 0.0;
	int recessive,all1,all2,x;
	int **genotype = NULL;
	
	min = 2147483647L;
	max = -2147483647L;
	for (i = 0; i < 2; ++i) {
		for (j = 0; j < deme[i].n; ++j) {
			if (deme[i].tree[j].allele < min) {
				min = deme[i].tree[j].allele;
			}
			if (deme[i].tree[j].allele > max) {
				max = deme[i].tree[j].allele;
			}
		}
	}
	if (min < 0) {
		for (i = 0; i < deme[0].n; ++i) {
			deme[0].tree[i].allele -= min;
		}
		for (i = 0; i < deme[1].n; ++i) {
			deme[1].tree[i].allele -= min;
		}
		max += 1 - min;                                                              // new 14-12-2006
	}
	else {
		max += 1;                                                                    // new 14-12-2006
	}
	tmp = (int **) malloc(2 * sizeof (int *));
	for (i = 0; i < 2; ++i) {
		tmp[i] = (int *) malloc(max * sizeof (int));
	}
	for(i = 0; i < 2; ++i) {
		for(j = 0; j < max; ++j) {
			tmp[i][j] = 0;
		}
	}
	for (i = 0; i < deme[0].n; ++i) {
		++tmp[0][deme[0].tree[i].allele];
	}
	for (i = 0; i < deme[1].n; ++i) {
		++tmp[1][deme[1].tree[i].allele];
	}
	S -> K[0] = 0;                                                                 // K[2] gives the total number of alleles in the pooled sample
	S -> K[1] = 0;
	S -> K[2] = 0;
	for (i = 0; i < max; ++i) {
		if (tmp[0][i] + tmp[1][i] != 0) (S -> K[2]) += 1;
		if (tmp[0][i] != 0) (S -> K[0]) += 1;
		if (tmp[1][i] != 0) (S -> K[1]) += 1;
	}
	n = (int **) malloc(2 * sizeof (int *));                                     // n is the [2 x cpt] array that gives allele counts in both pops
	for (i = 0; i < 2; ++i) {
		n[i] = (int *) malloc((S -> K[2]) * sizeof (int));
	}
	k = 0;
	for (i = 0; i < max; ++i) {
		if (tmp[0][i] + tmp[1][i] != 0) {
			n[0][k] = tmp[0][i];
			n[1][k] = tmp[1][i];
			++k;
		}
	}
	S -> maf = -2147483647L;
	p = (double **) malloc(3 * sizeof(double *));
	for (i = 0; i < 3; ++i) {
		p[i] = (double *) malloc((S -> K[2]) * sizeof(double));
	}
	if (dominance && ((S -> K[2]) > 1)) {							// 08-03-2011 : problem in the following when 1 allele fixed
		genotype = (int **) malloc(2 * sizeof (int *));
		for (i = 0; i < 2; ++i) {
			genotype[i] = (int *) malloc(2 * sizeof (int));
		}
		recessive = (int) (2 * unif_rand());                                     // Draw which allele is recessive
		for (i = 0; i < 2; ++i) {
			genotype[i][0] = genotype[i][1] = 0;
			while (n[i][0] > 0 || n[i][1] > 0) {
				x = (int) ((n[i][0] + n[i][1]) * unif_rand()) + 1;
				if(x <= n[i][0]) {
					--n[i][0];
					all1 = 0;
				}
				else {
					--n[i][1];
					all1 = 1;
				}
				x = (int) ((n[i][0] + n[i][1]) * unif_rand()) + 1;
				if(x <= n[i][0]) {
					--n[i][0];
					all2 = 0;
				}
				else {
					--n[i][1];
					all2 = 1;
				}
				if (all1 == recessive && all2 == recessive) {
					++genotype[i][0];
				}
				else {
					++genotype[i][1];
				}
			}
		}
		for (u = 0; u < S -> K[2]; ++u) {
			p[2][u] = (double) (genotype[0][u] + genotype[1][u]) / (D.n1 + D.n2);
			if (p[2][u] > S -> maf) {
				S -> maf = p[2][u];
			}
/*			if ((1 - p[2][u]) > S -> maf) {			BUG 17-09-2012
				S -> maf = (1 - p[2][u]);
			}*/			
		}
	} else {
		for (u = 0; u < S -> K[2]; ++u) {
			p[2][u] = (double) (n[0][u] + n[1][u]) / (D.n1 + D.n2);
			if (p[2][u] > S -> maf) {
				S -> maf = p[2][u];
			}
/*			if ((1 - p[2][u]) > S -> maf) {			BUG 17-09-2012
				S -> maf = (1 - p[2][u]);
			}*/			
		}
	}
	if (S -> maf <= MAF) {
		if (dominance) {
			p[0][0] = exp(lgammafn(genotype[0][0] + A + 0.5) + lgammafn(D.n1 + A + B) - lgammafn(genotype[0][0] + A) - lgammafn(D.n1 + A + B + 0.5)); // Zhivotovsky 1999 eq. (10)
			p[1][0] = exp(lgammafn(genotype[1][0] + A + 0.5) + lgammafn(D.n2 + A + B) - lgammafn(genotype[1][0] + A) - lgammafn(D.n2 + A + B + 0.5));
			p[2][0] = exp(lgammafn(genotype[0][0] + genotype[1][0] + A + 0.5) + lgammafn(D.n1 + D.n2 + A + B) - lgammafn(genotype[0][0] + genotype[1][0] + A) - lgammafn(D.n1 + D.n2 + A + B + 0.5));
			for (k = 0; k < 3; ++k) {
				p[k][1] = 1.0 - p[k][0];
			}
			var = exp(lgammafn(genotype[0][0] + genotype[1][0] + A + 1.0) + lgammafn(D.n1 + D.n2 + A + B) - lgammafn(genotype[0][0] + genotype[1][0] + A) - lgammafn(D.n1 + D.n2 + A + B + 1.0)) - p[2][0] * p[2][0];			
		} else {
			for (u = 0; u < S -> K[2]; ++u) {
				p[0][u] = (double) n[0][u] / D.n1;
				p[1][u] = (double) n[1][u] / D.n2;
			}
		}
		S -> num_1 = S -> num_2 = S -> den = 0;
		nsum = (double) D.n1 + D.n2;
		nsum2 = pow(D.n1,2) + pow(D.n2,2);
		n_c = nsum - nsum2 / nsum;
		for (u = 0; u < S -> K[2]; ++u) {
		    S -> num_1 += p[0][u] * (1 - p[0][u]);
			S -> num_2 += p[1][u] * (1 - p[1][u]);
			S -> den += (double) D.n1 * pow((p[0][u] - p[2][u]),2) + (D.n1 - pow(D.n1,2) / (D.n1 + D.n2)) * p[0][u] * (1 - p[0][u]);
			S -> den += (double) D.n2 * pow((p[1][u] - p[2][u]),2) + (D.n2 - pow(D.n2,2) / (D.n1 + D.n2)) * p[1][u] * (1 - p[1][u]);
		}
		S -> num_1 *= n_c * D.n1 / (D.n1 - 1);
		S -> num_2 *= n_c * D.n2 / (D.n2 - 1);
		if (S -> den > 0) S -> F[1] = 1 - (S -> num_1 / S -> den); else S -> F[1] = 9.0;
		if (S -> den > 0) S -> F[2] = 1 - (S -> num_2 / S -> den); else S -> F[2] = 9.0;
		if (dominance) {
			S -> H = 2 * p[2][0] * (1 - p[2][0]) + 2 * var;
		} else {
			S -> H = 0;
			for (u = 0; u < S -> K[2]; ++u) {
				S -> H += p[2][u] * (1 - p[2][u]);
			}
			S -> H *= (double) (D.n1 + D.n2) / (D.n1 + D.n2 - 1);
		}
		nbar = nsum / 2;
		xx = yy = 0.0;
		for (u = 0; u < S -> K[2]; ++u) {
		    xx += D.n1 * pow(p[0][u],2);
			xx += D.n2 * pow(p[1][u],2);
			yy += pow((p[0][u] * D.n1 + p[1][u] * D.n2),2);
		}
		q2 = (xx - 2) / (2 * (nbar - 1.0));
		q3 = 1.0 / (2 * nbar * n_c) * (yy - nbar * (n_c - 1.0) / (nbar - 1.0) * xx) + (nbar - n_c) / (n_c * (nbar - 1.0)) * (1.0 - xx);
		S -> num_FST = (1 - q2);
		S -> den_FST = (1 - q3);
		if ((1 - q3) > 0) S -> F[ST] = 1.0 - (1 - q2) / (1 - q3); else S -> F[ST] = 9.0;
	}	
	if (dominance && ((S -> K[2]) > 1)) {
		for (i = 0; i < 2; ++i) {
			free(genotype[i]);
		}
		free(genotype);
	}	
	for (k = 0; k < 3; ++k) {
		free(p[k]);
	}
	free(p);		
	for (i = 0; i < 2; ++i) {
		free(tmp[i]);
		free(n[i]);
	}
	free(n);	
	free(tmp);
}
