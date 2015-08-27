#include "STeM.h"

#define MAX_NUM_GENES 100
#define MAX_NUM_CHROM 5000
#define MAX_GEN_LENGTH 32
#define MAX_NUM_FLEX_RES 100
#define IM1 2147483563 
#define IM2 2147483399 
#define AM (1.0/IM1) 
#define IMM1 (IM1-1) 
#define IA1 40014 
#define IA2 40692 
#define IQ1 53668 
#define IQ2 52774 
#define IR1 12211 
#define IR2 3791 
#define NTAB 32 
#define NDIV (1+IMM1/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS) 

void load_eigen(gsl_vector *v,gsl_matrix *m,char filename[100],int atom)
{
	FILE *file;
 	file = fopen(filename,"r");
 	char line[10000];
 	float temp;
 	int i,j;
 	while(fgets(line,10000,file))
	{
 		if (3 == sscanf(line,"%d\t%d\t%f",&i,&j,&temp))
		{
 			if (i>atom || j>atom)
			{
 				//printf("Index to big: (%d,%d) %d\n",i,j,atom); 
 			
 			} 
 			else
			{
 				gsl_matrix_set(m,i-1,j-1,temp);
			}
		}
		
 		if (2 == sscanf(line,"%dth Eigenvalue:%f",&i,&temp))
		{
 			if (i>atom) 
			{
 				printf("Index to big: (%d) %d\n",i,atom); 
 			}
 			else
			{
 				gsl_vector_set(v,i-1,temp);
 			}
 		}
 	}
 	
 	fclose(file);
}
  void write_eigen(char filename[100], gsl_matrix *m,gsl_vector *v,int nb_atom) {
 	FILE *out_file;
 	out_file = fopen(filename,"w");
	 int i,j;
	 for(i=0; i < nb_atom; ++i) {
	 	fprintf (out_file,"%dth Eigenvalue:%14.10f\n",i+1,gsl_vector_get (v, i));
	 	for(j=0; j < nb_atom; ++j) {
	 		fprintf (out_file,"%d\t%d\t%10.6f\n",j+1,i+1,gsl_matrix_get (m, j, i));
	 	}
	 	fprintf (out_file,"\n\n");
	 }
	 fclose(out_file);
 }
 void write_matrix(char filename[100], gsl_matrix *m,int nb_atom,int nb_atom_2) {
 	FILE *out_file;
 	out_file = fopen(filename,"w");
	int i,j;
	 for(i=0; i < nb_atom; ++i) {
	 	for(j=0; j < nb_atom_2; ++j) {
	 		fprintf (out_file,"%12.6g ",gsl_matrix_get (m, i, j));
	 	}
	 fprintf (out_file,"\n");
	 }
	 fclose(out_file);
 }
 
 
 float ran2(long *idum){
  int j; 
  long k; 
  static long idum2=123456789; 
  static long iy=0; 
  static long iv[NTAB]; 
  float temp; 

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum); 
    idum2=(*idum); 
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1; 
      *idum=IA1*(*idum-k*IQ1)-k*IR1; 
      if (*idum < 0) *idum += IM1; 
      if (j < NTAB) iv[j] = *idum; 
    } 
    iy=iv[0]; 
  } 
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1; 
  k=idum2/IQ2; 
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2; 
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum; 
  if (iy < 1) iy += IMM1; 
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp; 
}

int load_matrix(gsl_matrix *newstrc,char filename[100])
{
	FILE *file;
 	file = fopen(filename,"r");
 	char line[40000];
 	float temp;
 	char *ptr;
 	int count;
 	int i = 0,j;
	int flag = 0;
 	while(fgets(line,40000,file))
 	{
		ptr = line;
		j = 0;
		sscanf(line,"%g",&temp);
		//printf("I STORE:%d %d %f\n",i,j,temp);
		gsl_matrix_set(newstrc,i,j,temp);
		count = 0;
		flag =0;
		while(1)
		{
			if (flag == 0 && count == 0 && (line[count] == ' ' || line[count] == '\t')) {flag = 2;continue;} // Espace en commenceant
			if (flag == 2 && (line[count] != ' ' || line[count] != '\t')) {flag = 1;} 
			++count;
			//printf("Count:%d Flag:%d\n",count,flag);
			if (line[count] == '\0' || line[count] == '\n') {break;}
			//printf("%s\n",ptr+count);
			
			if (flag == 1) {if (line[count] == ' ' || line[count] == '\t') {continue;} else {flag = 0;}}
			
			if ((line[count] == ' ' || line[count] == '\t' ) && flag == 0)
			{
				++j;
				
				if(line[count+1] == '\0' || line[count+1] == '\n') {break;}
					
				if (1 == sscanf(ptr+count+1,"%g",&temp))
				{
					//printf("I STORE:%d %d %f\n",i,j,temp);
					gsl_matrix_set(newstrc,i,j,temp);
					flag = 1;
				}
			}
		}
		
		++i;
	}
	return(i);
}
void load_eigen_grid(gsl_vector *v,gsl_matrix *m,char filename[100],int atom, int mode) {
	FILE *file;
 	file = fopen(filename,"r");
 	char line[10000];
 	float temp;
 	int i,j;
 	while(fgets(line,10000,file)) {
 		if (3 == sscanf(line,"%d\t%d\t%f",&i,&j,&temp)) {gsl_matrix_set(m,i-1,j-1,temp);}
 		if (2 == sscanf(line,"%dth Eigenvalue:%f",&i,&temp)) {gsl_vector_set(v,i-1,temp);if (i > mode+1) {fclose(file);return;}}
 	}
 	fclose(file);
}

void write_grid_mat(char filename[500], gsl_matrix *m,int nb_atom,int nb_atom_1) {
 	FILE *out_file;
 	out_file = fopen(filename,"w"); 	
	int i,j;
	for(i=0; i < nb_atom; ++i) {
	 for(j=0; j < nb_atom_1; ++j) {

	 		fprintf (out_file,"%8.3f ",gsl_matrix_get (m, i, j));
	 	}
	 fprintf (out_file,"\n");	
	 }
	 fclose(out_file);
 
 }
 
 void write_eigen_mat(char filename[500], gsl_matrix *m,int nb_atom,int nb_atom_1,int mode, int nbr_mode) {
 	FILE *out_file;
 	out_file = fopen(filename,"w"); 	
	 int i,j;
for(i=0; i < nb_atom_1; ++i) {
	 for(j=mode; j < nbr_mode+mode; ++j) {

	 		fprintf (out_file,"%8.3f ",gsl_matrix_get (m, i, j));
	 	}
	 fprintf (out_file,"\n");	
	 }
	 fclose(out_file);
 
 }
  double dot_p(double a[3], double b[3]) {
 	double dot;
 	dot = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
 	return(dot);
 }
 

 float write_strc(char filename[100], struct pdb_atom *newstrc,int nb_atom,float factor)
 {
 	FILE *out_file;
	
	int k;
	
	out_file = fopen(filename,"w");
	
	float avg = 0;
	
	float tot = 1;
	
	for (k = 0; k<nb_atom;k++)
	{
		if (newstrc[k].atom_type == 1) {fprintf(out_file,"ATOM  ");}
		if (newstrc[k].atom_type == 4) {fprintf(out_file,"ATOM  ");}
		if (newstrc[k].atom_type == 5) {fprintf(out_file,"HETATM");}
	 	if (newstrc[k].atom_type == 2) {fprintf(out_file,"HETATM");}
	 	if (newstrc[k].atom_type == 3) {fprintf(out_file,"HETATM");}
	 	if (newstrc[k].atom_type == 4) {fprintf(out_file,"ATOM  ");}
	 	if (newstrc[k].atom_type == 5) {fprintf(out_file,"HETATM");}
		fprintf(out_file,"%5.d %s%s %s%4d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
			newstrc[k].atom_number,
 			newstrc[k].atom_prot_type,
 			newstrc[k].res_type,
 			newstrc[k].chain,
 			newstrc[k].res_number,
 			newstrc[k].x_cord,
 			newstrc[k].y_cord,
 			newstrc[k].z_cord,
 			newstrc[k].b_factor
 			);
		
		if (newstrc[k].main_vars[0] > 0)
		{
			// Rebuild matrix
			
			gsl_matrix *temp = gsl_matrix_alloc(3,3);
			gsl_matrix *temp2 = gsl_matrix_alloc(3,3);
			gsl_matrix *evec = gsl_matrix_alloc(3,3);
			gsl_matrix *tevec = gsl_matrix_alloc(3,3);
			gsl_matrix *eval = gsl_matrix_alloc(3,3);
			gsl_matrix_set_all(temp,0);
			gsl_matrix_set_all(temp2,0);
			gsl_matrix_set_all(eval,0);
			
			
			int i,j;
			
			for(i=0;i<3;++i)
			{
				for(j=0;j<3;++j)
				{
					gsl_matrix_set(evec,i,j,newstrc[k].global_evecs[i][j]);
					gsl_matrix_set(tevec,i,j,newstrc[k].global_evecs[j][i]);	
				}
				
				gsl_matrix_set(eval,i,i,newstrc[k].main_vars[i]);
			}
			
			multiplie_matrix(evec,3,3,eval,3,3,temp2);
			multiplie_matrix(temp2,3,3,tevec,3,3,temp);
			
			fprintf(out_file,"ANISOU");
			fprintf(out_file,"%5.d %s%s %s%4d  %7d%7d%7d%7d%7d%7d\n",
				newstrc[k].atom_number,
				newstrc[k].atom_prot_type,
				newstrc[k].res_type,
				newstrc[k].chain,
				newstrc[k].res_number,
				(int) (gsl_matrix_get(temp,0,0)*factor),
				(int) (gsl_matrix_get(temp,1,1)*factor),
				(int) (gsl_matrix_get(temp,2,2)*factor),
				(int) (gsl_matrix_get(temp,0,1)*factor),
				(int) (gsl_matrix_get(temp,0,2)*factor),
				(int) (gsl_matrix_get(temp,1,2)*factor)
				);
			
			tot += 6.0;
			avg += gsl_matrix_get(temp,0,0)*factor;
			avg += gsl_matrix_get(temp,1,1)*factor;
			avg += gsl_matrix_get(temp,2,2)*factor;
			avg += gsl_matrix_get(temp,0,1)*factor;
			avg += gsl_matrix_get(temp,0,2)*factor;
			avg += gsl_matrix_get(temp,1,2)*factor;
			
			gsl_matrix_free(temp);
			gsl_matrix_free(temp2);
			gsl_matrix_free(eval);
			gsl_matrix_free(evec);
			gsl_matrix_free(tevec);
		}
 	}
 	 fclose(out_file);
 	 return(avg/tot);
 }
 
 void write_anisou_file(char filename[100], struct pdb_atom *newstrc,int nb_atom, int ihess)
 {
	 FILE *out_file;
	 
	 int k;
	 
	 out_file = fopen(filename,"w");
	 
	 for (k = 0; k<nb_atom;k++)
	 {
		 if (newstrc[k].atom_type == 1) {fprintf(out_file,"ATOM  ");}
		 if (newstrc[k].atom_type == 2) {fprintf(out_file,"HETATM");}
		 if (newstrc[k].atom_type == 3) {fprintf(out_file,"HETATM");}
		 fprintf(out_file,"%5.d %s%s %s%4d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
			 newstrc[k].atom_number,
			 newstrc[k].atom_prot_type,
			 newstrc[k].res_type,
			 newstrc[k].chain,
			 newstrc[k].res_number,
			 newstrc[k].x_cord,
			 newstrc[k].y_cord,
			 newstrc[k].z_cord,
			 newstrc[k].b_factor
			 );
			 
		fprintf(out_file,"ANISOU");
		fprintf(out_file,"%5.d %s%s %s%4d  %7d%7d%7d%7d%7d%7d\n",
			newstrc[k].atom_number,
			newstrc[k].atom_prot_type,
			newstrc[k].res_type,
			newstrc[k].chain,
			newstrc[k].res_number,
			(int) (newstrc[k].covar[0][0]),
			(int) (newstrc[k].covar[1][1]),
			(int) (newstrc[k].covar[2][2]),
			(int) (newstrc[k].covar[0][1]),
			(int) (newstrc[k].covar[0][2]),
			(int) (newstrc[k].covar[1][2])
			);
	 }
	 fclose(out_file);
 }
 
 
 void write_strc_b(char filename[100], struct pdb_atom *old,int nb_atom,gsl_matrix *m, int node) {
 	int k, i = 0;
 	FILE *out_file;
 	out_file = fopen(filename,"w");
 	
 	
 	for(k=0;k<nb_atom;++k) {
 		
 		if ((old[k].node != old[k+1].node) && (k != nb_atom-1)) {
	 			++i;
		}
		if (old[k].atom_type == 1) {fprintf(out_file,"ATOM  ");}
		if (old[k].atom_type == 2) {fprintf(out_file,"HETATM");}
		if (old[k].atom_type == 3) {fprintf(out_file,"HETATM");}
			fprintf(out_file,"%5.d %s%s %s%4d%12.3f%8.3f%8.3f  1.00  %2.2f\n",
				old[k].atom_number,
				old[k].atom_prot_type,
				old[k].res_type,
				old[k].chain,
				old[k].res_number,
				old[k].x_cord,
				old[k].y_cord,
				old[k].z_cord,
				gsl_matrix_get(m,i,i)*2
		);
	}
 }
 void assignArray(gsl_matrix *m,double **a,int count,int count_1) {
	int i,j;
	for (i=0;i<count;++i) {
		for (j=0;j<count_1;++j) {
			gsl_matrix_set(m,i,j,a[i][j]);
		}
	}
}
 void correlate_sp(gsl_matrix *m,struct pdb_atom *strc, int atom) {
 	int i,k;
 	double moy_x = 0,moy_y = 0;
 	double std_x = 0,std_y=0;
 	double prod_tot=0;
 	for (k=0;k<3;++k) {
	 	for (i=0;i<atom/3;++i) {
	 	//printf("%f\t%f\n",gsl_matrix_get(m,i*3+k,i*3+k),strc[i*3+k].b_factor);
	 		moy_x += gsl_matrix_get(m,i*3+k,i*3+k);
	 		moy_y += strc[i*3+k].b_factor;
	 	}
	 	moy_x /= (atom/3);
	 	moy_y /= (atom/3);
	 	//printf("Moy X:%g\tMoy Y:%g\n",moy_x,moy_y);

	 	for (i=0;i<atom/3;++i) {
	 		std_x += (gsl_matrix_get(m,i*3+k,i*3+k)-moy_x)*(gsl_matrix_get(m,i*3+k,i*3+k)-moy_x);
	 		std_y += (strc[i*3+k].b_factor-moy_y)*(strc[i*3+k].b_factor-moy_y);
	 	}
	 	//printf("Std X:%g\tStd Y:%g\n",std_x,std_y);
	 	if (std_y == 0) {break;}


	 	for (i=0;i<atom/3;++i) {
	 		prod_tot += (gsl_matrix_get(m,i*3+k,i*3+k)-moy_x)*(strc[i*3+k].b_factor-moy_y);
	 	}
	 	//printf("Product Tot:%g\n",prod_tot);
	 	printf("%d:%g\n",k,prod_tot/sqrt(std_y*std_x));
 	}
 }
 
 int count_atom_CA_n(struct pdb_atom *all,int atom, int node, int ligand) {
 	int i,k=-1;
 	for (i=0;i<atom;++i) {
 		//printf("Type:-%s-\n",all[i].atom_prot_type);
 		if (strncmp(all[i].res_type,"HOH",3) == 0) {;continue;}
 		if (all[i].atom_type == 4 && (strncmp(all[i].atom_prot_type," P  ",4) == 0)) {++k;}
 		if ((strncmp(all[i].atom_prot_type," CA ",4) == 0) ||
 			((node == 3) && ((strncmp(all[i].atom_prot_type," N  ",4) == 0) || (strncmp(all[i].atom_prot_type," C  ",4) == 0))) ||
 			((ligand == 1) && (all[i].atom_type == 3))) {
 			++k;
	 	}
 	
 	}
 	return(k+1);
 }
 
 float correlate(gsl_matrix *m,struct pdb_atom *strc, int atom) {
 	int i;
 	double moy_x = 0,moy_y = 0;
 	double std_x = 0,std_y=0;
 	double prod_tot=0;
 	int c = 0;
 	float diff = 0;
 	for (i=0;i<atom;++i) {
 		if (strc[i].atom_type == 3){continue;}
 		if (strc[i].atom_type == 4 && strncmp(strc[i].atom_prot_type," P  ",4) != 0) {continue;} // added for bfact
 		if (strc[i].b_factor < 0.00001) {continue;} 																						 // added for bfact

 		++c;
 		moy_x += gsl_matrix_get(m,i,i);
 		moy_y += strc[i].b_factor;
 		diff += (gsl_matrix_get(m,i,i)-strc[i].b_factor)*(gsl_matrix_get(m,i,i)-strc[i].b_factor);
 	}
 	moy_x /= c;
 	moy_y /= c;
 	//printf("Moy X:%g\tMoy Y:%g\n",moy_x,moy_y);
 	
 	 for (i=0;i<atom;++i) {
 	 	if (strc[i].atom_type == 3){continue;}
 		if (strc[i].atom_type == 4 && strncmp(strc[i].atom_prot_type," P  ",4) != 0) {continue;} // added for bfact
 		if (strc[i].b_factor < 0.00001) {continue;} 																						 // added for bfact
 		std_x += (gsl_matrix_get(m,i,i)-moy_x)*(gsl_matrix_get(m,i,i)-moy_x);
 		std_y += (strc[i].b_factor-moy_y)*(strc[i].b_factor-moy_y);
 	}
 	//printf("Std X:%g\tStd Y:%g\n",std_x,std_y);
 	if (std_y == 0) {return(0);}

 	
 	 for (i=0;i<atom;++i) {
 	 	if (strc[i].atom_type == 3){continue;}
 		if (strc[i].atom_type == 4 && strncmp(strc[i].atom_prot_type," P  ",4) != 0) {continue;} // added for bfact
 		if (strc[i].b_factor < 0.00001) {continue;} 																						 // added for bfact
 		prod_tot += (gsl_matrix_get(m,i,i)-moy_x)*(strc[i].b_factor-moy_y);
 	}
 	//printf("Product Tot:%g\n",prod_tot);
 	//printf("Absolute Difference:%10.6f\n",diff);
 	return(prod_tot/sqrt(std_y*std_x));
 	
 }
 
 void k_inverse_matrix_stem(gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc,int mode,int nm) 
 {
 	//gsl_matrix *buffer = gsl_matrix_alloc(nb_atom, nb_atom); /*Matrix buffer a additionner*/
 	gsl_matrix_set_all (m, 0);
 	int i,j,k,l;
 	
 	for (k=mode;k<mode+nm;++k) 
	{
 		if (k > int (evl->size-1)) {break;}
 		
 		if  (gsl_vector_get (evl, k) < 0.000001) 
		{
 			printf("K = %d -> Eval to small I next:%f\n",k,gsl_vector_get (evl, k));
 			continue;
 		}
 		
 		for (i=0;i<nb_atom;++i) 
		{
 			for (j=0;j<nb_atom;++j) 
			{
				//if (i != j) {continue;}
		 		for (l=0;l<3;++l) 
				{	
			 			gsl_matrix_set(m,i,j,
			 				gsl_matrix_get(evc,i*3+l,k)*gsl_matrix_get(evc,j*3+l,k)/gsl_vector_get(evl, k) + gsl_matrix_get(m,i,j)
			 			);
				}
		 	}
	 	}
	 	//break;
 	}
 }
 
 void k_tot_inv_matrix_stem(gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc,int mode,int nm) 
 {
 	printf("In function\n");
	 //gsl_matrix *buffer = gsl_matrix_alloc(nb_atom, nb_atom); /*Matrix buffer a additionner*/
	 gsl_matrix_set_all (m, 0);
	 int i,j,k;
	 
	 for (k=mode;k<mode+nm;++k)
	 {
		 if (k > int (evl->size-1)) {break;}
		 if  (gsl_vector_get (evl, k) < 0.000001) 
		 {
			 printf("K = %d -> Eval to small I next:%f\n",k,gsl_vector_get (evl, k));
			 continue;
		 }
		 
		printf("%6.10f\n", gsl_vector_get (evl, k));
		 
		 for (i=0;i<3*nb_atom;++i)
		 {
			 for (j=0;j<3*nb_atom;++j)
			 {
					 gsl_matrix_set(m, i, j, gsl_matrix_get(evc, i, k)*gsl_matrix_get(evc, j, k)/gsl_vector_get(evl, k) + gsl_matrix_get(m, i, j));
			 }
		 }
		 //break;
	 }
 }
 
 void k_cov_inv_matrix_stem(gsl_matrix *m,int nb_atom, gsl_vector *evl,gsl_matrix *evc,int mode,int nm) 
 {
	 //gsl_matrix *buffer = gsl_matrix_alloc(nb_atom, nb_atom); /*Matrix buffer a additionner*/
	 gsl_matrix_set_all (m, 0);
	 int i,j,k,l;
	 
	 for (k=mode;k<mode+nm;++k)
	 {
		 if (k > int (evl->size-1)) {break;}
		 if  (gsl_vector_get (evl, k) < 0.000001) 
		 {
			 //printf("K = %d -> Eval to small I next:%g\n",k,gsl_vector_get (evl, k));
			 continue;
		 }
		 
 		 //printf("%6.10f\n", gsl_vector_get (evl, k));
		 
		 for (i=0;i<nb_atom*3;++i)
		 {
			for (j=0;j<nb_atom*3;++j) 
			{
				//	printf("I:%d et J:%d -> %f * %f -> %f\n",i,j,gsl_matrix_get(evc,i,k),gsl_matrix_get(evc, j, k),gsl_matrix_get(evc,i,k)*gsl_matrix_get(evc, j, k));

				 	
					gsl_matrix_set(m, i,j,
						gsl_matrix_get(evc,i,k)*gsl_matrix_get(evc, j, k)/gsl_vector_get(evl, k) 
						+ gsl_matrix_get(m, i,j)
					);
				
			 }
		 }
	//	break;
	 }
 }
 
 void buid_pre(gsl_matrix *m,gsl_matrix *evc,int nb_atom) {
 	gsl_matrix_set_all (m, 0);
 	int i,k,l;
 	double temp;
	for (k=0;k < nb_atom*3; ++k) {
		for (i = 0;i<nb_atom;++i) {
	 		for (l=0;l<3;++l) {
	 			temp = gsl_matrix_get (evc, (i+1)*3-l, k)*gsl_matrix_get (evc, (i+1)*3-l, k);
	 			//printf("I:%d K:%d Temp:%f\n",i,k,temp);
	 			gsl_matrix_set(m, i, k, temp+gsl_matrix_get(m,i,k));
			}
	 	}
	}
 
 }

 
long time_seed(){
  long s;
  //long s_max;
  time_t timer;
  struct tm *t;
  int val[8],i;


  //s_max=2147483647; // largest integer


   /*
     struct tm is a structure used to hold the time and date. 
     Its members are as follows:

     int tm_sec;      seconds after the minute (0 to 61)  
     int tm_min;      minutes after the hour (0 to 59)  
     int tm_hour;     hours since midnight (0 to 23)  
     int tm_mday;     day of the month (1 to 31)  
     int tm_mon;      months since January (0 to 11)  
     int tm_year;     years since 1900  
     int tm_wday;     days since Sunday (0 to 6 Sunday=0)  
     int tm_yday;     days since January 1 (0 to 365)  
     int tm_isdst;    Daylight Savings Time  
   */


  timer=time(NULL);
  t=localtime(&timer);
  
  
  //printf("The current time is %s.\n",asctime(localtime(&timer)));


  val[0]=t->tm_sec+1;
  val[1]=t->tm_min+1;
  val[2]=t->tm_hour+1;
  val[3]=t->tm_mday+1;
  val[4]=t->tm_mon+1; 
  val[5]=t->tm_year+1;
  val[6]=t->tm_wday+1;
  val[7]=t->tm_yday+1;

  s=1;
  for(i=0;i<=7;i++){
    s *= val[i];
    //printf("%d ",val[i]);
  }
  //printf("\n");
  //printf("s=%lu\n",s);

  //s /= 2;
  //if (s < 0){s *= -1;}

  return s;
}

float calc_energy(int atom,gsl_vector *eval,float t,int m) {
	int i;
	double ental=0.00000;
	double entro=0.00000;
	double Na = 6.022 * pow(10,23);
	double h = 6.626068 * pow(10,-34);
	double v;
	double e =2.7182;
	double power;

	double k = 1.3806503 *pow(10,-23);
	double R= 8.31447;
	double sum = 0.000000000;
	for (i = m;i<atom*3;++i) {
		v = gsl_vector_get(eval,i);
		//printf("V:%f\n",v);
		power = pow(e,-h*v/(k*t));
		ental += 0.5*Na*h*v+(Na*h*v*power)/(1-power);
		float rond = h*v/k;
		//entro += -R * log(1-power)+(Na*h*v*power)/(t*(1-power));
		entro += (rond/t)/(pow(e,rond/t)-1)-log(1-pow(e,-rond/t));
		sum += log(1.0/v);
	}
	
	//printf("Ental contribution:%f\n",ental/1000);
	//printf("Temperature:%f\n",t);
	//printf("Entro contribution:%f\n",t*entro/1000);
	//printf("Ln Summation:%f\n",sum);		
	return(sum);
	//return((ental-t*entro)/1000);
}

void assign_anisou_all(struct pdb_atom *init,int atom,struct pdb_atom *strc,int all, int ihess)
{
	int i,j,k,l;
	
	for(i = 0;i<atom;++i) 
	{
		for(j=0;j<all;++j) 
		{
			if (init[i].node == strc[j].node) 
			{
				for(k=0;k<3;++k) 
				{
					for(l=0;l<3;++l)
					{
						strc[j].covar[k][l] = init[i].covar[k][l];
					}
				}
			}
		}
	}
}
