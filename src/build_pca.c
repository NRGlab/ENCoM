#include "STeM.h"

float covar_vince(float *a,float *b,int atom);

double vector_product(float *a,gsl_matrix *evec,int m,int atom) {
	int i;
	float sum = 0;
	for (i = 0;i<atom;++i) {
//		printf("I:%d -> %f += %f * %f\n",i,sum,a[i],gsl_matrix_get(evec,i,m));
	if (a[i] < -9999.0) {continue;}
		sum += a[i] * gsl_matrix_get(evec,i,m);
	}
	return(sum);
}



int main(int argc, char *argv[]) {
	int i,j,k;
	//int prot = 0;
	printf("In function\n");

	int nbrAT[5000];
	
	int nbrCA[5000];
	//int* align[argc-1];
	
	struct pdb_atom* allpdb[5000];
	struct pdb_atom* allall[5000];
	int toread = 0;
	for (i = 0;i < argc;i++) {
		if (strcmp("-il",argv[i]) == 0) {++toread;continue;}
		if (toread == 0) {continue;}
    printf("I:%d/%d %s\n",i,argc,argv[i]);
    
    // Build Strc principale, only CA
		
		nbrAT[toread-1] = count_atom(argv[i]);
		int nconn = count_connect(argv[i]);
		
		int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(k=0;k<nconn;k++) { connect_h[k]=(int *)malloc(6*sizeof(int));}
   
    assign_connect(argv[i],connect_h);
    
    struct pdb_atom* strc_all;
		
		strc_all = (pdb_atom*)malloc(nbrAT[toread-1]*sizeof(pdb_atom));
		if (strc_all == NULL) {printf("Cannot assign strc\n");return(1);}
    
    nbrCA[toread-1] = build_all_strc(argv[i],strc_all);
		check_lig(strc_all,connect_h,nconn,nbrAT[toread-1]);
		
		struct pdb_atom* strc_node;
		strc_node = (pdb_atom*)malloc(nbrCA[toread-1]*sizeof(pdb_atom));
		if (strc_node == NULL) {return(1);}
		//printf("In function:%d \n",nconn[toread-1]);
											
		nbrCA[toread-1] = build_cord_CA(strc_all, strc_node,nbrAT[toread-1],0,connect_h,nconn);
    
    allpdb[toread-1] = strc_node;
		
		allall[toread-1] = strc_all;
		
    
    for(k=0;k<nconn;k++) { free(connect_h[k]);}
    free(connect_h);
    
    ++toread;
   }
   
  printf("ALIGNING\n");

	int **align=(int **)malloc(toread*sizeof(int *)); 
  for(k=0;k<toread-1;k++) { align[k]=(int *)malloc(nbrCA[0]*sizeof(int));}

	printf("In function\n");
	for(i=1;i<toread-1;++i) {
		
		
		int score = node_align(allpdb[0],nbrCA[0],allpdb[i],nbrCA[i],align[i]);
		int align_temp[5000];
		node_align(allpdb[i],nbrCA[i],allpdb[0],nbrCA[0],align_temp);
	
 		printf("I:%d RMSD:%8.5f Score: %d/%d\n",i,sqrt(rmsd_yes(allpdb[i],allpdb[0],nbrCA[i], align_temp,allall[i],nbrAT[i])),score,nbrCA[0]);
 	/*	char temp[50];
 		sprintf(temp,"strc_%d.pdb",i);
		write_strc(temp, allall[i],nbrAT[i],1);
		write_strc("strc_0.pdb", allall[0],nbrAT[0],1);*/
		
	}
	
	printf("Master ALIGNING\n");
	
	
	int master_align[nbrCA[0]];
	int index = 0;
	

	float **array=(float **)malloc(5000*sizeof(float *)); 
  for(k=0;k<5000;k++) { array[k]=(float *)malloc(toread*sizeof(float));}
  
  float **array_rev=(float **)malloc(5000*sizeof(float *)); 
  for(k=0;k<5000;k++) { array_rev[k]=(float *)malloc(toread*sizeof(float));}
	
	
	
for(i=0;i<nbrCA[0];++i) {
	//	printf("I:%d %3d%s",i,allpdb[0][i].res_number,allpdb[0][i].res_type);
		k=0;
		for (j=0;j<toread-1;++j) {
			//printf("::%3d%s ",allpdb[j+1][align[j][i]].res_number,allpdb[j+1][align[j][i]].res_type);
			if (align[j][i] == -1) {++k;}
		}	
		//printf("K:%d\n",k);
		// Hard coding the caca pour prendre data IVET
		/*if ( (allpdb[0][i].res_number > 4    && allpdb[0][i].res_number < 32 ) ||
		     (allpdb[0][i].res_number  > 35  && allpdb[0][i].res_number < 117) ||
		     (allpdb[0][i].res_number  > 120 && allpdb[0][i].res_number < 169) ||
		     (allpdb[0][i].res_number  > 184 && allpdb[0][i].res_number < 353)) {
		     	k = 0;
		     } else {
		     	k = 100;
		     }*/
		
		master_align[i] = k;
		
		

		     
		if (k < (toread-1) * 0.1) {
			for(j=0;j<toread-1;++j) {
				//printf("J:%d -> %f %f %f\n",j,allpdb[j][i].x_cord,allpdb[j][i].y_cord,allpdb[j][i].z_cord);
				if (j == 0) {
					array[index+0][j] = allpdb[0][i].x_cord;
					array[index+1][j] = allpdb[0][i].y_cord;
					array[index+2][j] = allpdb[0][i].z_cord;
					
					array_rev[j][index+0] = allpdb[0][i].x_cord;
					array_rev[j][index+1] = allpdb[0][i].y_cord;
					array_rev[j][index+2] = allpdb[0][i].z_cord;
					
				} else {
					if (align[j][i] == -1) {
						//printf("OFF\n");
						array[index+0][j] = -9999.9;
						array[index+1][j] = -9999.9;
						array[index+2][j] = -9999.9;
						
						array_rev[j][index+0] = -9999.9;
						array_rev[j][index+1] = -9999.9;
						array_rev[j][index+2] = -9999.9;
						
					} else {
						
						array[index+0][j] = allpdb[j][align[j][i]].x_cord;
						array[index+1][j] = allpdb[j][align[j][i]].y_cord;
						array[index+2][j] = allpdb[j][align[j][i]].z_cord;
						
						array_rev[j][index+0] = allpdb[j][align[j][i]].x_cord;
						array_rev[j][index+1] = allpdb[j][align[j][i]].y_cord;
						array_rev[j][index+2] = allpdb[j][align[j][i]].z_cord;
					}
				}
			  //	printf("ALign:%d:: %d ==  align[%d][%d] ---- array[%d][%d] = %f\n",i,align[j][i],j,i,index,j,array[index+0][j]);
			}
			
			index += 3;
		}
		
	}
	
	
	printf("Build Covariance matrix with index of:%d\n",index);
	
	gsl_matrix *cova = gsl_matrix_alloc(index,index);
	for(i=0;i<index;++i) {
		for (j=0;j<index;++j) {
			// Storer dans covariance matrix au point i,j  la covariance de array[i] qui contient autant d élément que de prot
			float myVal = covar_vince(array[i],array[j],toread-1);
			gsl_matrix_set(cova,i,j,myVal);
			//printf("I:%d J:%d Covar:%f\n",i,j,covar(array[i],array[j],argc-1));
						
		}
	}
	
	
	printf("Diagonalyse Matrix\n");
	
	//Diagonalyse the matrix
	
	gsl_vector *eval_temp = gsl_vector_alloc(index); /*Déclare un vector qui va contenir les eigenvalue */	
	gsl_matrix *evec_temp = gsl_matrix_alloc (index,index); /*Déclare une matrice qui va contenir les eigenvector correspondant à eval*/
	diagonalyse_matrix(cova,index,eval_temp,evec_temp); /*Diagonalyse la matrice*/
	gsl_eigen_symmv_sort (eval_temp, evec_temp,GSL_EIGEN_SORT_ABS_DESC);
	// Rebuild the eigenvector for strc one... will put 0 to Node that were not match
	
	gsl_vector *eval = gsl_vector_alloc(3*nbrCA[0]); /*Déclare un vector qui va contenir les eigenvalue */	
	gsl_matrix *evec = gsl_matrix_alloc (3*nbrCA[0],3*nbrCA[0]);
	gsl_matrix_set_all(evec,0);
	gsl_vector_set_all(eval,0);
	// J veux copien un eigen genre voir
	printf("Rebuild Eigenvector\n");
	
	float tot = 0.0000000000000000;
	
	for (i=0;i<index;++i) {
		tot += gsl_vector_get(eval_temp,i);
	
	}
	
 	FILE *file; /*Pointe le file défini par la fonction*/
 	file = fopen("project.dat","w"); /*Ouvre le fichier*/
	for (i=0;i<10;++i) {
		printf("I:%d %5.3f \n",i,gsl_vector_get(eval_temp,i)/tot*100);
		
		// Project on eigenvector
		
		for (k = 0 ;k<toread-1;++k) {
			fprintf(file,"	I:%d K:%d %f\n",i,k,vector_product(array_rev[k],evec_temp,i,index));
		}
		
	
	}
	fclose(file);
	for(j=0;j<index;++j) {
		if (j+6 > nbrCA[0]*3-2) {break;}
		k=0;
		gsl_vector_set(eval,j+6,gsl_vector_get(eval_temp,j));
		for(i=0;i<nbrCA[0];++i) {
			//printf("I:%d Master Align:%d < %f ?\n",i,master_align[i] , (argc-1) * 0.1);
			if (master_align[i] < (toread-1) * 0.1) {
				gsl_matrix_set(evec,i*3+0,j+6,gsl_matrix_get(evec_temp,k*3+0,j));
				gsl_matrix_set(evec,i*3+1,j+6,gsl_matrix_get(evec_temp,k*3+1,j));
				gsl_matrix_set(evec,i*3+2,j+6,gsl_matrix_get(evec_temp,k*3+2,j));
				++k;
			}
		}
	}
	write_eigen("pca_eigen.dat",evec,eval,3*nbrCA[0]);
	
	return(0);
}









float covar_vince(float *a,float *b,int atom) {
	
	float moya = 0.0;
	float moyb = 0.0;
	
	int i;
	int count = 0;
	// Moyenne de a et b
	
	for(i=0;i<atom;++i) {
		
		if (a[i] < -9999.0 || b[i] < -9999.0) {continue;}
		++count;
		moya += a[i];
		moyb += b[i];
    //printf("	%f::%f\n",a[i],b[i]);
	}
	moya /= count;
	moyb /= count;
  //printf("Moy:%f %f\n",moya,moyb);
	float total = 0.0;
	float sum = 0.0;
	for(i=0;i<atom;++i) {
		float temp = (a[i]-moya)*(b[i]-moyb);
		//printf("TEMP:%f\n",temp);
		sum = sum + temp;
		
	}
	//printf("THIS IS SUM:%f\n",sum);
	total = sum;
	return(total);
}
