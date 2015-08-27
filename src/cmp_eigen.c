#include "STeM.h"

int count_eigen(char filename[100]) {
	int i = 0;
	int j,max=0;
	FILE *file;
 	file = fopen(filename,"r");
 	if (file == NULL) {return(1);}
 	char line[100];
	
 	while(fgets(line,100,file)) {
 	
 		
 		if (strncmp("2th",line,3) == 0) {break;}
 		sscanf(line,"%d\t%d",&i,&j);
 		if (i>max) {max=i;}
	}
	
	fclose(file);
	
	
	
	return(max);
}

float corr_eigen(gsl_matrix *evec,gsl_matrix *pca,int atom,int a,int b) {
	int i;
	double moy_x = 0.00000,moy_y = 0.00000;
 	double std_x = 0.00000,std_y=0.0000000;
 	double prod_tot=0.0000000000;
 	int count = 0;
	for (i=0;i<atom;++i) {
		//printf("I:%d	%f	%f\n",i,gsl_matrix_get(evec,i,a),gsl_matrix_get(pca,i,b));
		//printf("I:%d %8.5f::%8.5f\n",i,gsl_matrix_get(evec,i,b),gsl_matrix_get(pca,i,b));
		if (gsl_matrix_get(pca,i,b) == 0) {continue;}
		++count;
		moy_x += gsl_matrix_get(evec,i,a);
		moy_y += gsl_matrix_get(pca ,i,b);
	}
	moy_x /= count;
	moy_y /= count;
	
	for (i=0;i<atom;++i) {
		if (gsl_matrix_get(pca,i,b) == 0) {continue;}
		std_x += (moy_x - gsl_matrix_get(evec,i,a))*(moy_x - gsl_matrix_get(evec,i,a));
		std_y += (moy_y - gsl_matrix_get(pca ,i,b))*(moy_y - gsl_matrix_get(pca ,i,b));
		prod_tot += (gsl_matrix_get(evec,i,a)-moy_x)*(gsl_matrix_get(pca,i,b)-moy_y);
	}
	
	
	
	return(prod_tot/sqrt(std_y*std_x));

}

float over_eigen(gsl_matrix *evec,gsl_matrix *pca,int atom,int aa,int bb, struct pdb_atom *strc_node) {
	int i;
	int i_node;
	float a = 0.00000;
 	float x = 0.000000,z=0.00000;
 	float la=0.0000000000,lb=0.000000000;
 	// Caculate lenght (because of the vector reject)
 	
 	for(i=0;i<atom;++i) {
 		if (i % 3 == 0) {
 			i_node = i/3; // 3 eigenvectors par node
 			if (strc_node[i_node].atom_type == 4 && strncmp(strc_node[i_node].atom_prot_type," P  ",4) !=0) {
 				i += 2;
 				continue;
 			}
 		} 
 		if (gsl_matrix_get(pca,i,aa) == 0) {continue;}
 		//printf("I COMPARE: %f et %f\n",gsl_matrix_get(evec,i,aa),gsl_matrix_get(pca,i,bb));
 		la += gsl_matrix_get(evec,i,aa)*gsl_matrix_get(evec,i,aa);
 		lb += gsl_matrix_get(pca ,i,bb)*gsl_matrix_get(pca ,i,bb);
 		
 	
 	}
 	la = sqrt(la);
 	lb = sqrt(lb);
 	//printf("La:%f Lb:%f\n",la,lb);
 	for(i=0;i<atom;++i) {
 		if (i % 3 == 0) {
 			i_node = i/3; // 3 eigenvectors par node
 			if (strc_node[i_node].atom_type == 4 && strncmp(strc_node[i_node].atom_prot_type," P  ",4) !=0) {
 				i += 2;
 				continue;
 			}
 		} 
 		
// <<<<<<< HEAD
//	 		x = gsl_matrix_get(evec,i,aa)*la;
//	 		z = gsl_matrix_get(pca ,i,bb)*lb;
// =======
	 		 x = gsl_matrix_get(evec,i,aa);
	 		 z = gsl_matrix_get(pca ,i,bb);
// >>>>>>> master
	 		if (z == 0) {continue;}
			a += x*z;
	 	
 	}
 	
 	
 	return(a/(lb*la));

}

int main(int argc, char *argv[]) {
	int i,j;
	char eigen_name[100] = "eigen.dat";
	char pca_name[100] = "pca_eigen.dat";
	int atom[2];
	int mode = 6;
	int nm = 10;
	char filename[100] = "undefined";
	int verbose =1 ;
	int lig = 0;
	for (i = 1;i < argc;i++) {
		if (strcmp("-i",argv[i]) == 0) {strcpy(filename,argv[i+1]);}
		if (strcmp("-lig",argv[i]) == 0) {++lig;}
 		if (strcmp("-ieig1",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
		if (strcmp("-ieig2",argv[i]) == 0)  {strcpy(pca_name,argv[i+1]);}
		if (strcmp("-m",argv[i]) == 0)  {sscanf(argv[i+1],"%d",&mode);}
		if (strcmp("-nm",argv[i]) == 0)  {sscanf(argv[i+1],"%d",&nm);}
	}
	printf("-ieig1 %s -ieig2 %s -i %s\n",eigen_name,pca_name,filename);
	if (strcmp("undefined",filename) == 0) {
		printf("You need to give the WT form pdb\n");
		return(0);
	}
	
	int all = count_atom(filename);
 	int nconn = count_connect(filename);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
    
    assign_connect(filename,connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	int atom_strc = build_all_strc(filename,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom_strc,all);}

	check_lig(strc_all,connect_h,nconn,all);

	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	//atom = count_atom_CA_n(strc_all,all,super_node,lig);
	if (verbose == 1) {printf("	Node:%d\n",atom_strc);}
	struct pdb_atom strc_node[atom_strc];
	atom_strc = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_strc);}
	
	if (verbose == 1) {printf("Open eigen 1:%s\n",eigen_name);}
	
	atom[0] = count_eigen(eigen_name);
	
	gsl_vector *eval = gsl_vector_alloc(atom[0]);
	gsl_matrix *evec = gsl_matrix_alloc (atom[0],atom[0]);
	
	load_eigen(eval,evec,eigen_name,atom[0]);
	if (verbose == 1) {printf("Open eigen 2:%s\n",pca_name);}
	atom[1] = count_eigen(pca_name);
	
	gsl_vector *evalpca = gsl_vector_alloc(atom[1]);
	gsl_matrix *evecpca = gsl_matrix_alloc (atom[1],atom[1]);
	
	load_eigen(evalpca,evecpca,pca_name,atom[1]);
	
	printf("Atom:   %d::%d\n",atom[0],atom[1]);
	if (atom[0] != atom[1]) {printf("I exited, car pas même nombre de Ca entre les eigenvecoctors.... la première structure de build_pca, est la strc de référence pour comparer\n");return(1);}
	float corr;
	int k;
	printf("J:NMA I:PCA\n");
	//float overlap;
	for(i=mode-1;i<nm+mode-1;++i) {
		for(j=mode-1;j<nm+mode-1;++j) {
			corr = over_eigen(evec,evecpca,atom[0],i,j, strc_node);
			
			
				printf("I:%3d J:%3d %8.5f Val %8.5f :: %8.5f\n",i,j,corr*corr,gsl_vector_get(eval,i),gsl_vector_get(evalpca,j));
				
			
		}
	}
	return(1);
}
