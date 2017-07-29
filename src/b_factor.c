#include "STeM.h"

void assign_bfactor(struct pdb_atom *strc, gsl_matrix *m,int atom,int lig) {
	int i;
	float max = 0;
	float min = 100000;
	float sum = 0;
	float count = 0;
	for(i=0;i< (int) m->size1 ;++i) {
		sum += gsl_matrix_get(m,i,i);
		count += 1.0;
		if (gsl_matrix_get(m,i,i) < min) {min = gsl_matrix_get(m,i,i);}
		if (gsl_matrix_get(m,i,i) > max) {max = gsl_matrix_get(m,i,i);}
	}
	sum /= count;
	float ratio = (max-min)/90.0;
	
	printf("Average:%g Min:%g Max:%g\n",sum,min,max);
	for(i=0;i<atom;++i) {
		int node = strc[i].node;
		strc[i].b_factor = (gsl_matrix_get(m,node,node)-min)/ratio+5.0;
	}

}




int main(int argc, char *argv[]) {

	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
 	int help_flag = 1;
 	char file_name[500];
 	char eigen_name[500] = "eigen.dat";
 	char out_name[500] = "b_factor.pdb";
 	int verbose = 0;
	int mode = 6;
	int i;
	int toprint = 0;
	int nconn;
	int lig = 0;

 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);}
 		if (strcmp("-ieig",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}   
 		if (strcmp("-print",argv[i]) == 0) {toprint= 1;}   
 		if (strcmp("-m",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);mode = temp;} 
 	}
 	
 	if (help_flag == 0) { } else {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-o\tOutput Name Motion\n-ieig\tFile Name Eigen\n-v\tVerbose\n-sp\tSuper Node Mode (CA, N, C)\n-m\tMode\n-nm\tNombre de mode\n-lig\tTient compte duligand (sauf HOH)\n****************************\n");
 		return(0); 
 	} 

 	
 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
    
    assign_connect(file_name,connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom,all);}

	check_lig(strc_all,connect_h,nconn,all);

	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	//atom = count_atom_CA_n(strc_all,all,super_node,lig);
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	
	
	//***************************************************
 	//*													*
 	//* Load eigenvector et eigenvalue					*
 	//*													*
 	//***************************************************
 	
 	if (verbose == 1) {printf("Loading Eigenvector\n");}
 	
 	gsl_vector *eval = gsl_vector_alloc(3*atom);
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom);
	load_eigen(eval,evec,eigen_name,3*atom);
	
	if (verbose == 1) {printf("Inversing Matrix\n");}
	gsl_matrix *k_inverse = gsl_matrix_alloc(atom, atom); /*Déclare et crée une matrice qui va être le pseudo inverse*/
	if (mode < 0) {
		mode = atom*3+mode;
	}
	printf("Mode:%d\n",mode);
	k_inverse_matrix_stem(k_inverse,atom,eval,evec,mode,atom*3);

	assign_bfactor(strc_all,k_inverse,all,lig);
	
	write_strc(out_name, strc_all,all,1.0);
		
	gsl_matrix_free(k_inverse);
	
}
	
	
	
