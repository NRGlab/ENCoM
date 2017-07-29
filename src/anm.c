#include "STeM.h"

int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/

	int i,j;
	int help_flag = 1;
 	char file_name[500];

 	int lig = 0;

 	int verbose = 0;
	int nconn;

	char eigen_name[500] = "eigen.dat";
	char hessian_name[500] = "hessian.dat";
	int hessian_flag = 0;
	float temperature = 310;
	float cutoff = 18.0;
	int pfanm = 0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {++help_flag;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
 		if (strcmp("-pfanm",argv[i]) == 0) {pfanm= 1;}
 		if (strcmp("-cut",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);cutoff = temp;} 
 		
 		if (strcmp("-o",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-hes",argv[i]) == 0) {strcpy(hessian_name,argv[i+1]);++hessian_flag;}
 		
 		
 	}
 	
 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
  
  if (verbose == 1) {printf("Filename:%s\n",file_name);} 	
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
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
		
	// Pour les besoins... on limite à 800 atomes
	
	if (atom > 2000) {
		printf("Too much atom, if you want to remove the limite.... vincent.frappier@usherbrooke.ca\n");
		return(1);
	}
	
	
	//***************************************************
 	//*													*
 	//*Build Hessian									*
 	//*													*
 	//***************************************************


    double **hessian=(double **)malloc(3*atom*sizeof(double *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
    for(i=0;i<3*atom;i++) { hessian[i]=(double *)malloc(3*atom*sizeof(double));}
    for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
        
	gsl_matrix *h_matrix = gsl_matrix_alloc(3*atom, 3*atom);	/*Déclare une matrice hessian matrix de 3n par 3n*/
 	
 	//***************************************************
 	//*													*
 	//*Build template									*
 	//*													*
 	//***************************************************
	
	if (verbose == 1) {printf("Building Hessian\n");}
	

	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
	gsl_matrix *templaate = gsl_matrix_alloc(3*atom, 3*atom);
	gsl_matrix_set_all(templaate,1);
	if (false) {load_matrix(templaate, "template_name");} // matrice n X n séparée par des tabs
	build_enm(strc_node,hessian,atom,templaate,cutoff,pfanm);	// modidfier templaate, mettre cutoff énorme, pfanm à 0
 	if (verbose == 1) {printf("\tAssigning Array\n");}	
	assignArray(h_matrix,hessian,3*atom,3*atom);



	//***************************************************
	//*													*
	//*Diagonalyse the matrix							*
	//*													*
	//***************************************************
	

	if (verbose == 1) {printf("Diagonalizing Hessian\n");}
	gsl_vector *eval = gsl_vector_alloc(3*atom); /*Déclare un vector qui va contenir les eigenvalue */	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom); /*Déclare une matrice qui va contenir les eigenvector correspondant à eval*/
	diagonalyse_matrix(h_matrix,3*atom,eval,evec); /*Diagonalyse la matrice*/

	if (verbose == 1) {
		printf("First eigenvalues\n");
		for (i=0;i<10;++i) {printf("I:%d %.10f\n",i,gsl_vector_get(eval,i));}
	}
 	
 	write_eigen(eigen_name,evec,eval,3*atom);


		 	
 	float ener = calc_energy(atom,eval,temperature,6);
 	printf("Energy:%10.8f\n",ener);
 	printf("Energy/node:%10.8f\n",ener/(float)(atom*3));
 	// addition pour printer la correlation
 	gsl_matrix *k_inverse = gsl_matrix_alloc(atom, atom);
	k_inverse_matrix_stem(k_inverse,atom,eval,evec,6,atom*3-6);

	printf("Correlation:%f\n",correlate(k_inverse,strc_node, atom));
	gsl_matrix_free(k_inverse);
 
 	
 	gsl_matrix_free(templaate);
 	
 	gsl_matrix_free(h_matrix);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
 	
 	free(connect_h);

 	
 	
 	
 	
 	return(1);
 	
}
