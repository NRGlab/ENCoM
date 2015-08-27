#include "STeM.h"

void add_model(char filename[500],gsl_matrix *m,char matrix_name[500],float init_templaate,float vinit,float bond_factor,float angle_factor,float kp_factor,int lig,char inputname[500]);
void write_template(char filename[500],gsl_matrix *m,int nb_atom, struct pdb_atom *old);
void load_default_matrix(gsl_matrix *m);


int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/

	int i,j;
	int help_flag = 1;
 	char file_name[500][500];

 	int lig = 0;
 	char matrix_name[500] = "undef";
 	char param_name[500] = "param_m.dat";
 	char inputname[500] ="none";
 	int verbose = 0;
	int nconn;
	float vinit = 0.001; // Valeur de base
	float bond_factor = 1000;		// Facteur pour poid des bond strechcing
	float angle_factor = 10000;		// Facteur pour poid des angles
	double K_phi1 = 1;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5;
	float init_templaate = 1;
	float kp_factor = 10000;					// Facteur pour poid des angles dièdres
	int weight_factor = 0;
	char eigen_name[500] = "eigen.dat";
	char hessian_name[500] = "hessian.dat";
	char covariance_name[500] = "structure.cov";
	int hessian_flag = 0;
	float temperature = 310;
	int print_template = 0;
	int invert_hessian = 0;
	char ihess_name[500];
	int no_write = 0;
	int covariance_flag = 0;
	int total_model = 0;
	float noligand = -1;
	int starting_mode = 6;
	int only_bfact = 0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-il",argv[i]) == 0) {
 			while (strncmp(argv[i+1+total_model],".pdb",4)>0) {
 				//printf("%d I read:%s ,%d\n",total_model,argv[i+1+total_model],strncmp(argv[i+1+total_model],".pdb",4));
 				strcpy(file_name[total_model],argv[i+1+total_model]);
 				++total_model;
 				if (i+1+total_model == argc) {break;}
 				//if (total_model == 2) {break;}
 			}
 			--help_flag;
 		}
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name[total_model],argv[i+1]);++total_model;}
 		if (strcmp("-inp",argv[i]) == 0) {strcpy(inputname,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {++help_flag;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}
 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;}
 		if (strcmp("-m",argv[i]) == 0) {strcpy(matrix_name,argv[i+1]);}
 		if (strcmp("-sm",argv[i]) == 0) {int temp;sscanf(argv[i+1],"%d",&temp);starting_mode = temp;}
 		if (strcmp("-param",argv[i]) == 0) {strcpy(param_name,argv[i+1]);}
 		if (strcmp("-init",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);vinit = temp;}
 		if (strcmp("-t",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);init_templaate = temp;} 
 		if (strcmp("-kr",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);bond_factor = temp;}
 		if (strcmp("-kt",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);angle_factor = temp;}
 		if (strcmp("-kpf",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); kp_factor = temp;}
 		if (strcmp("-temp",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp); temperature = temp;}
 		if (strcmp("-w",argv[i]) == 0) {weight_factor = 1;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(eigen_name,argv[i+1]);}
 		if (strcmp("-hes",argv[i]) == 0) {strcpy(hessian_name,argv[i+1]);++hessian_flag;}
 		if (strcmp("-ihes",argv[i]) == 0) {strcpy(ihess_name,argv[i+1]);++invert_hessian;}
 		if (strcmp("-pt",argv[i]) == 0) {print_template = 1;}
 		if (strcmp("-b",argv[i]) == 0) {only_bfact = 1;}
 		if (strcmp("-no",argv[i]) == 0) {no_write = 1;}
 		if (strcmp("-cov",argv[i]) == 0)  {covariance_flag = 1;strcpy(covariance_name,argv[i+1]);}
 		if (strcmp("-fcov",argv[i]) == 0) {covariance_flag = 2;strcpy(covariance_name,argv[i+1]);}
 		if (strcmp("-nolig",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);noligand = temp;}
 		
 	}
 	int count = 0;
 	printf("\n***********ARNMajor*********\nFile:%s\nAlpha1:%f\nAlpha2:%f\nAlpha3:%f\nAlpha4:%f\nAlpha5:%f\n********************\n",
 	file_name[0],bond_factor,angle_factor,kp_factor,init_templaate,vinit);
 	
 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	if (verbose == 1) {printf("I do this file:%s\n",file_name[0]);}
 	
 	all = count_atom(file_name[0]);
 	
 	nconn = count_connect(file_name[0]);
 	
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
	for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
	
	assign_connect(file_name[0],connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name[0],strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom,all);}

	check_lig(strc_all,connect_h,nconn,all);

	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
		
	// Pour les besoins... on limite à 800 atomes
	
	if (atom > 10000) {
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
	
	/*Déclare une matrice hessian matrix de 3n par 3n*/
	gsl_matrix *h_matrix = gsl_matrix_alloc(3*atom, 3*atom);
 	//***************************************************
 	//*													*
 	//*Build template									*
 	//*													*
 	//***************************************************
	 	
	assign_atom_type(strc_all, all);
	if (strcmp(inputname,"none") == 0) {} else {assign_lig_type(strc_all, all, inputname);}
	gsl_matrix *vcon = gsl_matrix_alloc(all,all);
	gsl_matrix *inter_m = gsl_matrix_alloc(8,8);
	gsl_matrix *templaate = gsl_matrix_alloc(atom*3, atom*3);
	gsl_matrix_set_all(templaate,vinit);
	gsl_matrix_set_all(vcon,0);
	
	if (verbose == 1) {printf("Do Vcon !!!\n");}
	 	
	vcon_file_dom(strc_all,vcon,all);
	 	
	if (verbose == 1) {printf("Reading Interaction Matrix %s\n",matrix_name);}
	if (strcmp(matrix_name,"undef") == 0) {
		load_default_matrix(inter_m);
	} else {
		load_matrix(inter_m,matrix_name);
	}
	//print_matrix(inter_m);
	//write_matrix("vcon_vince.dat", vcon,all,all);
	if (verbose == 1) {printf("Building templaate\n");}
	 	
	all_interaction(strc_all,all, atom, templaate,lig,vcon,inter_m,strc_node);
	gsl_matrix_scale (templaate, init_templaate);
	
	if (print_template == 1) {write_template("template.dat", templaate,atom,strc_node);}
	
	if (verbose == 1) {printf("Building Hessian\n");}
	
	if (verbose == 1) {printf("\tCovalent Bond Potential\n");}		
	build_1st_matrix(strc_node,hessian,atom,bond_factor);
	
	if (verbose == 1) {printf("\tAngle Potential\n");}	
	build_2_matrix(strc_node,hessian,atom,angle_factor);
	
	if (verbose == 1) {printf("\tDihedral Potential\n");}	
	build_3_matrix(strc_node, hessian,atom,K_phi1/2+K_phi3*9/2,kp_factor);
	
	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
	build_4h_matrix(strc_node,hessian,atom,1.0,templaate);
	
 	if (verbose == 1) {printf("\tAssigning Array\n");}	
	assignArray(h_matrix,hessian,3*atom,3*atom);
	
	// S'il y a d'autres structure, les loader et et addiotionner à la matrix
	// Il faut passer touuuut les arguments !
	
	if (total_model != 1) {	
		printf("Total model:%d\n",total_model);
		for(i=1;i<total_model;++i) {
			printf("	I load:%s\n",file_name[i]);
			add_model(file_name[i],h_matrix,matrix_name,init_templaate,vinit,bond_factor,angle_factor,kp_factor,lig,inputname);
		//	printf("0,0 -> %f\n",gsl_matrix_get(h_matrix,0,0));
		}
	}
	
	gsl_matrix_scale(h_matrix, 1.0/(float)total_model);
	////////////////
	
	
	// Si on veut juste les valeurs de la prot sans le ligand
	
	if (noligand > 0) {
		// Combien de ligand node
		int atom_no_lig = 0;
		for(i=0;i<atom;++i) {
			if (strc_node[i].atom_type != 3 && strc_node[i].atom_type != 4) {
				++atom_no_lig;
				
				// On veut enlever un cutoff autour du ligand
				
				
				
				
			}
		}
		printf("Atom no lig:%d\n",atom_no_lig);
		
		gsl_matrix *h_matrix_no_lig = gsl_matrix_alloc(3*atom_no_lig, 3*atom_no_lig);
		int i_count = 0;
		int j_count = 0;
		int k,h;
		for (i=0;i<atom;++i) {
			j_count = 0;
			if (strc_node[i].atom_type == 3 || strc_node[i].atom_type == 4) {continue;}
			for (j=0;j<atom;++j) {
				if (strc_node[j].atom_type == 3 || strc_node[j].atom_type == 4) {continue;}
				//printf("I:%d (%d) J:%d (%d) -> Type:\n",i_count,i,j_count,j.);
				for(k=0;k<3;++k) {
					for (h=0;h<3;++h) {
						//printf("%d et %d -> (%d,%d) = (%d,%d)\n",i_count,j_count,3*i_count+k,3*j_count+h,3*i+k,3*j+h);
						gsl_matrix_set(h_matrix_no_lig,3*i_count+k,3*j_count+h,gsl_matrix_get(h_matrix,3*i+k,3*j+h));
					}
				}
				++j_count;
			}
			++i_count;
		}
		gsl_matrix_free(h_matrix);
		gsl_matrix *h_matrix = gsl_matrix_alloc(3*atom_no_lig, 3*atom_no_lig);
		gsl_matrix_memcpy (h_matrix, h_matrix_no_lig);
		gsl_matrix_free(h_matrix_no_lig);
		atom = atom_no_lig;
		if (verbose == 1) {printf("	Node:%d\n",atom);}
	}
	if (weight_factor == 1) {mass_weight_hessian(h_matrix,atom,strc_node);}
	if (hessian_flag != 0) {
		printf("Writing Hessian\n");
		write_matrix(hessian_name,h_matrix,3*atom,3*atom);
		write_matrix("template.dat",templaate,atom,atom);
	}

	


	//***************************************************
	//*													*
	//*Diagonalyse the matrix							*
	//*													*
	//***************************************************
	

	if (verbose == 1) {printf("Diagonalizing Hessian\n");}
	
	gsl_vector *eval = gsl_vector_alloc(3*atom); /*Déclare un vector qui va contenir les eigenvalue */	
	
	gsl_matrix *evec= gsl_matrix_alloc (3*atom,3*atom); /*Déclare une matrice qui va contenir les eigenvector correspondant à eval*/
	
	diagonalyse_matrix(h_matrix,3*atom,eval,evec); /*Diagonalyse la matrice*/
	
	if (weight_factor == 1) {adjust_weight_evec(evec,atom,strc_node);}
	
	if (verbose == 1)
	{
		printf("First eigenvalues\n");
		for (i=0;i<10;++i) {printf("I:%d %.10f\n",i,gsl_vector_get(eval,i));}
	}
	
	if (no_write == 0) {write_eigen(eigen_name,evec,eval,3*atom);}
	
	float ener  = 0.0;
	if (temperature  < 0) {
		float factor;
		for(factor = -2000.0;factor < 2000.0;++factor) {
			calc_energy(atom,eval,300.0-factor/1000.0,starting_mode);
		}
		for(factor = -300.0;factor < 300;++factor) {
			calc_energy(atom,eval,300.0-factor/10.0,starting_mode);
		}
	
	} else {
			ener = calc_energy(atom,eval,temperature,starting_mode);
			printf("Energy:%10.8f\n",ener);
			printf("Energy/node:%10.8f\n",ener/(float)(atom*3));
	}
	if(invert_hessian == 1)
	{
		gsl_matrix *k_totinv = gsl_matrix_alloc(3*atom, 3*atom);
		
		k_tot_inv_matrix_stem(k_totinv,atom,eval,evec,6,atom*3-6);
		
		printf("Writing matrix...\n");
		
		write_matrix(ihess_name,k_totinv,3*(atom-lig),3*(atom-lig));
	}

	if (covariance_flag > 0) 
	{
		gsl_matrix *k_inverse = gsl_matrix_alloc(atom, atom);
		k_inverse_matrix_stem(k_inverse,atom,eval,evec,starting_mode,atom*3-6);
	//	k_tot_inv_matrix_stem(k_inverse,atom,eval,evec,6,atom);
		FILE *file;
 		file = fopen(covariance_name,"w");
 			fprintf(file,"File:%s\n",file_name[0]);
			fprintf(file,"Energy:%f\n",ener);
   		for (i=0;i<atom;++i)
   		{
			for (j=i;j<atom;++j)
			{
				if (i != j && covariance_flag == 1) {continue;}
				fprintf(file,"COV: %d %d %s%d%s %s%d%s %f\n",i,j,strc_node[i].res_type,strc_node[i].res_number,strc_node[i].chain,
																	  strc_node[j].res_type,strc_node[j].res_number,strc_node[j].chain
				,gsl_matrix_get(k_inverse,i,j)*1000);
			}	
		}
		gsl_matrix_scale(k_inverse, 1000.0);
		printf("Correlation:%f\n",correlate(k_inverse,strc_node, atom));
		gsl_permutation * p = gsl_permutation_alloc (atom);
		int signum;
		gsl_linalg_LU_decomp (k_inverse, p,&signum);
		double det = gsl_linalg_LU_det (k_inverse, signum);
		printf("Det:%.10f\n",det);
		fprintf(file,"Det:%.10f\n",det);
		gsl_matrix_free(k_inverse);
		fclose(file);
	}
    if (only_bfact == 1) {
    
        gsl_matrix *k_inverse = gsl_matrix_alloc(atom, atom);
                k_inverse_matrix_stem(k_inverse,atom,eval,evec,6,atom*3-6);
                printf("Correlation:%f\n",correlate(k_inverse,strc_node, atom));
                        gsl_matrix_free(k_inverse);
                        
    }
	gsl_matrix_free(templaate);
	gsl_matrix_free(inter_m);
	gsl_matrix_free(vcon);
	
 	gsl_matrix_free(h_matrix);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
 	
 	free(connect_h);
	
 	return(1);
 	
}

void write_template(char filename[500],gsl_matrix *m,int nb_atom, struct pdb_atom *old) {
	FILE *out_file;
 	out_file = fopen(filename,"w");
	int i,j;
	for(i=0; i < nb_atom; ++i) {
		fprintf(out_file,"%s%d",old[i].res_type,old[i].res_number);
	 	for(j=0; j < nb_atom; ++j) {
	 		fprintf (out_file," %.6g",gsl_matrix_get (m, i, j));
	 	}
	 fprintf (out_file,"\n");
	 }
	 fclose(out_file);
}




void add_model(char filename[500],gsl_matrix *m,char matrix_name[500],float init_templaate,float vinit,float bond_factor,float angle_factor,float kp_factor,int lig,char inputname[500]) {
	int i,j;
	
	double K_phi1 = 1;				// Facteurs pour angles dièdres
	double K_phi3 = 0.5;
	int verbose = 1;
	
	int all = count_atom(filename);
	
 	int nconn = count_connect(filename);
 	printf("In function:%s\n",filename);
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
	for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(7*sizeof(int));}
	
	assign_connect(filename,connect_h);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	int atom = build_all_strc(filename,strc_all); // Retourne le nombre de Node
	
	if (verbose == 1) {printf("	Node:%d\n	Atom:%d\n",atom,all);}

	check_lig(strc_all,connect_h,nconn,all);

// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	double **hessian=(double **)malloc(3*atom*sizeof(double *)); // Matrix of the Hessian 1 2 3 (bond, angle, dihedral)
	for(i=0;i<3*atom;i++) { hessian[i]=(double *)malloc(3*atom*sizeof(double));}
	for(i=0;i<3*atom;i++)for(j=0;j<(3*atom);j++){hessian[i][j]=0;}
	
	/*Déclare une matrice hessian matrix de 3n par 3n*/
	gsl_matrix *h_matrix = gsl_matrix_alloc(3*atom, 3*atom);
	
	assign_atom_type(strc_all, all);
	if (strcmp(inputname,"none") == 0) {} else {assign_lig_type(strc_all, all, inputname);}
	gsl_matrix *vcon = gsl_matrix_alloc(all,all);
	gsl_matrix *inter_m = gsl_matrix_alloc(8,8);
	gsl_matrix *templaate = gsl_matrix_alloc(atom*3, atom*3);
	gsl_matrix_set_all(templaate,vinit);
	gsl_matrix_set_all(vcon,0);
	
	if (verbose == 1) {printf("Do Vcon !!!\n");}
	 	
	vcon_file_dom(strc_all,vcon,all);
	 	
	if (verbose == 1) {printf("Reading Interaction Matrix %s\n",matrix_name);}
	load_matrix(inter_m,matrix_name);
	//write_matrix("vcon_vince.dat", vcon,all,all);
	if (verbose == 1) {printf("Building templaate\n");}
	 	
	all_interaction(strc_all,all, atom, templaate,lig,vcon,inter_m,strc_node);
	gsl_matrix_scale (templaate, init_templaate);
	
	
	if (verbose == 1) {printf("Building Hessian\n");}
	
	if (verbose == 1) {printf("\tCovalent Bond Potential\n");}		
	build_1st_matrix(strc_node,hessian,atom,bond_factor);
	
	if (verbose == 1) {printf("\tAngle Potential\n");}	
	build_2_matrix(strc_node,hessian,atom,angle_factor);
	
	if (verbose == 1) {printf("\tDihedral Potential\n");}	
	build_3_matrix(strc_node, hessian,atom,K_phi1/2+K_phi3*9/2,kp_factor);
	
	if (verbose == 1) {printf("\tNon Local Interaction Potential\n");}	
	build_4h_matrix(strc_node,hessian,atom,1.0,templaate);
	
 	if (verbose == 1) {printf("\tAssigning Array\n");}	
	assignArray(h_matrix,hessian,3*atom,3*atom);
	gsl_matrix_add (m ,h_matrix);
}

void load_default_matrix(gsl_matrix *m) {
	gsl_matrix_set_all(m,3.0);
	
	gsl_matrix_set(m,3,0,1.0);
	
	gsl_matrix_set(m,1,1,1.0);
	gsl_matrix_set(m,3,1,1.0);
	gsl_matrix_set(m,7,1,1.0);


	gsl_matrix_set(m,2,2,1.0);
	gsl_matrix_set(m,3,2,1.0);
	gsl_matrix_set(m,6,2,1.0);
	
	gsl_matrix_set(m,0,3,1.0);
	gsl_matrix_set(m,1,3,1.0);
	gsl_matrix_set(m,2,3,1.0);
	
	
	gsl_matrix_set(m,2,6,1.0);
	gsl_matrix_set(m,6,6,1.0);
	
	gsl_matrix_set(m,1,7,1.0);
	gsl_matrix_set(m,7,7,1.0);

}









