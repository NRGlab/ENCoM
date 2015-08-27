#include "STeM.h"

float rmsd_hack(struct pdb_atom *init,struct pdb_atom *targ,int atom, int *align,struct pdb_atom *all_init,int all) {
 	// Fonction qui suimperpose deux structures en utilisant single value decomposition
 	// Centre les deux structures autour des atomes dans align et rotate init et init all pour fitter !
 	
 	int i,j,k; // Dummy
 	int t_atom;	
 	double cen_init[3],cen_targ[3];
 	double rmsd;
	int d;
	int temp;
	
 	gsl_matrix *init_v = gsl_matrix_alloc(atom,3);
 	gsl_matrix *targ_v = gsl_matrix_alloc(3,atom);
 	gsl_matrix *corr   = gsl_matrix_alloc(3,3);
 	gsl_vector *vec_s  = gsl_vector_alloc(3);
 	gsl_vector *vec_w  = gsl_vector_alloc(3);
 	gsl_matrix *mat_v  = gsl_matrix_alloc(3,3);
 	gsl_matrix *rota   = gsl_matrix_alloc(3,3);
	
	// Construit structure qui comprend seulement les atomes aligné
	
	struct pdb_atom t_init[atom];
	struct pdb_atom t_targ[atom];
	k=-1;
	for(i=0;i<atom;++i) {
 		d = align[i];
 		//printf("D:%d	%d\n",d,align[i]);
 		if(d == -1) {continue;}

 		//printf("D:%d	%d\n",i,align[i]);
 		//printf("(%8.5f,%8.5f,%8.5f) (%8.5f,%8.5f,%8.5f)\n",init[i].x_cord,init[i].y_cord,init[i].z_cord,targ[d].x_cord,targ[d].y_cord,targ[d].z_cord);
 		++k;
 		t_init[k].x_cord = init[i].x_cord;
 		t_init[k].y_cord = init[i].y_cord;
 		t_init[k].z_cord = init[i].z_cord;
 		
 		t_targ[k].x_cord = targ[d].x_cord;
 		t_targ[k].y_cord = targ[d].y_cord;
 		t_targ[k].z_cord = targ[d].z_cord;
 	}
 	t_atom =k+1;
 	// Centre les deux structures
 	if (k == -1) {return(-1);}
 	for(i=0;i<3;++i) {
 		cen_init[i] = 0;
 		cen_targ[i] = 0;
 	}
 	rmsd = 0;
 	// Trouve le centre de masse pour chaque structure
 	for(i=0;i<t_atom;++i) {
 		
 		cen_init[0] += t_init[i].x_cord;
 		cen_init[1] += t_init[i].y_cord;
 		cen_init[2] += t_init[i].z_cord;
 		
 		cen_targ[0] += t_targ[i].x_cord;
 		cen_targ[1] += t_targ[i].y_cord;
 		cen_targ[2] += t_targ[i].z_cord;
 		rmsd += ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 	}
 	
 	
 	for(i=0;i<3;++i) {
 		cen_init[i] /= t_atom;
 		cen_targ[i] /= t_atom;
 	}

	// Translationne et crée des matrice à multiplier (on retournera au centre tantot)
		
 	rmsd = 0;
 	for(i=0;i<t_atom;++i) {
 	 		

 		
 		t_init[i].x_cord -= cen_init[0];
 		t_init[i].y_cord -= cen_init[1];
 		t_init[i].z_cord -= cen_init[2];
 		
 		
 		t_targ[i].x_cord -= cen_targ[0];
 		t_targ[i].y_cord -= cen_targ[1];
 		t_targ[i].z_cord -= cen_targ[2];

 		gsl_matrix_set(init_v,i,0,t_init[i].x_cord);
 		gsl_matrix_set(init_v,i,1,t_init[i].y_cord);
 		gsl_matrix_set(init_v,i,2,t_init[i].z_cord);

 		gsl_matrix_set(targ_v,0,i,t_targ[i].x_cord);
 		gsl_matrix_set(targ_v,1,i,t_targ[i].y_cord);
 		gsl_matrix_set(targ_v,2,i,t_targ[i].z_cord);
 
 		rmsd += ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 	}
 	//printf("	Trans RMSD:%f\n",rmsd/t_atom);
 	
 	// Multiply les deux vecteurs qui comprennent les valeurs
 	 	
 	multiplie_matrix(targ_v,3,t_atom,init_v,t_atom,3,corr);


 	gsl_linalg_SV_decomp (corr, mat_v,vec_s,vec_w);
	//gsl_linalg_SV_decomp_jacobi (corr, mat_v, vec_s);
 	
 	// Multiplie les matrix
 	
	gsl_matrix_transpose (corr);
	
	multiplie_matrix(mat_v,3,3,corr,3,3,rota);
	//print_matrix(rota);
	
	// Centre avant de rotationer, les pdb qui n'ont pas été copier dans t_init
	
	for(i=0;i<atom;++i) {
		init[i].x_cord -= cen_init[0];
 		init[i].y_cord -= cen_init[1];
 		init[i].z_cord -= cen_init[2];
	}
	
	for(i=0;i<all;++i) {
		all_init[i].x_cord -= cen_init[0];
 		all_init[i].y_cord -= cen_init[1];
 		all_init[i].z_cord -= cen_init[2];
	}
	
	
	
 	rotate_all(rota,init,atom);
 	rotate_all(rota,t_init,t_atom);
 	rotate_all(rota,all_init,all);
 	
 	
 	rmsd = 0;
 	for(i=0;i<t_atom;++i) {
 		rmsd += ((t_init[i].x_cord-t_targ[i].x_cord)*(t_init[i].x_cord-t_targ[i].x_cord)+
 				 (t_init[i].y_cord-t_targ[i].y_cord)*(t_init[i].y_cord-t_targ[i].y_cord)+
 				 (t_init[i].z_cord-t_targ[i].z_cord)*(t_init[i].z_cord-t_targ[i].z_cord));
 	}
 		

 	
 	gsl_matrix_free(corr);
 	gsl_matrix_free(init_v);
 	gsl_matrix_free(targ_v);
 	gsl_matrix_free(mat_v);
    gsl_vector_free(vec_s);
  	gsl_vector_free(vec_w);
 	gsl_matrix_free(rota);

 	
 	return(rmsd/t_atom);
 }





int main(int argc, char *argv[]) {
	int all; /*Nombre d'atomes dans pdb*/
	int atom; /*Nombre de carbone CA*/
	int all_t; /*Nombre d'atomes dans pdb*/
	int atom_t; /*Nombre de carbone CA*/
 	int help_flag = 1;
	int aln_flag = 0;
 	char file_name[500];
 	char check_name[500];
 	char out_name[500];
 	int verbose = 0;
	int i,j;
	int lig = 0;
	int resnumc_flag = 0;
	int nconn;
	int print_flag = 0;
	int surf_flag = 0;
	float ligalign = 0; // Flag/valeur pour aligner seulement les résidus dans un cutoff du ligand, 0, one le fait pas... > 0... le cutoff
	int flag_ligalign = 0;
 	for (i = 1;i < argc;i++) {
 		if (strcmp("-i",argv[i]) == 0) {strcpy(file_name,argv[i+1]);--help_flag;}
 		if (strcmp("-h",argv[i]) == 0) {help_flag = 1;}
 		if (strcmp("-v",argv[i]) == 0) {verbose = 1;}

 		if (strcmp("-lig",argv[i]) == 0) {lig= 1;} 
 		if (strcmp("-t",argv[i]) == 0) {strcpy(check_name,argv[i+1]);help_flag = 0;}
 		if (strcmp("-o",argv[i]) == 0) {strcpy(out_name,argv[i+1]);++print_flag;}
 		if (strcmp("-ligc",argv[i]) == 0) {float temp;sscanf(argv[i+1],"%f",&temp);ligalign = temp;flag_ligalign = 1;lig = 1;}
 		if (strcmp("-resnumc",argv[i]) == 0) {resnumc_flag = 1;}
 		if (strcmp("-num",argv[i]) == 0) {resnumc_flag = 2;} 
 		if (strcmp("-surf",argv[i]) == 0) {surf_flag = 1;}
 		if (strcmp("-a",argv[i]) == 0) {aln_flag = 1;}
 	}
	 	
 	if (help_flag == 1) {
 		printf("****************************\nHelp Section\n-i\tFile Input (PDB)\n-v\tVerbose\n-w\tWeight Vector\n-t\tInitial Value of template (negative value for random)\n\tIf Load Template, multiply the template\n-lt\tLoad input template\n-sp\tSuper Node Mode (CA, N, C)\n-kt\tPoid de l'angle entre les nodes (1)\n-kr\tPoid de la distance entre les nodes (1)\n-f\tFile to fit\n****************************\n");
 		return(0); 
 	} 

 	//***************************************************
 	//*													*
 	//*Build a structure contaning information on the pdb
 	//*													*
 	//***************************************************
 	
 	// Première strucutre
 	
 	all = count_atom(file_name);
 	nconn = count_connect(file_name);
 	if (verbose == 1) {printf("First file:%s\n",file_name);}
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_h=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_h[i]=(int *)malloc(6*sizeof(int));}
    
    assign_connect(file_name,connect_h);

	// Assign tous les atoms
	
	struct pdb_atom strc_all[all];
	atom = build_all_strc(file_name,strc_all); // Retourne le nombre de Node
	if (atom > 2000) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	if (verbose == 1) {printf("	Atom:%d\n",all);}
	check_lig(strc_all,connect_h,nconn,all);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom);}
	struct pdb_atom strc_node[atom];
	atom = build_cord_CA(strc_all, strc_node,all,lig,connect_h,nconn);
	
	if (verbose == 1) {printf("	Assign Node:%d\n",atom);}
	
	// Free Connect
		
	//for(i=0;i<nconn;i++) {printf("I:%d\n",i);free(connect_h[i]);}
	//free(connect_h);
	
	//Construit la structure a comparer
 	nconn = 0;
 	all_t = count_atom(check_name);
 	nconn = count_connect(check_name);
 	if (verbose == 1) {printf("Sec file:%s\n",check_name);}
 	if (verbose == 1) {printf("Connect:%d\n",nconn);}
 	
	if (verbose == 1) {printf("Assigning Structure\n\tAll Atom\n");}
	
	// Array qui comprend tous les connects
	
	int **connect_t=(int **)malloc(nconn*sizeof(int *)); 
    for(i=0;i<nconn;i++) { connect_t[i]=(int *)malloc(6*sizeof(int));}

    assign_connect(check_name,connect_t);
	
	// Assign tous les atoms
	
	struct pdb_atom strc_all_t[all_t];
	atom_t = build_all_strc(check_name,strc_all_t); // Retourne le nombre de Node
	if (atom_t > 2000) {printf("Too much node.... To fix, ask vincent.frappier@usherbrooke.ca\n");return(1);}
	if (verbose == 1) {printf("	Atom:%d\n",all_t);}
	check_lig(strc_all_t,connect_t,nconn,all_t);
	
	// Assign les Nodes
	
	if (verbose == 1) {printf("	CA Structure\n");}
	
	if (verbose == 1) {printf("	Node:%d\n",atom_t);}
	struct pdb_atom strc_node_t[atom_t];

	atom_t = build_cord_CA(strc_all_t, strc_node_t,all_t,lig,connect_t,nconn);
	if (verbose == 1) {printf("	Assign Node:%d\n",atom_t);}
 
 	int align[atom];
 	int score = node_align(strc_node,atom,strc_node_t,atom_t,align);
 	printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	if (resnumc_flag == 2) {
 		// Align only based on number (pour GPCR renumber)
 		for(i=0;i<atom;++i) {align[i] = -1;}
 		score = 0;
 		for(i=0;i<atom;++i) {
 			if (strc_node[i].atom_type == 3) {continue;}
 			int j;
 			for(j = 0;j<atom_t;++j) {
 				if (strc_node_t[j].atom_type == 3) {continue;}
 				if (strc_node[i].res_number == strc_node_t[j].res_number) {align[i] = j;++score;break;}
 				
 			}
 		}
 		printf("Only num RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	}
 	if ((float)score/(float)atom < 0.8 && resnumc_flag == 0)
	{
 		printf("Low Score... Will try an homemade alignement !!!\n");
 		score = node_align_onechain(strc_node,atom,strc_node_t,atom_t,align);
 		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
 	}
 	
 	if ((float)score/(float)atom < 1.8 && resnumc_flag == 0)
	{
		printf("Low Score... Will try an homemade alignement !!!\n");
		score = node_align_low(strc_node,atom,strc_node_t,atom_t,align);
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
 	
 	if (flag_ligalign == 1) {
 		
		score = node_align_lig(strc_node,atom,strc_node_t,atom_t,align,strc_all,all,strc_all_t,all_t,ligalign);
		
		printf("RMSD:%8.5f Score: %d/%d\n",sqrt(rmsd_no(strc_node,strc_node_t,atom, align)),score,atom);
	}
	
	rmsd_yes(strc_node,strc_node_t,atom, align,strc_all,all);
	//rmsd_hack(strc_node_t,strc_node_t,atom, align,strc_all,all);
	if (print_flag != 0)
	{
		write_strc(out_name,strc_all,all,1.0);
	}
	
	if(aln_flag == 1)
	{
		for(i = 0; i < atom; i++)
		{
			if(align[i] != -1)
			{
				printf("%1i_%s %1i_%s\n", strc_node[i].res_number, strc_node[i].chain, strc_node_t[align[i]].res_number, strc_node_t[align[i]].chain);
			}
		}
	}
	
	
	if (surf_flag == 1) {
	
		// Regarde les surfaces en contact conserv
	
		gsl_matrix *inter_m = gsl_matrix_alloc(8,8);
		gsl_matrix_set_all(inter_m,1.0);
	
	
		// Init
	
		gsl_matrix *vcon = gsl_matrix_alloc(all,all);
		gsl_matrix *templaate = gsl_matrix_alloc(atom*3, atom*3);

		gsl_matrix_set_all(templaate,0.0);
		gsl_matrix_set_all(vcon,0.0);

		assign_atom_type(strc_all, all);
		vcon_file_dom(strc_all,vcon,all);
		all_interaction(strc_all,all, atom, templaate,lig,vcon,inter_m,strc_node);

		// Target
	
		gsl_matrix *vcon_t = gsl_matrix_alloc(all_t,all_t);
		gsl_matrix *templaate_t = gsl_matrix_alloc(atom_t*3, atom_t*3);
	
		gsl_matrix_set_all(templaate_t,0.0);
		gsl_matrix_set_all(vcon_t,0.0);

		assign_atom_type(strc_all_t, all_t);
		vcon_file_dom(strc_all_t,vcon_t,all_t);
		all_interaction(strc_all_t,all_t, atom_t, templaate_t,lig,vcon_t,inter_m,strc_node_t);
		float sum = 0.0;
		int count = 0;
	
		int tam[4];
		tam[0] = 0; // TP
		tam[1] = 0; // FP
		tam[2] = 0; // TN
		tam[3] = 0; // FN
	
		for(i=0;i<atom;++i) {
			if (align[i] == -1) {continue;}
			for(j=i;j<atom;++j) {
				if (abs(i-j) < 5) {continue;}
				if (align[j] == -1) {continue;}
				sum += pow(gsl_matrix_get(templaate,i,j) - gsl_matrix_get(templaate_t,align[i],align[j]),2);
				++count;
			
				if (gsl_matrix_get(templaate,i,j) > 10 && gsl_matrix_get(templaate_t,align[i],align[j]) > 10) {++tam[0];}
				if (gsl_matrix_get(templaate,i,j) < 10 && gsl_matrix_get(templaate_t,align[i],align[j]) > 10) {++tam[1];}
				if (gsl_matrix_get(templaate,i,j) < 10 && gsl_matrix_get(templaate_t,align[i],align[j]) < 10) {++tam[2];}
				if (gsl_matrix_get(templaate,i,j) > 10 && gsl_matrix_get(templaate_t,align[i],align[j]) < 10) {++tam[3];}
			
			}
		}
		printf("Vcon_RMSD=%f\n",sqrt(sum/count));
		printf("TP:%d\nFP:%d\nTN:%d\nFN:%d\n",tam[0],tam[1],tam[2],tam[3]);
	}
	return(1);
 	
 	
}
