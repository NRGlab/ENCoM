#include "STeM.h"

void read_vcon(char filename[100],int atom, struct pdb_atom *strc, gsl_matrix *m) {
	FILE *file;
 	file = fopen(filename,"r");
 	int actual;
 	int second;
 	float value;
 	int i = 0;
 	int k;
 	char line[82],buffer[82];
 	int pair[atom];
 	char *line_ptr;
 	gsl_matrix_set_all(m,0);
 	while(fgets(line,82,file)) {
 		strcpy(buffer,line);
 		buffer[6] = '\0';
 		if (sscanf(buffer,"%d",&actual) == 1) {pair[i] = actual;++i;continue;}	
	}
	printf("V con atom:%d\n",i);
	fclose(file);
	i=-1;
	file = fopen(filename,"r");
	while(fgets(line,82,file)) {
 		strcpy(buffer,line);
 		buffer[6] = '\0';
 		if (sscanf(buffer,"%d",&actual) == 1) {++i;continue;}
 		strcpy(buffer,line);
 		buffer[28] = '\0';
 		if (sscanf(buffer,"%d",&second) == 1) {
 			//printf("I:%d  Actual:%d Second:%d\t",i,actual,second);
 			for (k=0;k<atom;++k) {
 				if (pair[k] == second) {break;}
 			}
 			line_ptr = line;
 			//printf("%s",line_ptr+46);
 			sscanf(line_ptr+46,"%f",&value);
 			//if (strc[i].atom_type == 3 && strc[k].atom_type == 3) {printf("I:%d K:%d Value:%f \n",i,k,value);}
 			gsl_matrix_set(m,i,k,value);
 			
 		} else {continue;} 		
	}
	fclose(file);
	
}
void print_templaate(struct pdb_atom *newstrc,int atom,gsl_matrix *m,char filename[100],float min,float max) {
	FILE *out_file;
 	out_file = fopen(filename,"w"); 	
	int i,j,k;
	int temp;
	int count = 1;
	for (i=0;i<atom;++i) {
		temp = 0;
		for (j=i;j<atom;++j) {
			if (i==j) {continue;}
			temp = 0;
			for (k=0;k<6;++k) {
				if (newstrc[i].node_c[k] == newstrc[j].node) {++temp;} 
			}
			//printf("I:%d J:%d TEMP:%d\n",i,j,temp);
			if (temp != 0) {continue;}
			if ((gsl_matrix_get(m,i,j) < max) &&  (gsl_matrix_get(m,i,j) > min)) {
				fprintf(out_file,"HETATM");
				fprintf(out_file,"%5.d %s%s%6.d%12.3f%8.3f%8.3f  1.00 %5.2f\n",
		 			count,
		 			newstrc[i].atom_prot_type,
		 			newstrc[i].res_type,
		 			newstrc[i].res_number,
		 			newstrc[i].x_cord,
		 			newstrc[i].y_cord,
		 			newstrc[i].z_cord,
		 			gsl_matrix_get(m,i,j)/10
	 			);
	 			++count;
	 			fprintf(out_file,"HETATM");
				fprintf(out_file,"%5.d %s%s%6.d%12.3f%8.3f%8.3f  1.00 %5.2f\n",
		 			count,
		 			newstrc[j].atom_prot_type,
		 			newstrc[j].res_type,
		 			newstrc[j].res_number,
		 			newstrc[j].x_cord,
		 			newstrc[j].y_cord,
		 			newstrc[j].z_cord,
		 			gsl_matrix_get(m,j,i)/10
	 			);
	 			++count;
	 			
			}
		}
	}
	for (i=1;i<(count+1)/2;++i) {
		
		fprintf(out_file,"CONECT %4d %4d\n",2*i-1,2*i);
		fprintf(out_file,"CONECT %4d %4d\n",2*i,2*i-1);
	}
}
 double templaate_average(gsl_matrix *m,int atom) {
 	 int i,j;
	double sum = 0;
	 for(i=0; i < atom; ++i) {
	 	for(j=0; j < atom; ++j) {
	 		sum += gsl_matrix_get(m,i,j);
	 	}
	 }
 	return (sum/(i*j));
 }
void all_interaction(struct pdb_atom *strc,int atom,int res_n, gsl_matrix *ma, int flag,gsl_matrix *mvcon,gsl_matrix *mint,struct pdb_atom *ca_strc) {
 	int k,l,i,j,o;
 	float temp;
 	int nod[2][2];
 	float max_inter = 0;
 	float all_score = 0.0;
 	char strong[100];
 	// On veut scorer chaque node, une par appaort à l'autre, pour ce faire on chaque atom et fais la somme du score
 	// ma est la dite matrice
	//gsl_matrix_set_all(ma,0);
	
	//printf("In all interaction\n");
	
	gsl_matrix *conta = gsl_matrix_alloc(8,8);
	gsl_matrix_set_all(conta,0);
	
	/*for(k=0;k<8;++k) {
		for (l=0;l<8;++l) {
			printf("%10.5f ",gsl_matrix_get(mint,k,l));
		}
		printf("\n");
	}*/
 	for (k=0;k<atom;++k) {
 		for (l=0;l<atom;++l) { 

 			if (strc[k].type == 0 || strc[l].type == 0) {continue;}
 			if (gsl_matrix_get(mvcon,l,k) == 0) {continue;} // Si pas de contact on next	
 			temp = 0;
 			if (strc[k].node == strc[l].node) {continue;} // On next si même node !
 			//printf("L:%4d K:%4d Type:%d :: %d Surface:%f\n",l,k,strc[k].type,strc[l].type,gsl_matrix_get(mvcon,l,k));

 			// On next si same node or connect node connect to node
 			
 			if (strc[k].node > res_n-1 || strc[l].node > res_n-1) {
				
				continue;
			}
			// On veut nexter si dans coller
 			for (i=0;i<6;++i) {
 			//	printf("I:%d\n",i);
 				if (ca_strc[strc[k].node].node_c[i] == strc[l].node) {++temp;break;}
 				if (ca_strc[strc[l].node].node_c[i] == strc[k].node) {++temp;break;}
 				for (j=0;j<6;++j) {
 					//printf("Node:%d %d\n",strc[k].node,strc[l].node);
 					nod[0][0] = ca_strc[strc[k].node].node_c[i];
 					nod[0][1] = ca_strc[strc[l].node].node_c[i];
 					if (nod[0][0] < 0){} else  {if (ca_strc[nod[0][0]].node_c[j] == strc[l].node) {++temp;}}				
 					if (nod[0][1] < 0){} else  {if (ca_strc[nod[0][1]].node_c[j] == strc[k].node) {++temp;}}
 					if (nod[0][0] > res_n - 1) {
 					   // printf("0 0 Strange con:%d > %d\n",nod[0][0],res_n-1);
 					    nod[0][0] = -1;
 					}
					if (nod[0][1] > res_n - 1) {
					  //  printf("0 1 Strange con:%d > %d\n",nod[0][1],res_n-1);
					    nod[0][1] = -1;
					}
 					for (o=0;o<6;++o) {
 						if (nod[0][0] < 0 ) {nod[1][0] = -1;} else { nod[1][0] = ca_strc[nod[0][0]].node_c[j];}
 						if (nod[1][1] < 0) {nod[1][1] = -1;} else { nod[1][1] = ca_strc[nod[0][1]].node_c[j];}
 						if (nod[1][0] > res_n - 1) {
 						  //  printf("1 0 Strange con:%d > %d\n",nod[1][0],res_n);
 						    nod[1][0] = -1;
 						}
						if (nod[1][1] > res_n - 1) {
						  //  printf("1 1 Strange con:%d > %d\n",nod[1][1],res_n);
						    nod[1][1] = -1;
						}
						if (nod[1][0] < 0){} else  {if (ca_strc[nod[1][0]].node_c[j] == strc[l].node) {++temp;}}				
 						if (nod[1][1] < 0){} else  {if (ca_strc[nod[1][1]].node_c[j] == strc[k].node) {++temp;}}
 					
 					}
 				}
 			}
 			if (temp != 0) {continue;} // On next si node connecté
 			
 			
 			// Le score pour a adder
 			
			temp = gsl_matrix_get(mvcon,l,k)*gsl_matrix_get(mint,strc[k].type-1,strc[l].type-1);
			all_score += temp;
			gsl_matrix_set(ma,strc[k].node,strc[l].node,gsl_matrix_get(ma,strc[k].node,strc[l].node)+temp);
			//gsl_matrix_set(ma,strc[l].node,strc[k].node,gsl_matrix_get(ma,strc[k].node,strc[l].node));
			gsl_matrix_set(ma,strc[l].node,strc[k].node,gsl_matrix_get(ma,strc[l].node,strc[k].node)+temp);
		
			// Regarde strongest interaction !
			if (gsl_matrix_get(ma,strc[k].node,strc[l].node) > max_inter) {
				max_inter = gsl_matrix_get(ma,strc[k].node,strc[l].node);
				sprintf(strong,"Strongest interaction:%f Node:%d et %d, Amino:%s %d et %s %d\n",max_inter,strc[k].node,strc[l].node,strc[k].res_type,strc[k].res_number,strc[l].res_type,strc[l].res_number);
			}
		
			//printf("Node:(%4d,%4d) ResNumC: %s%d%s %s%d%sTemp:%10.5f Total:%10.5f %10.5f\n",strc[k].node,strc[l].node,strc[k].res_type,strc[k].res_number,strc[k].chain,strc[l].res_type,strc[l].res_number,strc[l].chain,temp,gsl_matrix_get(ma,strc[k].node,strc[l].node),gsl_matrix_get(ma,strc[l].node,strc[k].node));
			
			// Regarde les sufaces total en contact par type
			
			gsl_matrix_set(conta,strc[k].type-1,strc[l].type-1,gsl_matrix_get(conta,strc[k].type-1,strc[l].type-1)+gsl_matrix_get(mvcon,l,k));
			gsl_matrix_set(conta,strc[l].type-1,strc[k].type-1,gsl_matrix_get(conta,strc[l].type-1,strc[k].type-1)+gsl_matrix_get(mvcon,l,k));
 		}
 	}
 	printf("	%s",strong);
 /*	for(k=0;k<8;++k) {
 		for(l=0;l<8;++l) {
 			printf("%10.2f ",gsl_matrix_get(conta,k,l));
 		} 	
 		printf("\n");
 	}*/
 	
 	printf("	CF score:%f\n",all_score);
 	
 	gsl_matrix_free(conta);
 	
 }
 void all_interaction_leStatium(struct pdb_atom *strc,int atom,int res_n, gsl_matrix *ma, int flag,gsl_matrix *mvcon,struct pdb_atom *ca_strc) {
 	int k,l;
 	float temp;

 	// On veut scorer chaque node, une par appaort à l'autre, pour ce faire on chaque atom et fais la somme du score
 	// ma est la dite matrice
	//gsl_matrix_set_all(ma,0);
	
	printf("In all interaction\n");

 	for (k=0;k<atom;++k) {
 		for (l=0;l<atom;++l) { 

 			if (gsl_matrix_get(mvcon,l,k) == 0) {continue;} // Si pas de contact on next	
 			temp = 0;
 			//printf("L:%4d K:%4d Type:%d :: %d Surface:%f\n",l,k,strc[k].type,strc[l].type,gsl_matrix_get(mvcon,l,k));
 		 			
 			// Le score pour a adder
 			
			temp = gsl_matrix_get(mvcon,l,k);
			gsl_matrix_set(ma,strc[k].node,strc[l].node,gsl_matrix_get(ma,strc[k].node,strc[l].node)+temp);
			gsl_matrix_set(ma,strc[l].node,strc[k].node,gsl_matrix_get(ma,strc[l].node,strc[k].node)+temp);
		
 		}
 	}
 	
 	
 }
 void assign_lig_type(struct pdb_atom *strc, int atom, char inp[500]) {
 	FILE *file;
 	file = fopen(inp,"r");
 	char line[100];
 	int type[100];
 	int i;
 	for(i=0;i<100;++i) {type[i] = -1;}
 	i = 0;
 	while(fgets(line,82,file)) {
 		if (strncmp(line,"HET",3) == 0) {
 			int temp;
 			//printf("Line:%s",line);
 			sscanf(line+12,"%d",&temp);
 			type[i] = temp;
 			++i;
 		}
 		
	}
	i = 0;
	int k;
 	for (k =0;k<atom;++k) {
 		if (strc[k].atom_type != 1) {
 			if (type[i] == -1) {printf("Problem with atom type of ligand\n");continue;}
 			strc[k].type = type[i];
 			++i;
 		}
 	}
 	
 }
 void assign_atom_type(struct pdb_atom *strc, int atom) {
 	int k;
 	for (k =0;k<atom;++k) {
 		strc[k].type = 0;
 		if (strncmp(strc[k].atom_prot_type," H",2) == 0)	{continue;}
 		if (strncmp(strc[k].atom_prot_type,"H",1) == 0)	    {continue;}
 		if (strncmp(strc[k].atom_prot_type," OXT",4) == 0)	{strc[k].type = 2;}
 		if (strc[k].atom_type != 1) {
 			strc[k].type = 6;
 			if (strncmp(strc[k].atom_prot_type," N",2) == 0)	{strc[k].type = 3;}
 			if (strncmp(strc[k].atom_prot_type," O",2) == 0)	{strc[k].type = 1;}
 			if (strncmp(strc[k].atom_prot_type," HG",3) == 0)	{strc[k].type = 1;}
 			if (strncmp(strc[k].atom_prot_type," CL",3) == 0)	{strc[k].type = 4;}
 			if (strncmp(strc[k].atom_prot_type," BR",3) == 0)	{strc[k].type = 4;}
 			if (strncmp(strc[k].atom_prot_type," FE",3) == 0)	{strc[k].type = 1;}
 			if (strncmp(strc[k].atom_prot_type," II",3) == 0)	{strc[k].type = 4;}
 			//printf("Res:%s\tAtom Type:-%s-\tType:%d\n",strc[k].res_type,strc[k].atom_prot_type,strc[k].type);
 			
 			
 		}
		if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
		if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
		if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
		if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
		if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			
	 	if (strncmp(strc[k].res_type,"ALA",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
		}

		if (strncmp(strc[k].res_type,"ARG",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," NE",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CZ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," NH",3) == 0)	{strc[k].type = 3;}
		}

		if (strncmp(strc[k].res_type,"ASN",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," OD",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," ND",3) == 0)	{strc[k].type = 3;}
		}

		if (strncmp(strc[k].res_type,"ASP",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," OD",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," OD",3) == 0)	{strc[k].type = 2;}
		}

		if (strncmp(strc[k].res_type,"CYS",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," SG",3) == 0)	{strc[k].type = 6;}
		}

		if (strncmp(strc[k].res_type,"GLN",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," OE",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," NE",3) == 0)	{strc[k].type = 3;}
		}

		if (strncmp(strc[k].res_type,"GLU",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," OE",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," OE",3) == 0)	{strc[k].type = 2;}
		}

		if (strncmp(strc[k].res_type,"GLY",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
		}

		if (strncmp(strc[k].res_type,"HIS",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," ND",3) == 0)	{strc[k].type = 1;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CE",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," NE",3) == 0)	{strc[k].type = 1;}
		}

		if (strncmp(strc[k].res_type,"ILE",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 4;}
		}

		if (strncmp(strc[k].res_type,"LEU",3) == 0) {//printf("	Res:%s\tAtom Type:-%s-\tRes Num:%d\n",strc[l].res_type,strc[l].atom_prot_type,strc[l].res_number);
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 4;}
		}

		if (strncmp(strc[k].res_type,"LYS",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CE",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," NZ",3) == 0)	{strc[k].type = 3;}
		}

		if (strncmp(strc[k].res_type,"MET",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," SD",3) == 0)	{strc[k].type = 8;}
			if (strncmp(strc[k].atom_prot_type," CE",3) == 0)	{strc[k].type = 4;}
		}

		if (strncmp(strc[k].res_type,"PHE",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CE",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CZ",3) == 0)	{strc[k].type = 5;}
			
		}

		if (strncmp(strc[k].res_type,"PRO",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 4;}
		}

		if (strncmp(strc[k].res_type,"SER",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," OG",3) == 0)	{strc[k].type = 1;}
		}

		if (strncmp(strc[k].res_type,"THR",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," OG",3) == 0)	{strc[k].type = 1;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
		}

		if (strncmp(strc[k].res_type,"TRP",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," NE",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CE",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CZ",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CH",3) == 0)	{strc[k].type = 5;}
		}

		if (strncmp(strc[k].res_type,"TYR",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CD",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CE",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," CZ",3) == 0)	{strc[k].type = 5;}
			if (strncmp(strc[k].atom_prot_type," OH",3) == 0)	{strc[k].type = 1;}
		}

		if (strncmp(strc[k].res_type,"VAL",3) == 0) {
			if (strncmp(strc[k].atom_prot_type," N ",3) == 0)	{strc[k].type = 3;}
			if (strncmp(strc[k].atom_prot_type," CA",3) == 0)	{strc[k].type = 7;}
			if (strncmp(strc[k].atom_prot_type," C ",3) == 0)	{strc[k].type = 6;}
			if (strncmp(strc[k].atom_prot_type," O ",3) == 0)	{strc[k].type = 2;}
			if (strncmp(strc[k].atom_prot_type," CB",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
			if (strncmp(strc[k].atom_prot_type," CG",3) == 0)	{strc[k].type = 4;}
	 	}


	 	if (strc[k].atom_type == 4) { // Assignation des atom types pour l'ARN.
	 		//basés sur SOBOLEV ET AL., p.122, LIGIN classes

	 		//types spécifiques aux bases
	 		if (strncmp(strc[k].res_type,"  A",4) == 0) {
				if (strncmp(strc[k].atom_prot_type," N1 ",4) == 0)	{strc[k].type = 2;}
	 		}
	 		if (strncmp(strc[k].res_type,"  G",4) == 0) {
				if (strncmp(strc[k].atom_prot_type," N1 ",4) == 0)	{strc[k].type = 3;}
	 		}
	 		if (strncmp(strc[k].res_type,"  C",4) == 0) {
				if (strncmp(strc[k].atom_prot_type," N1 ",4) == 0)	{strc[k].type = 6;}
				if (strncmp(strc[k].atom_prot_type," N5 ",4) == 0)	{strc[k].type = 2;}
	 		}
	 		if (strncmp(strc[k].res_type,"  U",4) == 0) {
				if (strncmp(strc[k].atom_prot_type," N1 ",4) == 0)	{strc[k].type = 6;}
				if (strncmp(strc[k].atom_prot_type," N5 ",4) == 0)	{strc[k].type = 3;}
	 		}

	 		//Types généraux
	 		//Atomes de la base
	 		//carbones
	 		if (strncmp(strc[k].atom_prot_type," C2 ",4) == 0)	{strc[k].type = 5;}
	 		if (strncmp(strc[k].atom_prot_type," C3 ",4) == 0)	{strc[k].type = 5;}
	 		if (strncmp(strc[k].atom_prot_type," C4 ",4) == 0)	{strc[k].type = 5;}
	 		if (strncmp(strc[k].atom_prot_type," C5 ",4) == 0)	{strc[k].type = 5;}
	 		if (strncmp(strc[k].atom_prot_type," C6 ",4) == 0)	{strc[k].type = 5;}
	 		if (strncmp(strc[k].atom_prot_type," C8 ",4) == 0)	{strc[k].type = 5;}
	 		//azotes
	 		if (strncmp(strc[k].atom_prot_type," N2 ",4) == 0)	{strc[k].type = 3;}
	 		if (strncmp(strc[k].atom_prot_type," N3 ",4) == 0)	{strc[k].type = 2;}
	 		if (strncmp(strc[k].atom_prot_type," N4 ",4) == 0)	{strc[k].type = 3;}
	 		if (strncmp(strc[k].atom_prot_type," N6 ",4) == 0)	{strc[k].type = 3;}
	 		if (strncmp(strc[k].atom_prot_type," N7 ",4) == 0)	{strc[k].type = 2;}
	 		if (strncmp(strc[k].atom_prot_type," N9 ",4) == 0)	{strc[k].type = 6;}

	 		//Atomes du sucre
	 		//carbones
	 		if (strncmp(strc[k].atom_prot_type," C1'",4) == 0)	{strc[k].type = 8;}
	 		if (strncmp(strc[k].atom_prot_type," C2'",4) == 0)	{strc[k].type = 6;}
	 		if (strncmp(strc[k].atom_prot_type," C3'",4) == 0)	{strc[k].type = 8;}
	 		if (strncmp(strc[k].atom_prot_type," C4'",4) == 0)	{strc[k].type = 8;}
	 		if (strncmp(strc[k].atom_prot_type," C5'",4) == 0)	{strc[k].type = 8;}
	 		//oxygenes
	 		if (strncmp(strc[k].atom_prot_type," O4'",4) == 0)	{strc[k].type = 2;}
	 		if (strncmp(strc[k].atom_prot_type," O2'",4) == 0)	{strc[k].type = 1;}

	 		//Atomes du phosphate
	 		if (strncmp(strc[k].atom_prot_type," P  ",4) == 0)	{strc[k].type = 8;}
	 		if (strncmp(strc[k].atom_prot_type," OP1",4) == 0)	{strc[k].type = 2;}
	 		if (strncmp(strc[k].atom_prot_type," OP2",4) == 0)	{strc[k].type = 2;}
	 		if (strncmp(strc[k].atom_prot_type," OP3",4) == 0)	{strc[k].type = 2;}
	 		if (strncmp(strc[k].atom_prot_type," O5'",4) == 0)	{strc[k].type = 2;}
	 		if (strncmp(strc[k].atom_prot_type," O3'",4) == 0)	{strc[k].type = 2;}

	 	}
	 	//if (strc[k].type == 0) {printf("Res:%s\tAtom Type:-%s-\tType:%d\tK:%d\tAtom:%d\n",strc[k].res_type,strc[k].atom_prot_type,strc[k].type,k,strc[k].atom_number);}

	 	if (strc[k].type == 0) {printf("Res:%s\tAtom Type:-%s-\tType:%d\tK:%d\tAtom:%d\n",strc[k].res_type,strc[k].atom_prot_type,strc[k].type,k,strc[k].atom_number);}

	 }
 }
 
 
