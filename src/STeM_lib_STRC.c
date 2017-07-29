#include "STeM.h"
 
 
 int covalent_bond(int i, int j,int **con, int ncon) {
	int m,l;
	int con_flag = 0;
 	for (m=0;m<ncon;++m) {
 		if (j == con[m][0]) {
 			for(l=1;l<5;++l) {
 				
 				if (i == con[m][l]) { con_flag = 1; }
 			}
 		}
 	}
 	return(con_flag);
 }
int count_connect(char filename[100]) {
	FILE *file;
 	file = fopen(filename,"r");
 	if (file == NULL) {return(1);}
 	int count = 0; /*Count les lignes dans pdb*/
 	char line[100];

 	while(fgets(line,100,file)) {
 		if (strncmp("CONECT",line,4) == 0) {++count;}
	}
	
	fclose(file);
 	return(count);
}
int assign_connect(char filename[100],int **con) {

 	
	FILE *file;
 	file = fopen(filename,"r");
 	if (file == NULL) {return(1);}

 	char line[100];
 	int temp1,temp2,temp3,temp4,temp5,temp6;
 	int j=-1;
 	while(fgets(line,100,file)) {
 		if (strncmp("CONECT",line,6) == 0){
 			++j;
 			temp1 =-1;
 			temp2 =-1;
 			temp3 =-1;
 			temp4 =-1;
 			temp5 =-1;
 			temp6 =-1;
 			//printf("1:%d\t2:%d\t3:%d\t4:%d\n",temp1,temp2,temp3,temp4);
 			sscanf(line, "CONECT %d %d %d %d %d %d",&temp1,&temp2,&temp3,&temp4,&temp5,&temp6);
 			con[j][0] =temp1;
 			con[j][1] =temp2;
 			con[j][2] =temp3;
 			con[j][3] =temp4;
 			con[j][4] =temp5;
 			//con[j][5] =temp6;
 			//printf("J:%d Value:%d %d -- %d %d %d %d %d %d\n",j,temp6,con[j][5],temp1,temp2,temp3,temp4,temp5,temp6);
 			//printf("%s",line);
 		}
	}
	fclose(file);
	//printf("FCLOSE:%d\n",);

	return(j+1);
}

int check_covalent_CA(struct pdb_atom *CA,struct pdb_atom *strc,int atom,int all,int a,int b) {
 	int i,j;
 	float dist;
 	for(i=0;i<all;++i) {
 		if (strc[i].node != a) {continue;}
 		if ((strcmp(strc[i].atom_prot_type," N   ") == 0) || (strcmp(strc[i].atom_prot_type," C   ") == 0)) {} else {continue;}
 		for(j=0;j<all;++j) {
 			if (strc[j].node != b) {continue;}
 			if ((strcmp(strc[j].atom_prot_type," N   ") == 0) || (strcmp(strc[j].atom_prot_type," C   ") == 0)) {} else {continue;}
 			dist = (strc[i].x_cord - strc[j].x_cord)*(strc[i].x_cord - strc[j].x_cord)+(strc[i].y_cord - strc[j].y_cord)*(strc[i].y_cord - strc[j].y_cord)+(strc[i].z_cord - strc[j].z_cord)*(strc[i].z_cord - strc[j].z_cord);
 			//	printf("Node:%d %dAtom:-%s-	-%s-	Dist:%f	Res:%d	%d Node:%d %d\n",a,b,strc[i].atom_prot_type,strc[j].atom_prot_type,dist,strc[i].res_number,strc[j].res_number);
 			

 			if(2.9 > dist) {
 				
 				return(1);
 			}
 		}
 	}
 	return(0);
 }
 
 int check_covalent_DNA(struct pdb_atom *CA,struct pdb_atom *strc,int atom,int all,int a,int b) {
 	int i,j;
 	float dist = 0.0;
 	for(i=0;i<all;++i) {
 		//printf("Node:%d	Atom:-%s- Res:-%s-\n",strc[i].node,strc[i].atom_prot_type,strc[i].res_type);
 		if (strc[i].node != a) {continue;}
 		//printf("	Node:%d	Atom:-%s- Res:-%s-\n",strc[i].node,strc[i].atom_prot_type,strc[i].res_type);
 		//if ((strcmp(strc[i].atom_prot_type," P   ") == 0) || (strncmp(strc[i].atom_prot_type," O3'",4) == 0)) {} else {continue;}
 		//printf("		Node:%d	Atom:-%s- Res:-%s-\n",strc[i].node,strc[i].atom_prot_type,strc[i].res_type);
 		for(j=0;j<all;++j) {
 			if (strc[j].node != b) {continue;}

 			// Commenté pour pas skipper des pareils 
 			//if(strncmp(strc[j].atom_prot_type,strc[i].atom_prot_type,4) == 0) {continue;} 
 			//if ((strcmp(strc[j].atom_prot_type," P   ") == 0) || (strncmp(strc[j].atom_prot_type," O3'",4) == 0)) {} else {continue;}
 			dist = (strc[i].x_cord - strc[j].x_cord)*(strc[i].x_cord - strc[j].x_cord)+(strc[i].y_cord - strc[j].y_cord)*(strc[i].y_cord - strc[j].y_cord)+(strc[i].z_cord - strc[j].z_cord)*(strc[i].z_cord - strc[j].z_cord);
 			

 			if(6 > dist) {
 				//printf("Node:%d %dAtom:-%s-	-%s-	Dist:%f	Res:%d	%d\n",a,b,strc[i].atom_prot_type,strc[j].atom_prot_type,dist,strc[i].res_number,strc[j].res_number);
 				
 				return(1);
 			}
 		}

 	}
 	return(0);
 }
 
void check_lig(struct pdb_atom *strc,int **con,int ncon, int atom) {
	// strc => Structure de tous les atoms
	// con => Array qui contient tous les connections entre toutes les connects
	// ncon => Nombre de connect dans TOUS le pdb (problématique dans les nmr)
	// atom => Nombre d'atom dans la structure
	
	int i,k,l,n;
	int flag=0;
	int index;
	
	// On recommence, on cherche les atoms des CONNECT et vérifie s'il s'agit d'un atom type 3 (HETATM), si oui, ben on modifie pour 2 ou 5
		
	// On loop plusieurs fois, pour ne rien oublier, mais on ne passe jamais à travers car il y a un flag de quittage
	for (n=0;n<ncon;++n) {
		//printf("N:%d\n",n);
		int quit_f = 0; // Le flag pour quitter
		for (i=0;i<ncon;++i) {
			
			// Tous les atomes qui font des liens avec un Het ou qui sont un HET on leur CONNECT d'assigner, le premier chiffre est l'atome
			//printf("CON[%d][0] = %d\n",i,con[i][0]);
			flag = 0;
			for (l=0;l<atom;++l) {
				if (strc[l].atom_number == con[i][0]) {
					index = l;
					if (strc[l].atom_type == 1 || strc[l].atom_type == 4 || strc[l].atom_type == 2) {
						flag = 1;
						break;
					}
					continue;
				}			
			}
			if (flag == 1) {continue;}
			//printf("	Flag:%d\n",flag);
			// On regarde si ce connect qui est lié à un ligand est lié a d'autre atoms
			flag = 0;
			for(l=1;l<6;++l) {
			
				// Si égale à -1, c'est qu'il en a plus après
				if (con[i][l] == -1) {break;}
				
				//On cherche l'atome connecté et regarde son type
				for(k=0;k<atom;++k) {
					if(strc[k].atom_number == con[i][l]) {
						//printf("		I found the atom connect, is type is:%d et num:%d et connect num:%d\n",strc[k].atom_type,strc[k].atom_number,l);
						if (strc[k].atom_type == 2 || strc[k].atom_type == 1 || strc[k].atom_type == 4) {
							//printf("			Atom num %d will be change to cov and was type:%d\n",strc[index].atom_number,strc[index].atom_type);
							strc[index].atom_type = 2;
							flag = 1;
							++quit_f;
							break;
						}
					}
				}
				if (flag == 1) {break;}
			}
				
		}
		//printf("	Iteration:%d	Flag:%d\n",n,quit_f);
		//Quit après un loop si pas changement
		if (quit_f == 0) {
		//	printf("My flag:%d	It:%d/%d\n",quit_f,n,ncon);
			break;
		}
	}
	int prot_flag = 0;
	int dna_flag = 0;
	
	//Regarde Node Type (2,3 ou 5)
	
	
	for(i=0;i<atom;++i) {
		
		// Si jamais pas atom type 2... On skip	
		if (strc[i].atom_type == 1 || strc[i].atom_type == 3 || strc[i].atom_type == 4) {continue;}
		
		// Atom ds le backbone des prots
		if (strncmp(strc[i].atom_prot_type," CA ",4) == 0) {++prot_flag;}
		if (strncmp(strc[i].atom_prot_type," C  ",4) == 0) {++prot_flag;}
		if (strncmp(strc[i].atom_prot_type," N  ",4) == 0) {++prot_flag;}
		
		// Atom ds le backbone des ADN/RNA
		if (strncmp(strc[i].atom_prot_type," C1'",4) == 0) {++dna_flag;}
		if (strncmp(strc[i].atom_prot_type," C2 ",4) == 0) {++dna_flag;}
		if (strncmp(strc[i].atom_prot_type," P  ",4) == 0) {++dna_flag;}
		
		//printf("L:%d	Node:%d	Res:%s	Type:%d	Atom:-%s- et Prot:%d DNA:%d\n",i,strc[i].node,strc[i].res_type,strc[i].atom_type,strc[i].atom_prot_type,prot_flag,dna_flag);
		
		// Si On change de Node
		if (strc[i].node != strc[i+1].node || i == atom-1) {
			
			// Si on avait pas le Backbone de Prot, il s'agissait d'un ligand lié covalent
			if (prot_flag != 3 && dna_flag != 3) {
				for(l=0;l<atom;++l) {
					if(strc[l].node == strc[i].node) {strc[l].atom_type = 3;}
				}
				printf("Warning %s et %d et node %d was atom type 2, now is 3\n",strc[i].res_type,strc[i].res_number,strc[i].node);
			} else {
				//printf("Le Prot Flag était de:%d\n",prot_flag);
			}
			
			// Si on avait un backbone de DNA et qu'on pensait que c'était une prot
			
			if (dna_flag == 3) {
				for(l=0;l<atom;++l) {
					if(strc[l].node == strc[i].node) {strc[l].atom_type = 5;}
				}
				printf("Warning %s et %d et node %d was atom type 2, now is 5\n",strc[i].res_type,strc[i].res_number,strc[i].node);
			} else {
				//if (dna_flag > 0) {printf("Le flag DNA était de %d\n",dna_flag);} 
			}
			
			prot_flag = 0; dna_flag = 0;
		}
		
	}

	
}

void copy_strc(struct pdb_atom *target, struct pdb_atom *initial, int atom) {
        int i, j,k;
        for (i = 0; i < atom; i++) {
                target[i].atom_number = initial[i].atom_number;
                target[i].x_cord = initial[i].x_cord;
                target[i].y_cord = initial[i].y_cord;
                target[i].z_cord = initial[i].z_cord;
                strcpy(target[i].atom_prot_type, initial[i].atom_prot_type);
                target[i].res_number = initial[i].res_number;
                target[i].atom_type = initial[i].atom_type;
                strcpy(target[i].res_type, initial[i].res_type);
                target[i].b_factor = initial[i].b_factor;
                strcpy(target[i].chain, initial[i].chain);
                target[i].type = initial[i].type;
                target[i].node = initial[i].node;
                for (j = 0; j < 6; ++j) {
                        target[i].node_c[j] = initial[i].node_c[j];
                }
                for(j=0;j<3;++j) {
                	target[i].main_vars[j] = initial[i].main_vars[j];
                	for(k=0;k<3;++k) {
                		target[i].global_evecs[j][k] = initial[i].global_evecs[j][k];
                	}
                }
        }
}

int count_atom(char filename[100]) {
 	FILE *file;
 	file = fopen(filename,"r");
 	if (file == NULL) {return(0);}
 	int count = 0; /*Count les lignes dans pdb*/

 	char line[82];
 	while(fgets(line,82,file)) {
 				
 		// Il arrive parfois qu'il existe d'autre conformation d'acides aminées... On ne tient compte que des conformations A
 		
 		if ((line[16] == ' ') || (line[16] == 'A')) { } else {continue;}
 		
 		// On ignore les Hydrogènes (pas tenu compte ds Vcon)
 		
 		if (line[13] == 'H' || line[12] == 'H') { continue;}
 		
 		// On veut pas l'eau non plus (p-e a rajouter)
 		
 		if (line[17] == 'H' && 	line[18] == 'O' && line[19] == 'H') {continue;}
 		// On ne veut que des atomes ou hetatomes
 		
 		if ((strncmp("ATOM",line,4) == 0) ||(strncmp("HETATM",line,6) == 0) ) {} else {continue;}
 		++count;

	}
	fclose(file);
 return(count);
 }
 int build_cord_CA(struct pdb_atom *all, struct pdb_atom *CA, int atom,int ligand,int **con, int ncon) {
	int i,j,k,l,n;
	int lig_flag = -1;
	float cord[3];
	int lig_atom;
	int lig_tot;
	float dist;
	float min;
	int off = 0;
 	cord[0] = 0.0;
 	cord[1] = 0.0;
 	cord[2] = 0.0;
 	lig_tot = 0;
 	lig_atom = 0;
 	lig_flag = 0;
 	k = -1;
 	// On veut assigner les nodes (CA, P, atome centrale du node d'un ligand)
 	// On garde les nodes préablement assigner
 	
 	// Pour chaque atom on cherche les nodes
 	
 	for (i=0;i<atom;++i) {
 		
 		//printf("I:%d Node:%d	Atom:%d	Type:%s	Res num:%d	Res Type:%s Type:%d Flag:%d\n",i,all[i].node,all[i].atom_number,all[i].atom_prot_type,all[i].res_number,all[i].res_type,all[i].atom_type,lig_flag);		
 		
 		// Cas spécial des ligands (on prend la position de l'atome le plus au centre)
 		
 		if (lig_flag == 2) {
 			dist = (cord[0]- all[i].x_cord)*(cord[0]- all[i].x_cord)+(cord[1]- all[i].y_cord)*(cord[1]- all[i].y_cord)+(cord[2]- all[i].z_cord)+(cord[2]- all[i].z_cord);
 			if (i == lig_atom) {min = dist;}
 			if (dist < min) {min = dist;lig_atom = i;}
 			//printf("I:%d Dist:%8.5f Min:%8.5f Index:%d Node Now:%d Next:%d\n",i,dist,min,lig_atom,all[i].node,all[i+1].node);
 			
 		}
 		
 		
 		if(all[i].atom_type == 3 && ligand == 1 && lig_flag == 0) {lig_atom = i; lig_flag = 1;}
 		
 		if(all[i].atom_type == 3 && ligand == 1 && lig_flag == 1) {
 			
 			cord[0] += all[i].x_cord;
 			cord[1] += all[i].x_cord;
 			cord[2] += all[i].x_cord;
 			++lig_tot;
 		}
 		
 		if (lig_flag == 1 && (all[i].node != all[i+1].node || i == atom-1)) {
 			//printf("I CHANGE I\n");
 			cord[0] /= lig_tot;
 			cord[1] /= lig_tot;
 			cord[2] /= lig_tot;
 			lig_flag = 2;
 			i = lig_atom-1;
 			continue;
 		}
 		
 		if ((all[i].node != all[i+1].node || i == atom-1) && lig_flag == 2) {
 		//	printf("I CHANGE I ET FLAG\n");
 			lig_flag = 3;
 			lig_tot = i;
 			i = lig_atom;
 		}

 		// Si il s'agit d'un atome principal de Node ou bien atom centrale d'un ligand
 		if  (
 				(strncmp(all[i].atom_prot_type," CA ",4) == 0 && (all[i].atom_type == 2 || all[i].atom_type == 1)) || 
 				(strncmp(all[i].atom_prot_type," P  ",4) == 0 && (all[i].atom_type == 4 || all[i].atom_type == 5)) ||
 				(strncmp(all[i].atom_prot_type," C1'",4) == 0 && (all[i].atom_type == 4 || all[i].atom_type == 5)) ||
 				(strncmp(all[i].atom_prot_type," C2 ",4) == 0 && (all[i].atom_type == 4 || all[i].atom_type == 5)) ||

 				lig_flag == 3
 			) {} else {
 			//if (all[i].node < 40) {	printf("I:%d Node:%d	Atom:%d	Type:%s	Res num:%d	Res Type:%s Type:%d Flag:%d\n",i,all[i].node,all[i].atom_number,all[i].atom_prot_type,all[i].res_number,all[i].res_type,all[i].atom_type,lig_flag);}
 				continue;
 				
 				}
 		if (ligand == 0 && all[i].atom_type == 3) {continue;}

 		if (k > -1) {
 			if (CA[k].node ==  all[i].node) {continue;}
 		}
 		++k;
 		
 		//printf("K:%d Node:%d	Atom:%d	Type:%s	Res num:%d	Res Type:%s Con:%d %d %d Type:%d Lig:%d\n",k,all[i].node,all[i].atom_number,all[i].atom_prot_type,all[i].res_number,all[i].res_type,all[i].node_c[0],all[i].node_c[1],all[i].node_c[2],all[i].atom_type,ligand);
 		
 		
 		
 		// On copie l'information
 		for (j=0;j<6;++j) {CA[k].node_c[j] = -1;}
 		CA[k].node = k;
 		CA[k].mass = 0.0000;
 		CA[k].atom_number = all[i].atom_number;
	 	CA[k].res_number = all[i].res_number;
	 	CA[k].atom_type = all[i].atom_type;
	 	CA[k].b_factor = all[i].b_factor;

	 			
	 	CA[k].x_cord = all[i].x_cord;
	 	CA[k].y_cord = all[i].y_cord;
	 	CA[k].z_cord = all[i].z_cord;
	 		
	 	strcpy(CA[k].atom_prot_type,all[i].atom_prot_type);
	 	strcpy(CA[k].res_type,all[i].res_type);
	 	strcpy(CA[k].chain,all[i].chain);
	 	if (k+off != all[i].node) {
	 		printf("Node was off -> K:%d Node:%d	Atom:%d	Type:%s	Res num:%d	Res Type:%s Con:%d %d %d Type:%d Lig:%d\n",k,all[i].node,all[i].atom_number,all[i].atom_prot_type,all[i].res_number,all[i].res_type,all[i].node_c[0],all[i].node_c[1],all[i].node_c[2],all[i].atom_type,ligand);
	 		++off;
	 	}
	 	if (lig_flag == 3) {
	 		//printf("I JUST ADD LIG:%d\n",i);
	 		i = lig_tot;
	 		cord[0] = 0.0;
 			cord[1] = 0.0;
 			cord[2] = 0.0;
 			lig_tot = 0;
 			lig_atom = 0;
 			lig_flag = 0;
	 	}
 	}

 	// Check For Connection
 	int flag;
 	for (i=0;i<k+1;++i) {
 		
 		for (l=i+1;l<k+1;++l) {
 	
 		if (i == l) {continue;}
 		
		if (CA[i].atom_type == 1 || CA[i].atom_type == 2 || CA[i].atom_type == 4 || CA[i].atom_type == 5) {
			if (i+20 < l) {continue;}
			
			if ( CA[i].atom_type == 4 || CA[i].atom_type == 5) {
			    flag = check_covalent_DNA(CA,all,k,atom,CA[i].node,CA[l].node);
			} else {
			    flag = check_covalent_CA(CA,all,k,atom,CA[i].node,CA[l].node);
			    
			}
			if (flag != 1) {
					if (covalent_bond(CA[i].atom_number,CA[l].atom_number,con,ncon) == 1) {flag = 1;}
					
				
			}
			//printf("I TEST:%d   %d  Type:%d %d -> Flag:%d\n",i,l,CA[i].atom_type,CA[l].atom_type,flag);
			if (flag == 1) {
				for(n=0;n<6;++n) {
					if (CA[i].node_c[n] == -1) {CA[i].node_c[n] = CA[l].node;break;}
		 		}
		 		for(n=0;n<6;++n) {
					if (CA[l].node_c[n] == -1) {CA[l].node_c[n] = CA[i].node;break;}
		 		}
		 		continue; 
			 }
		}
 		if (CA[i].atom_type == 3) {
	 		if (covalent_bond(CA[i].atom_number,CA[l].atom_number,con,ncon) != 0) {
	 			for(n=0;n<6;++n) {
					if (CA[i].node_c[n] == -1) {CA[i].node_c[n] = CA[l].node;break;}
		 		}
		 		for(n=0;n<6;++n) {
					if (CA[l].node_c[n] == -1) {CA[l].node_c[n] = CA[i].node;break;}
		 		}
		 		continue;
	 		}
 		}
 	
 	}
 	}
 	
 	// Assign les masses aux nodes
 	
 	for (i=0;i<atom;++i) {
 		CA[all[i].node].mass += all[i].mass;
 	}
 	// for (i=0;i<k+1;++i) {
 	// 	printf("CA I:%d Type:%d Node:%d Atom:%d Type:%s Res num:%d %s Res Type:%s Con:%d %d %d %d %d %d Cord:%f,%f,%f Mass:%f\n",i,CA[i].atom_type,CA[i].node,CA[i].atom_number,CA[i].atom_prot_type,CA[i].res_number,CA[i].chain,CA[i].res_type,CA[i].node_c[0],CA[i].node_c[1],CA[i].node_c[2],CA[i].node_c[3],CA[i].node_c[4],CA[i].node_c[5],CA[i].x_cord,CA[i].y_cord,CA[i].z_cord,CA[i].mass);
 	// 	if (i != CA[i].node) {
 	// 		printf(" Node don't fit ID, I last\n");
 	// 		break;
 	// 	}
 	// }
 	return(k+1);
 
 }
 int assign_cord_all_atom(char filename[100],struct pdb_atom *structure) {
 	
 	float x,y,z,b_factor; 				/*Coordonné temporaire pour storer les coordonné*/
 	int at_nb,res; 						/*Temporaire Numéro d'atome et numéro de résidus*/
 	char line[82];						/*Store l'input d'une ligne de fichiers*/
 	unsigned int index = 0; 			/*Index pour storer les différente atomes*/
 	char *line_ptr; 					/*Pointeur pour la fonction sscanf des floating point*/
 	char buffer[6]; 					/*Buffer temporaire pour storer les string*/
 	int temp_count; 					/*Index pour storer les différente string*/
 	int n = 0; 								// Counter to assign node
 	//int node_flag = 0;
 	
 	// Ouvre le FILE
 	
 	FILE *file; /*Pointe le file défini par la fonction*/
 	file = fopen(filename,"r"); /*Ouvre le fichier*/
 	
 	
 	while(fgets(line,82,file)) {
 		line_ptr = line;
 		
 		// Il arrive parfois qu'il existe d'autre conformation d'acides aminées... On ne tient compte que des conformations A
 		
 		if ((line[16] == ' ') || (line[16] == 'A')) { } else {continue;}
 		
 		// On ignore les Hydrogènes (pas tenu compte ds Vcon)
 		
 		if (line[13] == 'H') { continue;}
 		
 		// On veut pas l'eau non plus (p-e a rajouter)
 		
 		if (line[17] == 'H' && 	line[18] == 'O' && line[19] == 'H') {continue;}
 		// On ne veut que des atomes ou hetatomes
 		
 		if ((strncmp("ATOM",line,4) == 0) ||(strncmp("HETATM",line,6) == 0) ) {} else {continue;}
			
 		// Assign l'atom Type, 1 = prot, 2 = prot unusual (mod aa), 3 = Lig, 4 = ADN/ARN, 5 = ADN,ARN modified
 			
 		if (strncmp("ATOM",line,4) == 0) {structure[index].atom_type = 1;}
 		if (strncmp("HETA",line,4) == 0) {structure[index].atom_type = 3;}
 		
 		// Si DNA/ARN, l'ACIDE AMINÉ VA ETRE SOIS -A  - ou -  A-
 		
 		if ((line[17] == 'A' || line[17] == 'U'|| line[17] == 'G'||line[17] == 'T' ||line[17] == 'C') && line[18] == ' ' && line[19] == ' ') {structure[index].atom_type = 4;}
 		if ((line[19] == 'A' || line[19] == 'U'|| line[19] == 'G'||line[19] == 'T' ||line[19] == 'C') && line[18] == ' ' && line[17] == ' ') {structure[index].atom_type = 4;}
 		
 		// Scan pour des chiffres de residus, b-factor,coord, atom 
			// Reset les valeurs	 		
 		
		x =-1;
		y=-1;
		z=-1;
		at_nb = -1;
		res = -1;
		b_factor = -1;
			
	 	sscanf((line_ptr+30),"%f %f %f",&x,&y,&z);
	 	sscanf((line_ptr+6),"%d",&at_nb);
	 	sscanf((line_ptr+23),"%d",&res);
	 	sscanf((line_ptr+60),"%f",&b_factor);	
	 	
	 	structure[index].x_cord = x;
	 	structure[index].y_cord = y;
	 	structure[index].z_cord = z;
	 	structure[index].atom_number = at_nb;
	 	structure[index].res_number = res;
	 	structure[index].b_factor = b_factor;
	 	structure[index].main_vars[0] = -1;
	 	structure[index].main_vars[1] = -1;
	 	structure[index].main_vars[2] = -1;
	 		
		// Copy le Residu Type

 		for (temp_count = 17; temp_count < 20;++temp_count) {
 			buffer[temp_count-17] = line[temp_count];
 			buffer[temp_count-16] = '\0';
 		}
 		strcpy(structure[index].res_type,buffer);
	 	
	 	//Copy la CHAIN
	 		
	 	for (temp_count = 21; temp_count < 22;++temp_count) {
	 		buffer[temp_count-21] = line[temp_count];
	 		buffer[temp_count-20] = '\0';
	 	}
	 	strcpy(structure[index].chain,buffer);
	 	
	 	//Copy le nom de l'atom
	 		
 		for (temp_count = 12; temp_count < 17;++temp_count) {
 			buffer[temp_count-12] = line[temp_count];
 			buffer[temp_count-11] = '\0';
 		}	 		
 		strcpy(structure[index].atom_prot_type,buffer);

		// On veut assigner des nodes (pour chaque residus)
  		if(index > 0) {
			if (structure[index-1].res_number != structure[index].res_number || strcmp(structure[index-1].res_type,structure[index-1].res_type) != 0) {++n;}
		}
		if (index > 0) {
		    if (structure[index-1].res_number == structure[index].res_number && line[29] != ' ')  {printf("%s\n",line);++n;}
		}
		structure[index].node = n;	
		/*if (strncmp(structure[index].atom_prot_type," P  ",4) == 0 && structure[index].atom_type == 4) {++node_flag;}
		if (strncmp(structure[index].atom_prot_type," CA ",4) == 0 && structure[index].atom_type != 4) {++node_flag;}
 			
 			 		
			if (node_flag != 0) {
			
				if 	((strcmp(structure[index].chain, structure[index-1].chain) != 0) || 
					 (structure[index].res_number != structure[index-1].res_number)  ||
					
					 ((structure[index].res_number == structure[index-1].res_number) && (strcmp(structure[index].res_type,structure[index-1].res_type) != 0))
					 )
				{
					if (node == 1 ||  structure[index].atom_type == 3) {++n;}
					if (node == 3 && structure[index].atom_type != 3) {++n;++n;++n;}
				}			
			
			}
 			structure[index].node = n;
 			if ((node == 3) && ((strncmp(structure[index].atom_prot_type," N  ",4) == 0))) {structure[index].node = n-1;}
 			if ((node == 3) && ((strncmp(structure[index].atom_prot_type," C  ",4) == 0))) {structure[index].node = n+1;}
 			if ((node == 3) && ((strncmp(structure[index].atom_prot_type," O  ",4) == 0))) {structure[index].node = n+1;}
 			//printf("Node:%d	Atom:-%s- Res:-%s- Num:%d\n",structure[index].node,structure[index].atom_prot_type,structure[index].res_type,structure[index].res_number);*/
 			printf("Index:%d Type:-%s- Res Num:%d	Node:%d Res:-%s- Atom type:%d Node:%4d\n",index,structure[index].atom_prot_type,structure[index].res_number,n,structure[index].res_type,structure[index].atom_type,structure[index].node);
 			++index;
			
	 	
 	}
 	fclose(file);
 	return(index);
 }
 
 int build_all_strc(char filename[100],struct pdb_atom *structure) {
 	
 	float x,y,z,b_factor; 				/*Coordonné temporaire pour storer les coordonné*/
 	int at_nb,res; 						/*Temporaire Numéro d'atome et numéro de résidus*/
 	char line[82];						/*Store l'input d'une ligne de fichiers*/
 	unsigned int index = 0; 			/*Index pour storer les différente atomes*/
 	char *line_ptr; 					/*Pointeur pour la fonction sscanf des floating point*/
 	char buffer[6]; 					/*Buffer temporaire pour storer les string*/
 	int temp_count; 					/*Index pour storer les différente string*/
 	int n = 0; 								// Counter to assign node
 	char pline[82]; // Previous Line
 	// Ouvre le FILE
 	
 	FILE *file; /*Pointe le file défini par la fonction*/
 	file = fopen(filename,"r"); /*Ouvre le fichier*/
 	if (file == NULL) {return(0);}
 	
 	while(fgets(line,82,file)) {
 		line_ptr = line;
 		
 		// Il arrive parfois qu'il existe d'autre conformation d'acides aminées... On ne tient compte que des conformations A
 		
 		if ((line[16] == ' ') || (line[16] == 'A')) { } else {continue;}
 		
 		// On ignore les Hydrogènes (pas tenu compte ds Vcon)
 		
 		
 		if (line[13] == 'H' || line[12] == 'H') { continue;}
 		
 		// On veut pas l'eau non plus (p-e a rajouter)
 		
 		if (line[17] == 'H' && 	line[18] == 'O' && line[19] == 'H') {continue;}
 		// On ne veut que des atomes ou hetatomes
 		
 		if ((strncmp("ATOM",line,4) == 0) ||(strncmp("HETATM",line,6) == 0) ) {} else {continue;}
		//printf("%s",line);
 		// Assign l'atom Type, 1 = prot, 2 = prot unusual (mod aa), 3 = Lig, 4 = ADN/ARN, 5 = ADN,ARN modified
 			
 		if (strncmp("ATOM",line,4) == 0) {structure[index].atom_type = 1;}
 		if (strncmp("HETA",line,4) == 0) {structure[index].atom_type = 3;}
 		
 		// Si DNA/ARN, l'ACIDE AMINÉ VA ETRE SOIS -A  - ou -  A-
 		
 		if ((line[17] == 'A' || line[17] == 'U'|| line[17] == 'G'||line[17] == 'T' ||line[17] == 'C') && line[18] == ' ' && line[19] == ' ') {structure[index].atom_type = 4;}
 		if ((line[19] == 'A' || line[19] == 'U'|| line[19] == 'G'||line[19] == 'T' ||line[19] == 'C') && line[18] == ' ' && line[17] == ' ') {structure[index].atom_type = 4;}
 		if ((line[19] == 'A' || line[19] == 'U'|| line[19] == 'G'||line[19] == 'T' ||line[19] == 'C') && line[18] == 'D' && line[17] == ' ') {structure[index].atom_type = 4;}
 		
 		// Scan pour des chiffres de residus, b-factor,coord, atom 
			// Reset les valeurs	 		
 		
		x =-1;
		y=-1;
		z=-1;
		at_nb = -1;
		res = -1;
		b_factor = 0;
			
	 	sscanf((line_ptr+30),"%f %f %f",&x,&y,&z);
	 	sscanf((line_ptr+6),"%d",&at_nb);
	 	sscanf((line_ptr+22),"%d",&res);
	 	sscanf((line_ptr+60),"%f",&b_factor);	
	 	
	 	structure[index].x_cord = x;
	 	structure[index].y_cord = y;
	 	structure[index].z_cord = z;
	 	structure[index].atom_number = at_nb;
	 	structure[index].res_number = res;
	 	structure[index].b_factor = b_factor;
	 		
	 	structure[index].main_vars[0] = -1;	
	 	structure[index].main_vars[1] = -1;	
	 	structure[index].main_vars[2] = -1;	
	 		
		// Copy le Residu Type

 		for (temp_count = 17; temp_count < 20;++temp_count) {
 			buffer[temp_count-17] = line[temp_count];
 			buffer[temp_count-16] = '\0';
 		}
 		strcpy(structure[index].res_type,buffer);
	 	
	 	//Copy la CHAIN
	 		
	 	for (temp_count = 21; temp_count < 22;++temp_count) {
	 		buffer[temp_count-21] = line[temp_count];
	 		buffer[temp_count-20] = '\0';
	 	}
	 	strcpy(structure[index].chain,buffer);
	 	
	 	//Copy le nom de l'atom
	 		
 		for (temp_count = 12; temp_count < 17;++temp_count) {
 			buffer[temp_count-12] = line[temp_count];
 			buffer[temp_count-11] = '\0';
 		}
 		buffer[4] = 	' ';	
 		strcpy(structure[index].atom_prot_type,buffer);
 		
 		// Assign la masse à l'atome
 		
 		structure[index].mass = 0.00;
		if (structure[index].atom_prot_type[1] == 'C') {structure[index].mass = 12.0107;}
  		if (structure[index].atom_prot_type[1] == 'N') {structure[index].mass = 14.0067;}
  		if (structure[index].atom_prot_type[1] == 'O') {structure[index].mass = 15.9994;}
  		if (structure[index].atom_prot_type[1] == 'S') {structure[index].mass = 32.065;}
  		if (structure[index].atom_prot_type[1] == 'E') {structure[index].mass = 55.845;}
  		
  		
  		
		// On veut assigner des nodes (pour chaque residus)
		
		if(index > 0) {
			// Si le numero de residus est différent
			if (structure[index-1].res_number != structure[index].res_number || strcmp(structure[index-1].res_type,structure[index].res_type) != 0 || strcmp(structure[index-1].chain,structure[index].chain) != 0) {++n;}
		  if (structure[index-1].res_number == structure[index].res_number && line[26] != ' ' && line[26] != pline[26])  {++n;}
		  
		}
		structure[index].node = n;
		if (structure[index].mass == 0.00) {
  		//	printf("Type:-%s- Res Num:%d	Node:%d Res:-%s- Atom type:%d Node:%4d Mass:%8f\n",structure[index].atom_prot_type,structure[index].res_number,n,structure[index].res_type,structure[index].atom_type,structure[index].node,structure[index].mass);	
  		}
  	//printf("Index:%d Type:-%s- Res Num:%d	Node:%d Res:-%s- Atom type:%d Node:%4d Mass:%8f Cord:%f %f %f\n",index,structure[index].atom_prot_type,structure[index].res_number,n,structure[index].res_type,structure[index].atom_type,structure[index].node,structure[index].mass,structure[index].x_cord,structure[index].y_cord,structure[index].z_cord);	
		//printf("%s N:%d\n",line,n);
		strcpy(pline,line);
 		++index;
			
	 	
 	}
 	fclose(file);
 	return(n+1); // Retourne le nombre de node
 }
 
 

