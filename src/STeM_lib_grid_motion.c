#include "STeM.h"

void write_matrix_pdb(char filename[100], gsl_matrix *m,int nb_atom,int nb_atom_1) {
 	FILE *out_file;
 	out_file = fopen(filename,"w"); 	
	 int i,j;

	 for(i=0; i < nb_atom; ++i) {
	 fprintf(out_file,"HETATM %4d  C   GRD A   1      ",i+1);
	 	for(j=0; j < nb_atom_1; ++j) {
	 //	printf("-%f-\t",gsl_matrix_get(m,i,j));
	 	//HETATM 1776  C   GRD A   1      31.163  35.014   1.520  1.00  1.00 
	 		fprintf (out_file,"%6.3f  ",gsl_matrix_get (m, i, j));
	 	}
	 fprintf (out_file,"1.00  1.00\n");	
	// printf("\n");
	 }
	 fclose(out_file);
 }
 
 void load_eigen_motion(gsl_vector *v,gsl_matrix *m,char filename[100],int atom, int mode) {
	FILE *file;
 	file = fopen(filename,"r");
 	char line[10000];
 	float temp;
 	int i,j;
 	int actual_mode;
 	while(fgets(line,10000,file)) {
 		if (2 == sscanf(line,"%dth Eigenvalue:%f",&i,&temp)) {
 			gsl_vector_set(v,i-1,temp);
 			actual_mode = i;
 			//printf("Actual mode:%d\n",i);
 			if (mode+1 == i) {return;}
 		}
 		if (actual_mode == mode) {
 			if (3 == sscanf(line,"%d\t%d\t%f",&i,&j,&temp)) {
 				//printf("I:%d J:%d Temp:%f\n",i,j,temp);
 				gsl_matrix_set(m,i-1,j-1,temp);
 			}
 		} 		
 	}
 	fclose(file);
}
 
 void print_image_motion(struct pdb_atom *old,gsl_matrix *m,int mode,float max,float min,int atom,char out_name[100], int lig_f) {
 	FILE *out_file;
 	out_file = fopen(out_name,"w");
 	int l,k;
 	float pi = 3.14159265358979323846264338327950288419716939937510;
 	float temp;
 	float amplitude;
 	int node = 0;
 	struct pdb_atom newstrc[atom];
 	
 	for (l = 0; l<31;++l) {
 		temp = l;
 		amplitude = sin(temp*pi/15)*(max-min)/2+(max+min)/2;
 		//printf("Frame:%2d\tAmplitude:%f\n",l+1,amplitude);
 		fprintf(out_file,"MODEL     %4d\n",l+1);
		// Nous avons une nouvelle strucuture qui va être eigen bouger
		
	 	for (k=0;k<atom;++k) {
	 			 			
	 		newstrc[k].x_cord = old[k].x_cord;
	 		newstrc[k].y_cord = old[k].y_cord;
	 		newstrc[k].z_cord = old[k].z_cord;

	 		if ((old[k].atom_type == 3) && (lig_f == 0)) {continue;}
			
			node = old[k].node;
	 		newstrc[k].x_cord = old[k].x_cord+(amplitude*gsl_matrix_get(m,node*3,mode-1));
	 		newstrc[k].y_cord = old[k].y_cord+(amplitude*gsl_matrix_get(m,node*3+1,mode-1));
	 		newstrc[k].z_cord = old[k].z_cord+(amplitude*gsl_matrix_get(m,node*3+2,mode-1));
	 	}	 	
	 	
	 	for (k = 0;k<atom;++k) {
	 		
	 		if (strcmp(newstrc[k].res_type,"HOH") == 0) {continue;}
	 		if (old[k].atom_type == 1 || old[k].atom_type == 4 )  {fprintf(out_file,"ATOM  ");}
	 		if (old[k].atom_type == 3 || old[k].atom_type == 2 || old[k].atom_type == 5) {fprintf(out_file,"HETATM");}
	 		fprintf(out_file,"%5.d %s%s %s%4d%12.3f%8.3f%8.3f  1.00  %3.2f             \n",
	 			old[k].atom_number,
	 			old[k].atom_prot_type,
	 			old[k].res_type,
	 			old[k].chain,
	 			old[k].res_number,
	 			newstrc[k].x_cord,
	 			newstrc[k].y_cord,
	 			newstrc[k].z_cord,
	 			old[k].b_factor
	 			
	 			);
	 		//ATOM      1  N   ALA     2      -0.677  -1.230  -0.491  1.00  0.00           N //
 		}
 	fprintf(out_file,"TER\nENDMDL\n\n");
 	}
 	fclose(out_file);
 }
 
void print_grid_motion(struct pdb_atom *old,gsl_matrix *m,int mode,int atom,char out_name[100],gsl_matrix *g,int line_l,int nbr,int lig) {
 	FILE *out_file;
 	out_file = fopen(out_name,"w");
 	int h;
 	struct pdb_atom newstrc[atom];
 	int k,l;

 	
 	for (h=0;h<line_l;++h) {
 		fprintf(out_file,"Model %d\n",h+1);
 		//printf("Model %d\n",h+1);
	 	for (k=0;k<atom;++k) {
	 		
	 		newstrc[k].atom_number = old[k].atom_number;
	 		newstrc[k].res_number = old[k].res_number;
	 		newstrc[k].atom_type = old[k].atom_type;
	 		newstrc[k].b_factor = old[k].b_factor;

	 			
	 		newstrc[k].x_cord = old[k].x_cord;
	 		newstrc[k].y_cord = old[k].y_cord;
	 		newstrc[k].z_cord = old[k].z_cord;
	 		
	 		strcpy(newstrc[k].atom_prot_type,old[k].atom_prot_type);
	 		strcpy(newstrc[k].res_type,old[k].res_type);
	 		strcpy(newstrc[k].chain,old[k].chain);
	 		if (strcmp(old[k].res_type,"HOH") == 0) {continue;}
	 		if ((lig == 0) && (old[k].atom_type == 3)) {continue;}
	 		//printf("Atom:%d	Res Num:%d	Res:%s	Type:%s\n",old[k].atom_number,old[k].res_number,newstrc[k].res_type,newstrc[k].atom_prot_type);
			for (l=0;l<nbr;++l) {
				
		 		double amplitude = gsl_matrix_get(g,h,l);
		 		//printf("%d > %d ?\n",old[k].node*3,m->size1);
		 		if (old[k].node*3 > m->size1-1) {continue;}
		 	//	if (k == 10) {printf("L:%d H:%d K:%d Node:%d\n",l,h,k,old[k].node*3+2);}
		 		newstrc[k].x_cord += (amplitude*gsl_matrix_get(m,old[k].node*3,mode+l));
		 		newstrc[k].y_cord += (amplitude*gsl_matrix_get(m,old[k].node*3+1,mode+l));
		 		newstrc[k].z_cord += (amplitude*gsl_matrix_get(m,old[k].node*3+2,mode+l));
	 		
	 		}	 		
	 	}
	 	
	 	
	 	
	 	for (k = 0;k<atom;++k) {
	 		//printf("K:%d Res_num:%d et %4d\n",k,newstrc[k].res_number,newstrc[k].res_number);
	 		if (newstrc[k].atom_type == 1 || newstrc[k].atom_type == 4) {fprintf(out_file,"ATOM  ");}
	 		if (newstrc[k].atom_type == 2 || newstrc[k].atom_type == 3 || newstrc[k].atom_type == 5) {fprintf(out_file,"HETATM");}
	 		fprintf(out_file,"%5.d %s%s %s%4d%12.3f%8.3f%8.3f  1.00  0.00\n",
	 			newstrc[k].atom_number,
	 			newstrc[k].atom_prot_type,
	 			newstrc[k].res_type,
	 			newstrc[k].chain,
	 			newstrc[k].res_number,
	 			newstrc[k].x_cord,
	 			newstrc[k].y_cord,
	 			newstrc[k].z_cord
	 			
	 			);
	 		//ATOM      1  N   ALA     2      -0.677  -1.230  -0.491  1.00  0.00           N //
 		}
 		fprintf(out_file,"TER\nENDMDL\n\n");
 	}
 	fclose(out_file);
 }
  
 int test_grid(double **grid,int nbr_mode,int mode,gsl_matrix *matrix,struct pdb_atom *strc,int atom,float max_dist, int point) {
 	int i,j;
 	int l = 0;
 	
 	double actual[nbr_mode];
 	for(i=0;i<point;++i) {
 		for(j=0;j<nbr_mode;++j) {actual[j] = grid[i][j];}
 		
 		int ok = test_point_two(strc, actual, max_dist, atom,mode, nbr_mode,matrix);
 		if (ok == 1) {
 			++l; 
 			for(j=0;j<nbr_mode;++j) { grid[i][j] = actual[j];}
 		}
 	}
 	return(l); 
 }
 
 void find_max_ampli(struct pdb_atom *strc,gsl_matrix *m,int mode,int atom,float max,float *max_value,float *min_value) {
 	int i,j = mode;
 	double bx,by,bz;
 	double vx,vy,vz;
 	double a,b,c;
 	float max_temp;
 	float min_temp;
 	double dist;
	
 	max_temp = 99999;
 	min_temp = -99999;
 	for (i = 0;i<atom-1;++i) {
 		
 		int connect_flag = 0;
 		int t;
 		for (t=0;t<6;++t) {
 			if (strc[i].node_c[t] == strc[1+i].node) {connect_flag = 1;}
 		}
 		
 		if (connect_flag == 0) {printf("Covalent Break between node %d et %d\n",strc[i].node,strc[i+1].node);continue;}

 		
		bx = strc[i].x_cord - strc[i+1].x_cord;
		by = strc[i].y_cord - strc[i+1].y_cord;
		bz = strc[i].z_cord - strc[i+1].z_cord;
		dist = sqrt(bx*bx+by*by+bz*bz);

		vx = gsl_matrix_get(m,i*3,j-1) - gsl_matrix_get(m,(i+1)*3,j-1);
		vy = gsl_matrix_get(m,i*3+1,j-1) - gsl_matrix_get(m,(i+1)*3+1,j-1);
		vz = gsl_matrix_get(m,i*3+2,j-1) - gsl_matrix_get(m,(i+1)*3+2,j-1);
			
		a = vx*vx+vy*vy+vz*vz;
	 	b = 2*(bx*vx+by*vy+bz*vz);
	 	c = bx*bx+by*by+bz*bz-(dist*dist*max*max);
	 		 	
	 	if (max_temp > (-b + sqrt(b*b-4*a*c))/(2*a)) {max_temp =  (-b + sqrt(b*b-4*a*c))/(2*a);}
	 	if (min_temp < (-b - sqrt(b*b-4*a*c))/(2*a)) {min_temp =  (-b - sqrt(b*b-4*a*c))/(2*a);}
	 	//printf("Atom pair:%d\tMax:%g\tMin:%g\tDist:%f\n",i+1,(-b + sqrt(b*b-4*a*c))/(2*a),(-b - sqrt(b*b-4*a*c))/(2*a),dist);
	 }
	*max_value = max_temp; *min_value = min_temp;	 	
}
 
 
void find_max_ampli_two(struct pdb_atom *strc,gsl_matrix *m,int mode,int atom,float max,double *max_value,double *min_value,struct pdb_atom *ori) {
 	int i,j = mode,l;
 	double bx,by,bz;
 	double bxi,byi,bzi;
 	double vx,vy,vz;
 	double a,b,c;
 	double max_temp;
 	double min_temp;
 	double dist,disti;
	
 	max_temp = 99999;
 	min_temp = -99999;
 	for (i = 0;i<atom;++i) {
 		for (l = 0;l<atom;++l) {
	 		int connect_flag = 0;
	 		int t;
	 		for (t=0;t<6;++t) {
	 			if (strc[i].node_c[t] == strc[l].node) {connect_flag = 1;}
	 		}
	 		
	 		if (connect_flag == 0) {continue;}

	 		
			bx = strc[i].x_cord - strc[l].x_cord;
			by = strc[i].y_cord - strc[l].y_cord;
			bz = strc[i].z_cord - strc[l].z_cord;
		
			bxi = ori[i].x_cord - ori[l].x_cord;
			byi = ori[i].y_cord - ori[l].y_cord;
			bzi = ori[i].z_cord - ori[l].z_cord;
		
			disti = sqrt(bxi*bxi+byi*byi+bzi*bzi);
			dist = sqrt(bx*bx+by*by+bz*bz);
		
			vx = gsl_matrix_get(m,i*3,j-1) - gsl_matrix_get(m,l*3,j-1);
			vy = gsl_matrix_get(m,i*3+1,j-1) - gsl_matrix_get(m,l*3+1,j-1);
			vz = gsl_matrix_get(m,i*3+2,j-1) - gsl_matrix_get(m,l*3+2,j-1);
			
			a = vx*vx+vy*vy+vz*vz;
		 	b = 2*(bx*vx+by*vy+bz*vz);
		 	c = bx*bx+by*by+bz*bz-(disti*disti*max*max);
		 	//printf("sqrt(%f)\n",b*b-4*a*c);
		 	if (b*b-4*a*c < 0) {max_temp = 0;min_temp = 0;break;}	
		 	if (max_temp > (-b + sqrt(b*b-4*a*c))/(2*a)) {max_temp =  (-b + sqrt(b*b-4*a*c))/(2*a);}
		 	if (min_temp < (-b - sqrt(b*b-4*a*c))/(2*a)) {min_temp =  (-b - sqrt(b*b-4*a*c))/(2*a);}
		 	
		 	//printf("	Atom pair:%d\tMax:%g\tMin:%g\tDist:%f\tDisti:%f\n",i+1,(-b + sqrt(b*b-4*a*c))/(2*a),(-b - sqrt(b*b-4*a*c))/(2*a),dist,disti);
		 }
	 }
	*max_value = max_temp; *min_value = min_temp;	 	
}
  int test_point(struct pdb_atom *strc, double *actual, float limit, int atom,int mode, int nb_mode, gsl_matrix *matrix) {
 	
 	int j,k,l;
 	double bx,by,bz;
 	double dist,newstrc_dist,bxi,byi,bzi;
 	//printf("Max Dist:%f\n",limit);
 	for (k=0;k<atom;++k) {
 		for(l=k+1;l<atom;++l) {
			int connect_flag = 0;
	 		int t;
			for (t=0;t<6;++t) {
	 			if (strc[k].node_c[t] == strc[l].node) {++connect_flag;}
	 		}
	 		if (connect_flag != 0) {} else {continue;}
	 		
			bx = strc[k].x_cord - strc[l].x_cord;
			by = strc[k].y_cord - strc[l].y_cord;
			bz = strc[k].z_cord - strc[l].z_cord;
			dist = sqrt(bx*bx+by*by+bz*bz);
			//printf("\n	Atom:%d\t",k);
			bxi = bx;
			byi = by;
			bzi = bz;
				
			for (j=0;j<nb_mode;++j) {
				bxi += (actual[j]*gsl_matrix_get(matrix,k*3,j-1+mode))   - (actual[j]*gsl_matrix_get(matrix,l*3,j-1+mode));
				byi += (actual[j]*gsl_matrix_get(matrix,k*3+1,j-1+mode)) - (actual[j]*gsl_matrix_get(matrix,l*3+1,j-1+mode));
				bzi += (actual[j]*gsl_matrix_get(matrix,k*3+2,j-1+mode)) - (actual[j]*gsl_matrix_get(matrix,l*3+2,j-1+mode));			
			}
			newstrc_dist = sqrt(bxi*bxi+byi*byi+bzi*bzi);
			//printf("	Limit:%f	Dist:%f",dist*limit,newstrc_dist);		
			if (dist*limit >= newstrc_dist) {} else {
			//	printf("	Atom:%d %d Limit:[%f]	newstrc Dist:%f Ori:%f\n",k,l,dist*limit,newstrc_dist,dist);	
				return(0);}
			}
		}
		return(1);		
}

 int test_point_two(struct pdb_atom *strc, double *actual, float limit, int atom,int mode, int nb_mode, gsl_matrix *matrix) {
 	
 	int j,k,l,m;
 	double bx,by,bz;
 	double dist,newstrc_dist,bxi,byi,bzi;
 	//printf("Max Dist:%f\n",limit);
 //	for(i=0;i<nb_mode;++i) {printf("%.3f\t",actual[i]);}
 	//printf("\n");
 	for (k=0;k<atom;++k) {
 		for(l=k+1;l<atom;++l) {
			int connect_flag = 0;
	 		int t;
			for (t=0;t<6;++t) {
	 			if (strc[k].node_c[t] == strc[l].node) {++connect_flag;}
	 		}
	 		if (connect_flag != 0) {} else {continue;}
	 		
			bx = strc[k].x_cord - strc[l].x_cord;
			by = strc[k].y_cord - strc[l].y_cord;
			bz = strc[k].z_cord - strc[l].z_cord;
			dist = sqrt(bx*bx+by*by+bz*bz);
			//printf("\n	Atom:%d\t",k);
			bxi = bx;
			byi = by;
			bzi = bz;
				
			for (j=0;j<nb_mode;++j) {
				bxi += (actual[j]*gsl_matrix_get(matrix,k*3,j-1+mode))   - (actual[j]*gsl_matrix_get(matrix,l*3,j-1+mode));
				byi += (actual[j]*gsl_matrix_get(matrix,k*3+1,j-1+mode)) - (actual[j]*gsl_matrix_get(matrix,l*3+1,j-1+mode));
				bzi += (actual[j]*gsl_matrix_get(matrix,k*3+2,j-1+mode)) - (actual[j]*gsl_matrix_get(matrix,l*3+2,j-1+mode));			
			}
			newstrc_dist = sqrt(bxi*bxi+byi*byi+bzi*bzi);
			//printf("	Limit:%f	Dist:%f",dist*limit,newstrc_dist);		
			if (dist*limit >= newstrc_dist && dist*(limit-1) <= newstrc_dist) {} else {
				//printf("	Atom:%d %d Limit:[%f]	newstrc Dist:%f Ori:%f\n",k,l,dist*limit,newstrc_dist,dist);	
				return(0);}
			}
		}
		
		for (k=0;k<atom;++k) {
 			for(l=k+1;l<atom;++l) {
 				for(m=l+1;m<atom;++m) {
 					int connect_flag = 0;
	 				int t;
 					for (t=0;t<6;++t) {
	 					if (strc[k].node_c[t] == strc[l].node) {++connect_flag;}
	 					if (strc[l].node_c[t] == strc[m].node) {++connect_flag;}
	 				}
	 				if (connect_flag != 0 && connect_flag != 1) {} else {continue;}
	 				
	 				gsl_vector *u = gsl_vector_alloc(3);
	 				gsl_vector *v = gsl_vector_alloc(3);
	 				
	 				gsl_vector_set(u,0,strc[k].x_cord - strc[l].x_cord);
	 				gsl_vector_set(u,1,strc[k].y_cord - strc[l].y_cord);
	 				gsl_vector_set(u,2,strc[k].z_cord - strc[l].z_cord);
	 				
	 				gsl_vector_set(v,0,strc[l].x_cord - strc[m].x_cord);
	 				gsl_vector_set(v,1,strc[l].y_cord - strc[m].y_cord);
	 				gsl_vector_set(v,2,strc[l].z_cord - strc[m].z_cord);
	 				
	 				double scalar = gsl_vector_get(u,0)*gsl_vector_get(v,0)+gsl_vector_get(u,1)*gsl_vector_get(v,1)+gsl_vector_get(u,2)*gsl_vector_get(v,2);
	 				double lenght_u = gsl_vector_get(u,0)*gsl_vector_get(u,0)+gsl_vector_get(u,1)*gsl_vector_get(u,1)+gsl_vector_get(u,2)*gsl_vector_get(u,2);
	 				double lenght_v = gsl_vector_get(v,0)*gsl_vector_get(v,0)+gsl_vector_get(v,1)*gsl_vector_get(v,1)+gsl_vector_get(v,2)*gsl_vector_get(v,2);
	 				//printf("Scalar:%g 
	 				double angle = acos(scalar/(lenght_u*lenght_v))*180/3.1416;
	 				//printf("Angle between %d %d %d => %f\n",k,l,m,angle);
	 				
	 				struct pdb_atom newstrc[atom];
	 				
	 							
			 		newstrc[k].x_cord = strc[k].x_cord;
			 		newstrc[k].y_cord = strc[k].y_cord;
			 		newstrc[k].z_cord = strc[k].z_cord;
			 		
			 		newstrc[l].x_cord = strc[l].x_cord;
			 		newstrc[l].y_cord = strc[l].y_cord;
			 		newstrc[l].z_cord = strc[l].z_cord;
			 		
			 		newstrc[m].x_cord = strc[m].x_cord;
			 		newstrc[m].y_cord = strc[m].y_cord;
			 		newstrc[m].z_cord = strc[m].z_cord;
	 				
	 				for (j=0;j<nb_mode;++j) {			
		 				newstrc[k].x_cord += (actual[j]*gsl_matrix_get(matrix,strc[k].node*3,  j+mode-1));
		 				newstrc[k].y_cord += (actual[j]*gsl_matrix_get(matrix,strc[k].node*3+1,j+mode-1));
		 				newstrc[k].z_cord += (actual[j]*gsl_matrix_get(matrix,strc[k].node*3+2,j+mode-1));
		 			}
	 				
				
	 				gsl_vector_set(u,0,newstrc[k].x_cord - newstrc[l].x_cord);
	 				gsl_vector_set(u,1,newstrc[k].y_cord - newstrc[l].y_cord);
	 				gsl_vector_set(u,2,newstrc[k].z_cord - newstrc[l].z_cord);
	 				
	 				gsl_vector_set(v,0,newstrc[l].x_cord - newstrc[m].x_cord);
	 				gsl_vector_set(v,1,newstrc[l].y_cord - newstrc[m].y_cord);
	 				gsl_vector_set(v,2,newstrc[l].z_cord - newstrc[m].z_cord);
	 				
	 				double scalar_x = gsl_vector_get(u,0)*gsl_vector_get(v,0)+gsl_vector_get(u,1)*gsl_vector_get(v,1)+gsl_vector_get(u,2)*gsl_vector_get(v,2);
	 				double lenght_ux = gsl_vector_get(u,0)*gsl_vector_get(u,0)+gsl_vector_get(u,1)*gsl_vector_get(u,1)+gsl_vector_get(u,2)*gsl_vector_get(u,2);
	 				double lenght_vx = gsl_vector_get(v,0)*gsl_vector_get(v,0)+gsl_vector_get(v,1)*gsl_vector_get(v,1)+gsl_vector_get(v,2)*gsl_vector_get(v,2);
	 				//printf("Scalar:%g 
	 				double anglex = acos(scalar_x/(lenght_ux*lenght_vx))*180/3.1416;
	 					 				
	 				if (angle*limit >= anglex && angle*(limit-1) <= anglex) {} else {
						printf("Angle between %d %d %d => %f et %f\n",k,l,m,angle,anglex);
					return(0);}
						 				
	 				gsl_vector_free(u);
	 				gsl_vector_free(v);				
	 				
 				}
 			}
 		}
		
		return(1);		
}
 int build_grid(double **grid,double *resolution,int nbr_mode,int atom,gsl_matrix *matrix,struct pdb_atom *strc,float limit,int mode, int now, double *actual) {
	
	int i,j,l=now;
	double passing[nbr_mode];
	for(i=0;i<nbr_mode;++i) {passing[i] = actual[i];}
	
	//printf("L:%d	",l);
	//for(i=0;i<nbr_mode;++i) {printf("%4f\t",actual[i]);}
	
	int ok = test_point(strc, actual, limit, atom,mode, nbr_mode,matrix);
	//if (l > 50) {printf("More then 50\n");return(l);}
		
	for(j=0;j<l;++j) {
		int temp = 0;
		for(i=0;i<nbr_mode;++i) {
			if (actual[i] == grid[j][i]) {++temp;}
			//printf("%f =? %f Temp:%d\n",actual[i],grid[j][i],temp);
		}
		if (temp == nbr_mode) {return(l);}
	}
	
	if (ok == 0) {return(l);}
	if (ok == 1) {
		for(i=0;i<nbr_mode;++i) {
			grid[l][i] = actual[i];
			//printf("%.3f\t",actual[i]);
		}
		//printf("\n");
		++l;
	}
	
	for (i=0;i<nbr_mode;++i) {
		passing[i] += resolution[i];
		l = build_grid(grid,resolution,nbr_mode,atom,matrix,strc,limit,mode,l,passing);
	}
	for(i=0;i<nbr_mode;++i) {passing[i] = actual[i];}
	for (i=0;i<nbr_mode;++i) {
		passing[i] -= resolution[i];
		l = build_grid(grid,resolution,nbr_mode,atom,matrix,strc,limit,mode,l,passing);
	}
	//printf("\nThis IS THE END\n");
	return(l);
}

void print_image_torsion(struct pdb_atom *old,gsl_matrix *evec,int mode,float max,float min,int atom,char out_name[100],int nbnode,struct pdb_atom *strc_node) {
 	FILE *out_file;
 	out_file = fopen(out_name,"w");
 	int l,k;
 	float pi = 3.14159265358979323846264338327950288419716939937510;
 	float amplitude;

 	struct pdb_atom newstrc[atom];
 	struct pdb_atom newnode[nbnode];
 	int align[nbnode];
 	for(k=0;k<nbnode;++k) {
 		align[k] = k;
 	
 	}
 	gsl_vector *vect = gsl_vector_alloc(3);
 	gsl_matrix *rota = gsl_matrix_alloc(3,3);
 	for (l = 0; l<31;++l) {
 		
 		amplitude = sin(l*pi/15)*(max-min)/2+(max+min)/2;
 		
 		printf("Frame:%2d\tAmplitude:%f\n",l+1,amplitude);
 		fprintf(out_file,"Model %d\n",l+1);
		// Nous avons une nouvelle strucuture qui va être eigen bouger
		
		// Copie les coordoné dans newstrc
	 	copy_strc(newstrc,old,atom);
	 	copy_strc(newnode,strc_node,nbnode);
	   // On va appliquer toutes les rotations
		
	   for (k=0;k<nbnode;++k) {	
			
			center_node(atom,newstrc,k);
			center_node(nbnode,newnode,k);
	 		
	 		assign_vector(newstrc,k," CA  ",atom,newstrc,k," C   ",atom,vect);
			rotate_matrix(rota,gsl_matrix_get(evec,k*2,mode-1)*amplitude,vect);
		 	rotate_phi(rota,newstrc,atom,k);
		 	rotate_phi(rota,newnode,nbnode,k);
		 	
		 	assign_vector(newstrc,k," CA  ",atom,newstrc,k," N   ",atom,vect);
			rotate_matrix(rota,gsl_matrix_get(evec,k*2+1,mode-1)*amplitude,vect);
		 	rotate_psy(rota,newstrc,atom,k);
		 	rotate_psy(rota,newnode,nbnode,k);
		 	
		 	
		 	rmsd_yes(newnode,strc_node,nbnode, align,newstrc,atom);
		 	
		// 	write_strc("nod_after_rmsd.pdb",newnode,nbnode);
		 	//break;
	 	}	 	
	 	//printf("RMSD IS:%f\n",);
	write_movie(out_file,newstrc,atom,k); 
 	

 	}
 	fclose(out_file);
 }

