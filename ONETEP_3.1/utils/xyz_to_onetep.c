/*

!===============================================================!
! This program converts a .xyz file to a ONETEP .dat input file.!
!---------------------------------------------------------------!
! Originally  written by Chris-Kriton Skylaris on 06/10/2001    !
! for the ONES code.                                            !
! Modified for ONETEP and improved by Chris-Kriton Skylaris on  !
! 09/02/2005.                                                   !
!===============================================================!

*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define paragon 1.889726313


/* cks C function prototypes */

int find_atomic_number(char atomo[]);

void test_arguments(int argc, char *argv[]);




main(int argc, char *argv[])
{

  FILE *fp_in,*fp_out;
  

int nat,nuc;

typedef struct {
  char name[3];
  int atomic_number;
  double x;
  double y;
  double z;
  double charge;
} particle;


int count;

 particle *nuclei;



/* cks: test that an input and output file name are provided as arguments and 
if not print a help message */ 
test_arguments( argc, argv);


/* cks: open .xyz file of cartesian atomic coordinates in Angstrom */
fp_in=fopen(argv[1],"r");  

/* cks: read number of atoms "nat"  */
fscanf(fp_in,"%d",&nat); 

/* cks: allocate memory for particle structure pointer  */
nuclei=(particle*)malloc(nat*sizeof(particle));


/* cks: skip the title line of the .xyz file */ 
count=0;
while( getc(fp_in) !='\n' )
count++;
while( getc(fp_in) !='\n' )
count++;





/* cks: open .dat output file for writing */ 
fp_out=fopen(argv[2],"w");


/* cks: output simulation cell, ppd and grid information  */
fprintf(fp_out,"! This file was created with the xyz to ONETEP format converter \n");
fprintf(fp_out,"! Version 1.0, (c) Chris-Kriton Skylaris 2001 and 2005  \n \n \n");




fprintf(fp_out,"! << LATTICE VECTORS (Cartesian coordinates in a0) >>  \n");
fprintf(fp_out,"%%block   lattice_cart            \n");
fprintf(fp_out,"  20.00     0.00     0.00        \n");
fprintf(fp_out,"   0.00    20.00     0.00        \n");
fprintf(fp_out,"   0.00     0.00    20.00        \n");
fprintf(fp_out,"%%endblock  lattice_cart          \n\n\n");


/* print calculation parameters */
fprintf(fp_out,"! << CALCULATION PARAMETERS >>   \n");
fprintf(fp_out,"maxit_ngwf_cg   35    ! Max number of NGWF iterations \n");
fprintf(fp_out,"cutoff_energy  600 eV ! K.E. cut-off \n");
fprintf(fp_out,"kernel_cutoff   20.0  ! Density matrix kernel cut-off (a0) \n");
fprintf(fp_out,"use_SPAM2        T    ! Use sparse matrix algebra  \n\n\n");


   
/* cks: read line by line the symbol and coordinates of each atom and output them
line by line in the format required by ONETEP */ 
fprintf(fp_out,"! << COORDINATES OF ATOMS IN SIMULATION CELL (in a0) >> \n");
fprintf(fp_out,"%%block  geometry    \n");
 
for( nuc=0; nuc< nat ; nuc++){

     fscanf(fp_in,"%s %lf %lf %lf",nuclei[nuc].name
     ,&nuclei[nuc].x,&nuclei[nuc].y,&nuclei[nuc].z  );

     /* cks: Angstrom to Bohr conversion */
     nuclei[nuc].x *= paragon;
     nuclei[nuc].y *= paragon;
     nuclei[nuc].z *= paragon;

     /* cks: find the atomic number of the chemical symbol of the current atom */
     nuclei[nuc].atomic_number =find_atomic_number(nuclei[nuc].name);


     
     /* cks: output line with element name, at. number, pseudo name, coordinates
	, nfunctions, support radius */

     switch (nuclei[nuc].atomic_number) {
       
     case 1:
       
       fprintf( fp_out,"%3s   %3i   %2s_00.recpot    %10.6f  %10.6f  %10.6f  1  6.00 \n"
		,nuclei[nuc].name
		,nuclei[nuc].atomic_number
		,nuclei[nuc].name
		,nuclei[nuc].x
		,nuclei[nuc].y
		,nuclei[nuc].z); 
                break ; 
       
     case 6:
       
       fprintf( fp_out,"%3s   %3i   %2s_01.recpot    %10.6f  %10.6f  %10.6f  4  6.30 \n" 
		,nuclei[nuc].name
		,nuclei[nuc].atomic_number
		,nuclei[nuc].name
		,nuclei[nuc].x
		,nuclei[nuc].y
		,nuclei[nuc].z);  
                break ;
       
     case 8:
       
       fprintf( fp_out,"%3s   %3i   %2s_02.recpot    %10.6f  %10.6f  %10.6f  4  6.30 \n"
		,nuclei[nuc].name
		,nuclei[nuc].atomic_number
		,nuclei[nuc].name
		,nuclei[nuc].x
		,nuclei[nuc].y
		,nuclei[nuc].z); 
                break ;

     case 15:

       fprintf( fp_out,"%3s   %3i   %2s_00.recpot    %10.6f  %10.6f  %10.6f  9  6.30 \n"
                ,nuclei[nuc].name
                ,nuclei[nuc].atomic_number
                ,nuclei[nuc].name
                ,nuclei[nuc].x
                ,nuclei[nuc].y
                ,nuclei[nuc].z);
                break ;

     case 20:
       fprintf( fp_out,"%3s   %3i   %2s_00.recpot    %10.6f  %10.6f  %10.6f  5  6.30 \n"
		,nuclei[nuc].name
		,nuclei[nuc].atomic_number
		,nuclei[nuc].name
		,nuclei[nuc].x
		,nuclei[nuc].y
		,nuclei[nuc].z); 
                break ;

       
     default:
       
       fprintf( fp_out,"%3s   %3i   %2s_00.recpot    %10.6f  %10.6f  %10.6f  4  6.30 \n"
		,nuclei[nuc].name
		,nuclei[nuc].atomic_number
		,nuclei[nuc].name
		,nuclei[nuc].x
		,nuclei[nuc].y
		,nuclei[nuc].z);  
       
     }     
     
}

 
 
fprintf(fp_out,"%%endblock  geometry    \n\n\n");

/* cks: close input .xyz file */
fclose(fp_in);



/* cks: one last loop to output atomic set info */
fprintf(fp_out,"! << ATOMIC SETS >>   \n");
fprintf(fp_out,"%%block atomic_sets       \n");
for( nuc=0; nuc< nat ; nuc++){




  switch (nuclei[nuc].atomic_number) {
       
  case 1:
    
    fprintf( fp_out,"  H_SZ_6.0_recpot00_600eV.fbl \n");   
    break;

  case 6:
        
    fprintf( fp_out,"  C_SZ_6.0_recpot00_600eV.fbl \n"); 
    break;
    
  case 7:

    fprintf( fp_out,"  N_SZ_6.0_recpot01_600eV.fbl \n"); 
    break;

  case 8:
    
    fprintf( fp_out,"  O_SZ_6.0_recpot01_600eV.fbl \n"); 
    break;

  case 15:
        
    fprintf( fp_out,"  P_STO-3Gs_val.fbl \n"); 
    break;

  case 20:

    fprintf( fp_out,"  Ca_SZ_6.0_recpot00_600eV.fbl \n");
    break;
       
  default:
       
  fprintf( fp_out,"  %s_SZ_6.0_fireball.fbl   \n"
	   ,nuclei[nuc].name);
  }



}
fprintf(fp_out,"%%endblock atomic_sets      \n\n\n");





/* cks:  close the output file */
fclose(fp_out);


free(nuclei);

return 0;

}






int find_atomic_number(char atomo[]){


    typedef struct {
        double number;
        char   *symbol; 
    } atom;

    int number=0;

#include "elements.h"


/* Loop until meet element with given name, or until you reach element with at. number 110  */
    do {
        number++;  
        if ( strcmp(atomo, element[number].symbol) == 0 ) break;    
    } while( number < 110  );

    if (number == 110) {
      printf("\n Unable to find the atomic number of element with symbol: %s \n",atomo);
      printf("\n ONETEP input construction aborted \n");
      exit(1);
    }

    return number;
       


}






void test_arguments(int argc, char *argv[]){


#include <ctype.h>

    if ( (argc == 1) || (argc == 2) ){
      printf("\nConvert an xmol xyz format file to a ONETEP input file:\n"
	     "xyz_to_ONETEP file_name.xyz file_name.dat \n"
	     "\nCopyright Chris-Kriton Skylaris, 2001 and 2005.\n" );
      exit(2);
    }

   
    if ( (isalpha(*argv[1])==0) || (isalpha(*argv[2])==0)) {
      printf("\n File names must start with a letter \n");
      exit(4);
    }

}









