/*constructs a polymer and creates a LAMMPS datafile for said polymer
 * this polymer is to be let loose in a constrained sphere so making this
 *  starting configuration small is probably a good idea. How small? not quite
 * clear yet, how much
 * may they overlap? Also not quite clear yet. The soft potentials permit these
 * initial configurations to be disgustingly overlapped. We'll probably setup
 * something that passes for a self avoiding walk that doesn't overlap or tie
 * any funny knots.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define in(X, Y, Z, W) X * W * W + Y * W + Z
#define Pi 3.14159265359

struct monomer {
        int id;
        double x;
        double y;
        double z;
};


int i=0,j=0,k=0, l=0, dj=1,  dk=1;

/*arguments:
 * Filename
 * number of monomers on polymer
 * number of monomers on tether
 */
int main(int argc, char *argv[]) {
        int N = atoi(argv[1]);  //total number of monomers for chain
        int T = atoi(argv[2]); //total number of monomers for tether
        double angle = (90-atoi(argv[3]))*Pi/180.0; //angle of original branches in degrees
        int NT=N+T;
        double width = fmax(1.0*N, 1.0*(T+1)); //xy dimensions
        double zlo = -width/2-N*cos(angle);
        //can't mess around with width on any uncessary spots since it will cause problems  
        //e.g. remember the average extension likes to know how far stretched you are from the top?
        double lb=0.95; //bondlength

        double lb_x = lb*cos(angle);
        double lb_z = lb*sin(angle);
        //spacing so that diagonal of cube fills 90% of sphere's length along that line.
        //printf("%d monomers => cube width %d\n", N, n);
	double zpoly =-T*lb+width/2 -1.0;//z coordinate of polymer
	
	//1.0IS A MAGIC NUMBER, EXTRA SHIFT BECAUSE OF WALL SHENANEGANS 
	
	double midshift = (width-(N+1)*lb_x)/2.0;
	if(N%2==0){
	midshift+=lb_x/2.0;
	}
	 //shifts things to the middle
	 
	//printf("AAGh %f %f\n", midshift,);
	if (N%2==0){ midshift-=0.5*lb;
	}
        struct monomer *listomonomers=(struct monomer*)calloc(N+T+1,sizeof(struct monomer));
        int MM=N/2+1; //ID of tethering point (middle Monomer)
        //first entry is empty cuz indexing from 1 is more convenient
        for(i=1; i<=N; i++) {
        	
                listomonomers[i].id = i;
                listomonomers[i].x=lb_x*i-width/2.0 + midshift;
                listomonomers[i].y=0.0;
                listomonomers[i].z=zpoly+-lb_z*abs(i-MM);
        }

        
        //double xtether=lb*floor(N/2)-width/2.0;
        int id = 0;
        //using this x coordinate (fixed at 0) does mean that the starting position of the attachment
        //might not be as orthogonal as we'd like since it's shifted one to the side, just have a look at a N13 T7 one
        for(i=1; i<=T; i++)
        {
          id=i+N; //unshifted ID... why do we not start at 0 again?
          listomonomers[id].id=id;
          listomonomers[id].x=0.0;
          listomonomers[id].y=0.0;
          listomonomers[id].z=zpoly+lb*i;
        }
	
	
	
        //write to tha file
        FILE *configuration = fopen(argv[4],"w");

        fprintf(configuration, "#Arrangement of polymer tethered to a walln\n");
        fprintf(configuration, "%8d atoms\n",NT );
        fprintf(configuration, "%8d bonds\n",NT-1 );
        fprintf(configuration, "%8d angles\n\n",NT-1 );
        fprintf(configuration, "%8d atom types\n",2 );
        fprintf(configuration, "%8d bond types\n",1 );
        fprintf(configuration, "%8d angle types\n\n",2 );
        fprintf(configuration, "%8lf %8lf xlo xhi\n",-width/2,width/2 );
        fprintf(configuration, "%8lf %8lf ylo yhi\n",-width/2,width/2 ); //the 
        fprintf(configuration, "%8lf %8lf zlo zhi\n\n",zlo,width/2);
        fprintf(configuration, "Masses\n\n 1 1\n  2 1\n\n");
        fprintf(configuration, "Atoms\n\n");
        //positions of the first N long chain type particles 
        for ( i = 1; i <= N; i++) {
                fprintf(configuration,"%8d %8d %8d %8lf %8lf %8lf\n",listomonomers[i].id,1,
                        1,listomonomers[i].x,listomonomers[i].y,listomonomers[i].z );
        }
        //positions of the tethering particles (LESS STIFF AND SEPARATE FOR THE PURPOSES OF EXTRACTING DATA)
        for ( i = N+1; i <= NT; i++) {
       	     fprintf(configuration,"%8d %8d %8d %8lf %8lf %8lf\n",listomonomers[i].id, 1,
                        2,listomonomers[i].x,listomonomers[i].y,listomonomers[i].z );
        }
        //bonds
        fprintf(configuration, "\n Bonds\n\n");
        for ( i = 1; i <= N-1; i++) {
                fprintf(configuration,"%8d %8d %8d %8d\n",i, 1,i,i+1);
        }
        //tether point bond
        fprintf(configuration,"%8d %8d %8d %8d\n",N, 1, MM, N+1);
        //bonds relating to the tether
        for ( i = 1; i <= T-1; i++) {
                fprintf(configuration,"%8d %8d %8d %8d\n",N+i, 1, N+i, N+i+1);
        }

        fprintf(configuration, "\n Angles\n\n");
        for ( i = 1; i <= N-2; i++) {
                fprintf(configuration,"%8d %8d %8d %8d %8d\n",i, 1,i,i+1,i+2);
        }
        //tether point angles
        fprintf(configuration,"%8d %8d %8d %8d %8d\n",N-1, 2, MM-1, MM, N+1); //bottom L shape angle
        fprintf(configuration,"%8d %8d %8d %8d %8d\n",N, 2, MM+1, MM, N+1); //top L shape angle
        fprintf(configuration,"%8d %8d %8d %8d %8d\n",N+1,   2, MM, N+1, N+2);
        //Angle connecting two of those on the branch

        //tether's angles with no anchor point involvement
        for ( i = 1; i <= T-2; i++) {
                fprintf(configuration,"%8d %8d %8d %8d %8d\n",N+i+1, 2, N+i, N+i+1, N+i+2);
        }

        free(listomonomers);
        return 0;
}
