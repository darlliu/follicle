#ifndef FOLLICLE_SYSTEM_H
#define FOLLICLE_SYSTEM_H
#include "headers.h"
#include "follicle.h"
#include "pathway.h"

namespace fol
{
	class fol_sys
	{
	public:
		fol_sys(): 
			paramset(1)
		{
			fname="parameters.txt";
			ant=NULL;
			lig=NULL;
			l1=l2=l3=40;
			//default dimension of grid

		};
		void setfile(std::string fnamein){fname=fnamein;};
		void read_param();
		//routines for reading files

		/*------preparation routines -----*/

		void allocate_follicles();
		void allocate_pathways();
		void initialize();
		//apply the IC

		/*------main control flow routines -----*/
		void run();
		//run once to solve for next round
		void runcycles(const unsigned cycles)
		{
			unsigned cycle=0;
			double err=100;
			do
			{
				run();
				cycle++;
			} while (cycle<cycles);
		};


		/*------regular solving routines -----*/
		void makelaplacian_abbd2();
		//make laplacian according to l1 l2 l3
		//and individual follicle placements
		void diffuse();
		//calculate the diffused V vector V=Del*U+F(U)
		//note here we do not scale by dt yet.
		void bigf();
		//aggregate parts of F vectors from elim and from follicles
		void unext();
		//get unext=v*dt+u or some equivalent


		/*------regular output routines -----*/
		void writeout();
		std::stringstream& writestates();
		//sstream out from follicles 
		std::stringstream& writeconc();
		std::stringstream& writepos();
		//sstream out from receptors
		
		/*------other routines -----*/
        unsigned int getidc(unsigned d1, unsigned d2, unsigned d3){return d1+l1*d2+l1*l2*d3;};
        unsigned int getidr(unsigned d1, unsigned d2, unsigned d3){return d3+l3*d2+l3*l2*d1;};
		

	protected:
		std::string fname;
		std::vector<follicle> follicles;
		std::vector<follicle>::iterator fol;
		std::vector<pathway> pathways;
		std::vector<follicle>::iterator path;
		// holders for follicles and pathways

		double *ant, *lig, *ant2, *lig2, *ant3, *lig3, *lig4, *ant4;
		// this is actually the F vectors, actual ligand and 
		// antagonists are stored in individual pathways
		// we expect each to hold about 100m of data at any point, which means
		// for 2 pathways it will be 800m at least. AB4BD4 may be difficult
		unsigned l1, l2, l3, num_pathway, offset_for_growth;
		// basic parameter
		double dt, dx, coefx, coefy, coefz;
		// time factor
		bool paramset;
		std::vector<unsigned> cols, rows;
		std::vector<double> nabla;

		// the flow goes like this:
		/* 1, choose a pathway
		   1.5, diffuse into F, eliminates into F
		   2, for each follicle, touch and gen multiple times to F
		   3, repeat for >2nd order methods
		   4, get next lig and ant
		   5, bind and gen
		   6, evolve and grow
		   7, write.
		   8, repeat for other pathways
		   */



	};

};

#endif