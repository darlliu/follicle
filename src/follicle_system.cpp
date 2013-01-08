#include "follicle_system.h"

namespace fol
{

	void fol_sys::allocate_follicles()
	{
		//spacing is adjusted here
		if (!paramset) throw ("Parameters unset!");
		unsigned counter=0;
		follicle *f1=new follicle();
		f1->l1=l1; f1->l2=l2; f1->l3=l3;
		for (unsigned i=0; i<l1; i+=4)
		for (unsigned j=0; j<l2; j+=4)
		{
			follicles.push_back(*f1);
			follicles[counter].posit(i,j,0,counter);
			
			counter++;
		}
		delete f1;
		return;
	}
	void fol_sys::allocate_pathways()
	{
		if (!paramset) throw("Parameters unset!");
		pathway *p1=new pathway(), *p2=new pathway();
		double rxn[7]={1e-3, 1e-3, 5e-4, 5e-4, 1e-7, 1e-7, 1e-2}, 
			rxn2[7]={1e-3, 1e-3, 5e-4, 5e-4, 1e-7, 1e-7, 1e-2};
		p1->rxn_rates.assign(rxn, rxn+7);
		p2->rxn_rates.assign(rxn2,rxn2+7);
		pathways.push_back(*p1);
		pathways.push_back(*p2);
		delete p1, p2;
		num_pathway=2;
		return;
	}
	void fol_sys::read_param()
	{

		paramset=1;
	}

	void fol_sys::makelaplacian_abbd2()
	{
		if (!paramset) throw("Error: no parameters set!");
		//we ask that the size of grid is divisible by 2 4 8 and 10
		if (40%l1||40%l2||40%l3) throw("Error: Grid size incorrect!");

		for (unsigned i = 0; i < l1; ++i){
		for (unsigned j = 0; j < l2; ++j){
		for (unsigned k = 0; k < l3; ++k){
			rows.push_back(getidc(i,j,k));
			cols.push_back(getidc(i,j,k));
			nabla.push_back(6);
			if (i==0)
			{
				rows.push_back(getidc(i,j,k));
				cols.push_back(getidc(i+1,j,k));
				nabla.push_back(-2);
			}
			else if (i==l1-1)
			{
				rows.push_back(getidc(i,j,k));
				cols.push_back(getidc(i-1,j,k));
				nabla.push_back(-2);
			}
			else
			{
				rows.push_back(getidc(i,j,k));
				cols.push_back(getidc(i-1,j,k));
				nabla.push_back(-1);
				rows.push_back(getidc(i,j,k));
				cols.push_back(getidc(i+1,j,k));
				nabla.push_back(-1);
			}
			if (j==0)
			{				
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,l2-1,k));
				nabla.push_back(-1);
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j+1,k));
				nabla.push_back(-1);
			}
			else if (j==l2-1)
			{
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,0,k));
				nabla.push_back(-1);
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j-1,k));
				nabla.push_back(-1);
			}
			else 
			{
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j-1,k));
				nabla.push_back(-1);
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j+1,k));
				nabla.push_back(-1);
			}
			if (k==0)
			{
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j,l3-1));
				nabla.push_back(-1);
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j,k+1));
				nabla.push_back(-1);
			}
			else if (k==l3-1)
			{
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j,0));
				nabla.push_back(-1);
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j,k-1));
				nabla.push_back(-1);
			}
			else 
			{
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j,k-1));
				nabla.push_back(-1);
				cols.push_back(getidc(i,j,k));
				rows.push_back(getidc(i,j,k+1));
				nabla.push_back(-1);
			}
		}}}
		;;//nested for loops
		//note: this is a sparse container for fortran routines
		//and only the lower triagle is represented.
	}

	
}