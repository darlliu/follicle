#include"classes.h"

int main()
{
	try
	{
		fol::follicle *f1=new fol::follicle();
		(*f1).init(20);
		f1->posit(1,1,1,1);
		fol::pathway *p1=new fol::pathway();
		fol::fol_sys *fsys=new fol::fol_sys();
		fsys->allocate_follicles();
		fsys->allocate_pathways();
		fsys->makelaplacian_abbd2();
		//size_t sz=sizeof(*fsys);
		f1->grow();
		delete f1, p1, fsys;
	}
	catch (char* err)
	{
		std::cout<<err;
		return 1;
	}
	return 0;
}