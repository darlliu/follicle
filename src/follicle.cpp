#include "follicle.h"

namespace fol
{
    follicle::follicle(unsigned int cycles_in)
    {
        cycles=cycles_in;
        cycle=0;
        num_path = 0;
        index=rowl=coll=0;
        major=1; // default to column major
        for (int i = 0; i<2; i++)
        {
            I_n[i]=D_n[i]=ligands_n[i]=antagonists_n[i]\
                 = receptors_n[i]=active_receptors_n[i]=0;
        }// initialize to zero noise;
        states=new int [cycles];
        for (int i = 0; i<cycles; i++)
        {
            states[i]=UNKNOWN;
        }// initialize the states
    }

    void follicle::cur_ligands (unsigned int pathway_num, double * io)
    {
        if ((sizeof (io)/sizeof (double))!=num_lig[pathway_num] \
                || pathway_num>num_path)
            throw ("access error: pathway mismatch")
        catch (char* info)
        {
            std::cout<<info<<std::endl;
            return;
        }


    }
}
