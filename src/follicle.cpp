#include "follicle.h"

namespace fol
{
    follicle::follicle(unsigned int cycles_in)
    {
        cycles=cycles_in;
        cycle=0;
        num_path = 0;
        index=0;
        major=1; // default to column major
        /*for (int i = 0; i<2; i++)
        {
            I_n[i]=D_n[i]=ligands_n[i]=antagonists_n[i]\
                 = receptors_n[i]=active_receptors_n[i]=0;
        }// initialize to zero noise;
        */
        states=new int [cycles];
        for (int i = 0; i<cycles; i++)
        {
            states[i]=unknown;
        }// initialize the states
        top.resize(cycles);
        dp.resize(cycles);
        bulge.resize(cycles);
        prc.resize(cycles);

    }
/*  
    void follicle::cur_ligands (unsigned int pathway_num, double * io)
    {
        if ((sizeof (io)/sizeof (double))!=num_lig[pathway_num] \
                || pathway_num>num_path)
            throw ("access error: pathway size incorrect");
        if (ligands[pathway_num]==NULL)
            throw ("local storage array size incorrect");
        catch (char* info)
        {
            std::cout<<info<<std::endl;
            return;
        }
        // the output format is [c1_ligand1, c1_ligand2 ... c2_ligand1...]

        for (int i = 0; i<num_lig[pathway_num]; i++)
            io[i]=ligands[pathway_num][i+num_lig[pathway_num]*cycle];
    }

    void follicle::cur_antagonists (unsigned int pathway_num, double * io)
    {
        if ((sizeof (io)/sizeof (double))!=num_antagonists[pathway_num] \
                || pathway_num>num_path)
            throw ("access error: pathway size incorrect");
        try
        {
            if (antagonists[pathway_num]==NULL)
            throw ("local storage array size incorrect");
        }
        catch (char* info)
        {
            std::cout<<info<<std::endl;
            return;
        }
        // the output format is [c1_ligand1, c1_ligand2 ... c2_ligand1...]

        for (int i = 0; i<num_antagonists[pathway_num]; i++)
            io[i]=antagonists[pathway_num][i+num_antagonists[pathway_num]*cycle];
    }


    void follicle::cur_active_receptors (unsigned int pathway_num, double * io)
    {
        if ((sizeof (io)/sizeof (double))!=num_receptor[pathway_num] \
                || pathway_num>num_path)
            throw ("access error: pathway size incorrect");
        if (num_receptor[pathway_num]==NULL)
            throw ("local storage array size incorrect");
        catch (char* info)
        {
            std::cout<<info<<std::endl;
            return;
        }
        // the output format is [c1_ligand1, c1_ligand2 ... c2_ligand1...]

        for (int i = 0; i<num_receptor[pathway_num]; i++)
            io[i]=active_receptors[pathway_num][i+num_receptor[pathway_num]*cycle];
    }
*/
    void follicle::set_next_ligands (unsigned int pathway_num, double * io)
    {
        if ((sizeof (io)/sizeof (double))!=num_lig[pathway_num] \
                || pathway_num>num_path)
            throw ("access error: pathway size incorrect");
        if (ligands[pathway_num]==NULL)
            throw ("local storage array size incorrect");
        if (ligands[pathway_num][cycle]<0)
            throw ("Wrong cycle number");

        catch (char* info)
        {
            std::cout<<info<<std::endl;
            return;
        }
        // the output format is [c1_ligand1, c1_ligand2 ... c2_ligand1...]

        for (int i = 0; i<num_lig[pathway_num]; i++)
            ligands[pathway_num][i+num_lig[pathway_num]*cycle]=io[i];
    }

    void follicle::set_next_ligands (unsigned int pathway_num, double * io)
    {
        if ((sizeof (io)/sizeof (double))!=num_ang[pathway_num] \
                || pathway_num>num_path)
            throw ("access error: pathway size incorrect");
        if (antagonists[pathway_num]==NULL)
            throw ("local storage array size incorrect");
        if (antagonists[pathway_num][cycle]<0)
            throw ("Wrong cycle number");

        catch (char* info)
        {
            std::cout<<info<<std::endl;
            return;
        }
        // the output format is [c1_ligand1, c1_ligand2 ... c2_ligand1...]

        for (int i = 0; i<num_ang[pathway_num]; i++)
            antagonists[pathway_num][i+num_ang[pathway_num]*cycle]=io[i];
    }
}
