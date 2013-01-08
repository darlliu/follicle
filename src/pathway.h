#ifndef PATHWAY_H
#define PATHWAY_H

#include"headers.h"
namespace fol
{
    /*
     * =====================================================================================
     *        Class:  pathway
     *  Description:  a class of pathways to contain first parameters of the system and then
	 *				  also perform ligand pass-through routines to each follicle
     * =====================================================================================
     */
    class pathway
    {
        public:
            /* ====================  LIFECYCLE     ======================================= */
            pathway ()
            {
				rxn_rates.resize(7);
				gen_rates.resize(2);
				receptor_consts.resize(2);
				thresholds.resize(2);
                lig=ant=0;
            };                             /* constructor */
            ~pathway ()
            {
                delete [] lig;
                delete [] ant; 
            };
            void init ( const unsigned int l11, const unsigned int l22, const unsigned int l33)
            {
                l1=l11;
                l2=l22;
                l3=l33;
                lig=new double[l1*l2*l3];
                ant=new double[l1*l2*l3];
                for (unsigned int i = 0; i < l1*l2*l3; i++) 
                {
                    lig[i]=0;
                    ant[i]=0;
                }
            };


            /* ====================  ACCESSORS     ======================================= */
            double gen_lig(unsigned int in)
            {
                return gen_rates[0];
            };
            double gen_ant(unsigned int in)
            {
                return gen_rates[1];
            };
            double lig_aff() {return rxn_rates[0];};
            double ant_aff() {return rxn_rates[1];};
            double lig_off() {return rxn_rates[2];};
            double ant_off() {return rxn_rates[3];};
            
            bool lig_thr(double in)
            {
                if (in>0) {
                    return 1;
                }
                else return 0;
            };
			bool lig_thc(double in)
            {
                if (in>0) {
                    return 1;
                }
                else return 0;
            };
            bool ant_thr(double in)
            {
                if (in>0) {
                    return 1;
                }
                else return 0;
            };
			double lig_deg()
			{
				return rxn_rates[4];
			};
			double ant_deg()
			{
				return rxn_rates[5];
			};



            /* ====================  MUTATORS      ======================================= */
            void anil()
				//antagonist anniliation part
			{
				double temp;
				for (unsigned i = 0; i < l1*l2*l3; i++)
				{
					temp=lig[i]*ant[i]*rxn_rates[6];
					ant[i]-=temp;
					lig[i]-=temp;
				}
			};


            /* ====================  OPERATORS     ======================================= */
			double *lig, *ant;
            // data members
            /* ====================  DATA MEMBERS  ======================================= */
            static unsigned int l1, l2, l3;
			// grid size, note: only these are statics
            std::vector<double> rxn_rates,gen_rates,receptor_consts,thresholds;
			//double aff_lig, aff_ant, off_lig, off_ant,  deg_lig, deg_ant, elimination rate;
            // reaction rates
            //double g_lig, g_ant;
            // generation rates
            //double rtot0_lig, rtot0_ant;
            // receptor constants
            //double thr1, thr2;
            // thresholds -- 0 means no check.

            
            



    }; /* -----  end of class Pathway  ----- */

}
#endif