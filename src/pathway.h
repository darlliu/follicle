#include"headers.h"
namespace fol
{
    /*
     * =====================================================================================
     *        Class:  pathway
     *  Description:  
     * =====================================================================================
     */
    class pathway
    {
        public:
            /* ====================  LIFECYCLE     ======================================= */
            pathway ()
            {
                aff_lig=aff_ant=0;
                g_lig=g_ant=0;
                lig=ant=NULL;
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
            double * lig(){return lig};
            double * ant(){return ant};
            double gen_lig(unsigned int in)
            {
                return g_lig;
            };
            double ant_lig(unsigned int in)
            {
                return g_ant;
            };
            double lig_aff() {return aff_lig};
            double ant_aff() {return aff_ant};
            double lig_off() {return off_lig};
            double ant_off() {return off_ant};
            
            bool lig_thr(double in)
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
            /* ====================  MUTATORS      ======================================= */
            
            /* ====================  OPERATORS     ======================================= */

        protected:
            /* ====================  DATA MEMBERS  ======================================= */
            const unsigned int l1, l2, l3;
            double aff_lig, aff_ant, off_lig, off_ant, eli, deg_lig, deg_ant;
            // reaction rates
            double g_lig, g_ant;
            // generation rates
            double rtot0_lig, rtot0_ant;
            // receptor constants
            double thr1, thr2;
            // thresholds -- 0 means no check.
            double *lig, *ant;
            // data members
            
            



    }; /* -----  end of class Pathway  ----- */

}
