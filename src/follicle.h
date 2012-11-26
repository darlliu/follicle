#ifndef GRANDFATHER_H
#define GRANDFATHER_H

#include "headers.h"
#include "pathway.h"

namespace fol
{
    typedef enum 
    {
        unknown = -1,
        r_telo=0,
        c_telo,
        p_ana,
        a_ana,
        cata,
    } state;

    
 
    /*
     * =====================================================================================
     *        Class:  follicle
     *  Description:  a class of follicles with predefined pathways and data members.
     * =====================================================================================
     */
    class follicle
    {
        public:
            /* ====================  LIFECYCLE     ======================================= */
            follicle (unsigned int cycles_in) ;/* constructor */
            follicle ();
            ~follicle();
            void init(unsigned int cycles_in) ;
            /* ====================  ACCESSORS     ======================================= */
            state cur_states()
            {
                return states[cycle];
            };
            
			std::vector<state> all_states(){ return states; };

            //double* get_active_receptors(unsigned int pathway_num) 
                //{return active_receptors[pathway_num];};

            /* ====================  MUTATORS      ======================================= */
            void evolve (); 
            // do the following:
            // 0, make sure new ligand/antagonists info are loaded from receptors and at current cycle.
            // 1, set the new state depending on our thresholds
            // algo: touch multiple times
            // bind -> write
            // evolve, grow
            
            void grow();
            // do the following:
            // 1. depending on the current state, grow once or not
            // 2. set new growth states, inform receptors.

            void set_path (unsigned int pathway_num)
            {
                if (pathway_num>num_path) throw ("Critical error: pathway_num wrong!");
                paths=pathways[pathway_num];  
            };
            void bind(unsigned int pathway_num);
            // for each receptor with binding affinity x
            // do the binding given the correct amount of gradient
            
            void write(unsigned int pathway_num);

            void touch(unsigned int pathway_num);

            void set_next_ligands ( double *in)
            {
                lig=in;
            };
            void set_next_antagonists (double *in)
            {
                ant=in;
            };
			void
			follicle::gen ( unsigned int pathway_num);
            void posit ( const unsigned int l11, const unsigned int l22, const unsigned int l33, unsigned int idx );
            void add_path ( pathway path );
            // set next data arrays
            void generate_noise (double *mean_n, double *var_n,unsigned int pathway_num);
            // given the guassian mean and variance of noise,
            // and the pathway number
            /* ====================  OPERATORS     ======================================= */

        protected:
            /* ====================  DATA MEMBERS  ======================================= */
            unsigned int cycles, l1,l2,l3, t_factor;
                                       // number of cycles, follicle index and global bounds
                                       // note: does not change
            unsigned int cycle, cnt, t;    // current cycle
            unsigned int num_path; //total number of pathways
            //std::vector<unsigned int> num_lig, num_ang, num_receptor;
                               // number of ligands, antagonists, receptors as indexed by pathway number
                               // these are 1 by default
                               // they are unused as of now!
            bool major;         // whether to output data row major or column major
            unsigned int index; // index of follicle, convertible to (d1, d2, d3);
                                // note: this shall become the element number (reordered) if meshed
            std::vector<grid> top, bulge, dp, prc;
            // position vectors. top can always be determined by index
            // since the follicle grid is 75 by 75, top = l1/75*index+l2/75*(index-75), l3';
            // where l1 l2 always divisible by 75 (or sth else) and l3' is a const.
            
            std::vector<state> states; // states if each cycle of this follicle
            //std::vector<double*> ligands, antagonists;  
            std::vector<receptors> lig_r, ant_r;
            // the second receptor vector is optional
            // for example in wnt dkk may be used with receptor as well

            pathway paths;
            std::vector<pathway> pathways;
            // handles to pathways
            double *lig, *ant;
            // collapsed data array -- big F vector
            //
            // the algorithm for follicle should be this:
            // 1. loop through pathways and touch() to generate R vector
            // 1.5 repeat 1 and get next round ligands/antagonists
            // 2. when next round of pathway ligands and antagonists are determined 
            //    bind() to determine receptors bound
            // 3. rinse and repeat
        private:
            /* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class Follicle  ----- */

}
#endif