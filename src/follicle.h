#include "headers.h"


namespace fol
{
    typedef enum 
    {
        unknown = -1;
        r_telo=0;
        c_telo;
        p_ana;
        a_ana;
    } states;

    
 
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
                
            ~follicle();
            
            /* ====================  ACCESSORS     ======================================= */
            states cur_state()
            {
                return state[cycle];
            };
            
            states * all_state(){ return state; };

            double* get_active_receptors(unsigned int pathway_num) 
                {return active_receptors[pathway_num];};

            /* ====================  MUTATORS      ======================================= */
            void evolve (); 
            // do the following:
            // 0, make sure new ligand/antagonists info are loaded.
            // 1, set the new state depending on our thresholds
            // 2, prepare the generation vector G.
            // 3, calculate the recepor-binding event, receptor decay event, get the R vector
            // 4, combine into F vector and prepare for output
            
            void grow();
            // do the following:
            // 1. depending on the current state, grow once or not
            // 2. set new growth states, inform receptors.
            // 3. 


            void bind();
            // for each receptor with binding affinity x
            // do the binding given the correct amount of gradient

            void set_next_ligands (unsigned int pathway_num, double *in);
            void set_next_antagonists (unsigned int pathway_num, double *in);
            
            void generate_noise (double *mean_n, double *var_n, pathway_num);
            // given the guassian mean and variance of noise,
            // and the pathway number
            /* ====================  OPERATORS     ======================================= */

        protected:
            /* ====================  DATA MEMBERS  ======================================= */
            const unsigned int cycles; // number of cycles
                                       // note: does not change
            unsigned int cycle;    // current cycle
            unsigned int num_path; //total number of pathways
            std::vector<unsigned int> num_lig, num_ang, num_receptor;
                               // number of ligands, antagonists, receptors as indexed by pathway number
                               // these are 1 by default
                               // they are unused as of now!

            int* path_type;     // pathway types, positive has greater
            bool major;         // whether to output data row major or column major
            unsigned int index; // index of follicle, convertible to (d1, d2, d3);
                                // note: this shall become the element number (reordered) if meshed
            std::vector<std::vector<grid>> top, bulge, dp, prc; 
            states * state; // states if each cycle of this follicle
            std::vector<double*> ligands, antagonists;  
            receptors lig_r, ant_r;
            // if more than one species per pathway then the species will be collapsed columnwise
            // otherwise the length of each entries will be the number of cycles, the index will be pathway number

        private:
            /* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class Follicle  ----- */

}
