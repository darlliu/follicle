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

            void cur_ligands (unsigned int pathway_num, double * io);
            void cur_antagonists (unsigned int pathway_num, double * io);
            void cur_receptors (unsigned int pathway_num, double * io );
            void cur_active_receptors (unsigned int pathway_num, double * io);
            

            double * all_ligands (unsigned int pathway_num) {return ligands[pathway_num];};
            double * all_antagonists (unsigned int pathway_num) {return antagonists[pathway_num];};
            double * all_receptors (unsigned int pathway_num){return receptors[pathway_num];};
            double * all_active_receptors (unsigned int pathway_num){return active_receptors[pathway_num];};

            /* ====================  MUTATORS      ======================================= */
            void evolve (); 
            // set the new state depending on our metric
            
            void set_next_ligands (unsigned int pathway_num, double *in);
            void set_next_antagonists (unsigned int pathway_num, double *in);
            void set_next_receptors (unsigned int pathway_num, double *in);
            void set_next_active_receptors (unsigned int pathway_num, double *in);
            
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
                               // number of ligands, antagonists, receptors
            double* path_type;     // pathway types, positive has greater than zero
            bool major;         // whether to output data row major or column major
            unsigned int index; // index of follicle, convertible to (x, y)
                                // note: this shall become the element number (reordered) if meshed

            unsigned int rowl, coll; // row and column length of the grid
            states * state; // states if each cycle of this follicle
            std::vector<double*> ligands, ligands_noise, antagonists, antagonists_noise, receptors, active_receptors;  
            // synthesis, decay, ligand concentration (final), activator, receptor and receptor bound conc.
            // if more than one species per pathway then the species will be collapsed columnwise
            

        private:
            /* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class Follicle  ----- */

}
