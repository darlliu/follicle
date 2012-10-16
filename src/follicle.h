#include "headers.h"

namespace fol
{
    
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
            follicle ();                             /* constructor */

            /* ====================  ACCESSORS     ======================================= */
            to_path(); // adaptor(s) to get information to pathways

            /* ====================  MUTATORS      ======================================= */
            evolve (); // all the local operations (determine state, birth and decay, noise)

            /* ====================  OPERATORS     ======================================= */

        protected:
            /* ====================  DATA MEMBERS  ======================================= */
            unsigned int num_path; //number of pathways
            unsigned int *num_lig, *num_ang, *num_receptor;// number of ligands, antagonists, receptors
            int *path_type;     // pathway type
            bool major; // whether to output data row major or column major
            unsigned int index; // index of follicle, convertible to (x, y)
                                // note: this shall become the element number (reordered) if meshed

            unsigned int rowl, coll; // row and column length of the grid
            double **I, **D, **ligands, **activators, **receptors, **active_receptors;  
            // synthesis, decay, ligand concentration (final), activator, receptor and receptor bound conc.
            double I_n[2], D_n[2], ligands_n[2], activators_n[2], receptors_n[2], active_receptors_n[2];  
            // the guassian mean and var of noise in each term

            

        private:
            /* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class Follicle  ----- */

}
