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
     *        Class:  Grid
     *  Description:  a simple class of affined grids that allow simple enumeration
     * =====================================================================================
     */
    class grid
    {
        public:
            /* ====================  LIFECYCLE     ======================================= */
            grid (direction_in, offset_in,l11,l22,l33)
            {
                direction=direction_in;
                offset=offset_in;
                d1=0;
                d2=0;
                d3=0;
                l1=l11;
                l2=l22;
                l3=l33;
            };                             /* constructor */

            /* ====================  ACCESSORS     ======================================= */
            unsigned int getx() {return d1};
            unsigned int gety() {return d2};
            unsigned int getz() {return d3};
            void get(unsigned int* io)
            {
                if (size(io)/size(int)<3)
                    throw("Error: io parser size incorrect");
                io[0]=d1;
                io[1]=d2;
                io[2]=d3;
            };
            /* ====================  MUTATORS      ======================================= */
            void set(unsigned int d11, unsigned int d22, unsigned int d33)
            {
                d1=d11;
                d2=d22;
                d3=d33;
            };

            /* ====================  OPERATORS     ======================================= */

            grid& operator++ ()
            {
                if(cur_offset==0)
                {
                    d3++;
                    cur_offset=offset;
                    test1(d3,l3);
                }
                else
                {
                    cur_offset--;
                    if (direction==0)
                    {
                        d1++;
                        test1(d1,l1);
                    }
                    else
                    {
                        d2++;
                        test1(d2,l2);
                    }
                }
                return this;
            };
            grid& operator-- ()
            {
                if(cur_offset==0)
                {
                    d3--;
                    cur_offset=offset;
                    test2(d3,l3);
                }
                else
                {
                    if (direction==0)
                    {
                        d1--;
                        test2(d1,l1);
                    }
                    else
                    {
                        d2--;
                        test2(d2,l2);
                    }
                    cur_offset--;
                }
                return this;
            };
            grid& operator+ (int in)
            {
                return plus(in);
            };
            grid& operator- (int in)
            {
                return minus (in);
            };
        protected:
            /* ====================  METHODS       ======================================= */
            grid& plus (int in)
            {
                if (in<0)
                    return minus (-in);
                else
                {
                    for (int i = 0; i< in; i++)
                        (*this)++;
                    return this;
                }
            };

            grid& minus (int in)
            {
                if (in<0)
                    return plus (-in);
                else
                {
                    for (int i = 0; i< in; i++)
                        (*this)--;
                    return this;
                }
            };

            void test1(d&,l)
            {
                // does not protect against right boundary overflow on non-periodic z
                if (d<l)
                    return;
                else
                    if (l==0)
                        continue;
                    else
                    {
                        d=0;
                    }
            };

            void test2(d&,l)
            {
                if (d>0)
                    return;
                else
                    if (l==0)
                        throw("Critical error with index: negative");
                    else
                    {
                        d=l;
                    }
            };
            /* ====================  DATA MEMBERS  ======================================= */

        private:
            /* ====================  METHODS       ======================================= */
            unsigned int offset, cur_offset, d1, d2, d3, l1, l2, l3;
            bool direction;
            /* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class Grid  ----- */



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
            // 2. 
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

            int* path_type;     // pathway types, positive has greater
            bool major;         // whether to output data row major or column major
            unsigned int index; // index of follicle, convertible to (d1, d2, d3);
                                // note: this shall become the element number (reordered) if meshed
            std::vector<grid> top, bulge, dp, prc; 
            states * state; // states if each cycle of this follicle
            std::vector<double*> ligands, antagonists, active_receptors;  
            // if more than one species per pathway then the species will be collapsed columnwise
            // otherwise the length of each entries will be the number of cycles, the index will be pathway number

        private:
            /* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class Follicle  ----- */

}
