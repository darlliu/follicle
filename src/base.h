#include "headers.h"

namespace fol
{
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
            gird()
            {
                d1=0;
                d2=0;
                d3=0;
            };
            void init (direction_in, offset_in,l11,l22,l33)
            {
                direction=direction_in;
                offset=offset_in;
                l1=l11;
                l2=l22;
                l3=l33;
            };                             /* constructor */

            /* ====================  ACCESSORS     ======================================= */
            unsigned int getx() {return d1};
            unsigned int gety() {return d2};
            unsigned int getz() {return d3};
            unsigned int getidr(){return d3+l3*d2+l3*l2*d1};
            unsigned int getidc(){return d1+l1*d2+l1*l2*d3};
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

        protected:
            /* ====================  METHODS       ======================================= */
            unsigned int offset, cur_offset, d1, d2, d3, l1, l2, l3;
            bool direction;
            /* ====================  DATA MEMBERS  ======================================= */

    }; /* -----  end of class Grid  ----- */





    /*
     * =====================================================================================
     *        Class:  receptors
     *  Description:  a simple class of receptors with basic data io functionality 
     * =====================================================================================
     */
    class receptors
    {
        public:
            /* ====================  LIFECYCLE     ======================================= */
            receptors ();                             /* constructor */
            void init(unsigned int cycles, double rtot0_in)
            {
                data.resize(cycles);
                rtot0=rtot0_in;
                rtot=rtot0;
            };
            ~receptors();
            {
                for (it=data.begin(); it<data.end(); it++)
                    delete [] (*it);
            };
            /* ====================  ACCESSORS     ======================================= */
            double rtot() {return rtot;};
            double* data(int cycle)
            {
                return data[cycle];
            };
            double* data_now()
            {
                return data[now];
            };
            unsigned int now(){return now;};
            double touch(unsigned int i)
            {
                if (i>size(data_now)/size(double)) throw "Critical error: data vector overflow";
                return rtot0-data_now()[i]; 
            };
            /* ====================  MUTATORS      ======================================= */
            void newsrc(grid in)
            {
                src.push_back(in);
                //done! only need to inform not to copy src.
            };
            void write(unsigned int i,double d , double kdeg )
            {
                data_now()[i]=d-kdeg*data[now-1][i];
            };

            /* ====================  OPERATORS     ======================================= */
            int operator++()
            {
                now++;
                if (now>cycles) throw("Overflowing cycles!");
                if (data[now]!=NULL) throw ("Uninitialized!")
                data[now]=new double [src.size()]; 
                for (int i=0; i<src.size(); i++) data[now][i]=data[now-1][i];
                // Here each entry of src is a set of  
                // sources (array of grid elements) which are copied by vector methods
                // data is also copied accordingly. only need to del data for each
                return now;
            };

            double operator[] (int cycle)
            {
                double sum=0;
                for (unsigned int j = 0; j < src.size(); j++) 
                {
                    sum+=data_now()[j];
                }
                return sum;
            };

        protected:
            /* ====================  METHODS       ======================================= */

            /* ====================  DATA MEMBERS  ======================================= */
            std::vector<double*>::iterator it;
            unsigned int cycles, now;
            double rtot0, rtot;    
            std::vector<grid> src;
            std::vector<double*> data;


    }; /* -----  end of class Receptors  ----- */


}
