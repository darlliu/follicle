#include "follicle.h"

namespace fol
{
	    unsigned grid::offset=2, grid::l1=0, grid::l2=0, grid::l3=0;
	
		unsigned int follicle::cycles=0,
        follicle::t_factor=20,
		follicle::l1=15,follicle::l2=15,follicle::l3=15,
        follicle::num_path = 0;

        //index=0;
        bool follicle::major=1; // default to column major
	

    follicle::follicle(unsigned int cycles_in)
    {
		cnt=0;
		cycle=0;
		t=0;
        init(cycles_in);
    }

	follicle::~follicle()
	{
	}


    follicle::follicle ()
    {
		t=0;
		cnt=0;
		cycle=0;
    }

    void
        follicle::init (unsigned int cycles_in)
    {
        cycles=cycles_in;
        states.resize(cycles);
        for (unsigned i = 0; i<cycles; i++)
        {
            states[i]=unknown;
        }// initialize the states
        /*
		top.resize(cycles);
        dp.resize(cycles);
        bulge.resize(cycles);
        prc.resize(cycles);
		*/
        return;
    }

    void
        follicle::posit ( const unsigned i, const unsigned j, const unsigned k, unsigned int idx )
        {
            index=idx;
            grid* g1= new grid();
            g1->init(1,2,l1,l2,l3);
            g1->set(i,j,k);
            //for test only
            top.push_back(*g1);
            bulge.push_back((*g1)+5);
            bulge.push_back((*g1)+6);
            dp.push_back(*g1+9);
            dp.push_back(*g1+10);
			prc.push_back((*g1)+7);
			prc.push_back((*g1)+8);
			delete g1;
			//std::cout<<prc[0].gety()<<"\t"<<top[0].gety()<<std::endl;

            return ;
        }		/* -----  end of function follicle::posit  ----- */
    void
        follicle::add_path ( pathway *path )
        {
            num_path++;
            pathways.push_back(path);
            return ;
        }		/* -----  end of function follicle::add_path  ----- */
    void
        follicle::grow (  )
        {
            switch (states[cycle])
            {
                case unknown:
                    throw("Currrent cycle not initialized, try evolve?");
                    return;
                case r_telo:
                    return;
                case c_telo:
                    return;
                case p_ana:
                    {
                        for (unsigned int i = 0; i < dp.size(); i++) 
                        {
                            ++dp[i];
                        }
                        return;
                    }
                case a_ana:
                    return;
                case cata:
                    {
                        for (unsigned int i = 0; i < dp.size(); i++) 
                        {
                            --dp[i];
                        }
                        return;
                    }
                default:
                    return;
            }
            // what more to do:
            // inform receptors if they are on the growing parts
            // now there's no need to do so
            return ;
        }		/* -----  end of function follicle::grow  ----- */


    void
        follicle::evolve (  )
        {
            cycle++;
            switch (states[cycle])
            {
                case unknown:
                    states[cycle]=r_telo;
                    break;
                case r_telo:
                    {
                        for (unsigned int i = 0; i < num_path; i++) 
                        {
                            if (!pathways[i]->lig_thr(lig_r[i][cycle]))
                            {
                                states[cycle]=r_telo;
                                break;
                            }
                        }
                        states[cycle]=c_telo;
                        break;
                    }
                case c_telo:
                    {
                        for (unsigned int i = 0; i < num_path; i++) 
                        {
                            if (!pathways[i]->lig_thc(lig_r[i][cycle]))
                            {
                                states[cycle]=c_telo;
                                break;
                            }
                        }
                        states[cycle]=p_ana;
                    }
                case p_ana:
                    {
                        if (dp[0].getz()-top[0].getz()>40)
                            states[cycle]=a_ana;
                        else
                            states[cycle]=p_ana;
                        break;
                    }
                case a_ana:
                    {
                        if (cnt>40)
                        {
                            cnt=0;
                            states[cycle]=cata;
                        }
                        else
                            states[cycle]=a_ana;
                    }
                case cata:
                    {
                        //here we assume originally dp is 10 z distance from
                        //the top --> will reset once done one cycle 
                        if (dp[0].getz()-top[0].getz()<=10) {
                            states[cycle]=r_telo;
                        }
                        else
                            states[cycle]=cata;
                    }
                default:
                    break;
            }
            if(t>t_factor)
            {
                cnt++;
                grow();
                t=0;
            }
            return ;
        }		/* -----  end of function follicle::evolve  ----- */

    void
        follicle::touch ( unsigned int pathway_num, double* F_lig, double * F_ant)
        {
            set_path(pathway_num);
            double templig=0, tempant=0;
            for (unsigned int i = 0; i < bulge.size(); i++) 
            {

                templig=(lig_r[pathway_num].touch(i))*\
                        (paths->lig_aff()*paths->lig[bulge[i].getidr()])+\
                        lig_r[pathway_num].data_now()[i]*paths->lig_off();
                F_lig[bulge[i].getidr()]=-templig;
                tempant=(ant_r[pathway_num].touch(i))*\
                        (paths->ant_aff()*paths->ant[bulge[i].getidr()])+\
                        ant_r[pathway_num].data_now()[i]*paths->ant_off();
                F_ant[bulge[i].getidr()]=-tempant;
           }

           // get bulge positions and convert them to global index
           // then get receptor affinity and receptor conc.
            return ;
        }		/* -----  end of function follicle::bind  ----- */
    void
        follicle::gen ( unsigned int pathway_num)
        {
            set_path(pathway_num);
            double templig=0, tempant=0;
            for (unsigned int i = 0; i < bulge.size(); i++) 
            {

                templig=paths->gen_lig(1);
                lig[bulge[i].getidr()]=templig;
                tempant=paths->gen_ant(1);
                ant[bulge[i].getidr()]=tempant;
           }
            for (unsigned int i = 0; i < top.size(); i++) 
            {

                templig=paths->gen_lig(0);
                lig[top[i].getidr()]=templig;
                tempant=paths->gen_ant(0);
                ant[top[i].getidr()]=tempant;
            }
            if (states[cycle]==p_ana || states[cycle]==a_ana)
				for (unsigned int i = 0; i < prc.size(); i++) 
            {

                templig=paths->gen_lig(2);
                lig[prc[i].getidr()]=templig;
                tempant=paths->gen_ant(2);
                ant[prc[i].getidr()]=tempant;
            }
            for (unsigned int i = 0; i < dp.size(); i++) 
            {

                templig=paths->gen_lig(3);
                lig[dp[i].getidr()]=templig;
                tempant=paths->gen_ant(3);
                ant[dp[i].getidr()]=tempant;
            }
            return ;
        }		/* -----  end of function follicle::bind  ----- */

    void
        follicle::bind ( unsigned int pathway_num )
        {
            set_path(pathway_num);
            double templig=0, tempant=0;
            ++lig_r[pathway_num];
            ++ant_r[pathway_num];
            for (unsigned int i = 0; i < bulge.size(); i++) 
            {
               
                templig=(lig_r[pathway_num].touch(i))*\
                        (paths->lig_aff()*paths->lig[bulge[i].getidr()])+\
                        lig_r[pathway_num].data_now()[i]*paths->lig_off();
                lig[bulge[i].getidr()]=-templig;
                tempant=(ant_r[pathway_num].touch(i))*\
                        (paths->ant_aff()*paths->ant[bulge[i].getidr()])+\
                        ant_r[pathway_num].data_now()[i]*paths->ant_off();
                ant[bulge[i].getidr()]=-tempant;

                lig_r[pathway_num].write(i, templig,paths->lig_deg());
                ant_r[pathway_num].write(i, tempant,paths->ant_deg());
            }

            // get bulge positions and convert them to global index
            // then get receptor affinity and receptor conc.
            return ;
        }		/* -----  end of function follicle::bind  ----- */
}
