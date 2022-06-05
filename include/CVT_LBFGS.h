#pragma once
#include "Non_Linear_Optimization.h"
#include "cvt.h"
using namespace std;
class CVT_LBFGS : public Non_Linear_Optimization
{
protected:
	CVT& m_cvt;
	double m_epsilon;
public:
    CVT_LBFGS(CVT& cvt, double epsilon = 1.0e-16, bool fDisplay = true)
		: m_cvt(cvt), m_epsilon(epsilon)
    {
		m_printprogress = fDisplay;
		m_variables.resize(m_cvt.m_new_site.size() * 2); // positions
		for (int i = 0; i < m_cvt.m_new_site.size(); ++i)
		{
			m_variables[i] = m_cvt.m_new_site[i].x();
			m_variables[i + m_cvt.m_new_site.size()] = m_cvt.m_new_site[i].y();
		}
		m_parameters.epsilon = 1e-16;
    }
	double GetTerminationAccuracy() const { return m_epsilon; }
	virtual int  run()
	{
		m_energysequence.clear();
		m_normsequence.clear();
		lbfgsfloatval_t fx = 0;
		int numOfIterations = 0;
		int ret = lbfgs(m_variables.size(), &m_variables[0], &fx, _evaluate, _progress, this, &m_parameters, &numOfIterations);
		return numOfIterations;
	}
protected:
   virtual lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    { 
	   m_cvt.m_old_site = m_cvt.m_new_site;
	   int size = m_cvt.m_new_site.size();
	   m_cvt.m_new_site.clear();
	    for (int i = 0; i <size; ++i)
	   {
			m_cvt.m_new_site.push_back(Point_2(x[i], x[i + size]));
		   m_variables[i] = m_cvt.m_new_site.at(i).hx();
		   m_variables[i + size] = m_cvt.m_new_site.at(i).hy();
		
	   }
		m_cvt.update_partition();
		pair<double, vector<double>> energy_gradient = m_cvt.ComputeEnergyAndGradient_CVT();
	
		double energy = energy_gradient.first;
	
		memcpy(g, &energy_gradient.second[0], sizeof(lbfgsfloatval_t) * n);
		return energy;
    }
 	virtual int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
		m_energysequence.push_back(fx);
		m_normsequence.push_back(Norm_Inf(g, n));
		if (m_printprogress)
		{
			//cerr << " Iteration #" << k << ": ";
			//cerr << " Energy = " << setprecision(10) << m_energysequence.back() << ", ";
			//cerr << " Error = " << setprecision(10) << m_normsequence.back() << ", ";
			//cerr << " step = " << step << "\n";
			//cerr << "-----------------------------------------\n";
		}

		if (m_normsequence.back() < m_epsilon)
			return 10000;
		return 0;
	}
};