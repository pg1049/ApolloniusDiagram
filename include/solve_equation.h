#include "non_linear_optimization.h"

class Solve_Equation : public Non_Linear_Optimization
{

public:
	vector<vector<double>>  m_point;
	int dimension;
	void initial(int dimension1, const vector<double>& initialVals, const vector<vector<double>>&  m_points)
	{
		dimension = dimension1;
		m_variables.resize(dimension1);
		m_point = m_points;
		for (int i = 0; i < dimension1; ++i)
		{
			m_variables[i] = initialVals[i];
		}
		m_parameters.epsilon = 1e-6;
	}
	virtual lbfgsfloatval_t evaluate(
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step
	)
	{
		double fx = 0;
		double size = m_point.size();
		vector<double> F;
		vector<vector<double>>G;
		for (int i = 0; i < size; ++i)
		{
			double sum = 0;
			for (int j = 0; j < dimension; ++j)
			{
				sum += (x[j] - m_point.at(i).at(j))* (x[j] - m_point.at(i).at(j));
			}
			double f = sqrt(sum) - m_point.at(i).at(dimension);
			F.push_back(f);
		}

		double Fx = 0;
		for (int i = 0; i < F.size()-1; ++i)
		{
			for (int j = i + 1; j < F.size(); j++)
			{
				Fx += (F[i] - F[j]) * (F[i] - F[j]);
			}
		}
		fx += Fx;
		for (int i = 0; i < dimension; ++i)
		{
			vector<double> g;
			for (int j = 0; j <size; ++j)
			{
				if (F.at(j) + m_point.at(j).at(dimension) == 0)
				{
					double g_result = (x[i] - m_point.at(j).at(i)) / ((F.at(j) + m_point.at(j).at(dimension)) + 1e-6);
					g.push_back(g_result);
				}
				else
				{
					double g_result = (x[i] - m_point.at(j).at(i)) / (F.at(j) + m_point.at(j).at(dimension));
					g.push_back(g_result);
				}
			}
			G.push_back(g);
		}

		for (int i = 0; i < G.size(); ++i)
		{
			vector<double> G_one = G.at(i);
			double result_g = 0;
			for (int j = 0; j < G_one.size()-1; ++j)
			{
				for (int k = j+1; k < G_one.size(); ++k)
				result_g += 2 * (F[j] - F[k])* (G_one[j] - G_one[k]);
			}
			g[i] = result_g;
		}

		return fx;
	}

	vector<double> solve(int dimen, vector<double> &X, vector<vector<double>>& m_point)
	{
		initial(dimen, X, m_point);
		run();
		vector<double> variables = GetVariables();
		vector<double> variablesnull;
		double result = GetFx();
	//	cout << "result" << result << endl;
		if (result < 1e-8)
			return variables;
		else
			return variablesnull;
	}

	vector<double> solve(int dimen, vector<double> &X, vector<vector<double>> &m_point, vector<double> &range)
	{
		initial(dimen, X, m_point);
		run();
		vector<double> variables = GetVariables();
		vector<double> variablesnull;
		double result = GetFx();
		
		//cout << variables.at(0) << "   " << variables.at(1) << "   " << variables.at(2) << endl;
		if (result <= 1e-8)
		{
			int num = 0;
			for (int i = 0; i < variables.size(); ++i)
			{
				double a = range.at(2 * i);
				double b = range.at(2 * i + 1);
				if ((variables.at(i) >= a - 1e-8) && (variables.at(i) <= b + 1e-8))
					num++;
			}
			if (num == variables.size())
				return  variables;
			else return  variablesnull;
		}
		else return  variablesnull;
	}
};

