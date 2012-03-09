#ifndef HMM_H
#define HMM_H

#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>

gsl_rng *r;
double epsilon = 1e-10;

using namespace::std;

class hmm_recomb
{
	public:
//		basic objects - constructors, declarations, et c.
		hmm_recomb() {emit_probs.clear(); trans_probs.clear(); init_probs.clear(); data_progeny.clear(); data_parents.clear(); num_states = -1; num_transitions = -1;
				r = gsl_rng_alloc(gsl_rng_rand48); long seed; seed = rand()*time (NULL) * getpid();  gsl_rng_set (r, seed); rho = 0; g = 0;};
		~hmm_recomb() {};


//		basic variables	
		int num_states;
		int num_transitions;
		int num_progeny;
		double rho; 
		double g;

		vector < double > init_probs;
		vector < vector < vector < double > > > emit_probs;
		vector < vector < vector < double > > > trans_probs;

		vector < vector < double > > data_progeny;
		vector < vector < double > > data_parents;

		vector < double > snp_pos;

//		basic variables used in viterbi
		vector < vector < double > > phi;
		vector < vector < double > > delta;	

// 		basic hmm - viterbi, initialization, backward/forward, viterbi recursion
		double backward(int progeny);
		double backward(int progeny,int start,int stop);
		void viterbi(int);
		void mle(double max_rho, double max_g, int max_iter, double tol);
		vector<int> recurse_viterbi();

	
//		basic initialisations		
		void set_emit_probs(double);
		void set_num_states(int);
		void set_num_trans(int);
		void set_num_progeny();
		void set_trans_probs(double);
		void set_init_probs();

//		read data functions
		void read_snp_positions(string);
		void read_observed(string,int);

//		private:
		double gen_ran_val(double, double, double, int);
		int converged(vector<double>, int , int,double);
		int which_min( vector < vector < double > >, int, double);
		double find_min(vector < vector < double > >, int);
};


void hmm_recomb::set_num_progeny()
{
	num_progeny = data_progeny[0].size();
}

void hmm_recomb::set_init_probs()
{
	if (num_states == -1)
	{
		cout<<"No states set!"<<endl;
		exit(-1);
	}

	init_probs.clear();
	init_probs.resize(num_states,0);

	int i = 0;

	for (i = 0; i < num_states; i++)
	{
		init_probs[i] = (double)(1.0)/((double)num_states);
	}	
}

void hmm_recomb::set_num_states(int num)
{
	num_states = num;
}

void hmm_recomb::set_num_trans(int num)
{
	num_transitions = num;
}


double hmm_recomb::backward(int progeny, int start, int stop)
{

	int t, size = stop-start, i, j;
	double back = 0.0;

	vector < double > observed, dummyValues;

	vector < vector < double > > backward_values;

	observed.clear();
	observed.resize(size,0);
	dummyValues.clear();
	dummyValues.resize(size,0);
	backward_values.clear();
	backward_values.resize(num_states,dummyValues);

	for (i = 0; i < num_states; i++)
	{
		backward_values[i][size-1]= 1;
	}


	for (t =start; t<stop; t++)
	{
		observed[t-start] = data_progeny[t][progeny];	
	}

	for ( t = (stop-2); t >= start; t--)
	{
		for (i = 0; i < num_states; i++)
		{
			for (j = 0; j < num_states; j++)
			{
				backward_values[i][t-start] = backward_values[i][t-start] + trans_probs[i][j][t]*backward_values[j][t+1-start]*emit_probs[(int)observed[t+1]][j][t+1];
			}
		}

	}

	back = 0.0;

	for (i =0; i < num_states; i++)
	{
		back = back + init_probs[i]*emit_probs[(int)observed[0]][i][0]*backward_values[i][0];
	}

	return(back);
}


double hmm_recomb::backward(int progeny)
{

	int t, size = data_progeny.size(), i, j;
	double back = 0.0;

	vector < double > observed, dummyValues;

	vector < vector < double > > backward_values;
	
	observed.clear();
	observed.resize(size,0);
	dummyValues.clear();
	dummyValues.resize(size,0);
	backward_values.clear();
	backward_values.resize(num_states,dummyValues);

	for (i = 0; i < num_states; i++)
	{
		backward_values[i][size-1]= 1;
	}


	for (t =0; t<size; t++)
	{
		observed[t] = data_progeny[t][progeny];	
	}

	for ( t = (size-2); t >= 0; t--)
	{
		for (i = 0; i < num_states; i++)
		{
			for (j = 0; j < num_states; j++)
			{
				backward_values[i][t] = backward_values[i][t] + trans_probs[i][j][t]*backward_values[j][t+1]*emit_probs[(int)observed[t+1]][j][t+1];
			}
		}

	}

	back = 0.0;


	for (i =0; i < num_states; i++)
	{

		back = back + init_probs[i]*emit_probs[(int)observed[0]][i][0]*backward_values[i][0];
	}

	return(back);
}

double hmm_recomb::gen_ran_val(double current_val, double max_val, double min_val, int value)
{
	double test = gsl_ran_flat(r,0,1), x = current_val, y = 0, ran_val;

	if (test > 0.5)
	{
		y = x + gsl_ran_beta(r,1,value)*(max_val - x - epsilon);
	}
	else
	{
		y = x - gsl_ran_beta(r,1,value)*(x - min_val + epsilon);
	}
	return(y);
}

int hmm_recomb::converged(vector<double> vals, int counter, int depth,double tol)
{
	int i = 0,total = 0;
	double val;
	if (counter > depth)
	{
		val = vals[counter-1];
		for (i = 1; i < (depth-1); i++)
		{
			if (fabs(val - vals[counter - i])<tol)
			{	total = total + 1;	}
		}
		if (total == (depth - 2))
		{	return(1);	}
		else
		{	return(0);	}
	}
	else
	{
		return(0);
	}
}


void hmm_recomb::mle(double max_rho, double max_g, int max_iter, double tol)
{

	double current_rho = max_rho/10, last_rho, last_g, current_g = max_g/10, last_sum = 0, back_sum = 0, new_val,test, new_rho, new_g, value;
	int k, counter;

	vector<double> sums, mle_vals;

	sums.clear();
	sums.resize(max_iter,0);

	last_rho = current_rho;
	last_g 	 = current_g;

	new_rho = current_rho;
	new_g = current_g;

	set_trans_probs(1.0/current_rho);
	set_emit_probs(current_g);

	last_sum = 0;

	cout<<new_rho<<" "<<new_g<<" "<<k<<endl;		
	
	for (counter =1; counter < max_iter; counter++)
	{

		if (converged(sums,counter,50,0.00001)==1)
		{	
			counter = max_iter + 1;	
			break;
		}
		new_rho = current_rho;
		new_g = current_g;

		if (counter > 100)
		{
			test = gsl_ran_flat(r,0,1);
			if (test <0.5)
			{	new_rho  = gen_ran_val(current_rho, max_rho,0, 5);	}
			else
			{	new_g  = gen_ran_val(current_g, max_g,0, 5);	}
		}
		else
		if (counter <= 100)
		{
			test = gsl_ran_flat(r,0,1);
			if (test <0.5)
			{	new_rho  = gen_ran_val(current_rho, max_rho,0,5);	}
			else
			{	new_g  = gen_ran_val(current_g, max_g,0,100);	}

		}

		back_sum = 0;

		set_trans_probs(1.0/new_rho);
		set_emit_probs(new_g);

		for (k = 0; k < num_progeny; k++)
		{
			value = backward(k);
			if (value < 1e-300)
			{
				back_sum += -671;
			}
			else
			{
				back_sum += gsl_sf_log(backward(k));
			}
		}

		if (back_sum >= last_sum)
		{	
			current_rho = new_rho;
			current_g = new_g;
			last_sum = back_sum;	
		}

		sums[counter]=last_sum;
	}

	counter = 0;
	
	rho = current_rho;
	g = current_g;
}


vector<int> hmm_recomb::recurse_viterbi()
{
	int T = num_transitions+1, i;
	vector < int > q;
	q.clear(); q.resize(T,-1);

//	cout<<delta.size()<<endl;
	if (delta[0][T-1] <= delta[1][T-1])
	{	q[T-1] = 0;}
	else
	{	q[T-1] = 1;}

	for ( i = (T-2); i >= 0; i--)
	{
		q[i] = phi[q[i+1]][i+1];	
//		cout<<i<<" "<<q[i]<<endl;
	}
	return(q);
}

void hmm_recomb::viterbi(int progeny)
{
	// check function - is everything as expected?

	int T = num_transitions+1, t, i = 0, size = data_progeny.size(),j;

	double min_path = 0;
	vector<int> observed;
	vector<double> dummy;
	vector < vector < double > > path_cost, position_cost;

	observed.clear();
	observed.resize(size,0);
	dummy.clear(); dummy.resize(num_states,0);
	path_cost.clear(); path_cost.resize(num_states,dummy);

	position_cost = path_cost;

	for (t =0; t<size; t++)
	{
		observed[t] = (int)data_progeny[t][progeny];	
	}

	vector < double > dummyValues;

	dummyValues.clear();
	dummyValues.resize(T,0);

	delta.clear();
	delta.resize(num_states,dummyValues);
	
	phi.clear();
	phi.resize(num_states,dummyValues);
	
//	cout<<phi.size()<<" "<<phi[0].size()<<endl;
	for (i = 0; i < num_states; i++)
	{
		delta[i][0] = -gsl_sf_log(init_probs[i])-gsl_sf_log(emit_probs[observed[0]][i][0]);
		phi[i][0]   = 0;
	}

	for (t = 1; t<T; t++)
	{
		position_cost.clear(); position_cost.resize(num_states,dummy);	
		for (j = 0; j < num_states; j++)
		{
			for (i= 0; i < num_states; i++)
			{
				position_cost[i][j] = delta[i][t-1]-gsl_sf_log(trans_probs[i][j][t-1]) - gsl_sf_log(emit_probs[observed[t]][j][t]);
				path_cost[i][j] = delta[i][t-1] - gsl_sf_log(trans_probs[i][j][t-1]);
			}
			delta[j][t] = find_min(position_cost,j);
			min_path = find_min(path_cost,j);
			phi[j][t] = which_min(path_cost,j,min_path);
		}	
	}

}


double hmm_recomb::find_min(vector < vector < double > > position, int j)
{
	double min = 1e10;
	int i;

	for (i=0; i<position.size(); i++)
	{
		if (position[i][j] < min)
		{
			min = position[i][j];
		}
	}
	return(min);
}

int hmm_recomb::which_min(vector < vector < double > > position, int j, double min)
{
	int min_k=-1, i= 0;

	for (i=0; i<position.size(); i++)
	{
		if (position[i][j] == min)
		{
			min_k = i;
		}
	}
	return(min_k);
}


void hmm_recomb::set_emit_probs(double g)
{
	int i;
	vector < double > dummyValues;
	vector < vector < double > > dummyVector;

	if (data_parents.size()==0)
	{
		cerr<<"There are no parents..."<<endl;
		exit(-1);
	}
	if (num_states == -1)
	{
		cerr<<"There are negative number of states..."<<endl;
		exit(-1);		
	}
	
	// emit_probs:
	// color // num_states // position

	dummyValues.clear();
	dummyValues.resize(data_parents.size(),0);

	dummyVector.clear();
	dummyVector.resize(num_states,dummyValues);	

	emit_probs.clear();
	emit_probs.resize(num_states,dummyVector);

	for (i = 0; i < data_parents.size(); i++)
	{

		if (data_parents[i][0] == 0)
		{	
			emit_probs[0][0][i] = 1-g;
			emit_probs[1][0][i] = g;
		}
		else
		{
			emit_probs[0][0][i] = g;
			emit_probs[1][0][i] = 1-g;
		}
		if (data_parents[i][1] == 0)
		{	
			emit_probs[0][1][i] = 1-g;
			emit_probs[1][1][i] = g;
		}
		else
		{
			emit_probs[0][1][i] = g;
			emit_probs[1][1][i] = 1-g;
		}			
	}	
}

void hmm_recomb::read_snp_positions(string data_file)
{

	ifstream dataStream;
	dataStream.open(data_file.c_str());

	char buffer[256];
	string sample;
	int i, j,counter, check;


	string line;
	snp_pos.clear();
	cout<<"Reading position data..."<<endl;
	if (dataStream.is_open())
	{
		while(dataStream.good())
		{	
			getline(dataStream,line);
			snp_pos.push_back(atoi(line.c_str()));
		}		
	}
	else
	{
		cout<<"No data file with that name"<<endl;
		exit(-1);
	}
	snp_pos.pop_back();
	cout<<"SNP position data read."<<endl;
	set_num_trans(snp_pos.size()-1);
	dataStream.close();
	
}


void hmm_recomb::set_trans_probs(double rho)
{

	int i,j,k, num_trans = snp_pos.size()-1;
	vector < double > dummyValues;
	vector < vector < double > > dummyVector;

	if (num_states == -1)
	{
		cerr<<"There are negative number of states..."<<endl;
		exit(-1);		
	}
	
	// emit_probs:
	// num states // color // position

	dummyValues.clear();
	dummyValues.resize(snp_pos.size()-1,0);

	dummyVector.clear();
	dummyVector.resize(num_states,dummyValues);	

	trans_probs.clear();
	trans_probs.resize(num_states,dummyVector);

	for (i=0; i<num_states;i++)
	{
		for (j=0; j<num_states; j++)
		{
			for (k = 0; k < num_trans; k++)
			{
				if (i == j)
				{
					trans_probs[i][j][k] = 1 - gsl_cdf_exponential_P(double(snp_pos[k+1]-snp_pos[k]), rho);
				}
				else
				{
					trans_probs[i][j][k] = gsl_cdf_exponential_P(double(snp_pos[k+1]-snp_pos[k]), rho);
				}
			}
		}
	}
	
}

void hmm_recomb::read_observed(string data_file, int state)
{
	ifstream dataStream;
	dataStream.open(data_file.c_str());

	char buffer[256];
	string sample;
	int i, j,counter, check;

	if (state == 0)  // then it's the parents
	{
		data_parents.clear();
	}
	else
	{
		data_progeny.clear();
	}
	
	vector<double> dummyValues;	

	string line;
	cout<<"Reading data..."<<endl;
	if (dataStream.is_open())
	{
		while(dataStream.good())
		{	
			getline(dataStream,line);
			dummyValues.clear();
			for (i = 0; i < line.size(); i++)
			{

				switch(line[i])
				{
					case 48: 
						dummyValues.push_back(0);
					break;		
					case 49:
						dummyValues.push_back(1);
					break;
					case 44:
					break;		
					default:
						cout<<"OMG!!!"<<endl;
						exit(-1);
					break;
				}
			}		
			if (state == 0)
			{
				data_parents.push_back(dummyValues);
			}
			else
			{
				data_progeny.push_back(dummyValues);
			}
		}
//		
	}
	else
	{
		cout<<"No data file with that name"<<endl;
		exit(-1);
	}

	if (state ==0)
	{	data_parents.pop_back();	}
	else
	{	data_progeny.pop_back();	}
	cout<<"Data read."<<endl;
//	data_progeny.pop_back();

	if (state ==0)
	{
		set_num_trans(data_parents.size()-1);
	}
	else
	{
		set_num_trans(data_progeny.size()-1);
	}
	dataStream.close();
	//return &data_progeny;
	
}


#endif /* LINEAGEMCMC_H */
