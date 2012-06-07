#include <math.h>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include "gsl/gsl_sf.h"
#include "hmm.h"

string ToString(int value)
{
	stringstream ss;
	ss << value;
	return ss.str();
}

int main(int argc, char* argv[])
{

	if (argc != 4)
	{
		cerr<<"HMM requires 3 input variables: path, recombination rate, and genotype rate!"<<endl;
		exit(-1);
	}
	hmm_recomb test;
	int j = 0;
	std::string out_path_file, out_para_file, progeny_file, position_file, parent_file, in_1, in_2;
	ofstream fout, gout;
	
	int i = 0, k = 0;
	vector<int> this_path;

	//out_para_file = in_1+"para.txt";
	cout<<out_para_file<<endl;
	gout.open(out_para_file.c_str());
	for (j = 1; j <= 14; j++)
	{
		in_1 = argv[1];
		parent_file = in_1 + "parents.MAL"+ ToString(j);
		cout<<parent_file <<endl;
		progeny_file = in_1 + "progeny.MAL"+ ToString(j);
		cout<<progeny_file <<endl;
		position_file = in_1 + "positions.MAL"+ ToString(j);
		cout<<position_file <<endl;
		out_path_file = in_1+"paths.MAL"+ToString(j);

//		cout<<parent_file<<endl;
		test.read_observed(parent_file,0);
		test.read_observed(progeny_file,1);
		test.read_snp_positions(position_file);	
	
		fout.open(out_path_file.c_str());


		cout<<test.data_progeny.size()<<" "<<test.data_progeny[0].size()<<endl;

		test.set_num_states(2);

		test.set_init_probs();
		test.set_num_progeny();

		cout<<test.num_progeny<<endl;

	//	test.mle(0.01,0.1, 1000,1e-7);

		gout<<test.rho<<","<<test.g<<endl;

		test.set_trans_probs(1/3.e-6);
		test.set_emit_probs(0.003);
		for (k = 0; k < test.num_progeny;k++)
		{
			test.viterbi(k);
			this_path = test.recurse_viterbi();

			for (i = 0; i < this_path.size(); i++)
			{
				fout<<this_path[i];
				if (i < this_path.size()-1)
				{	fout<<",";	}
			}
			fout<<endl;
		}
		fout.close();
	}
	gout.close();
	return(0);

}
