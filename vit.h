#ifndef VIT_H
#define VIT_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

class vit {
private: 
	int state_size, observables_size; 
	std::vector<int> states;
	std::vector<std::string> observables; 
	std::vector<double> pi; //initProbs
	std::vector<std::vector<double> > A, phi;// transProbs and emProbs
	std::string original; 
	double pb(int input, int state); 
	double p(int input, int state, std::vector<double> prev_calc); 
	double pb_log(int input, int state); 
	double p_log(int input, int state, std::vector<double> prev_calc); 
public:
	vit(std::string f);
	~vit();
	std::string backtrack(std::vector<std::vector<double> > input);
	std::string log_backtrack(std::vector<std::vector<double> > input);
	std::vector<std::vector<double> > viterbi(std::string str);
	std::vector<std::vector<double> > log_viterbi(std::string str);
};

#endif