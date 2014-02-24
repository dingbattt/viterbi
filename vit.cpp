#include "vit.h"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
using namespace std;

// g++ -I /home/frederik/Boost/ vit.h -std=c++11 vit.cpp -o vit.o && ./vit.o
map<char, int> obs_map;

vit::vit(string f){
	ifstream mf(f.c_str());
	if(mf.is_open()){
		string line;
		while(getline(mf, line)){
			if (( line.compare("states") ) == 0){
				getline(mf, line);
				state_size = stoi(line.c_str());
				getline(mf, line);
				vector<string> string_of_states;
				boost::algorithm::split(string_of_states, line, boost::algorithm::is_any_of(" "));
				if (string_of_states.size()==state_size){
					for (int r = 0; r<string_of_states.size(); r++){
						states.push_back(atoi(string_of_states[r].c_str()));
					}
				}
			}
			if (( line.compare("observables") ) == 0){
				getline(mf, line);
				observables_size = stoi(line.c_str());
				getline(mf, line);
				boost::algorithm::split(observables, line, boost::algorithm::is_any_of(" "));
				int idx = 0;
				for(auto obs : observables){
					obs_map[obs[0]] = idx;
					idx++;
				}
			}
			if ((line.compare("initProbs")) == 0){
				getline(mf, line);
				vector<string> init_vector; 
				boost::algorithm::split(init_vector, line, boost::algorithm::is_any_of(" "));
				for (int r = 0; r<init_vector.size(); r++){
					pi.push_back(atof(init_vector[r].c_str()));
				}
			}
			if ((line.compare("transProbs")) == 0){
				for(int s=0; s<state_size; s++){
					vector<string> trans_vector; 
					getline(mf, line);
					boost::algorithm::split(trans_vector, line, boost::algorithm::is_any_of(" "));
					vector<double> temp; 
					A.push_back(temp);
					for (int r = 0; r<trans_vector.size(); r++){
						A[s].push_back(atof(trans_vector[r].c_str()));
					}
				}
			}
			if ((line.compare("emProbs")) ==0){
				for(int s=0; s<state_size; s++){
					vector<string> em_vector; 
					getline(mf, line);
					boost::algorithm::split(em_vector, line, boost::algorithm::is_any_of(" "));
					vector<double> temp; 
					phi.push_back(temp);
					for (int r = 0; r<em_vector.size(); r++){
						phi[s].push_back(atof(em_vector[r].c_str()));
					}
				}
			}
		}
	}
}

vit::~vit(){

}

vector<vector<double> > vit::viterbi(string str){
	original = str; 
	int X = str.size();
	char input[X];
	strcpy(input, str.c_str());
	vector< vector<double> > res;
	vector<double> temp (X);
	for(int i=0; i<state_size; i++){
		res.push_back(temp);
	}
	int viter[X]; 
	for(int xs=0; xs<X; xs++){
		pair<int, double> tmp (0,0);
		for(int st=0; st<state_size; st++){
			if(xs>0){
				vector<double> pc;
				for(int i=0; i<state_size; i++){
					pc.push_back(res[i][(xs-1)]);
				}
				res[st][xs]=p(obs_map[input[xs]], st, pc);
			} else {
				res[st][xs]=pb(obs_map[input[xs]], st);
			}
			if(res[st][xs] > tmp.second){
				tmp.second = res[st][xs];
				tmp.first = (st+1);
			}
		}
		viter[xs]=tmp.first;
	}
	return res;
}

vector<vector<double> > vit::log_viterbi(string str){
	original = str; 
	int X = str.size();
	char input[X];
	strcpy(input, str.c_str());
	vector< vector<double> > res;
	vector<double> temp (X);
	for(int i=0; i<state_size; i++){
		res.push_back(temp);
	}
	int viter[X]; 
	for(int xs=0; xs<X; xs++){
		pair<int, double> tmp (0,0);
		for(int st=0; st<state_size; st++){
			if(xs>0){
				vector<double> pc;
				for(int i=0; i<state_size; i++){
					pc.push_back(res[i][(xs-1)]);
				}
				res[st][xs]=p(obs_map[input[xs]], st, pc);
			} else {
				res[st][xs]=pb(obs_map[input[xs]], st);
			}
			if(res[st][xs] > tmp.second){
				tmp.second = res[st][xs];
				tmp.first = (st+1);
			}
		}
		viter[xs]=tmp.first;
	}
	return res;
}

double vit::pb(int input, int state){
	double init_p = pi[state]; 
	double emission_p = phi[state][input];
	return (init_p * emission_p);
}

double vit::pb_log(int input, int state){
	double tmp = 0;
	for(int k = 1; k<=state_skze; k++){
		tmp = A[k][1] * ()
	}
	double init_p = pi[state]; 
	double emission_p = phi[state][input];
	return (init_p * emission_p);
}

double vit::p(int xn, int j, vector<double> prev_calc){
	double max_so_far = 0;
	for(int k=0; k<prev_calc.size(); k++){
		double init_p = prev_calc[k];
		double emission_p = phi[j][xn];
		double trans_p = A[k][j];
		double res = (init_p * emission_p * trans_p);
		if (res>max_so_far) {max_so_far = res;} 
	}
	return max_so_far;
}

double vit::p_log(int xn, int j, vector<double> prev_calc){
	double max_so_far = 0;
	for(int k=0; k<prev_calc.size(); k++){
		double init_p = prev_calc[k];
		double emission_p = phi[j][xn];
		double trans_p = A[k][j];
		double res = (init_p * emission_p * trans_p);
		if (res>max_so_far) {max_so_far = res;} 
	}
	return max_so_far;
}

string vit::backtrack(vector<vector<double> > input){
	int z_size = input[0].size();
	vector<pair<int, double> > z(z_size);
	pair<int, double> tmp (0,-1);
	for (int i=0; i<state_size; i++){
		if (tmp.second<input[i][z_size-1]){
			tmp.second=input[i][z_size-1];
			tmp.first = i+1;
		}
	}
	z[z_size-1] = tmp;
	for(int i=z_size-2; i>=0; i--){
		tmp.second = -1;
		tmp.first = 0;
		for(int k=0; k<state_size; k++){
			for(int t=0; t<state_size; t++){
				double prob = phi[z[i+1].first][obs_map[original[i+1]]] * input[t][i] * A[k][z[i+1].first];
				if (tmp.second < prob){
					tmp.second = prob;
					tmp.first = t+1;
				}
			}
		}
		z[i] = tmp; 
	}
	for(pair<int, double> foo : z){
		cout << foo.first; 
	}
	cout << endl;
	string res; 
	for (pair<int, double> foo : z){
		res = res + to_string(foo.first);
	}
	return res; 
}

string vit::log_backtrack(vector<vector<double> > input){
	int z_size = input[0].size();
	vector<pair<int, double> > z(z_size);
	pair<int, double> tmp (0,-1);
	for (int i=0; i<state_size; i++){
		if (tmp.second<input[i][z_size-1]){
			tmp.second=input[i][z_size-1];
			tmp.first = i+1;
		}
	}
	z[z_size-1] = tmp;
	for(int i=z_size-2; i>=0; i--){
		tmp.second = -1;
		tmp.first = 0;
		for(int k=0; k<state_size; k++){
			for(int t=0; t<state_size; t++){
				double prob = log(phi[z[i+1].first][obs_map[original[i+1]]]) * input[t][i] * log(A[k][z[i+1].first]);
				if (tmp.second < prob){
					tmp.second = prob;
					tmp.first = t+1;
				}
			}
		}
		z[i] = tmp; 
	}
	for(pair<int, double> foo : z){
		cout << foo.first; 
	}
	cout << endl;
	string res; 
	for (pair<int, double> foo : z){
		res = res + to_string(foo.first);
	}
	return res; 
}

int main(){
	auto v = new vit("vit_input.txt");
	auto f = v->viterbi("GTTTCCCAGTGTATATCGAGGGATACTACGTGCATAGTAACATCGGCCAA");
	v->backtrack(f);
}