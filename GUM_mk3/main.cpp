// test commit 
#include <iostream>
#include <string>
#include "atom.h"
#include "rule.h"
#include "monte_carlo.h"
using namespace std;

int main(void) {
	vector<Atom> atom_list;
	vector<Rule> cluster_rules;
	vector<Rule> spin_rules;
	int shape[3] = { 10, 10, 20 };//{ 8,8,16 }; //{ 20,20,40 }; //{ 14,14,18 }; //{ 2,2,4 }; //{ 4,4,8 };
	int species[3] = { 1000, 1000, 0 };//{ 512, 384, 128 }; //{ 8000,8000,0 };   //{ 2744,2744,0 }; //{8,6,2}; //{ 64,48,16 };
	cout << shape[0] << ',' << shape[1] << ',' << shape[2] << '\n';
	fillAtomList(atom_list, shape, species, "MART", "RAND", "RAND");
	cout << species[0] << ',' << species[1] << ',' << species[2] << '\n';
	fillRuleList(cluster_rules, "CLUSTER_RULES.txt", "FIT.txt",0);
	cout << "Testing" << '\n';
	fillRuleList(spin_rules, "SPIN_RULES.txt", "FIT.txt",cluster_rules.size());
	cout << "begining MC" << '\n';
	runMetropolis7(100,1,500,1, atom_list, cluster_rules, spin_rules);
	int exit;
	std::cin >> exit;
}