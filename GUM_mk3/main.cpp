// test commit 
#include <iostream>
#include <string>
#include "atom.h"
#include "rule.h"
#include "monte_carlo.h"
using namespace std;

int main(void) {
	vector<Atom*> atom_list;
	vector<Rule*> cluster_rules;
	vector<Rule*> spin_rules;
	int shape[3] = { 4,4,8 };//{ 20,20,40 }; //{ 14,14,18 }; //{ 2,2,4 };//{ 4,4,8 };
	int species[3] = { 64,64,0 };//{ 8000,8000,0 }; //{ 2744,2744,0 }; //{8,6,2};//{ 64,48,16 };
	fillAtomList(atom_list, shape, species, "MART", "FM", "RAND");
	fillRuleList(cluster_rules, "CLUSTER_RULES.txt", "FIT.txt",0);
	fillRuleList(spin_rules, "SPIN_RULES.txt", "FIT.txt",cluster_rules.size());
	const int length = 1;
	runMetropolis(10, 300, 10000, 5, atom_list, cluster_rules, spin_rules);
	int exit;
	std::cin >> exit;
}