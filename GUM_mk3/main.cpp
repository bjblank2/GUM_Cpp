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
	int shape[3] = {2,2,4 };//{ 8,8,16 }; //{ 20,20,40 }; //{ 14,14,18 }; //{ 2,2,4 }; //{ 4,4,8 };
	int species[3] = {8,8,0};//{ 512, 384, 128 }; //{ 8000,8000,0 };   //{ 2744,2744,0 }; //{8,6,2}; //{ 64,48,16 };
	fillAtomList(atom_list, shape, species, "MART", "STRIPED", "RAND");
	fillRuleList(cluster_rules, "CLUSTER_RULES.txt", "FIT.txt",0);
	fillRuleList(spin_rules, "SPIN_RULES.txt", "FIT.txt",cluster_rules.size());
	runMetropolis(500, 0,401, 5, atom_list, cluster_rules, spin_rules);
	int exit;
	std::cin >> exit;
}