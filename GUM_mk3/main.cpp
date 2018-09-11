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
	int shape[3] = { 4,4,8 };
	int species[3] = { 8,7,1 };
	fillAtomList(atom_list, shape, species, "AUST", "FM", "RAND");
	fillRuleList(cluster_rules, "CLUSTER_RULES.txt", "FIT.txt",0);
	fillRuleList(spin_rules, "SPIN_RULES.txt", "FIT.txt",cluster_rules.size());

	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	cout << "gen rand numbs";
	cout << '\n';
	uint64_t timestart = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	for (int i = 0; i < 10; i++) {
		double spin_rand = unif(rng);
	}
	uint64_t timestop = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	cout << timestop - timestart;
	cout << '\n';
	cout << '\n';
	cout << "clac site energy";
	cout << '\n';
	timestart = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	cout << evalSiteEnergy2(100, 0, atom_list, cluster_rules, spin_rules);
	timestop = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	cout << '\n';
	cout << '\n';
	cout << timestop - timestart;
	int exit;
	std::cin >> exit;
}