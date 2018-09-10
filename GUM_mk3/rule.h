#pragma once
#ifndef rule_h
#define rule_h
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;

class Rule {
private:
	int order;
	string plain;
	string phase;
	int coordination;
	string neighbor_arrangment;
	float energy_contribution;
	string name;

public:
	vector<int> home_species;
	vector<int> neighbor_species;
	Rule(void);
	Rule(string _name, float energy_contribution, int _order, string _plain, string _phase, int _coordination, string _neighbor_arrangment, vector<int> _home_species, vector<int> _neighbor_species);
	float getEnergyContribution(void);
	void setEnergyContribution(float input);
	string getPlain();
	int getOrder();
	int getPhase();
	string getNeighborArrangment();
	// need to add atom first... float applyRule()
};

#endif