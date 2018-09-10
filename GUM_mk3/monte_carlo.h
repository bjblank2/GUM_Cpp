#pragma once
#ifndef mc_h
#define mc_h
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include "rule.h"
#include "atom.h"

using namespace std;

int applyBC(int i, int inc, int limit);
void fillRuleList(vector<Rule> &list, const char * rule_file, const char * fit_file);
void fillAtomList(vector<Atom> &atom_list, int shape[3], int numb_species[3], string phase_init, string spin_init, string species_init);
float evalSiteEnergy(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void calcBEGParams(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, float BEG_params[]);
void applyRules(int site, int neighbor, int order, vector<Atom> &atom_list, vector<Rule> & cluster_rules, vector<Rule> &spin_rules, float BEG_params[]);
#endif // !mc_h