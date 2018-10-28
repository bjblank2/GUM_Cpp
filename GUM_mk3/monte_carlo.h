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
void fillRuleList(vector<Rule*> &list, const char * rule_file, const char * fit_file, int offset);
void fillAtomList(vector<Atom*> &atom_list, int shape[3], int numb_species[3], string phase_init, string spin_init, string species_init);
//float evalSiteEnergy(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
//void calcBEGParams(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, float BEG_params[]);
vector<float> calcBEGParams(int site, vector<Atom*> &atom_list, vector<Rule*> &cluster_rules, vector<Rule*> &spin_rules);
vector<float> calcBEGParams(void);
void clacBEGParams(vector<float> &J_K, vector<Atom*> &atom_list);
void clacBEGParams(int site, vector<Atom*> &atom_list, vector<Rule*> &cluster_rules, vector<Rule*> &spin_rules, vector<float> &J_K);
//void applyRules(int site, int neighbor, int order, vector<Atom> &atom_list, vector<Rule> & cluster_rules, vector<Rule> &spin_rules, float BEG_params[]);
float evalSiteEnergy2(float temp, int site, vector<Atom*> &atom_list, vector<Rule*> &cluster_rules, vector<Rule*> &spin_rules);
float evalSiteEnergy3(float temp, int site, vector<Atom*> &atom_list, vector<Rule*> &cluster_rules, vector<Rule*> &spin_rules, vector<float> &J_K);
float evalLattice(float temp, vector<Atom*> &atom_list, vector<Rule*> &cluster_rules, vector<Rule*> &spin_rules);
void runMetropolis(float passes, float temp1, float temp2, float temp_inc, vector<Atom*> &atom_list, vector<Rule*> &cluster_rules, vector<Rule*> &spin_rules);
#endif // !mc_h