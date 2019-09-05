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
void fillRuleList(vector<Rule> &list, const char * rule_file, const char * fit_file, int offset);
void fillAtomList(vector<Atom> &atom_list, int shape[3], int numb_species[3], string phase_init, string spin_init, string species_init);
//float evalSiteEnergy(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
//void calcBEGParams(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, float BEG_params[]);
void init_calcJK(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void pair_calcJK(int site, int neighbor, vector<Atom> &atom_list, vector<float> J_K);
void re_calcJK(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void clacBEGParams(vector<float> &J_K);
void clacBEGParams(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K);
void clacBEGParamsNEW(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K);
float evalSiteEnergy3(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K);
float evalSiteEnergy4(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K);
float evalSiteEnergySINGLE(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K);
float evalSiteEnergyTOTAL(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K);
float evalLattice(float temp, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules,vector<float> &J_K);
float evalSiteEnergy5(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
float evalSiteEnergy6(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K);
void eval_flip(float temp, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, int site, int new_state[3], vector<float> &flip_enrgs);
void growCluster(int site, float temp, int seed_phase, int new_phase, vector<int> &links, vector<int> &cluster, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void flipCluster(int seed_phase, int new_phase, vector<Atom> &atom_list, vector<int> &cluster);
float evalCluster(vector<Atom> &atom_list, vector<int> &cluster, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K, float temp);
void runMetropolis1(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void runMetropolis2(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void runMetropolis3(int spin_passes, int cluster_passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void runMetropolis4(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void runMetropolis5(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void runMetropolis6(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);
void runMetropolis7(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules);

#endif // !mc_h