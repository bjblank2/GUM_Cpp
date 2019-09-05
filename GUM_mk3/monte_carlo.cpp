#include "monte_carlo.h"
#include <cmath>
#include <random>

int applyBC(int i, int inc, int limit) {
	int new_i;
	if (i + inc >= limit) { new_i = i + inc - limit; }
	else if (i + inc < 0) { new_i = i + inc + limit; }
	else { new_i = i + inc; }
	return new_i;
}

void fillRuleList(vector<Rule> &list, const char * rule_file, const char * fit_file, int offset) {
	int order = 0;
	string plain;
	string phase;
	int coordination = 0;
	string neighbor_arrangment;
	float energy_contribution;
	string name;
	string rule_line;
	string fit_line;
	vector<string> rule_lines;
	vector<string> fit_lines;
	vector<float> fit_vals;
	ifstream rule_list;
	rule_list.open(rule_file);
	ifstream fit_list;
	fit_list.open(fit_file);
	if (rule_list.is_open()) {
		if (fit_list.is_open()) {
			while (getline(rule_list, rule_line))
			{
				rule_lines.push_back(rule_line);
			}
			while (getline(fit_list, fit_line))
			{
				fit_lines.push_back(fit_line);
			}
			rule_list.close();
			fit_list.close();

			for (int i = 0; i < fit_lines.size(); i++) {
				size_t found = fit_lines[i].find('=');
				if (found != std::string::npos) {
					istringstream iss0(fit_lines[i]);
					vector<string> results0(istream_iterator<string>{iss0}, istream_iterator<string>());
					string fit_val = results0[2];
					stringstream buffer(fit_val);
					buffer >> energy_contribution;
					fit_vals.push_back(energy_contribution);
				}
			}
			int rule_itter = 0;
			for (int i = 0; i < rule_lines.size(); i++) {
				vector<int> home_species;
				vector<int> neighbor_species;
				size_t found = rule_lines[i].find('#');
				if (found != std::string::npos) {
					name = rule_lines[i].erase(found, 2);
					string order_string = rule_lines[i + 1];
					stringstream buffer(order_string);
					buffer >> order;
					neighbor_arrangment = rule_lines[i + 2];
					istringstream iss1(rule_lines[i + 3]);
					vector<string> results1(istream_iterator<string>{iss1}, istream_iterator<string>());
					for (int i = 0; i < results1.size(); i++) {
						string home_atoms = results1[i];
						stringstream buffer(home_atoms);
						int home_atom;
						buffer >> home_atom;
						home_species.push_back(home_atom);
					}

					istringstream iss2(rule_lines[i + 4]);
					vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
					for (int i = 0; i < results2.size(); i++) {
						string neighbor_atoms = results2[i];
						stringstream buffer(neighbor_atoms);
						int neighbor_atom;
						buffer >> neighbor_atom;
						neighbor_species.push_back(neighbor_atom);
					}
					phase = rule_lines[i + 5];
					plain = rule_lines[i + 6];
					list.push_back(Rule(name, fit_vals[rule_itter + offset], order, plain, phase, 0, neighbor_arrangment, home_species, neighbor_species));
					rule_itter += 1;
				}
			}
		}
	}
	else cout << "Unable to open file";
}

void fillAtomList(vector<Atom> &atom_list, int shape[3], int numb_species[3], string phase_init, string spin_init, string species_init) {
	int numb_atoms = shape[0] * shape[1] * shape[2];
	int atom_index = 0;
	int spin;
	int phase;
	int species;
	int pos[3];
	double spin_rand;
	double phase_rand;
	int index_rand;
	bool use_rand = false;

	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	std::mt19937_64 rng_int;
	uint64_t timeSeed_int = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss_int{ uint32_t(timeSeed_int & 0xffffffff), uint32_t(timeSeed_int >> 32) };
	rng_int.seed(ss_int);
	std::uniform_int_distribution<int> unif_int(0, numb_atoms - 1);
	int id = 0;
	for (int i = 0; i < shape[0]; i++) {
		for (int j = 0; j < shape[1]; j++) {
			for (int k = 0; k < shape[2]; k++) {
				spin_rand = unif(rng);
				phase_rand = unif(rng);
				// Set atom spin
				if (spin_init == "FM") {
					spin = 1;
				}
				else if (spin_init == "RAND") {
					if (spin_rand >= (.6666666666666666)) { spin = -1; }
					if (spin_rand >= (.3333333333333333) and spin_rand < (.6666666666666666)) { spin = 0; }
					if (spin_rand < (.3333333333333333)) { spin = 1; }
				}
				else if (spin_init == "AFM") {
					if (k % 2 == 0) {
						if ((i + j) % 2 == 0) { spin = 1; }
						else { spin = -1; }
					}
					else {
						if (((i + j + (k + 1) / 2) % 2) == 0) { spin = 1; }
						else { spin = -1; }
					}
				}
				else { spin = 1; }
				// Set atom phase
				if (phase_init == "AUST") { phase = 0; }
				else if (phase_init == "MART") { phase = 1; }
				else if (phase_init == "RAND") {
					if (phase_rand >= (.6666666666666666)) { phase = -1; }
					if (phase_rand >= (.3333333333333333) and phase_rand < (.6666666666666666)) { phase = 0; }
					if (phase_rand < (.3333333333333333)) { phase = 1; }
				}
				else { phase = 1; }
				// Set atom species
				if (k % 2 == 0) { species = 1; }
				else { species = 0; }
				if (numb_species[2] != 0) {
					if (species_init != "RAND") {
						if (numb_species[1] / numb_species[2] == 1) {
							if (k % 2 == 0) {
								if ((i + j + k / 2) % 2 == 0) { species = 2; }
								else { species = 1; }
							}
							else { species = 0; }
						}
						else {
							if (species_init == "ORDERED") {
								if ((numb_species[0] + numb_species[1] + numb_species[2]) / numb_species[2] == 8) {
									if (k % 2 == 0) {
										if ((i / 2 + (j + 1) / 2 + k / 2 % 2) == 0) { species = 2; }
										else { species = 1; }
									}
									else { species = 0; }
								}
							}
						}
					}
				}
				else { use_rand = true; }
				pos[0] = i;
				pos[1] = j;
				pos[2] = k;
				atom_list.push_back(Atom(atom_index, species, spin, phase, pos));
				atom_index += 1;
				id += 1;

			}
		}
	}
	if (species_init == "RAND" || use_rand == true) {
		int numb_comp = 0;
		while (numb_comp < numb_species[2]) {
			index_rand = unif_int(rng_int);
			if (atom_list[index_rand].getSpecies() != 0) {
				if (atom_list[index_rand].getSpecies() != 2) {
					atom_list[index_rand].setSpecies(2);
					numb_comp += 1;
				}
			}
		}
	}
	// Set neighbors
	atom_index = 0;
	for (int i = 0; i < shape[0]; i++) {
		for (int j = 0; j < shape[1]; j++) {
			for (int k = 0; k < shape[2]; k++) {
				int neighbor_index = 0;
				for (int n_i = 0; n_i < shape[0]; n_i++) {
					for (int n_j = 0; n_j < shape[1]; n_j++) {
						for (int n_k = 0; n_k < shape[2]; n_k++) {
							if (n_i == i && n_j == j && n_k == applyBC(k, 1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == j && n_k == applyBC(k, -1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == applyBC(k, 1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == applyBC(k, -1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j, -1, shape[1]) && n_k == applyBC(k, 1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j, -1, shape[1]) && n_k == applyBC(k, -1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == applyBC(j, -1, shape[1]) && n_k == applyBC(k, 1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == applyBC(j, -1, shape[1]) && n_k == applyBC(k, -1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, 1, shape[0]) && n_j == j && n_k == k) {
								atom_list[atom_index].setNeighbor(2, "IN", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == k) {
								atom_list[atom_index].setNeighbor(2, "IN", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j, 1, shape[1]) && n_k == k) {
								atom_list[atom_index].setNeighbor(2, "IN", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j, -1, shape[1]) && n_k == k) {
								atom_list[atom_index].setNeighbor(2, "IN", neighbor_index);
							}
							if (n_i == i && n_j == j && n_k == applyBC(k, 2, shape[2])) {
								atom_list[atom_index].setNeighbor(2, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == j && n_k == applyBC(k, -2, shape[2])) {
								atom_list[atom_index].setNeighbor(2, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, 1, shape[0]) && n_j == j && n_k == applyBC(k, 2, shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, 1, shape[0]) && n_j == j && n_k == applyBC(k, -2, shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == applyBC(k, 2, shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == applyBC(k, -2, shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j, 1, shape[1]) && n_k == applyBC(k, 2, shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j, 1, shape[1]) && n_k == applyBC(k, -2, shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j, -1, shape[1]) && n_k == applyBC(k, 2, shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j, -1, shape[1]) && n_k == applyBC(k, -2, shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, 1, shape[0]) && n_j == applyBC(j, 1, shape[1]) && n_k == k) {
								atom_list[atom_index].setNeighbor(3, "IN", neighbor_index);
							}
							if (n_i == applyBC(i, 1, shape[0]) && n_j == applyBC(j, -1, shape[1]) && n_k == k) {
								atom_list[atom_index].setNeighbor(3, "IN", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == applyBC(j, 1, shape[1]) && n_k == k) {
								atom_list[atom_index].setNeighbor(3, "IN", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == applyBC(-j, 1, shape[1]) && n_k == k) {
								atom_list[atom_index].setNeighbor(3, "IN", neighbor_index);
							}
							neighbor_index += 1;
						}
					}
				}
				atom_index += 1;
			}
		}
	}
}

void init_calcJK(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	int home_spin = atom_list[site].getSpin();
	int home_phase = atom_list[site].getPhase();
	int home_species = atom_list[site].getSpecies();
	int neighbor_spin;
	int neighbor_phase;
	int neighbor_species;
	int neighbor_order;
	string neighbor_plain;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		for (int i = 0; i < cluster_rules.size(); i++) {
			if (neighbor_order == cluster_rules[i].getOrder()) {
				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
						if (neighbor_plain == cluster_rules[i].getPlain() || cluster_rules[i].getPlain() == "ALL") {
							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (cluster_rules[i].getPhase() == 1) {
										atom_list[site].incJK(cluster_rules[i].getEnergyContribution(),0);
									}
									if (cluster_rules[i].getPhase() == 0) {
										atom_list[site].incJK(0,cluster_rules[i].getEnergyContribution());
									}
								}
							}
							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
								if (cluster_rules[i].getPhase() == 1) {
									atom_list[site].incJK(cluster_rules[i].getEnergyContribution(), 0);
								}
								if (cluster_rules[i].getPhase() == 0) {
									atom_list[site].incJK(0, cluster_rules[i].getEnergyContribution());
								}
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < spin_rules.size(); i++) {
			if (neighbor_order == spin_rules[i].getOrder()) {
				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
							if (spin_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (spin_rules[i].getPhase() == 1) {
										atom_list[site].incJK(cluster_rules[i].getEnergyContribution(), 0);
									}
									if (spin_rules[i].getPhase() == 0) {
										atom_list[site].incJK(0, cluster_rules[i].getEnergyContribution());
									}
								}
							}
							if (spin_rules[i].getNeighborArrangment() == "COMB") {
								if (spin_rules[i].getPhase() == 1) {
									atom_list[site].incJK(cluster_rules[i].getEnergyContribution(), 0);
								}
								if (spin_rules[i].getPhase() == 0) {
									atom_list[site].incJK(0, cluster_rules[i].getEnergyContribution());
								}
							}
						}
					}
				}
			}
		}
	}
}
void pair_calcJK(int site, int neighbor, vector<Atom> &atom_list, vector<float> J_K) {
	J_K[0] = (atom_list[site].J + atom_list[neighbor].J) / 2;
	J_K[1] = (atom_list[site].K + atom_list[neighbor].K) / 2;
}
void re_calcJK(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	atom_list[site].J = 0;
	atom_list[site].K = 0;
	int home_spin = atom_list[site].getSpin();
	int home_phase = atom_list[site].getPhase();
	int home_species = atom_list[site].getSpecies();
	int neighbor_index;
	int neighbor_spin;
	int neighbor_phase;
	int neighbor_species;
	int neighbor_order;
	string neighbor_plain;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		neighbor_index = atom_list[site].neighbors[neighbor];
		for (int i = 0; i < cluster_rules.size(); i++) {
			if (neighbor_order == cluster_rules[i].getOrder()) {
				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
						if (neighbor_plain == cluster_rules[i].getPlain() || cluster_rules[i].getPlain() == "ALL") {
							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (cluster_rules[i].getPhase() == 1) {
										atom_list[site].incJK(cluster_rules[i].getEnergyContribution(), 0);
									}
									if (cluster_rules[i].getPhase() == 0) {
										atom_list[site].incJK(0, cluster_rules[i].getEnergyContribution());
									}
								}
							}
							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
								if (cluster_rules[i].getPhase() == 1) {
									atom_list[site].incJK(cluster_rules[i].getEnergyContribution(), 0);
								}
								if (cluster_rules[i].getPhase() == 0) {
									atom_list[site].incJK(0, cluster_rules[i].getEnergyContribution());
								}
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < spin_rules.size(); i++) {
			if (neighbor_order == spin_rules[i].getOrder()) {
				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
							if (spin_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (spin_rules[i].getPhase() == 1) {
										atom_list[site].incJK(cluster_rules[i].getEnergyContribution()*home_spin*neighbor_spin, 0);
										atom_list[neighbor_index].incJK(2 * cluster_rules[i].getEnergyContribution()*home_spin*neighbor_spin, 0);
									}
									if (spin_rules[i].getPhase() == 0) {
										atom_list[site].incJK(0, cluster_rules[i].getEnergyContribution()*home_spin*neighbor_spin);
										atom_list[neighbor_index].incJK(0, 2 * cluster_rules[i].getEnergyContribution()*home_spin*neighbor_spin);
									}
								}
							}
							if (spin_rules[i].getNeighborArrangment() == "COMB") {
								if (spin_rules[i].getPhase() == 1) {
									atom_list[site].incJK(cluster_rules[i].getEnergyContribution()*home_spin*neighbor_spin, 0);
									atom_list[neighbor_index].incJK(2 * cluster_rules[i].getEnergyContribution()*home_spin*neighbor_spin, 0);
								}
								if (spin_rules[i].getPhase() == 0) {
									atom_list[site].incJK(0, cluster_rules[i].getEnergyContribution()*home_spin*neighbor_spin);
									atom_list[neighbor_index].incJK(0, 2 * cluster_rules[i].getEnergyContribution()*home_spin*neighbor_spin);
								}
							}
						}
					}
				}
			}
		}
	}
}
float evalSiteEnergy6(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	float Kb = 0.000086173324;
	float uB = .000057883818012;
	float H = 0;
	float site_energy = 0;
	int site_phase = atom_list[site].getPhase();
	int neighbor_phase;
	int neighbor_site;
	int site_spin = atom_list[site].getSpin();
	int sig1;
	int sig2;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_site = atom_list[site].neighbors[neighbor];
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		//pair_calcJK(site, neighbor_site, atom_list, J_K);
		J_K[0] = (atom_list[site].J + atom_list[neighbor_site].J) / 4 / atom_list[site].getNumbNeighbors();
		J_K[1] = (atom_list[site].K + atom_list[neighbor_site].K) / 4 / atom_list[site].getNumbNeighbors();
		sig1 = 1 - pow(site_phase, 2);
		sig2 = 1 - pow(neighbor_phase, 2);
		site_energy += J_K[0] * site_phase*neighbor_phase + J_K[1] * sig1*sig2;
	}
	//site_energy /= 8; ////////////////////////////////////////////////////////////////////////// AAAAAAAAAAAAAAAHHHHHHHHH !!!!!!!!!! ////////////
	site_energy -= Kb * temp * log(2)*(1 - pow(site_phase, 2));
	site_energy -= 3 * uB*H*site_spin;
	// add mag contribution
	return site_energy;
}

void clacBEGParams(vector<float> &J_K) {
	J_K[0] = -0.1; //-0.822868 , -0.771559
	J_K[1] = -0.035;
}

void clacBEGParams(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	int home_spin = atom_list[site].getSpin();
	int home_phase = atom_list[site].getPhase();
	int home_species = atom_list[site].getSpecies();
	int neighbor_spin;
	int neighbor_phase;
	int neighbor_species;
	int neighbor_order;
	string neighbor_plain;
	J_K[0] = 0;
	J_K[1] = 0;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		for (int i = 0; i < cluster_rules.size(); i++) {
			if (neighbor_order == cluster_rules[i].getOrder()) {
				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
						if (neighbor_plain == cluster_rules[i].getPlain() || cluster_rules[i].getPlain() == "ALL") {
							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (cluster_rules[i].getPhase() == 1) {
										J_K[0] += cluster_rules[i].getEnergyContribution();
									}
									if (cluster_rules[i].getPhase() == 0) {
										J_K[1] += cluster_rules[i].getEnergyContribution();
									}
								}
							}
							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
								if (cluster_rules[i].getPhase() == 1) {
									J_K[0] += cluster_rules[i].getEnergyContribution();
								}
								if (cluster_rules[i].getPhase() == 0) {
									J_K[1] += cluster_rules[i].getEnergyContribution();
								}
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < spin_rules.size(); i++) {
			if (neighbor_order == spin_rules[i].getOrder()) {
				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
							if (spin_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (spin_rules[i].getPhase() == 1) {
										J_K[0] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
									}
									if (spin_rules[i].getPhase() == 0) {
										J_K[1] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
									}
								}
							}
							if (spin_rules[i].getNeighborArrangment() == "COMB") {
								if (spin_rules[i].getPhase() == 1) {
									J_K[0] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
								if (spin_rules[i].getPhase() == 0) {
									J_K[1] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
							}
						}
					}
				}
			}
		}
	}
	//for (int i = 0; i < cluster_rules.size(); i++) {
	//	if (cluster_rules[i].getOrder() == 0) {
	//		if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
	//			if (cluster_rules[i].getPhase() == 1) {
	//				J_K[0] += cluster_rules[i].getEnergyContribution();
	//			}
	//			if (cluster_rules[i].getPhase() == 0) {
	//				J_K[1] += cluster_rules[i].getEnergyContribution();
	//			}
	//		}
	//	}
	//}
	J_K[0] -= .0;
	J_K[1] -= .0;
	//J_K[0] /= 200;
	//J_K[1] /= 200;
}

void clacBEGParamsNEW(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	int home_phase = atom_list[site].getPhase();
	int home_species = atom_list[site].getSpecies();
	int neighbor_spin;
	int neighbor_phase;
	int neighbor_species;
	int neighbor_order;
	string neighbor_plain;
	J_K[0] = 0;
	J_K[1] = 0;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		for (int i = 0; i < cluster_rules.size(); i++) {
			if (neighbor_order == cluster_rules[i].getOrder()) {
				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
						if (neighbor_plain == cluster_rules[i].getPlain() || cluster_rules[i].getPlain() == "ALL") {
							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (cluster_rules[i].getPhase() == 1) {
										J_K[0] += cluster_rules[i].getEnergyContribution();
									}
									if (cluster_rules[i].getPhase() == 0) {
										J_K[1] += cluster_rules[i].getEnergyContribution();
									}
								}
							}
							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
								if (cluster_rules[i].getPhase() == 1) {
									J_K[0] += cluster_rules[i].getEnergyContribution();
								}
								if (cluster_rules[i].getPhase() == 0) {
									J_K[1] += cluster_rules[i].getEnergyContribution();
								}
							}
						}
					}
				}
			}
		}
	}
	J_K[0] -= .0;
	J_K[1] -= .0;
}

float evalSiteEnergy3(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	float Kb = 0.000086173324;
	float uB = .000057883818012;
	float H = 0;
	float site_energy = 0;
	int site_phase = atom_list[site].getPhase();
	int neighbor_phase;
	int site_spin = atom_list[site].getSpin();
	int sig1;
	int sig2;
	// select wether to use fixed or on the fly J-K calcuations 
	//clacBEGParams(J_K); // Fixed J-K
	clacBEGParams(site, atom_list, cluster_rules, spin_rules, J_K);  // on the fly J-K
	for (int neighbor = 0; neighbor < 8; neighbor++) {
		neighbor_phase = atom_list[site].getNeighborPhase(1, neighbor, atom_list);
		sig1 = 1 - pow(site_phase, 2);
		sig2 = 1 - pow(neighbor_phase, 2);
		site_energy += J_K[0] * site_phase*neighbor_phase + J_K[1] * sig1*sig2;
	}
	site_energy /= 8; ////////////////////////////////////////////////////////////////////////// AAAAAAAAAAAAAAAHHHHHHHHH !!!!!!!!!! ////////////
	site_energy -= Kb * temp * log(2)*(1 - pow(site_phase, 2));
	site_energy -= 3 * uB*H*site_spin;
	// add mag contribution
	return site_energy;
}

float evalSiteEnergy4(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	float site_energy = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
	for (int neighbor = 0; neighbor < 26; neighbor++) {
		int n_site = atom_list[site].neighbors[neighbor];
		site_energy += evalSiteEnergy3(temp, n_site, atom_list, cluster_rules, spin_rules, J_K);
	}
	return site_energy/27;
}

float evalSiteEnergySINGLE(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	float Kb = 0.000086173324;
	float uB = .000057883818012;
	float H = 0;
	float site_energy = 0;
	int site_phase = atom_list[site].getPhase();
	int neighbor_phase;
	int site_spin = atom_list[site].getSpin();
	int sig1;
	int sig2;
	// select wether to use fixed or on the fly J-K calcuations 
	//calcBEGParams(J_K); // Fixed J-K
	clacBEGParamsNEW(site, atom_list, cluster_rules, spin_rules, J_K);  // on the fly J-K
	for (int neighbor = 0; neighbor < 8; neighbor++) {
		neighbor_phase = atom_list[site].getNeighborPhase(1, neighbor, atom_list);
		sig1 = 1 - pow(site_phase, 2);
		sig2 = 1 - pow(neighbor_phase, 2);
		site_energy += J_K[0] * site_phase*neighbor_phase + J_K[1] * sig1*sig2;
	}
	site_energy /= 8; ////////////////////////////////////////////////////////////////////////// AAAAAAAAAAAAAAAHHHHHHHHH !!!!!!!!!! ////////////
	site_energy -= Kb * temp * log(8)*(1 - pow(site_phase, 2));
	site_energy -= 3 * uB*H*site_spin;
	// add mag contribution
	int home_spin = atom_list[site].getSpin();
	int home_phase = atom_list[site].getPhase();
	int home_species = atom_list[site].getSpecies();
	int neighbor_spin;
	int neighbor_species;
	int neighbor_order;
	string neighbor_plain;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		for (int i = 0; i < spin_rules.size(); i++) {
			if (neighbor_order == spin_rules[i].getOrder()) {
				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
							if (spin_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (abs(home_phase) == 1) {
										if (abs(spin_rules[i].getPhase()) == 1){
											site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
										}
									}
									else if (home_phase == 0) {
										if (spin_rules[i].getPhase() == 0) {
											site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
										}
									}
								}
							}
							else if (spin_rules[i].getNeighborArrangment() == "COMB") {
								if (spin_rules[i].getPhase() == 1 == abs(home_phase)) {
									site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
								else if (spin_rules[i].getPhase() == 0 == home_phase) {
									site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
							}
						}
					}
				}
			}
		}
	}
	return site_energy;
}

float evalSiteEnergyTOTAL(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	float site_energy = evalSiteEnergySINGLE(temp, site, atom_list, cluster_rules, spin_rules, J_K);
	for (int neighbor = 0; neighbor < 26; neighbor++) {
		int n_site = atom_list[site].neighbors[neighbor];
		site_energy += evalSiteEnergySINGLE(temp, n_site, atom_list, cluster_rules, spin_rules, J_K);
	}
	int home_spin = atom_list[site].getSpin();
	int home_phase = atom_list[site].getPhase();
	int home_species = atom_list[site].getSpecies();
	int neighbor_spin;
	int neighbor_species;
	int neighbor_order;
	int neighbor_phase;
	string neighbor_plain;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		for (int i = 0; i < spin_rules.size(); i++) {
			if (neighbor_order == spin_rules[i].getOrder()) {
				if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
					if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
						if (neighbor_plain == spin_rules[i].getPlain() || spin_rules[i].getPlain() == "ALL") {
							if (spin_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (spin_rules[i].getPhase() == 1 == abs(home_phase)) {
										site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
									}
									else if (spin_rules[i].getPhase() == 0 == home_phase) {
										site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
									}
								}
							}
							if (spin_rules[i].getNeighborArrangment() == "COMB") {
								if (spin_rules[i].getPhase() == 1 == abs(home_phase)) {
									site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
								else if (spin_rules[i].getPhase() == 0 == home_phase) {
									site_energy += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
							}
						}
					}
				}
			}
		}
	}
	return site_energy;
}

float evalLattice(float temp, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	float e_total = 0;
	for (int site = 0; site < atom_list.size(); site++) {
		e_total += evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
	}
	return e_total / atom_list.size() * 16;
}

void runMetropolis1(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old = 0;
	float e_site_new = 0;
	float spin_rand = 0;
	float phase_rand = 0;
	float keep_rand = 0;
	int old_phase = 0;
	int new_phase = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	int current_phase = 0;
	bool phase_same;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float phase_total = 0;
	float phase_avg = 0;
	float phase_total_abs = 0;
	float phase_avg_abs = 0;
	float keep_prob = 0;
	int numb_atoms = size(atom_list);
	float current_J;
	float current_K;
	float new_J;
	float new_K;
	float atom_avg_J;
	float atom_avg_K;
	float pass_avg_J;
	float pass_avg_K;
	int flip_count = 0;
	int flip_count2 = 0;
	vector<float> J_K = { 0,0 };
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	cout << evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K);
	cout << "\n";
	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		phase_avg = 0;
		phase_avg_abs = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		pass_avg_J = 0;
		pass_avg_K = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			phase_total = 0;
			phase_total_abs = 0;

			spin_total = 0;
			spin_total2 = 0;
			atom_avg_J = 0;
			atom_avg_K = 0;
			for (int site = 0; site < atom_list.size(); site++) {
				// Flip Phase
				bool keep = false;
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_old = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);////////////////////////          //////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				old_phase = atom_list[site].getPhase();
				phase_same = true;
				while (phase_same == true) {
					phase_rand = unif(rng);
					if (phase_rand <= 0.3333333333333333) {
						new_phase = -1;
					}
					else if (phase_rand <= 0.6666666666666666) {
						new_phase = 0;
					}
					else {
						new_phase = 1;
					}
					if (new_phase != old_phase) { phase_same = false; }
				}
				atom_list[site].setPhase(new_phase);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_new = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);////////////////////////         ///////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (e_site_new <= e_site_old) {
					flip_count2 += 1;
					keep = true;
					e_total += e_site_new;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand <= keep_prob) {
						keep = true;
						e_total += e_site_new;
						flip_count += 1;
					}
					else {
						atom_list[site].setPhase(old_phase);
						keep = false;
						e_total += e_site_old;
					}
				}
				current_phase = atom_list[site].getPhase();
				phase_total += current_phase;
				phase_total_abs += abs(current_phase);

				atom_avg_J += J_K[0];
				atom_avg_K += J_K[1];
				// Flip Spin
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//e_site_old = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);                   ////////////////
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//current_J = J_K[0];
				//current_K = J_K[1];
				//old_spin = atom_list[site].getSpin();
				//spin_same = true;
				//while (spin_same == true) {
				//	spin_rand = unif(rng);
				//	if (spin_rand <= 0.3333333333333333) {
				//		new_spin = -1;
				//	}
				//	else if (spin_rand <= 0.6666666666666666) {
				//		new_spin = 0;
				//	}
				//	else {
				//		new_spin = 1;
				//	}
				//	if (new_spin != old_spin) { spin_same = false; }
				//}
				//atom_list[site].setSpin(new_spin);
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//e_site_new = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);                     //////////////
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//new_J = J_K[0];
				//new_K = J_K[1];
				//if (e_site_new <= e_site_old) {
				//	e_total += e_site_new;
				//	flip_count2 += 1;
				//	keep = true;
				//	atom_avg_J += new_J;
				//	atom_avg_K += new_K;
				//}
				//else {
				//	keep_rand = unif(rng);
				//	keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
				//	if (keep_rand < keep_prob) {
				//		e_total += e_site_new;
				//		flip_count += 1;
				//		keep = true;
				//		atom_avg_J += new_J;
				//		atom_avg_K += new_K;
				//	}
				//	else {
				//		atom_list[site].setSpin(old_spin);
				//		e_total += e_site_old;
				//		keep = false;
				//		atom_avg_J += current_J;
				//		atom_avg_K += current_K;
				//	}
				//}
				//current_spin = atom_list[site].getSpin();
				//spin_total2 += current_spin;
				//if (atom_list[site].getSpecies() != 0) {
				//	for (int neighbors = 0; neighbors < 6; neighbors++) {
				//		spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
				//	}
				//}
			}
			phase_avg += phase_total;
			phase_avg_abs += phase_total_abs;
			spin_avg += spin_total;
			spin_avg2 += spin_total2;
			e_avg += e_total;
			pass_avg_J += atom_avg_J;
			pass_avg_K += atom_avg_K;
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / passes / numb_atoms * 16;
		cout << " , ";
		cout << phase_avg / passes / numb_atoms;
		cout << " , ";
		cout << phase_avg_abs / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms / 6 * 2;
		cout << " , ";
		cout << spin_avg2 / passes / numb_atoms;
		cout << " , ";
		cout << pass_avg_J / passes / numb_atoms;
		cout << " , ";
		cout << pass_avg_K / passes / numb_atoms;
		cout << " , ";
		cout << pass_avg_K / pass_avg_J;
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}

void runMetropolis2(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old = 0;
	float e_site_new = 0;
	float spin_rand = 0;
	float phase_rand = 0;
	float keep_rand = 0;
	int old_phase = 0;
	int new_phase = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	int current_phase = 0;
	bool both_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float phase_total = 0;
	float phase_avg = 0;
	float keep_prob = 0;
	int numb_atoms = size(atom_list);
	float current_J;
	float current_K;
	float new_J;
	float new_K;
	float atom_avg_J;
	float atom_avg_K;
	float pass_avg_J;
	float pass_avg_K;
	int flip_count = 0;
	int flip_count2 = 0;
	vector<float> J_K = { 0,0 };
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	float e_lattice = evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K); //only for debug
	cout << evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K);
	cout << "\n";
	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		phase_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		pass_avg_J = 0;
		pass_avg_K = 0;
		flip_count = 0;
		flip_count2 = 0;
		//e_total = evalLattice(temp, atom_list, cluster_rules, spin_rules);
		for (int i = 0; i < passes; i++) {
			//e_total = evalLattice(temp, atom_list, cluster_rules, spin_rules);
			e_total = 0;
			phase_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			atom_avg_J = 0;
			atom_avg_K = 0;
			for (int site = 0; site < atom_list.size(); site++) {
				//cout << atom_list[site].getSpecies();
				// Flip Phase and Spin
				bool keep = false;
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_old = evalSiteEnergy4(temp, site, atom_list, cluster_rules, spin_rules, J_K);/////////////////////////////MADE IT 4!!!!!!!
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				current_J = J_K[0];
				current_K = J_K[1];
				old_phase = atom_list[site].getPhase();
				old_spin = atom_list[site].getSpin();
				both_same = true;
				while (both_same == true) {
					phase_rand = unif(rng);
					spin_rand = unif(rng);
					if (phase_rand <= 0.3333333333333333) {
						new_phase = -1;
					}
					else if (phase_rand <= 0.6666666666666666) {
						new_phase = 0;
					}
					else {
						new_phase = 1;
					}
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					if ((old_phase == new_phase) and (old_spin == new_spin)) { both_same = true; }
					else { both_same = false; }
				}
				atom_list[site].setSpin(new_spin);
				atom_list[site].setPhase(new_phase);
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_new = evalSiteEnergy4(temp, site, atom_list, cluster_rules, spin_rules, J_K);//MADE IT 4!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				new_J = J_K[0];
				new_K = J_K[1];
				if (e_site_new <= e_site_old) {
					flip_count2 += 1;
					keep = true;
					e_total += e_site_new;
					atom_avg_J += new_J;
					atom_avg_K += new_K;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand <= keep_prob) {
						keep = true;
						flip_count += 1;
						e_total += e_site_new;
						atom_avg_J += new_J;
						atom_avg_K += new_K;
					}
					else {
						atom_list[site].setPhase(old_phase);
						atom_list[site].setSpin(old_spin);
						keep = false;
						e_total += e_site_old;
						atom_avg_J += current_J;
						atom_avg_K += current_K;
					}
				}
				current_phase = atom_list[site].getPhase();
				phase_total += abs(current_phase);
				current_spin = atom_list[site].getSpin();
				if (atom_list[site].getSpecies() != 0) {
					for (int neighbors = 0; neighbors < 6; neighbors++) {
						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
					}
				}
			}
			phase_avg += phase_total;
			spin_avg += spin_total;
			spin_avg2 += spin_total2;
			e_avg += e_total;
			pass_avg_J += atom_avg_J;
			pass_avg_K += atom_avg_K;
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / passes / numb_atoms * 16;
		cout << " , ";
		cout << phase_avg / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms / 6 * 2;
		cout << " , ";
		cout << spin_avg2 / passes / numb_atoms;
		cout << " , ";
		cout << pass_avg_J / passes / numb_atoms;
		cout << " , ";
		cout << pass_avg_K / passes / numb_atoms;
		cout << " , ";
		cout << pass_avg_K / pass_avg_J;
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}

void runMetropolis3(int spin_passes, int cluster_passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old = 0;
	float e_site_new = 0;
	float spin_rand = 0;
	float phase_rand = 0;
	float keep_rand = 0;
	int old_phase = 0;
	int new_phase = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	int current_phase = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float phase_total = 0;
	float phase_avg = 0;
	float keep_prob = 0;
	int numb_atoms = size(atom_list);
	float atom_avg_J;
	float atom_avg_K;
	float pass_avg_J;
	float pass_avg_K;
	int flip_count = 0;
	int flip_count2 = 0;
	float phase_avg_abs = 0;
	float phase_total_abs = 0;
	vector<int> cluster;
	vector<int> links;
	float rand;
	float H_cluster_old;
	float H_cluster_new;
	float prob;
	int seed_site;
	int seed_phase;
	bool phase_same;
	vector<float> J_K = { 0,0 };
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng_i(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_int_distribution<int> uni(0, atom_list.size()); // guaranteed unbiased
	cout << evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K);
	cout << "\n";
	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		phase_avg = 0;
		phase_avg_abs = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		pass_avg_J = 0;
		pass_avg_K = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < spin_passes; i++) {
			e_total = 0;
			phase_total = 0;
			phase_total_abs = 0;

			spin_total = 0;
			spin_total2 = 0;
			atom_avg_J = 0;
			atom_avg_K = 0;
			for (int site = 0; site < atom_list.size(); site++) {
				// Flip Spin
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_old = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);                   ////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				old_spin = atom_list[site].getSpin();
				spin_same = true;
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_list[site].setSpin(new_spin);
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_new = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);                     //////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (e_site_new <= e_site_old) {
					e_total += e_site_new;
					flip_count2 += 1;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand < keep_prob) {
						e_total += e_site_new;
						flip_count += 1;
					}
					else {
						atom_list[site].setSpin(old_spin);
						e_total += e_site_old;
					}
				}
				current_spin = atom_list[site].getSpin();
				spin_total2 += current_spin;
				if (atom_list[site].getSpecies() != 0) {
					for (int neighbors = 0; neighbors < 6; neighbors++) {
						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
					}
				}
			}
			phase_avg += phase_total;
			phase_avg_abs += phase_total_abs;
			spin_avg += spin_total;
			spin_avg2 += spin_total2;
			e_avg += e_total;
			pass_avg_J += atom_avg_J;
			pass_avg_K += atom_avg_K;
		}
		
		seed_site = uni(rng_i);
		seed_phase = atom_list[seed_site].getPhase();
		phase_rand = unif(rng);
		phase_same = true;
		while (phase_same == true)
		{
			if (phase_rand <= 0.3333333333333333) {
				new_phase = -1;
			}
			else if (phase_rand <= 0.6666666666666666) {
				new_phase = 0;
			}
			else {
				new_phase = 1;
			}
			if (new_phase == seed_phase) { phase_same = true; }
			else { phase_same = false; }
		}
		growCluster(seed_site, temp, seed_phase, new_phase, links, cluster, atom_list, cluster_rules, spin_rules);
		if (seed_phase*new_phase == -1) {
			flipCluster(seed_phase, new_phase, atom_list, cluster);
		}
		else {
			H_cluster_old = evalCluster(atom_list, cluster, cluster_rules, spin_rules, J_K, temp);
			flipCluster(seed_phase, new_phase, atom_list, cluster);
			H_cluster_new = evalCluster(atom_list, cluster, cluster_rules, spin_rules, J_K, temp);
			if (H_cluster_new <= H_cluster_old) {
				cout << "\accepting MC cluster flip: new energy < old energy";
			}
			else {
				rand = unif(rng);
				prob = exp(-1 / (Kb*temp)*(H_cluster_new - H_cluster_old));
				if (rand < prob) {
					cout << "\accepting MC cluster flip";
				}
				else {
					cout << "rejecting MC cluster flip";
					flipCluster(new_phase, seed_phase, atom_list, cluster);
				}
			}
		}
		e_total = evalLattice(temp, atom_list, cluster_rules, spin_rules, J_K);
	}
}

void growCluster(int site, float temp, int seed_phase, int new_phase, vector<int> &links, vector<int> &cluster, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	vector<float> J_K = { 0,0 };
	float Kb = .000086173324;
	float B = 1 / (Kb*temp);
	int site_phase = atom_list[site].getPhase();
	int proposed_link;
	int new_site;
	clacBEGParams(site, atom_list, cluster_rules, spin_rules, J_K);
	float BEG_K = -2 * B*J_K[0];
	float BEG_M = -2 * B*J_K[1];
	float rand;
	float prob;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// // Wolff Algorithm
	if (new_phase*seed_phase == -1) {
		for (int neighbor = 0; neighbor < 8; neighbor++) {
			if (atom_list[site].getNeighborPhase(1, neighbor, atom_list) == seed_phase) {
				proposed_link = atom_list[site].getNeighbor(1, neighbor, atom_list);
				if (std::find(links.begin(), links.end(), proposed_link) == links.end()) {
					rand = unif(rng);
					prob = 1 - exp(-2 * BEG_K);
					if (rand <= prob) {
						new_site = atom_list[site].getNeighbor(1,neighbor,atom_list);
						growCluster(new_site, temp, seed_phase, new_phase, links, cluster, atom_list, cluster_rules, spin_rules);
					}
				}
			}
		}
	}
	// // Mixed Cluster Algorithm
	if ((seed_phase == 1 and new_phase == 0) or (seed_phase == 0 and new_phase == -1)) {
		for (int neighbor = 0; neighbor < 8; neighbor++) {
			proposed_link = atom_list[site].getNeighbor(1, neighbor, atom_list);
			if (std::find(links.begin(), links.end(), proposed_link) == links.end()) {
				if (atom_list[site].getNeighborPhase(1, neighbor, atom_list) == 1 or atom_list[site].getNeighborPhase(1, neighbor, atom_list) == 0) {
					if (atom_list[site].getNeighborPhase(1, neighbor, atom_list) == site_phase) {
						rand = unif(rng);
						prob = 1 - exp(-BEG_K - BEG_M / 3);
						if (rand < prob) {
							new_site = atom_list[site].getNeighbor(1, neighbor, atom_list);
							growCluster(new_site, temp, seed_phase, new_phase, links, cluster, atom_list, cluster_rules, spin_rules);
						}
					}
					else {
						rand = unif(rng);
						prob = 1 - exp(-BEG_K + BEG_M / 3);
						if (rand < prob) {
							new_site = atom_list[site].getNeighbor(1, neighbor, atom_list);
							growCluster(new_site, temp, seed_phase, new_phase, links, cluster, atom_list, cluster_rules, spin_rules);
						}
					}
				}
			}
		}
	}
	if ((seed_phase == -1 and new_phase == 0) or (seed_phase == 0 and new_phase == 0, 1)) {
		for (int neighbor = 0; neighbor < 8; neighbor++) {
			proposed_link = atom_list[site].getNeighbor(1, neighbor, atom_list);
			if (std::find(links.begin(), links.end(), proposed_link) == links.end()) {
				if (atom_list[site].getNeighborPhase(1, neighbor, atom_list) == -1 or atom_list[site].getNeighborPhase(1, neighbor, atom_list) == 0) {
					if (atom_list[site].getNeighborPhase(1, neighbor, atom_list) == site_phase) {
						rand = unif(rng);
						prob = 1 - exp(-BEG_K - BEG_M / 3);
						if (rand < prob) {
							new_site = atom_list[site].getNeighbor(1, neighbor, atom_list);
							growCluster(new_site, temp, seed_phase, new_phase, links, cluster, atom_list, cluster_rules, spin_rules);
						}
					}
					else {
						rand = unif(rng);
						prob = 1 - exp(-BEG_K + BEG_M / 3);
						if (rand < prob) {
							new_site = atom_list[site].getNeighbor(1, neighbor, atom_list);
							growCluster(new_site, temp, seed_phase, new_phase, links, cluster, atom_list, cluster_rules, spin_rules);
						}
					}
				}
			}
		}
	}
}

float evalCluster(vector<Atom> &atom_list, vector<int> &cluster, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K, float temp) {
	float Kb = .000086173324;
	float total_H = 0;
	int site;
	float total_H_inc;
	float inc_count = 0;
	int site_phase;
	int neighbor_phase;
	int proposed_link = 0;
	if (cluster.size() == 1) {
		site = cluster[0];
		total_H = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
	}
	else {
		for (int i = 0; i < cluster.size(); i++) {
			site = cluster[i];
			site_phase = atom_list[site].getPhase();
			clacBEGParams(site, atom_list, cluster_rules, spin_rules, J_K);
			total_H_inc = 0;
			inc_count = 0;
			for (int neighbor = 0; neighbor < 8; neighbor++) {
				proposed_link = atom_list[site].getNeighbor(1, neighbor, atom_list);
				if (std::find(cluster.begin(), cluster.end(), proposed_link) == cluster.end()) {
					neighbor_phase = atom_list[site].getNeighborPhase(1, neighbor, atom_list);
					total_H_inc += J_K[0] * site_phase*neighbor_phase + J_K[1] * (1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
					inc_count += 1;
				}
			}
			total_H += (total_H_inc / inc_count + Kb * temp*log(2)*pow(site_phase,2));
		}
	}
	return total_H;
}

void flipCluster(int seed_phase, int new_phase, vector<Atom> &atom_list, vector<int> &cluster) {
	int site;
	int old_phase;
	if (seed_phase*new_phase == -1) {
		for (int i = 0; i < cluster.size(); i++) {
			site = cluster[i];
			atom_list[site].setPhase(new_phase);
		}
	}
	else {
		for (int i = 0; i <= cluster.size(); i++) {
			site = cluster[i];
			if ((seed_phase == 1 and new_phase == 0) or (seed_phase == 0 and new_phase == -1)) {
				old_phase = atom_list[site].getPhase();
				if (old_phase == 1) {
					atom_list[site].setPhase(0);
				}
				else if (old_phase == 0) {
					atom_list[site].setPhase(-1);
				}
				if ((seed_phase == -1 and new_phase == 0) or (seed_phase == 0 and new_phase == 1)) {
					old_phase = atom_list[site].getPhase();
					if (old_phase == -1) {
						atom_list[site].setPhase(0);
					}
					else if (old_phase == 0) {
						atom_list[site].setPhase(1);
					}
				}
			}
		}
	}
}

void runMetropolis4(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old = 0;
	float e_site_new = 0;
	float spin_rand = 0;
	float phase_rand = 0;
	float keep_rand = 0;
	int old_phase = 0;
	int new_phase = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	int current_phase = 0;
	bool both_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float phase_total = 0;
	float phase_avg = 0;
	float keep_prob = 0;
	int numb_atoms = size(atom_list);
	float current_J;
	float current_K;
	float new_J;
	float new_K;
	float atom_avg_J;
	float atom_avg_K;
	float pass_avg_J;
	float pass_avg_K;
	int flip_count = 0;
	int flip_count2 = 0;
	vector<float> J_K = { 0,0 };
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	float e_lattice = evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K); //only for debug
	cout << evalLattice(temp1, atom_list, cluster_rules, spin_rules, J_K);
	cout << "\n";
	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		phase_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		pass_avg_J = 0;
		pass_avg_K = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			phase_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			atom_avg_J = 0;
			atom_avg_K = 0;
			for (int site = 0; site < atom_list.size(); site++) {
				// Flip Phase and Spin
				bool keep = false;
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_old = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);/////////////MADE IT 4!!!!!!!
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				current_J = J_K[0];
				current_K = J_K[1];
				old_phase = atom_list[site].getPhase();
				old_spin = atom_list[site].getSpin();
				both_same = true;
				while (both_same == true) {
					phase_rand = unif(rng);
					spin_rand = unif(rng);
					if (phase_rand <= 0.3333333333333333) {
						new_phase = -1;
					}
					else if (phase_rand <= 0.6666666666666666) {
						new_phase = 0;
					}
					else {
						new_phase = 1;
					}
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					if ((old_phase == new_phase) and (old_spin == new_spin)) { both_same = true; }
					else { both_same = false; }
				}
				atom_list[site].setSpin(new_spin);
				atom_list[site].setPhase(new_phase);
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_new = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);/////////////MADE IT 4!!!!!!!
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (e_site_new <= e_site_old) {
					flip_count2 += 1;
					keep = true;
					e_total += e_site_new;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand <= keep_prob) {
						keep = true;
						flip_count += 1;
						e_total += e_site_new;
					}
					else {
						atom_list[site].setPhase(old_phase);
						atom_list[site].setSpin(old_spin);
						keep = false;
						e_total += e_site_old;
					}
				}
				current_phase = atom_list[site].getPhase();
				phase_total += current_phase;
				current_spin = atom_list[site].getSpin();
				if (atom_list[site].getSpecies() != 0) {
					for (int neighbors = 0; neighbors < 6; neighbors++) {
						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
					}
				}
			}
			phase_avg += phase_total;
			spin_avg += spin_total;
			spin_avg2 += spin_total2;
			e_avg += e_total;
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / passes / numb_atoms * 16;
		cout << " , ";
		cout << phase_avg / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms / 6 * 2;
		cout << " , ";
		cout << spin_avg2 / passes / numb_atoms;
		cout << " , ";
		cout << J_K[0];
		cout << " , ";
		cout << J_K[1];
		cout << " , ";
		cout << J_K[1] / J_K[0];
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}

void eval_flip(float temp, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, int site, int new_state[3], vector<float> &flip_enrgs) {
	float Kb = .0000861733035;
	int site_phase = atom_list[site].getPhase();
	int site_spin = atom_list[site].getSpin();
	int site_species = atom_list[site].getSpecies();
	int new_phase = new_state[0];
	int new_spin = new_state[1];
	int new_species = new_state[2];
	int neighbor_phase;
	int neighbor_spin;
	int neighbor_species;
	int neighbor_order;
	string neighbor_plain;

	float current_enrg = 0;
	float new_enrg = 0;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		for (int clust_rule = 0; clust_rule < cluster_rules.size(); clust_rule++) {
			if (neighbor_order == cluster_rules[clust_rule].getOrder()) {
				if (find(cluster_rules[clust_rule].home_species.begin(), cluster_rules[clust_rule].home_species.end(), site_species) != cluster_rules[clust_rule].home_species.end()) {
					if (find(cluster_rules[clust_rule].neighbor_species.begin(), cluster_rules[clust_rule].neighbor_species.end(), neighbor_species) != cluster_rules[clust_rule].neighbor_species.end()) {
						if (neighbor_plain == cluster_rules[clust_rule].getPlain() || cluster_rules[clust_rule].getPlain() == "ALL") {
							if (cluster_rules[clust_rule].getNeighborArrangment() == "PERM") {
								if (site_species != neighbor_species) {
									if (cluster_rules[clust_rule].getPhase() == 1) {
										current_enrg += cluster_rules[clust_rule].getEnergyContribution()*site_phase*neighbor_phase;
										new_enrg += cluster_rules[clust_rule].getEnergyContribution()*new_phase*neighbor_phase;
									}
									else if (cluster_rules[clust_rule].getPhase() == 0) {
										current_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
										new_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(new_phase, 2))*(1 - pow(neighbor_phase, 2));
									}
								}
							}
							else if (cluster_rules[clust_rule].getNeighborArrangment() == "COMB") {
								if (cluster_rules[clust_rule].getPhase() == 1) {
									current_enrg += cluster_rules[clust_rule].getEnergyContribution()*site_phase*neighbor_phase;
									new_enrg += cluster_rules[clust_rule].getEnergyContribution()*new_phase*neighbor_phase;
								}
								else if (cluster_rules[clust_rule].getPhase() == 0) {
									current_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
									new_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(new_phase, 2))*(1 - pow(neighbor_phase, 2));
								}
							}
						}
					}
				}
			}
		}
		for (int spin_rule = 0; spin_rule < spin_rules.size(); spin_rule++) {
			if (neighbor_order == spin_rules[spin_rule].getOrder()) {
				if (find(spin_rules[spin_rule].home_species.begin(), spin_rules[spin_rule].home_species.end(), site_species) != spin_rules[spin_rule].home_species.end()) {
					if (find(spin_rules[spin_rule].neighbor_species.begin(), spin_rules[spin_rule].neighbor_species.end(), neighbor_species) != spin_rules[spin_rule].neighbor_species.end()) {
						if (neighbor_plain == spin_rules[spin_rule].getPlain() || spin_rules[spin_rule].getPlain() == "ALL") {
							if (spin_rules[spin_rule].getNeighborArrangment() == "PERM") {
								if (site_species != neighbor_species) {
									if (spin_rules[spin_rule].getPhase() == 1 and abs(site_phase) == 1) {
										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
									}
									if (spin_rules[spin_rule].getPhase() == 1 and abs(new_phase) == 1) {
										new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin/2;
									}
									if (spin_rules[spin_rule].getPhase() == 0 and site_phase == 0) {
										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
									}
									if (spin_rules[spin_rule].getPhase() == 0 and new_phase == 0) {
										new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin/2;
									}
									if (spin_rules[spin_rule].getPhase() == 1 and abs(neighbor_phase) == 1) {
										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin/2;
										new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin/2;
									}
									if (spin_rules[spin_rule].getPhase() == 0 and neighbor_phase == 0) {
										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin/2;
										new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin/2;
									}
								}
							}
							else if (spin_rules[spin_rule].getNeighborArrangment() == "COMB") {
								if (spin_rules[spin_rule].getPhase() == 1 and abs(site_phase) == 1) {
									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
								}
								if (spin_rules[spin_rule].getPhase() == 1 and abs(new_phase) == 1) {
									new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin/2;
								}
								if (spin_rules[spin_rule].getPhase() == 0 and site_phase == 0) {
									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
								}
								if (spin_rules[spin_rule].getPhase() == 0 and new_phase == 0) {
									new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin/2;
								}
								if (spin_rules[spin_rule].getPhase() == 1 and abs(neighbor_phase) == 1) {
									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin/2;
									new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin/2;
								}
								if (spin_rules[spin_rule].getPhase() == 0 and neighbor_phase == 0) {
									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin/2;
									new_enrg += spin_rules[spin_rule].getEnergyContribution()*new_spin*neighbor_spin/2;
								}
							}
						}
					}
				}
			}
		}
	}
	current_enrg += Kb * temp*log(8)*(1 - pow(site_phase, 2));
	new_enrg += Kb * temp*log(8)*(1 - pow(site_phase, 2));
	flip_enrgs[0] = current_enrg;
	flip_enrgs[1] = new_enrg;
}

float evalSiteEnergy5(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .0000861733035;
	int site_phase = atom_list[site].getPhase();
	int site_spin = atom_list[site].getSpin();
	int site_species = atom_list[site].getSpecies();
	int neighbor_phase;
	int neighbor_spin;
	int neighbor_species;
	int neighbor_order;
	string neighbor_plain;

	float current_enrg = 0;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		for (int clust_rule = 0; clust_rule < cluster_rules.size(); clust_rule++) {
			if (neighbor_order == cluster_rules[clust_rule].getOrder()) {
				if (find(cluster_rules[clust_rule].home_species.begin(), cluster_rules[clust_rule].home_species.end(), site_species) != cluster_rules[clust_rule].home_species.end()) {
					if (find(cluster_rules[clust_rule].neighbor_species.begin(), cluster_rules[clust_rule].neighbor_species.end(), neighbor_species) != cluster_rules[clust_rule].neighbor_species.end()) {
						if (neighbor_plain == cluster_rules[clust_rule].getPlain() || cluster_rules[clust_rule].getPlain() == "ALL") {
							if (cluster_rules[clust_rule].getNeighborArrangment() == "PERM") {
								if (site_species != neighbor_species) {
									if (cluster_rules[clust_rule].getPhase() == 1) {
										current_enrg += cluster_rules[clust_rule].getEnergyContribution()*site_phase*neighbor_phase;
									}
									else if (cluster_rules[clust_rule].getPhase() == 0) {
										current_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
									}
								}
							}
							else if (cluster_rules[clust_rule].getNeighborArrangment() == "COMB") {
								if (cluster_rules[clust_rule].getPhase() == 1) {
									current_enrg += cluster_rules[clust_rule].getEnergyContribution()*site_phase*neighbor_phase;
								}
								else if (cluster_rules[clust_rule].getPhase() == 0) {
									current_enrg += cluster_rules[clust_rule].getEnergyContribution()*(1 - pow(site_phase, 2))*(1 - pow(neighbor_phase, 2));
								}
							}
						}
					}
				}
			}
		}
		for (int spin_rule = 0; spin_rule < spin_rules.size(); spin_rule++) {
			if (neighbor_order == spin_rules[spin_rule].getOrder()) {
				if (find(spin_rules[spin_rule].home_species.begin(), spin_rules[spin_rule].home_species.end(), site_species) != spin_rules[spin_rule].home_species.end()) {
					if (find(spin_rules[spin_rule].neighbor_species.begin(), spin_rules[spin_rule].neighbor_species.end(), neighbor_species) != spin_rules[spin_rule].neighbor_species.end()) {
						if (neighbor_plain == spin_rules[spin_rule].getPlain() || spin_rules[spin_rule].getPlain() == "ALL") {
							if (spin_rules[spin_rule].getNeighborArrangment() == "PERM") {
								if (site_species != neighbor_species) {
									if (spin_rules[spin_rule].getPhase() == 1 and abs(site_phase) == 1) {
										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
									}
									if (spin_rules[spin_rule].getPhase() == 0 and site_phase == 0) {
										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
									}
									if (spin_rules[spin_rule].getPhase() == 1 and abs(neighbor_phase) == 1) {
										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
									}
									if (spin_rules[spin_rule].getPhase() == 0 and neighbor_phase == 0) {
										current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
									}
								}
							}
							else if (spin_rules[spin_rule].getNeighborArrangment() == "COMB") {
								if (spin_rules[spin_rule].getPhase() == 1 and abs(site_phase) == 1) {
									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
								}
								if (spin_rules[spin_rule].getPhase() == 0 and site_phase == 0) {
									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
								}
								if (spin_rules[spin_rule].getPhase() == 1 and abs(neighbor_phase) == 1) {
									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
								}
								if (spin_rules[spin_rule].getPhase() == 0 and neighbor_phase == 0) {
									current_enrg += spin_rules[spin_rule].getEnergyContribution()*site_spin*neighbor_spin / 2;
								}
							}
						}
					}
				}
			}
		}
	}
	current_enrg += Kb * temp*log(8)*(1 - pow(site_phase, 2));
	return current_enrg;
}

void runMetropolis5(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old;
	float e_site_new;
	vector<float> flip_enrgs = { 0, 0 };
	float spin_rand = 0;
	float phase_rand = 0;
	float keep_rand = 0;
	int old_phase = 0;
	int new_phase = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	int current_phase = 0;
	int new_state[3];
	bool both_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float phase_total = 0;
	float phase_avg = 0;
	float phase_total_abs = 0;
	float phase_avg_abs = 0;
	float keep_prob = 0;
	int numb_atoms = size(atom_list);
	int flip_count = 0;
	int flip_count2 = 0;

	int check_phase = 0;
	int check_spin = 0;
	int check_species = 0;
	bool keep;

	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	cout << "\n";
	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		phase_avg = 0;
		phase_avg_abs = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			phase_total = 0;
			phase_total_abs = 0;
			spin_total = 0;
			spin_total2 = 0;
			for (int site = 0; site < atom_list.size(); site++) {
				// Flip Phase and Spin
				keep = false;
				old_phase = atom_list[site].getPhase();
				old_spin = atom_list[site].getSpin();
				both_same = true;
				while (both_same == true) {
					phase_rand = unif(rng);
					spin_rand = unif(rng);
					if (phase_rand <= 0.3333333333333333) {
						new_phase = -1;
					}
					else if (phase_rand <= 0.6666666666666666) {
						new_phase = 0;
					}
					else {
						new_phase = 1;
					}
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					new_state[0] = new_phase;
					new_state[1] = new_spin;
					if ((old_phase == new_phase) and (old_spin == new_spin)) { both_same = true; }
					else { both_same = false; }
				}
				new_state[2] = 9; // just a standin dummy in case I want to add species at some point. Not acctually used
				eval_flip(temp, atom_list, cluster_rules, spin_rules, site, new_state, flip_enrgs);
				e_site_old = float(flip_enrgs[0]);
				e_site_new = float(flip_enrgs[1]);
				/*cout << e_site_old;
				cout << ' ';
				cout << e_site_new;
				cout << " : ";*/
				check_phase = atom_list[site].getPhase();
				check_spin = atom_list[site].getSpin();
				check_species = atom_list[site].getSpecies();
				if (e_site_new <= e_site_old) {
					flip_count2 += 1;
					keep = true;
					atom_list[site].setSpin(new_spin);
					atom_list[site].setPhase(new_phase);
					e_total +=e_site_new;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand <= keep_prob) {
						keep = true;
						flip_count += 1;
						atom_list[site].setSpin(new_spin);
						atom_list[site].setPhase(new_phase);
						e_total += e_site_new;
					}
					else {
						//atom_list[site].setPhase(old_phase);
						//atom_list[site].setSpin(old_spin);
						keep = false;
						e_total += e_site_old;
					}
				}
				/*cout << e_site_old;
				cout << ' ';
				cout << e_site_new;
				cout << "\n";*/
 				current_phase = atom_list[site].getPhase();
				phase_total += current_phase;
				phase_total_abs += abs(current_phase);
				current_spin = atom_list[site].getSpin();
				if (atom_list[site].getSpecies() != 0) {
					for (int neighbors = 0; neighbors < 6; neighbors++) {
						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
					}
				}
			}
			phase_avg += phase_total;
			phase_avg_abs += phase_total_abs;
			spin_avg += spin_total;
			spin_avg2 += spin_total2;
			e_avg += e_total;
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / passes / numb_atoms * 16;
		cout << " , ";
		cout << phase_avg / passes / numb_atoms;
		cout << " , ";
		cout << phase_avg_abs / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms / 6 * 2;
		cout << " , ";
		cout << spin_avg2 / passes / numb_atoms;
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}

void runMetropolis6(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	////////////////////////////////////////////////////////
	// MC with two site "spiral" moves
	////////////////////////////////////////////////////////
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old = 0;
	float e_site_new = 0;
	float spin_rand = 0;
	float phase_rand = 0;
	float keep_rand = 0;
	int old_phase = 0;
	int new_phase = 0;
	int old_spin = 0;
	int new_spin = 0;
	int old_phase2 = 0;
	int new_phase2 = 0;
	int old_spin2 = 0;
	int new_spin2 = 0;
	int current_spin = 0;
	int current_phase = 0;
	bool all_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float phase_total = 0;
	float phase_avg = 0;
	float keep_prob = 0;
	int numb_atoms = size(atom_list);
	float atom_avg_J;
	float atom_avg_K;
	float pass_avg_J;
	float pass_avg_K;
	int flip_count = 0;
	int flip_count2 = 0;
	int site2;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	cout << "\n";
	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		phase_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		pass_avg_J = 0;
		pass_avg_K = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < int(passes/8); i++) {
			e_total = 0;
			phase_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			atom_avg_J = 0;
			atom_avg_K = 0;
			for (int site = 0; site < atom_list.size(); site++) {
				for (int neig = 0; neig < 8; neig++) {
					// Flip Phase and Spin
					site2 = atom_list[site].neighbors_1[neig];
					bool keep = false;
					old_phase = atom_list[site].getPhase();
					old_spin = atom_list[site].getSpin();
					old_phase2 = atom_list[site2].getPhase();
					old_spin2 = atom_list[site2].getSpin();
					all_same = true;
					while (all_same == true) {
						phase_rand = unif(rng);
						spin_rand = unif(rng);
						if (phase_rand <= 0.3333333333333333) {
							new_phase = -1;
						}
						else if (phase_rand <= 0.6666666666666666) {
							new_phase = 0;
						}
						else {
							new_phase = 1;
						}
						if (spin_rand <= 0.3333333333333333) {
							new_spin = -1;
						}
						else if (spin_rand <= 0.6666666666666666) {
							new_spin = 0;
						}
						else {
							new_spin = 1;
						}

						new_phase2 = new_phase;
						spin_rand = unif(rng);
						if (spin_rand <= 0.3333333333333333) {
							new_spin2 = -1;
						}
						else if (spin_rand <= 0.6666666666666666) {
							new_spin2 = 0;
						}
						else {
							new_spin2 = 1;
						}

						if ((old_phase == new_phase) and (old_spin == new_spin) and (old_spin2 == new_spin2)) { all_same = true; }
						else { all_same = false; }
					}
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					e_site_old = evalSiteEnergy5(temp, site, atom_list, cluster_rules, spin_rules);                        //////////
					e_site_old += evalSiteEnergy5(temp, site2, atom_list, cluster_rules, spin_rules);                      //////////
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					atom_list[site].setSpin(new_spin);
					atom_list[site].setPhase(new_phase);
					atom_list[site2].setSpin(new_spin2);
					atom_list[site2].setPhase(new_phase2);
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					e_site_new = evalSiteEnergy5(temp, site, atom_list, cluster_rules, spin_rules);                        //////////
					e_site_new += evalSiteEnergy5(temp, site2, atom_list, cluster_rules, spin_rules);                      //////////
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					if (e_site_new <= e_site_old) {
						flip_count2 += 1;
						keep = true;
						e_total += e_site_new;
					}
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
						if (keep_rand <= keep_prob) {
							keep = true;
							flip_count += 1;
							e_total += e_site_new;
						}
						else {
							atom_list[site].setPhase(old_phase);
							atom_list[site].setSpin(old_spin);
							atom_list[site2].setPhase(old_phase2);
							atom_list[site2].setSpin(old_phase2);
							keep = false;
							e_total += e_site_old;
						}
					}
					current_phase = atom_list[site].getPhase()+atom_list[site2].getPhase();
					phase_total += current_phase;
					current_spin = atom_list[site].getSpin()+atom_list[site2].getSpin();
					if (atom_list[site].getSpecies() != 0) {
						for (int neighbors = 0; neighbors < 6; neighbors++) {
							spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
						}
					}
					/*cout << site;
					cout << " ,";
					cout << atom_list[site].getSpecies();
					cout << " ,";
					cout << atom_list[site].getPhase();
					cout << " ,";
					cout << atom_list[site].getSpin();
					cout << "\n";*/
				}
			}
			phase_avg += phase_total/8;
			spin_avg += spin_total;
			spin_avg2 += spin_total2/8;
			e_avg += e_total/8;
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / passes / numb_atoms * 16;
		cout << " , ";
		cout << phase_avg / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms / 6 * 2;
		cout << " , ";
		cout << spin_avg2 / passes / numb_atoms;
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}

void runMetropolis7(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .0000861733035;
	float e_total = 0;
	float e_site_old = 0;
	float e_site_new = 0;
	float spin_rand = 0;
	float phase_rand = 0;
	float keep_rand = 0;
	int old_phase = 0;
	int new_phase = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	int current_phase = 0;
	bool both_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float phase_total = 0;
	float phase_avg = 0;
	float keep_prob = 0;
	int numb_atoms = size(atom_list);
	float current_J;
	float current_K;
	float new_J;
	float new_K;
	float atom_avg_J;
	float atom_avg_K;
	float pass_avg_J;
	float pass_avg_K;
	int flip_count = 0;
	int flip_count2 = 0;
	vector<float> J_K = { 0,0 };
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	cout << "\n";

	for (int site = 0; site < atom_list.size(); site++) {
		init_calcJK(site, atom_list, cluster_rules, spin_rules);
	}

	for (float temp = temp1; temp < temp2; temp += temp_inc) {
		e_avg = 0;
		phase_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		pass_avg_J = 0;
		pass_avg_K = 0;
		flip_count = 0;
		flip_count2 = 0;
		//e_total = evalLattice(temp, atom_list, cluster_rules, spin_rules);
		for (int i = 0; i < passes; i++) {
			//e_total = evalLattice(temp, atom_list, cluster_rules, spin_rules);
			e_total = 0;
			phase_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			atom_avg_J = 0;
			atom_avg_K = 0;
			pass_avg_J = 0;
			pass_avg_K = 0;
			for (int site = 0; site < atom_list.size(); site++) {
				//cout << atom_list[site].getSpecies();
				// Flip Phase and Spin
				bool keep = false;
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_old = evalSiteEnergy6(temp, site, atom_list, cluster_rules, spin_rules, J_K);                                           //
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				current_J = J_K[0];
				current_K = J_K[1];
				old_phase = atom_list[site].getPhase();
				old_spin = atom_list[site].getSpin();
				both_same = true;
				while (both_same == true) {
					phase_rand = unif(rng);
					spin_rand = unif(rng);
					if (phase_rand <= 0.3333333333333333) {
						new_phase = -1;
					}
					else if (phase_rand <= 0.6666666666666666) {
						new_phase = 0;
					}
					else {
						new_phase = 1;
					}
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					if ((old_phase == new_phase) and (old_spin == new_spin)) { both_same = true; }
					else { both_same = false; }
				}
				atom_list[site].setSpin(new_spin);
				atom_list[site].setPhase(new_phase);
				if (old_spin != new_spin) {
					re_calcJK(site, atom_list, cluster_rules, spin_rules);
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				e_site_new = evalSiteEnergy6(temp, site, atom_list, cluster_rules, spin_rules, J_K);                                     //
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				new_J = J_K[0];
				new_K = J_K[1];
				if (e_site_new <= e_site_old) {
					flip_count2 += 1;
					keep = true;
					e_total += e_site_new;
					atom_avg_J += new_J;
					atom_avg_K += new_K;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand <= keep_prob) {
						keep = true;
						flip_count += 1;
						e_total += e_site_new;
						atom_avg_J += new_J;
						atom_avg_K += new_K;
					}
					else {
						atom_list[site].setPhase(old_phase);
						atom_list[site].setSpin(old_spin);
						if (old_spin != new_spin) {
							re_calcJK(site, atom_list, cluster_rules, spin_rules);
						}
						keep = false;
						e_total += e_site_old;
						atom_avg_J += current_J;
						atom_avg_K += current_K;
					}
				}
				current_phase = atom_list[site].getPhase();
				phase_total += abs(current_phase);
				current_spin = atom_list[site].getSpin();
				if (atom_list[site].getSpecies() != 0) {
					for (int neighbors = 0; neighbors < 6; neighbors++) {
						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
					}
				}
			}
			phase_avg += phase_total;
			spin_avg += spin_total;
			spin_avg2 += spin_total2;
			e_avg += e_total;
			pass_avg_J += atom_avg_J;
			pass_avg_K += atom_avg_K;
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / passes / numb_atoms * 16;
		cout << " , ";
		cout << phase_avg / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms / 6 * 2;
		cout << " , ";
		cout << spin_avg2 / passes / numb_atoms;
		cout << " , ";
		cout << pass_avg_J / passes / numb_atoms;
		cout << " , ";
		cout << pass_avg_K / passes / numb_atoms;
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}