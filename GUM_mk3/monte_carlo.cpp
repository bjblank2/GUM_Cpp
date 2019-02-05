#include "monte_carlo.h"
#include <cmath>

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

void calcBEGParams(vector<float> &J_K) {
	J_K[0] = -.5;
	J_K[1] = -.4;
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
	for (int i = 0; i < cluster_rules.size(); i++) {
		if (cluster_rules[i].getOrder() == 0) {
			if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
				if (cluster_rules[i].getPhase() == 1) {
					J_K[0] += cluster_rules[i].getEnergyContribution();
				}
				if (cluster_rules[i].getPhase() == 0) {
					J_K[1] += cluster_rules[i].getEnergyContribution();
				}
			}
		}
	}
	//J_K[0] += 7;
	//J_K[1] += 7;
	//J_K[0] /= 200;
	//J_K[1] /= 200;
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
	//calcBEGParams(J_K); // Fixed J-K
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

float evalLattice(float temp, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, vector<float> &J_K) {
	float e_total = 0;
	for (int site = 0; site < atom_list.size(); site++) {
		e_total += evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
	}
	return e_total;
}

void runMetropolis(float passes, float temp1, float temp2, float temp_inc, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
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
				// Flip Phase
				bool keep = false;
				e_site_old = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
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
				e_site_new = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
				if (e_site_new <= e_site_old) {
					flip_count2 += 1;
					keep = true;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand <= keep_prob) {
						keep = true;
						flip_count += 1;
					}
					else {
						atom_list[site].setPhase(old_phase);
						keep = false;
					}
				}
				current_phase = atom_list[site].getPhase();
				phase_total += current_phase;

				// Flip Spin
				e_site_old = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
				current_J = J_K[0];
				current_K = J_K[1];
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
				e_site_new = evalSiteEnergy3(temp, site, atom_list, cluster_rules, spin_rules, J_K);
				new_J = J_K[0];
				new_K = J_K[1];
				if (e_site_new <= e_site_old) {
					e_total += e_site_new;
					flip_count2 += 1;
					keep = true;
					atom_avg_J += new_J;
					atom_avg_K += new_K;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_site_new - e_site_old));
					if (keep_rand < keep_prob) {
						e_total += e_site_new;
						flip_count += 1;
						keep = true;
						atom_avg_J += new_J;
						atom_avg_K += new_K;
					}
					else {
						atom_list[site].setSpin(old_spin);
						e_total += e_site_old;
						keep = false;
						atom_avg_J += current_J;
						atom_avg_K += current_K;
					}
				}
 				current_spin = atom_list[site].getSpin();
				//cout << atom_list[site].getSpecies();
				//cout << "[";
				//cout << current_spin;
				//cout << "]";
				spin_total2 += current_spin;
				if (atom_list[site].getSpecies() != 3) {
					for (int neighbors = 0; neighbors < 6; neighbors++) {
						spin_total += atom_list[site].getSpin() * atom_list[site].getNeighborSpin(2, neighbors, atom_list);
						//cout << " , ";
						//cout << atom_list[site].getNeighborSpecies(2, neighbors, atom_list);
					}
				}
				//cout << "\n";
				//current_spin = abs(atom_list[site].getSpin());
				//spin_total += current_spin;
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
		cout << e_avg / passes / numb_atoms * 16;// +-104.33856560606262;
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