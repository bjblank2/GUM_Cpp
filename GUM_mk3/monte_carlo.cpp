#include "monte_carlo.h"
//double randNumb(void) {
//	std::mt19937_64 rng;
//	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//	rng.seed(ss);
//	std::uniform_real_distribution<double> unif(0, 1);
//	double random_number = unif(rng);
//	return random_number;
//}
//double randNumb(mt19937_64 rng, uniform_real_distribution<double> unif) {
//	return unif(rng);
//}
int applyBC(int i, int inc, int limit) {
	int new_i;
	if (i + inc >= limit) { new_i = i + inc - limit; }
	else if (i + inc < 0) { new_i = i + inc + limit; }
	else { new_i = i + inc; }
	return new_i;
}
void fillRuleList(vector<Rule> &list, const char * rule_file, const char * fit_file) {
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
	if (rule_list.is_open()){
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
			for (int i = 0; i < rule_lines.size(); i++){
				vector<int> home_species;
				vector<int> neighbor_species;
				size_t found = rule_lines[i].find('#');
				if (found != std::string::npos) {
					name = rule_lines[i].erase(found, 2);
					string order_string = rule_lines[i + 1];
					stringstream buffer(order_string);
					buffer >> order;
					neighbor_arrangment = rule_lines[i + 2];
					istringstream iss1(rule_lines[i+3]);
					vector<string> results1(istream_iterator<string>{iss1},istream_iterator<string>());
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
					list.push_back(Rule(name,fit_vals[rule_itter],order,plain,phase,0,neighbor_arrangment,home_species,neighbor_species ));
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
	std::uniform_int_distribution<int> unif_int(0, numb_atoms-1);

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
					if (spin_rand <= 1 / 3) { spin = -1; }
					else if (spin_rand <= 2 / 3) { spin = 0; }
					else { spin = 1; }
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
					if (phase_rand <= 1 / 3) { phase = -1; }
					else if (phase_rand <= 2 / 3) { phase = 0; }
					else { phase = 1; }
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
							if (n_i == i && n_j == j && n_k == applyBC(k, 1, shape[2])){
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == j && n_k == applyBC(k, -1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i,-1,shape[0]) && n_j == j && n_k == applyBC(k, 1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == applyBC(k, -1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j,-1,shape[1]) && n_k == applyBC(k, 1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j,-1,shape[1]) && n_k == applyBC(k, -1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == applyBC(j,-1,shape[1]) && n_k == applyBC(k, 1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == applyBC(j,-1,shape[1]) && n_k == applyBC(k, -1, shape[2])) {
								atom_list[atom_index].setNeighbor(1, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, 1, shape[0]) && n_j == j && n_k == k) {
								atom_list[atom_index].setNeighbor(2, "IN", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == k) {
								atom_list[atom_index].setNeighbor(2, "IN", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j,1,shape[1]) && n_k == k) {
								atom_list[atom_index].setNeighbor(2, "IN", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j,-1,shape[1]) && n_k == k) {
								atom_list[atom_index].setNeighbor(2, "IN", neighbor_index);
							}
							if (n_i == i && n_j == j && n_k == applyBC(k,2,shape[2])) {
								atom_list[atom_index].setNeighbor(2, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == j && n_k == applyBC(k,-2,shape[2])) {
								atom_list[atom_index].setNeighbor(2, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, 1, shape[0]) && n_j == j && n_k == applyBC(k,2,shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, 1, shape[0]) && n_j == j && n_k == applyBC(k,-2,shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == applyBC(k,2,shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i, -1, shape[0]) && n_j == j && n_k == applyBC(k,-2,shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j,1,shape[1]) && n_k == applyBC(k,2,shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j,1,shape[1]) && n_k == applyBC(k,-2,shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j,-1,shape[1]) && n_k == applyBC(k,2,shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == i && n_j == applyBC(j,-1,shape[1]) && n_k == applyBC(k,-2,shape[2])) {
								atom_list[atom_index].setNeighbor(3, "OUT", neighbor_index);
							}
							if (n_i == applyBC(i,1,shape[0]) && n_j == applyBC(j, 1, shape[1]) && n_k == k) {
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
float evalSiteEnergy(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .000086173324;
	float site_energy = 0;
	int site_phase;
	int neighbor_phase;
	float BEG_params[2] = { 0, 0 };
	cout << '\n';
	calcBEGParams(site, atom_list, cluster_rules, spin_rules, BEG_params);
	cout << BEG_params[0];
	cout << BEG_params[1];
	cout << '\n';
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(1); neighbor++) {
		site_phase = atom_list[site].getPhase();
		neighbor_phase = atom_list[site].getNeighborPhase(1,neighbor,atom_list);
		site_energy += BEG_params[0] * site_phase*neighbor_phase + BEG_params[1] * (1 - site_phase ^ 2)*(1 - neighbor_phase ^ 2);
		cout << neighbor;
		cout << '\n';
	}
	//site_energy += Kb * temp * log(8)*(site_phase ^ 2);
	// add mag contribution
	return site_energy;
}
void applyRules(int site, int neighbor, int order, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, float BEG_params[]){
	int home_spin = atom_list[site].getSpin();
	int home_phase = atom_list[site].getPhase();
	int home_species = atom_list[site].getSpecies();
	int neighbor_spin = atom_list[site].getNeighborSpin(order, neighbor, atom_list);
	int neighbor_phase = atom_list[site].getNeighborPhase(order, neighbor, atom_list);
	int neighbor_species = atom_list[site].getNeighborSpecies(order, neighbor, atom_list);
	string neighbor_plain = atom_list[site].getNeighborPlain(order, neighbor);
	//cout << BEG_params[0];
	//cout << BEG_params[1];
	//cout << '\n';
	for (int i = 0; i < cluster_rules.size(); i++) {
		if (order == cluster_rules[i].getOrder()) {
			if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
				if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
					if (neighbor_plain == cluster_rules[i].getPlain()) {
						if (cluster_rules[i].getNeighborArrangment() == "PERM") {
							if (home_species != neighbor_species) {
								if (cluster_rules[i].getPhase() == 1) {
									BEG_params[0] += cluster_rules[i].getEnergyContribution();
								}
								if (cluster_rules[i].getPhase() == 0) {
									BEG_params[1] += cluster_rules[i].getEnergyContribution();
								}
							}
						}
						if (cluster_rules[i].getNeighborArrangment() == "COMB") {
							if (cluster_rules[i].getPhase() == 1) {
								BEG_params[0] += cluster_rules[i].getEnergyContribution();
							}
							if (cluster_rules[i].getPhase() == 0) {
								BEG_params[1] += cluster_rules[i].getEnergyContribution();
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < spin_rules.size(); i++) {
		if (order == spin_rules[i].getOrder()) {
			if (find(spin_rules[i].home_species.begin(), spin_rules[i].home_species.end(), home_species) != spin_rules[i].home_species.end()) {
				if (find(spin_rules[i].neighbor_species.begin(), spin_rules[i].neighbor_species.end(), neighbor_species) != spin_rules[i].neighbor_species.end()) {
					if (neighbor_plain== spin_rules[i].getPlain()) {
						string bla = spin_rules[i].getNeighborArrangment();
						if (spin_rules[i].getNeighborArrangment() == "PERM") {
							if (atom_list[site].getSpecies() != atom_list[neighbor].getSpecies()) {
								if (spin_rules[i].getPhase() == 1) {
									BEG_params[0] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
								if (spin_rules[i].getPhase() == 0) {
									BEG_params[1] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
							}
						}
						if (spin_rules[i].getNeighborArrangment() == "COMB") {
							if (spin_rules[i].getPhase() == 1) {
								BEG_params[0] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
							}
							if (spin_rules[i].getPhase() == 0) {
								BEG_params[1] += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
							}
						}
					}
				}
			}
		}
	}
}
void calcBEGParams(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules, float BEG_params[]) {
	for (int i = 0; i < atom_list[site].getNumbNeighbors(1); i++) {
		applyRules(site, i, 1, atom_list, cluster_rules, spin_rules, BEG_params);
	}
	for (int i = 0; i < atom_list[site].getNumbNeighbors(2); i++) {
		applyRules(site, i, 2, atom_list, cluster_rules, spin_rules, BEG_params);
	}
	for (int i = 0; i < atom_list[site].getNumbNeighbors(3); i++) {
		applyRules(site, i, 3, atom_list, cluster_rules, spin_rules, BEG_params);
	}
	BEG_params[0] /= 8;
	BEG_params[1] /= 8;
	cout << '\n';
	cout << BEG_params[0];
	cout << BEG_params[1];
}
vector<float> calcBEGParams(int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	int home_spin = atom_list[site].getSpin();
	int home_phase = atom_list[site].getPhase();
	int home_species = atom_list[site].getSpecies();
	float BEG_J = 0;
	float BEG_K = 0;
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(); neighbor++) {
		int neighbor_spin = atom_list[site].getNeighborSpin(neighbor, atom_list);
		int neighbor_phase = atom_list[site].getNeighborPhase(neighbor, atom_list);
		int neighbor_species = atom_list[site].getNeighborSpecies(neighbor, atom_list);
		int neighbor_order = atom_list[site].getNeighborOrder(neighbor, atom_list);
		string neighbor_plain = atom_list[site].getNeighborPlain(neighbor);
		for (int i = 0; i < cluster_rules.size(); i++) {
			if (neighbor_order == cluster_rules[i].getOrder()) {
				if (find(cluster_rules[i].home_species.begin(), cluster_rules[i].home_species.end(), home_species) != cluster_rules[i].home_species.end()) {
					if (find(cluster_rules[i].neighbor_species.begin(), cluster_rules[i].neighbor_species.end(), neighbor_species) != cluster_rules[i].neighbor_species.end()) {
						if (neighbor_plain == cluster_rules[i].getPlain()) {
							if (cluster_rules[i].getNeighborArrangment() == "PERM") {
								if (home_species != neighbor_species) {
									if (cluster_rules[i].getPhase() == 1) {
										BEG_K += cluster_rules[i].getEnergyContribution();
									}
									if (cluster_rules[i].getPhase() == 0) {
										BEG_K += cluster_rules[i].getEnergyContribution();
									}
								}
							}
							if (cluster_rules[i].getNeighborArrangment() == "COMB") {
								if (cluster_rules[i].getPhase() == 1) {
									BEG_J += cluster_rules[i].getEnergyContribution();
								}
								if (cluster_rules[i].getPhase() == 0) {
									BEG_K += cluster_rules[i].getEnergyContribution();
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
						if (neighbor_plain == spin_rules[i].getPlain()) {
							if (spin_rules[i].getNeighborArrangment() == "PERM") {
								if (atom_list[site].getSpecies() != atom_list[neighbor].getSpecies()) {
									if (spin_rules[i].getPhase() == 1) {
										BEG_J += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
									}
									if (spin_rules[i].getPhase() == 0) {
										BEG_K += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
									}
								}
							}
							if (spin_rules[i].getNeighborArrangment() == "COMB") {
								if (spin_rules[i].getPhase() == 1) {
									BEG_J += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
								if (spin_rules[i].getPhase() == 0) {
									BEG_K += spin_rules[i].getEnergyContribution()*home_spin*neighbor_spin;
								}
							}
						}
					}
				}
			}
		}
	}
	BEG_J /= 8;
	BEG_K /= 8;
	vector<float> BEG_params;
	BEG_params.push_back(BEG_J);
	BEG_params.push_back(BEG_K);
	return BEG_params;
}
float evalSiteEnergy2(float temp, int site, vector<Atom> &atom_list, vector<Rule> &cluster_rules, vector<Rule> &spin_rules) {
	float Kb = .000086173324;
	float site_energy = 0;
	int site_phase;
	int neighbor_phase;
	vector<float> BEG_params = calcBEGParams(site, atom_list, cluster_rules, spin_rules);
	for (int neighbor = 0; neighbor < atom_list[site].getNumbNeighbors(1); neighbor++) {
		site_phase = atom_list[site].getPhase();
		neighbor_phase = atom_list[site].getNeighborPhase(1, neighbor, atom_list);
		site_energy += BEG_params[0] * site_phase*neighbor_phase + BEG_params[1] * (1 - site_phase ^ 2)*(1 - neighbor_phase ^ 2);
	}
	//site_energy += Kb * temp * log(8)*(site_phase ^ 2);
	// add mag contribution
	return site_energy;
}