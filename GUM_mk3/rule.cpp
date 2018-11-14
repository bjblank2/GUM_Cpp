#include "rule.h"
 
Rule::Rule(void) {
	order = 10;
}
Rule::Rule(string _name, float _energy_contribution, int _order, string _plain, string _phase, int _coordination, string _neighbor_arrangment, vector<int> _home_species, vector<int> _neighbor_species) {
	name = _name;
	energy_contribution = _energy_contribution;
	order = _order;
	plain = _plain;
	phase = _phase;
	coordination = _coordination;
	neighbor_arrangment = _neighbor_arrangment;
	home_species = _home_species;
	neighbor_species = _neighbor_species;
}
float Rule::getEnergyContribution() {
	return energy_contribution;
}
int Rule::getOrder() {
	return order;
}
int Rule::getPhase() {
	if (phase == "mart" || phase == "MART") { return 1; }
	if (phase == "aust" || phase == "AUST") { return 0; }
}
string Rule::getNeighborArrangment() {
	return neighbor_arrangment;
}
string Rule::getPlain() {
	return plain;
}