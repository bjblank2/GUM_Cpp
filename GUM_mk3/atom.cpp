#include "atom.h"

Atom::Atom(void) {
	species = 10;
	spin = 10;
	phase = 10;
	index = 10;
}

Atom::Atom(int _index, int _species, int _spin, int _phase, int _pos[3]) {
	index = _index;
	species = _species;
	spin = _spin;
	phase = _phase;
	pos[0] = _pos[0];
	pos[1] = _pos[1];
	pos[2] = _pos[2];
}

void Atom::setNeighbor(int _order, string _plain, int _index) {
	if (_order == 1) {
		neighbors_1.push_back(_index);
		neighbor_plain_1.push_back(_plain);
	}
	if (_order == 2) {
		neighbors_2.push_back(_index);
		neighbor_plain_2.push_back(_plain);
	}
	if (_order == 3) {
		neighbors_3.push_back(_index);
		neighbor_plain_3.push_back(_plain);
	}
}

int Atom::getNeighborSpin(int _order, int _neighbor, vector<Atom> atom_list) {
	int neighbor_index;
	int neighbor_spin;
	if (_order == 1) {
		neighbor_index = neighbors_1[_neighbor];
		neighbor_spin = atom_list[neighbor_index].getSpin();
	}
	else if (_order == 2) {
		neighbor_index = neighbors_2[_neighbor];
		neighbor_spin = atom_list[neighbor_index].getSpin();
	}
	else if (_order == 3) {
		neighbor_index = neighbors_3[_neighbor];
		neighbor_spin = atom_list[neighbor_index].getSpin();
	}
	else {
		cout << "error";
		cout << '\n';
		neighbor_spin = 1000000000;
	}
	return neighbor_spin;
}
int Atom::getNeighborSpecies(int _order, int _neighbor, vector<Atom> atom_list) {
	int neighbor_index;
	int neighbor_species;
	if (_order == 1) {
		neighbor_index = neighbors_1[_neighbor];
		neighbor_species = atom_list[neighbor_index].getSpecies();
	}
	else if (_order == 2) {
		neighbor_index = neighbors_2[_neighbor];
		neighbor_species = atom_list[neighbor_index].getSpecies();
	}
	else if (_order == 3) {
		neighbor_index = neighbors_3[_neighbor];
		neighbor_species = atom_list[neighbor_index].getSpecies();
	}
	else {
		cout << "error";
		neighbor_species = 1000000000;
	}
	return neighbor_species;
}
int Atom::getNumbNeighbors(int _order) {
	int numb_neighbors = 0;
	if (_order == 1) { numb_neighbors = neighbors_1.size(); }
	else if (_order == 2) { numb_neighbors = neighbors_2.size(); }
	else if (_order == 3) { numb_neighbors = neighbors_3.size(); }
	return numb_neighbors;
}
int Atom::getNeighborPhase(int _order, int _neighbor, vector<Atom> atom_list) {
	int neighbor_index;
	int neighbor_phase;
	if (_order == 1) {
		neighbor_index = neighbors_1[_neighbor];
		neighbor_phase = atom_list[neighbor_index].getPhase();
	}
	else if (_order == 2) {
		neighbor_index = neighbors_2[_neighbor];
		neighbor_phase = atom_list[neighbor_index].getPhase();
	}
	else if (_order == 3) {
		neighbor_index = neighbors_3[_neighbor];
		neighbor_phase = atom_list[neighbor_index].getPhase();
	}
	else {
		cout << "error";
		neighbor_phase = 1000000000;
	}
	return neighbor_phase;
}
string Atom::getNeighborPlain(int _order, int _neighbor) {
	string neighbor_plain;
	if (_order == 1) {
		neighbor_plain = neighbor_plain_1 [_neighbor];
	}
	else if (_order == 2) {
		neighbor_plain = neighbor_plain_2[_neighbor];
	}
	else if (_order == 3) {
		neighbor_plain = neighbor_plain_3[_neighbor];
	}
	else {
		cout << "error";
		neighbor_plain = "ERROR";
	}
	return neighbor_plain;
}

void Atom::setSpin(int _spin) {
	spin = _spin;
}
void Atom::setSpecies(int _species) {
	species = _species;
}
void Atom::setPhase(int _phase) {
	phase = _phase;
}

int Atom::getSpin(void) {
	return spin;
}
int Atom::getSpecies(void) {
	return species;
}
int Atom::getPhase(void) {
	return phase;
}