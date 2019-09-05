#pragma once
#ifndef atom_h
#define atom_h
#include <vector>
#include <string>
#include <iostream>
using namespace std;

class Atom {
private:
	vector<int> neighbors_2;
	vector<int> neighbors_3;
	vector<int> neighbor_orders;
	vector<string> neighbor_plain_1;
	vector<string> neighbor_plain_2;
	vector<string> neighbor_plain_3;
	vector<string> neighbor_plains;
	int species;
	int spin;
	int phase;

public:
	vector<int> neighbors_1;
	vector<int> neighbors;
	float J;
	float K;
	int index;
	int pos[3];
	Atom(void);
	Atom(int _index, int _species, int _spin, int _phase, int _pos[3]);
	void setNeighbor(int _order, string _plain, int _index);
	void setSpin(int _spin);
	void setSpecies(int _species);
	void setPhase(int _phase);
	void incJK(float J_inc, float K_inc);
	int getSpin(void);
	int getSpecies(void);
	int getPhase(void);
	int getNeighbor(int _order, int _neighbor, vector<Atom> &atom_list);
	int getNeighborSpin(int _order, int _neighbor, vector<Atom> &atom_list);
	int getNeighborSpecies(int _order, int _neighbor, vector<Atom> &atom_list);
	int getNeighborPhase(int _order, int _neighbor, vector<Atom> &atom_list);
	int getNumbNeighbors(int _order);
	string getNeighborPlain(int _order, int _neighbor);

	int getNeighborSpin(int _neighbor, vector<Atom> &atom_list);
	int getNeighborSpecies(int _neighbor, vector<Atom> &atom_list);
	int getNeighborPhase(int _neighbor, vector<Atom> &atom_list);
	int getNumbNeighbors();
	int getNeighborOrder(int neighbor, vector<Atom> &atom_list);
	string getNeighborPlain(int _neighbor);
};
#endif // !atom_h
