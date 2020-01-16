#ifndef manhattan_H
#define manhattan_H

int manhattanProject(dVector &p,vector<vector<double>> &v, vector<double> &t);

vector<int> manhattanGenerateG(dVector &p, int k, vector<vector<double>> &v, vector<double> &t);

int fHushFunction(vector<int> &g, vector<int> &r, int n);

double manhattanDistance(vector<double> &x, vector<double> &y);

int cubeGenerateKey(dVector &p, int k, vector<vector<double>> &v, vector<double> &t);

#endif
