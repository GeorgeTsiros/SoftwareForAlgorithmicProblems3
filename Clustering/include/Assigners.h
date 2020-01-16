#include <string>
#include <vector>
#include "Library.h"
#ifndef DATABASE
#define DATABASE

#include "Database.h"

#endif

using namespace std;

template <class Point>
class Assigner {
protected:
    int K;
    int Grids;
    int L;
    int k;
    double w;
public:
    Assigner(){}
    virtual vector<int>** assign(vector<vector<Point>>*, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*) {return NULL;}
    virtual string get_name() {}
    virtual int get_K() {return K;}
};

template <class Point>
class Lloyd_assignment : public Assigner<Point> {
private:
    string name = "Lloyd's Assignment";
public:
    Lloyd_assignment(int K, int Grids, int L, int k){this->K = K; this->Grids = Grids; this->L = L; this->k = k;}
    vector<int>** assign(vector<vector<Point>>*, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*);
    string get_name();
};

template <class Point>
class Inverse_assignment : public Assigner<Point> {
private:
    string name = "Inverse Assignment";
public:
    Inverse_assignment(int K, int Grids, int L, int k, double w){this->K = K; this->Grids = Grids; this->L = L; this->k = k; this->w = w;}
    vector<int>** assign(vector<vector<Point>>*, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*);
    string get_name();
};