#include <string>
#include <vector>
#include "Library.h"
#ifndef DATABASE
#define DATABASE

#include "Database.h"

#endif

using namespace std;

template <class Point>
class Updater {
protected:
    int K;
public:
    Updater(){}
    virtual int update(vector<vector<Point>>*, vector<int>**, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*) {return 0;}
    virtual string get_name() {}
    virtual int get_K() {return K;}
};

template <class Point>
class PAM : public Updater<Point> {
private:
    string name = "Partitioning Around Medoids (PAM)";
public:
    PAM(int K) {this->K = K;};
    int update(vector<vector<Point>>*, vector<int>**, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*);
    string get_name();
    ~PAM();
};

template <class Point>
class MV_DTW : public Updater<Point> {
private:
    string name = "Mean Vector - DTW centroid Curve";
public:
    MV_DTW(int K) {this->K = K;};
    int update(vector<vector<Point>>*, vector<int>**, vector<pair<vector<Point>*, int>>*, DistanceDatabase<Point>*);
    string get_name();
    ~MV_DTW();
};

int mv_dtw_datatype(vector<vector<double>>*, vector<int>**, vector<pair<vector<double>*, int>>*);
int mv_dtw_datatype(vector<vector<double*>>*, vector<int>**, vector<pair<vector<double*>*, int>>*);