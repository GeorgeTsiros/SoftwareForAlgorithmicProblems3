#include "HashTable.h"
#ifndef DATABASE
#define DATABASE

#include "Database.h"

#endif

using namespace std;

template <class Point>
class LSH {
private:
    /* dataset */
    vector<vector<Point>>* dataset;
    /* Vector containing shifts of size(l,k,d) */
    vector<vector<vector<double>>> s;
    /* Amplified hash for dataset*/
    vector<vector<int>> data_amplified_g;
    /* Vector containing projections of data */
    vector<vector<int>> a_projects;
    /* L Hash Tables */
    HashTable <Point> **MyHashTable;
    /* Size of hash Table */
    int TableSize;
    /* array for hash computations */
    int * power;
    /* parameters LSH */
    int k;
    int L;
    double w;
    int d;
    /* end of params */
public:
    LSH(int, int, double);
    void fit(vector<vector<Point>>*);
    void evaluate(vector<vector<Point>>*, double, vector<vector<int>>*, Point**, double**, int**);
    void evaluate_clusters(vector<vector<Point>>*, vector<int>*, Point**, int**, int*);
    ~LSH();
};