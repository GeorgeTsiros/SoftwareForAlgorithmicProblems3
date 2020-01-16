#include "Initializers.h"
#include "Assigners.h"
#include "Updaters.h"

using namespace std;

template <class Point>
class Cluster {
private:
    Initializer<Point>* initializer;
    Assigner<Point>* assigner;
    Updater<Point>* updater;
    int K;
    int Grids;
    int L;
    int k;
    double w;
    vector<pair<vector<Point>*, int>> centroids;
    vector<int>** clusters;
public:
    Cluster(int*, string, string, string, double);
    void fit(vector<vector<Point>>*, DistanceDatabase<Point>*);
    vector<double> silhouette(vector<vector<Point>>*, DistanceDatabase<Point>*);
    double average_distance(int, vector<int>*, DistanceDatabase<Point>*);
    int find_closest_centroid(pair<vector<Point>*, int>, DistanceDatabase<Point>*);
    vector<pair<vector<Point>*, int>>* get_centroids() {return &this->centroids;}
    vector<int>** get_clusters() {return this->clusters;}
    ~Cluster();
};

void clear_centroids(vector<pair<vector<double>*, int>>*);
void clear_centroids(vector<pair<vector<double*>*, int>>*);