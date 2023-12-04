

//General includes
/*#include <stdio.h>      //I/O
#include <stdlib.h>
#include <time.h>       //for code timing purposes
#include <math.h>*/
#include <bits/stdc++.h>
#include <omp.h>


//Parameters;
#define VERTICES 1000            //number of vertices
#define DENSITY 100             //minimum number of edges per vertex. DO NOT SET TO >= VERTICES/2
#define MAX_WEIGHT 100000         //max edge length + 1
#define INF_DIST 1000000000     //"infinity" initial value of each node
#define IMPLEMENTATIONS 2       //number of Dijkstra implementations
#define THREADS 8               //number of OMP threads
#define RAND_SEED 1234          //random seed
typedef float data_t;           //data type
using namespace std;
#define CPG 3.611
#define GIG 1000000000

mutex mtx;

int main() {

    srand(RAND_SEED);

    //functions
    void setIntArrayValue(int* in_array, int array_size, int value);
    void setDataArrayValue(data_t* in_array, int array_size, data_t init_value);
    void initializeGraphZero(data_t* graph, int num_vertices);
    void constructGraphEdge(data_t* graph, int* edge_count, int num_vertices , vector<pair<pair<int , int> , data_t>> &edge_list);


    // Bellman-ford implementations
    void bellmanfordCPUSerial(data_t* graph, data_t* node_dist, int* parent_node, int* visited_node, int num_vertices, int v_start , vector<pair<pair<int , int> , data_t>> &edge_list);
    void bellmanfordCPUParallel(data_t* graph, data_t* node_dist, int* parent_node, int* visited_node, int num_vertices, int v_start , vector<pair<pair<int , int> , data_t>> &edge_list);

    //timing
    struct timespec diff(struct timespec start, struct timespec end);
    struct timespec start, end;                     //timespec
    struct timespec time_stamp[IMPLEMENTATIONS];


    //declare variables and allocate memory
    int graph_size      = VERTICES*VERTICES*sizeof(data_t);     //adjacency matrix representation of graph
    int int_array       = VERTICES*sizeof(int);                 //array of vertex IDs. Vertices have int IDs.
    int data_array      = VERTICES*sizeof(data_t);              //array of vertex distances (depends on type of data used)
    data_t* graph       = (data_t*)malloc(graph_size);                  //graph itself
    data_t* node_dist   = (data_t*)malloc(data_array);                  //distances from source indexed by node ID
    int* parent_node    = (int*)malloc(int_array);                      //previous nodes on SP indexed by node ID
    int* edge_count     = (int*)malloc(int_array);                      //number of edges per node indexed by node ID
    int* visited_node   = (int*)malloc(int_array);                      //pseudo-bool if node has been visited indexed by node ID
    int *pn_matrix      = (int*)malloc(IMPLEMENTATIONS*int_array);      //matrix of parent_node arrays (one per each implementation)
    data_t* dist_matrix = (data_t*)malloc(IMPLEMENTATIONS*data_array);
    vector <pair<pair<int , int> , data_t>> edge_list;

    //initialize arrays and graph
    setIntArrayValue(edge_count, VERTICES, 0);
    initializeGraphZero(graph, VERTICES);
    constructGraphEdge(graph, edge_count, VERTICES , edge_list);
    free(edge_count);



    int i;
    int origin = (rand() % VERTICES);               //starting vertex
    printf("Origin vertex: %d\n\n", origin);

    /*  SERIAL Bellman-ford  */
    int version = 0;
    printf("Running serial...");
    clock_gettime(CLOCK_REALTIME, &start);
    bellmanfordCPUSerial(graph, node_dist, parent_node, visited_node, VERTICES, origin , edge_list);
    clock_gettime(CLOCK_REALTIME, &end);
    time_stamp[version] = diff(start, end);               //record time
    for (i = 0; i < VERTICES; i++) {
        pn_matrix[version*VERTICES + i] = parent_node[i];
        dist_matrix[version*VERTICES + i] = node_dist[i];
    }
    printf("Done!\n");

    /*  PARALLEL (OPENMP) Bellman-ford  */
    version++;
    printf("Running OpenMP...");
    clock_gettime(CLOCK_REALTIME, &start);
    bellmanfordCPUParallel(graph, node_dist, parent_node, visited_node, VERTICES, origin , edge_list);
    clock_gettime(CLOCK_REALTIME, &end);
    time_stamp[version] = diff(start, end);               //record time
    for (i = 0; i < VERTICES; i++) {
        pn_matrix[version*VERTICES + i] = parent_node[i];
        dist_matrix[version*VERTICES + i] = node_dist[i];
    }
    printf("Done!\n");

    printf("\nVertices: %d", VERTICES);
    printf("\nDensity: %d", DENSITY);
    printf("\nMax Weight: %d", MAX_WEIGHT);
    printf("\n\nTime (seconds):\nSerial,OpenMP\n");
    for (i = 0; i < IMPLEMENTATIONS; i++) {
        /*printf("%ld,", (long long)( (double)(CPG)*(double)
            (GIG * time_stamp[i].tv_sec + time_stamp[i].tv_nsec)));*/
        printf("%.8f,", (double)( (double)
            (time_stamp[i].tv_sec + time_stamp[i].tv_nsec / (GIG * 1.0))));
    }


}



/*  Construct graph with randomized edges.  */
void constructGraphEdge(data_t* graph, int* edge_count, int num_vertices , vector<pair<pair<int , int> , data_t>> &edge_list) {
    int i;                  //iterator
    int rand_vertex;        //random previous vertex
    int curr_num_edges;     //current number of edges per vertex
    int num_edges;          //edges per vertex
    data_t edge, weight;    //edge chance and weight

    //initialize a connected graph
    for (i = 1; i < num_vertices; i++) {
        rand_vertex = (rand() % i);                     //select a random previous vertex to create a connected graph
        weight = (rand() % MAX_WEIGHT) + 1;             //random (non-zero) weight
        graph[rand_vertex*num_vertices + i] = weight;   //set edge weights
        graph[i*num_vertices + rand_vertex] = weight;
        edge_list.push_back({{i , rand_vertex} , weight});
        edge_list.push_back({{rand_vertex , i} , weight});
        edge_count[i] += 1;                             //increment edge counts for each vertex
        edge_count[rand_vertex] += 1;
    }

    //add additional edges until DENSITY reached for all vertices
    for (i = 0; i < num_vertices; i++) {    //for each vertex
        curr_num_edges = edge_count[i];         //current number of edges (degree) of vertex
        while (curr_num_edges < DENSITY) {      //add additional edges if number of edges < DENSITY
            rand_vertex = (rand() % num_vertices);  //choose any random vertex
            weight = (rand() % MAX_WEIGHT) + 1;     //choose a random (non-zero) weight
            if ((rand_vertex != i) && (graph[i*num_vertices + rand_vertex] == 0)) { //add edge if not trying to connect to itself and no edge currently exists
                graph[i*num_vertices + rand_vertex] = weight;
                graph[rand_vertex*num_vertices + i] = weight;
                edge_list.push_back({{i , rand_vertex} , weight});
                edge_list.push_back({{rand_vertex , i} , weight});
                edge_count[i] += 1;
                curr_num_edges++;               //one additional edge constructed
            }
        }
    }
}

/*  Initialize elements of a 1D int array with an initial value   */
void setIntArrayValue(int* in_array, int array_size, int init_value) {
    int i;
    for (i = 0; i < array_size; i++) {
        in_array[i] = init_value;
    }
}

/*  Initialize elements of a 1D data_t array with an initial value   */
void setDataArrayValue(data_t* in_array, int array_size, data_t init_value) {
    int i;
    for (i = 0; i < array_size; i++) {
        in_array[i] = init_value;
    }
}

/*  Construct graph with no edges or weights     */
void initializeGraphZero(data_t* graph, int num_vertices) {
    int i, j;

    for (i = 0; i < num_vertices; i++) {
        for (j = 0; j < num_vertices; j++) {           //weight of all edges initialized to 0
            graph[i*num_vertices + j] = (data_t)0;
        }
    }
}

/*  Print int array elements    */
void checkArray(int* a, int length) {
    int i;
    printf("Proof: ");
    for (i = 0; i < length; i++) {
        printf("%d, ", a[i]);
    }
    printf("\n");
}

/*   Difference in two timespec objects   */
struct timespec diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec - start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    }
    else {
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    return temp;
}


/****************Bellman-ford'S ALGORITHM IMPLEMENTATIONS****************/


void bellmanfordCPUSerial(data_t* graph, data_t* node_dist, int* parent_node, int* visited_node, int num_vertices, int v_start , vector<pair<pair<int , int> , data_t>> &edge_list) {

    //functions
    void setIntArrayValue(int* in_array, int array_size, int init_value);
    void setDataArrayValue(data_t* in_array, int array_size, data_t init_value);
    int closestNode(data_t* node_dist, int* visited_node, int num_vertices);
    void checkArray(int* a, int length);

    //reset/clear data from previous runs
    setDataArrayValue(node_dist, VERTICES, INF_DIST);     //all node distances are infinity
    node_dist[v_start] = 0;                     //start distance is zero; ensures it will be first pulled out

    for (int i = 0; i < num_vertices; i++) {

        for(int j = 0 ; j < edge_list.size() ; j++)
        {
            int u = edge_list[j].first.first;
            int v = edge_list[j].first.second;
            data_t weight = edge_list[j].second;

            if(node_dist[u] + weight < node_dist[v]){
                node_dist[v] = node_dist[u] + weight;
            }
        }
    }
}

void bellmanfordCPUParallel(data_t* graph, data_t* node_dist, int* parent_node, int* visited_node, int num_vertices, int v_start , vector<pair<pair<int , int> , data_t>> &edge_list) {

    //functions
    void setIntArrayValue(int* in_array, int array_size, int init_value);
    void setDataArrayValue(data_t* in_array, int array_size, data_t init_value);
    int closestNode(data_t* node_dist, int* visited_node, int num_vertices);
    void checkArray(int* a, int length);

    //reset/clear data from previous runs
    setDataArrayValue(node_dist, VERTICES, INF_DIST);     //all node distances are infinity
    node_dist[v_start] = 0;                     //start distance is zero; ensures it will be first pulled out

    for (int i = 0; i < num_vertices; i++) {
        omp_set_num_threads(THREADS);
        #pragma omp parallel for
        for(int j = 0 ; j < edge_list.size() ; j++)
        {
            int u = edge_list[j].first.first;
            int v = edge_list[j].first.second;
            data_t weight = edge_list[j].second;

            if(node_dist[u] + weight < node_dist[v]){
                mtx.lock();
                node_dist[v] = min(node_dist[u] + weight , node_dist[v]);
                mtx.unlock();
            }
        }
    }
}
