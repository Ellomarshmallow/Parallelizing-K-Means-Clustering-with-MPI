#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h> 

/*smaller values of EPS will result in more clusters, while larger values of EPS will result in fewer clusters. 
Similarly, smaller values of MIN_PTS will result in more clusters, while larger values of MIN_PTS will result in 
fewer clusters. The optimal values of EPS and MIN_PTS will depend on the characteristics of the dataset and the desired 
number of clusters.*/

#define MAX_LINE_LEN 1024
#define EPS 50
#define MIN_PTS 14

// Data structure to store a point in the 7-dimensional space
typedef struct {
  double x[7];
} Point;

// Data structure to store a cluster
typedef struct {
  int size;
  int *points;
} Cluster;

// Function prototypes
void dbscan(int rank, int size, Point *points, int num_points, int *cluster_id);
double euclidean_distance(Point p, Point q);
int find_neighbors(Point *points, int num_points, Point p, int *neighbors);
void expand_cluster(Point *points, int num_points, int *cluster_id, int p, int c, int *neighbors, int num_neighbors);
void write_output(int *cluster_id, int num_points, const char* filename);
void split_dataset(int rank, int size, Point *points, int num_points, Point **local_points, int *local_num_points);







int main(int argc, char *argv[]) {

  // Initialize the MPI
  MPI_Init(&argc, &argv);

  // Rank and Size of the MPI Cluster
  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD , &rank);
  MPI_Comm_size( MPI_COMM_WORLD , &size);

  // Check command line arguments
  if (argc != 2) {
    if (rank == 0) {
      fprintf(stderr, "Usage: %s filename\n", argv[0]);
    }
    MPI_Finalize();
    return 1;
  }

  // Open the input file
  FILE *fp = fopen(argv[1], "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: Can't open file %s\n", argv[1]);
    MPI_Finalize();
    return 1;
  }

  // Read the points from the input file
  Point points[MAX_LINE_LEN];
  int num_points = 0;
  char line[MAX_LINE_LEN];
  while (fgets(line, MAX_LINE_LEN, fp) != NULL) {
    // Parse??
    char *token =strtok(line, ",");
    int i = 0;
    while (token != NULL) {
      points[num_points].x[i++]=atof(token);
      token=strtok(NULL, ",");
    }
    num_points++;
  }
/*
    if (num_points == MAX_LINE_LEN) {
      fprintf(stderr, "Error: Too many points\n");
      return 1;
    }
    sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &points[num_points].x[0], &points[num_points].x[1], &points[num_points].x[2], &points[num_points].x[3], &points[num_points].x[4], &points[num_points].x[5], &points[num_points].x[6]);
    num_points++;
  }*/

  // Close the input file
fclose(fp);

  // Allocate memory
  int cluster_id[num_points];
  for (int i = 0; i < num_points; i++){
    cluster_id[i] = -1;
  }




dbscan(rank, size, points, num_points, cluster_id);
    //write_output(cluster_id, num_points, "output.csv"); //this just outputs all the points in the general cluster 

    //

   // write_output(cluster_id, local_num_points, "output2.csv"); // this only provides output of the local points
    

  MPI_Finalize();

  return 0;
}



/// Base dbscan function 
void dbscan(int rank, int size, Point *points, int num_points, int *cluster_id) {
// Allocate memory for the local dataset
Point *local_points;
int local_num_points;
split_dataset(rank, size, points, num_points, &local_points, &local_num_points);
// Broadcast the eps and min_pts to all the processes

double eps = EPS;
int min_pts = MIN_PTS;
MPI_Bcast(&eps, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&min_pts, 2, MPI_INT, 0, MPI_COMM_WORLD);

// Perform the DBSCAN algorithm on the local dataset
// Initialize cluster ID
for (int i = 0; i < local_num_points; i++) {
cluster_id[i] = -1;
}

int c = -1;
for (int i = 0; i < local_num_points; i++) {
if (cluster_id[i] != -1) {
continue;
}
int neighbors[local_num_points];
int num_neighbors = find_neighbors(local_points, local_num_points, local_points[i], neighbors);
if (num_neighbors < MIN_PTS) {
// Noise
cluster_id[i] = -1;
} else {
// Expand the cluster
c++;
expand_cluster(local_points, local_num_points, cluster_id, i, c, neighbors, num_neighbors);
}
}
MPI_Barrier(MPI_COMM_WORLD);

// Gather the cluster ID of each process
int *cluster_id_global = NULL;
if (rank == 0) {
    cluster_id_global = (int *)malloc(num_points * sizeof(int));
}
int recvcounts[size], displs[size];
for (int i = 0; i < size; i++) {
    recvcounts[i] = num_points / size;
if (i < num_points % size) {
    recvcounts[i]++;
}
if (i == 0) {
    displs[i] = 0;
} else {
    displs[i] = displs[i - 1] + recvcounts[i - 1];
}
}

MPI_Gatherv(cluster_id, local_num_points, MPI_INT, cluster_id_global, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

// Write the output to a file
if (rank == 0) {
write_output(cluster_id_global, num_points, "output6.csv");
free(cluster_id_global);
}

// Free the memory
free(local_points);
}



// Function to compute the Euclidean distance between two points
double euclidean_distance(Point p, Point q) {
  double sum = 0;
  for (int i = 0; i < 7; i++) {
    sum += (p.x[i] - q.x[i]) * (p.x[i] - q.x[i]);
    //printf("The sum is %f \n", sum);
  }
  return sqrt(sum);
}

// Function to find the neighbors of a point within a certain radius
int find_neighbors(Point *points, int num_points, Point p, int *neighbors) {
  int num_neighbors = 0;
  for (int i = 0; i < num_points; i++) {
    if (euclidean_distance(points[i], p) <= EPS) {
      neighbors[num_neighbors++] = i;
    }
  }
  return num_neighbors;
}


void expand_cluster(Point *points, int num_points, int *cluster_id, int p, int c, int *neighbors, int num_neighbors) {
  cluster_id[p] = c;
  for (int i = 0; i < num_neighbors; i++) {
    int neighbor = neighbors[i];
    if (cluster_id[neighbor] == -1) {
      cluster_id[neighbor] = c;
      int new_neighbors[num_points];
      int num_new_neighbors = find_neighbors(points, num_points, points[neighbor], new_neighbors);
      if (num_new_neighbors >= MIN_PTS) {
        expand_cluster(points, num_points, cluster_id, neighbor, c, new_neighbors, num_new_neighbors);
      }
    }
  }
}

void write_output(int *cluster_id, int num_points, const char* filename) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: Can't open output file\n");
        return;
    }
    for (int i = 0; i < num_points; i++) {
        fprintf(fp, "Point %d belongs to cluster %d\n", i, cluster_id[i]);
    }
    fclose(fp);
}

void split_dataset(int rank, int size, Point *points, int num_points, Point **local_points, int *local_num_points) {
    // Calculate the number of points assigned to each process
    int points_per_process = num_points / size;
    int extra_points = num_points % size;

    // Calculate the starting and ending indices for the local dataset
    int start_index = rank * points_per_process + (rank < extra_points ? rank : extra_points);
    int end_index = start_index + points_per_process + (rank < extra_points ? 1 : 0);
    if(rank == size-1) {
        end_index = num_points;
    }

    // Allocate memory for the local dataset
    *local_num_points = end_index - start_index;
    *local_points = (Point *)malloc(*local_num_points * sizeof(Point));

    // Copy the local dataset from the global dataset
    memcpy(*local_points, &points[start_index], *local_num_points * sizeof(Point));
}



