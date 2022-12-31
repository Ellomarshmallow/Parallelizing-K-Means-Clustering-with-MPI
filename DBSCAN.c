#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_LINE_LEN 1024
#define EPS 20
#define MIN_PTS 2

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
void dbscan(Point *points, int num_points, int *cluster_id);
double euclidean_distance(Point p, Point q);
int find_neighbors(Point *points, int num_points, Point p, int *neighbors);

int main(int argc, char *argv[]) {
  // Check command line arguments
  if (argc != 2) {
    fprintf(stderr, "Usage: %s filename\n", argv[0]);
    return 1;
  }

  // Open the input file
  FILE *fp = fopen(argv[1], "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: Can't open file %s\n", argv[1]);
    return 1;
  }

  // Read the points from the input file
  Point points[MAX_LINE_LEN];
  int num_points = 0;
  char line[MAX_LINE_LEN];
  while (fgets(line, MAX_LINE_LEN, fp) != NULL) {
    if (num_points == MAX_LINE_LEN) {
      fprintf(stderr, "Error: Too many points\n");
      return 1;
    }
    sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &points[num_points].x[0], &points[num_points].x[1], &points[num_points].x[2], &points[num_points].x[3], &points[num_points].x[4], &points[num_points].x[5], &points[num_points].x[6]);
    num_points++;
  }

  // Close the input file
  fclose(fp);

  // Run DBSCAN on the points
  int cluster_id[num_points];
  dbscan(points, num_points, cluster_id);

  // Print the results
  for (int i = 0; i < num_points; i++) {
    printf("Point %d belongs to cluster %d\n", i, cluster_id[i]);
  }

  return 0;
}

// Function to run DBSCAN on the points
void dbscan(Point *points, int num_points, int *cluster_id) {
  // Initialize the cluster IDs
  for (int i = 0; i < num_points; i++) {
    cluster_id[i] = -1;
  }

  // Initialize the cluster index
  int c = 0;

  // Iterate through each point
  for (int i = 0; i < num_points; i++) {
    // Skip already processed points
    if (cluster_id[i] != -1) {
      continue;
    }

    // Find the neighbors of point i
    int neighbors[num_points];
    int num_neighbors = find_neighbors(points, num_points, points[i], neighbors);
    printf("the number of neighbors is %i", num_neighbors);
   // printf("The minimum number of points is %i", MIN_PTS);


    // If the point has fewer than MIN_PTS neighbors, mark it as noise
    if (num_neighbors < MIN_PTS) { 
        // for the cluster id to produce only -2 then it is clear to say that the number neighbors is less than min point at every step
        // But why is this? 
        // there is one true reason; num < min, however why ? and how? what is min point and what is num neighbors.

        // Start of Min_PTS  : Input as 2 
        //  at the point before the conditional statement we find that the point is still 2 
        // num neighbors is set to 1... why is this 1, did i set it to one? 

        // Num neighbors taks (points, num_points, points[i], neighbros ) and gives the num of neighbors
        // the num of neigbors is initialized as 0
        // the for loop then states, starts at 0, the with condition that i stays less than the num_points, it increments by 1  num_points is 1000
        // conditional if statement: if euclidean distance (p , points[i]) <= Epsilon ; epsilon in this case is .99
        // t
        // THE euclidean distance
        // init the sum as 0
        // for int i starting at 0 and as long as i is less than 7 and increments by 1. 
        //10,1,32,38,47,32,20
        //42,31,20,1,16,45,9

        /*
        // Function to compute the Euclidean distance between two points
            double euclidean_distance(Point p, Point q) {
            double sum = 0;
                for (int i = 0; i < 7; i++) {
                sum += (p.x[i] - q.x[i]) * (p.x[i] - q.x[i]);
  }
  return sqrt(sum);
}
        
        
         */



        /* 
        // Function to find the neighbors of a point within a certain radius
        int find_neighbors(Point *points, int num_points, Point p, int *neighbors) {
        int num_neighbors = 0;
        for (int i = 0; i < num_points; i++) {
            if (euclidean_distance(p, points[i]) <= EPS) {
            neighbors[num_neighbors++] = i;
    }
  }
  return num_neighbors;
}
        */
      cluster_id[i] = -2;
    }
    // Otherwise, create a new cluster with point i as its seed
    else {
      Cluster cluster;
      cluster.size = 0;
      cluster.points = malloc(num_points * sizeof(int));
      cluster.points[cluster.size++] = i;
      cluster_id[i] = c; // this takes c which is  the cluster at a certain index 

      // Iterate through the neighbors
      for (int j = 0; j < num_neighbors; j++) {
        // Skip already processed points
        if (cluster_id[neighbors[j]] != -1 ) {
          continue;
        }

        // Mark the point as belonging to cluster c
        cluster_id[neighbors[j]] = c;
        cluster.points[cluster.size++] = neighbors[j];

        // Find the neighbors of the point
        int neighbors2[num_points];
        int num_neighbors2 = find_neighbors(points, num_points, points[neighbors[j]], neighbors2);

        // If the point has at least MIN_PTS neighbors, add them to the cluster
        if (num_neighbors2 >= MIN_PTS) {
          for (int k = 0; k < num_neighbors2; k++) {
            // Skip already processed points
            if (cluster_id[neighbors2[k]] != -1) {
              continue;
            }

            // Add the point to the cluster
            cluster.points[cluster.size++] = neighbors2[k];
            cluster_id[neighbors2[k]] = c;
          }
        }
      }

      // Increment the cluster index
      c++;
    }
  }
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
    if (euclidean_distance(p, points[i]) <= EPS) {
      neighbors[num_neighbors++] = i;
      printf("neighbors are %i \n", num_neighbors);
      
    }
  }
  return num_neighbors;
}