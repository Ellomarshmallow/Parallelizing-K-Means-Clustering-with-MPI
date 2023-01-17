#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <stddef.h>

#define MAX_LINE_LENGTH 18000
#define UNASSIGNED -1
#define NOISE -2
#define EPS  612
#define MIN_POINTS  3

//dbscan
//write
//read
//distance
//region
//expand
// split


/*smaller values of EPS will result in more clusters, while larger values of EPS will result in fewer clusters. 
Similarly, smaller values of MIN_PTS will result in more clusters, while larger values of MIN_PTS will result in 
fewer clusters. The optimal values of EPS and MIN_PTS will depend on the characteristics of the dataset and the desired 
number of clusters.*/

typedef struct {
    double x[6];
    int cluster;
} Point;



int read_csv(const char *file_name, Point *points, int *n_points) {
    FILE *file = fopen(file_name, "r");
    if (!file) {
        printf("Error: Could not open file.\n");
        return 1;
    }

    
    int i = 0;
/
    while (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf", &points[i].x[0], &points[i].x[1], &points[i].x[2], &points[i].x[3], &points[i].x[4], &points[i].x[5]) == 6) {
        points[i].cluster = UNASSIGNED;
        i++;
    }
    *n_points = i;

    fclose(file);
    return 0;
}

double distance(Point a, Point b) {
    double dx = a.x[0] - b.x[0];
    double dy = a.x[1] - b.x[1];
    double dz = a.x[2] - b.x[2];
    double dw = a.x[3] - b.x[3];
    double dv = a.x[4] - b.x[4];
    double du = a.x[5] - b.x[5];
    return dx * dx + dy * dy + dz * dz + dw * dw + dv * dv + du * du;
}

void region_query(Point point, Point *points, int n_points, int *neighbors, int *n_neighbors) {
    int i;
    for (i = 0; i < n_points; i++) {
        //printf("%d", n_points);
        if (distance(point, points[i]) <= EPS * EPS) {
            neighbors[(*n_neighbors)++] = i;
        }
    }
}

void expand_cluster(Point point, Point *points, int n_points, int *neighbors, int n_neighbors, int c) {
    int i;
    for (i = 0; i < n_neighbors; i++) {
        int j = neighbors[i];
        if (points[j].cluster == UNASSIGNED || points[j].cluster == NOISE) {
            points[j].cluster = c;
            int new_neighbors[n_points];
            int n_new_neighbors = 0;
            region_query(points[j], points, n_points, new_neighbors, &n_new_neighbors);
            if (n_new_neighbors >= MIN_POINTS) {
                expand_cluster(points[j], points, n_points, new_neighbors, n_new_neighbors, c);
            }
        }
    }
    }
    

int write_csv(const char *file_name, Point *points, int n_points) {
FILE *file = fopen(file_name, "w");
if (!file) {
printf("Error: Could not write file.\n");
return 1;
}
int i;
for (i = 0; i < n_points; i++) {
    fprintf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%d\n", points[i].x[0], points[i].x[1], points[i].x[2], points[i].x[3], points[i].x[4], points[i].x[5], points[i].cluster);
}

fclose(file);
return 0;
}

Point* split_dataset(const Point *original_points, int n_points, Point *local_points, int *local_n_points, int rank, int n_processes) {
    int i;
    int points_per_process = n_points / n_processes;
    int remaining_points = n_points % n_processes;
    int starting_index = rank * points_per_process + (rank < remaining_points ? rank : remaining_points);
    int ending_index = starting_index + points_per_process;
    if (rank < remaining_points) {
        ending_index++;
    }
    *local_n_points = 0;
    for (i = starting_index; i < ending_index; i++) {
        local_points[(*local_n_points)++] = original_points[i];
    }

    return local_points;
}



void dbscan(Point *points, int n_points) {
int c = 0;
int i;
for (i = 0; i < n_points; i++) {
if (points[i].cluster != UNASSIGNED) {
continue;
}
int neighbors[n_points];
int n_neighbors = 0;
region_query(points[i], points, n_points, neighbors, &n_neighbors);
if (n_neighbors < MIN_POINTS) {
points[i].cluster = NOISE;
} else {
c++;
points[i].cluster = c;
expand_cluster(points[i], points, n_points, neighbors, n_neighbors, c);
}
}
}



void gather_dataset(Point *local_points, int local_n_points, Point *clustered_points, int chunk_size, int rank, int size) {
    if (rank == 0) {
        memcpy(clustered_points, local_points, local_n_points * sizeof(Point));
        int i;
        for (i = 1; i < size; i++) {
            //MPI_Recv(clustered_points + i * chunk_size, chunk_size * sizeof(Point), MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    } else {
        MPI_Send(local_points, chunk_size * sizeof(Point), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
}






int main(int argc, char **argv) {
int rank, size;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

char *nnodes = argv[1];
char *ncpus = argv[2];
char *filename = argv[3];

Point *points = malloc(MAX_LINE_LENGTH * sizeof(Point));
int n_points;

read_csv(filename, points, &n_points);

Point *local_points = malloc(n_points * sizeof(Point));
int local_n_points;


// Broadcast the number of points to all ranks
MPI_Bcast(&n_points, 1, MPI_INT, 0, MPI_COMM_WORLD);

// Split the dataset into chunks for each rank
int chunk_size = n_points / size;

//split_dataset(points, n_points, local_points, &local_n_points, rank, size);
local_n_points = 0;
if (rank == size - 1) {
        local_n_points += n_points % size;
        
    }

// Define the MPI datatype for the Point struct
MPI_Datatype MPI_Point;
const int nitems = 6;
int blocklengths[6] = {1, 1, 1, 1, 1, 1};
MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
MPI_Aint offsets[6];


offsets[0] = offsetof(Point, x[1]);
offsets[1] = offsetof(Point, x[2]);
offsets[2] = offsetof(Point, x[3]);
offsets[3] = offsetof(Point, x[4]);
offsets[4] = offsetof(Point, x[5]);
offsets[5] = offsetof(Point, x[6]);

MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_Point);
MPI_Type_commit(&MPI_Point);


    /// measure the start time
double start_time = MPI_Wtime();
MPI_Scatter(points, chunk_size, MPI_Point, local_points, chunk_size, MPI_Point, 0, MPI_COMM_WORLD);

    // Each rank performs DBSCAN on its local dataset
dbscan(local_points, local_n_points);


    // Gather the clustered points from all ranks
Point *clustered_points = malloc(n_points * sizeof(Point));


//gather_dataset(local_points, local_n_points, clustered_points, chunk_size, rank, size);
MPI_Gather(local_points, local_n_points, MPI_Point, clustered_points, local_n_points, MPI_Point, 0, MPI_COMM_WORLD);
double end_time = MPI_Wtime();

    if (rank == 0) {
        
        int recv_count;
        MPI_Status status;
        MPI_Get_count(&status, MPI_BYTE, &recv_count);
        //printf("Rank %d received %d bytes in total\n", rank, recv_count);
        //printf("Cluster_points %Point \n", clustered_points);
        write_csv("/home/taylor.lucero/HPC4DS/output4.csv", clustered_points, n_points);
            

        // Calculate the total time taken
        double total_time = end_time - start_time;

        // Calculate the speedup
        //double speedup = (1.0 * n_points) / (total_time * size);
        

        //Calulate number of bytes 
        long int num_B = sizeof(Point) * n_points;
        


// Print the results
        //printf("%i", n_points); // returned nothing
       // printf("Total time taken: %f seconds\n", total_time);
        //printf("Start time: %f\n", start_time);
       //  printf("End time: %f\n", end_time);
       // printf("Speedup: %f\n", speedup);
        //printf("Number of nodes: %s\n", nnodes);
        //printf("Numbes of cpus: %s\n", ncpus);
        //printf("NNumber of points:%i\n", n_points);
        printf("%s,%s,%d,%ld,%lf\n", nnodes, ncpus, size, num_B, total_time);


//
    }
    MPI_Type_free(&MPI_Point);
    MPI_Finalize();
    return 0;
}

