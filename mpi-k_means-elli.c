#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <unistd.h>

typedef struct point
{
    double x, y; // data points
    int cluster; // cluster association
} point;

void save_data_sheet(point *pts, int k, int nptsincluster)
{
    int i;
    FILE *fptr = fopen("data.csv", "w");
    fprintf(fptr, "%s,%s,%s\n", "x", "y", "cluster");

    for (i = 0; i < k * nptsincluster; i++)
    {
        fprintf(fptr, "%lf,%lf,%d\n", pts[i].x, pts[i].y, pts[i].cluster);
    }
    fclose(fptr);
}

point *generate_data(int k, int nptsincluster) // XXX: Not ideal data generation for clustering
{
    int i;
    int j;
    int spread = 3;

    point *pts;

    pts = calloc(k * nptsincluster, sizeof(point));

    for (i = 0; i < k; i++)
    {

        for (j = 0; j < nptsincluster; j++)
        {
            double u1, u2, z1, z2;
            u1 = 1.0 * random() / RAND_MAX;
            u2 = 1.0 * random() / RAND_MAX;
            z1 = spread * i + sqrt(-2 * log2(u1)) * cos(2 * M_PI * u2);
            z2 = spread * i + sqrt(-2 * log2(u1)) * sin(2 * M_PI * u2);
            int n = i * nptsincluster + j;

            pts[n].x = z1;
            pts[n].y = z2;

            // printf("%lf %lf\n", pts[n].x, pts[n].y);
        }
    }
    if (1)
    {
        save_data_sheet(pts, k, nptsincluster);
    }
    return pts;
}

point *initial_centroids(int k, int nptsincluster, point *pts)
{
    int r;
    int j;
    point *init_centroids;

    init_centroids = calloc(k, sizeof(point));
    // printf("Random initial centroids:\n", r);

    for (j = 0; j < k; j++)
    {
        r = rand() % (k * nptsincluster); // XXX: Doesn't ensure exclusivity when generating the numbers
        init_centroids[j] = pts[r];
        // printf("Random Num %d\n", r);
        // printf("Init centroid[%d]\t%g\t%g\n", j, init_centroids[j].x, init_centroids[j].y);
    }
    return init_centroids;
}

point *calc_new_centroids(int k, int nptsincluster, point *pts)
{
    int i, j;
    point *new_centroids;
    new_centroids = calloc(k, sizeof(point));
    float new_x, new_y;

    for (i = 0; i < k; i++)
    {
        point sum;
        int count = 0;
        sum.x = sum.y = 0.0;

        for (j = 0; j < k * nptsincluster; j++)
        {
            if (pts[j].cluster == i)
            {
                sum.x += pts[j].x;
                sum.y += pts[j].y;
                count += 1;
            }
        }

        new_x = sum.x / count;
        new_y = sum.y / count;

        new_centroids[i].x = new_x;
        new_centroids[i].y = new_y;
    }
    return new_centroids;
}

int assign_cluster(point pt, point *current_centroids, int k)
{
    int i;
    float distance;
    float min_distance = 1000000; // XXX: Very hacky
    int cluster_assignment;

    for (i = 0; i < k; i++)
    {
        distance = sqrt(pow(pt.x - current_centroids[i].x, 2) + pow(pt.y - current_centroids[i].y, 2));
        if (distance < min_distance) // XXX: what if there is equal distance between two different centroids?
        {
            min_distance = distance;
            cluster_assignment = i;
            current_centroids[i].cluster = i;
        }
    }
    return cluster_assignment;
}

int compare_centroids(point *current_centroids, point *new_centroid, int k)
{
    int diff = 0;
    int i;

    for (i = 0; i < k; i++)
    {
        diff = current_centroids[i].x - new_centroid[i].x;
        if (diff != 0)
            return diff;

        diff = current_centroids[i].y - new_centroid[i].y;
        if (diff != 0)
            return diff;
    }
    return diff;
}

int main(int argc, char **argv)
{

    // For debugging
    // {
    //     int i=0;
    //     while(0==i){
    //         sleep(5);
    //     }
    // }

    MPI_Init(&argc, &argv);

    // Get the number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // XXX: Split variable declarations into for all processes and only root process?

    int nptsincluster = 30;
    int k = 3;
    const int root = 0;
    int elements_per_proc = ceil((k * nptsincluster) / size);
    double start_time, stop_time, elapsed_time;

    point *pts;
    point *init_centroids;
    point *current_centroids;
    point *final_centroids;
    point *sub_pts;

    current_centroids = calloc(k, sizeof(point));
    final_centroids = calloc(k * size, sizeof(point));
    sub_pts = calloc(elements_per_proc, sizeof(point));

    // Define own MPI Datatype for the struct
    MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Datatype custom_type;
    int array_of_blocklengths[] = {1, 1, 1};

    MPI_Aint array_of_displacements[] = {offsetof(point, x),
                                         offsetof(point, y),
                                         offsetof(point, cluster)};

    MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements, types, &custom_type);
    MPI_Type_commit(&custom_type);

    srand(1337);

    if (rank == 0) // Data needs to sit on root process to be scattered later on
    {
        // 0. Generate random 2d data for now
        pts = generate_data(k, nptsincluster);

        // Benchmark with data in memory already
        start_time = MPI_Wtime();

        // 1. Randomly choose k initial centroids
        init_centroids = initial_centroids(k, nptsincluster, pts);
        current_centroids = init_centroids;
    }

    // Step 2. Send data (split by n nodes) and inital k centroids to the different nodes from node 0
    MPI_Scatter(pts, elements_per_proc, custom_type, sub_pts, elements_per_proc, custom_type, root, MPI_COMM_WORLD);
    // XXX: think about the effects elements_per_proc has when it is not evenly distributed
    MPI_Bcast(current_centroids, k, custom_type, root, MPI_COMM_WORLD);

    // // In each node: --> removed to have process 0 work too XXX: think about if this is correct
    // if (rank != 0)
    // {
    int i, l;
    int max_iterations = 1000;
    point *new_centroids;
    new_centroids = calloc(k, sizeof(point));

    for (i = 0; i <= max_iterations; i++)
    {
        int moved = 0;
        // 3. Find the euclidean distance between all data points in our set with the k centroids.
        // 3b. Assign cluster based on distance
        for (l = 0; l < elements_per_proc; l++) // XXX: this might lead to one centroid being 0 because in non even elements_per_proc distributions?
        {
            sub_pts[l].cluster = assign_cluster(sub_pts[l], current_centroids, k);
            // printf("%lf %lf %d\n", sub_pts[l].x, sub_pts[l].y, sub_pts[l].cluster); //making sure the data looks good
        }

        // 4. Calculate new centroids
        new_centroids = calc_new_centroids(k, elements_per_proc, sub_pts);

        // 5. See if centroids have changed. If not, break loop.
        moved = compare_centroids(current_centroids, new_centroids, k);

        if (moved == 0)
        {
            // printf("Final centroids from rank %d:\n", rank);
            // int i;
            // for (i = 0; i < k; i++)
            // {
            //     printf("%lf %lf %d\n", current_centroids[i].x, current_centroids[i].y, current_centroids[i].cluster); // XXX: Seems like one centroid is incorrect? (0.0 0.0)
            // }
            break;
        }
        else
        {
            current_centroids = new_centroids;
        }
    }
    // }

    // XXX: Think about adding a barrier?

    MPI_Gather(current_centroids, k, custom_type, final_centroids, k, custom_type, root, MPI_COMM_WORLD);
    // XXX: Think about adding a barrier?

    if (rank == root)
    {
        // printf("Received centroids in rank %d:\n", rank);
        // int i;
        // for (i = 0; i < k * size; i++)
        // {
        //     printf("%lf %lf %d\n", final_centroids[i].x, final_centroids[i].y, final_centroids[i].cluster); // making sure the received centroids look good
        // }

        // Step 6 calc final centroids as average of all received ones TODO: this is wrong
        int i;
        for (i = 0; i < k; i++)
        {
            int x = 0;
            int y = 0;
            int j;
            for (j = 0; j < k * size; j++)
            {
                if (final_centroids[i].cluster == i)
                {
                    x += final_centroids[j].x;
                    y += final_centroids[j].y;
                }
            }
            final_centroids[i].x = x / size;
            final_centroids[i].y = y / size;
        }

        // int q;
        // for (q = 0; q < k; q++)
        // {
        //     printf("%lf %lf %d\n", final_centroids[q].x, final_centroids[q].y, final_centroids[q].cluster); // making sure the received centroids look good
        // }

        // Step 7. Assign final clusters to all points
        for (l = 0; l < k * nptsincluster; l++)
        {
            pts[l].cluster = assign_cluster(pts[l], final_centroids, k);
        }
    }
    // Think about adding barrier?

    long int num_B = k * nptsincluster * sizeof(point); // data size XXX: but not in bytes??
    long int B_in_GB = 1 << 30;
    double num_GB = (double)num_B / (double)B_in_GB;

    if ((rank == root) && 1)
    {
        stop_time = MPI_Wtime();
        elapsed_time = stop_time - start_time;

        // Store data sheet with cluster assignments
        save_data_sheet(pts, k, nptsincluster);
        printf("Transfer size (B): %10li, Transfer Time (s): %15.9f, Bandwidth (GB/s): %15.9f\n", num_B, elapsed_time, num_GB / elapsed_time);
    }

    // Free memory?

    MPI_Finalize();
    return 0;
}