#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// XXX I AM RUNNING THIS LOCALLY AS IT IS SHOULD NOT BE RUN ON THE LOGIN NODE AND I HAVEN't MADE A JOB SUBMIT SCRIPTXXX

typedef struct point
{
    double x, y; // data points
    int cluster; // cluster association
} point;

void save_data_sheet(point *pts, int k, int nptsincluster)
{
    int i;
    FILE *fptr = fopen("data.csv", "w");
    fprintf(fptr, "%c,%c,%c\n", "x", "y", "cluster");

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
    if (0)
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
        if (distance < min_distance) //XXX: what if there is equal distance between two different centroids?
        {
            min_distance = distance;
            cluster_assignment = i;
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

int main(void)
{

    int nptsincluster = 30;
    int k = 3;
    int i, l;
    int max_iterations = 1000;

    point *pts;
    point *init_centroids;
    point *current_centroids;

    current_centroids = calloc(k, sizeof(point));

    srand(1337);

    // 0. Generate random 2d data for now
    pts = generate_data(k, nptsincluster);

    // 1. Randomly choose k initial centroids
    init_centroids = initial_centroids(k, nptsincluster, pts);
    current_centroids = init_centroids;

    for (i = 0; i <= max_iterations; i++)
    {
        int moved = 0;
        // 2. Find the euclidean distance between all data points in our set with the k centroids.
        // 3. Assign cluster based on distance
        for (l = 0; l < k * nptsincluster; l++)
        {
            pts[l].cluster = assign_cluster(pts[l], current_centroids, k);
        }

        // 4. Calculate new centroids
        point *new_centroids;
        new_centroids = calloc(k, sizeof(point));
        new_centroids = calc_new_centroids(k, nptsincluster, pts);

        // 5. See if centroids have changed. If not, break loop.
        moved = compare_centroids(current_centroids, new_centroids, k);
        if (moved == 0)
        {
            break;
        }
    }

    if (1) // Store data sheet with cluster assignments
    {
        save_data_sheet(pts, k, nptsincluster);
    }

    return 0;
}