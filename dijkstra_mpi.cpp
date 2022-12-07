/*----------------------------------------------------
 * File:    mpi_Dijkstra.c
 *
 * Purpose: implement Dijkstra's shortest path algorithm
 *          for a weighted directed graph using mpi
 *
 * Compile: mpicc -g -Wall -o mpi_Dijkstra mpi_Dijkstra.c
 * Run:     mpiexec -n <number of processes> mpi_Dijkstra
 *
 * Input:   n: the number of vertices
 *          mat: the matrix where mat[i][j] is the length
 *               from vertex i to j
 *
 * Output:  length of the shortest path from vertex 0 to vertex v
 *          Shortest path to each vertex v from vertex 0
 *
 * Algorithm: the matrix mat is partioned by columns so that each
 *            process gets n / p columns. In each iteration each
 *            process finds its local vertex with the shortest distance
 *            from the source vertex 0. A global minimum vertex u of the found
 *            shortest distances is computed and then each process updates
 *            its local distance array if there's a shorter path that goes through u
 *
 * Note:    1. This program assumes n is evenly divisible by
 *             p (number of processes)
 *          2. Edge weights should be nonnegative
 *          3. If there is no edge between any two vertices the weight is the constant
 *             INFINITY
 *          4. The cost of traveling to itself is 0
 *          5. The adjacency matrix is stored as an 1-dimensional array and subscripts
 *             are computed using A[n * i + j] to get A[i][j] in the 2-dimensional case
 *
 * Author: Henrik Lehmann
 *-----------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <iostream>


using namespace std;


#define ROWMJR(R,C,NR,NC) (R*NC+C)
#define COLMJR(R,C,NR,NC) (C*NR+R)
/* define access directions for matrices */
#define a(R,C) a[ROWMJR(R,C,ln,n)]
#define b(R,C) b[ROWMJR(R,C,nn,n)]

int Read_n(int my_rank, MPI_Comm comm, int n);
MPI_Datatype Build_blk_col_type(int n, int loc_n);
void Read_matrix(float loc_mat[], int n, int loc_n,
                 MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm,const char * const filename);
void Dijkstra(float loc_mat[], float loc_dist[], float loc_pred[], int loc_n, int n,
              MPI_Comm comm);
void Dijkstra_Init(float loc_mat[], float loc_pred[], float loc_dist[], float loc_known[],int my_rank, int loc_n);
int Find_min_dist(float loc_dist[], float loc_known[], int loc_n);
void Print_matrix(float global_mat[], int rows, int cols);
void Print_dists(float global_dist[], int n);
void Print_paths(float global_pred[], int n);


static void
load(
  const char * const filename,
  int * const np,
  float ** const ap
)
{
  int i, j, n, ret;
  FILE * fp=NULL;
  float *a, *mat;

  /* open the file */
  fp = fopen(filename, "r");
  assert(fp);

  /* get the number of nodes in the graph */
  ret = fscanf(fp, "%d", &n);
  assert(1 == ret);

  /* allocate memory for local values */
  a = (float*) malloc(n*n*sizeof(*a));
  assert(a);

  /* read in roots local values */
  mat = (float*) malloc(n*n*sizeof(float));
  
  for (i=0; i<n; ++i) {
    for (j=0; j<n; ++j) {
      ret = fscanf(fp, "%f",&a(i,j));
      mat[i*n+j]=a[i * n + j];
      //cout<<a[i * n + j]<<" ";
      assert(1 == ret);
    }
  }

  /* close file */
  ret = fclose(fp);
  assert(!ret);
  Print_matrix(mat,n,n);
  /* record output values */
  *np = n;
  *ap = a;
}


static void
print_numbers(
  const char * const filename,
  const int n,
  const float * const numbers)
{
  int i;
  FILE * fout;

  /* open file */
  if(NULL == (fout = fopen(filename, "w"))) {
    fprintf(stderr, "error opening '%s'\n", filename);
    abort();
  }

  /* write numbers to fout */
  for(i=0; i<n; ++i) {
    fprintf(fout, "%10.4f\n", numbers[i]);
  }

  fclose(fout);
}


int main(int argc, char **argv) {
    float *loc_mat, *loc_dist, *loc_pred, *global_dist = NULL, *global_pred = NULL;
    int my_rank, p, loc_n, n;
    MPI_Comm comm;
    MPI_Datatype blk_col_mpi_t;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &p);
    n = Read_n(my_rank, comm, atoi(argv[2]));
    loc_n = n / p;
    loc_mat = (float*)malloc(n * loc_n * sizeof(float));
    loc_dist = (float*)malloc(loc_n * sizeof(float));
    loc_pred = (float*)malloc(loc_n * sizeof(float));
    blk_col_mpi_t = Build_blk_col_type(n, loc_n);

    if (my_rank == 0) {
        global_dist = (float*)malloc(n * sizeof(float));
        global_pred = (float*)malloc(n * sizeof(float));
    }
    Read_matrix(loc_mat, n, loc_n, blk_col_mpi_t, my_rank, comm,argv[1]);
    Dijkstra(loc_mat, loc_dist, loc_pred, loc_n, n, comm);

    /* Gather the results from Dijkstra */
    MPI_Gather(loc_dist, loc_n, MPI_INT, global_dist, loc_n, MPI_INT, 0, comm);
    MPI_Gather(loc_pred, loc_n, MPI_INT, global_pred, loc_n, MPI_INT, 0, comm);

    /* Print results */
    if (my_rank == 0) {
        Print_dists(global_dist, n);
        Print_paths(global_pred, n);
        free(global_dist);
        free(global_pred);
    }

    printf("ennnnnddd fuck offf");


//     if(argc >= 4){
//     printf("Computing result for source 0.\n");
//     dijkstra(0, n, a, &l);
//     printf("Writing result to %s.\n", argv[3]);
//     print_numbers(argv[3], n, l);
//    }

//     free(a);
//     free(l);

    free(loc_mat);
    free(loc_pred);
    free(loc_dist);
    MPI_Type_free(&blk_col_mpi_t);
    MPI_Finalize();
    return 0;
}






/*---------------------------------------------------------------------
 * Function:  Read_n
 * Purpose:   Read in the number of rows in the matrix on process 0
 *            and broadcast this value to the other processes
 * In args:   my_rank:  the calling process' rank
 *            comm:  Communicator containing all calling processes
 * Ret val:   n:  the number of rows in the matrix
 */
int Read_n(int my_rank, MPI_Comm comm, int n) {


    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    return n;
}






/*---------------------------------------------------------------------
 * Function:  Build_blk_col_type
 * Purpose:   Build an MPI_Datatype that represents a block column of
 *            a matrix
 * In args:   n:  number of rows in the matrix and the block column
 *            loc_n = n/p:  number cols in the block column
 * Ret val:   blk_col_mpi_t:  MPI_Datatype that represents a block
 *            column
 */
MPI_Datatype Build_blk_col_type(int n, int loc_n) {
    MPI_Aint lb, extent;
    MPI_Datatype block_mpi_t;
    MPI_Datatype first_bc_mpi_t;
    MPI_Datatype blk_col_mpi_t;

    MPI_Type_contiguous(loc_n, MPI_INT, &block_mpi_t);
    MPI_Type_get_extent(block_mpi_t, &lb, &extent);

    /* MPI_Type_vector(numblocks, elts_per_block, stride, oldtype, *newtype) */
    MPI_Type_vector(n, loc_n, n, MPI_INT, &first_bc_mpi_t);

    /* This call is needed to get the right extent of the new datatype */
    MPI_Type_create_resized(first_bc_mpi_t, lb, extent, &blk_col_mpi_t);

    MPI_Type_commit(&blk_col_mpi_t);

    MPI_Type_free(&block_mpi_t);
    MPI_Type_free(&first_bc_mpi_t);

    return blk_col_mpi_t;
}






/*---------------------------------------------------------------------
 * Function:  Read_matrix
 * Purpose:   Read in an nxn matrix of ints on process 0, and
 *            distribute it among the processes so that each
 *            process gets a block column with n rows and n/p
 *            columns
 * In args:   n:  the number of rows/cols in the matrix and the submatrices
 *            loc_n = n/p:  the number of columns in the submatrices
 *            blk_col_mpi_t:  the MPI_Datatype used on process 0
 *            my_rank:  the caller's rank in comm
 *            comm:  Communicator consisting of all the processes
 * Out arg:   loc_mat:  the calling process' submatrix (needs to be
 *               allocated by the caller)
 */
void Read_matrix(float loc_mat[], int n, int loc_n,
                 MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm,const char * const filename) {
    float * mat ;
    int i, j;
    if (my_rank == 0) {
        // mat = (int*)malloc(n * n * sizeof(int));
        // for (i = 0; i < n; i++)
        //     for (j = 0; j < n; j++)
        //         scanf("%d", &mat[i * n + j]);

        /* load data */
        printf("Loading graph from %s.\n", filename);
        load(filename, &n, &mat);




    }
    //Print_matrix(mat,n,n);

    MPI_Scatter(mat, 1, blk_col_mpi_t, loc_mat, n * loc_n, MPI_INT, 0, comm);

    if (my_rank == 0) free(mat);
}






/*-------------------------------------------------------------------
 * Function:   Dijkstra_Init
 * Purpose:    Initialize all the matrices so that Dijkstras shortest path
 *             can be run
 *
 * In args:    loc_n:    local number of vertices
 *             my_rank:  the process rank
 *
 * Out args:   loc_mat:  local matrix containing edge costs between vertices
 *             loc_dist: loc_dist[v] = shortest distance from the source to each vertex v
 *             loc_pred: loc_pred[v] = predecessor of v on a shortest path from source to v
 *             loc_known: loc_known[v] = 1 if vertex has been visited, 0 else
 *
 *
 */
void Dijkstra_Init(float loc_mat[], float loc_pred[], float loc_dist[], float loc_known[],
                   int my_rank, int loc_n) {
    int loc_v;

    if (my_rank == 0)
        loc_known[0] = 1;
    else
        loc_known[0] = 0;

    for (loc_v = 1; loc_v < loc_n; loc_v++)
        loc_known[loc_v] = 0;

    for (loc_v = 0; loc_v < loc_n; loc_v++) {
        loc_dist[loc_v] = loc_mat[0 * loc_n + loc_v];
        loc_pred[loc_v] = 0;
    }
}






/*-------------------------------------------------------------------
 * Function:   Dijkstra
 * Purpose:    compute all the shortest paths from the source vertex 0
 *             to all vertices v
 *
 *
 * In args:    loc_mat:  local matrix containing edge costs between vertices
 *             loc_n:    local number of vertices
 *             n:        total number of vertices (globally)
 *             comm:     the communicator
 *
 * Out args:   loc_dist: loc_dist[v] = shortest distance from the source to each vertex v
 *             loc_pred: loc_pred[v] = predecessor of v on a shortest path from source to v
 *
 */
void Dijkstra(float loc_mat[], float loc_dist[], float loc_pred[], int loc_n, int n,
              MPI_Comm comm) {

     
    float new_dist, dist_glbl_u;
    int i,my_rank,loc_u,loc_v,glbl_u;
    float *loc_known;
    float my_min[2];
    float glbl_min[2];

    MPI_Comm_rank(comm, &my_rank);
    loc_known = (float*)malloc(loc_n * sizeof(float));

    Dijkstra_Init(loc_mat, loc_pred, loc_dist, loc_known, my_rank, loc_n);

    /* Run loop n - 1 times since we already know the shortest path to global
       vertex 0 from global vertex 0 */
    for (i = 0; i < n - 1; i++) {
        loc_u = Find_min_dist(loc_dist, loc_known, loc_n);

        if (loc_u != -1) {
            my_min[0] = loc_dist[loc_u];
            my_min[1] = loc_u + my_rank * loc_n;
        }
        else {
            my_min[0] = INFINITY;
            my_min[1] = -1;
        }

        /* Get the minimum distance found by the processes and store that
           distance and the global vertex in glbl_min
        */
        MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC, comm);

        dist_glbl_u = glbl_min[0];
        glbl_u = glbl_min[1];

        /* This test is to assure that loc_known is not accessed with -1 */
        if (glbl_u == -1)
            break;

        /* Check if global u belongs to process, and if so update loc_known */
        if ((glbl_u / loc_n) == my_rank) {
            loc_u = glbl_u % loc_n;
            loc_known[loc_u] = 1;
        }

        /* For each local vertex (global vertex = loc_v + my_rank * loc_n)
           Update the distances from source vertex (0) to loc_v. If vertex
           is unmarked check if the distance from source to the global u + the
           distance from global u to local v is smaller than the distance
           from the source to local v
         */
        for (loc_v = 0; loc_v < loc_n; loc_v++) {
            if (!loc_known[loc_v]) {
                new_dist = dist_glbl_u + loc_mat[glbl_u * loc_n + loc_v];
                if (new_dist < loc_dist[loc_v]) {
                    loc_dist[loc_v] = new_dist;
                    loc_pred[loc_v] = glbl_u;
                }
            }
        }
    }
    free(loc_known);
}






/*-------------------------------------------------------------------
 * Function:   Find_min_dist
 * Purpose:    find the minimum local distance from the source to the
 *             assigned vertices of the process that calls the method
 *
 *
 * In args:    loc_dist:  array with distances from source 0
 *             loc_known: array with values 1 if the vertex has been visited
 *                        0 if not
 *             loc_n:     local number of vertices
 *
 * Return val: loc_u: the vertex with the smallest value in loc_dist,
 *                    -1 if all vertices are already known
 *
 * Note:       loc_u = -1 is not supposed to be used when this function returns
 *
 */
int Find_min_dist(float loc_dist[], float loc_known[], int loc_n) {
    int loc_u = -1, loc_v;
    int shortest_dist = INFINITY;

    for (loc_v = 0; loc_v < loc_n; loc_v++) {
        if (!loc_known[loc_v]) {
            if (loc_dist[loc_v] < shortest_dist) {
                shortest_dist = loc_dist[loc_v];
                loc_u = loc_v;
            }
        }
    }
    return loc_u;
}






/*-------------------------------------------------------------------
 * Function:  Print_matrix
 * Purpose:   Print the contents of the matrix
 * In args:   mat, rows, cols
 *
 *
 */
void Print_matrix(float mat[], int rows, int cols) {
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++)
            if (mat[i * cols + j] == INFINITY)
                printf("i ");
            else
                printf("%d ", mat[i * cols + j]);
        printf("\n");
    }

    printf("\n");
}






/*-------------------------------------------------------------------
 * Function:    Print_dists
 * Purpose:     Print the length of the shortest path from 0 to each
 *              vertex
 * In args:     n:  the number of vertices
 *              dist:  distances from 0 to each vertex v:  dist[v]
 *                 is the length of the shortest path 0->v
 */
void Print_dists(float global_dist[], int n) {
    int v;

    printf("  v    dist 0->v\n");
    printf("----   ---------\n");

    for (v = 1; v < n; v++) {
        if (global_dist[v] == INFINITY) {
            printf("%3d       %5s\n", v, "inf");
        }
        else
            printf("%3d       %4d\n", v, global_dist[v]);
        }
    printf("\n");
}






/*-------------------------------------------------------------------
 * Function:    Print_paths
 * Purpose:     Print the shortest path from 0 to each vertex
 * In args:     n:  the number of vertices
 *              pred:  list of predecessors:  pred[v] = u if
 *                 u precedes v on the shortest path 0->v
 */
void Print_paths(float global_pred[], int n) {
    int v, w, count, i;
    float *path;

    path =  (float*)malloc(n * sizeof(float));

    printf("  v     Path 0->v\n");
    printf("----    ---------\n");
    for (v = 1; v < n; v++) {
        printf("%3d:    ", v);
        count = 0;
        w = v;
        while (w != 0) {
            path[count] = w;
            count++;
            w = global_pred[w];
        }
        printf("0 ");
        for (i = count-1; i >= 0; i--)
            printf("%d ", path[i]);
        printf("\n");
    }

    free(path);
}