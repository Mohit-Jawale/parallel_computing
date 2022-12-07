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
      assert(1 == ret);
    }
  }

  /* close file */
  ret = fclose(fp);
  assert(!ret);
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

    MPI_Gather(loc_dist, loc_n, MPI_INT, global_dist, loc_n, MPI_INT, 0, comm);
    MPI_Gather(loc_pred, loc_n, MPI_INT, global_pred, loc_n, MPI_INT, 0, comm);


    if(argc >= 4 and my_rank==0){
    cout<<"Computing result for source 0"<<endl;
    printf("Writing result to %s.\n", argv[3]);
    print_numbers(argv[3], n, global_dist);
   }

    free(global_dist);
    free(global_pred);
    free(loc_mat);
    free(loc_pred);
    free(loc_dist);
    MPI_Type_free(&blk_col_mpi_t);
    MPI_Finalize();
    return 0;
}





int Read_n(int my_rank, MPI_Comm comm, int n) {


    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    return n;
}




MPI_Datatype Build_blk_col_type(int n, int loc_n) {
    MPI_Aint lb, extent;
    MPI_Datatype block_mpi_t;
    MPI_Datatype first_bc_mpi_t;
    MPI_Datatype blk_col_mpi_t;

    MPI_Type_contiguous(loc_n, MPI_INT, &block_mpi_t);
    MPI_Type_get_extent(block_mpi_t, &lb, &extent);

    MPI_Type_vector(n, loc_n, n, MPI_INT, &first_bc_mpi_t);


    MPI_Type_create_resized(first_bc_mpi_t, lb, extent, &blk_col_mpi_t);

    MPI_Type_commit(&blk_col_mpi_t);

    MPI_Type_free(&block_mpi_t);
    MPI_Type_free(&first_bc_mpi_t);

    return blk_col_mpi_t;
}



void Read_matrix(float loc_mat[], int n, int loc_n,
                 MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm,const char * const filename) {
    float * mat ;
    int i, j;
    if (my_rank == 0) {
        printf("Loading graph from %s.\n", filename);
        load(filename, &n, &mat);

    }

    MPI_Scatter(mat, 1, blk_col_mpi_t, loc_mat, n * loc_n, MPI_INT, 0, comm);

    if (my_rank == 0) free(mat);
}





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


        MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC, comm);

        dist_glbl_u = glbl_min[0];
        glbl_u = glbl_min[1];


        if (glbl_u == -1)
            break;

        if ((glbl_u / loc_n) == my_rank) {
            loc_u = glbl_u % loc_n;
            loc_known[loc_u] = 1;
        }

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
