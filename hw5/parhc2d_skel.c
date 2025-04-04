#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define DEBUG_RANK 0
#define DEBUG_INDEX_I 10
#define DEBUG_INDEX_J 10

void grid(int nx, int nxglob, int istglob, int ienglob, double xstglob, double xenglob, double *x, double *dx)
{
  int i, iglob;
  *dx = (xenglob - xstglob)/(double)(nxglob-1);

  for(i=0; i<nx; i++)
  {
    iglob = istglob + i;
    x[i] = xstglob + (double)iglob * (*dx);
  }
}

void enforce_bcs(int nx, int ny, int istglob, int ienglob, int jstglob, int jenglob, int nxglob, int nyglob, double *x, double *y, double **T)
{
  int i, j;

  if (istglob == 0)
    for(j=0; j<ny; j++)
      T[0][j] = 0.0;

  if (ienglob == nxglob - 1)
    for(j=0; j<ny; j++)
      T[nx-1][j] = 0.0;

  if (jstglob == 0)
    for(i=0; i<nx; i++)
      T[i][0] = 0.0;

  if (jenglob == nyglob - 1)
    for(i=0; i<nx; i++)
      T[i][ny-1] = 0.0;
}

void set_initial_condition(int nx, int ny, int istglob, int ienglob, int jstglob, int jenglob, int nxglob, int nyglob, double *x, double *y, double **T, double dx, double dy)
{
  int i, j;
  double del = 1.0;

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      T[i][j] = 0.25 * (tanh((x[i]-0.4)/(del*dx)) - tanh((x[i]-0.6)/(del*dx))) 
                     * (tanh((y[j]-0.4)/(del*dy)) - tanh((y[j]-0.6)/(del*dy)));

  enforce_bcs(nx, ny, istglob, ienglob, jstglob, jenglob, nxglob, nyglob, x, y, T);
}

void halo_exchange_2d_x(int rank, int rank_x, int rank_y, int size, int px, int py, int nx, int ny, int nxglob, int nyglob, double *x, double *y, double **T, double *xleftghost, double *xrightghost, double *sendbuf_x, double *recvbuf_x)
{
  MPI_Status status;
  int left_nb = (rank_x == 0) ? MPI_PROC_NULL : rank - 1;
  int right_nb = (rank_x == px - 1) ? MPI_PROC_NULL : rank + 1;
  int j;

  for (j = 0; j < ny; j++) sendbuf_x[j] = T[0][j];
  MPI_Recv(recvbuf_x, ny, MPI_DOUBLE, right_nb, 0, MPI_COMM_WORLD, &status);
  MPI_Send(sendbuf_x, ny, MPI_DOUBLE, left_nb, 0, MPI_COMM_WORLD);
  for (j = 0; j < ny; j++) xrightghost[j] = recvbuf_x[j];

  for (j = 0; j < ny; j++) sendbuf_x[j] = T[nx - 1][j];
  MPI_Recv(recvbuf_x, ny, MPI_DOUBLE, left_nb, 1, MPI_COMM_WORLD, &status);
  MPI_Send(sendbuf_x, ny, MPI_DOUBLE, right_nb, 1, MPI_COMM_WORLD);
  for (j = 0; j < ny; j++) xleftghost[j] = recvbuf_x[j];
}

void halo_exchange_2d_y(int rank, int rank_x, int rank_y, int size, int px, int py, int nx, int ny, int nxglob, int nyglob, double *x, double *y, double **T, double *ybotghost, double *ytopghost, double *sendbuf_y, double *recvbuf_y)
{
  MPI_Status status;
  int bot_nb = (rank_y == 0) ? MPI_PROC_NULL : rank - px;
  int top_nb = (rank_y == py - 1) ? MPI_PROC_NULL : rank + px;
  int i;

  for (i = 0; i < nx; i++) sendbuf_y[i] = T[i][0];
  MPI_Recv(recvbuf_y, nx, MPI_DOUBLE, top_nb, 2, MPI_COMM_WORLD, &status);
  MPI_Send(sendbuf_y, nx, MPI_DOUBLE, bot_nb, 2, MPI_COMM_WORLD);
  for (i = 0; i < nx; i++) ytopghost[i] = recvbuf_y[i];

  for (i = 0; i < nx; i++) sendbuf_y[i] = T[i][ny - 1];
  MPI_Recv(recvbuf_y, nx, MPI_DOUBLE, bot_nb, 3, MPI_COMM_WORLD, &status);
  MPI_Send(sendbuf_y, nx, MPI_DOUBLE, top_nb, 3, MPI_COMM_WORLD);
  for (i = 0; i < nx; i++) ybotghost[i] = recvbuf_y[i];
}

void get_rhs(int nx, int nxglob, int ny, int nyglob, int istglob, int ienglob, int jstglob, int jenglob, double dx, double dy, double *xleftghost, double *xrightghost, double *ybotghost, double *ytopghost, double kdiff, double *x, double *y, double **T, double **rhs, int rank)
{
  int i, j;
  double dxsq = dx*dx, dysq = dy*dy;

  for(i=1; i<nx-1; i++)
    for(j=1; j<ny-1; j++)
      rhs[i][j] = kdiff*(T[i+1][j]+T[i-1][j]-2.0*T[i][j])/dxsq +
                  kdiff*(T[i][j+1]+T[i][j-1]-2.0*T[i][j])/dysq;

  i = 0;
  for(j=1; j<ny-1; j++)
    rhs[i][j] = (istglob==0) ? 0.0 :
                kdiff*(T[i+1][j]+xleftghost[j]-2.0*T[i][j])/dxsq +
                kdiff*(T[i][j+1]+T[i][j-1]-2.0*T[i][j])/dysq;

  i = nx-1;
  for(j=1; j<ny-1; j++)
    rhs[i][j] = (ienglob==nxglob-1) ? 0.0 :
                kdiff*(xrightghost[j]+T[i-1][j]-2.0*T[i][j])/dxsq +
                kdiff*(T[i][j+1]+T[i][j-1]-2.0*T[i][j])/dysq;

  j = 0;
  for(i=1; i<nx-1; i++)
    rhs[i][j] = (jstglob==0) ? 0.0 :
                kdiff*(T[i+1][j]+T[i-1][j]-2.0*T[i][j])/dxsq +
                kdiff*(T[i][j+1]+ybotghost[i]-2.0*T[i][j])/dysq;

  j = ny-1;
  for(i=1; i<nx-1; i++)
    rhs[i][j] = (jenglob==nyglob-1) ? 0.0 :
                kdiff*(T[i+1][j]+T[i-1][j]-2.0*T[i][j])/dxsq +
                kdiff*(ytopghost[i]+T[i][j-1]-2.0*T[i][j])/dysq;

  i = 0; j = 0;
  rhs[i][j] = (istglob==0 || jstglob==0) ? 0.0 :
              kdiff*(T[i+1][j]+xleftghost[j]-2.0*T[i][j])/dxsq +
              kdiff*(T[i][j+1]+ybotghost[i]-2.0*T[i][j])/dysq;

  i = nx-1; j = 0;
  rhs[i][j] = (ienglob==nxglob-1 || jstglob==0) ? 0.0 :
              kdiff*(xrightghost[j]+T[i-1][j]-2.0*T[i][j])/dxsq +
              kdiff*(T[i][j+1]+ybotghost[i]-2.0*T[i][j])/dysq;

  i = 0; j = ny-1;
  rhs[i][j] = (istglob==0 || jenglob==nyglob-1) ? 0.0 :
              kdiff*(T[i+1][j]+xleftghost[j]-2.0*T[i][j])/dxsq +
              kdiff*(ytopghost[i]+T[i][j-1]-2.0*T[i][j])/dysq;

  i = nx-1; j = ny-1;
  rhs[i][j] = (ienglob==nxglob-1 || jenglob==nyglob-1) ? 0.0 :
              kdiff*(xrightghost[j]+T[i-1][j]-2.0*T[i][j])/dxsq +
              kdiff*(ytopghost[i]+T[i][j-1]-2.0*T[i][j])/dysq;
}

void timestep_FwdEuler(int rank, int size, int rank_x, int rank_y, int px, int py, int nx, int nxglob, int ny, int nyglob, int istglob, int ienglob, int jstglob, int jenglob, double dt, double dx, double dy, double *xleftghost, double *xrightghost, double *ybotghost, double *ytopghost, double kdiff, double *x, double *y, double **T, double **rhs, double *sendbuf_x, double *recvbuf_x, double *sendbuf_y, double *recvbuf_y)
{
  int i, j;

  halo_exchange_2d_x(rank, rank_x, rank_y, size, px, py, nx, ny, nxglob, nyglob, x, y, T, xleftghost, xrightghost, sendbuf_x, recvbuf_x);
  halo_exchange_2d_y(rank, rank_x, rank_y, size, px, py, nx, ny, nxglob, nyglob, x, y, T, ybotghost, ytopghost, sendbuf_y, recvbuf_y);

  get_rhs(nx, nxglob, ny, nyglob, istglob, ienglob, jstglob, jenglob, dx, dy, xleftghost, xrightghost, ybotghost, ytopghost, kdiff, x, y, T, rhs, rank);

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      T[i][j] = T[i][j] + dt * rhs[i][j];

  enforce_bcs(nx, ny, istglob, ienglob, jstglob, jenglob, nxglob, nyglob, x, y, T);
}
void output_soln(int rank, int nx, int ny, int it, double tcurr, double *x, double *y, double **T, int size)
{
  int i,j;
  FILE* fp;
  char fname[100];

  sprintf(fname, "./output_file/T_x_y_%06d_%04d_par_%d.dat", it, rank, size);
  fp = fopen(fname, "w");
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      fprintf(fp, "%lf %lf %lf\n", x[i], y[j], T[i][j]);
  fclose(fp);
}

void get_processor_grid_ranks(int rank, int size, int px, int py, int *rank_x, int *rank_y)
{
  *rank_y = rank / px;
  *rank_x = rank % px;
}
int main(int argc, char** argv)
{
  int nx, ny, nxglob, nyglob, rank, size, px, py, rank_x, rank_y;
  double *x, *y, **T, **rhs, start_time, end_time, time_step, tcurr, kdiff;
  double xstglob, xenglob, ystglob, yenglob, t_print, xlen, ylen, dx, dy, **Tnew;
  double *xleftghost, *xrightghost, *ybotghost, *ytopghost;
  double *sendbuf_x, *recvbuf_x, *sendbuf_y, *recvbuf_y;
  int i, j, it, num_time_steps, it_print;
  int istglob, ienglob, jstglob, jenglob;
  FILE* fid;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank == 0)
  {
    fid = fopen("input2d.in", "r");
    fscanf(fid, "%d %d\n", &nxglob, &nyglob);
    fscanf(fid, "%lf %lf %lf %lf\n", &xstglob, &xenglob, &ystglob, &yenglob);
    fscanf(fid, "%lf %lf %lf %lf\n", &start_time, &end_time, &time_step, &t_print);
    fscanf(fid, "%lf\n", &kdiff);
    fscanf(fid, "%d %d\n", &px, &py);
    fclose(fid);

    nx = nxglob / px;
    ny = nyglob / py;
    xlen = (xenglob - xstglob) / (double)px;
    ylen = (yenglob - ystglob) / (double)py;

    num_time_steps = (int)((end_time - start_time) / time_step) + 1;
    it_print = (int)(t_print / time_step);

    if(px * py != size)
    {
      printf("Processor grid does not match total processors. Exiting...\n");
      exit(1);
    }
  }

  int *arr_int = malloc(8 * sizeof(int));
  if(rank == 0) {
    arr_int[0] = nxglob; arr_int[1] = nx;
    arr_int[2] = nyglob; arr_int[3] = ny;
    arr_int[4] = num_time_steps; arr_int[5] = it_print;
    arr_int[6] = px; arr_int[7] = py;
  }
  MPI_Bcast(arr_int, 8, MPI_INT, 0, MPI_COMM_WORLD);
  if(rank != 0) {
    nxglob = arr_int[0]; nx = arr_int[1];
    nyglob = arr_int[2]; ny = arr_int[3];
    num_time_steps = arr_int[4]; it_print = arr_int[5];
    px = arr_int[6]; py = arr_int[7];
  }
  free(arr_int);

  get_processor_grid_ranks(rank, size, px, py, &rank_x, &rank_y);

  double *arr_dbl = malloc(11 * sizeof(double));
  if(rank == 0) {
    arr_dbl[0] = start_time; arr_dbl[1] = end_time; arr_dbl[2] = time_step; arr_dbl[3] = t_print;
    arr_dbl[4] = xlen; arr_dbl[5] = xstglob; arr_dbl[6] = xenglob;
    arr_dbl[7] = ylen; arr_dbl[8] = ystglob; arr_dbl[9] = yenglob; arr_dbl[10] = kdiff;
  }
  MPI_Bcast(arr_dbl, 11, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(rank != 0) {
    start_time = arr_dbl[0]; end_time = arr_dbl[1]; time_step = arr_dbl[2]; t_print = arr_dbl[3];
    xlen = arr_dbl[4]; xstglob = arr_dbl[5]; xenglob = arr_dbl[6];
    ylen = arr_dbl[7]; ystglob = arr_dbl[8]; yenglob = arr_dbl[9]; kdiff = arr_dbl[10];
  }
  free(arr_dbl);

  istglob = rank_x * (nxglob / px);
  ienglob = (rank_x + 1) * (nxglob / px) - 1;
  jstglob = rank_y * (nyglob / py);
  jenglob = (rank_y + 1) * (nyglob / py) - 1;

  x = (double *)malloc(nx * sizeof(double));
  y = (double *)malloc(ny * sizeof(double));
  T = (double **)malloc(nx * sizeof(double *));
  rhs = (double **)malloc(nx * sizeof(double *));
  Tnew = (double **)malloc(nx * sizeof(double *));
  for(i = 0; i < nx; i++) {
    T[i] = (double *)malloc(ny * sizeof(double));
    rhs[i] = (double *)malloc(ny * sizeof(double));
    Tnew[i] = (double *)malloc(ny * sizeof(double));
  }

  xleftghost = (double *)malloc(ny * sizeof(double));
  xrightghost = (double *)malloc(ny * sizeof(double));
  ybotghost = (double *)malloc(nx * sizeof(double));
  ytopghost = (double *)malloc(nx * sizeof(double));
  sendbuf_x = (double *)malloc(ny * sizeof(double));
  recvbuf_x = (double *)malloc(ny * sizeof(double));
  sendbuf_y = (double *)malloc(nx * sizeof(double));
  recvbuf_y = (double *)malloc(nx * sizeof(double));

  grid(nx, nxglob, istglob, ienglob, xstglob, xenglob, x, &dx);
  grid(ny, nyglob, jstglob, jenglob, ystglob, yenglob, y, &dy);

  set_initial_condition(nx, ny, istglob, ienglob, jstglob, jenglob, nxglob, nyglob, x, y, T, dx, dy);
  output_soln(rank, nx, ny, 0, start_time, x, y, T, size);

  for(it = 0; it < num_time_steps; it++)
  {
    tcurr = start_time + (it + 1) * time_step;

    timestep_FwdEuler(rank, size, rank_x, rank_y, px, py, nx, nxglob, ny, nyglob,
                      istglob, ienglob, jstglob, jenglob, time_step, dx, dy,
                      xleftghost, xrightghost, ybotghost, ytopghost, kdiff, x, y, T, rhs,
                      sendbuf_x, recvbuf_x, sendbuf_y, recvbuf_y);

    if(it % it_print == 0)
      output_soln(rank, nx, ny, it, tcurr, x, y, T, size);
  }

  for(i = 0; i < nx; i++) {
    free(T[i]);
    free(rhs[i]);
    free(Tnew[i]);
  }
  free(T); free(rhs); free(Tnew);
  free(x); free(y);
  free(xleftghost); free(xrightghost); free(ybotghost); free(ytopghost);
  free(sendbuf_x); free(recvbuf_x); free(sendbuf_y); free(recvbuf_y);

  MPI_Finalize();
  return 0;
}
