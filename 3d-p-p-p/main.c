// --------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
// --------------------------------------------------------------------------
#include <silo.h>
// --------------------------------------------------------------------------
#ifndef SQR
#define SQR(x) ((x) * (x))
#endif
// --------------------------------------------------------------------------
#ifndef PI
#define PI 3.14159265358979323846264338327950288419716940
#endif
// --------------------------------------------------------------------------
/*
 * The following macros addresses the problem of indexing arrays when
 * numerically solving mathematical physics problems.
 *
 * Boundary points may be included in the array, or may not be included.
 * If you don't include left boundary point in the array (i.e., you don't
 * allocate memory for it), then you should refer to the first internal
 * point with index `0', but it is more convinient to reference it with
 * index `1'. In multidimensional case you also need to specify the size of
 * the array.
 *
 * Thus, `IDX1(i0, b0, n0)' refers to the point with convinient index `i0'
 * in case of the whole array is referenced from `b0' and has `n0' elements
 * in total.
 *
 * Samples:
 *
 * 0 1           N
 * x-o-o-o-...-o-x  -  IDX1(i, 1, N-1)
 * x-o-o-o-...-o-o  -  IDX1(i, 1, N)
 * o-o-o-o-...-o-x  -  IDX1(i, 0, N)
 * o-o-o-o-...-o-o  -  IDX1(i, 0, N+1)
 *
 * (`x' denotes points that are not included in the array, `o' denotes
 * points that are included in the array.)
 *
 * Actually, the last argument of the macro (`n0') is not used, but it was
 * added for consistency with multidimensional cases.
 */
#ifndef IDX1
#define IDX1(i0, b0, n0) ((i0)-(b0))
#endif
// --------------------------------------------------------------------------
/*
 * We use the C programming language (not Fortran), so it is natural to use
 * the column-major order, i.e. the last index is the fastest.
 *
 *   0 1           N 
 * 0 x-x-x-x-...-x-x -> i1
 *   | | | |     | |
 * 1 x-o-o-o-...-o-x
 *   | | | |     | |
 *   . . . .     . .
 *   | | | |     | |
 *   x-o-o-o-...-o-x
 *   | | | |     | |
 * N x-x-x-x-...-x-x
 */
#ifndef IDX2
#define IDX2(i0, i1, b0, b1, n0, n1) (((i0)-(b0))*(n1) + (i1)-(b1))
#endif
// --------------------------------------------------------------------------
#ifndef IDX3
#define IDX3(i0, i1, i2, b0, b1, b2, n0, n1, n2) (((i0)-(b0))*(n1)*(n2) + ((i1)-(b1))*(n2) + (i2)-(b2))
#endif
// --------------------------------------------------------------------------
typedef struct XDraw_Record_s
{
    int m1;
    float x;
    float y;
    int m2;
} XDraw_Record;
// --------------------------------------------------------------------------
const struct
{
    int m1;
    int m2;
} XDraw_EOL = { 0, 0 };
// --------------------------------------------------------------------------
/*
 * It is natural to store number of segments, because it gives an ability to
 * immidiately calculate grid step as h = Ns / L.
 *
 * Number of actual points that are stored in an array is easily obtained
 * from these variables.
 */
typedef struct Grid_Info_s
{
    int Ns0, Ns1, Ns2;   // # of grid segments in each spatial dir'n
    double L0, L1, L2;  // # of internal grid points in each spatial dir'n
    double h0, h1, h2;
} Grid_Info;
// --------------------------------------------------------------------------
typedef double fun_t(double, double, double, double, Grid_Info);
// --------------------------------------------------------------------------
double U_basic(double x0, double x1, double x2, double t, Grid_Info g)
{
    assert(t >= 0.0);
    assert(g.L0 > 0.0);
    assert(g.L1 > 0.0);
    assert(g.L2 > 0.0);

    /*
    return
        sin(2*PI * x0 / g.L0) *
        sin(2*PI * x1 / g.L1) *
        sin(2*PI * x2 / g.L2) *
        cos(2*PI * t * sqrt(
                    1.0 / SQR(g.L0) +
                    1.0 / SQR(g.L1) +
                    1.0 / SQR(g.L2)
        ));
        */

//    return sin(PI * x0 / g.L0) * sin(PI * x1 / g.L1) * sin(PI * x2 / g.L2) *
//        cos(PI * t * sqrt(1 / SQR(g.L0) + 1 / SQR(g.L1) + 1 / SQR(g.L2)));

    const double x0_star = 0.35;
    const double x1_star = 0.25;
    const double x2_star = 0.15;

    if (fabs(x0 - x0_star) > g.L0/2) x0 -= g.L0;
    if (fabs(x1 - x1_star) > g.L1/2) x1 -= g.L1;
    if (fabs(x2 - x2_star) > g.L2/2) x2 -= g.L2;

    return exp(-(SQR(x0 - x0_star) + SQR(x1 - x1_star) + SQR(x2 - x2_star)) / 0.01);

    /*
    if (x0 == 0.5 && x1 == 0.5 && x2 == 0.5) return 10.0;
    else return 0.0;
    */
}
// --------------------------------------------------------------------------
static void Save_Graph(
        const double *y, XDraw_Record *b, int N,
        const char *fname, const char *mode)
{
    assert(y);
    assert(b);
    assert(N > 0);
    assert(fname);
    assert(mode);

    int i, res, num;
    FILE *fout;

    // Fill buffer
    for (i = 0; i <= N; ++i)
    {
        int idx = IDX1(i, 0, N+1);
        b[idx].y = (float)y[idx];
    }

    // Save file
    fout = fopen(fname, mode);
    if (!fout)
    {
        fprintf(stderr, "Error: Unable to open file in a specified mode: "
                "%s, %s\n", fname, mode);
        exit(EXIT_FAILURE);
    }

    num = fwrite(b, sizeof(*b), N+1, fout);
    if (num != N+1)
    {
        fprintf(stderr, "Error: Unable to write to file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    num = fwrite(&XDraw_EOL, sizeof(XDraw_EOL), 1, fout);
    if (num != 1)
    {
        fprintf(stderr, "Error: Unable to write to file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    res = fclose(fout);
    if (res != 0)
        fprintf(stderr, "Warning: Unable to close file %s\n", fname);
}
// --------------------------------------------------------------------------
static double Calculate_Error(
        const double *u, const double *y, double *z, Grid_Info g)
{
    assert(u);
    assert(y);
    assert(z);
    assert(g.Ns0 > 0);
    assert(g.Ns1 > 0);
    assert(g.Ns2 > 0);

    int i0, i1, i2;

    double z_Ch = 0.0;

    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g.Ns1; ++i1)
        {
            for (i2 = 0; i2 <= g.Ns2; ++i2)
            {
                int idx = IDX3(i0, i1, i2, 0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);

                z[idx] = y[idx] - u[idx];

                if (fabs(z[idx]) > z_Ch) z_Ch = fabs(z[idx]);
            }
        }
    }

    return z_Ch;
}
// --------------------------------------------------------------------------
static void Prepare_Initial_Data(
        double *x[], double *y[], Grid_Info g, double tau, fun_t fun)
{
    assert(x);
    assert(x[0]);
    assert(x[1]);
    assert(x[2]);
    assert(y);
    assert(y[0]);
    assert(y[1]);
    assert(g.Ns0 > 0);
    assert(g.Ns1 > 0);
    assert(g.Ns2 > 0);
    assert(g.h0 > 0.0);
    assert(g.h1 > 0.0);
    assert(g.h2 > 0.0);
    assert(tau > 0.0);
    assert(fun);

    double gamma0 = SQR(tau / g.h0);
    double gamma1 = SQR(tau / g.h1);
    double gamma2 = SQR(tau / g.h2);
    const int mid = 1;
    const int bot = 0;
    double t = 0.0;

    int i0, i1, i2, idx;
    int idx0l, idx0r, idx1l, idx1r, idx2l, idx2r;

    // Grid
    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        idx = IDX1(i0, 0, g.Ns0+1);
        x[0][idx] = i0 * g.h0;
    }

    for (i1 = 0; i1 <= g.Ns1; ++i1)
    {
        idx = IDX1(i1, 0, g.Ns1+1);
        x[1][idx] = i1 * g.h1;
    }

    for (i2 = 0; i2 <= g.Ns2; ++i2)
    {
        idx = IDX1(i2, 0, g.Ns2+1);
        x[2][idx] = i2 * g.h2;
    }

    // t^0
    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g.Ns1; ++i1)
        {
            for (i2 = 0; i2 <= g.Ns2; ++i2)
            {
                idx = IDX3(i0, i1, i2, 0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);

                y[bot][idx] = fun(x[0][i0], x[1][i1], x[2][i2], t, g);
            }
        }
    }

    //memcpy(y[mid], y[bot], (g.Ns0+1)*(g.Ns1+1)*(g.Ns2+1)*sizeof(y[mid][0]));

    // t^1
    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g.Ns1; ++i1)
        {
            for (i2 = 0; i2 <= g.Ns2; ++i2)
            {
                idx   = IDX3(i0,   i1,   i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx0l = IDX3((i0-1+g.Ns0)%g.Ns0, i1,   i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx0r = IDX3((i0+1)%g.Ns0, i1,   i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx1l = IDX3(i0,   (i1-1+g.Ns1)%g.Ns1, i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx1r = IDX3(i0,   (i1+1)%g.Ns1, i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx2l = IDX3(i0,   i1,   (i2-1+g.Ns2)%g.Ns2, 0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx2r = IDX3(i0,   i1,   (i2+1)%g.Ns2, 0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);

                y[mid][idx] = y[bot][idx] +
                    gamma0 / 2 * (y[bot][idx0l] - 2 * y[bot][idx] + y[bot][idx0r]) +
                    gamma1 / 2 * (y[bot][idx1l] - 2 * y[bot][idx] + y[bot][idx1r]) +
                    gamma2 / 2 * (y[bot][idx2l] - 2 * y[bot][idx] + y[bot][idx2r]);
            }
        }
    }
}
// --------------------------------------------------------------------------
static void Get_Analytical_Solution(
        double *x[], double *u, Grid_Info g, int n, double tau, fun_t fun)
{
    assert(x);
    assert(x[0]);
    assert(x[1]);
    assert(x[2]);
    assert(u);
    assert(g.Ns0 > 0);
    assert(g.Ns1 > 0);
    assert(g.Ns2 > 0);
    assert(n >= 0);
    assert(tau > 0.0);
    assert(fun);

    double t = tau * n;
    int i0, i1, i2;

    // t^n
    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g.Ns1; ++i1)
        {
            for (i2 = 0; i2 <= g.Ns2; ++i2)
            {
                int idx = IDX3(i0, i1, i2, 0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                u[idx] = fun(x[0][i0], x[1][i1], x[2][i2], t, g);
            }
        }
    }
}
// --------------------------------------------------------------------------
static void Solve(double *y[], Grid_Info g, double tau, int top)
{
    assert(y);
    assert(y[0]);
    assert(y[1]);
    assert(g.Ns0 > 0);
    assert(g.Ns1 > 0);
    assert(g.Ns2 > 0);
    assert(g.h0 > 0);
    assert(g.h1 > 0);
    assert(g.h2 > 0);
    assert(tau > 0.0);
    assert(top == 0 || top == 1);

    double gamma0 = SQR(tau / g.h0);
    double gamma1 = SQR(tau / g.h1);
    double gamma2 = SQR(tau / g.h2);

    int mid = (top + 1) % 2;
    int bot = (mid + 1) % 2;

    int i0, i1, i2, idx;
    int idx0l, idx0r, idx1l, idx1r, idx2l, idx2r;

    // Internal points
    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g.Ns1; ++i1)
        {
            for (i2 = 0; i2 <= g.Ns2; ++i2)
            {
                idx   = IDX3(i0,   i1,   i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx0l = IDX3((i0-1+g.Ns0)%g.Ns0, i1,   i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx0r = IDX3((i0+1)%g.Ns0, i1,   i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx1l = IDX3(i0,   (i1-1+g.Ns1)%g.Ns1, i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx1r = IDX3(i0,   (i1+1)%g.Ns1, i2,   0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx2l = IDX3(i0,   i1,   (i2-1+g.Ns2)%g.Ns2, 0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);
                idx2r = IDX3(i0,   i1,   (i2+1)%g.Ns2, 0, 0, 0, g.Ns0+1, g.Ns1+1, g.Ns2+1);

                y[top][idx] = 2 * y[mid][idx] - y[bot][idx] +
                    gamma0 * (y[mid][idx0l] - 2 * y[mid][idx] + y[mid][idx0r]) +
                    gamma1 * (y[mid][idx1l] - 2 * y[mid][idx] + y[mid][idx1r]) +
                    gamma2 * (y[mid][idx2l] - 2 * y[mid][idx] + y[mid][idx2r]);
            }
        }
    }
}
// --------------------------------------------------------------------------
void Write_Silo(int n, double tau, double *x[], Grid_Info g,
        double *y, double *u, double *z)
{
    assert(n >= 0);
    assert(tau > 0.0);
    assert(n <= 9999);
    assert(x);
    assert(x[0]);
    assert(x[1]);
    assert(x[2]);
    assert(g.Ns0 > 0);
    assert(g.Ns1 > 0);
    assert(g.Ns2 > 0);
    assert(y);
    assert(u);
    assert(z);

    char *fname;
    const char *pattern = "%04d.silo";

    double dtime = n * tau;

    DBfile *dbfile;
    int dims[] = { g.Ns0+1, g.Ns1+1, g.Ns2+1 };
    int ndims = sizeof(dims) / sizeof(*dims);
    DBoptlist *optlist = DBMakeOptlist(3);
    int major_order = 1; 
    // 1. VisIt Silo plugin seems to always use the column-major order...
    // 2. The previouse statemant seems to be wrong. The stuff works,
    //    but I don't know why
    // - set major_order to 1 (col-maj)
    // - use row-maj
    // - add optlist to all operation
    
    // Create file
    // TODO: replace magic `4' with max number of decimal digits in n
    fname = (char *)malloc((4 + 5 + 1) * sizeof(*fname));
    sprintf(fname, pattern, n);
    dbfile = DBCreate(fname, DB_CLOBBER, DB_LOCAL,
            "Comment about the data", DB_PDB); // DB_PDB, DB_HDF5

    if (dbfile == NULL)
    {
        fprintf(stderr, "Could not create Silo file: %s\n", fname);
        exit(EXIT_FAILURE);
    }

    free(fname);

    // Save cycle and time values to an option list
    DBAddOption(optlist, DBOPT_DTIME, &dtime);
    DBAddOption(optlist, DBOPT_CYCLE, &n);
    DBAddOption(optlist, DBOPT_MAJORORDER, &major_order);

    // Write a rectilinear mesh
    DBPutQuadmesh(dbfile, "quadmesh", NULL, x, dims, ndims,
                DB_DOUBLE, DB_COLLINEAR, optlist);

    // Write scalar values
    DBPutQuadvar1(dbfile, "y", "quadmesh", y, dims,
        ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
    /*
    DBPutQuadvar1(dbfile, "u", "quadmesh", u, dims,
        ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
    DBPutQuadvar1(dbfile, "z", "quadmesh", z, dims,
        ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
        */

    DBFreeOptlist(optlist);

    DBClose(dbfile);
}
// --------------------------------------------------------------------------
int main()
{
    double *x[3];   // Coordinates
    double *u;      // Analytical solution at a particular time-step
    double *z;      // Error at each point at a particular time-step
    double *y[2];   // Numerical solution at two consecutive time-steps
    double *z_Ch;   // Error

    XDraw_Record *b;    // Buffer to save data in XDraw format

    Grid_Info g = { 30, 30, 30, 1.0, 1.0, 1.0, -1, -1, -1 };

    const double T = 1*g.L0;
    const int K = 2*g.Ns0;
    const double tau = T / (1*K);   // TODO: what is stability condition?

    int n;  // Counter through time-steps
    const char *mode = NULL;

    g.h0 = g.L0 / g.Ns0;
    g.h1 = g.L1 / g.Ns1;
    g.h2 = g.L2 / g.Ns2;

    // Allocate memory
    puts("Allocate memory");
    x[0] = malloc((g.Ns0+1) * sizeof(*(x[0])));
    x[1] = malloc((g.Ns1+1) * sizeof(*(x[1])));
    x[2] = malloc((g.Ns2+1) * sizeof(*(x[2])));
    u = malloc((g.Ns0+1) * (g.Ns1+1) * (g.Ns2+1) * sizeof(*u));
    z = malloc((g.Ns0+1) * (g.Ns1+1) * (g.Ns2+1) * sizeof(*z));
    y[0] = malloc((g.Ns0+1) * (g.Ns1+1) * (g.Ns2+1) * sizeof(*(y[0])));
    y[1] = malloc((g.Ns0+1) * (g.Ns1+1) * (g.Ns2+1) * sizeof(*(y[1])));
    z_Ch = malloc((K+1) * sizeof(*z_Ch));

    // Prepare initial data
    puts("Prepare initial data");
    Prepare_Initial_Data(x, y, g, tau, U_basic);
    
    // t^0
    puts("t^0");
    n = 0;
    mode = "wb";
    Get_Analytical_Solution(x, u, g, n, tau, U_basic);
    z_Ch[n] = Calculate_Error(u, y[n], z, g);
    Write_Silo(n, tau, x, g, y[n], u, z);

    // t^1
    puts("t^1");
    n = 1;
    mode = "ab";
    Get_Analytical_Solution(x, u, g, n, tau, U_basic);
    z_Ch[n] = Calculate_Error(u, y[n], z, g);
    Write_Silo(n, tau, x, g, y[n], u, z);

    // Solve equation, calculate error and save data
    puts("Solve equation, calculate error and save data");
    for (n = 2; n <= K; ++n)
    {
        int top = n % 2;

        //puts("Solve equation");
        Solve(y, g, tau, top);

        //puts("Calculate error");
        Get_Analytical_Solution(x, u, g, n, tau, U_basic);
        z_Ch[n] = Calculate_Error(u, y[top], z, g);

        //puts("Save data");
        Write_Silo(n, tau, x, g, y[top], u, z);
    }

    // Save error
    puts("Save error");
    b = malloc((K+1) * sizeof(*b));
    for (n = 0; n <= K; ++n)
    {
        b[n].m1 = b[n].m2 = sizeof(b[n].x) + sizeof(b[n].y);
        b[n].x = (float)n;
        printf("%e\n", z_Ch[n]);
    }
    Save_Graph(z_Ch, b, K, "z_Ch.bin", "wb");

    free(b);

    // Free memory
    //puts("Free memory");
    free(z_Ch);
    free(y[1]);
    free(y[0]);
    free(z);
    free(u);
    free(x[2]);
    free(x[1]);
    free(x[0]);

    return EXIT_SUCCESS;
}
// --------------------------------------------------------------------------
