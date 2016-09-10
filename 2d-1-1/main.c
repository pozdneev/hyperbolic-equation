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
    int Ns0, Ns1;   // # of grid segments in each spatial dir'n
    double L0, L1;  // # of internal grid points in each spatial dir'n
    double h0, h1;
} Grid_Info;
// --------------------------------------------------------------------------
/*
static void Print(double *y[2], int layer, Grid_Info *g)
{
    int i0, i1;

    for (i0 = 0; i0 <= g->Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g->Ns1; ++i1)
        {
            int idx = IDX2(i0, i1, 0, 0, g->Ns0+1, g->Ns1+1);
            printf("%+9.2e ", y[layer][idx]);
        }
        putchar('\n');
    }
}
*/
// --------------------------------------------------------------------------
/*
 * NB: In multi-multi-multidimensional case we should think carefully about
 * type of storage for point data. Should we use smth like
 *      typedef Point_s {
 *          double x[Ndims];
 *      } Point;
 *      Point x[Np];
 * or smth like
 *      typedef Point_s {
 *          double x[Np];
 *      } Point;
 *      Point x[Ndims];
 * ? 
 */
/*
typedef Point_s
{
    double x0, x1;
} Point;
*/
// --------------------------------------------------------------------------
typedef double fun_t(double, double, double, Grid_Info);
// --------------------------------------------------------------------------
double U_basic(double x0, double x1, double t, Grid_Info g)
{
    assert(t >= 0.0);
    assert(g.L0 > 0.0);
    assert(g.L1 > 0.0);

    return sin(PI * x0 / g.L0) * sin(PI * x1 / g.L1) * 
        cos(PI * t * sqrt(1 / SQR(g.L0) + 1 / SQR(g.L1)));
//    if (x0 >= 0.3 && x0 <= 0.4 && x1 >= 0.4 && x1 <= 0.5) return 1; else return 0;
//    return exp(-(SQR(x0 - 0.35) + SQR(x1 - 0.45)) / 0.01);
}
// --------------------------------------------------------------------------
/*
static void Init_Buffer(const double *x, XDraw_Record *b, int N)
{
    assert(x);
    assert(b);
    assert(N > 0);

    int i;

    for (i = 0; i <= N; ++i)
    {
        int idx = IDX1(i, 0, N+1);
        b[idx].m1 = b[idx].m2 = sizeof(b[idx].x) + sizeof(b[idx].y);
        b[idx].x = (float)x[idx];
    }
}
*/
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
    /*
    assert(g.L0 > 0.0);
    assert(g.L1 > 0.0);
    assert(g.h0 > 0.0);
    assert(g.h1 > 0.0);
    */

    int i0, i1;

    double z_Ch = 0.0;

    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g.Ns1; ++i1)
        {
            int idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);

            z[idx] = y[idx] - u[idx];

            /*
            if (fabs(u[idx]) > u_Ch) u_Ch = fabs(u[idx]);
            if (fabs(y[idx]) > y_Ch) y_Ch = fabs(y[idx]);
            */
            if (fabs(z[idx]) > z_Ch) z_Ch = fabs(z[idx]);

            /*
            u_L1h += SQR(u[idx]);
            y_L1h += SQR(y[idx]);
            z_L1h += SQR(z[idx]);
            */
        }
    }

    /*
    u_L1h = sqrt(u_L1h * h);
    y_L1h = sqrt(y_L1h * h);
    z_L1h = sqrt(z_L1h * h);

    fprintf(fout, "||u||_{C^h} = %e\n", u_Ch);
    fprintf(fout, "||y||_{C^h} = %e\n", y_Ch);
    fprintf(fout, "||z||_{C^h} = %e\n", z_Ch);

    fprintf(fout, "||u||_{L_2^h} = %e\n", u_L1h);
    fprintf(fout, "||y||_{L_2^h} = %e\n", y_L1h);
    fprintf(fout, "||z||_{L_2^h} = %e\n", z_L1h);
    */

    return z_Ch;
}
// --------------------------------------------------------------------------
static void Prepare_Initial_Data(
        double *x[], double *y[], Grid_Info g, double tau, fun_t fun)
{
    assert(x);
    assert(x[0]);
    assert(x[1]);
    assert(y);
    assert(y[0]);
    assert(y[1]);
    assert(g.Ns0 > 0);
    assert(g.Ns1 > 0);
    assert(g.h0 > 0.0);
    assert(g.h1 > 0.0);
    assert(tau > 0.0);
    assert(fun);

    double gamma0 = SQR(tau / g.h0);
    double gamma1 = SQR(tau / g.h1);
    const int mid = 1;
    const int bot = 0;
    double t = 0.0;

    int i0, i1, idx;

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

    // t^0
    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g.Ns1; ++i1)
        {
            idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);

            y[bot][idx] = fun(x[0][i0], x[1][i1], t, g);
        }
    }

    /*
    puts("=== 00 ===");
    Print(y, bot, &g);
    */

    // t^1
    for (i0 = 1; i0 <= g.Ns0-1; ++i0)
    {
        for (i1 = 1; i1 <= g.Ns1-1; ++i1)
        {
            int idx0l, idx0r, idx1l, idx1r;

            idx   = IDX2(i0, i1,   0, 0, g.Ns0+1, g.Ns1+1);
            idx0l = IDX2(i0-1, i1, 0, 0, g.Ns0+1, g.Ns1+1);
            idx0r = IDX2(i0+1, i1, 0, 0, g.Ns0+1, g.Ns1+1);
            idx1l = IDX2(i0, i1-1, 0, 0, g.Ns0+1, g.Ns1+1);
            idx1r = IDX2(i0, i1+1, 0, 0, g.Ns0+1, g.Ns1+1);

            y[mid][idx] = y[bot][idx] +
                gamma0 / 2 * (y[bot][idx0l] - 2 * y[bot][idx] + y[bot][idx0r]) +
                gamma1 / 2 * (y[bot][idx1l] - 2 * y[bot][idx] + y[bot][idx1r]);
        }
    }
    
    // Boundary points
    for (i1 = 0; i1 <= g.Ns1; ++i1)
    {
        i0 = 0;
        idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
        y[mid][idx] = y[bot][idx];

        i0 = g.Ns0;
        idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
        y[mid][idx] = y[bot][idx];
    }

    for (i0 = 1; i0 <= g.Ns0-1; ++i0)
    {
        i1 = 0;
        idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
        y[mid][idx] = y[bot][idx];

        i1 = g.Ns1;
        idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
        y[mid][idx] = y[bot][idx];
    }

    /*
    puts("=== 01 ===");
    Print(y, mid, &g);
    */
}
// --------------------------------------------------------------------------
static void Get_Analytical_Solution(
        double *x[], double *u, Grid_Info g, int n, double tau, fun_t fun)
{
    assert(x);
    assert(x[0]);
    assert(x[1]);
    assert(u);
    assert(g.Ns0 > 0);
    assert(g.Ns1 > 0);
    assert(n >= 0);
    assert(tau > 0.0);
    assert(fun);

    double t = tau * n;
    int i0, i1;

    // t^n
    for (i0 = 0; i0 <= g.Ns0; ++i0)
    {
        for (i1 = 0; i1 <= g.Ns1; ++i1)
        {
            int idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
            u[idx] = fun(x[0][i0], x[1][i1], t, g);
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
    assert(g.h0 > 0);
    assert(g.h1 > 0);
    assert(tau > 0.0);
    assert(top == 0 || top == 1);

    double gamma0 = SQR(tau / g.h0);
    double gamma1 = SQR(tau / g.h1);

    int mid = (top + 1) % 2;
    int bot = (mid + 1) % 2;

    int i0, i1, idx;

    // Internal points
    for (i0 = 1; i0 <= g.Ns0-1; ++i0)
    {
        for (i1 = 1; i1 <= g.Ns1-1; ++i1)
        {
            int idx0l = IDX2(i0-1, i1, 0, 0, g.Ns0+1, g.Ns1+1);
            int idx0r = IDX2(i0+1, i1, 0, 0, g.Ns0+1, g.Ns1+1);
            int idx1l = IDX2(i0, i1-1, 0, 0, g.Ns0+1, g.Ns1+1);
            int idx1r = IDX2(i0, i1+1, 0, 0, g.Ns0+1, g.Ns1+1);
                idx   = IDX2(i0, i1,   0, 0, g.Ns0+1, g.Ns1+1);

            y[top][idx] = 2 * y[mid][idx] - y[bot][idx] +
                gamma0 * (y[mid][idx0l] - 2 * y[mid][idx] + y[mid][idx0r]) +
                gamma1 * (y[mid][idx1l] - 2 * y[mid][idx] + y[mid][idx1r]);
        }
    }

    // Boundary points
    for (i1 = 0; i1 <= g.Ns1; ++i1)
    {
        i0 = 0;
        idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
        y[top][idx] = y[mid][idx];

        i0 = g.Ns0;
        idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
        y[top][idx] = y[mid][idx];
    }

    for (i0 = 1; i0 <= g.Ns0-1; ++i0)
    {
        i1 = 0;
        idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
        y[top][idx] = y[mid][idx];

        i1 = g.Ns1;
        idx = IDX2(i0, i1, 0, 0, g.Ns0+1, g.Ns1+1);
        y[top][idx] = y[mid][idx];
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
    assert(g.Ns0 > 0);
    assert(g.Ns1 > 0);
    assert(y);
    assert(u);
    assert(z);

    char *fname;
    const char *pattern = "%04d.silo";

    double dtime = n * tau;

    DBfile *dbfile = NULL;
    DBoptlist *optlist = DBMakeOptlist(2);
    int dims[] = { g.Ns0+1, g.Ns1+1 };
    int ndims = 2;
    
    // Create file
    // TODO: replace magic `4' with max number of decimal digits in n
    fname = (char *)malloc((4 + 5 + 1) * sizeof(*fname));
    sprintf(fname, pattern, n);
    dbfile = DBCreate(fname, DB_CLOBBER, DB_LOCAL,
            "Comment about the data", DB_PDB); // DB_PDB, DB_HDF5

    if(dbfile == NULL)
    {
        fprintf(stderr, "Could not create Silo file!\n");
        exit(EXIT_FAILURE);
    }

    free(fname);

    // Save cycle and time values to an option list
    DBAddOption(optlist, DBOPT_DTIME, &dtime);
    DBAddOption(optlist, DBOPT_CYCLE, &n);

    // Write a rectilinear mesh
    DBPutQuadmesh(dbfile, "quadmesh", NULL, x, dims, ndims,
                DB_DOUBLE, DB_COLLINEAR, optlist);

    DBFreeOptlist(optlist);

    // Write scalar values
    DBPutQuadvar1(dbfile, "y", "quadmesh", y, dims,
        ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
    DBPutQuadvar1(dbfile, "u", "quadmesh", u, dims,
        ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
    DBPutQuadvar1(dbfile, "z", "quadmesh", z, dims,
        ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);

    DBClose(dbfile);
}
// --------------------------------------------------------------------------
int main()
{
    double *x[2];   // Coordinates
    double *u;      // Analytical solution at a particular time-step
    double *z;      // Error at each point at a particular time-step
    double *y[2];   // Numerical solution at two consecutive time-steps
    double *z_Ch;   // Error

    XDraw_Record *b;    // Buffer to save data in XDraw format

    Grid_Info g = { 10, 10, 1.0, 1.0, -1, -1 };

    const double T = 1*g.L0;
    const int K = 2*g.Ns0;
    const double tau = T / (1*K);   // TODO: what is stability condition?

    int n;  // Counter through time-steps
    const char *mode = NULL;

    g.h0 = g.L0 / g.Ns0;
    g.h1 = g.L1 / g.Ns1;

    // Allocate memory
//    puts("Allocate memory");
    x[0] = malloc((g.Ns0+1) * sizeof(*(x[0])));
    x[1] = malloc((g.Ns1+1) * sizeof(*(x[1])));
    u = malloc((g.Ns0+1) * (g.Ns1+1) * sizeof(*u));
    z = malloc((g.Ns0+1) * (g.Ns1+1) * sizeof(*z));
    y[0] = malloc((g.Ns0+1) * (g.Ns1+1) * sizeof(*(y[0])));
    y[1] = malloc((g.Ns0+1) * (g.Ns1+1) * sizeof(*(y[1])));
    z_Ch = malloc((K+1) * sizeof(*z_Ch));

    // Prepare initial data
//    puts("Prepare initial data");
    Prepare_Initial_Data(x, y, g, tau, U_basic);
    
    /*
    // Prepare buffer
    b = malloc((N+1) * sizeof(*b));
    Init_Buffer(x, b, N);
    */

    // t^0
//    puts("t^0");
    n = 0;
    mode = "wb";
    Get_Analytical_Solution(x, u, g, n, tau, U_basic);
    z_Ch[n] = Calculate_Error(u, y[n], z, g);
    /*
    Save_Graph(u, b, N, "u.bin", mode);
    Save_Graph(y[n], b, N, "y.bin", mode);
    Save_Graph(z, b, N, "z.bin", mode);
    */
    Write_Silo(n, tau, x, g, y[n], u, z);

    // t^1
//    puts("t^1");
    n = 1;
    mode = "ab";
    Get_Analytical_Solution(x, u, g, n, tau, U_basic);
    z_Ch[n] = Calculate_Error(u, y[n], z, g);
    /*
    Save_Graph(u, b, N, "u.bin", mode);
    Save_Graph(y[n], b, N, "y.bin", mode);
    Save_Graph(z, b, N, "z.bin", mode);
    */
    Write_Silo(n, tau, x, g, y[n], u, z);

    // Solve equation, calculate error and save data
//    puts("Solve equation, calculate error and save data");
    for (n = 2; n <= K; ++n)
    {
        int top = n % 2;

        //puts("Solve equation");
        Solve(y, g, tau, top);

        /*
        printf("=== %02d ===\n", n);
        Print(y, top, &g);
        */

        //puts("Calculate error");
        Get_Analytical_Solution(x, u, g, n, tau, U_basic);
        z_Ch[n] = Calculate_Error(u, y[top], z, g);
//        printf("%e\n", z_Ch[n]);

        //puts("Save data");
        /*
        Save_Graph(u, b, N, "u.bin", "ab");
        Save_Graph(y[top], b, N, "y.bin", "ab");
        Save_Graph(z, b, N, "z.bin", "ab");
        */
        Write_Silo(n, tau, x, g, y[top], u, z);
    }

    /*
    free(b);
    */

    // Save error
//    puts("Save error");
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
//    puts("Free memory");
    free(z_Ch);
    free(y[1]);
    free(y[0]);
    free(z);
    free(u);
    free(x[1]);
    free(x[0]);

    return EXIT_SUCCESS;
}
// --------------------------------------------------------------------------
