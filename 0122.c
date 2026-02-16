#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ===========================================
   設定
   =========================================== */
#define MAX_ITER 600
#define EPS 1e-12
#define ALPHA 10.0
#define BETA 0.1
#define RHO 2
#define DT 0.1
#define External_Force 0.12

/* ノードの座標 */
typedef struct
{
    double x, y, z;
} Node;

typedef struct
{
    double x, y, z;
} MshNode;

/* ===========================================
   グラフ構造体
   =========================================== */
typedef struct
{
    int n;           /* ノード数 */
    int **adj;       /* 隣接行列 (n×n) */
    Node *nodes;     /* nodes[i].x, nodes[i].y, nodes[i].z */
    double **length; /* length[i][j] エッジ長。エッジ無しなら0 */
} Graph;

/* エラー処理 */
void die(const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    exit(EXIT_FAILURE);
}

FILE *xfopen(const char *fname, const char *mode)
{
    FILE *fp = fopen(fname, mode);
    if (!fp)
    {
        perror(fname);
        exit(EXIT_FAILURE);
    }
    return fp;
}

Graph create_graph_from_csv(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
        die("cannot open csv file");

    char line[1024];

    /* ----------------------------
       1. ノード読み込み
       ---------------------------- */
    int capacity = 128;
    int num_nodes = 0;

    int *node_ids = malloc(capacity * sizeof(int));
    Node *nodes_tmp = malloc(capacity * sizeof(Node));
    if (!node_ids || !nodes_tmp)
        die("memory allocation failed");

    /* "# ===== Nodes =====" を探す */
    while (fgets(line, sizeof(line), fp))
    {
        if (strstr(line, "Nodes") != NULL)
            break;
    }

    /* ヘッダ行 node_id,x,y,z */
    if (!fgets(line, sizeof(line), fp))
        die("node header missing");

    /* ノード本体 */
    while (fgets(line, sizeof(line), fp))
    {
        /* 空行 or コメントでノード部終了 */
        if (line[0] == '#' || strlen(line) < 3)
            break;

        int id;
        double x, y, z;

        if (sscanf(line, "%d,%lf,%lf,%lf", &id, &x, &y, &z) != 4)
            continue;

        if (num_nodes >= capacity)
        {
            capacity *= 2;
            node_ids = realloc(node_ids, capacity * sizeof(int));
            nodes_tmp = realloc(nodes_tmp, capacity * sizeof(Node));
            if (!node_ids || !nodes_tmp)
                die("realloc failed");
        }

        node_ids[num_nodes] = id;
        nodes_tmp[num_nodes].x = x;
        nodes_tmp[num_nodes].y = y;
        nodes_tmp[num_nodes].z = z;
        num_nodes++;
    }

    if (num_nodes == 0)
        die("no nodes found");

    /* ----------------------------
       2. ID → index 対応
       ---------------------------- */
    int max_id = 0;
    for (int i = 0; i < num_nodes; i++)
        if (node_ids[i] > max_id)
            max_id = node_ids[i];

    int *id_to_index = malloc((max_id + 1) * sizeof(int));
    if (!id_to_index)
        die("allocation failed");

    for (int i = 0; i <= max_id; i++)
        id_to_index[i] = -1;

    for (int i = 0; i < num_nodes; i++)
        id_to_index[node_ids[i]] = i;

    /* ----------------------------
       3. Graph 初期化
       ---------------------------- */
    Graph g;
    g.n = num_nodes;
    g.nodes = malloc(sizeof(Node) * g.n);
    g.adj = malloc(sizeof(int *) * g.n);
    g.length = malloc(sizeof(double *) * g.n);

    for (int i = 0; i < g.n; i++)
    {
        g.nodes[i] = nodes_tmp[i];
        g.adj[i] = calloc(g.n, sizeof(int));
        g.length[i] = calloc(g.n, sizeof(double));
    }

    /* ----------------------------
       4. エッジ読み込み
       ---------------------------- */
    /* エッジヘッダまでスキップ */
    while (fgets(line, sizeof(line), fp))
    {
        if (strstr(line, "edge_id") != NULL)
            break;
    }

    while (fgets(line, sizeof(line), fp))
    {
        int eid, n1, n2;
        double x1, y1, z1, x2, y2, z2;

        if (sscanf(line,
                   "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf",
                   &eid, &n1, &n2,
                   &x1, &y1, &z1,
                   &x2, &y2, &z2) != 9)
            continue;

        if (n1 > max_id || n2 > max_id)
            continue;

        int i = id_to_index[n1];
        int j = id_to_index[n2];
        if (i < 0 || j < 0)
            continue;

        g.adj[i][j] = g.adj[j][i] = 1;

        double dx = g.nodes[i].x - g.nodes[j].x;
        double dy = g.nodes[i].y - g.nodes[j].y;
        double dz = g.nodes[i].z - g.nodes[j].z;
        double L = sqrt(dx * dx + dy * dy + dz * dz);

        g.length[i][j] = g.length[j][i] = L;
    }

    fclose(fp);
    free(node_ids);
    free(nodes_tmp);
    free(id_to_index);

    return g;
}

void apply_external_force_to_top_nodes(
    const Graph *g,
    double *b,
    double force_z)
{
    if (!g || !b)
    {
        fprintf(stderr, "apply_external_force_to_top_nodes: NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    int n = g->n;

    /* -------------------------------
       z 最大値を探す
       ------------------------------- */
    double z_max = g->nodes[0].z;
    for (int i = 1; i < n; i++)
    {
        if (g->nodes[i].z > z_max)
            z_max = g->nodes[i].z;
    }

    /* 許容誤差（メッシュサイズ依存を避けるため相対） */
    double EPS_Z = 1e-8 * (fabs(z_max) + 1.0);

    /* -------------------------------
       上面ノードに外力を加える
       ------------------------------- */
    for (int i = 0; i < n; i++)
    {
        if (fabs(g->nodes[i].z - z_max) < EPS_Z)
        {
            int iz = 3 * i + 2;
            b[iz] += force_z;
        }
    }
}

void apply_external_force_to_top_nodes_xyz(
    const Graph *g,
    double *b,
    double fx,
    double fy,
    double fz)
{
    if (!g || !b)
        die("apply_external_force_to_top_nodes_xyz: NULL pointer");

    int n = g->n;

    /* z 最大値を探す */
    double z_max = g->nodes[0].z;
    for (int i = 1; i < n; i++)
        if (g->nodes[i].z > z_max)
            z_max = g->nodes[i].z;

    double EPS_Z = 1e-8 * (fabs(z_max) + 1.0);

    /* 上面ノードに外力付加 */
    for (int i = 0; i < n; i++)
    {
        if (fabs(g->nodes[i].z - z_max) < EPS_Z)
        {
            b[3 * i + 0] += fx;
            b[3 * i + 1] += fy;
            b[3 * i + 2] += fz;
        }
    }
}

/* 隣接行列の解放 */
void free_graph(Graph *g)
{
    for (int i = 0; i < g->n; i++)
    {
        free(g->adj[i]);
        free(g->length[i]);
    }
    free(g->adj);
    free(g->length);
    free(g->nodes);
}

/* === 分子計算関数の種類 === */
typedef double (*NumeratorFunc)(double s, double x, double y, double z);

static inline double num_A(double s, double x, double y, double z) { return s * x * x; }
static inline double num_B(double s, double x, double y, double z) { return s * x * y; }
static inline double num_C(double s, double x, double y, double z) { return s * x * z; }
static inline double num_D(double s, double x, double y, double z) { return s * y * y; }
static inline double num_E(double s, double x, double y, double z) { return s * y * z; }
static inline double num_F(double s, double x, double y, double z) { return s * z * z; }

/* === 汎用 Calc 関数 === */
static inline double generic_calc(double s,
                                  double x1, double x2,
                                  double y1, double y2,
                                  double z1, double z2,
                                  NumeratorFunc num_func)
{
    double x = x1 - x2;
    double y = y1 - y2;
    double z = z1 - z2;

    double sx = x * x;
    double sy = y * y;
    double sz = z * z;

    double denom2 = sx + sy + sz;
    if (denom2 < 1e-12)
        return 0.0;   

    double denom = pow(denom2, 1.5);
    return num_func(s, x, y, z) / denom;

}

double *compute_matrix(const Graph *g, double s, NumeratorFunc num_func)
{
    if (!g || g->n <= 0)
        return NULL;

    int n = g->n;
    double *M = calloc((size_t)n * n, sizeof(double));
    if (!M)
        die("calloc failed in compute_matrix");

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (g->adj[i][j] == 1)
            {

                double xi = g->nodes[i].x;
                double yi = g->nodes[i].y;
                double zi = g->nodes[i].z;
                double xj = g->nodes[j].x;
                double yj = g->nodes[j].y;
                double zj = g->nodes[j].z;

                double val = generic_calc(s, xi, xj, yi, yj, zi, zj, num_func);

                M[i * n + j] = val;
                M[j * n + i] = val;
            }
        }
        M[i * n + i] = 0.0;
    }

    return M;
}

double *compute_A_matrix(const Graph *g, double s) { return compute_matrix(g, s, num_A); }
double *compute_B_matrix(const Graph *g, double s) { return compute_matrix(g, s, num_B); }
double *compute_C_matrix(const Graph *g, double s) { return compute_matrix(g, s, num_C); }
double *compute_D_matrix(const Graph *g, double s) { return compute_matrix(g, s, num_D); }
double *compute_E_matrix(const Graph *g, double s) { return compute_matrix(g, s, num_E); }
double *compute_F_matrix(const Graph *g, double s) { return compute_matrix(g, s, num_F); }

/*２つのノード間の内力計算*/
void compute_f_ij(
    const Graph *g,
    const double *Acoef,
    const double *Bcoef,
    const double *Ccoef,
    const double *Dcoef,
    const double *Ecoef,
    const double *Fcoef,
    const double *v,
    const double *b,
    int i,
    int j,
    double f[3] /* f[0]=fx, f[1]=fy, f[2]=fz */
)
{
    if (!g || !Acoef || !Bcoef || !Ccoef ||
        !Dcoef || !Ecoef || !Fcoef || !v || !f)
    {
        fprintf(stderr, "compute_f_ij: NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    int n = g->n;

    if (i < 0 || i >= n || j < 0 || j >= n)
    {
        fprintf(stderr, "compute_f_ij: invalid index i=%d j=%d\n", i, j);
        exit(EXIT_FAILURE);
    }

    if (i == j || g->adj[i][j] == 0)
    {
        f[0] = f[1] = f[2] = 0.0;
        return;
    }

    int ix = 3 * i + 0;
    int iy = 3 * i + 1;
    int iz = 3 * i + 2;

    int jx = 3 * j + 0;
    int jy = 3 * j + 1;
    int jz = 3 * j + 2;

    double a = Acoef[i * n + j];
    double b_ = Bcoef[i * n + j];
    double c = Ccoef[i * n + j];
    double d = Dcoef[i * n + j];
    double e = Ecoef[i * n + j];
    double fzeta = Fcoef[i * n + j];

    /* fx */
    f[0] =
        -a * v[ix] - b_ * v[iy] - c * v[iz] + a * v[jx] + b_ * v[jy] + c * v[jz];

    /* fy */
    f[1] =
        -b_ * v[ix] - d * v[iy] - e * v[iz] + b_ * v[jx] + d * v[jy] + e * v[jz];

    /* fz */
    f[2] =
        -c * v[ix] - e * v[iy] - fzeta * v[iz] + c * v[jx] + e * v[jy] + fzeta * v[jz];
}

void compute_internal_force_i(
    const Graph *g,
    const double *Acoef,
    const double *Bcoef,
    const double *Ccoef,
    const double *Dcoef,
    const double *Ecoef,
    const double *Fcoef,
    const double *v,
    const double *b,
    int i,
    double f_i[3] /* 出力: f_i[x,y,z] */
)
{
    if (!g || !Acoef || !Bcoef || !Ccoef ||
        !Dcoef || !Ecoef || !Fcoef || !v || !f_i)
    {
        fprintf(stderr, "compute_internal_force_i: NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    int n = g->n;

    if (i < 0 || i >= n)
    {
        fprintf(stderr, "compute_internal_force_i: invalid node index %d\n", i);
        exit(EXIT_FAILURE);
    }

    /* 初期化 */
    f_i[0] = 0.0;
    f_i[1] = 0.0;
    f_i[2] = 0.0;

    /* 隣接ノードとの力を加算 */
    for (int j = 0; j < n; j++)
    {
        if (j == i)
            continue;
        if (g->adj[i][j] == 0)
            continue;

        double fij[3];
        compute_f_ij(g,
                     Acoef, Bcoef, Ccoef,
                     Dcoef, Ecoef, Fcoef,
                     v, b, i, j, fij);

        f_i[0] += fij[0];
        f_i[1] += fij[1];
        f_i[2] += fij[2];
    }
}

void compute_internal_edge_force(
    const Graph *g,
    const double *Acoef,
    const double *Bcoef,
    const double *Ccoef,
    const double *Dcoef,
    const double *Ecoef,
    const double *Fcoef,
    const double *v,
    int i,
    int j,
    double f_ij[3])
{
    /* 既存の内力計算をそのまま使用 */
    compute_f_ij(g,
                 Acoef, Bcoef, Ccoef,
                 Dcoef, Ecoef, Fcoef,
                 v, NULL, i, j, f_ij);
}

/*内力の大きさ|f|を計算する関数*/
double internal_force_magnitude(const double f_i[3])
{
    if (!f_i)
    {
        fprintf(stderr, "internal_force_magnitude: NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return sqrt(
        f_i[0] * f_i[0] +
        f_i[1] * f_i[1] +
        f_i[2] * f_i[2]);
}

int node_degree(const Graph *g, int i)
{
    int deg = 0;
    for (int j = 0; j < g->n; j++)
    {
        if (g->adj[i][j])
            deg++;
    }
    return deg;
}

void apply_external_force_to_edge_force(
    const Graph *g,
    const double *F_ext, /* 外力 [3n] */
    int i,
    int j,
    double f_ij[3] /* 入出力 */
)
{
    if (!g || !F_ext || !f_ij)
    {
        fprintf(stderr, "apply_external_force_to_edge_force: NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (i < 0 || i >= g->n || j < 0 || j >= g->n)
        return;

    if (g->adj[i][j] == 0 || i == j)
        return;

    int deg = node_degree(g, i);
    if (deg == 0)
        return;

    /* 外力成分 */
    double Fx = F_ext[3 * i + 0];
    double Fy = F_ext[3 * i + 1];
    double Fz = F_ext[3 * i + 2];

    f_ij[0] = Fx / deg;
    f_ij[1] = Fy / deg;
    f_ij[2] = Fz / deg;
}

/*エッジの成長式*/
double edge_evo(double in_force, double r)
{
    double evo_result;

    double x = fabs(in_force);
    evo_result = ((ALPHA * x) / (M_PI * pow(r, RHO))) - (BETA * r);

    return evo_result;
}

/*ノードiのrを使って成長させる*/
double compute_time_evo_i(
    const Graph *g,
    const double *Acoef,
    const double *Bcoef,
    const double *Ccoef,
    const double *Dcoef,
    const double *Ecoef,
    const double *Fcoef,
    const double *v,
    const double *b,
    int i,
    double r)
{
    double f_i[3];

    compute_internal_force_i(g,
                             Acoef, Bcoef, Ccoef,
                             Dcoef, Ecoef, Fcoef,
                             v, b, i, f_i);

    double force_mag = internal_force_magnitude(f_i);

    return edge_evo(force_mag, r);
}

/* ===========================================
   対角ブロックの和を計算
   =========================================== */
void compute_row_sums(const Graph *g,
                      const double *A, const double *B, const double *C,
                      const double *D, const double *E, const double *F,
                      int i, double sums[6])
{
    int n = g->n;
    for (int k = 0; k < 6; k++)
        sums[k] = 0.0;

    for (int j = 0; j < n; j++)
    {
        if (i == j)
            continue;
        if (g->adj[i][j] == 0)
            continue;

        int idx = i * n + j;

        sums[0] += A[idx];
        sums[1] += B[idx];
        sums[2] += C[idx];
        sums[3] += D[idx];
        sums[4] += E[idx];
        sums[5] += F[idx];
    }
}

/* ===========================================
   ブロック行列 M[n][n][3][3] を構築
   =========================================== */
double ****build_jacobi_block_matrix(
    const Graph *g,
    const double *A,
    const double *B,
    const double *C,
    const double *D,
    const double *E,
    const double *F)
{
    int n = g->n;

    double ****M = malloc(n * sizeof(double ***));
    for (int i = 0; i < n; ++i)
    {
        M[i] = malloc(n * sizeof(double **));
        for (int j = 0; j < n; ++j)
        {
            M[i][j] = malloc(3 * sizeof(double *));
            for (int a = 0; a < 3; ++a)
            {
                M[i][j][a] = malloc(3 * sizeof(double));
            }
        }
    }

    /* --- 非対角ブロック --- */
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {

            int idx = i * n + j;

            if (i != j)
            {
                if (g->adj[i][j] == 0)
                {
                    for (int a = 0; a < 3; a++)
                        for (int b = 0; b < 3; b++)
                            M[i][j][a][b] = 0.0;
                    continue;
                }

                M[i][j][0][0] = A[idx];
                M[i][j][0][1] = B[idx];
                M[i][j][0][2] = C[idx];
                M[i][j][1][0] = B[idx];
                M[i][j][1][1] = D[idx];
                M[i][j][1][2] = E[idx];
                M[i][j][2][0] = C[idx];
                M[i][j][2][1] = E[idx];
                M[i][j][2][2] = F[idx];
            }
        }
    }

    /* --- 対角ブロック（負の総和） --- */
    for (int i = 0; i < n; ++i)
    {

        double sums[6];
        compute_row_sums(g, A, B, C, D, E, F, i, sums);

        M[i][i][0][0] = -sums[0];
        M[i][i][0][1] = -sums[1];
        M[i][i][0][2] = -sums[2];
        M[i][i][1][0] = -sums[1];
        M[i][i][1][1] = -sums[3];
        M[i][i][1][2] = -sums[4];
        M[i][i][2][0] = -sums[2];
        M[i][i][2][1] = -sums[4];
        M[i][i][2][2] = -sums[5];
    }

    return M;
}

/* 解放 */
void free_jacobi_block_matrix(double ****M, int n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int a = 0; a < 3; ++a)
                free(M[i][j][a]);
            free(M[i][j]);
        }
        free(M[i]);
    }
    free(M);
}

double **alloc_matrix(int N)
{
    double **m = malloc(sizeof(double *) * N);
    for (int i = 0; i < N; i++)
        m[i] = malloc(sizeof(double) * N);
    return m;
}

void free_matrix(double **m, int N)
{
    for (int i = 0; i < N; i++)
        free(m[i]);
    free(m);
}

void flatten_block_matrix(int n, double ****M, double **A)
{
    int N = 3 * n;

    for (int bi = 0; bi < n; bi++)
    {
        for (int bj = 0; bj < n; bj++)
        {

            int row0 = bi * 3;
            int col0 = bj * 3;

            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    A[row0 + a][col0 + b] = M[bi][bj][a][b];
        }
    }
}

void zero_row(double **A, int N, int row)
{
    if (row < 0 || row >= N)
    {
        fprintf(stderr, "zero_row: invalid row index %d\n", row);
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < N; j++)
    {
        A[row][j] = 0.0;
    }
}

void zero_col(double **A, int N, int col)
{
    if (col < 0 || col >= N)
    {
        fprintf(stderr, "zero_col: invalid col index %d\n", col);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++)
    {
        A[i][col] = 0.0;
    }
}

void apply_dirichlet_bc(double **A, double *b, int N, int idx, double value)
{
    if (idx < 0 || idx >= N)
    {
        fprintf(stderr, "apply_dirichlet_bc: invalid index %d\n", idx);
        exit(EXIT_FAILURE);
    }

    /* 行を 0 */
    for (int j = 0; j < N; j++)
        A[idx][j] = 0.0;

    /* 列を 0 */
    for (int i = 0; i < N; i++)
        A[i][idx] = 0.0;

    /* 対角成分を 1 */
    A[idx][idx] = 1.0;

    /* 右辺を指定値に */
    b[idx] = value;
}

void apply_dirichlet_bc_to_bottom_nodes(
    const Graph *g,
    double **A,
    double *b,
    int N)
{
    if (!g || !A || !b)
        die("apply_dirichlet_bc_to_bottom_nodes: NULL pointer");

    int n = g->n;

    /* -------------------------------
       z 最小値を探す
       ------------------------------- */
    double z_min = g->nodes[0].z;
    for (int i = 1; i < n; i++)
        if (g->nodes[i].z < z_min)
            z_min = g->nodes[i].z;

    /* 相対誤差 */
    double EPS_Z = 1e-8 * (fabs(z_min) + 1.0);

    /* -------------------------------
       bottom ノードに Dirichlet BC
       ------------------------------- */
    for (int i = 0; i < n; i++)
    {
        if (fabs(g->nodes[i].z - z_min) < EPS_Z)
        {
            /* x, y, z を全て固定 */
            apply_dirichlet_bc(A, b, N, 3 * i + 0, 0.0);
            apply_dirichlet_bc(A, b, N, 3 * i + 1, 0.0);
            apply_dirichlet_bc(A, b, N, 3 * i + 2, 0.0);
        }
    }
}

void gaussian_elimination(int N, double **A, double *b, double *x)
{
    /* 拡大係数行列を作る */
    double **M = alloc_matrix(N);
    double *bb = malloc(sizeof(double) * N);

    for (int i = 0; i < N; i++)
    {
        bb[i] = b[i];
        for (int j = 0; j < N; j++)
            M[i][j] = A[i][j];
    }

    /* 前進消去（部分ピボット） */
    for (int k = 0; k < N; k++)
    {

        /* ピボット選択 */
        int pivot = k;
        double max = fabs(M[k][k]);
        for (int i = k + 1; i < N; i++)
        {
            if (fabs(M[i][k]) > max)
            {
                max = fabs(M[i][k]);
                pivot = i;
            }
        }

        if (max < 1e-14)
        {
            M[k][k] = (M[k][k] >= 0 ? 1 : -1) * 1e-14;
        }

        /* 行交換 */
        if (pivot != k)
        {
            double *tmp = M[k];
            M[k] = M[pivot];
            M[pivot] = tmp;

            double tb = bb[k];
            bb[k] = bb[pivot];
            bb[pivot] = tb;
        }

        /* 消去 */
        for (int i = k + 1; i < N; i++)
        {
            double factor = M[i][k] / M[k][k];
            for (int j = k; j < N; j++)
                M[i][j] -= factor * M[k][j];
            bb[i] -= factor * bb[k];
        }
    }

    /* 後退代入 */
    for (int i = N - 1; i >= 0; i--)
    {
        double sum = bb[i];
        for (int j = i + 1; j < N; j++)
            sum -= M[i][j] * x[j];
        x[i] = sum / M[i][i];
    }

    free_matrix(M, N);
    free(bb);
}

/*ここからcsvファイル出力用関数*/
void write_v_csv(const char *filename, const double *v, int N)
{
    FILE *fp = xfopen(filename, "w");

    fprintf(fp, "dof,value\n");
    for (int i = 0; i < N; i++)
    {
        fprintf(fp, "%d,%.10f\n", i, v[i]);
    }

    fclose(fp);
}

void write_edge_radius_csv(
    const char *filename,
    const Graph *g,
    double r)
{
    FILE *fp = xfopen(filename, "w");

    fprintf(fp, "node_i,node_j,r\n");

    for (int i = 0; i < g->n; i++)
    {
        for (int j = i + 1; j < g->n; j++)
        {
            if (g->adj[i][j] == 0)
                continue;

            fprintf(fp, "%d,%d,%.10f\n", i, j, r);
        }
    }

    fclose(fp);
}

void write_edge_force_csv(
    const char *filename,
    const Graph *g,
    const double *Acoef,
    const double *Bcoef,
    const double *Ccoef,
    const double *Dcoef,
    const double *Ecoef,
    const double *Fcoef,
    const double *v,
    const double *b)
{
    FILE *fp = xfopen(filename, "w");

    fprintf(fp, "node_i,node_j,fx,fy,fz,fmag\n");

    for (int i = 0; i < g->n; i++)
    {
        for (int j = i + 1; j < g->n; j++)
        {
            if (g->adj[i][j] == 0)
                continue;

            double f_ij[3];
            compute_f_ij(g,
                         Acoef, Bcoef, Ccoef,
                         Dcoef, Ecoef, Fcoef,
                         v, b, i, j, f_ij);

            double fmag = internal_force_magnitude(f_ij);

            fprintf(fp, "%d,%d,%.10f,%.10f,%.10f,%.10f\n",
                    i, j,
                    f_ij[0], f_ij[1], f_ij[2],
                    fmag);
        }
    }

    fclose(fp);
}

double **alloc_edge_radius(const Graph *g, double r0)
{
    int n = g->n;
    double **r = malloc(n * sizeof(double *));
    if (!r)
        die("alloc_edge_radius failed");

    for (int i = 0; i < n; i++)
    {
        r[i] = malloc(n * sizeof(double));
        if (!r[i])
            die("alloc_edge_radius failed");

        for (int j = 0; j < n; j++)
        {
            r[i][j] = (g->adj[i][j]) ? r0 : 0.0;
        }
    }
    return r;
}

void free_edge_radius(double **r, int n)
{
    for (int i = 0; i < n; i++)
        free(r[i]);
    free(r);
}

void update_edge_radius(
    const Graph *g,
    const double *Acoef,
    const double *Bcoef,
    const double *Ccoef,
    const double *Dcoef,
    const double *Ecoef,
    const double *Fcoef,
    const double *v,
    double **edge_r,
    double *max_dr 
)
{
    int n = g->n;
    *max_dr = 0.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {

            if (!g->adj[i][j])
                continue;

            double f_ij[3];
            compute_internal_edge_force(
                g,
                Acoef, Bcoef, Ccoef,
                Dcoef, Ecoef, Fcoef,
                v,
                i, j,
                f_ij);

            double fmag = sqrt(
                f_ij[0] * f_ij[0] +
                f_ij[1] * f_ij[1] +
                f_ij[2] * f_ij[2]);

            double r_old = edge_r[i][j];

            double dr =
                (ALPHA * fabs(fmag)) /
                    (M_PI * pow(r_old, RHO)) -
                (BETA * r_old);

            double r_new = r_old + DT * dr;

            if (r_new < 1e-6)
                r_new = 1e-6;

            double diff = fabs(r_new - r_old);
            if (diff > *max_dr)
                *max_dr = diff;

            edge_r[i][j] = r_new;
            edge_r[j][i] = r_new;
        }
    }
}

void write_edge_radius_csv_updated(
    const char *filename,
    const Graph *g,
    double **edge_r)
{
    FILE *fp = xfopen(filename, "w");
    fprintf(fp, "node_i,node_j,r\n");

    for (int i = 0; i < g->n; i++)
    {
        for (int j = i + 1; j < g->n; j++)
        {
            if (!g->adj[i][j])
                continue;
            fprintf(fp, "%d,%d,%.10f\n",
                    i, j, edge_r[i][j]);
        }
    }
    fclose(fp);
}

void update_stiffness_by_radius(
    const Graph *g,
    const double *A0, const double *B0, const double *C0,
    const double *D0, const double *E0, const double *F0,
    double **edge_r,
    double r0,
    double *A, double *B, double *C,
    double *D, double *E, double *F)
{
    int n = g->n;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {

            int idx = i * n + j;

            if (!g->adj[i][j])
            {
                A[idx] = B[idx] = C[idx] =
                    D[idx] = E[idx] = F[idx] = 0.0;
                continue;
            }

            double scale = pow(edge_r[i][j] / r0, 2.0);

            A[idx] = A0[idx] * scale;
            B[idx] = B0[idx] * scale;
            C[idx] = C0[idx] * scale;
            D[idx] = D0[idx] * scale;
            E[idx] = E0[idx] * scale;
            F[idx] = F0[idx] * scale;
        }
    }
}

void write_deformed_nodes_and_edges_csv(
    const char *filename,
    const Graph *g,
    const double *v,
    double **edge_r)
{
    if (!g || !v || !edge_r)
        die("write_deformed_nodes_and_edges_csv: NULL pointer");

    FILE *fp = xfopen(filename, "w");

    /* 添付CSVと完全一致するヘッダ */
    fprintf(fp, "type,id,i,j,x,y,z,r\n");

    /* ======================
       ノード出力
       ====================== */
    for (int i = 0; i < g->n; i++)
    {
        double x_new = g->nodes[i].x - v[3 * i + 0];
        double y_new = g->nodes[i].y - v[3 * i + 1];
        double z_new = g->nodes[i].z - v[3 * i + 2];

        fprintf(fp,
                "node,%d,,,% .10f,% .10f,% .10f,\n",
                i, x_new, y_new, z_new);
    }

    /* ======================
       エッジ出力
       ====================== */
    int edge_id = 0;

    for (int i = 0; i < g->n; i++)
    {
        for (int j = i + 1; j < g->n; j++)
        {
            if (!g->adj[i][j])
                continue;

            /* ノード i の変形後座標 */
            double xi = g->nodes[i].x - v[3 * i + 0];
            double yi = g->nodes[i].y - v[3 * i + 1];
            double zi = g->nodes[i].z - v[3 * i + 2];

            fprintf(fp,
                    "edge,%d,%d,%d,% .10f,% .10f,% .10f,% .10f\n",
                    edge_id,
                    i, j,
                    xi, yi, zi,
                    edge_r[i][j]);

            edge_id++;
        }
    }

    fclose(fp);
}

void write_mesh_nodes_edges_with_radius_csv(
    const char *filename,
    const Graph *g,
    const double *v,
    double **edge_r)
{
    if (!g || !v || !edge_r)
        die("write_mesh_nodes_edges_with_radius_csv: NULL pointer");

    FILE *fp = xfopen(filename, "w");

    /* ======================
       Nodes
       ====================== */
    fprintf(fp, "# ===== Nodes =====\n");
    fprintf(fp, "node_id,x,y,z\n");

    for (int i = 0; i < g->n; i++)
    {
        double x_new = g->nodes[i].x - v[3 * i + 0];
        double y_new = g->nodes[i].y - v[3 * i + 1];
        double z_new = g->nodes[i].z - v[3 * i + 2];

        fprintf(fp, "%d,%.10f,%.10f,%.10f\n",
                i, x_new, y_new, z_new);
    }

    /* ======================
       Edges (+ radius r)
       ====================== */
    fprintf(fp, "\n# ===== Edges =====\n");
    fprintf(fp,
            "edge_id,n1_id,n2_id,x1,y1,z1,x2,y2,z2,r\n");

    int edge_id = 0;

    for (int i = 0; i < g->n; i++)
    {
        for (int j = i + 1; j < g->n; j++)
        {
            if (!g->adj[i][j])
                continue;

            double x1 = g->nodes[i].x - v[3 * i + 0];
            double y1 = g->nodes[i].y - v[3 * i + 1];
            double z1 = g->nodes[i].z - v[3 * i + 2];

            double x2 = g->nodes[j].x - v[3 * j + 0];
            double y2 = g->nodes[j].y - v[3 * j + 1];
            double z2 = g->nodes[j].z - v[3 * j + 2];

            fprintf(fp,
                    "%d,%d,%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                    edge_id,
                    i, j,
                    x1, y1, z1,
                    x2, y2, z2,
                    edge_r[i][j]);

            edge_id++;
        }
    }

    fclose(fp);
}

void write_deformed_snapshot_every_100(
    int step,
    const Graph *g,
    const double *v,
    double **edge_r)
{
    char filename[256];

    /* 100ステップごとだけ出力 */
    if (step % 100 != 0)
        return;

    snprintf(
        filename,
        sizeof(filename),
        "mesh_3_deformed_step_%04d.csv",
        step);

    write_mesh_nodes_edges_with_radius_csv(
        filename,
        g,
        v,
        edge_r);

    printf("  Snapshot written: %s\n", filename);
}


/* ===========================================
   main
   =========================================== */
int main(void)
{
    double s_young = 0.22 * M_PI;
    double r0 = 0.5;
    int nstep = MAX_ITER;
    double tol = 1e-8;

    double *v_final = NULL;

    Graph g = create_graph_from_csv("mesh.csv");
    int n = g.n;
    int N = 3 * n;

    FILE *fp_vhist = xfopen("v_history.csv", "w");
    FILE *fp_rhist = xfopen("edge_radius_history.csv", "w");

    double *A0 = compute_A_matrix(&g, s_young);
    double *B0 = compute_B_matrix(&g, s_young);
    double *C0 = compute_C_matrix(&g, s_young);
    double *D0 = compute_D_matrix(&g, s_young);
    double *E0 = compute_E_matrix(&g, s_young);
    double *F0 = compute_F_matrix(&g, s_young);

    double *Acoef = calloc(n * n, sizeof(double));
    double *Bcoef = calloc(n * n, sizeof(double));
    double *Ccoef = calloc(n * n, sizeof(double));
    double *Dcoef = calloc(n * n, sizeof(double));
    double *Ecoef = calloc(n * n, sizeof(double));
    double *Fcoef = calloc(n * n, sizeof(double));

    double **edge_r = alloc_edge_radius(&g, r0);

    int **edge_id = malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++)
    {
        edge_id[i] = malloc(n * sizeof(int));
        for (int j = 0; j < n; j++)
            edge_id[i][j] = -1;
    }

    /* edge_id 割り当て */
    int eid = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (g.adj[i][j])
            {
                edge_id[i][j] = edge_id[j][i] = eid++;
            }
        }
    }
    int num_edges = eid;

    /* ヘッダ */
    fprintf(fp_vhist, "step");
    for (int i = 0; i < N; i++)
        fprintf(fp_vhist, ",v%d", i);
    fprintf(fp_vhist, "\n");

    /* ヘッダ */
    fprintf(fp_rhist, "step");
    for (int e = 0; e < num_edges; e++)
        fprintf(fp_rhist, ",e%d", e);
    fprintf(fp_rhist, "\n");

    for (int step = 0; step < nstep; step++)
    {

        printf("time step %d\n", step);

        /* (1) 半径 → 剛性 */
        update_stiffness_by_radius(
            &g,
            A0, B0, C0,
            D0, E0, F0,
            edge_r, r0,
            Acoef, Bcoef, Ccoef,
            Dcoef, Ecoef, Fcoef);

        /* (2) 行列構築 */
        double ****Mblock =
            build_jacobi_block_matrix(
                &g,
                Acoef, Bcoef, Ccoef,
                Dcoef, Ecoef, Fcoef);

        double **A = alloc_matrix(N);
        flatten_block_matrix(n, Mblock, A);

        double *b = calloc(N, sizeof(double));
        double *v = calloc(N, sizeof(double));

        /* z 方向の外力（常に） */
apply_external_force_to_top_nodes(&g, b, External_Force);

/* step が 5 の倍数 → x 方向外力 */
if (step % 5 == 0)
{
    apply_external_force_to_top_nodes_xyz(
        &g, b,
        0.06,   /* fx */
        0.0,    /* fy */
        0.0     /* fz */
    );
}

/* step が 7 の倍数 → y 方向外力 */
if (step % 7 == 0)
{
    apply_external_force_to_top_nodes_xyz(
        &g, b,
        0.0,    /* fx */
        0.06,   /* fy */
        0.0     /* fz */
    );
}
        /* bottom ノードを完全固定 */
        apply_dirichlet_bc_to_bottom_nodes(&g, A, b, N);

        /* (3) 変位 */
        gaussian_elimination(N, A, b, v);

        write_deformed_snapshot_every_100(
    step,
    &g,
    v,
    edge_r);

        /* 各ステップ */
        fprintf(fp_vhist, "%d", step);
        for (int i = 0; i < N; i++)
            fprintf(fp_vhist, ",%.10f", v[i]);
        fprintf(fp_vhist, "\n");

        /* (4) 半径更新＋収束量 */
        double max_dr;
        update_edge_radius(
            &g,
            Acoef, Bcoef, Ccoef,
            Dcoef, Ecoef, Fcoef,
            v,
            edge_r,
            &max_dr);

        /* 各ステップ */
        fprintf(fp_rhist, "%d", step);

        double *rbuf = calloc(num_edges, sizeof(double));

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (!g.adj[i][j])
                    continue;
                int eid = edge_id[i][j];
                rbuf[eid] = edge_r[i][j];
            }
        }

        for (int e = 0; e < num_edges; e++)
            fprintf(fp_rhist, ",%.10f", rbuf[e]);

        fprintf(fp_rhist, "\n");
        free(rbuf);

        if (step == nstep - 1 || max_dr < tol)
        {
            v_final = malloc(sizeof(double) * N);
            if (!v_final)
                die("malloc failed for v_final");
            memcpy(v_final, v, sizeof(double) * N);
        }

        if (max_dr < tol)
        {
            printf("Converged at step %d\n", step);
            free_jacobi_block_matrix(Mblock, n);
            free_matrix(A, N);
            free(b);
            free(v);
            break;
        }

        /* 後始末 */
        free_jacobi_block_matrix(Mblock, n);
        free_matrix(A, N);
        free(b);
        free(v);
    }

    write_edge_radius_csv_updated("edge_radius_final.csv", &g, edge_r);

    if (v_final)
    {
        write_mesh_nodes_edges_with_radius_csv(
    "mesh_nodes_edges_3_deformed.csv",
    &g,
    v_final,
    edge_r);
    }

    FILE *gp = popen("gnuplot", "w");
    if (!gp)
        die("gnuplot failed");

    fprintf(gp,
            "set terminal pngcairo size 1200,800\n"
            "set output 'v_history.png'\n"
            "set datafile separator ','\n"
            "set xlabel 'step'\n"
            "set ylabel 'displacement'\n"
            "set key off\n"
            "plot for [i=2:%d] 'v_history.csv' using 1:i with lines\n",
            N + 1);

    pclose(gp);

    gp = popen("gnuplot", "w");
    if (!gp)
        die("gnuplot failed");

    fprintf(gp,
            "set terminal pngcairo size 1200,800\n"
            "set output 'edge_radius_history.png'\n"
            "set datafile separator ','\n"
            "set xlabel 'step'\n"
            "set ylabel 'radius'\n"
            "set key off\n"
            "plot for [i=2:%d] 'edge_radius_history.csv' using 1:i with lines\n",
            num_edges + 1);

    pclose(gp);

    /* 解放 */
    free_edge_radius(edge_r, n);
    free(A0);
    free(B0);
    free(C0);
    free(D0);
    free(E0);
    free(F0);
    free(Acoef);
    free(Bcoef);
    free(Ccoef);
    free(Dcoef);
    free(Ecoef);
    free(Fcoef);
    free_graph(&g);

    fclose(fp_vhist);
    fclose(fp_rhist);

    for (int i = 0; i < n; i++)
        free(edge_id[i]);
    free(edge_id);

    return 0;
}