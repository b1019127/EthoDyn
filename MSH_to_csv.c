#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LINE_LEN 1024
#define EDGE_LENGTH_LIMIT 5.0
#define MAX_EDGES 2000000

typedef struct {
    double x, y, z;
} Node;

typedef struct {
    int n1, n2;
} Edge;

/* ================= 判定関数 ================= */

double edge_length(int n1, int n2, Node *nodes)
{
    double dx = nodes[n2].x - nodes[n1].x;
    double dy = nodes[n2].y - nodes[n1].y;
    double dz = nodes[n2].z - nodes[n1].z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

int edge_exists(Edge *edges, int count, int a, int b)
{
    for (int i = 0; i < count; i++) {
        if (edges[i].n1 == a && edges[i].n2 == b)
            return 1;
    }
    return 0;
}

int is_removed_node(int nid, Node *nodes)
{
    if (nodes[nid].z == 0.0 &&
        nodes[nid].x >= -1.3 &&
        nodes[nid].x <=  1.3) {
        return 1;
    }
    return 0;
}

/* ================= main ================= */

int main(void)
{
    FILE *fp = fopen("bridge.msh", "r");
    FILE *out = fopen("mesh.csv", "w");
    if (!fp || !out) return 1;

    char line[LINE_LEN];
    int maxNodeTag = 0;
    Node *nodes = NULL;

    /* ---------- Nodes ---------- */
    while (fgets(line, LINE_LEN, fp)) {
        if (strncmp(line, "$Nodes", 6) == 0) {
            int nb, nn, minTag;
            fscanf(fp, "%d %d %d %d\n", &nb, &nn, &minTag, &maxNodeTag);
            nodes = calloc(maxNodeTag + 1, sizeof(Node));

            for (int b = 0; b < nb; b++) {
                int dim, tag, param, n;
                fscanf(fp, "%d %d %d %d\n", &dim, &tag, &param, &n);

                int *ids = malloc(sizeof(int) * n);
                for (int i = 0; i < n; i++) fscanf(fp, "%d\n", &ids[i]);

                for (int i = 0; i < n; i++)
                    fscanf(fp, "%lf %lf %lf\n",
                           &nodes[ids[i]].x,
                           &nodes[ids[i]].y,
                           &nodes[ids[i]].z);
                free(ids);
            }
        }
    }

    printf("Removed nodes (z=0 && -1.7<=x<=1.7):\n");
    for (int i = 1; i <= maxNodeTag; i++) {
        if (nodes[i].z == 0.0 &&
            nodes[i].x >= -1.7 &&
            nodes[i].x <= 1.7) {
            printf("node_id = %d (x=%.6f, z=%.6f)\n",
                   i, nodes[i].x, nodes[i].z);
        }
    }

    /* ---------- Output Nodes ---------- */
    fprintf(out, "# ===== Nodes =====\n");
    fprintf(out, "node_id,x,y,z\n");
    for (int i = 1; i <= maxNodeTag; i++)
        fprintf(out, "%d,%.10f,%.10f,%.10f\n",
                i, nodes[i].x, nodes[i].y, nodes[i].z);

    /* ---------- Elements → Edges ---------- */
    rewind(fp);
    fprintf(out, "\n# ===== Edges =====\n");
    fprintf(out, "edge_id,n1_id,n2_id,x1,y1,z1,x2,y2,z2\n");

    Edge *edges = malloc(sizeof(Edge) * MAX_EDGES);
    int edge_count = 0;
    int edge_id = 1;

    while (fgets(line, LINE_LEN, fp)) {
        if (strncmp(line, "$Elements", 9) == 0) {
            int nb, ne, mi, ma;
            fscanf(fp, "%d %d %d %d\n", &nb, &ne, &mi, &ma);

            for (int b = 0; b < nb; b++) {
                int dim, tag, type, n;
                fscanf(fp, "%d %d %d %d\n", &dim, &tag, &type, &n);

                for (int i = 0; i < n; i++) {
                    int eid;
                    fscanf(fp, "%d", &eid);

                    int v[4], nv = 0;
                    if (type == 1) { fscanf(fp, "%d %d\n", &v[0], &v[1]); nv = 2; }
                    else if (type == 2) { fscanf(fp, "%d %d %d\n", &v[0], &v[1], &v[2]); nv = 3; }
                    else if (type == 4) { fscanf(fp, "%d %d %d %d\n", &v[0], &v[1], &v[2], &v[3]); nv = 4; }
                    else { fgets(line, LINE_LEN, fp); continue; }

                    for (int a = 0; a < nv; a++) {
                        for (int b2 = a + 1; b2 < nv; b2++) {
                            int n1 = v[a], n2 = v[b2];

                            if (n1 == n2) continue;
                            if (is_removed_node(n1, nodes) || is_removed_node(n2, nodes)) continue;

                            if (n1 > n2) { int t = n1; n1 = n2; n2 = t; }
                            if (edge_exists(edges, edge_count, n1, n2)) continue;
                            if (edge_length(n1, n2, nodes) >= EDGE_LENGTH_LIMIT) continue;

                            edges[edge_count++] = (Edge){n1, n2};

                            fprintf(out,
                                "%d,%d,%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                                edge_id++, n1, n2,
                                nodes[n1].x, nodes[n1].y, nodes[n1].z,
                                nodes[n2].x, nodes[n2].y, nodes[n2].z);
                        }
                    }
                }
            }
        }
    }

    fclose(fp);
    fclose(out);
    free(nodes);
    free(edges);

    return 0;
}