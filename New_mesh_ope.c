#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LINE_LEN 2048

typedef struct {
    int id;
    double x, y, z;
} Node;

typedef struct {
    int id;
    int n1, n2;
    double r;   /* 元の半径 r_o */
    double L;
} Edge;

int main(void)
{
    FILE *fin = fopen("mesh_3_deformed_step_0100.csv", "r");
    FILE *fout_raw = fopen("mesh_raw_2.csv", "w");
    FILE *fout_uniform = fopen("mesh_uniform_2.csv", "w");

    if (!fin || !fout_raw || !fout_uniform) {
        fprintf(stderr, "File open error\n");
        return 1;
    }

    char line[LINE_LEN];
    int mode = 0;

    Node *nodes = NULL;
    int node_count = 0, node_cap = 0;

    Edge *edges = NULL;
    int edge_count = 0, edge_cap = 0;

    /* ===== 読み込み ===== */
    while (fgets(line, LINE_LEN, fin)) {

        if (strstr(line, "# ===== Nodes")) {
            mode = 1; fgets(line, LINE_LEN, fin); continue;
        }
        if (strstr(line, "# ===== Edges")) {
            mode = 2; fgets(line, LINE_LEN, fin); continue;
        }
        if (line[0] == '#' || strlen(line) < 3) continue;

        if (mode == 1) {
            Node n;
            if (sscanf(line, "%d,%lf,%lf,%lf",
                       &n.id, &n.x, &n.y, &n.z) != 4) continue;

            if (node_count >= node_cap) {
                node_cap = node_cap ? node_cap * 2 : 256;
                nodes = realloc(nodes, node_cap * sizeof(Node));
            }
            nodes[node_count++] = n;
        }

        if (mode == 2) {
            Edge e;
            double x1,y1,z1,x2,y2,z2;
            if (sscanf(line,
                "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                &e.id, &e.n1, &e.n2,
                &x1,&y1,&z1,&x2,&y2,&z2,&e.r) != 10) continue;

            double dx=x2-x1, dy=y2-y1, dz=z2-z1;
            e.L = sqrt(dx*dx + dy*dy + dz*dz);

            if (edge_count >= edge_cap) {
                edge_cap = edge_cap ? edge_cap * 2 : 256;
                edges = realloc(edges, edge_cap * sizeof(Edge));
            }
            edges[edge_count++] = e;
        }
    }
    fclose(fin);

    /* ===== r_min, r_max 計算 ===== */
    double r_max = edges[0].r;
    double r_min = 0.4;

    for (int i=1; i<edge_count; i++) {
        if (edges[i].r > r_max) r_max = edges[i].r;
        r_min = r_max/2;
    }

    /* ===== 使用ノード判定 ===== */
    int *used = calloc(node_count, sizeof(int));
    for (int i=0;i<edge_count;i++) {
        used[edges[i].n1] = 1;
        used[edges[i].n2] = 1;
    }

    /* ===== ノードID再割当 ===== */
    int *id_map = malloc(node_count * sizeof(int));
    int new_node_count = 0;

    for (int i=0;i<node_count;i++) {
        if (used[nodes[i].id])
            id_map[nodes[i].id] = new_node_count++;
        else
            id_map[nodes[i].id] = -1;
    }

    Node *new_nodes = malloc(new_node_count * sizeof(Node));
    for (int i=0;i<node_count;i++) {
        if (id_map[nodes[i].id] >= 0) {
            int nid = id_map[nodes[i].id];
            new_nodes[nid].id = nid;
            new_nodes[nid].x = nodes[i].x;
            new_nodes[nid].y = nodes[i].y;
            new_nodes[nid].z = nodes[i].z;
        }
    }

    for (int i=0;i<edge_count;i++) {
        edges[i].n1 = id_map[edges[i].n1];
        edges[i].n2 = id_map[edges[i].n2];
        edges[i].id = i;
    }

    /* ===== 体積保存用計算（変換後半径使用） ===== */
    double V=0.0, Lsum=0.0;
    for (int i=0;i<edge_count;i++) {
        double r_f = (r_max - r_min) / r_max * edges[i].r + r_min;
        V += M_PI * r_f * r_f * edges[i].L;
        Lsum += edges[i].L;
    }
    double r_uniform = sqrt(V / (M_PI * Lsum));

/* ===== 出力マクロ ===== */
#define OUTPUT(f, radius_expr) \
    fprintf(f,"# ===== Nodes =====\nnode_id,x,y,z\n"); \
    for (int i=0;i<new_node_count;i++) \
        fprintf(f,"%d,%.10f,%.10f,%.10f\n", \
            new_nodes[i].id,new_nodes[i].x,new_nodes[i].y,new_nodes[i].z); \
    fprintf(f,"\n# ===== Edges =====\n"); \
    fprintf(f,"edge_id,n1_id,n2_id,x1,y1,z1,x2,y2,z2,r\n"); \
    for (int i=0;i<edge_count;i++) { \
        double r_f = radius_expr; \
        fprintf(f,"%d,%d,%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n", \
            edges[i].id,edges[i].n1,edges[i].n2, \
            new_nodes[edges[i].n1].x,new_nodes[edges[i].n1].y,new_nodes[edges[i].n1].z, \
            new_nodes[edges[i].n2].x,new_nodes[edges[i].n2].y,new_nodes[edges[i].n2].z, \
            r_f); \
    }

    OUTPUT(fout_raw,
        (r_max - r_min) / r_max * edges[i].r + r_min)

    OUTPUT(fout_uniform,
        r_uniform)

#undef OUTPUT

    fclose(fout_raw);
    fclose(fout_uniform);

    free(nodes);
    free(edges);
    free(new_nodes);
    free(used);
    free(id_map);

    return 0;
}