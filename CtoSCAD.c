#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LINE_LEN 1024

#define INVALID_VAL (-11999999999997.0)
#define REPLACE_VAL (3.1313558375)
#define TOL 1e-6

static inline double sanitize(double v)
{
    if (fabs(v - INVALID_VAL) < TOL)
        return REPLACE_VAL;
    return v;
}

int main(void)
{
    FILE *fin  = fopen("mesh_3_deformed_step_0500.csv", "r");
    FILE *fout = fopen("mesh_s500.txt", "w");

    if (!fin || !fout) {
        perror("File open error");
        return 1;
    }

    /* =========================================
       ノードごとの最大エッジ半径を保存
       ========================================= */
    double *max_r = NULL;
    int max_node_id = 0;

    char line[LINE_LEN];
    int mode = 0;

    /* ---------- 1st pass : edge 半径収集 ---------- */
    while (fgets(line, LINE_LEN, fin)) {

        line[strcspn(line, "\r\n")] = '\0';

        if (line[0] == '\0' || line[0] == '#')
            continue;

        if (strcmp(line, "node_id,x,y,z") == 0) {
            mode = 1;
            continue;
        }

        if (strcmp(line, "edge_id,n1_id,n2_id,x1,y1,z1,x2,y2,z2,r") == 0) {
            mode = 2;
            continue;
        }

        if (mode == 2) {
            int edge_id, n1_id, n2_id;
            double x1, y1, z1, x2, y2, z2, r;

            if (sscanf(line,
                       "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                       &edge_id, &n1_id, &n2_id,
                       &x1, &y1, &z1,
                       &x2, &y2, &z2, &r) == 10) {

                r = sanitize(r);

                int max_id = (n1_id > n2_id) ? n1_id : n2_id;
                if (max_id > max_node_id) {
                    max_r = realloc(max_r, (max_id + 1) * sizeof(double));
                    for (int i = max_node_id + 1; i <= max_id; i++)
                        max_r[i] = 0.0;
                    max_node_id = max_id;
                }

                if (r > max_r[n1_id]) max_r[n1_id] = r;
                if (r > max_r[n2_id]) max_r[n2_id] = r;
            }
        }
    }

    rewind(fin);

    /* =========================================================
       OpenSCAD module 定義 + union 開始
       ========================================================= */
    fprintf(fout,
        "module atom(r,x0,y0,z0)\n"
        "{\n"
        "  translate(v=[x0,y0,z0])\n"
        "  sphere(r=r*0.025,$fn=10);\n"
        "}\n\n"
        "module bond(x2,y2,z2,x1,y1,z1,rc)\n"
        "{\n"
        "  tx = (x2 + x1)/2;\n"
        "  ty = (y2 + y1)/2;\n"
        "  tz = (z2 + z1)/2;\n"
        "  ax = x2 - x1;\n"
        "  ay = y2 - y1;\n"
        "  az = z2 - z1;\n\n"
        "  translate(v=[tx,ty,tz])\n"
        "  rotate(a = [-acos(az/sqrt(ax*ax+ay*ay+az*az)), 0, -atan2(ax, ay)])\n"
        "  cylinder(r=rc*0.025,h=sqrt(ax*ax+ay*ay+az*az),center=true,$fn=10);\n"
        "}\n\n"
        "union(){\n\n"
    );

    mode = 0;

    /* ---------- 2nd pass : 出力 ---------- */
    while (fgets(line, LINE_LEN, fin)) {

        line[strcspn(line, "\r\n")] = '\0';

        if (line[0] == '\0' || line[0] == '#')
            continue;

        if (strcmp(line, "node_id,x,y,z") == 0) {
            mode = 1;
            continue;
        }

        if (strcmp(line, "edge_id,n1_id,n2_id,x1,y1,z1,x2,y2,z2,r") == 0) {
            mode = 2;
            continue;
        }

        if (mode == 1) {
            int node_id;
            double x, y, z;

            if (sscanf(line, "%d,%lf,%lf,%lf",
                       &node_id, &x, &y, &z) == 4) {

                x = sanitize(x);
                y = sanitize(y);
                z = sanitize(z);

                double r = (node_id <= max_node_id && max_r[node_id] > 0.0)
                           ? max_r[node_id]
                           : 0.15;

                fprintf(fout,
                        "  atom(%.10f, %.10f, %.10f, %.10f); // %d\n",
                        r, x, y, z, node_id);
            }
            continue;
        }

        if (mode == 2) {
            int edge_id, n1_id, n2_id;
            double x1, y1, z1, x2, y2, z2, r;

            if (sscanf(line,
                       "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                       &edge_id, &n1_id, &n2_id,
                       &x1, &y1, &z1,
                       &x2, &y2, &z2, &r) == 10) {

                x1 = sanitize(x1);
                y1 = sanitize(y1);
                z1 = sanitize(z1);
                x2 = sanitize(x2);
                y2 = sanitize(y2);
                z2 = sanitize(z2);
                r  = sanitize(r);

                fprintf(fout,
                        "  bond(%.10f, %.10f, %.10f, "
                        "%.10f, %.10f, %.10f, %.10f); // %d - %d\n",
                        x2, y2, z2,
                        x1, y1, z1, r,
                        n1_id, n2_id);
            }
            continue;
        }
    }

    fprintf(fout, "\n}\n");

    free(max_r);
    fclose(fin);
    fclose(fout);

    printf("Output written to mesh_du.txt\n");
    return 0;
}