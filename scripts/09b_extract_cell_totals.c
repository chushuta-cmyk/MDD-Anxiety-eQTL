#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
    if (argc < 4) return 1;
    FILE *fm = fopen(argv[1], "r");
    static char buf[8 * 1024 * 1024];
    setvbuf(fm, buf, _IOFBF, sizeof(buf));
    char line[512];
    while (fgets(line, sizeof(line), fm)) {
        if (line[0] != '%') break;
    }
    long n_genes = 0, n_cells = 0, n_nnz = 0;
    sscanf(line, "%ld %ld %ld", &n_genes, &n_cells, &n_nnz);
    int64_t *tot = calloc(n_cells, sizeof(int64_t));
    long r, c, v;
    while (fscanf(fm, "%ld %ld %ld", &r, &c, &v) == 3) {
        if (c >= 1 && c <= n_cells) tot[c - 1] += v;
    }
    fclose(fm);

    FILE *fq = fopen(argv[2], "rb");
    FILE *fo = fopen(argv[3], "wb");
    int32_t idx;
    while (fread(&idx, sizeof(int32_t), 1, fq) == 1) {
        fwrite(&tot[idx], sizeof(int64_t), 1, fo);
    }
    fclose(fq);
    fclose(fo);
    return 0;
}
