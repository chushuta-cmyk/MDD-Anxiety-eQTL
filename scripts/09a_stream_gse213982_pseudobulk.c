#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// Reads MatrixMarket header + entries fast with custom buffer
// Usage: stream_mtx <mtx_path> <n_cells> <n_genes> <cell_qc_in.bin> <pass1_or_2> ...
// Actually simpler: two passes directly inside one C binary!
// Input map: cell_map.tsv (160711 lines: donor_id_idx \t subtype_idx \t is_mix)
// Outputs:
//  1) cell_qc_summary.tsv (total_cells, qc_pass_cells, per-sample counts)
//  2) donor_subtype_cell_counts.tsv (71 x 41 matrix)
//  3) pb_counts.bin + pb_groups.tsv (all non-empty donor x subtype pseudobulk sums for QC-pass cells)
//  4) lr_counts.mtx (sparse MatrixMarket for genes in lr_gene_mask.bin across QC-pass cells)

int main(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr, "Usage: %s <mtx> <cell_map.bin> <lr_mask.bin> <n_donors> <n_subtypes> <out_dir>\n", argv[0]);
        return 1;
    }
    const char *mtx_path = argv[1];
    const char *cell_map_path = argv[2];
    const char *lr_mask_path = argv[3];
    int n_donors = atoi(argv[4]);
    int n_subtypes = atoi(argv[5]);
    const char *out_dir = argv[6];

    FILE *fm = fopen(mtx_path, "r");
    if (!fm) { perror("fopen mtx"); return 1; }
    static char io_buf[8 * 1024 * 1024];
    setvbuf(fm, io_buf, _IOFBF, sizeof(io_buf));

    char line[512];
    // skip comment lines
    while (fgets(line, sizeof(line), fm)) {
        if (line[0] != '%') break;
    }
    long n_genes = 0, n_cells = 0, n_nnz = 0;
    if (sscanf(line, "%ld %ld %ld", &n_genes, &n_cells, &n_nnz) != 3) {
        fprintf(stderr, "Bad header: %s\n", line);
        return 1;
    }
    fprintf(stderr, "[C] MTX dims: %ld genes x %ld cells, %ld nnz\n", n_genes, n_cells, n_nnz);

    // Load cell_map.bin: 3 int32 per cell (donor_idx, subtype_idx, is_mix)
    int32_t *cell_donor = malloc(n_cells * sizeof(int32_t));
    int32_t *cell_sub   = malloc(n_cells * sizeof(int32_t));
    uint8_t *cell_mix   = malloc(n_cells * sizeof(uint8_t));
    FILE *fc = fopen(cell_map_path, "rb");
    if (!fc) { perror("cell_map"); return 1; }
    for (long i = 0; i < n_cells; i++) {
        int32_t buf[3];
        if (fread(buf, sizeof(int32_t), 3, fc) != 3) { fprintf(stderr, "short read cell_map\n"); return 1; }
        cell_donor[i] = buf[0];
        cell_sub[i]   = buf[1];
        cell_mix[i]   = (uint8_t)buf[2];
    }
    fclose(fc);

    uint8_t *lr_mask = calloc(n_genes, sizeof(uint8_t));
    FILE *fl = fopen(lr_mask_path, "rb");
    if (!fl) { perror("lr_mask"); return 1; }
    if (fread(lr_mask, sizeof(uint8_t), n_genes, fl) != (size_t)n_genes) {
        fprintf(stderr, "short read lr_mask\n"); return 1;
    }
    fclose(fl);

    // Pass 1: count n_genes_by_counts & total_counts per cell
    int32_t *n_genes_cell = calloc(n_cells, sizeof(int32_t));
    int64_t *tot_cnt_cell = calloc(n_cells, sizeof(int64_t));

    long data_offset = ftell(fm);
    long r, c, v;
    long lines_read = 0;
    while (fscanf(fm, "%ld %ld %ld", &r, &c, &v) == 3) {
        long ci = c - 1;
        if (ci >= 0 && ci < n_cells) {
            n_genes_cell[ci]++;
            tot_cnt_cell[ci] += v;
        }
        lines_read++;
    }
    fprintf(stderr, "[C] Pass 1 done (%ld entries).\n", lines_read);

    // Apply QC matching Notebook 03: n_genes > 500 && n_genes < 8000 && tot_cnt > 1000 && !is_mix
    uint8_t *qc_pass = calloc(n_cells, sizeof(uint8_t));
    int32_t *new_cell_idx = malloc(n_cells * sizeof(int32_t));
    int32_t n_pass = 0;
    int32_t *ds_cell_cnt = calloc((size_t)n_donors * n_subtypes, sizeof(int32_t));

    for (long i = 0; i < n_cells; i++) {
        if (n_genes_cell[i] > 500 && n_genes_cell[i] < 8000 && tot_cnt_cell[i] > 1000 && cell_mix[i] == 0) {
            qc_pass[i] = 1;
            new_cell_idx[i] = n_pass++;
            ds_cell_cnt[(size_t)cell_donor[i] * n_subtypes + cell_sub[i]]++;
        } else {
            new_cell_idx[i] = -1;
        }
    }
    fprintf(stderr, "[C] QC pass cells (pre-doublet): %d / %ld\n", n_pass, n_cells);

    // Assign compact group IDs to non-empty (donor, subtype) pairs
    int32_t *ds_to_grp = malloc((size_t)n_donors * n_subtypes * sizeof(int32_t));
    int n_groups = 0;
    char grp_path[1024];
    snprintf(grp_path, sizeof(grp_path), "%s/pb_groups.tsv", out_dir);
    FILE *fg = fopen(grp_path, "w");
    fprintf(fg, "grp_idx\tdonor_idx\tsubtype_idx\tn_cells\n");
    for (int d = 0; d < n_donors; d++) {
        for (int s = 0; s < n_subtypes; s++) {
            int cnt = ds_cell_cnt[(size_t)d * n_subtypes + s];
            if (cnt > 0) {
                ds_to_grp[(size_t)d * n_subtypes + s] = n_groups;
                fprintf(fg, "%d\t%d\t%d\t%d\n", n_groups, d, s, cnt);
                n_groups++;
            } else {
                ds_to_grp[(size_t)d * n_subtypes + s] = -1;
            }
        }
    }
    fclose(fg);
    fprintf(stderr, "[C] Non-empty (donor, subtype) pseudobulks: %d (PB array: %.1f MB)\n",
            n_groups, (double)n_groups * n_genes * sizeof(int32_t) / (1024.0 * 1024.0));

    // Allocate pseudobulk array: [n_groups x n_genes] int32
    int32_t *pb = calloc((size_t)n_groups * n_genes, sizeof(int32_t));
    if (!pb) { fprintf(stderr, "Failed to alloc pb matrix\n"); return 1; }

    // Write sparse LR matrix & QC-passed cell indices
    char qc_path[1024], lr_path[1024];
    snprintf(qc_path, sizeof(qc_path), "%s/qc_pass_orig_idx.bin", out_dir);
    FILE *fq = fopen(qc_path, "wb");
    for (long i = 0; i < n_cells; i++) {
        if (qc_pass[i]) {
            int32_t idx32 = (int32_t)i;
            fwrite(&idx32, sizeof(int32_t), 1, fq);
        }
    }
    fclose(fq);

    snprintf(lr_path, sizeof(lr_path), "%s/lr_triplets.bin", out_dir);
    FILE *flr = fopen(lr_path, "wb");
    static char lr_buf[4 * 1024 * 1024];
    setvbuf(flr, lr_buf, _IOFBF, sizeof(lr_buf));

    // Pass 2: aggregate pseudobulk + extract LR triplets
    fseek(fm, data_offset, SEEK_SET);
    long lr_nnz = 0;
    while (fscanf(fm, "%ld %ld %ld", &r, &c, &v) == 3) {
        long ci = c - 1;
        long gi = r - 1;
        if (ci >= 0 && ci < n_cells && gi >= 0 && gi < n_genes && qc_pass[ci]) {
            int gidx = ds_to_grp[(size_t)cell_donor[ci] * n_subtypes + cell_sub[ci]];
            pb[(size_t)gidx * n_genes + gi] += (int32_t)v;
            if (lr_mask[gi]) {
                int32_t tri[3] = {(int32_t)gi, new_cell_idx[ci], (int32_t)v};
                fwrite(tri, sizeof(int32_t), 3, flr);
                lr_nnz++;
            }
        }
    }
    fclose(fm);
    fclose(flr);
    fprintf(stderr, "[C] Pass 2 done. LR nonzeros extracted: %ld\n", lr_nnz);

    char pb_path[1024];
    snprintf(pb_path, sizeof(pb_path), "%s/pb_counts.bin", out_dir);
    FILE *fpb = fopen(pb_path, "wb");
    fwrite(pb, sizeof(int32_t), (size_t)n_groups * n_genes, fpb);
    fclose(fpb);

    fprintf(stderr, "[C] ALL DONE cleanly!\n");
    return 0;
}
