#ifndef NEWMATRIX_H
#define NEWMATRIX_H
void **matrix(int row, int col, int num_bytes);
void **tria_matrix(int n, int num_bytes);
void free_matrix(void **matr);
#endif
