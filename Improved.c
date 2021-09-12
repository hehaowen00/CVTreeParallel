#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

long M, M1, M2;
short code[26] = {
    0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3
};

#define encode(ch)  code[ch-'A']
#define LEN         6
#define AA_NUMBER   20
#define	EPSILON	    1e-010

void init()
{
    M2 = 1;
    for (int i=0; i<LEN-2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
        M2 *= AA_NUMBER; 
    M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
    M  = M1 *AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
}

typedef struct  {
    long* vector;
    long* second;
    long* one_l;
    long indexs;
    long total;
    long total_l;
    long complement;
} Bacteria;

Bacteria * bacteria_new();
void bacteria_delete();
int bacteria_load(Bacteria *, char *);
void bacteria_init_buffer(Bacteria*, char*);
void bacteria_cont_buffer(Bacteria*, char);
double bacteria_compare(Bacteria *, Bacteria*);

Bacteria * bacteria_new()
{
    Bacteria * b = malloc(sizeof(Bacteria));
    b->vector = malloc(sizeof(long) * M);
    b->second = malloc(sizeof(long) * M1);
    b->one_l = malloc(sizeof(long) * AA_NUMBER);
    memset(b->vector, 0, M * sizeof(long));
    memset(b->second, 0, M1 * sizeof(long));
    memset(b->one_l, 0, AA_NUMBER * sizeof(long));
    b->total = 0;
    b->total_l = 0;
    b->complement = 0;
    return b;
}

void bacteria_delete(Bacteria * b)
{
    free(b->vector);
    free(b->second);
    free(b->one_l);
    free(b);
}

int bacteria_load(Bacteria * bacteria, char * filename)
{
    FILE * file =  fopen(filename, "r");
    if (file == NULL)
    {
        fprintf(stderr, "error: failed to open file %s\n", filename);
        exit(1);
    }

    char ch;
    while ((ch = fgetc(file)) != EOF)
    {
        if (ch == '>')
        {
            while (fgetc(file) != '\n');
            char buffer[LEN-1];
            fread(buffer, sizeof(char), LEN-1, file);
            bacteria_init_buffer(bacteria, buffer);
        }
        else if (ch != '\n')
        {
            bacteria_cont_buffer(bacteria, ch);
        }
    }

/*     fseek(file, 0, SEEK_END); */
/*     long long length = ftell(file); */
/*     fseek(file, 0, SEEK_SET); */
/*     char * file_buf = malloc(sizeof(char) * length); */
/*     fread(file_buf, 1, length, file); */

/*     long idx = 0; */
/*     char temp = file_buf[idx]; */

/*     do { */
/*         if (temp == '>') */
/*         { */
/*             while (file_buf[idx++] != '\n'); */
/*             char buffer[LEN-1]; */
/*             memcpy(buffer, file_buf + idx, LEN-1); */
/*             bacteria_init_buffer(bacteria, buffer); */
/*         } */
/*         else if (temp != '\n') */
/*         { */
/*             bacteria_cont_buffer(bacteria, temp); */
/*         } */

/*         if (idx + 1 >= length) */
/*         { */
/*             break; */
/*         } */
/*         idx++; */
/*         temp = file_buf[idx]; */
/*     } while (true); */

    fclose(file);
    return 0;
}

void bacteria_init_buffer(Bacteria * b, char * buf)
{
    b->complement++;
    b->indexs = 0;
    for (int i = 0; i < LEN-1; i++)
    {
        short enc = encode(buf[i]);
        b->one_l[enc]++;
        b->total_l++;
        b->indexs = b->indexs * AA_NUMBER + enc;
    }
    b->second[b->indexs]++;
}

void bacteria_cont_buffer(Bacteria * b, char c)
{
    short enc = encode(c);
    b->one_l[enc]++;
    b->total_l++;
    long index = b->indexs * AA_NUMBER + enc;
    b->vector[index]++;
    b->total++;
    b->indexs = (b->indexs % M2) * AA_NUMBER + enc;
    b->second[b->indexs]++;
}

double bacteria_compare(Bacteria *b1, Bacteria *b2)
{
    return 0;
}

void find_pairs()
{
    for (int i = 0; i < 41 && i + 1 < 41; i++)
    {
        for (int j = i+1; j < 41; j++)
        {
            printf("%02d, %02d\n", i, j);
        }
    }
}

int main(int argc, char * argv[])
{
    init();

    Bacteria * b = bacteria_new();
    bacteria_load(b, "./data/AcMNPV.faa");
    printf("total: %ld, total_l: %ld\n", b->total, b->total_l);
    printf("one_l: ");
    for (int i = 0; i < AA_NUMBER; i++)
        printf("%ld, ", b->one_l[i]);
    printf("\n");
    printf("indexs: %ld, complement: %ld\n", b->indexs, b->complement);

    bacteria_delete(b);
    /* find_pairs(); */

    return 0;
}

// Assignment TODO:
// Port improved.cpp to Improved.c
// Pool memory
// Parallelize this using threadpool and SIMD
// Make a zero copy GPU accelerated version (CUDA Streams)

// Architecture:
// Load data from files into bacteria arrays

// TODO: 
// Contiguous memory allocations
// Keep most recently used in memory
// Work stealing threadpool?
// SIMD text processing
// Go from array of structs to struct of arrays to avoid repeated allocations

// Questions:
// 1. Is time_t suitable for benchmarking
// 2. Can we run the parallel code in a VM
// 3. Writing it in C because Rust won't work properly D:
// 4. Sad noises

/* total: 41036, total_l: 41816 */
/* one_l: 2178, 1052, 2518, 2231, 2083, 1365, 945, 2731, 3010, 3870, 1147, 3310, 1642, 1597, 1926, 2662, 2362, 2803, 319, 2065, */
/* indexs: 218967, complement: 156 */
