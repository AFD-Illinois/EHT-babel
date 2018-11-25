#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// 4d -> 1d index
#define IJKP(I,J,K,P,NY,NZ,NP) (P+(K+(J+I*NY)*NZ)*NP)

// string utils (parsing)
char *trimstr(char *s);

// memory
double ****malloc4(int nx, int ny, int nz, int np);
void free4(int nx, int ny, int nz, int np, double ****mem);

// linked list
typedef struct {
  char key[256];
  void *next;
} strlist_el;
void *push_strlist(char *s, void *head);
int strlist_count(char *s, void *head);
void delete_strlist(void *head);

#endif // UTILS_H