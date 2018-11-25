#include "utils.h"

char *trimstr(char *s) {
  // remove whitespace at both ends of string. returns pointer 
  // to modified version of memory
  char *e;
  while(isspace((unsigned char)*s)) s++;
  if(*s == 0)  return s;
  e = s + strlen(s) - 1;
  while(e > s && isspace((unsigned char)*e)) e--;
  e[1] = '\0';
  return s;
}

double ****malloc4(int nx, int ny, int nz, int np) {
  // mallocs a 4d array with dimensions (nx,ny,nz,np)
  double ****mem = malloc(sizeof(*mem)*nx);
  if (mem == NULL) { fprintf(stderr, " ! bad malloc.\n"); exit(-3); }
  for (int i=0; i<nx; ++i) {
    mem[i] = malloc(sizeof(**mem)*ny);
    if (mem[i] == NULL) { fprintf(stderr, " ! bad malloc.\n"); exit(-4); }
    for (int j=0; j<ny; ++j) {
      mem[i][j] = malloc(sizeof(***mem)*nz);
      if (mem[i][j] == NULL) { fprintf(stderr, " ! bad malloc.\n"); exit(-5); }
      for (int k=0; k<nz; ++k) {
        mem[i][j][k] = malloc(sizeof(****mem)*np);
        if (mem[i][j][k] == NULL) { fprintf(stderr, " ! bad malloc.\n"); exit(-6); }
      }
    }
  }
  return mem;
}

void free4(int nx, int ny, int nz, int np, double ****mem) {
  // free a 4d array with dimensions (nx,ny,nz,np)
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      for (int k=0; k<nz; ++k) {
        free(mem[i][j][k]);
      }
      free(mem[i][j]);
    }
    free(mem[i]);
  }
  free(mem);
}

void *push_strlist(char *s, void *head) {
  // push a new node to the beginning of linked list. returns new head.
  strlist_el *el = malloc(sizeof(*el));
  strncpy(el->key, s, 255);
  if (head != NULL) el->next = head;
  return el;
}

int strlist_count(char *s, void *head) {
  // count how many nodes have key value "s" for list at head.
  strlist_el *el = head;
  int count=0;
  while (el != NULL) {
    if (strncmp(el->key, s, 255) == 0) count++;
    el = el->next;
  }
  return count;
}

void delete_strlist(void *head) {
  // frees all nodes in a linked list. be sure to set to NULL!
  while (head != NULL) {
    void *n = ((strlist_el *)head)->next;
    free(head);
    head = n;
  }
}

