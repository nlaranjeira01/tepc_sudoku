#ifndef MINHEAP_H
#define MINHEAP_H
 
#include "dynarray.h"
 
struct entry
{
    void *item;
    unsigned int value;
};
typedef struct entry entry;
 
struct minheap
{
    dynarray *entries;
};
typedef struct minheap minheap;
 
typedef void(*minheap_forfn)(void*);
 
minheap *minheap_create(void);
void minheap_delete(minheap *heap);
void minheap_add(minheap *heap, void *item, unsigned int value);
void *minheap_remove_min(minheap *heap);
void minheap_for_each(const minheap *heap, minheap_forfn fun);
unsigned int minheap_get_count(const minheap *heap);
 
#endif /* MINHEAP_H */