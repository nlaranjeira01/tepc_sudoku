#include <stdlib.h>
 
#include "minheap.h"
 
static entry *entry_create(void *item, unsigned int value)
{
    entry *e = malloc(sizeof(entry));
    if (e) {
        e->item = item;
        e->value = value;
    }
    return e;
}
 
minheap *minheap_create(void)
{
    minheap *heap = malloc(sizeof(minheap));
    if (heap) {
        heap->entries = dynarray_create(13);
    }
    return heap;
}
 
void minheap_delete(minheap *heap)
{
    if (heap) {
        dynarray_for_each(heap->entries, free);
        dynarray_delete(heap->entries);
        free(heap);
    }
}
 
static void minheap_swap(minheap *heap, unsigned int index1, unsigned int index2)
{
    void *temp = dynarray_get(heap->entries, index1);
    dynarray_set(heap->entries, index1, dynarray_get(heap->entries, index2));
    dynarray_set(heap->entries, index2, temp);
}
 
static void minheap_bubble_up(minheap *heap, unsigned int index)
{
    entry *e = dynarray_get(heap->entries, index);
    unsigned int parent_index = (index - 1) / 2;
    entry *parent = dynarray_get(heap->entries, parent_index);
    if (e->value < parent->value) {
        minheap_swap(heap, index, parent_index);
        if (parent_index > 0) {
            minheap_bubble_up(heap, parent_index);
        }
    }
}
 
static void minheap_bubble_down(minheap *heap, unsigned int index)
{
    entry *e = dynarray_get(heap->entries, index);
    unsigned int left_child_index = (index * 2) + 1;
    unsigned int right_child_index = left_child_index + 1;
    unsigned int swapped = 0;
    unsigned int swapped_index;
    if (right_child_index < dynarray_get_count(heap->entries) /* There is a right child */
            && ((entry*)dynarray_get(heap->entries, right_child_index))->value 
                    < ((entry*)dynarray_get(heap->entries, left_child_index))->value 
                   /* And it's less than left child */
            && e->value > ((entry*)dynarray_get(heap->entries, right_child_index))->value) {
        minheap_swap(heap, index, right_child_index);
        swapped = 1;
        swapped_index = right_child_index;
    }
    else if (e->value > ((entry*)dynarray_get(heap->entries, left_child_index))->value) {
        minheap_swap(heap, index, left_child_index);
        swapped = 1;
        swapped_index = left_child_index;
    }
    if (swapped && (swapped_index * 2) + 1 < dynarray_get_count(heap->entries) - 1) {
        minheap_bubble_down(heap, swapped_index);
    }
}
 
void minheap_add(minheap *heap, void *item, unsigned int value)
{
    entry *e = entry_create(item, value);
    if (e) {
        dynarray_add_tail(heap->entries, e);
        unsigned int count = dynarray_get_count(heap->entries);
        if (count > 1) {
            minheap_bubble_up(heap, count - 1);
        }
    }
}
 
void *minheap_remove_min(minheap *heap)
{
    void *item = NULL;
    unsigned int count = dynarray_get_count(heap->entries);
    if (count > 1) {
        minheap_swap(heap, 0, count - 1);
    }
    if (count > 0) {
        entry *e = dynarray_remove_tail(heap->entries);
        item = e->item;
        free(e);
    }
    if (dynarray_get_count(heap->entries) > 1) {
        minheap_bubble_down(heap, 0);
    }
    return item;
}
 
void minheap_for_each(const minheap *heap, minheap_forfn fun)
{
    dynarray_for_each(heap->entries, (dynarray_forfn)fun);
}
 
unsigned int minheap_get_count(const minheap *heap)
{
    return dynarray_get_count(heap->entries);
}