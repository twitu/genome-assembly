#ifndef LLIST_H
#define LLIST_H

#include <stdlib.h>
#include <stdbool.h>

typedef struct ll_node {
    struct ll_node* next;
    union {
        int read_id;
        void* item;
    };
} ll_node;

typedef struct ll_node_queue {
    ll_node* head;
    ll_node* tail;
} ll_node_queue;

// node creator functions
ll_node* create_node_num(int id);
ll_node* create_node_item(void* item);

// queue operations
ll_node_queue* create_queue();
void enqueue_item(ll_node_queue* queue, void* item);
void* dequeue_item(ll_node_queue* queue);
bool queue_is_empty(ll_node_queue* queue);

// llist manipulator functions
ll_node* merge_sorted_list(ll_node* a, ll_node* b);
ll_node* duplicate_llist(ll_node* list);
void free_llist(ll_node* list);

#endif

