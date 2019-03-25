#ifndef LLIST_H
#define LLIST_H

#include <stdlib.h>

typedef struct ll_node {
    struct ll_node* next;
    union {
        int read_id;
        void* item;
    };
} ll_node;

ll_node* create_node_num(int id);
ll_node* create_node_item(void* item);
ll_node* merge_sorted_list(ll_node* a, ll_node* b);

#endif

