#include <stdlib.h>

#include "llist.h"

ll_node* create_node_num(int id) {
    ll_node* new_node = malloc(sizeof(ll_node));
    new_node->next = NULL;
    new_node->read_id = id;
    return new_node;
}

ll_node* create_node_item(void* item) {
    ll_node* new_node = malloc(sizeof(ll_node));
    new_node->next = NULL;
    new_node->item = item;
    return new_node;
}

ll_node* merge_sorted_list(ll_node* a, ll_node* b) {

    ll_node* sorted = NULL;
    ll_node** traverse = &sorted;

    while (a != NULL && b != NULL) {
        if (a->read_id > b->read_id) {
            *traverse = a;
            a = a->next;
        } else if (a->read_id < b->read_id) {
            *traverse = b;
            b = b->next;
        } else {
            // equality case free one of the nodes
            *traverse = a;
            a = a->next;
            ll_node* temp = b;
            b = b->next;
            free(temp);
        }

        // move traverse pointer to position of next node to be added
        traverse = &(*traverse)->next;
    }

    // if any list is not exhausted append to existing
    if (a != NULL) {
        *traverse = a;
    } else {
        *traverse = b;
    }

    return sorted;
}
