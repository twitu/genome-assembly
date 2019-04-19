#include <stdlib.h>
#include <stdbool.h>

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

ll_node_queue* create_queue() {
    return calloc(1, sizeof(ll_node_queue));
}

void enqueue_item(ll_node_queue* queue, void* item) {
    if (queue->head == NULL) {
        queue->head = create_node_item(item);
        queue->tail = queue->head;
    } else {
        queue->tail->next = create_node_item(item);
        queue->tail = queue->tail->next;
    }
}

void* dequeue_item(ll_node_queue* queue) {
    void* to_return = queue->head->item;
    ll_node* temp = queue->head;
    queue->head = queue->head->next;
    free(temp);
    return to_return;
}

bool queue_is_empty(ll_node_queue* queue) {
    return (queue->head == NULL) ? true : false;
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
    }

    if (b != NULL) {
        *traverse = b;
    }

    return sorted;
}

ll_node* duplicate_llist(ll_node* list) {
    ll_node* traverse = NULL, *new_list = NULL;

    while (list != NULL) {
        if (traverse == NULL) {
            traverse = create_node_num(list->read_id);
            new_list = traverse;
        } else {
            traverse->next = create_node_num(list->read_id);
            traverse = traverse->next;
        }
        
        list = list->next;
    }

    return new_list;
}

void free_llist(ll_node* list) {
    ll_node* traverse = NULL;
    while (list != NULL) {
        traverse = list;
        list = list->next;
        free(traverse);
    }
}
