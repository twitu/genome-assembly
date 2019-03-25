#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zhash.h"
#include "bucket_sort.h"
#include "llist.h"

#define MMER_SIZE 4 // efficient to keep mmer_size as powers of 2
#define KMER_SIZE 10

// defined constants for faster multiplication
const int power_val[] = {1, 4, 16, 64, 256, 1024, 4096, 16384};
// possible sizes for hash table; must be prime numbers
static const size_t hash_sizes[] = {
  53, 101, 211, 503, 1553, 3407, 6803, 12503, 25013, 50261,
  104729, 250007, 500009, 1000003, 2000029, 4000037, 10000019,
  25000009, 50000047, 104395301, 217645177, 512927357, 1000000007
};

// in line helper functions
#define MIN(A,B) \
   ({ __typeof__ (A) _A = (A); \
       __typeof__ (B) _B = (B); \
     _A < _B ? _A : _B; })

#define MAX(A,B) \
   ({ __typeof__ (A) _A = (A); \
       __typeof__ (B) _B = (B); \
     _A > _B ? _A : _B; })

// helper function for converting from numeric value to base pair
char getbp(int bp){
	switch(bp){
		case 0:
			return 'T';
		case 1:
			return 'G';
		case 2:
			return 'C';
		case 3:
			return 'A';
	}
}

// helper function for converting from base pair to numeric value
int getval(char c){
	switch(c){
		case 'T':
			return 0;
		case 'G':
			return 1;
		case 'C':
			return 2;
		case 'A':
			return 3;
		default :
			return -1;
	}
}

// helper function for generating canonical mmers
// returns next lexicographically smaller mmer than the one passed
// wraps around from AAAA to TTTT
char* next_smaller_mmer(char* mmer) {
    for (int i = MMER_SIZE - 1; i >= 0; i++) {
        if (mmer[i] == 'A') {
            mmer[i] == 'T';
        } else {
            mmer[i] = getbp(getval(mmer[i]) - 1);
            break;
        }
    }

    return mmer;
}

// merge lists in the forward direction from a to b
ll_node* merge_lists(int len, ll_node* a, ll_node* b, bool forward) {
    ll_node* temp = NULL;
    int i;

    // swap for backward direction
    if (!forward) {
        temp = a;
        a = b;
        b = temp;
    }


    // new_list points to starting list
    ll_node* new_list = a;


    // skip uptil last KMER_SIZE - 1 nodes of list
    for (i = 0; i < len - (KMER_SIZE - 1); i++) {
        a = a->next;
    }

    // merge entries of two kmers for KMER_SIZE - 1 nodes
    // nodes of b are freed as their values are transfered to a nodes
    for (i = 0; i < KMER_SIZE - 1; i++) {
        a->item = merge_sorted_list(a->item, b->item);

        temp = b;
        b = b->next;
        free(temp);

        if (i == KMER_SIZE - 1) {
            // if last node link rest of b nodes to new list
            a->next = b;
        } else {
            a = a->next;
        }
    }

    return new_list;
}

// iteration for level one hash using mmers
void* iterate_level_one_hash(struct ZHashTable* hash_table, bool indirection) {
    static struct ZHashTable* table = NULL;
    static struct ZHashEntry** entry;
    static int index;

    // reset variables table does not match hash_table
    if (table != hash_table) {
        table = hash_table;
        entry = NULL;
        index = 0;
    }

    // entry points to previously returned value in the same chain
    // move to next entry in the chain
    if (entry != NULL && *entry != NULL) {
        entry = &(*entry)->next;
    }

    // entry points to empty chain
    // move to next chain non empty chain
    if (entry == NULL || *entry == NULL) {
        while (index < hash_sizes[table->size_index]) {
            if (table->entries[index] != NULL) {
                entry = &table->entries[index];
                index++;
                break;
            }
            index++;
        }
    }
    
    // if iteration has ended keep returning NULL for subsequent calls
    if (entry == NULL) {
        return NULL;
    }

    if (indirection) {
        return entry;
    } else {
        return *entry;
    }
}

// iteration for level two hash using kmers
void* iterate_level_two_hash(struct ZHashTable* hash_table, bool indirection) {
    static struct ZHashTable* table = NULL;
    static struct ZHashEntry** entry;
    static int index;

    // reset variables table does not match hash_table
    if (table != hash_table) {
        table = hash_table;
        entry = NULL;
        index = 0;
    }

    // entry points to previously returned value in the same chain
    // move to next entry in the chain
    if (entry != NULL && *entry != NULL) {
        entry = &(*entry)->next;
    }

    // entry points to empty chain
    // move to next chain non empty chain
    if (entry == NULL || *entry == NULL) {
        while (index < hash_sizes[table->size_index]) {
            if (table->entries[index] != NULL) {
                entry = &table->entries[index];
                index++;
                break;
            }
            index++;
        }
    }
    
    // if iteration has ended keep returning NULL for subsequent calls
    if (entry == NULL) {
        return NULL;
    }

    if (indirection) {
        return entry;
    } else {
        return *entry;
    }
}

char* merge_keys(int a_len, int b_len, char* a, char* b, bool forward) {
    int len = a_len + b_len + 1 - (KMER_SIZE - 1);
    char* new_key = malloc(len);

    if (forward) {
        strncpy(new_key, a, a_len);
        strncpy(&new_key[a_len], &b[KMER_SIZE - 1], b_len - (KMER_SIZE - 1));
    } else {
        strncpy(new_key, b, b_len);
        strncpy(&new_key[b_len], &a[KMER_SIZE - 1], a_len - (KMER_SIZE - 1));
    }
    
    new_key[len] = '\0';
    return new_key;
}

// extend two kmer entries
void extend_kmers(struct ZHashTable* hash_table, struct ZHashEntry** a, struct ZHashEntry** b, bool forward) {
    struct ZHashEntry* temp;
    ll_node* new_read_ids;
    char* new_key = NULL;
    int a_len = strlen((*a)->key);
    int b_len = strlen((*b)->key);

    // merge read ids of both entries and concatenate keys
    new_read_ids = merge_lists(a_len, (ll_node*) (*a)->val, (ll_node*) (*b)->val, forward);
    new_key = merge_keys(a_len, b_len, (char*) (*a)->key, (char*) (*b)->key, forward);

    // free previous entries by non recursive method
    temp = (*a);
    a = &(*a)->next;
    zfree_entry(temp, false);
    temp = (*b);
    b = &(*b)->next;
    zfree_entry(temp, false);

    // store new entry
    zhash_set(hash_table, new_key, new_read_ids);
    free(new_key);
}

void print_kmers(struct ZHashTable* hash_table) {
    struct ZHashTable* kmer_hash;
    struct ZHashEntry* mmer_entry, *kmer_entry;
    ll_node* read_id, *traverse;
    
    while ((mmer_entry = (struct ZHashEntry*) iterate_level_one_hash(hash_table, false)) != NULL) {
        printf("%s\n", (mmer_entry)->key); // print mmer
        kmer_hash = (mmer_entry)->val;
        // iterate over all kmers of mmer
        while ((kmer_entry = (struct ZHashEntry*) iterate_level_two_hash(kmer_hash, false)) != NULL) {
            printf("%s\n", kmer_entry->key);
            read_id = (ll_node*) kmer_entry->val;
            // iterate over read id lists of each base pair
            while (read_id != NULL) {
                traverse = (ll_node*) read_id->item;
                // print each read id of base pair in same line
                while (traverse != NULL) {
                    printf("%d ", traverse->read_id);
                    traverse = traverse->next;
                }
                printf("\n");
                read_id = read_id->next;
            }
            printf("\n");
        }
    }
}

// stores all kmers of a read
// kmers are stored in 2 level hashing
// first level is hashed by signature of kmer
// second level is hashed by kmer itself
struct ZHashTable* process_read(struct ZHashTable* hash_table, char* read, int read_id) {
    int read_len = strlen(read);
    char* kmer = read;
    char* signature = NULL;
    int i, j;

    // initialize local variables for extracting signature
    char kmer_key[KMER_SIZE + 1];
    char mmer[MMER_SIZE + 1];
    char signature_cpy[MMER_SIZE + 1];
    int score = 0, rev_score = 0, max_score = 0;
    bool is_rev;
    int msb;

    // slide a window of KMER_SIZE over read
    for (i = 0; i < read_len - KMER_SIZE + 1; i++) {

        // iterate over the read and calculate the signature from scratch
        if (kmer > signature) {

            // store first MMER_SIZE characters in array
            for (j = 0; j < MMER_SIZE; j++) {
                mmer[j] = kmer[j];
                score = score*4 + getval(kmer[j]);
                rev_score = rev_score*4 + 3 - getval(kmer[j]);
            }
            mmer[MMER_SIZE] = '\0';

            // initialize max score
            if (score > rev_score) {
                max_score = score;
                is_rev = false;
            } else {
                max_score = rev_score;
                is_rev = true;
            }
            signature = kmer;
            msb = 0;

            // iterate over other characters
            j = MMER_SIZE;
            while (kmer[j] != '\0') {
                // incrementally update scores
                score = (score - getval(mmer[msb])*power_val[MMER_SIZE-1])*4;
                score += getval(kmer[j]);
                rev_score = (rev_score - ((3 - getval(mmer[msb]))*power_val[MMER_SIZE-1]))*4;
                rev_score += 3 - getval(kmer[j]);

                // rotate msb of mmer
                mmer[msb] = kmer[j];
                msb = (msb + 1)%MMER_SIZE;

                // increment pointer to next char
                j++;

                // change max score if required
                // set boolean if rev complement of mmer is signature
                if (MAX(score, rev_score) > max_score) {
                    if (score > rev_score) {
                        max_score = score;
                        is_rev = false;
                    } else {
                        max_score = rev_score;
                        is_rev = true;
                    }

                    // set new signature to point MMER_SIZE behind the next char to be added
                    signature = &kmer[j] - MMER_SIZE;
                }
            }
        } else {
            // compare current signature with new mmer 
            // new mmer is created by last letter added to the current kmer
            // store last MMER_SIZE characters of kmer in mmer array
            for (j = KMER_SIZE - MMER_SIZE; j < MMER_SIZE; j++) {
                mmer[j] = kmer[j];
                score = score*4 + getval(kmer[j]);
                rev_score = rev_score*4 + 3 - getval(kmer[j]);
            }
            mmer[MMER_SIZE] = '\0';

            if (MAX(score, rev_score) > max_score) {
                if (score > rev_score) {
                    max_score = score;
                    is_rev = false;
                } else {
                    max_score = rev_score;
                    is_rev = true;
                }

                // set signature as newly calculated mmer
                signature = &kmer[KMER_SIZE - MMER_SIZE];
            }
        }


        strncpy(signature_cpy, signature, MMER_SIZE);
        strncpy(kmer_key, kmer, KMER_SIZE);
        signature_cpy[MMER_SIZE] = '\0';
        kmer_key[KMER_SIZE] = '\0';

        // get reverse complement of signature if rev complement has higher score
        if (is_rev) {
            for (i = 0; i < MMER_SIZE; i++) {
                signature_cpy[i] = getbp(3 - getval(signature_cpy[i]));
            }

            for (i = 0; i < KMER_SIZE; i++) {
                kmer_key[i] = getbp(3 - getval(kmer_key[i]));
            }
        }

        // check if this mmer has been stored before
        // if not create a new hash table to store kmers for this signature
        struct ZHashTable* kmer_storage;
        if ((kmer_storage = zhash_get(hash_table, signature_cpy)) == NULL) {
            kmer_storage = zcreate_hash_table();
            zhash_set(hash_table, signature_cpy, kmer_storage);
        }
        
        // check if this kmer has been stored previously
        ll_node* read_id_list, *traverse;
        char temp;
        if ((read_id_list = zhash_get(kmer_storage, kmer_key)) == NULL) {
            
            // create linked list for all base pairs and store read id
            // creating list in reverse order so it matches order of bp in kmers
            for (j = 0; j < KMER_SIZE; j++) {
                traverse = (ll_node*) create_node_item((void*) create_node_num(read_id));
                traverse->next = read_id_list;
                read_id_list = traverse;
            }

            // store newly created linked list with kmer as key
            zhash_set(kmer_storage, kmer_key, read_id_list);
        } else {
            // kmer entry exists
            // add read id to read id list
            while (read_id_list != NULL) {
                traverse = (ll_node*) create_node_num(read_id);
                traverse->next = read_id_list->item;
                read_id_list->item = traverse;
                read_id_list = read_id_list->next;
            }
        }

        // increment kmer pointer
        kmer++;
    }

    return hash_table;
}

// prune kmers which have base pair occuring in only one read
// such kmers are highly likely to have been generated by errors
struct ZHashTable* prune_kmers(struct ZHashTable* hash_table) {
    struct ZHashEntry** traverse;
    struct ZHashEntry* temp;
    ll_node* read_id_list;
    for (int i = 0; i < hash_sizes[hash_table->size_index]; i++) {
        if (hash_table->entries[i] != NULL) {
            // entries exist for this hash value
            traverse = &hash_table->entries[i];
            while (*traverse != NULL) {
                read_id_list = (*traverse)->val;
                // check first entry of read id list to see if it has a next node
                if (((ll_node*) read_id_list->item)->next == NULL) {
                    // kmer has only one read id entry
                    temp = *traverse;
                    *traverse = ((struct ZHashEntry*) *traverse)->next;
                    zfree_entry(temp, false);
                } else {
                    traverse = &((struct ZHashEntry*) *traverse)->next;
                }
            }
        }
    }
}

// iterates over all entries of first level hash table
// applies prune_kmers function on all non null entries
struct ZHashTable* prune_data(struct ZHashTable* hash_table) {
    struct ZHashEntry* traverse;
    for (int i = 0; i < hash_sizes[hash_table->size_index]; i++) {
        if ((traverse = hash_table->entries[i]) != NULL) {
            // entries exist for this hash value
            while (traverse != NULL) {
                prune_kmers((struct ZHashTable*) traverse->val);
                traverse = traverse->next;
            }
        }
    }
}

int main() {
    // initialize file and structures
    FILE* file = fopen("input.txt", "r");
    struct ZHashTable* hash_table = zcreate_hash_table();
    
    // initialize variables
    char read[20];
    int read_id = 0;
    
    // get all the reads from file
    // expecting read of length less than 20
    while (fgets(read, 20, file) != NULL) {

        // pre process and store read
        int len = strlen(read);
        read[--len] = '\0';

        process_read(hash_table, read, read_id++);
    }

    // print kmers
    print_kmers(hash_table);

    // TODO: prune data and remove possibly erroneous kmers

    // TODO:
    // apply unitig extension to the data
    // extension can only occur 

}
