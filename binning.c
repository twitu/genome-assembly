#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zhash.h"
#include "bucket_sort.h"
#include "llist.h"

#define MMER_SIZE 4 // efficient to keep mmer_size as powers of 2
#define KMER_SIZE 31
#define ABUNDANCE_CUTOFF 1
#define READ_LENGTH 101

// defined constants for faster multiplication
const int power_val[] = {1, 4, 16, 64, 256, 1024, 4096, 16384};
// possible sizes for hash table; must be prime numbers
static const size_t hash_sizes[] = {
  53, 101, 211, 503, 1553, 3407, 6803, 12503, 25013, 50261,
  104729, 250007, 500009, 1000003, 2000029, 4000037, 10000019,
  25000009, 50000047, 104395301, 217645177, 512927357, 1000000007
};

typedef struct kmer_extension_node {
    struct ZHashEntry** extend_entry;
    struct ZHashTable* extend_table;
} kmer_extension_node;

typedef struct more_kmer_extension_node {
    char* key;
    ll_node* read_id_lists;
} more_kmer_extension_node;

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

// helper function to get score of string
int getscore(char* string) {
    int score = 0;
    while (*string != '\0') {
        score = score*4 + getval(*string);
        string++;
    }
    
    return score;
}

// helper function for generating canonical mmers
// returns next lexicographically smaller mmer than the one passed
// wraps around from AAAA to TTTT
int next_smaller_mmer(char* mmer, int mmer_score) {
    for (int i = MMER_SIZE - 1; i >= 0; i--) {
        if (mmer[i] == 'A') {
            mmer[i] = 'T';
        } else {
            mmer[i] = getbp(getval(mmer[i]) + 1);
            break;
        }
    }

    return ++mmer_score;
}

// merge lists in the forward direction from a to b
ll_node* merge_lists(int a_len, int b_len, ll_node* aa, ll_node* bb, bool forward) {
    ll_node* temp = NULL;
    int i, len = a_len;
ll_node* a = aa, *b = bb;
    // swap for backward direction
    if (!forward) {
        temp = a;
        a = b;
        b = temp;
        len = b_len;
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

        if (i == KMER_SIZE - 2) {
            // if last node link rest of b nodes to new list
            a->next = b;
        } else {
            a = a->next;
        }
    }

    return new_list;
}

bool compare_overlap(char* a, char* b, bool forward) {
    char* temp;
    int len;

    // swap to always keep a on left side of overlap
    // len contains length of string that will be left side of overlap
    if (!forward) {
        temp = a;
        a = b;
        b = temp;
    }

    len = strlen(a);

    for (int i = 0; i < KMER_SIZE - 1; i++) {
        if (a[len - (KMER_SIZE - 1) + i] != b[i]) return false;
    }

    return true;
}

// merge key strings in both backward or forward direction depending on parameter
char* merge_keys(int a_len, int b_len, char* a, char* b, bool forward) {
    int len = a_len + b_len + 1 - (KMER_SIZE - 1);
    // critical to use calloc, malloc can give unitialized string which can cause error
    char* new_key = calloc(len, sizeof(char));

    if (forward) {
        strncpy(new_key, a, a_len);
        strncpy(&new_key[a_len], &b[KMER_SIZE - 1], b_len - (KMER_SIZE - 1));
    } else {
        strncpy(new_key, b, b_len);
        strncpy(&new_key[b_len], &a[KMER_SIZE - 1], a_len - (KMER_SIZE - 1));
    }
    
    return new_key;
}

// iteration for level one hash using mmers
void* iterate_level_one_hash(struct ZHashTable* hash_table, bool indirection, bool remove_current) {
    static struct ZHashTable* table = NULL;
    static struct ZHashEntry** entry;
    static int index;
    static bool remove;

    // remove currently pointed node
    // when the next iterate is called
    if (hash_table == NULL && remove_current) {
        remove = remove_current;
        return NULL;
    }

    // reset variables table does not match hash_table
    if (table != hash_table) {
        table = hash_table;
        entry = NULL;
        index = 0;
    }

    // entry points to previously returned value in the same chain
    // move to next entry in the chain
    if (entry != NULL && *entry != NULL) {
        if (!remove) {
            entry = &(*entry)->next;
        } else {
            // if current node is marked for removal handle differently
            struct ZHashEntry* temp = *entry;
            *entry = (*entry)->next;
            zfree_entry(temp, false);
            remove = false;
        }
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
    if (entry == NULL || *entry == NULL) {
        table = NULL;
        return NULL;
    }

    if (indirection) {
        return entry;
    } else {
        return *entry;
    }
}

// iteration for level two hash using kmers
void* iterate_level_two_hash(struct ZHashTable* hash_table, bool indirection, bool remove_current) {
    static struct ZHashTable* table = NULL;
    static struct ZHashEntry** entry;
    static int index;
    static bool remove;

    // remove currently pointed node
    // when the next iterate is called
    if (hash_table == NULL && remove_current) {
        remove = remove_current;
        return NULL;
    }

    // reset variables table does not match hash_table
    if (table != hash_table) {
        table = hash_table;
        entry = NULL;
        index = 0;
    }

    // entry points to previously returned value in the same chain
    // move to next entry in the chain
    if (entry != NULL && *entry != NULL) {
        if (!remove) {
            entry = &(*entry)->next;
        } else {
            // if current node is marked for removal handle differently
            struct ZHashEntry* temp = *entry;
            *entry = (*entry)->next;
            zfree_entry(temp, false);
            remove = false;
        }
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
    if (entry == NULL || *entry == NULL) {
        table = NULL;
        return NULL;
    }

    if (indirection) {
        return entry;
    } else {
        return *entry;
    }
}

// extend two kmer entries
more_kmer_extension_node extend_kmers(struct ZHashEntry* a, struct ZHashEntry* b, bool forward) {
    int a_len = strlen(a->key);
    int b_len = strlen(b->key);

    // merge read ids of both entries and concatenate keys
    ll_node* new_read_ids = merge_lists(a_len, b_len, (ll_node*) a->val, (ll_node*) b->val, forward);
    char* new_key = merge_keys(a_len, b_len, (char*) a->key, (char*) b->key, forward);
    more_kmer_extension_node to_return;
    to_return.key = new_key;
    to_return.read_id_lists = new_read_ids;
    return to_return;
}

more_kmer_extension_node further_extend_kmers(more_kmer_extension_node a, struct ZHashEntry* b, bool forward) {
    int a_len = strlen(a.key);
    int b_len = strlen(b->key);

    // merge read ids of both entries and concatenate keys
    ll_node* new_read_ids = merge_lists(a_len, b_len, (ll_node*) a.read_id_lists, (ll_node*) b->val, forward);
    char* new_key = merge_keys(a_len, b_len, (char*) a.key, (char*) b->key, forward);

    free(a.key);
    a.key = new_key;
    a.read_id_lists = new_read_ids;
    return a;
}

// find kmer extension from existing table else return empty kmer_extension_node
kmer_extension_node find_kmer_extension(struct ZHashTable* hash_table, struct ZHashEntry* entry, int mmer_score, bool forward) {
    // initialize key structure
    char* key = entry->key;
    int key_len = strlen(key);

    // initialize signature mmer
    char compare_mmer[MMER_SIZE + 1];
    strncpy(compare_mmer, &key[key_len - (MMER_SIZE - 1)], MMER_SIZE - 1);

    bool multiple_extension = false;
    struct ZHashEntry** extend_entry = NULL, **compare_entry = NULL;
    struct ZHashTable* compare_mmer_hash = NULL, *extend_table = NULL;

    for (int i = 0; i < 4; i++) {
        compare_mmer[MMER_SIZE-1] = getbp(i);

        if (getscore(compare_mmer) > mmer_score) {
            // extension only with lexicographically larger mmers
            continue;
        }

        if ((compare_mmer_hash = (struct ZHashTable*) zhash_get(hash_table, compare_mmer)) == NULL) {
            // ignore if mmer does not have entry
            continue;
        }

        // compare signature is lexicographically greater than or equal
        while((compare_entry = iterate_level_two_hash(compare_mmer_hash, true, false)) != NULL){
            // handle equality case
            if (*compare_entry == entry) continue;

            // check if KMER_SIZE - 1 characters overlap
            if (!compare_overlap(key, (*compare_entry)->key, forward)) continue;

            // if extension entry already exists
            // there are multiple possible extensions
            // unitig extension is not possible
            if (extend_entry != NULL) {
                extend_entry = NULL;
                extend_table = NULL;
                multiple_extension = true;
                break;
            } else {
                extend_table = compare_mmer_hash;
                extend_entry = compare_entry;
            }
        }

        if (multiple_extension) {
            break;
        }
    }

    // deduct count from entry table
    kmer_extension_node to_return;
    to_return.extend_entry = extend_entry;
    to_return.extend_table = extend_table;
    return to_return;
}

// find further kmer extension from extended kmer else return empty further_kmer_extension_node
kmer_extension_node more_kmer_extension(struct ZHashTable* hash_table, char* key, int mmer_score, bool forward) {
    // initialize key structure
    int key_len = strlen(key);

    // initialize signature mmer
    char compare_mmer[MMER_SIZE + 1];
    strncpy(compare_mmer, &key[key_len - (MMER_SIZE - 1)], MMER_SIZE - 1);

    bool multiple_extension = false;
    struct ZHashEntry** extend_entry = NULL, **compare_entry = NULL;
    struct ZHashTable* compare_mmer_hash = NULL, *extend_table = NULL;

    for (int i = 0; i < 4; i++) {
        compare_mmer[MMER_SIZE-1] = getbp(i);

        if (getscore(compare_mmer) > mmer_score) {
            // extension only with lexicographically larger mmers
            continue;
        }

        if ((compare_mmer_hash = (struct ZHashTable*) zhash_get(hash_table, compare_mmer)) == NULL) {
            // ignore if mmer does not have entry
            continue;
        }

        // compare signature is lexicographically greater than or equal
        while((compare_entry = iterate_level_two_hash(compare_mmer_hash, true, false)) != NULL){
            // check if KMER_SIZE - 1 characters overlap
            if (!compare_overlap(key, (*compare_entry)->key, forward)) continue;

            // if extension entry already exists
            // there are multiple possible extensions
            // unitig extension is not possible
            if (extend_entry != NULL) {
                extend_entry = NULL;
                extend_table = NULL;
                multiple_extension = true;
                break;
            } else {
                extend_table = compare_mmer_hash;
                extend_entry = compare_entry;
            }
        }

        if (multiple_extension) {
            break;
        }
    }

    // deduct count from entry table
    kmer_extension_node to_return;
    to_return.extend_entry = extend_entry;
    to_return.extend_table = extend_table;
    return to_return;
}

void find_kmer_extensions(struct ZHashTable* hash_table, bool forward) {
    // initialize signature kmer
    char mmer[MMER_SIZE + 1];
    mmer[0] = 'C';
    mmer[MMER_SIZE] = '\0';
    memset(&mmer[1], 'T', MMER_SIZE-1);
    char compare_mmer[MMER_SIZE + 1];
    compare_mmer[MMER_SIZE] = '\0';
    int mmer_score = getscore(mmer);
    bool multiple_extension;
    char* a_key;
    int a_len;
    int score_limit = getbp('A')*MMER_SIZE;

    // iterate over all mmers from CTTT to AAAA..
    struct ZHashTable* mmer_hash, *compare_mmer_hash, *extend_table;
    struct ZHashEntry** kmer_entry, **compare_entry, **extend_entry = NULL;
    bool not_extended = false;
    while (mmer_score <= score_limit){
        // perform operation till mmer reaches AAAA..
        if ((mmer_hash = zhash_get(hash_table, mmer)) != NULL) {
            // iterate over all kmers of a particular mmer
            int array_index = 0;
            while (array_index < hash_sizes[mmer_hash->size_index]) {
                struct ZHashEntry** kmer_entry = &mmer_hash->entries[array_index];
                while (*kmer_entry != NULL) {
                    kmer_extension_node extension_node = find_kmer_extension(hash_table, *kmer_entry, mmer_score, forward);

                    if (extension_node.extend_entry != NULL) {
                        extend_entry = extension_node.extend_entry;
                        // create extension
                        more_kmer_extension_node further_extension =  extend_kmers(*kmer_entry, *extend_entry, forward);
                        // cannot delete both nodes directly as extend entry node points to kmer entry
                        if ((*extend_entry)->next == (*kmer_entry)) {
                            kmer_entry = extend_entry;
                            struct ZHashEntry* temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free extension node
                            temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free kmer node
                            mmer_hash->entry_count -= 2;
                        // cannot delete both nodes directly as kmer entry points to extend entry node
                        } else if ((*kmer_entry) == (*extend_entry)->next) {
                            struct ZHashEntry* temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free kmer node
                            temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free extension node
                            mmer_hash->entry_count -= 2;
                        // safe to delete
                        } else {
                            struct ZHashEntry* temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free kmer node
                            mmer_hash->entry_count--;
                            temp = *extend_entry;
                            *extend_entry = (*extend_entry)->next;
                            zfree_entry(temp, false); // free extension node
                            extension_node.extend_table->entry_count--;
                        }

                        // keep extending while possible
                        while (true) {
                            extension_node = more_kmer_extension(hash_table, further_extension.key, mmer_score, forward);
                            if (extension_node.extend_entry == NULL) {
                                break;
                            }

                            extend_entry = extension_node.extend_entry;
                            further_extension = further_extend_kmers(further_extension, *extend_entry, forward);
                            // extension node and kmer entry iterator are the same
                            if (*extend_entry == (*kmer_entry)) {
                                struct ZHashEntry* temp = *kmer_entry;
                                *kmer_entry = (*kmer_entry)->next;
                                zfree_entry(temp, false);
                            // extension node points to kmer entry iterator
                            } else if ((*extend_entry)->next == *kmer_entry) {
                                struct ZHashEntry* temp = *extend_entry;
                                kmer_entry = extend_entry;
                                *kmer_entry = (*extend_entry)->next;
                                zfree_entry(temp, false);
                            // kmer entry iterator points to extension node
                            } else {
                                struct ZHashEntry* temp = *extend_entry;
                                *extend_entry = (*extend_entry)->next;
                                free(temp);
                            }
                        }
                        // add further extended node to hash table
                        zhash_set(mmer_hash, further_extension.key, further_extension.read_id_lists);
                    } else {
                        kmer_entry = &(*kmer_entry)->next;
                    }   
                }
                array_index++;
            }
        }

        // get lexicographically next smallest mmer
        // mmer score is incremented
        mmer_score = next_smaller_mmer(mmer, mmer_score);
    }
}

void print_kmers(struct ZHashTable* hash_table) {
    struct ZHashTable* kmer_hash;
    struct ZHashEntry* mmer_entry, *kmer_entry;
    ll_node* read_id, *traverse;
    
    while ((mmer_entry = (struct ZHashEntry*) iterate_level_one_hash(hash_table, false, false)) != NULL) {
        kmer_hash = (mmer_entry)->val;
        printf("%s\n", mmer_entry->key); // print mmer
        // iterate over all kmers of mmer
        while ((kmer_entry = (struct ZHashEntry*) iterate_level_two_hash(kmer_hash, false, false)) != NULL) {
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
        }
        printf("\n");
    }
}

void expand_read_id_list(struct ZHashTable* hashtable) {
    struct ZHashEntry* mmer_entry = NULL, *kmer_entry = NULL;
    struct ZHashTable* kmer_hash = NULL;
    ll_node* read_id_list, *traverse = NULL;
    ll_node* read_id_lists;
    int kmer_len, i;
    while ((mmer_entry = iterate_level_one_hash(hashtable, false, false)) != NULL) {
        kmer_hash = mmer_entry->val;
        while ((kmer_entry = iterate_level_two_hash(kmer_hash, false, false)) != NULL) {
            traverse = NULL;
            read_id_list = kmer_entry->val;
            kmer_len = strlen(kmer_entry->key);
            for (i = 0; i < kmer_len; i++) {
                if (traverse == NULL) {
                    traverse = create_node_item(read_id_list);
                    read_id_lists = traverse;
                } else {
                    traverse->next = create_node_item(duplicate_llist(read_id_list));
                    traverse = traverse->next;
                }
            }
            kmer_entry->val = read_id_lists;
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
    int score, rev_score, max_score;
    bool is_rev;
    int msb;

    // slide a window of KMER_SIZE over read
    for (i = 0; i < read_len - KMER_SIZE + 1; i++) {

        // iterate over the read and calculate the signature from scratch
        if (kmer > signature) {

            // re intialize values each time for fresh calculation
            score = 0;
            rev_score = 0;
            max_score = 0;

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
            while (j < KMER_SIZE) {
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
            for (j = 0; j < MMER_SIZE; j++) {
                signature_cpy[j] = getbp(3 - getval(signature_cpy[j]));
            }

            for (j = 0; j < KMER_SIZE; j++) {
                kmer_key[j] = getbp(3 - getval(kmer_key[j]));
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
        if ((read_id_list = zhash_get(kmer_storage, kmer_key)) == NULL) {
            // create entry for the first time
            traverse = (ll_node*) create_node_num(read_id);
            zhash_set(kmer_storage, kmer_key, traverse);
        } else {
            // to make operation efficient and maintain descending order sorted linked list
            // shift read id of first node to second node and and put new read id in first node
            // all other nodes are untouched and there is no need to store the linked list again
            // as the pointer to first node has not changed
            traverse = (ll_node*) create_node_num(read_id_list->read_id);
            read_id_list->read_id = read_id;
            traverse->next = read_id_list->next;
            read_id_list->next = traverse;
        }

        // increment kmer pointer
        kmer++;
    }

    return hash_table;
}

// prune kmers which have base pair occuring in only one read
// such kmers are highly likely to have been generated by errors
struct ZHashTable* prune_kmers(struct ZHashTable* hash_table) {
    struct ZHashEntry** traverse, **to_remove = NULL;
    struct ZHashEntry* temp;
    ll_node* read_id_list, *temp_node;
    while ((traverse = iterate_level_two_hash(hash_table, true, false)) != NULL) {

        read_id_list = (ll_node*) (*traverse)->val;
        int count = 1;
        // check if number of reads exceeds cutoff
        while (read_id_list->next != NULL && count <= ABUNDANCE_CUTOFF) {
            count++;
            read_id_list = read_id_list->next;
        }

        if (count <= ABUNDANCE_CUTOFF) {
            // kmer has low occurence rate
            // free node and remove entry
            free((*traverse)->val);
            (*traverse)->val = NULL;
            // mark current node for removal
            iterate_level_two_hash(NULL, false, true);
        }
    }

    // if entire hash table is emptied free and return NULL
    if (hash_table->entry_count == 0) {
        free(hash_table);
        return NULL;
    } else {
        return hash_table;
    }
}

// iterates over all entries of first level hash table
// applies prune_kmers function on all non null entries
struct ZHashTable* prune_data(struct ZHashTable* hash_table) {
    struct ZHashEntry** traverse, *temp;
    while ((traverse = iterate_level_one_hash(hash_table, true, false)) != NULL) {
        // entries exist for this hash value
        if (prune_kmers((*traverse)->val) == NULL) {
            // hash table has been emptied remove entry
            // mark entry for removal
            (*traverse)->val = NULL;
            iterate_level_one_hash(NULL, false, true);
        }
    }
}

int main(int argc, char* argv[]) {
    // initialize file and structures
    FILE* file = fopen(argv[1], "r");
    struct ZHashTable* hash_table = zcreate_hash_table();
    
    // initialize variables
    char read[READ_LENGTH];
    int read_id = 0;
    
    // get all the reads from file
    // expecting read of length less than 20
    while (fgets(read, 50, file) != NULL) {

        // pre process and store read
        int len = strlen(read);
        read[--len] = '\0';

        process_read(hash_table, read, read_id++);
    }

    // prune stored values and remove possibly erroneous kmers
    prune_data(hash_table);
    // expand remaining entries
    expand_read_id_list(hash_table);

    // apply unitig extension to the data
    // first left to right directions
    // then in right to left direction
    find_kmer_extensions(hash_table, true);
    find_kmer_extensions(hash_table, false);

    // print kmers
    print_kmers(hash_table);
}
