// iterates over reads and stores kmers

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zhash.h"
#include "llist.h"

#define MMER_SIZE 4        // efficient to keep mmer_size as powers of 2
#define KMER_SIZE 31       // fixed size of initial kmer extracted from reads
#define ABUNDANCE_CUTOFF 1 // kmer should occur in more reads than cutoff to avoid deletion
#define READ_LENGTH 101    // maximum size of read supported

// defined constants for faster multiplication
// MMER_SIZE should not exceed length of power_val
const int power_val[] = {1, 4, 16, 64, 256, 1024, 4096, 16384};

// possible sizes for hash table
static const size_t hash_sizes[] = {
    53, 101, 211, 503, 1553, 3407, 6803, 12503, 25013, 50261,
    104729, 250007, 500009, 1000003, 2000029, 4000037, 10000019,
    25000009, 50000047, 104395301, 217645177, 512927357, 1000000007};

typedef struct kmer_extension_node
{
    struct ZHashEntry **extend_entry;
    struct ZHashTable *extend_table;
} kmer_extension_node;

typedef struct more_kmer_extension_node
{
    char *key;
    ll_node *read_id_lists;
} more_kmer_extension_node;

/*******************************************
 * Helper Macros
*******************************************/

// returns minimum of A and B
#define MIN(A, B) \
    ({ __typeof__ (A) _A = (A); \
       __typeof__ (B) _B = (B); \
     _A < _B ? _A : _B; })

// returns maximum of A and B
#define MAX(A, B) \
    ({ __typeof__ (A) _A = (A); \
       __typeof__ (B) _B = (B); \
     _A > _B ? _A : _B; })

// swap values in variables x and y
#define SWAP(x, y)          \
    do                      \
    {                       \
        typeof(x) SWAP = x; \
        x = y;              \
        y = SWAP;           \
    } while (0)

/*****************************************
 * Helper functions for conversion between numeric and ascii value of base pair
 * Helper functions for calculating score of string and for returning next lexically smaller string
 * Higher ascii values have lower scores so higher numeric score indicates lower in dictionary order
*****************************************/

// converts numeric value of bp to its ascii value
char getbp(int bp)
{
    switch (bp)
    {
    case 0:
        return 'T';

    case 1:
        return 'G';

    case 2:
        return 'C';

    case 3:
        return 'A';

    default:
        return 'A';
    }
}

// converts ascii character to its numeric value
int getval(char c)
{
    switch (c)
    {
    case 'T':
        return 0;

    case 'G':
        return 1;
        
    case 'C':
        return 2;

    case 'A':
        return 3;

    default:
        return 3;
        
    }
}

// calculates numeric score of "string" by summing numeric scores of all characters in the string
int getscore(char *string)
{
    int score = 0;
    while (*string != '\0')
    {
        score = score * 4 + getval(*string);
        string++;
    }

    return score;
}

// returns score of next smaller mmer in dictionary order
// converts passed "mmer" string to next smaller mmer representation in dictionary order
// wraps around from AAAA to TTTT
int next_smaller_mmer(char *mmer, int mmer_score)
{
    for (int i = MMER_SIZE - 1; i >= 0; i--)
    {
        if (mmer[i] == 'A')
        {
            mmer[i] = 'T';
        }
        else
        {
            mmer[i] = getbp(getval(mmer[i]) + 1);
            break;
        }
    }

    return ++mmer_score;
}

/*****************************************
 * Functions for merging read id lists, keys and strings and kmers
*****************************************/

// return new list of read id lists where continuous range of KMER_SIZE - 1 nodes of a_node and b_node are merged
// forward direction merges right end of a_node with left end of b_node
// backward direction merges right end of b_node with left end of a_node
ll_node *merge_lists(int a_len, int b_len, ll_node *a_node, ll_node *b_node, bool forward)
{
    // swap for merging in backward direction
    if (!forward)
    {
        SWAP(a_node, b_node);
        SWAP(a_len, b_len);
    }

    // new_list points to starting list
    ll_node *new_list = a_node;
    int len = a_len;

    // skip read id lists that don't overlap
    for (int i = 0; i < len - (KMER_SIZE - 1); i++)
    {
        a_node = a_node->next;
    }

    // merge read id lists of two kmers for KMER_SIZE - 1 nodes that overlap
    // nodes of b_node are freed as their values are transfered to a_node nodes
    for (int i = 0; i < KMER_SIZE - 1; i++)
    {
        a_node->item = merge_sorted_list(a_node->item, b_node->item);

        ll_node *temp = b_node;
        b_node = b_node->next;
        free(temp);

        if (i == KMER_SIZE - 2)
        {
            // after merging overlapping nodes link rest of b_node nodes to new list
            a_node->next = b_node;
        }
        else
        {
            a_node = a_node->next;
        }
    }

    return new_list;
}

// returns of true if a_string and b_string overlap at continuous KMER_SIZE - 1 base pairs
// forward direction compares right end of a_string with left end of b_string
// backward direction compares right end of b_string with left end of a_string
bool compare_overlap(char *a_string, char *b_string, bool forward)
{
    // swap variables a_string and b_string when comparing backwards overlap
    if (!forward)
    {
        SWAP(a_string, b_string);
    }

    int len = strlen(a_string);
    for (int i = 0; i < KMER_SIZE - 1; i++)
    {
        if (a_string[len - (KMER_SIZE - 1) + i] != b_string[i])
        {
            return false;
        }
    }

    return true;
}

// returns merged key of a_key and b_key which overlap at continuous KMER_SIZE - 1 base pairs
// forward direction merges right end of a_key with left end of b_key
// backward direction merges right end of b_key with left end of a_key
char *merge_keys(int a_len, int b_len, char *a_key, char *b_key, bool forward)
{
    int len = a_len + b_len + 1 - (KMER_SIZE - 1);
    // Note: critical to use calloc, malloc can give unitialized string which can cause error
    char *new_key = calloc(len, sizeof(char));

    if (forward)
    {
        strncpy(new_key, a_key, a_len);
        strncpy(&new_key[a_len], &b_key[KMER_SIZE - 1], b_len - (KMER_SIZE - 1));
    }
    else
    {
        strncpy(new_key, b_key, b_len);
        strncpy(&new_key[b_len], &a_key[KMER_SIZE - 1], a_len - (KMER_SIZE - 1));
    }

    return new_key;
}

// extends to kmer entries pointed to by given hash entries in given direction
// returns node containing pointer to merged kmer and read id list
// Note: does not free given a and b hash table entries
more_kmer_extension_node extend_kmers(struct ZHashEntry *a, struct ZHashEntry *b, bool forward)
{
    int a_len = strlen(a->key);
    int b_len = strlen(b->key);

    // merge read ids of both entries and concatenate keys
    ll_node *new_read_ids = merge_lists(a_len, b_len, (ll_node *)a->val, (ll_node *)b->val, forward);
    char *new_key = merge_keys(a_len, b_len, (char *)a->key, (char *)b->key, forward);
    more_kmer_extension_node to_return;
    to_return.key = new_key;
    to_return.read_id_lists = new_read_ids;
    return to_return;
}

// takes partially extended kmer and extends with kmer in hash entry b
// returns node containing pointer to merged kmer and read id list
// Note: manages de-allocation of memory for previous key
more_kmer_extension_node further_extend_kmers(more_kmer_extension_node a, struct ZHashEntry *b, bool forward)
{
    int a_len = strlen(a.key);
    int b_len = strlen(b->key);

    // merge read ids of both entries and concatenate keys
    ll_node *new_read_ids = merge_lists(a_len, b_len, (ll_node *)a.read_id_lists, (ll_node *)b->val, forward);
    char *new_key = merge_keys(a_len, b_len, (char *)a.key, (char *)b->key, forward);

    free(a.key);
    a.key = new_key;
    a.read_id_lists = new_read_ids;
    return a;
}

/*****************************************
 * Functions for easy and efficient iteration of records from hash table
 * Maintains static pointers to for efficiency and supports in place deletion
 * Note: level one and level two functions two support nested iteration
*****************************************/

/**
 * Usage: 
 * 1. Iteration: returns entries in hash_table one by one and NULL when no other entries left
 * hash_table: pass the same pointer while iterating it
 * indirection: false returns a pointer to the entry, true returns a pointer to the pointer of the entry
 * indirection is useful when modifying entries
 * remove_current: don't care
 * if deletion was called previously current entry is removed and next entry is returned
 * 2. Deletion: marks current entry for removal which is removed when next iteration in called
 * hash_table: pass NULL
 * indirection: don't care
 * remove_current: true
 * Note: be careful with deletion when performing nested iteration on the same hash_table
 */
void *iterate_level_one_hash(struct ZHashTable *hash_table, bool indirection, bool remove_current)
{
    static struct ZHashTable *table = NULL;
    static struct ZHashEntry **entry;
    static int index;
    static bool remove;

    // remove currently pointed node
    // when the next iterate is called
    if (hash_table == NULL && remove_current)
    {
        remove = remove_current;
        return NULL;
    }

    // reset variables table does not match hash_table
    if (table != hash_table)
    {
        table = hash_table;
        entry = NULL;
        index = 0;
    }

    // entry points to previously returned value in the same chain
    // move to next entry in the chain
    if (entry != NULL && *entry != NULL)
    {
        if (!remove)
        {
            entry = &(*entry)->next;
        }
        else
        {
            // if current node is marked for removal handle differently
            struct ZHashEntry *temp = *entry;
            *entry = (*entry)->next;
            zfree_entry(temp, false);
            table->entry_count--;
            remove = false;
        }
    }

    // entry points to empty chain
    // move to next chain non empty chain
    if (entry == NULL || *entry == NULL)
    {
        while (index < hash_sizes[table->size_index])
        {
            if (table->entries[index] != NULL)
            {
                entry = &table->entries[index];
                index++;
                break;
            }
            index++;
        }
    }

    // if iteration has ended keep returning NULL for subsequent calls
    if (entry == NULL || *entry == NULL)
    {
        table = NULL;
        return NULL;
    }

    if (indirection)
    {
        return entry;
    }
    else
    {
        return *entry;
    }
}

/**
 * Usage: useful for nested iteration like iterating the hash_table for an mmer
 * 1. Iteration: returns entries in hash_table one by one and NULL when no other entries left
 * hash_table: pass the same pointer while iterating it
 * indirection: false returns a pointer to the entry, true returns a pointer to the pointer of the entry
 * indirection is useful when modifying entries
 * remove_current: don't care
 * if deletion was called previously current entry is removed and next entry is returned
 * 2. Deletion: marks current entry for removal which is removed when next iteration in called
 * hash_table: pass NULL
 * indirection: don't care
 * remove_current: true
 * Note: be careful with deletion when performing nested iteration on the same hash_table
 */
void *iterate_level_two_hash(struct ZHashTable *hash_table, bool indirection, bool remove_current)
{
    static struct ZHashTable *table = NULL;
    static struct ZHashEntry **entry;
    static int index;
    static bool remove;

    // remove currently pointed node
    // when the next iterate is called
    if (hash_table == NULL && remove_current)
    {
        remove = remove_current;
        return NULL;
    }

    // reset variables table does not match hash_table
    if (table != hash_table)
    {
        table = hash_table;
        entry = NULL;
        index = 0;
    }

    // entry points to previously returned value in the same chain
    // move to next entry in the chain
    if (entry != NULL && *entry != NULL)
    {
        if (!remove)
        {
            entry = &(*entry)->next;
        }
        else
        {
            // if current node is marked for removal handle differently
            struct ZHashEntry *temp = *entry;
            *entry = (*entry)->next;
            zfree_entry(temp, false);
            table->entry_count--;
            remove = false;
        }
    }

    // entry points to empty chain
    // move to next chain non empty chain
    if (entry == NULL || *entry == NULL)
    {
        while (index < hash_sizes[table->size_index])
        {
            if (table->entries[index] != NULL)
            {
                entry = &table->entries[index];
                index++;
                break;
            }
            index++;
        }
    }

    // if iteration has ended keep returning NULL for subsequent calls
    if (entry == NULL || *entry == NULL)
    {
        table = NULL;
        return NULL;
    }

    if (indirection)
    {
        return entry;
    }
    else
    {
        return *entry;
    }
}

/*****************************************
 * Find possible kmer extensions and extend
*****************************************/

/**
 * Usage:
 * returns kmer information that overlaps at KMER_SIZE - 1 base pairs with given kmer entry key
 * returns entry and table of candidate kmer, NULL values of entry and table if no candidate is found
 * to be used when finding first extension for kmer
 * Arguments:
 * hash_table: pass mmer hashtable
 * entry: entry should contain kmer information that is to be extended
 * mmer_score: mmer_score of the kmer that is to be extended
 * forward: true to return candidate for right end extension and false to for left end extension
 */
kmer_extension_node find_kmer_extension(struct ZHashTable *hash_table, struct ZHashEntry *entry, int mmer_score, bool forward)
{
    // initialize key structure
    char *key = entry->key;
    int key_len = strlen(key);

    // initialize signature mmer
    char compare_mmer[MMER_SIZE + 1];
    compare_mmer[MMER_SIZE] = '\0';
    if (forward) {
        // forward direction generates possible mmers from right end of kmer
        strncpy(compare_mmer, &key[key_len - (MMER_SIZE - 1)], MMER_SIZE - 1);
    } else {
        // backward direction generated possible mmers from left end of kmer
        strncpy(&compare_mmer[1], key, MMER_SIZE - 1);
    }

    bool multiple_extension = false;
    struct ZHashEntry **extend_entry = NULL, **compare_entry = NULL;
    struct ZHashTable *compare_mmer_hash = NULL, *extend_table = NULL;

    // iterate over 4 possible compare mmers that can be created by changing one base pair
    // on the right end for forward direction or the left end for backward direction
    for (int i = 0; i < 4; i++)
    {
        if (forward) {
            compare_mmer[MMER_SIZE - 1] = getbp(i);
        } else {
            compare_mmer[0] = getbp(i);
        }

        if (getscore(compare_mmer) > mmer_score)
        {
            // extension only with lexicographically larger mmers
            continue;
        }

        if ((compare_mmer_hash = (struct ZHashTable *)zhash_get(hash_table, compare_mmer)) == NULL)
        {
            // ignore if mmer does not have entry
            continue;
        }

        // compare signature is lexicographically greater than or equal
        while ((compare_entry = iterate_level_two_hash(compare_mmer_hash, true, false)) != NULL)
        {
            // handle equality case
            if (*compare_entry == entry)
                continue;

            // check if KMER_SIZE - 1 characters overlap
            if (!compare_overlap(key, (*compare_entry)->key, forward))
                continue;

            // if extension entry already exists
            // there are multiple possible extensions
            // unitig extension is not possible
            if (extend_entry != NULL)
            {
                extend_entry = NULL;
                extend_table = NULL;
                multiple_extension = true;
                break;
            }
            else
            {
                extend_table = compare_mmer_hash;
                extend_entry = compare_entry;
            }
        }

        if (multiple_extension)
        {
            break;
        }
    }

    // deduct count from entry table
    kmer_extension_node to_return;
    to_return.extend_entry = extend_entry;
    to_return.extend_table = extend_table;
    return to_return;
}

/**
 * Usage: 
 * returns kmer infromation that overlaps at KMER_SIZE - 1 base pairs with given kmer string
 * returns entry and table of candidate kmer, NULL values of entry and table if no candidate is found
 * to be used when one extension has already been performed on a kmer
 * Arguments:
 * hash_table: pass mmer hashtable
 * key: pass kmer string that is to be extended
 * mmer_score: mmer_score of the kmer that is to be extended
 * forward: true to return candidate for right end extension and false to for left end extension
 */
kmer_extension_node more_kmer_extension(struct ZHashTable *hash_table, char *key, int mmer_score, bool forward)
{
    // initialize key structure
    int key_len = strlen(key);

    // initialize signature mmer
    char compare_mmer[MMER_SIZE + 1];
    compare_mmer[MMER_SIZE] = '\0';
    if (forward) {
        // forward direction generates possible mmers from right end of kmer
        strncpy(compare_mmer, &key[key_len - (MMER_SIZE - 1)], MMER_SIZE - 1);
    } else {
        // backward direction generated possible mmers from left end of kmer
        strncpy(&compare_mmer[1], key, MMER_SIZE - 1);
    }

    bool multiple_extension = false;
    struct ZHashEntry **extend_entry = NULL, **compare_entry = NULL;
    struct ZHashTable *compare_mmer_hash = NULL, *extend_table = NULL;

    // iterate over 4 possible compare mmers that can be created by changing one base pair
    // on the right end for forward direction or the left end for backward direction
    for (int i = 0; i < 4; i++)
    {
        if (forward) {
            compare_mmer[MMER_SIZE - 1] = getbp(i);
        } else {
            compare_mmer[0] = getbp(i);
        }

        if (getscore(compare_mmer) > mmer_score)
        {
            // extension only with lexicographically larger mmers
            continue;
        }

        if ((compare_mmer_hash = (struct ZHashTable *)zhash_get(hash_table, compare_mmer)) == NULL)
        {
            // ignore if mmer does not have entry
            continue;
        }

        // compare signature is lexicographically greater than or equal
        while ((compare_entry = iterate_level_two_hash(compare_mmer_hash, true, false)) != NULL)
        {
            // check if KMER_SIZE - 1 characters overlap
            if (!compare_overlap(key, (*compare_entry)->key, forward))
                continue;

            // if extension entry already exists
            // there are multiple possible extensions
            // unitig extension is not possible
            if (extend_entry != NULL)
            {
                extend_entry = NULL;
                extend_table = NULL;
                multiple_extension = true;
                break;
            }
            else
            {
                extend_table = compare_mmer_hash;
                extend_entry = compare_entry;
            }
        }

        if (multiple_extension)
        {
            break;
        }
    }

    // deduct count from entry table
    kmer_extension_node to_return;
    to_return.extend_entry = extend_entry;
    to_return.extend_table = extend_table;
    return to_return;
}

/**
 * Usage:
 * find extensions for all kmers in all hash tables and perform extension
 * Automatically hash table entries that have been merged and creates new ones
 * Arguments:
 * hash_table: pass mmer hash table
 * forward: true for right end extension and false for left end extension
 */
void find_kmer_extensions(struct ZHashTable *hash_table, bool forward)
{
    // initialize signature kmer
    char mmer[MMER_SIZE + 1];
    mmer[0] = 'C';
    mmer[MMER_SIZE] = '\0';
    memset(&mmer[1], 'T', MMER_SIZE - 1);
    char compare_mmer[MMER_SIZE + 1];
    compare_mmer[MMER_SIZE] = '\0';
    int mmer_score = getscore(mmer);
    bool multiple_extension;
    char *a_key;
    int a_len;
    int score_limit = getbp('A') * MMER_SIZE;

    // iterate over all mmers from CTTT to AAAA..
    struct ZHashTable *mmer_hash, *compare_mmer_hash, *extend_table;
    struct ZHashEntry **kmer_entry, **compare_entry, **extend_entry = NULL;
    bool not_extended = false;
    while (mmer_score <= score_limit)
    {
        // perform operation till mmer reaches AAAA..
        if ((mmer_hash = zhash_get(hash_table, mmer)) != NULL)
        {
            // iterate over all kmers of a particular mmer
            int array_index = 0;
            while (array_index < hash_sizes[mmer_hash->size_index])
            {
                struct ZHashEntry **kmer_entry = &mmer_hash->entries[array_index];
                while (*kmer_entry != NULL)
                {
                    kmer_extension_node extension_node = find_kmer_extension(hash_table, *kmer_entry, mmer_score, forward);

                    if (extension_node.extend_entry != NULL)
                    {
                        // create first extension
                        extend_entry = extension_node.extend_entry;
                        more_kmer_extension_node further_extension = extend_kmers(*kmer_entry, *extend_entry, forward);
                        // cannot delete both nodes directly as extend entry node points to kmer entry
                        if ((*extend_entry)->next == (*kmer_entry))
                        {
                            kmer_entry = extend_entry;
                            struct ZHashEntry *temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free extension node
                            temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free kmer node
                            mmer_hash->entry_count -= 2;
                        }
                        // cannot delete both nodes directly as kmer entry points to extend entry node
                        else if ((*kmer_entry) == (*extend_entry)->next)
                        {
                            struct ZHashEntry *temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free kmer node
                            temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free extension node
                            mmer_hash->entry_count -= 2;
                            // safe to delete
                        }
                        else
                        {
                            struct ZHashEntry *temp = *kmer_entry;
                            *kmer_entry = (*kmer_entry)->next;
                            zfree_entry(temp, false); // free kmer node
                            mmer_hash->entry_count--;
                            temp = *extend_entry;
                            *extend_entry = (*extend_entry)->next;
                            zfree_entry(temp, false); // free extension node
                            extension_node.extend_table->entry_count--;
                        }

                        // keep extending while possible
                        while (true)
                        {
                            extension_node = more_kmer_extension(hash_table, further_extension.key, mmer_score, forward);
                            if (extension_node.extend_entry == NULL)
                            {
                                break;
                            }

                            extend_entry = extension_node.extend_entry;
                            further_extension = further_extend_kmers(further_extension, *extend_entry, forward);
                            // extension node and kmer entry iterator are the same
                            if (*extend_entry == (*kmer_entry))
                            {
                                struct ZHashEntry *temp = *kmer_entry;
                                *kmer_entry = (*kmer_entry)->next;
                                zfree_entry(temp, false);
                            }
                            // extension node points to kmer entry iterator
                            else if ((*extend_entry)->next == *kmer_entry)
                            {
                                struct ZHashEntry *temp = *extend_entry;
                                kmer_entry = extend_entry;
                                *kmer_entry = (*extend_entry)->next;
                                zfree_entry(temp, false);
                            }
                            // kmer entry iterator points to extension node
                            else
                            {
                                struct ZHashEntry *temp = *extend_entry;
                                *extend_entry = (*extend_entry)->next;
                                free(temp);
                            }
                        }
                        // add further extended node to hash table
                        zhash_set(mmer_hash, further_extension.key, further_extension.read_id_lists);
                    }
                    else
                    {
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

/*****************************************
 * Debugger functions for printing information in mmer hash table
 * TODO: Add more functions to print more information
*****************************************/

// Usage: prints all mmers along with their kmers and reads ids
// Arguments: pass mmer hash table
void print_kmer_read_ids(struct ZHashTable *hash_table)
{
    struct ZHashTable *kmer_hash;
    struct ZHashEntry *mmer_entry, *kmer_entry;
    ll_node *read_id, *traverse;

    while ((mmer_entry = (struct ZHashEntry *)iterate_level_one_hash(hash_table, false, false)) != NULL)
    {
        kmer_hash = (mmer_entry)->val;
        printf("%s\n", mmer_entry->key); // print mmer
        // iterate over all kmers of mmer
        while ((kmer_entry = (struct ZHashEntry *)iterate_level_two_hash(kmer_hash, false, false)) != NULL)
        {
            printf("%s\n", kmer_entry->key);
            read_id = (ll_node *)kmer_entry->val;
            // iterate over read id lists of each base pair
            while (read_id != NULL)
            {
                traverse = (ll_node *)read_id->item;
                // print each read id of base pair in same line
                while (traverse != NULL)
                {
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

// Usage: prints all kmers
// Arguments: pass mmer hash table
void print_kmers(struct ZHashTable *hash_table)
{
    struct ZHashTable *kmer_hash;
    struct ZHashEntry *mmer_entry, *kmer_entry;
    ll_node *read_id, *traverse;

    // iterate over all mmers in hash table
    while ((mmer_entry = (struct ZHashEntry *)iterate_level_one_hash(hash_table, false, false)) != NULL)
    {
        kmer_hash = (mmer_entry)->val;
        // iterate over all kmers of mmer
        while ((kmer_entry = (struct ZHashEntry *)iterate_level_two_hash(kmer_hash, false, false)) != NULL)
        {
            printf("%s\n", kmer_entry->key);
        }
    }
}

/*****************************************
 * Process reads to populate kmer and mmer hash table
 * Perform pruning and deletion of low abundance kmers
*****************************************/

/**
 * Usage:
 * duplicate read id list for each base pair
 * to be called after pruning so that only abundant kmers have read ids expanded
 * Arguments:
 * pass mmer hash table
 */
void expand_read_id_list(struct ZHashTable *hashtable)
{
    struct ZHashEntry *mmer_entry = NULL, *kmer_entry = NULL;
    struct ZHashTable *kmer_hash = NULL;
    ll_node *read_id_list, *traverse = NULL;
    ll_node *read_id_lists;
    int kmer_len, i;
    while ((mmer_entry = iterate_level_one_hash(hashtable, false, false)) != NULL)
    {
        kmer_hash = mmer_entry->val;
        while ((kmer_entry = iterate_level_two_hash(kmer_hash, false, false)) != NULL)
        {
            traverse = NULL;
            read_id_list = kmer_entry->val;
            kmer_len = strlen(kmer_entry->key);
            for (i = 0; i < kmer_len; i++)
            {
                if (traverse == NULL)
                {
                    traverse = create_node_item(read_id_list);
                    read_id_lists = traverse;
                }
                else
                {
                    traverse->next = create_node_item(duplicate_llist(read_id_list));
                    traverse = traverse->next;
                }
            }
            kmer_entry->val = read_id_lists;
        }
    }
}

/**
 * Usage:
 * stores all kmers of a read
 * kmers are stored in 2 level hashing
 * first level is hashed by signature of kmer i.e. a mmer
 * second level is hashed by kmer itself
 * lexically smaller of kmer and its reverse complement is stored
 * Arguments:
 * hash_table: mmer hash table
 * read: read from which kmers are be parsed and stored
 * read_id: for debugging purposes
 */
struct ZHashTable *process_read(struct ZHashTable *hash_table, char *read, int read_id)
{
    int read_len = strlen(read);
    char *kmer = read;
    char *signature = NULL;
    int i, j;

    // initialize local variables for extracting signature
    char kmer_key[KMER_SIZE + 1];
    char mmer[MMER_SIZE + 1];
    char signature_cpy[MMER_SIZE + 1];
    int score, rev_score, max_score;
    bool is_rev;
    int msb;

    // slide a window of KMER_SIZE over read
    for (i = 0; i < read_len - KMER_SIZE + 1; i++)
    {

        // iterate over the read and calculate the signature from scratch
        if (kmer > signature)
        {

            // re intialize values each time for fresh calculation
            score = 0;
            rev_score = 0;
            max_score = 0;

            // store first MMER_SIZE characters in array
            for (j = 0; j < MMER_SIZE; j++)
            {
                mmer[j] = kmer[j];
                score = score * 4 + getval(kmer[j]);
                rev_score = rev_score * 4 + 3 - getval(kmer[j]);
            }
            mmer[MMER_SIZE] = '\0';

            // initialize max score
            if (score > rev_score)
            {
                max_score = score;
                is_rev = false;
            }
            else
            {
                max_score = rev_score;
                is_rev = true;
            }
            signature = kmer;
            msb = 0;

            // iterate over other characters
            j = MMER_SIZE;
            while (j < KMER_SIZE)
            {
                // incrementally update scores
                score = (score - getval(mmer[msb]) * power_val[MMER_SIZE - 1]) * 4;
                score += getval(kmer[j]);
                rev_score = (rev_score - ((3 - getval(mmer[msb])) * power_val[MMER_SIZE - 1])) * 4;
                rev_score += 3 - getval(kmer[j]);

                // rotate msb of mmer
                mmer[msb] = kmer[j];
                msb = (msb + 1) % MMER_SIZE;

                // increment pointer to next char
                j++;

                // change max score if required
                // set boolean if rev complement of mmer is signature
                if (MAX(score, rev_score) > max_score)
                {
                    if (score > rev_score)
                    {
                        max_score = score;
                        is_rev = false;
                    }
                    else
                    {
                        max_score = rev_score;
                        is_rev = true;
                    }

                    // set new signature to point MMER_SIZE behind the next char to be added
                    signature = &kmer[j] - MMER_SIZE;
                }
            }
        }
        // No need for calculating signature from scratch
        // check if new signatures are smaller than previous signature
        else
        {
            // compare current signature with new mmer
            // new mmer is created by last letter added to the current kmer
            // store last MMER_SIZE characters of kmer in mmer array
            for (j = KMER_SIZE - MMER_SIZE; j < MMER_SIZE; j++)
            {
                mmer[j] = kmer[j];
                score = score * 4 + getval(kmer[j]);
                rev_score = rev_score * 4 + 3 - getval(kmer[j]);
            }
            mmer[MMER_SIZE] = '\0';

            if (MAX(score, rev_score) > max_score)
            {
                if (score > rev_score)
                {
                    max_score = score;
                    is_rev = false;
                }
                else
                {
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
        if (is_rev)
        {
            for (j = 0; j < MMER_SIZE; j++)
            {
                signature_cpy[j] = getbp(3 - getval(signature_cpy[j]));
            }

            for (j = 0; j < KMER_SIZE; j++)
            {
                kmer_key[j] = getbp(3 - getval(kmer_key[j]));
            }
        }

        // check if this mmer has been stored before
        // if not create a new hash table to store kmers for this signature
        struct ZHashTable *kmer_storage;
        if ((kmer_storage = zhash_get(hash_table, signature_cpy)) == NULL)
        {
            kmer_storage = zcreate_hash_table();
            zhash_set(hash_table, signature_cpy, kmer_storage);
        }

        // check if this kmer has been stored previously
        ll_node *read_id_list, *traverse;
        if ((read_id_list = zhash_get(kmer_storage, kmer_key)) == NULL)
        {
            // create entry for the first time
            traverse = (ll_node *)create_node_num(read_id);
            zhash_set(kmer_storage, kmer_key, traverse);
        }
        else
        {
            // to make operation efficient and maintain descending order sorted linked list
            // shift read id of first node to second node and and put new read id in first node
            // all other nodes are untouched and there is no need to store the linked list again
            // as the pointer to first node has not changed
            traverse = (ll_node *)create_node_num(read_id_list->read_id);
            read_id_list->read_id = read_id;
            traverse->next = read_id_list->next;
            read_id_list->next = traverse;
        }

        // increment kmer pointer
        kmer++;
    }

    return hash_table;
}

/**
 * Usage:
 * delete kmers that don't occur in more than ABUNDANCE_CUTOFF number of reads
 * such kmers are highly likely to have been generated by errors
 * returns NULL if all kmers in the hash table are freed
 * Arguments: pass kmer hash table
 */
struct ZHashTable *prune_kmers(struct ZHashTable *hash_table)
{
    struct ZHashEntry **traverse, **to_remove = NULL;
    struct ZHashEntry *temp;
    ll_node *read_id_list, *temp_node;
    while ((traverse = iterate_level_two_hash(hash_table, true, false)) != NULL)
    {

        read_id_list = (ll_node *)(*traverse)->val;
        int count = 1;
        // check if number of reads exceeds cutoff
        while (read_id_list->next != NULL && count <= ABUNDANCE_CUTOFF)
        {
            count++;
            read_id_list = read_id_list->next;
        }

        if (count <= ABUNDANCE_CUTOFF)
        {
            // kmer has low occurence rate
            // free node and remove entry
            free((*traverse)->val);
            (*traverse)->val = NULL;
            // mark current node for removal
            iterate_level_two_hash(NULL, false, true);
        }
    }

    // if entire hash table is emptied free and return NULL
    if (hash_table->entry_count == 0)
    {
        free(hash_table);
        return NULL;
    }
    else
    {
        return hash_table;
    }
}

/**
 * Usage:
 * delete all kmers that don't occur in more than ABUNDANCE_CUTOFF number of reads
 * Arguments: pass mmer hash table
 */
struct ZHashTable *prune_data(struct ZHashTable *hash_table)
{
    struct ZHashEntry **traverse, *temp;
    while ((traverse = iterate_level_one_hash(hash_table, true, false)) != NULL)
    {
        // entries exist for this hash value
        if (prune_kmers((*traverse)->val) == NULL)
        {
            // hash table has been emptied remove entry
            // mark entry for removal
            (*traverse)->val = NULL;
            iterate_level_one_hash(NULL, false, true);
        }
    }
}

// pass file name containing reads
int main(int argc, char *argv[])
{
    // initialize file and structures
    FILE *file = fopen(argv[1], "r");
    struct ZHashTable *hash_table = zcreate_hash_table();

    // initialize variables
    char read[READ_LENGTH];
    int read_id = 0;

    // get all the reads from file
    while (fgets(read, READ_LENGTH, file) != NULL)
    {

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
