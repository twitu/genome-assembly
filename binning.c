#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zhash.h"
#include "bucket_sort.h"

#define MMER_SIZE 4 // efficient to keep mmer_size as powers of 2


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

// possible sizes for hash table; must be prime numbers
static const size_t hash_sizes[] = {
  53, 101, 211, 503, 1553, 3407, 6803, 12503, 25013, 50261,
  104729, 250007, 500009, 1000003, 2000029, 4000037, 10000019,
  25000009, 50000047, 104395301, 217645177, 512927357, 1000000007
};

int main() {
    FILE* file = fopen("input.txt", "r");
    struct ZHashTable *hash_table;
    hash_table = zcreate_hash_table();
    
    char read[20];
    char signature[MMER_SIZE + 1];
    char mmer[MMER_SIZE + 1];
    int i;
    int pow_del = 64;
    int insert_point = 0;
    // get reads from file
    while (fgets(read, 20, file) != NULL) {
        // find signature for read
        int len = strlen(read);
        read[--len] = '\0';
        char* store = (char*) malloc(sizeof(char)*len+1);
        strncpy(store, read, len + 1);

        // initializer mmer counter
        for (i=0; i < MMER_SIZE; i++) {
            mmer[i] = 'Z';
        }
        mmer[MMER_SIZE] = '\0';
        signature[MMER_SIZE] = '\0';

        unsigned int mer_count = 0;
        int score = -85;
        int max_score = score;
        for (i=0; i < len; i++) {
            score -= getval(mmer[mer_count])*pow_del;
            score *= MMER_SIZE;
            score += getval(read[i]);
            mmer[mer_count] = read[i];
            if (score > max_score) {
                max_score = score;
                insert_point = mer_count + 1;
                strncpy(signature, mmer, MMER_SIZE);
            }
            mer_count++;
            if (mer_count==MMER_SIZE) {
                mer_count = 0;
            }
        }
        // TODO fix mapping by sorting
        rotating_sort(signature, MMER_SIZE, insert_point);
        zhash_set(hash_table, signature, store);
    }

    // iterate over zhash and print all stored reads
    int size, ii;
    size = hash_sizes[hash_table->size_index];

    for (ii = 0; ii < size; ii++) {
        struct ZHashEntry *entry;
        if ((entry = hash_table->entries[ii])) {
            printf("values for key %s\n", entry->key);
            while (entry) {
                printf("%s\n", (char*) entry->val);
                entry = entry->next;
            }
        };
    }
}