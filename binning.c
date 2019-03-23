#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zhash.h"
#include "bucket_sort.h"
#include "llist.h"

#define MMER_SIZE 4 // efficient to keep mmer_size as powers of 2
#define KMER_SIZE 10

// in line helper functions
#define MIN(A,B) \
   ({ __typeof__ (A) _A = (A); \
       __typeof__ (B) _B = (B); \
     _A < _B ? _A : _B; })

#define MAX(A,B) \
   ({ __typeof__ (A) _A = (A); \
       __typeof__ (B) _B = (B); \
     _A > _B ? _A : _B; })

const int power_val[] = {1, 4, 16, 64, 256, 1024, 4096, 16384};


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

// stores all kmers of a read
// kmers are stored in 2 level hashing
// first level is hashed by signature of kmer
// second level is hashed by kmer itself
struct ZHashTable* process_read(struct ZHashTable* hash_table, char* read, int* read_id) {
    int read_len = strlen(read);
    char* kmer = read;
    char* signature = NULL;
    int i, j;

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
            // incrementally update scores
            j = MMER_SIZE;
            while (kmer[j] != '\0') {
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
            // store first MMER_SIZE characters in array
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
        signature_cpy[MMER_SIZE] = '\0';

        // get reverse complement of signature if rev complement has higher score
        if (is_rev) {
            for (i = 0; i < MMER_SIZE; i++) {
                signature_cpy[i] = getbp(3 - getval(signature_cpy[i]));
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
        // TODO: addition to list is O(n), has to be optimized
        ll_node* list;
        char temp;
        if ((list = zhash_get(kmer_storage, read)) == NULL) {
            list = add_ll_item(list, read_id);
            temp = kmer[KMER_SIZE];
            kmer[KMER_SIZE] = '\0';
            zhash_set(kmer_storage, kmer, list);
            kmer[KMER_SIZE] = temp;
        } else {
            add_ll_item(list, read_id);
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
        int* id = malloc(sizeof(int));
        *id = read_id++;

        process_read(hash_table, read, id);
    }
}
