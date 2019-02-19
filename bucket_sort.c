#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int get_char_val(char c){
	switch(c){
		case 'T':
			return 3;
		case 'G':
			return 2;
		case 'C':
			return 1;
		case 'A':
			return 0;
		default :
			return -1;
	}
}

char get_char_bp(__uint8_t bp){
	switch(bp){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
}

void bucket_sort(char string[], int len) {
    int buckets[4] = {0, 0, 0, 0};
    int i;
    for (i = 0; i < len; i++) {
        buckets[get_char_val
    (string[i])]++;
    }

    int j = 0, count = 0;
    for (i = 0; i < 4; i++) {
        char c = get_char_bp(i);
        for (j = 0; j < buckets[i]; j++) {
            string[count++] = c;
        }
    }
}

void rotating_sort(char string[], int len, int pointer) {
    int i;
    int count = 0;
    char* sorted = (char*) malloc(sizeof(char)*(len+1));
    sorted[len] = '\0';

    for (i = pointer; i < len; i++) {
        sorted[count++] = string[i];
    }

    for (i = 0; i < pointer; i++) {
        sorted[count++] = string[i];
    }

    strncpy(string, sorted, len);
    free(sorted);
}
