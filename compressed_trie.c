#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #include "zhash.h"
#include "cstring.h"


typedef struct node node;
typedef struct node* Node;
typedef struct leafnode leafnode;
typedef struct leafnode* Leafnode;
typedef union basepairs basepairs;

union basepairs {
    unsigned char u;
    struct {
        __uint8_t bp0:2;
        __uint8_t bp1:2;
        __uint8_t bp2:2;
        __uint8_t bp3:2;
    };
};

typedef union next_ptr {
    Node child[4];
    Leafnode leaf;
} next_ptr;

struct leafnode {
    // struct ZHashTable set;
    int count;
    char* mmer;
    int mmer_index;
};

struct node {
    unsigned int size:4;
    unsigned int isleaf:2;
    basepairs bp;
    next_ptr next;
};

Node init() {
    Node root = (Node) malloc(sizeof(node));
    root->size = 0;
    root->bp.u = 0;
    root->next.child[0] = NULL;
    root->next.child[1] = NULL;
    root->next.child[2] = NULL;
    root->next.child[3] = NULL;
    root->isleaf = 0;
    return root;
}

int getval(char c){
	int i;
	switch(c){
		case 'A':
			i=0;
			break;
		case 'T':
			i=1;
			break;
		case 'G':
			i=2;
			break;
		case 'C':
			i=3;
			break;
		default :
			i=-1;
	}
	return i;
}

char getbp(__uint8_t bp){
	switch(bp){
		case 0:
			return 'A';
		case 1:
			return 'T';
		case 2:
			return 'G';
		case 3:
			return 'C';
	}
}

int checkcommon(char* str, int str_len, basepairs bp){
	int i;
	int size = str_len;
	if (size>=1 && bp.bp0==getval(str[0])){
		if(size>=2 &&bp.bp1==getval(str[1])){
			if(size>=3 &&bp.bp2==getval(str[2])){
				if(size>=4 &&bp.bp3==getval(str[3])){
					i=4;
				}else{
					i=3;
				}
			}else{
				i=2;
			}
		}else{
			i=1;
		}
	}else{
		i=0;
		printf("ERROR NOT POSSIBLE in checkcommon\n");
	}
	return i;
}

basepairs setbps(char* str, int str_len){
	basepairs bp;
	if(str_len > 0){
		bp.bp0 = getval(str[0]);
	}
	if(str_len > 1){
		bp.bp1 = getval(str[1]);
	}
	if(str_len > 2){
		bp.bp2 = getval(str[2]);
	}
	if(str_len > 3){
		bp.bp3 = getval(str[3]);
	}

	return bp;
}

basepairs reducebps(basepairs bp, int match){
	if (match==1) {
		bp.bp0=bp.bp1;
		bp.bp1=bp.bp2;
		bp.bp2=bp.bp3;
	} else if (match==2) {
		bp.bp0=bp.bp2;
		bp.bp1=bp.bp3;
	} else if (match==3) {
		bp.bp0=bp.bp3;
	}

	return bp;
}

void insert(Node root, Cstring str){
	Node parent = NULL, current = root;
    int i, j, match;

	do{
		i=getval(str->ptr[0]);
		if(i == -1) {
			printf("ERROR NOT POSSIBLE in insert\n");
		}
		if(current->next.child[i]){
			// next[i] != NULL
			parent = current;
			current = current->next.child[i];
			match = checkcommon(str->ptr, str->new_len, current->bp);
			if(match >= current->size) { 
				// follow current node
				erase_str(str, current->size);
				continue;
			} else {
				// part current node
				Node temp = (Node) malloc(sizeof(node));
				temp->size = match;
				temp->bp = setbps(str->ptr, match);
				temp->isleaf = 0;
				erase_str(str, match);

				for(j=0; j<4; j++) {
					temp->next.child[j]=NULL;
				}

				parent->next.child[i] = temp;
				parent = temp;
				current->size -= match;
				current->bp = reducebps(current->bp,match);
				temp->next.child[current->bp.bp0] = current;
				current=temp;
				// str.erase(0,current->size);
				//divide current node based on match;
			}
		}else{ 
			// next[i] == NULL
			current->next.child[i] = (Node) malloc(sizeof(node));
			parent = current;
			current = current ->next.child[i];
			int size = str->new_len;

			while(size>0){

				if(size>4){
					current-> size = 4;
					current-> bp = setbps(str->ptr, 4);
					current-> isleaf = 0;
					erase_str(str, 4);

					for(j=0;j<4;j++) {
						current->next.child[j]=NULL;
					}

					size -= 4;
					i = getval(str->ptr[0]);

					current->next.child[i] = (Node) malloc(sizeof(node));
					parent = current;
					current = current->next.child[i];
					continue;
				}else{
					current->size = size;
					current->bp = setbps(str->ptr, size);
					current->isleaf = 1;
					erase_str(str, size);
					current->next.leaf = (Leafnode) malloc(sizeof(leafnode));
					size=0;
					break;
					//leaf node
				}
			}
			break;
		}
	} while(str->new_len > 0);
	current->next.leaf->count++;
	return;
}

void printbps(basepairs bp, int size){
	if(size>0){
		printf("%c",getbp(bp.bp0));
	}
	if(size>1){
		printf("%c",getbp(bp.bp1));
	}
	if(size>2){
		printf("%c",getbp(bp.bp2));
	}
	if(size>3){
		printf("%c",getbp(bp.bp3));
	}
}

void print_trie(Node root,int tabs){
	int i;
	if(root->isleaf!=1){
		for(i=0;i<tabs;i++){
			printf("\t");
		}
		printbps(root->bp,root->size);
		printf("\n");
		for(i=0;i<4;i++)
			if(root->next.child[i])print_trie(root->next.child[i],tabs+1);
	}else{
		for(i=0;i<tabs;i++){
			printf("\t");
		}
		printbps(root->bp,root->size);
		printf(":%d\n",root->next.leaf->count);
	}
}

int main() {
    Cstring str;
	FILE* file;
    Node root = init();
    char* reads[] = {
		"CAGCCGCTGGGTCCG",
		"CGCGACCGGCTGGTG",
		"CTGGGGCAGGTCGGG",
		"GGGGAGCAGATCCGG",
		"CTCGGCCTGCCGCCC",
		"GGCCGGTGCACGCCG",
		"TCGCGCGGGCCAGCC",
		"GCTGGGTCCGCGCGA",
		"CCGGCTGGTGCTGGG",
		"GCAGGTCGGGGGGGA",
		"GCAGATCCGGCTCGG",
		"CCTGCCGCCCGGCCG",
		"GTGCACGCCGTCGCG",
		"CGGGCCAGCCGCTGG",
		"GTCCGCGCGACCGGC",
		"TGGTGCTGGGGCAGG",
		"TCGGGGGGGAGCAGA",
		"TCCGGCTCGGCCTGC",
		"CGCCCGGCCGGTGCA",
		"CGCCGTCGCGCGGGC",
	};

    int i, n = 20;
    for (i = 0; i < n; i++) {
        str = init_str(reads[i]);
        insert(root, str);
    }

    print_trie(root, -1);
    return 0;
}
