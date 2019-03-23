CC=gcc
CFLAG=-g
compile: cstring.c cstring.h compressed_trie.c
	$(CC) $(CFLAG) cstring.c compressed_trie.c -o a.out
binning: zhash.h zhash.c binning.c llist.c llist.h
	$(CC) $(CFLAG) zhash.c binning.c llist.c -o a.out
clean:
	rm -rf *o a.out
