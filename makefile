CC=gcc
CFLAG=-g
compile: cstring.c cstring.h compressed_trie.c
	$(CC) $(CFLAG) cstring.c compressed_trie.c -o a.out
binning: zhash.h zhash.c binning.c bucket_sort.h bucket_sort.c
	$(CC) $(CFLAG) zhash.c binning.c bucket_sort.c -o a.out
clean:
	rm -rf *o a.out
