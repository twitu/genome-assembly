compile: cstring.c cstring.h compressed_trie.c
	gcc -g cstring.c compressed_trie.c -o a.out
binning: zhash.h zhash.c binning.c bucket_sort.h bucket_sort.c
	gcc -g zhash.c binning.c bucket_sort.c -o a.out