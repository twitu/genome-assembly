CC=gcc
CFLAG=-g

binning: zhash.h zhash.c binning.c llist.c llist.h
	$(CC) $(CFLAG) zhash.c binning.c llist.c -o a.out
clean:
	rm -rf *o a.out
