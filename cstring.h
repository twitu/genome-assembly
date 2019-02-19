#ifndef CSTRING_H
#define CSTRING_H

typedef struct cstring {
    char* str;
    char* ptr;
    int len;
    int index;
    int new_len;
} cstring;

typedef cstring* Cstring;

Cstring init_str (char* str);
Cstring erase_str (Cstring str, int offset);

#endif