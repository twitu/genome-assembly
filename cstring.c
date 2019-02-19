#include <stdlib.h>
#include <string.h>

#include "cstring.h"

Cstring init_str(char* str) {
    Cstring str_obj = (Cstring) malloc(sizeof(cstring));
    str_obj->str = str;
    str_obj->ptr = str;
    str_obj->index = 0;
    str_obj->len = strlen(str);
    str_obj->new_len = str_obj->len;
    return str_obj;
}

Cstring erase_str(Cstring str, int offset) {
    if (str->index + offset > str->len) return NULL;
    str->index += offset;
    str->str = &str->str[str->index];
    str->new_len -= offset;
    return str;
}
