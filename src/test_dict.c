#include "header.h"


double test_hash(){
    // Example usage
    /*pair array[] = {{3, 4}, {5, 4}, {2, 2}};
    double value_array[7]; // log perms
    value_array[0] = 30;
    value_array[1] = 20.3;
    value_array[2] = 12.3;


    int table_size = 7;
    int table[7];
    int used_len = 3;

    for (int i = 0; i < table_size; i++) {
        table[i] = -1;  // Initialize all entries to -1 (empty)
    }
    for (int i = 0; i < used_len; i++) {
        int h = hash(array[i], table_size);
        while (table[h] >= 0) {
            h = (h + 1) % table_size;  // Linear probing to find the next empty slot
        }
        table[h] = i;
    }



    dictionary dict = {.value_array = value_array, .array = array, .table=table, .used_len = used_len, .table_size = table_size,};
  
*/

    dictionary * dict = init_dictionary(20);

    pair pair_ex = { .x=3, .y = 4 };
    add_to_dictionary(dict, pair_ex, 30);

    pair pair_ex2 = { .x=5, .y = 4 };
    add_to_dictionary(dict, pair_ex2, 2.3);

    pair pair_ex3 = { .x=2, .y = 2 };
    add_to_dictionary(dict, pair_ex3, 12.3);

    for (int i = 0; i < (*dict).used_len; ++i)
    {
        Rprintf("(%d, %d) = %f\n", (*dict).array[i].x, (*dict).array[i].y, (*dict).value_array[i]);
    }
    return 0;
    pair pair_ex4 = { .x=1, .y = 1 };
    double result = lookup_logperm(dict, pair_ex4);
    Rprintf("RESULT = %f\n", result);

    pair new_pair_ex = { .x=3, .y = 4 };
    result = lookup_logperm(dict, new_pair_ex);
    Rprintf("RESULT = %f\n", result);

    pair new_pair_ex2 = { .x=5, .y = 4 };
    result = lookup_logperm(dict, new_pair_ex2);
    Rprintf("RESULT = %f\n", result);

    pair new_pair_ex3 = { .x=4, .y = 10 };
    Rprintf("updating (4,10) with val 2\n");
    update_dict(new_pair_ex3, 2, dict);
    result = lookup_logperm(dict, new_pair_ex3);
    Rprintf("RESULT = %f\n", result);

    
    pair new_pair_ex4 = { .x=5, .y = 4 };
    Rprintf("updating (5,4) with val 2\n");
    update_dict(new_pair_ex4, 2, dict);
    update_dict(new_pair_ex4, 3, dict);
    result = lookup_logperm(dict, new_pair_ex4);
    Rprintf("RESULT = %f\n", result);



/*
    int n = sizeof(arr) / sizeof(pair);
    int table_size = 2 * n;  // Choose a table size that's at least twice the number of elements
    int* table = malloc(table_size * sizeof(int));
    for (int i = 0; i < table_size; i++) {
        table[i] = -1;  // Initialize all entries to -1 (empty)
    }
    for (int i = 0; i < n; i++) {
        int h = hash(arr[i], table_size);
        while (table[h] >= 0) {
            h = (h + 1) % table_size;  // Linear probing to find the next empty slot
        }
        table[h] = i;
    }
    pair pair_ex = { .x=1, .y = 1 };
    int index = lookup(arr, table, table_size, pair_ex);
    Rprintf("%d\n", index);  // Output: 2
    pair_ex.x=5;
    pair_ex.y=4;
    index = lookup(arr, table, table_size, pair_ex);
    Rprintf("%d\n", index);  // Output: 2
    free(table);*/
    return 0;
}
/*int test_hash() {
    // Example usage
    pair arr[] = {{3, 4}, {5, 4}, {2, 2}};
    int n = sizeof(arr) / sizeof(pair);
    int table_size = 2 * n;  // Choose a table size that's at least twice the number of elements
    int* table = malloc(table_size * sizeof(int));
    for (int i = 0; i < table_size; i++) {
        table[i] = -1;  // Initialize all entries to -1 (empty)
    }
    for (int i = 0; i < n; i++) {
        int h = hash(arr[i], table_size);
        while (table[h] >= 0) {
            h = (h + 1) % table_size;  // Linear probing to find the next empty slot
        }
        table[h] = i;
    }
    pair pair_ex = { .x=1, .y = 1 };
    int index = lookup(arr, table, table_size, pair_ex);
    Rprintf("%d\n", index);  // Output: 2
    pair_ex.x=5;
    pair_ex.y=4;
    index = lookup(arr, table, table_size, pair_ex);
    Rprintf("%d\n", index);  // Output: 2
    free(table);
    return 0;
}*/