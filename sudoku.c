#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>
#include <ctype.h>
#include <errno.h>
#include "minheap.h"

#define INT_TYPE unsigned long long 
#define INT_TYPE_SIZE (sizeof(INT_TYPE) * 8)
#define CELL_VAL_SIZE 1
#define MAX_BDIM 8
#define TAG_SIZE 1
#define TAG_NEW_JOB 2
#define TAG_ADD_NODES 3
#define TAG_SOLVED 4
#define TAG_TERMINATED 5

enum SOLVE_STRATEGY {SUDOKU_SOLVE, SUDOKU_COUNT_SOLS};

typedef enum {
    STR2INT_SUCCESS,
    STR2INT_OVERFLOW,
    STR2INT_UNDERFLOW,
    STR2INT_INCONVERTIBLE
} str2int_errno;

#ifndef SUDOKU_SOLVE_STRATEGY
    #define SUDOKU_SOLVE_STRATEGY SUDOKU_SOLVE
#endif

#define BUILD_ERROR_IF(condition) ((void)sizeof(char[1 - 2*!!(condition)]))

int mpi_rank, mpi_size, max_workload, min_workload, max_nodes_per_message, coord_send_nodes_amount, n_threads, stop_threads = 0;
dynarray *node_stack;
MPI_Request solved_request;
pthread_mutex_t solved_lock;
struct timeval start_time, end_time;

void BUILD_TIME_CHECKS()
{
    BUILD_ERROR_IF(INT_TYPE_SIZE * CELL_VAL_SIZE < MAX_BDIM * MAX_BDIM);
}

typedef struct cellval
{
    INT_TYPE v[CELL_VAL_SIZE];
}cell_v;

typedef struct cell_coord
{
    int r, c;
}cell_coord;

typedef struct sudoku
{
    int bdim;
    int dim;
    int peers_size;
    char* grid;    
    int grid_size;
    cell_coord ****unit_list;
    cell_coord ***peers;
    cell_v **values;   
    unsigned long long sol_count;
}sudoku;

typedef struct Node
{
    cell_v **values;
    int dim;
    char sqrI;
    char sqrJ;
    char digits_count;
    char digit;
    unsigned int value;
}Node;

typedef struct ThreadArgs
{
    sudoku *s;
    int workload;
    int ret;
}ThreadArgs;

static int assign (sudoku *s, int i, int j, int d);

static inline int cell_v_get(cell_v *v, int p)
{
    return !!((*v).v[(p - 1) / INT_TYPE_SIZE] & (((INT_TYPE)1) << ((p - 1) % INT_TYPE_SIZE)));
}

static inline void cell_v_unset(cell_v *v, int p)
{
    (*v).v[(p - 1) / INT_TYPE_SIZE] &= ~(((INT_TYPE)1) << ((p - 1) % INT_TYPE_SIZE));
}

static inline void cell_v_set(cell_v *v, int p)
{
    (*v).v[(p - 1) / INT_TYPE_SIZE] |= ((INT_TYPE)1) << ((p -1) % INT_TYPE_SIZE);
}

static inline int cell_v_count(cell_v *v)
{
    int acc = 0;

    for(int i = 0; i < CELL_VAL_SIZE; i++)
    {
        acc += __builtin_popcountll((*v).v[i]);
    }

    return acc;
}

static inline int digit_get (cell_v *v)
{
    int count = cell_v_count(v);

    if(count != 1)
    {
        return -1;
    }

    for(int i = 0; i < CELL_VAL_SIZE; i++)
    { 
        if((*v).v[i])
        {
            return 1 + INT_TYPE_SIZE * i + __builtin_ctzll((*v).v[i]);
        }
    }

    return -1;
}

static void destroy_sudoku(sudoku *s)
{
    for(int i = 0; i < s->dim; i++)
    {
        for(int j = 0; j < s->dim; j++)
        {
            for(int k = 0; k < 3; k++)
            {
                free(s->unit_list[i][j][k]);
            }

            free(s->unit_list[i][j]);
        }

        free(s->unit_list[i]);
    }

    free(s->unit_list);
    
    for(int i = 0; i < s->dim; i++)
    {
        for(int j = 0; j < s->dim; j++)
        {
            free(s->peers[i][j]);
        }

        free(s->peers[i]);
    }

    free(s->peers);
    
    for(int i = 0; i < s->dim; i++)
    {
        free(s->values[i]);
    }

    free(s->values);
    free(s);
}

static void init(sudoku *s)
{
    int i, j, k, l, pos;
    
    for(i = 0; i < s->dim; i++)
    {
        int ibase = i / s->bdim * s->bdim;

        for(j = 0; j < s->dim; j++)
        {
            for(pos = 0; pos < s->dim; pos++) 
            {
                s->unit_list[i][j][0][pos].r = i;
                s->unit_list[i][j][0][pos].c = pos;
                s->unit_list[i][j][1][pos].r = pos;
                s->unit_list[i][j][1][pos].c = j;
            }

            int jbase = j / s->bdim * s->bdim;

            for(pos = 0, k = 0; k < s->bdim; k++)
            {
                for(l = 0; l < s->bdim; l++, pos++)
                {
                    s->unit_list[i][j][2][pos].r = ibase + k;
                    s->unit_list[i][j][2][pos].c = jbase + l;
                }
            }
        }
    }
    
    for(i = 0; i < s->dim; i++)
    {
        for(j = 0; j < s->dim; j++)
        {
            pos = 0;

            for(k = 0; k < s->dim; k++)
            {
                if(s->unit_list[i][j][0][k].c != j)
                {
                    s->peers[i][j][pos++] = s->unit_list[i][j][0][k]; 
                }
            }

            for(k = 0; k < s->dim; k++)
            { 
                cell_coord sq = s->unit_list[i][j][1][k];

                if(sq.r != i)
                {
                    s->peers[i][j][pos++] = sq; 
                }

                sq = s->unit_list[i][j][2][k];

                if (sq.r != i && sq.c != j)
                {
                    s->peers[i][j][pos++] = sq; 
                }
            }
        }
    }

    assert(pos == s->peers_size);
}

static int parse_grid(sudoku *s)
{
    int i, j, k;
    int ld_vals[s->dim][s->dim];

    for(k = 0, i = 0; i < s->dim; i++)
    {
        for(j = 0; j < s->dim; j++, k++)
        {
            ld_vals[i][j] = s->grid[k];
        }
    }

    for(i = 0; i < s->dim; i++)
    {
        for(j = 0; j < s->dim; j++)
        {
            for(k = 1; k <= s->dim; k++)
            {
                cell_v_set(&s->values[i][j], k);
            }
        }
    }
    
    for(i = 0; i < s->dim; i++)
    {
        for(j = 0; j < s->dim; j++)
        {
            if(ld_vals[i][j] > 0 && !assign(s, i, j, ld_vals[i][j]))
            {
                return 0;
            }
        }
    }

    return 1;
}

static sudoku *create_sudoku(int bdim, char *grid)
{
    assert(bdim <= MAX_BDIM);
    
    sudoku *r = malloc(sizeof(sudoku));
    r->bdim = bdim;
    int dim = bdim * bdim;
    r->dim = dim;
    r->grid_size = dim * dim;
    r->peers_size = 3 * dim - 2 * bdim - 1;
    r->grid = grid;
    r->sol_count = 0; 
    r->unit_list = malloc(sizeof(cell_coord***) * dim);

    assert(r->unit_list);

    for(int i = 0; i < dim; i++)
    {
        r->unit_list[i] = malloc(sizeof(cell_coord**) * dim);

        assert (r->unit_list[i]);

        for(int j = 0; j < dim; j++)
        {
            r->unit_list[i][j] = malloc(sizeof(cell_coord*) * 3);

            assert(r->unit_list[i][j]);

            for(int k = 0; k < 3; k++)
            {
                r->unit_list[i][j][k] = calloc(dim, sizeof(cell_coord));

                assert(r->unit_list[i][j][k]);
            }
        }
    }
    
    r->peers = malloc(sizeof(cell_coord**) * dim);

    assert(r->peers);

    for(int i = 0; i < dim; i++)
    {
        r->peers[i] = malloc(sizeof(cell_coord*) * dim);

        assert(r->peers[i]);

        for(int j = 0; j < dim; j++)
        {
            r->peers[i][j] = calloc(r->peers_size, sizeof(cell_coord));

            assert(r->peers[i][j]);
        }
    }
    
    r->values = malloc (sizeof(cell_v*) * dim);

    assert(r->values);

    for(int i = 0; i < dim; i++)
    {
        r->values[i] = calloc(dim, sizeof(cell_v));

        assert(r->values[i]);
    }
    
    init(r);

    if(grid)    //Se a grid é NULL, não fazemos o parse_grid
    {
        if(!parse_grid(r))
        {
            printf("Error parsing grid\n");
            destroy_sudoku(r);

            return 0;
        }        
    }
    
    return r;
}

static int eliminate (sudoku *s, int i, int j, int d)
{
    int k, ii, cont, pos;
    
    if (!cell_v_get(&s->values[i][j], d))
    {
        return 1;
    }

    cell_v_unset(&s->values[i][j], d);

    int count = cell_v_count(&s->values[i][j]);

    if(count == 0)
    {
        return 0;
    }
    else
    if(count == 1)
    {
        for(k = 0; k < s->peers_size; k++)
        {
            if(!eliminate(s, s->peers[i][j][k].r, s->peers[i][j][k].c, digit_get(&s->values[i][j])))
            {
                return 0;
            }
        }
    }

    for(k = 0; k < 3; k++)
    { 
        cont = 0;
        pos = 0;
        cell_coord* u = s->unit_list[i][j][k];

        for(ii = 0; ii < s->dim; ii++)
        {
            if(cell_v_get(&s->values[u[ii].r][u[ii].c], d))
            {
                cont++;
                pos = ii;
            }
        }

        if(cont == 0)
        {
            return 0;
        }
        else
        if(cont == 1)
        {
            if(!assign(s, u[pos].r, u[pos].c, d))
            {
                return 0;
            }
        }
    }

    return 1;
}

static int assign (sudoku *s, int i, int j, int d)
{
    for(int d2 = 1; d2 <= s->dim; d2++)
    {
        if(d2 != d)
        {
            if(!eliminate(s, i, j, d2))
            {
               return 0;
            }
        }
    }

    return 1;
}

static void display(sudoku *s)
{
    printf("%d\n", s->bdim);

    for(int i = 0; i < s->dim; i++)
    {
        for(int j = 0; j < s->dim; j++)
        {
            printf("%d ",  digit_get(&s->values[i][j]));
        }
        printf("\n");
    }
    printf("\n");
}

static int search (sudoku *s, int status, int *workload)
{
    if(stop_threads)
    {
        return 1;
    }

    int flag;
    MPI_Test(&solved_request, &flag, MPI_STATUS_IGNORE);

    if(flag)
    {
        stop_threads = 1;
        return 1;
    }

    if(!status)
    {
        return 0;
    }

    int solved = 1;

    for(int i = 0; solved && i < s->dim; i++)
    { 
        for(int j = 0; j < s->dim; j++)
        { 
            if(cell_v_count(&s->values[i][j]) != 1)
            {
                solved = 0;

                break;
            }
        }
    }

    if(solved)
    {
        pthread_mutex_lock(&solved_lock);

        if(stop_threads)
        {
            return 1;
        }
        stop_threads = 1;
        char solved_buffer[s->grid_size * sizeof(cell_v)];
        void *buffer_pos = solved_buffer;

        for(int i = 0; i < s->dim; i++, buffer_pos += s->dim * sizeof(cell_v))
        {
            memcpy(buffer_pos, s->values[i], s->dim * sizeof(cell_v));
        }

        MPI_Send(solved_buffer, s->grid_size * sizeof(cell_v), MPI_BYTE, 0, TAG_SOLVED, MPI_COMM_WORLD);
        pthread_mutex_unlock(&solved_lock);

        return 1;
    }

    int min = INT_MAX;
    int minI = -1;
    int minJ = -1;
    int ret = 0;
    
    cell_v **values_bkp = malloc(sizeof(cell_v *) * s->dim);

    for(int i = 0; i < s->dim; i++)
    {
        values_bkp[i] = malloc(sizeof(cell_v) * s->dim);
    }
    
    for(int i = 0; i < s->dim; i++)
    { 
        for(int j = 0; j < s->dim; j++)
        {
            int used = cell_v_count(&s->values[i][j]);

            if (used > 1 && used < min)
            {
                min = used;
                minI = i;
                minJ = j;
            }
        }
    }

    for(int i = 0; i < s->dim; i++)
    {
        for(int j = 0; j < s->dim; j++)
        {
            values_bkp[i][j] = s->values[i][j];
        }
    }

    Node *node = NULL;

    for(int k = 1; k <= s->dim; k++)
    {
        if(cell_v_get(&s->values[minI][minJ], k))
        {            
            if(*workload)
            {
                *workload = *workload - 1;
                if(search(s, assign(s, minI, minJ, k), workload))
                {
                    ret = 1;

                    goto FR_RT;
                }
                else
                {
                    for(int i = 0; i < s->dim; i++)
                    {
                        for (int j = 0; j < s->dim; j++)
                        {
                            s->values[i][j] = values_bkp[i][j];
                        }
                    }
                }                
            }
            else
            {
                node = malloc(sizeof(Node));
                node->values = malloc(sizeof(cell_v *) * s->dim);

                for(int i = 0; i < s->dim; i++)
                {
                    node->values[i] = malloc(sizeof(cell_v) * s->dim);
                }

                for(int i = 0; i < s->dim; i++)
                {
                    for (int j = 0; j < s->dim; j++)
                    {
                        node->values[i][j] = s->values[i][j];
                    }
                }

                node->digit = k;
                node->sqrI = minI;
                node->sqrJ = minJ;
                node->digits_count = min;
                node->dim = s->dim;

                dynarray_add_tail(node_stack, node);

                node = NULL;
            }
        }
    }
    
    FR_RT:
    for (int i = 0; i < s->dim; i++)
    {
        free(values_bkp[i]);
    }

    free (values_bkp);
    
    return ret;
}

void *thread_search(void *_args)
{
    ThreadArgs *args = (ThreadArgs*)_args;
    Node *node;
    int workload = args->workload;
    sudoku *s = args->s;
    args->ret = 0;  //assumimos primeiramente que não chegaremos à solução

    while(workload)
    {
        if(!dynarray_get_count(node_stack))
        {
            workload = 0;   //para poder pedir novos nós
            break;
        }
        
        node = dynarray_remove_tail(node_stack);

        if(!node)
        {
            printf("solve_sudoku: nó nulo encontrado na stack...\n");
            continue;
        }

        s->values = node->values;
        
        if(!assign(s, node->sqrI, node->sqrJ, node->digit))
        {
            for(int i = 0; i < s->dim; i++)
            {
                free(node->values[i]);
            }

            free(node->values);
            free(node);
            node = NULL;
            continue;
        }
        
        free(node);

        if(search(s, 1, &workload))
        {
            args->ret = 1;
            break;
        }
    }

    args->workload = workload;

    return NULL;
}

int solve_sudoku(sudoku *s)
{
    int workload_bkp = min_workload + (rand() % (max_workload - min_workload + 1));
    int workload = 0;
    MPI_Irecv(NULL, 0, MPI_BYTE, 0, TAG_SOLVED, MPI_COMM_WORLD, &solved_request);
    node_stack = dynarray_create(s->grid_size * s->grid_size);
    Node *node;

    if(n_threads < 1)
    {
        n_threads = 1;
    }

    pthread_t threads[n_threads];
    ThreadArgs threadArgs[n_threads];

    for(int i = 0; i < n_threads; i++)
    {
        threadArgs[i].s = malloc(sizeof(sudoku));
        memcpy(threadArgs[i].s, s, sizeof(sudoku));
    }

    if(s)
    {
        int values_size = s->grid_size * sizeof(cell_v);
        int node_size = values_size + 4; //+4: d + i + j + digits_count (4 bytes no total)
        int msg_size = node_size * max_nodes_per_message + 1; //+1: 1 byte para indicar quantos nós foram enviados/recebidos (nodes_amount)
        char msg_buffer[msg_size];
        void *buffer_pos;
        int nodes_amount;

        while(1)
        {
            if(!workload)
            {
                workload = workload_bkp;

                MPI_Send(NULL, 0, MPI_INT, 0, TAG_NEW_JOB, MPI_COMM_WORLD);
                MPI_Recv(msg_buffer, msg_size, MPI_BYTE, 0, TAG_NEW_JOB, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                nodes_amount = msg_buffer[msg_size - 1];

                for(int i = nodes_amount - 1; i >= 0; i--)
                {
                    buffer_pos = msg_buffer + i * node_size;
                    node = malloc(sizeof(Node));
                    node->values = malloc(sizeof(cell_v*) * s->dim);

                    for(int j = 0; j < s->dim; j++)
                    {
                        node->values[j] = malloc(sizeof(cell_v) * s->dim);
                    }

                    for(int j = 0; j < s->dim; j++, buffer_pos += s->dim * sizeof(cell_v))
                    {
                        memcpy(node->values[j], buffer_pos, s->dim * sizeof(cell_v));
                    }

                    node->digit = *(char*)buffer_pos++;
                    node->sqrI = *(char*)buffer_pos++;
                    node->sqrJ = *(char*)buffer_pos++;
                    node->digits_count = *(char*)buffer_pos++;

                    dynarray_add_tail(node_stack, node);

                    node = NULL;
                }
            }

            while(1)
            {

                for(int i = 0; i < n_threads; i++)
                {
                    threadArgs[i].workload = workload;
                    pthread_create(&threads[i], NULL, thread_search, &threadArgs[i]);
                }

                int min_wl_remaining = INT_MAX;

                for(int i = 0; i < n_threads; i++)
                {
                    pthread_join(threads[i], NULL);

                    if(threadArgs[i].ret)
                    {
                        return 1;
                    }

                    if(threadArgs[i].workload < min_wl_remaining)
                    {
                        min_wl_remaining = threadArgs[i].workload;
                    }
                }

                workload = min_wl_remaining;
                
                if(!workload)
                {
                    while(dynarray_get_count(node_stack))
                    {
                        nodes_amount = 0;
                        buffer_pos = msg_buffer;

                        for(int i = 0; i < max_nodes_per_message; i++)
                        {
                            node = dynarray_remove_tail(node_stack);

                            if(!node)
                            {
                                //printf("sudoku_solve: nó nulo encontrado na stack...\n");
                                break;
                            }

                            for(int j = 0; j < s->dim; j++, buffer_pos += s->dim * sizeof(cell_v))
                            {
                                memcpy(buffer_pos, node->values[j], s->dim * sizeof(cell_v));
                            }

                            *(char*)buffer_pos++ = node->digit;
                            *(char*)buffer_pos++ = node->sqrI;
                            *(char*)buffer_pos++ = node->sqrJ;
                            *(char*)buffer_pos++ = node->digits_count;                    
                            nodes_amount++;

                            for(int j = 0; j < s->dim; j++)
                            {
                                free(node->values[j]);
                            }

                            free(node->values);
                            free(node);

                            node = NULL;
                        }

                        msg_buffer[msg_size - 1] = nodes_amount;
                        
                        //printf("Processo %d\t envia %d\t nós ao processo 0\t\n", mpi_rank, nodes_amount);

                        MPI_Send(msg_buffer, msg_size, MPI_BYTE, 0, TAG_ADD_NODES, MPI_COMM_WORLD);
                    }
                    break;
                }

            }
        }        
    }
    else
    {
        printf("Puzzle == NULL\n");
        return 0;
    }

    return 1;
}

unsigned int sum_values(cell_v **values, int dim)
{
    unsigned int sum = 0;

    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            sum += cell_v_count(&values[i][j]);
        }
    }

    return sum;
}

unsigned int get_node_value(Node *node)
{
    unsigned int max_possible_values_sum = node->dim * node->dim * node->dim; //Existem dim*dim quadrados, cada um com no máximo dim possibilidades
    unsigned int current_possible_values_sum = sum_values(node->values, node->dim);
    double possib_digits_sqr_ratio = node->digits_count / (double) (node->dim);
    double possib_digits_board_ratio = 1.0 - current_possible_values_sum / (double) max_possible_values_sum;

    return (unsigned int)((0.8 * possib_digits_board_ratio + 0.2 * possib_digits_sqr_ratio) * UINT_MAX);
}

int coordinate(sudoku *s)
{
    minheap *nodes = minheap_create();

    Node *node;

    for(int i = 0; i < s->dim; i++)
    {
        for(int j = 0; j < s->dim; j++)
        {
            int digits_count = cell_v_count(&s->values[i][j]);

            if(digits_count <= 1)
            {
                continue;
            }

            for(char d = 1; d <= s->dim; d++)
            {
                if(cell_v_get(&s->values[i][j], d)) 
                {
                    node = malloc(sizeof(Node));
                    node->digit = d;
                    node->sqrI = i;
                    node->sqrJ = j;
                    node->digits_count = digits_count;
                    node->dim = s->dim;                    
                    node->values = malloc(sizeof(cell_v*) * s->dim);

                    for(int k = 0; k < s->dim; k++)
                    {
                        node->values[k] = malloc(sizeof(cell_v) * s->dim);
                    }

                    for(int k = 0; k < s->dim; k++)
                    {
                        for(int l = 0; l < s->dim; l++)
                        {
                            node->values[k][l] = s->values[k][l];
                        }
                    }

                    minheap_add(nodes, node, get_node_value(node));
                    node = NULL;
                }
            }
        }
    }

    MPI_Status status;
    int values_size = s->grid_size * sizeof(cell_v);
    int node_size = values_size + 4; //+4: d + i + j + digits_count (4 bytes no total)
    int msg_size = node_size * max_nodes_per_message + 1; //+1: 1 byte para indicar quantos nós foram enviados
    char msg_buffer[msg_size];
    void *buffer_pos;
    char solved_buffer[values_size];
    int nodes_amount;
    int already_solved = 0;
    int processes_terminated = 0;

    do
    {
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        switch(status.MPI_TAG)
        {
            case TAG_NEW_JOB:
                MPI_Recv(NULL, 0, MPI_INT, status.MPI_SOURCE, TAG_NEW_JOB, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                nodes_amount = 0;
                buffer_pos = msg_buffer;

                for(int i = 0; i < coord_send_nodes_amount; i++)
                {
                    node = minheap_remove_min(nodes);

                    if(!node)
                    {
                        printf("coordinate: nó nulo encontrado na heap...\n");
                        continue;
                    }

                    for(int j = 0; j < s->dim; j++, buffer_pos += s->dim * sizeof(cell_v))
                    {
                        memcpy(buffer_pos, node->values[j], s->dim * sizeof(cell_v));
                    }

                    *(char*)buffer_pos++ = node->digit;
                    *(char*)buffer_pos++ = node->sqrI;
                    *(char*)buffer_pos++ = node->sqrJ;
                    *(char*)buffer_pos++ = node->digits_count;                    
                    nodes_amount++;

                    for(int j = 0; j < s->dim; j++)
                    {
                        free(node->values[j]);
                    }

                    free(node->values);
                    free(node);

                    node = NULL;
                }

                msg_buffer[msg_size - 1] = nodes_amount;

                MPI_Send(msg_buffer, msg_size, MPI_BYTE, status.MPI_SOURCE, TAG_NEW_JOB, MPI_COMM_WORLD);

                break;

            
            case TAG_ADD_NODES:
                MPI_Recv(msg_buffer, msg_size, MPI_BYTE, status.MPI_SOURCE, TAG_ADD_NODES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                nodes_amount = msg_buffer[msg_size - 1];
                buffer_pos = msg_buffer;

                for(int i = 0; i < nodes_amount; i++)
                {
                    node = malloc(sizeof(Node));
                    node->values = malloc(sizeof(cell_v*) * s->dim);

                    for(int j = 0; j < s->dim; j++)
                    {
                        node->values[j] = malloc(sizeof(cell_v) * s->dim);
                    }

                    for(int j = 0; j < s->dim; j++, buffer_pos += s->dim * sizeof(cell_v))
                    {
                        memcpy(node->values[j], buffer_pos, s->dim * sizeof(cell_v));
                    }

                    node->digit = *(char*)buffer_pos++;
                    node->sqrI = *(char*)buffer_pos++;
                    node->sqrJ = *(char*)buffer_pos++;
                    node->digits_count = *(char*)buffer_pos++;
                    node->dim = s->dim;

                    minheap_add(nodes, node, get_node_value(node));

                    node = NULL;
                }
                break;
            
            case TAG_SOLVED:
                if(!already_solved)
                {
                    already_solved = 1;
                }    
                else
                {
                    break;
                }

                gettimeofday(&end_time, NULL);
                double elapsed_time = (((end_time.tv_sec * 1000000 + end_time.tv_usec) - (start_time.tv_sec * 1000000 + start_time.tv_usec)))/1000000.0;

                MPI_Recv(&solved_buffer, values_size, MPI_BYTE, status.MPI_SOURCE, TAG_SOLVED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                printf("Puzzle resolvido pelo processo %d:\n", status.MPI_SOURCE);

                buffer_pos = solved_buffer;

                FILE *results = fopen("results.txt", "a");
                fprintf(results, "%d,%d,%d,%d,%d,%d,%lf\n", max_workload, min_workload, max_nodes_per_message, coord_send_nodes_amount, mpi_size, n_threads, elapsed_time);
                fclose(results);

                MPI_Abort(MPI_COMM_WORLD, MPI_SUCCESS);

                for(int i = 0; i < s->dim; i++, buffer_pos += s->dim * sizeof(cell_v))
                {
                    memcpy(s->values[i], buffer_pos, s->dim * sizeof(cell_v));    
                }

                //display(s);

                for(int i = 1; i < mpi_size; i++)
                {
                    MPI_Send(NULL, 0, MPI_BYTE, i, TAG_SOLVED, MPI_COMM_WORLD);                        
                }

                break;

            case TAG_TERMINATED:
                MPI_Recv(NULL, 0, MPI_BYTE, status.MPI_SOURCE, TAG_TERMINATED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                processes_terminated++;
                printf("%d %d\n", processes_terminated, mpi_size - 1);
                if(processes_terminated == mpi_size - 1)
                {
                    return 1;
                }

                break;

            default:
                printf("switch case erro inesperado (MPI_TAG = %d)\n", status.MPI_TAG);
                break;
        }
    }
    while(minheap_get_count(nodes));

    printf("Puzzle não resolvido\n");
    MPI_Abort(MPI_COMM_WORLD, MPI_SUCCESS);

    return 0;
}

str2int_errno str2int(int *out, char *s, int base) {
    char *end;
    if (s[0] == '\0' || isspace(s[0]))
        return STR2INT_INCONVERTIBLE;
    errno = 0;
    long l = strtol(s, &end, base);

    if (l > INT_MAX || (errno == ERANGE && l == LONG_MAX))
        return STR2INT_OVERFLOW;
    if (l < INT_MIN || (errno == ERANGE && l == LONG_MIN))
        return STR2INT_UNDERFLOW;
    if (*end != '\0')
        return STR2INT_INCONVERTIBLE;
    *out = l;
    return STR2INT_SUCCESS;
}

int main (int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    assert(argc == 6);
    assert(str2int(&max_workload, argv[1], 10) == STR2INT_SUCCESS);
    assert(str2int(&min_workload, argv[2], 10) == STR2INT_SUCCESS);
    assert(str2int(&max_nodes_per_message, argv[3], 10) == STR2INT_SUCCESS);
    assert(str2int(&coord_send_nodes_amount, argv[4], 10) == STR2INT_SUCCESS);
    assert(str2int(&n_threads, argv[5], 10) == STR2INT_SUCCESS);

    if(mpi_size < 2)
    {
        printf("É preciso de no mínimo 2 processos\n");
        return 0;
    }
    gettimeofday(&start_time, NULL);
    pthread_mutex_init(&solved_lock, NULL);

    if(!mpi_rank)
    {
        int size;
        assert(scanf("%d", &size) == 1);
        assert(size <= MAX_BDIM);
        int board_size = size * size * size * size;
        char board[board_size];
    
        for (int i = 0; i < board_size; i++)
        {
            int temp;
            if (scanf("%d", &temp) != 1)
            {
                printf("Erro ao carregar arquivo (%d)\n", i);
                exit(1);
            }
            else
            {
                board[i] = temp;
            }
        }

        for(int i = 1; i < mpi_size; i++)
        {
            MPI_Send(&size, 1, MPI_INT, i, TAG_SIZE, MPI_COMM_WORLD);
        }

        sudoku *s = create_sudoku(size, board);
        coordinate(s);
    }
    else
    {
        srand(time(NULL));
        
        for(int i = 0; i < mpi_rank - 1; i++)
        {
            rand();
        }

        int size;

        MPI_Recv(&size, 1, MPI_INT, 0, TAG_SIZE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        sudoku *s = create_sudoku(size, NULL);
        
        solve_sudoku(s);
    }

    if(mpi_rank)
    {
        MPI_Send(NULL, 0, MPI_BYTE, 0, TAG_TERMINATED, MPI_COMM_WORLD);
    }

    pthread_mutex_destroy(&solved_lock);
    MPI_Finalize();

    return 0;
}
