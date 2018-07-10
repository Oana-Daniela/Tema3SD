#ifndef __SOLVE_H__
#define __SOLVE_H__
#include <stdlib.h>
#include <limits.h>

#define LINE 0
#define COLUMN 1
#define NO_ZERO -1

typedef struct _TGraphM {
  int nn; // numar noduri
  int **Ma; // matrice de adiacenta
} TGraphM;

TGraphM* copy(TGraphM graph)
{
  int i, j;

  TGraphM* my_copy = malloc(sizeof(TGraphM));
  my_copy->nn = graph.nn;
  my_copy->Ma = malloc(my_copy->nn * sizeof(int *));
  for (i = 0; i < graph.nn; i++)
  {
    my_copy->Ma[i] = malloc(my_copy->nn * sizeof(int));
    for (j = 0; j < graph.nn; j++)
    {
      my_copy->Ma[i][j] = graph.Ma[i][j];
    }
  }

  return my_copy;
}

TGraphM* read_file(char * testInputFileName)
{
  int i, j;

  TGraphM* graph = malloc(sizeof(TGraphM));
  if (graph == NULL)
  { // verific ca s-a putut aloca memorie
    return NULL;
  }

  FILE* file = fopen(testInputFileName, "r");
  if (file == NULL)
  { // verific faptul ca s-a putut deschide fisierul spre citire
    free(graph); // daca nu s-a putut fac rollback la toate operatiile realizate si intorc NULL
    return NULL;
  }
  // daca totul e in ordine incep sa citesc fisierul
  fscanf(file, "%d", &graph->nn);
  graph->Ma = malloc(graph->nn * sizeof(int *));
  for (i = 0; i < graph->nn; i++)
  {
    graph->Ma[i] = malloc(graph->nn * sizeof(int));
    for (j = 0; j < graph->nn; j++)
    {
      fscanf(file, "%d", &graph->Ma[i][j]);
    }
  }
  fclose(file);

  return graph;
}

/*
 * Functia returneaza minorantul unei linii sau a unei coloane din matricea Ma.
 * Daca isLine este LINE atunci se selecteaza linia, daca e COLUMN atunci se va selecta coloana.
 */
int lower_bound(TGraphM graph, int index, int isLine)
{
  int i, min = INT_MAX; // INT_MAX este o constanta din limits.h

  if (isLine == LINE)
  {
    for (i = 0; i < graph.nn; i++)
    {
      if (graph.Ma[index][i] < min)
      {
        min = graph.Ma[index][i];
      }
    }
  } 
  else
  {
    for (i = 0; i < graph.nn; i++)
    {
      if (graph.Ma[i][index] < min)
      {
        min = graph.Ma[i][index];
      }
    }
  }

  return min;
}

void substract_lower_bound(TGraphM *graph, int lower_bound, int index, int isLine)
{
  int i;

  if (isLine == LINE)
  {
    for (i = 0; i < graph->nn; i++)
    {
      graph->Ma[index][i] -= lower_bound;
    }
  } 
  else
  {
    for (i = 0; i < graph->nn; i++)
    {
      graph->Ma[i][index] -= lower_bound;
    }
  }
}

int first_elem(TGraphM graph, int index, int value, int isLine)
{
  int i;

  if (isLine == LINE)
  {
    for (i = 0; i < graph.nn; i++)
    {
      if (graph.Ma[index][i] == value)
      {
        return i;
      }
    }
  } 
  else
  {
    for (i = 0; i < graph.nn; i++)
    {
      if (graph.Ma[i][index] == value)
      {
        return i;
      }
    }
  }

  return NO_ZERO;
}

int min_zero_line(TGraphM graph)
{
  int i, j, line_count, min_value = INT_MAX, min_line = -1;

  for (i = 0; i < graph.nn; i++)
  {
    line_count = 0;
    for (j = 0; j < graph.nn; j++)
    {
      if (graph.Ma[i][j] == 0)
      {
        line_count++;
      }
      else 
        if (graph.Ma[i][j] == 1)
        {
          // daca gasim un zero incadrat atunci line_countul va fi maxim
          line_count = INT_MAX;
          break; // nu mai are rost sa mergem mai departe
        }
    }
    if (line_count < min_value && line_count > 0)
    {
      min_value = line_count;
      min_line = i;
    }
  }

  return min_line;
}

int is_finished(TGraphM graph)
{
  int i, j, line_counter;

  for (i = 0; i < graph.nn; i++)
  {
    line_counter = 0;
    for (j = 0; j < graph.nn; j++)
    {
      if (graph.Ma[i][j] == 0)
      {
        return 0; // o prima conditie de terminare e sa nu mai avem zero in matricea support
      }
      if (graph.Ma[i][j] == 1)
      {
        line_counter++;
      }
    }
    if (line_counter != 1)
    {
      return 0; // o a doua coditie sa avem un zero marcat pe fiecare linie (nici mai mult nici mai putin)
    }
  }

  // daca se ajunge aici inseamna ca totul e in ordine
  return 1;
}

void do_step_2(TGraphM *graph)
{
  int i;

  // mai intai obtinem cate un zero pe fiecare linie din matrice
  for (i = 0; i < graph->nn; i++)
  {
    substract_lower_bound(graph, lower_bound(*graph, i, LINE), i, LINE);
  }

  // verificam ca avem cel putin un zero si pe fiecare coloana
  for (i = 0; i < graph->nn; i++)
  {
    int zero_pos = first_elem(*graph, i, 0, COLUMN);
    if (zero_pos == NO_ZERO)
    {
      substract_lower_bound(graph, lower_bound(*graph, i, COLUMN), i, COLUMN);
    }
  }

  // in acest moment avem cate un 0 pe fiecare linie si coloana
}

void do_step_4(TGraphM *support, TGraphM *graph, int *marked_lines, int *marked_columns)
{
  int i, j, min = INT_MAX;

  // se marcheaza liniile care nu contin niciun zero incadrat
  // cu alte cuvinte, se marcheaza liniile care nu contin in matricea support niciun 1
  for (i = 0; i < support->nn; i++)
  {
    if (first_elem(*support, i, 1, LINE) == -1)
    {
      marked_lines[i] = 1;
    }
  }

  // cat timp mai avem ceva de marcat
  int marked_flag = 1;
  while (marked_flag)
  {
    marked_flag = 0;
    // se marcheaza coloanele care contin zerouri taiate pe liniile marcate
    for (i = 0; i < support->nn; i++)
    {
      if (marked_lines[i] == 1)
      {
        // obtine toate coloanele care au zero-uri taiate si marcheazale
        for (j = 0; j < support->nn; j++)
        {
          if (support->Ma[i][j] == 2 && marked_columns[j] == 0)
          {
            // se marcheaza coloana j
            marked_columns[j] = 1;
            marked_flag = 1; // algoritmul va continua
          }
        }
      }
    }
    if (marked_flag == 1)
    {
      // se marcheaza liniile care contin un zero incadrat pe coloanele deja marcate
      for (i = 0; i < support->nn; i++)
      {
        if (marked_lines[i] == 0)
        {
          // pentru liniile nemarcate, verificam conditia
          for (j = 0; j < support->nn; j++)
          {
            if (support->Ma[i][j] == 1 && marked_columns[j] == 1)
            {
              marked_lines[i] = 1;
              marked_flag = 1; // algoritmul va continua si in acest caz
            }
          }
        }
      }
    }
  }

  // se determina minorantul din liniile si coloanele care au ramas netaiate
  for (i = 0; i < support->nn; i++)
  {
    for (j = 0; j < support->nn; j++)
    {
      if (marked_lines[i] == 1 && marked_columns[j] == 0 && graph->Ma[i][j] < min)
      {
        min = graph->Ma[i][j];
      }
    }
  }

  // se fac schimbarile legate de minorant in matricea graph
  for (i = 0; i < support->nn; i++)
  {
    for (j = 0; j < support->nn; j++)
    {
      // avem 3 cazuri:
      // cazul 1: minorantul se adauga elementelor dublu taiate
      // cum determinam un element dublu taiat? foarte simplu:
      // > linia nu trebuie sa fie marcata, iar coloana trebuie sa fie marcata
      if (marked_lines[i] == 0 && marked_columns[j] == 1)
      {
        graph->Ma[i][j] += min;
      }
      // cazul 2: minorantul se scade din elementele netaiate
      if (marked_lines[i] == 1 && marked_columns[j] == 0)
      {
        graph->Ma[i][j] -= min;
      }
      // cazul 3: elementele simplu taiate sunt lasate neschimbate
    }
  }

  // in acest moment, pe noua matrice se poate repeta procedeul de la pasul 3
}

TGraphM* do_step_3(TGraphM *graph) 
{
  int i, j;

  int* marked_lines = malloc(graph->nn * sizeof(int));
  int* marked_columns = malloc(graph->nn * sizeof(int));

  while (1)
  {
    // avem nevoie de inca o copie a matricei de adiacenta pentru stocarea elementelor incadrate si barate
    TGraphM* support = copy(*graph);
    // toate elementele diferite de 0 din matricea support sunt marcate cu NO_ZERO
    for (i = 0; i < support->nn; i++)
    {
      for (j = 0; j < support->nn; j++)
      {
        if (support->Ma[i][j] != 0)
        {
          support->Ma[i][j] = NO_ZERO;
        }
      }
    }
    // se cauta un zero pe fiecare linie
    int min_index = -1, zero_pos;
    do {
      min_index = min_zero_line(*support);
      if (min_index > -1)
      {
        zero_pos = first_elem(*support, min_index, 0, LINE);
        // se marcheaza pozitia si se bareaza toate ceelalte zero-uri de pe linie sau coloana
        support->Ma[min_index][zero_pos] = 1;
        // bararea inseamna ca in matricea support vom pune 2
        for (i = 0; i < support->nn; i++)
        {
          // in acelasi for vom trata si linia si coloana
          if (support->Ma[min_index][i] == 0)
          {
            support->Ma[min_index][i] = 2;
          }
          if (support->Ma[i][zero_pos] == 0)
          {
            support->Ma[i][zero_pos] = 2;
          }
        }
      }
    } while (min_index != -1);

    if (is_finished(*support))
    {
      // am gasit solutia, eliberez memoria alocata dinamic pentru liniile si coloanele marcate
      free(marked_lines);
      free(marked_columns);
      // returnez solutia pe care am gasit-o
      return support;
    } 
    else
    {
      // toate liniile si coloanele devin nemarcate
      for (i = 0; i < graph->nn; i++)
      {
        marked_lines[i] = 0;
        marked_columns[i] = 0;
      }
      do_step_4(support, graph, marked_lines, marked_columns);
    }
  }
}

int solve(char *testInputFileName)
{
  int i, j, total_cost = 0;

  // algoritmul se gaseste la adresa: gabi.bizoi.home.ro/Algoritmul%20ungar_aplicatie.pdf

  // citeste fisierul de input
  TGraphM* graph = read_file(testInputFileName);
  if (graph == NULL)
  {
    printf("Fisierul nu a putut fi deschis sau nu mai exista memorie libera pe heap!\n");
    exit(1); // programul intoarce eroare
  }

  // pasul 1 este deja implementat - matricea fiind patratica
  // pana sa trecem la pasul 2, salvam matricea initiala de costuri
  TGraphM* copy_graph = copy(*graph);

  // trecem la pasul 2, aici va trebui sa formam cate un 0 pe fiecare linie si coloana
  do_step_2(graph);

  TGraphM* result = do_step_3(graph);
  if (result != NULL)
  {
    // daca avem o solutie optima, calculam costul final ca fiind suma costurilor de pe pozitiile marcate
    for (i = 0; i < graph->nn; i++)
    {
      for (j = 0; j < graph->nn; j++)
      {
        if (result->Ma[i][j] == 1)
        {
          total_cost += copy_graph->Ma[i][j]; 
        }
      }
    }
  }
  return total_cost;
}

#endif
