#ifndef _FIT_FORMULAS_H_
#define _FIT_FORMULAS_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mymatrices.h"

typedef struct {
	double **data; /*enthaelt bis zu 4 double arrays*/
	int len;
	int x, y, e_x, e_y; /* gibt die indizes der arrays in data an, ungenutztes auf -1 setzen */
} dataset;

typedef struct {
	int x;
	int y;
	int e_x;
	int e_y;
	int start;
	int len;
	double min_chng;
	int max_iter;
} fit_conf;

/*Mittelwert: Array, länge, startpunkt, 0=x 1=y*/
double average(dataset *dat, int start, int xy);

/*Varianz: Array, länge, startpunkt, mittelwert, index des datensatzes*/
double variance(dataset *dat, fit_conf *conf, double avg, int xy);

/*Initial Fit ohne Fehlerberuecksichtigung*/
double *init_fit(dataset *dat, fit_conf *conf);

/*Varianzgewichteter Fit (iterativ): Array, Fehlerarray, laenge, startpunkt, ausgangswert fuer steigung  */
/*Rueckgabe 0 - inc, 1 - inc_err, 2 - abs, 3 - abs_err, 4-8 - reserviert*/
/*Rueckgabe muss freigegeben werden mit free() !!! */
double *m_fit(dataset *dat, fit_conf *conf, double inclin, double absol);

/*Varianz fuer gewichtung: Fehlerarray, laenge, startpunkt, steigung*/
double *m_variance(dataset *dat, fit_conf *conf, double inclin);

/* Korellationskoeffizient -> sagt aus wie gut ein linearer zusammenhang besteht.
 *  Betrag gegen 1 -> nah an linearem Zusammenhang*/
double correlation(dataset *dat, fit_conf *conf);

void array_print(double *a, int n);

/*liest die Daten aus einer Datei aus und gibt Datenarray zurueck*/
double **read_data(char *);

/*Fittet die Daten und gibt das Fit-Array zurueck,
 *0 - Steigung, 1 - Fehler der Steigung, 2 - Absolutes Glied, 3 - Fehler Abs. Gl., 4 - Korellation,
 *5 - Art des Abbruchs, 6 Anzahl Iterationen, 7 letzte Differenz, 8 Guete (Chi^2 / d.o.f)
 *Arten des Abbruchs: 1 - max. Iterationen, 2 - aenderung zu klein, 3 - eingependelt, 4 - speicherallokation*/
double *linear_fit(dataset *dat, fit_conf *conf);

double *linear_file_fit(char *file, fit_conf *conf);

/* laedt daten aus datei aus, falls ein error nicht genutzt wird e_x bzw e_y auf -1 setzen*/
dataset *load_data(char *file, fit_conf *conf);

dataset *dataset_create(int sets, int len);

void dataset_free(dataset *dat);

fit_conf *load_config();

int write_config();

int fit_log_push(double *dat);

/* Gibt die Infos ueber den Fit aus */
int fit_print(double *dat);

#endif
