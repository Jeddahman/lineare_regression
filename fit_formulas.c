#include "fit_formulas.h"

/*Mittelwert: Array, länge, startpunkt, 0=x 1=y*/
double average(dataset *dat, int start, int xy) {
	int k;
	double res = 0.;
	if (xy == 0)
		xy = dat->x;
	else
		xy = dat->y;

	if (start >= dat->len - 1)
		return 0;
	for (k = start; k < dat->len + start; ++k) {
		res += dat->data[xy][k];
	}
	res /= (k - start);
	return res;
}

/*Varianz: Array, länge, startpunkt, mittelwert, index des datensatzes*/
double variance(dataset *dat, fit_conf *conf, double avg, int xy) {
	int k;
	double res = 0.;
	if (xy == 0)
		xy = dat->e_x;
	else
		xy = dat->e_y;

	if (conf->start >= conf->len - 1)
		return 0;
	for (k = conf->start; k < conf->len + conf->start; ++k) {
		res += pow((dat->data[xy][k] - avg), 2);
	}
	res /= (k - conf->start);
	return res;
}

double *init_fit(dataset *dat, fit_conf *conf) {
	double *res, v_sq_x, v_x, v_xy, v_y, var_y;
	int k;

	res = malloc(9 * sizeof(double));
	if (res == NULL) {
		return NULL;
	}

	/*ausrechnen der Summen */
	v_sq_x = 0; /*x^2 / v_sq*/
	v_x = 0; /*x / v_sq*/
	v_xy = 0; /*xy / v_sq*/
	v_y = 0; /*y / v_sq*/

	for (k = conf->start; k < conf->start + conf->len; ++k) {
		v_sq_x += pow(dat->data[dat->x][k], 2);
		v_x += dat->data[dat->x][k];
		v_xy += (dat->data[dat->x][k] * dat->data[dat->y][k]);
		v_y += dat->data[dat->y][k];
	}

	var_y = variance(dat, conf, v_x, conf->y);

	/*Steigungsberechnung */
	res[0] = v_xy - v_x * v_y;
	res[0] /= v_sq_x - pow(v_x, 2);

	/*Fehler der Steigung */
	res[1] = var_y/(conf->len);
	res[1] /= (v_sq_x - pow(v_x,2));

	/*Absolutes Glied*/
	res[2] = v_y - (res[0] * v_x);

	/*Fehler Absolutes Glied*/
	res[3] = res[1] * v_sq_x;
	/* hier muss evtl. noch gegausst werden! */

	/* Korellationskoeffizient -> wird in uebergeorndeter Funktion berechnet*/
	res[4] = 0;
	/* reserviert fuer Abbruchinformationen */
	res[5] = 0;
	/* reserviert fuer Anzahl Iterationen */
	res[6] = 0;
	/* reserviert fuer letzte Aenderung*/
	res[7] = 0;
	/* chi_sq Guete des Fits */
	res[8] = 0;

	return res;
}

/*Varianzgewichteter Fit (iterativ): Array, Fehlerarray, laenge, startpunkt, ausgangswert fuer steigung  */
/*Rueckgabe 0 - inc, 1 - inc_err, 2 - abs, 3 - abs_err, 4-8 reserviert*/
double *m_fit(dataset *dat, fit_conf *conf, double inclin, double absol) {
	double *vari, *res, vsum, v_sq_x, v_x, v_xy, v_sq_y, v_y, v_sq, rez_d,
			chi_sq;
	int k;

	/*Ansetzen der Varianzgewichtung fuer gegebenes inclin */

	vari = m_variance(dat, conf, inclin);
	if (vari == NULL)
		return NULL;
	res = malloc(9 * sizeof(double));
	if (res == NULL) {
		free(vari);
		return NULL;
	}

	/*ausrechnen der Summen */
	v_sq = 0; /* 1/sigma_i^2*/
	vsum = 0; /* Summe 1/sigma_i^2*/
	v_sq_x = 0; /*x^2 / v_sq*/
	v_x = 0; /*x / v_sq*/
	v_xy = 0; /*xy / v_sq*/
	v_sq_y = 0; /*y^2 / v_sq*/
	v_y = 0; /*y / v_sq*/
	chi_sq = 0; /* Guete des Fits */

	for (k = conf->start; k < conf->start + conf->len; ++k) {
		v_sq = vari[k];
		vsum += 1. / v_sq;
		v_sq_x += pow(dat->data[dat->x][k], 2) / v_sq;
		v_x += dat->data[dat->x][k] / v_sq;
		v_xy += (dat->data[dat->x][k] * dat->data[dat->y][k]) / v_sq;
		v_sq_y += pow(dat->data[dat->y][k], 2) / v_sq;
		v_y += dat->data[dat->y][k] / v_sq;
		chi_sq += pow(
				dat->data[dat->y][k] - inclin * dat->data[dat->x][k] - absol, 2)
				/ v_sq;
	}
	chi_sq /= conf->len - 2;
	/*Steigungsberechnung */
	rez_d = vsum * v_sq_x - pow(v_x, 2);
	rez_d = 1. / rez_d;
	res[0] = vsum * v_xy - v_x * v_y;
	res[0] *= rez_d;

	/*Fehler der Steigung */
	res[1] = rez_d * vsum;
	res[1] = sqrt(res[1]);

	/*Absolutes Glied*/
	res[2] = v_y - (res[0] * v_x);
	res[2] /= vsum;/*Varianzgemittelte Mittelwerte*/

	/*Fehler Absolutes Glied*/
	res[3] = res[1] * v_x / vsum;
	res[3] = sqrt(res[3]);
	/* hier muss evtl. noch gegausst werden! */

	/* Korellationskoeffizient -> wird in uebergeorndeter Funktion berechnet*/
	res[4] = 0;
	/* reserviert fuer Abbruchinformationen */
	res[5] = 0;
	/* reserviert fuer Anzahl Iterationen */
	res[6] = 0;
	/* reserviert fuer letzte Aenderung*/
	res[7] = 0;
	/* chi_sq Guete des Fits */
	res[8] = chi_sq;

	free(vari);
	return res;
}

/*Varianz fuer gewichtung: Fehlerarray, laenge, startpunkt, steigung*/
double *m_variance(dataset *dat, fit_conf *conf, double inclin) {
	double *res, hx, hy;
	int k;

	res = malloc(dat->len * sizeof(double));
	if (res == NULL)
		return NULL;

	for (k = conf->start; k < conf->start + conf->len; ++k) {
		if (dat->e_x != -1) {
			hx = inclin * dat->data[dat->e_x][k];
			hx = pow(hx, 2);
		} else
			hx = 0;
		if (dat->e_y != -1)
			hy = pow(dat->data[dat->e_y][k], 2);
		else
			hy = 0;
		res[k] = hx + hy;
	}
	return res;
}

double correlation(dataset *dat, fit_conf *conf) {
	double res, a_xy, a_x, a_y, a_sq_x, a_sq_y;
	int k;

	a_x = 0; /* arithmetisches mittel von x*/
	a_y = 0; /* arithmetisches mittel von y*/
	a_xy = 0; /* arithmetisches mittel von xy*/
	a_sq_x = 0; /* arithmetisches mittel von x^2*/
	a_sq_y = 0; /* arithmetisches mittel von y^2*/

	for (k = conf->start; k < conf->start + conf->len; ++k) {
		a_x += dat->data[dat->x][k];
		a_y += dat->data[dat->y][k];
		a_xy += dat->data[dat->x][k] * dat->data[dat->y][k];
		a_sq_x += pow(dat->data[dat->x][k], 2);
		a_sq_y += pow(dat->data[dat->y][k], 2);
	}
	a_x /= dat->len;
	a_y /= dat->len;
	a_xy /= dat->len;
	a_sq_x /= dat->len;
	a_sq_y /= dat->len;

	res = a_xy - a_x * a_y;
	res /= sqrt(a_sq_x - pow(a_x, 2)) * sqrt(a_sq_y - pow(a_y, 2));
	return res;
}

void array_print(double *a, int n) {
	int k;
	for (k = 0; k < n - 1; ++k) {
		if (k % 100 == 0 && k != 0)
			printf("\n\\");
		printf("%f, ", a[k]);
	}
	printf("%f\n", a[k]);
}

double *linear_fit(dataset *dat, fit_conf *conf) {
	double *fit, *fit_act;
	double chng = 1;
	double chng_o = 0;
	int k, l;

	fit = malloc(9 * sizeof(double));
	if (!fit)
		return NULL;

	if (!(fit_act = init_fit(dat, conf)))
		return NULL;
	for (k = 0;
			(k < conf->max_iter) && (chng >= conf->min_chng)
					&& (fabs(chng_o - chng) > 1e-20); ++k) {
		for (l = 0; l < 5; ++l)
			fit[l] = fit_act[l];
		free(fit_act);
		fit_act = m_fit(dat, conf, fit[0], fit[2]);
		if (fit_act != NULL) {
			chng_o = chng;
			chng = fabs(fit[0] - fit_act[0]);
		} else {
			fit_act = fit;
			fit_act[5] = 4;
			/*printf("Speicher konnte nicht zugewiesen werden,\n");*/
			break;
		}
	}
	if (k == conf->max_iter) {
		/*printf("Maximale Anzahl Iterationen erreicht.\n");*/
		fit_act[5] = 1;
	} else if (chng < conf->min_chng) {
		/*printf("Aenderung wurde zu klein.\n");*/
		fit_act[5] = 2;

	} else if (fabs(chng_o - chng) < 1e-20) {
		/*printf("Daten eingependelt.\n");*/
		fit_act[5] = 3;
	}
	/* 4 = speicher konnte nicht zugewiesen werden, iteration musste abgebrochen werden*/

	/*Rueckgabe der Anzahl an Iterationen*/
	fit_act[6] = k;
	/*Rueckgabe der letzten Aenderung*/
	fit_act[7] = chng;
	/*Berechnung des Korellationskoeffizienten */
	fit_act[4] = correlation(dat, conf);

	if (fit != fit_act)
		free(fit);
	return fit_act;
}

double *linear_file_fit(char *file, fit_conf *conf) {
	dataset *dts;
	double *res;
	dts = load_data(file, conf);
	if (!dts)
		return NULL;
	res = linear_fit(dts, conf);
	return res;
}

dataset *load_data(char *file, fit_conf *conf) {
	dataset *res;
	FILE *fp = fopen(file, "r");
	if (!fp)
		return NULL;
	char ch;
	int sets = 1, len = 1, k, l;

	while (!feof(fp)) {
		fscanf(fp, "%c", &ch);
		if (ch == '\n')
			++len;
		else if ((ch == ' ' || ch == '\t') && len == 1)
			++sets;
	}
	fseek(fp, 0, SEEK_SET);

	res = dataset_create(sets, len);

	for (k = 0; k < len; ++k) {
		for (l = 0; l < sets; ++l) {
			fscanf(fp, "%lf", &(res->data[l][k]));
		}
	}
	fclose(fp);
	if (conf->len == 0)
		conf->len = len;
	res->len = len;
	res->x = conf->x;
	res->y = conf->y;
	res->e_x = conf->e_x;
	res->e_y = conf->e_y;
	return res;
}

dataset *dataset_create(int sets, int len) {
	dataset *res;
	int k;

	res = malloc(sizeof(dataset));
	if (!res)
		return res;
	res->data = malloc(sets * sizeof(void *));
	if (!res->data) {
		free(res);
		return NULL;
	}
	for (k = 0; k < sets; ++k) {
		res->data[k] = malloc(len * sizeof(double));
		if (!res->data[k])
			break;
	}
	if (k < sets) {
		for (; k > 0; --k)
			free(res->data[k]);
		free(res->data);
		free(res);
		return NULL;
	}
	res->x = -1;
	res->y = -1;
	res->e_x = -1;
	res->e_y = -1;
	res->len = len;
	return res;
}

void dataset_free(dataset *dat) {
	free(dat->data[dat->x]);
	free(dat->data[dat->y]);
	if (dat->e_x != -1)
		free(dat->data[dat->e_x]);
	if (dat->e_y != -1)
		free(dat->data[dat->e_y]);
	free(dat->data);
	free(dat);
}

fit_conf *load_config() {
	FILE *fp = fopen("config.txt", "r");
	if (!fp)
		return NULL;
	fit_conf *res = malloc(sizeof(fit_conf));
	if (!res) {
		fclose(fp);
		return NULL;
	}
	fscanf(fp, "x = %i\n", &res->x);
	fscanf(fp, "y = %i\n", &res->y);
	fscanf(fp, "ex = %i\n", &res->e_x);
	fscanf(fp, "ey = %i\n", &res->e_y);
	fscanf(fp, "start = %i\n", &res->start);
	fscanf(fp, "laenge = %i\n", &res->len);
	fscanf(fp, "min_chng = %lf\n", &res->min_chng);
	fscanf(fp, "max_iter = %i\n", &res->max_iter);
	return res;
}

int write_config() {
	FILE *fp = fopen("config.txt", "w");
	if (!fp)
		return 1;
	fprintf(fp, "x = 0\n");
	fprintf(fp, "y = 1\n");
	fprintf(fp, "ex = 2\n");
	fprintf(fp, "ey = 3\n");
	fprintf(fp, "start = 0\n");
	fprintf(fp, "laenge = 0\n");
	fprintf(fp, "min_chng = 1e-20\n");
	fprintf(fp, "max_iter = 30\n");
	fclose(fp);
	return 1;
}


int fit_log_push(double *dat){
	 FILE *lgfile = fopen("fit_log.txt", "a");
	 	if(lgfile == NULL)return 1;
		fprintf(lgfile,"Fitdaten fuer Gerade f = m*x+n:\n");
		fprintf(lgfile,"m: %f +- %f\n",dat[0],dat[1]);
		fprintf(lgfile,"n: %f +- %f\n",dat[2],dat[3]);
		fprintf(lgfile,"Ausgewertet nach %.0f Iterationen, letzte Aenderung: %.5e\n",dat[6],dat[7]);
		fprintf(lgfile,"Guete des Fits (Chi square / r.o.f): %f\n",dat[8]);
		fprintf(lgfile,"Abbruchkriterium der Iteration: ");
		switch ((int)dat[5]){
		case 1:
			fprintf(lgfile,"Maximale Anzahl Iterationen erreicht.");
			break;
		case 2:
			fprintf(lgfile,"Letzte Aenderung klein genug.");
			break;
		case 3:
			fprintf(lgfile,"Daten haben sich eingependelt.");
			break;
		case 4:
			fprintf(lgfile,"Speicherzuweisung fehlgeschlagen.");
			break;
		}
		fprintf(lgfile,"\n\n");
		fclose(lgfile);
		return 0;
}

int fit_print(double *dat){
	printf("Fitdaten fuer Gerade f = m*x+n:\n");
	printf("m: %f +- %f\n",dat[0],dat[1]);
	printf("n: %f +- %f\n",dat[2],dat[3]);
	printf("Ausgewertet nach %.0f Iterationen, letzte Aenderung: %.5e\n",dat[6],dat[7]);
	printf("Guete des Fits (Chi square / r.o.f): %f\n",dat[8]);
	printf("Abbruchkriterium der Iteration: ");
	switch ((int)dat[5]){
	case 1:
		printf("Maximale Anzahl Iterationen erreicht.");
		break;
	case 2:
		printf("Letzte Aenderung klein genug.");
		break;
	case 3:
		printf("Daten haben sich eingependelt.");
		break;
	case 4:
		printf("Speicherzuweisung fehlgeschlagen.");
		break;
	}
	printf("\n");
	return 0;
}
