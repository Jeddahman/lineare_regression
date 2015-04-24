#include <stdio.h>
#include <stdlib.h>
#include "fit_formulas.h"

int main(int cargs, char **vargs) {
	/* Argumente in Reihenfolge: dateipfad, indikatoren (zulaessig sind x,y,ex,ey)
	 * indikator gibt reihenfolge und auftauchen der Daten - wenn nichts uebergeben wird geht man vom standard der cfg aus
	 * auch ein Ueberspringen von x fehler oder y fehler ist moeglich, genau so wie das komplette auslassen dieser
	 * Beispiel fuer standard: x y ex ey
	 */
	double *fit;
	char *file;
	int k;
	fit_conf *config;
	if (cargs < 2)
		file = "werte.txt";
	else
		file = vargs[1];
	config = load_config();
	if (!config)
		write_config();
	else
		config = load_config();
	if (cargs > 2)
		for (k = 2; k < cargs; ++k) {
			if (vargs[k][0] == 'x')
				config->x = k - 2;
			else if (vargs[k][0] == 'y')
				config->y = k - 2;
			else {
				if (vargs[k][0] == 'e') {
					if (vargs[k][1] == 'x')
						config->e_x = k - 2;
					else
						config->e_y = k - 2;
				}
			}
		}
	fit = linear_file_fit(file, config);
	fit_print(fit);
	fit_log_push(fit);
	free(fit);
	free(config);
	return 0;
}
