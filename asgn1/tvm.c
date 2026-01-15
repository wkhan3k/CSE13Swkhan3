/*
 * tvm:  Time Value of Money
 * Assignment 1 of CSE 13S, Winter 2026.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * The financial calculator's variables:
 */
double n, i, PV, PMT, FV;

static void die_line(int line_number, const char *msg) {
    fprintf(stderr, "line %d: %s\n", line_number, msg);
    exit(1);
}

static int is_integer_double(double x) {
    return isfinite(x) && floor(x) == x;
}

/*
 * Purpose:      Function f(n) used by the Newton-Raphson root finder.
 *
 * Parameter:    Number of periods n.
 *
 *               See equation (13) in the assignment PDF.
 *
 * Returns:      f(n)
 */
double fn(double n_param) {
    // Use the rearranged form that avoids negative exponents:
    // f(x) = (PV + PMT/i) * (1+i)^x - PMT/i + FV
    double a = 1.0 + i;
    double term = PV + (PMT / i);
    return term * pow(a, n_param) - (PMT / i) + FV;
}

/*
 * Purpose:      Function f(i) used by the Newton-Raphson root finder.
 *
 * Parameter:    Interest rate i.
 *
 *               See equation (13) in the assignment PDF.
 *
 * Returns:      f(i)
 */
double fi(double i_param) {
    double a = 1.0 + i_param;
    double term = PV + (PMT / i_param);
    return term * pow(a, n) - (PMT / i_param) + FV;
}

/*
 * Purpose:      Function f'(n) used by the Newton-Raphson root finder.
 *
 * Parameter:    Number of periods n.
 *
 *               See equation (14) in the assignment PDF.
 *
 * Returns:      f'(n)
 */
double fn_prime(double n_param) {
    // f'(n) = ln(1+i) * (PV + PMT/i) * (1+i)^n
    double a = 1.0 + i;
    return log(a) * (PV + (PMT / i)) * pow(a, n_param);
}

/*
 * Purpose:      Function f'(i) used by the Newton-Raphson root finder.
 *
 * Parameter:    Interest rate i.
 *
 *               See equation (15) in the assignment PDF.
 *
 * Returns:      f'(i)
 */
double fi_prime(double i_param) {
    // f'(i) = n*(PV+PMT/i)*(1+i)^(n-1) - PMT * ((1+i)^n + 1)/i^2
    double a = 1.0 + i_param;
    double a_to_n = pow(a, n);
    double first = n * (PV + (PMT / i_param)) * pow(a, n - 1.0);
    double second = PMT * (a_to_n + 1.0) / (i_param * i_param);
    return first - second;
}

/*
 * Purpose:      Find x such that fn(x) == 0 using a simple Newton-Raphson
 *               root finder.
 *
 * Parameters:   line_number for error messages.
 *
 * Returns:      x such that fn(x) == 0.
 */
double newton_raphson_n(int line_number) {
    // initial guess from PDF
    double x = 360.0;

    // Require i positive (division by i, log(1+i))
    if (!(i > 0.0)) {
        die_line(line_number, "i must be positive");
    }

    for (int iter = 0; iter < 100000; iter++) {
        double f = fn(x);
        double fp = fn_prime(x);

        if (!isfinite(f) || !isfinite(fp) || fp == 0.0) {
            die_line(line_number, "solver did not converge");
        }

        double delta = f / fp;
        x -= delta;

        if (!isfinite(x)) {
            die_line(line_number, "solver did not converge");
        }

        if (fabs(delta) < 1e-8) {
            if (x <= 0.0) {
                die_line(line_number, "solver did not converge");
            }
            return ceil(x);
        }
    }

    die_line(line_number, "solver did not converge");
    return 0.0;
}

/*
 * Purpose:      Find x such that fi(x) == 0 using a simple Newton-Raphson
 *               root finder.
 *
 * Parameters:   line_number for error messages.
 *
 * Returns:      x such that fi(x) == 0.
 */
double newton_raphson_i(int line_number) {
    // initial guess from PDF
    double x = 0.0025;

    for (int iter = 0; iter < 100000; iter++) {
        // x must stay positive because we divide by x and use (1+x)
        if (!(x > 0.0)) {
            die_line(line_number, "solver did not converge");
        }

        double f = fi(x);
        double fp = fi_prime(x);

        if (!isfinite(f) || !isfinite(fp) || fp == 0.0) {
            die_line(line_number, "solver did not converge");
        }

        double delta = f / fp;
        x -= delta;

        if (!isfinite(x)) {
            die_line(line_number, "solver did not converge");
        }

        if (fabs(delta) < 1e-8) {
            if (!(x > 0.0)) {
                die_line(line_number, "solver did not converge");
            }
            return x;
        }
    }

    die_line(line_number, "solver did not converge");
    return 0.0;
}

/*
 * Purpose:      Check whether i is 0.0.  If it is, report an error ("i must
 *               be positive"), and exit the program with an exit code of 1.
 */
void check_i(int line_number) {
    if (!(i > 0.0)) {
        die_line(line_number, "i must be positive");
    }
}

/*
 * Purpose:      Compute and print a variable's new value.
 */
void tvm_compute_variable(char *name, int line_number) {
    if (strcmp(name, "n") == 0) {
        n = newton_raphson_n(line_number);
        printf("n = %.0f\n", n);
        return;
    }

    if (strcmp(name, "i") == 0) {
        i = newton_raphson_i(line_number);
        printf("i = %.6f\n", i);
        return;
    }

    if (strcmp(name, "PV") == 0) {
        check_i(line_number);
        // PV = - PMT * (1 - (1+i)^(-n))/i - FV*(1+i)^(-n)
        double a = 1.0 + i;
        double disc = pow(a, -n);
        PV = -PMT * (1.0 - disc) / i - FV * disc;
        printf("PV = %.2f\n", PV);
        return;
    }

    if (strcmp(name, "PMT") == 0) {
        check_i(line_number);
        // PMT = i(PV(1+i)^n + FV) / (1 - (1+i)^n)
        double a = 1.0 + i;
        double a_to_n = pow(a, n);
        double denom = 1.0 - a_to_n;
        if (denom == 0.0) {
            die_line(line_number, "solver did not converge");
        }
        PMT = i * (PV * a_to_n + FV) / denom;
        printf("PMT = %.2f\n", PMT);
        return;
    }

    if (strcmp(name, "FV") == 0) {
        check_i(line_number);
        // FV = -PV(1+i)^n - PMT * ((1+i)^n - 1)/i
        double a = 1.0 + i;
        double a_to_n = pow(a, n);
        FV = -PV * a_to_n - PMT * (a_to_n - 1.0) / i;
        printf("FV = %.2f\n", FV);
        return;
    }

    die_line(line_number, "invalid variable name");
}

/*
 * Purpose:      Set a variable's new value.
 */
void tvm_set_variable(char *name, double value, int line_number) {
    if (strcmp(name, "n") == 0) {
        if (!(value > 0.0) || !is_integer_double(value)) {
            die_line(line_number, "n must be a positive integer");
        }
        n = value;
        return;
    }

    if (strcmp(name, "i") == 0) {
        if (!(value > 0.0)) {
            die_line(line_number, "i must be positive");
        }
        i = value;
        return;
    }

    if (strcmp(name, "PV") == 0) {
        PV = value;
        return;
    }

    if (strcmp(name, "PMT") == 0) {
        PMT = value;
        return;
    }

    if (strcmp(name, "FV") == 0) {
        FV = value;
        return;
    }

    die_line(line_number, "invalid variable name");
}

/*
 * Purpose:      Clear all of the financial variables:  n, i, PV, PMT, and FV.
 */
void tvm_clear(void) {
    n = 0.0;
    i = 0.0;
    PV = 0.0;
    PMT = 0.0;
    FV = 0.0;
}

/*
 * Purpose:      Process one of these comnands, or print an error and exit.
 *
 *                   "set VAR NUMBER"
 *                   "compute VAR"
 *                   "clear"
 */
void tvm_process_command(char *command, int line_number) {
    // Ignore empty lines
    if (command[0] == '\0') return;

    // Make a writable copy for strtok
    char buf[40];
    strncpy(buf, command, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';

    char *tok = strtok(buf, " \t");
    if (tok == NULL) return;

    if (strcmp(tok, "clear") == 0) {
        // no extra tokens allowed
        if (strtok(NULL, " \t") != NULL) {
            die_line(line_number, "invalid command");
        }
        tvm_clear();
        return;
    }

    if (strcmp(tok, "compute") == 0) {
        char *var = strtok(NULL, " \t");
        if (var == NULL) {
            die_line(line_number, "invalid command");
        }
        if (strtok(NULL, " \t") != NULL) {
            die_line(line_number, "invalid command");
        }
        tvm_compute_variable(var, line_number);
        return;
    }

    if (strcmp(tok, "set") == 0) {
        char *var = strtok(NULL, " \t");
        char *num = strtok(NULL, " \t");
        if (var == NULL || num == NULL) {
            die_line(line_number, "invalid command");
        }
        if (strtok(NULL, " \t") != NULL) {
            die_line(line_number, "invalid command");
        }

        char *end = NULL;
        double val = strtod(num, &end);
        if (end == num || *end != '\0') {
            die_line(line_number, "invalid command");
        }

        tvm_set_variable(var, val, line_number);
        return;
    }

    die_line(line_number, "invalid command");
}

/*
 * Purpose:      Truncate string s at the first occurance of '\n'.
 */
void truncate_at_newline(char *s) {
    int i = 0;

    while (s[i] != '\0' && s[i] != '\n') ++i;

    s[i] = '\0';
}

/*
 * Purpose:      Read commands from stdin, and process them.
 * Exit code:    1 on error, else 0
 */
int main(void) {
    char command[40];

    int line_number = 0;

    tvm_clear();

    while (1) {
        char *res = fgets(command, sizeof(command), stdin);

        if (res == NULL) break;

        ++line_number;
        truncate_at_newline(command);
        tvm_process_command(command, line_number);
    }

    if (ferror(stdin)) {
        fprintf(stderr, "tvm:  Error reading input\n");
        exit(1);
    }

    return 0;
}