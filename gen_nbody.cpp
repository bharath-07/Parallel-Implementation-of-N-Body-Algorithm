#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>

#define MAX_X 800.0
#define MAX_Y 800.0
#define MAX_Z 800.0
#define MAX_M 100.0

using namespace std;


int main (int argc, char **argv) {
    if (argc < 2) {
        printf ("Insufficient number of arguments. Please try again.\n");
        return 0;
    }

    int n = atoi (argv[1]);

    FILE *fp;
    char filename[15];
    sprintf(filename, "%d.txt", n);
    fp = fopen (filename, "w");

    srand (time(NULL));

    fprintf (fp, "%d\n", n);
    for (int i = 0; i < n; i++) {
        double x = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/MAX_X));
        double y = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/MAX_Y));
        double z = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/MAX_Z));
        double m = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/MAX_M));

        int num;
        num = rand() % 2 + 1;
        if (num == 1)
            x = -x;
        num = rand() % 2 + 1;
        if (num == 1)
            y = -y;
        num = rand() % 2 + 1;
        if (num == 1)
            z = -z;

        fprintf (fp, "%lf %lf %lf %lf\n", x, y, z, m);
    }

    fclose(fp);

    return 0;
}
