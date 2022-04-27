#include <iostream>
#include <stack>
#include <vector>
#include <math.h>
#include <algorithm>
#include <omp.h>

// 6.67 x 10 ^ (-11)
#define g 0.0667

using namespace std;

const double theta = 0.5;
double max_x = 0.0, max_y = 0.0, max_z = 0.0;

struct Node {
    Node* child[8];

    double x;
    double y;
    double z;

    double m;

    // Data node = true implies a leaf node
    bool isDataNode;

    vector<Node*> objectsAtSamePoint;

    omp_lock_t lock;
};

int whichOctant (double x, double y, double z, double originX, double originY, double originZ) {
    if (x - originX >= 0 && y - originY >= 0 && z - originZ >= 0)
        return 0;
    else if (x - originX < 0 && y - originY >= 0 && z - originZ >= 0)
        return 1;
    else if (x - originX < 0 && y - originY < 0 && z - originZ >= 0)
        return 2;
    else if (x - originX >= 0 && y - originY < 0 && z - originZ >= 0)
        return 3;
    else if (x - originX >= 0 && y - originY >= 0 && z - originZ < 0)
        return 4;
    else if (x - originX < 0 && y - originY >= 0 && z - originZ < 0)
        return 5;
    else if (x - originX < 0 && y - originY < 0 && z - originZ < 0)
        return 6;
    else if (x - originX >= 0 && y - originY < 0 && z - originZ < 0)
        return 7;
    return 0;
}

void setOriginAndRange (double& originX, double& originY, double& originZ,
                        double& rangeX, double& rangeY, double& rangeZ,
                        int octant) {
    if (octant == 0) {
        originX = originX + rangeX / 2;
        originY = originY + rangeY / 2;
        originZ = originZ + rangeZ / 2;
    } else if (octant == 1) {
        originX = originX - rangeX / 2;
        originY = originY + rangeY / 2;
        originZ = originZ + rangeZ / 2;
    } else if (octant == 2) {
        originX = originX - rangeX / 2;
        originY = originY - rangeY / 2;
        originZ = originZ + rangeZ / 2;
    } else if (octant == 3) {
        originX = originX + rangeX / 2;
        originY = originY - rangeY / 2;
        originZ = originZ + rangeZ / 2;
    } else if (octant == 4) {
        originX = originX + rangeX / 2;
        originY = originY + rangeY / 2;
        originZ = originZ - rangeZ / 2;
    } else if (octant == 5) {
        originX = originX - rangeX / 2;
        originY = originY + rangeY / 2;
        originZ = originZ - rangeZ / 2;
    } else if (octant == 6) {
        originX = originX - rangeX / 2;
        originY = originY - rangeY / 2;
        originZ = originZ - rangeZ / 2;
    } else {
        originX = originX + rangeX / 2;
        originY = originY - rangeY / 2;
        originZ = originZ - rangeZ / 2;
    }

    rangeX = rangeX / 2;
    rangeY = rangeY / 2;
    rangeZ = rangeZ / 2;
}

void getLock (Node* start) {
    // Get lock
    omp_set_lock(&start->lock);
}

void releaseLock (Node* start) {
    // Release lock
    omp_unset_lock(&start->lock);
}


void insert (Node* start, double x, double y, double z, double m) {
    double originX = 0.0, originY = 0.0, originZ = 0.0;
    double rangeX = max_x, rangeY = max_y, rangeZ = max_z;

    // [LOCK] : Obtain lock on start
    getLock (start);

    while (true) {
        int octant = whichOctant (x, y, z, originX, originY, originZ);

        // If child of the current node is not null
        // And the child is a data node (which means its a leaf node)
        // Then there is a collision between the 2 points and we have to
        // expand the tree further to resolve the collision
        if (start->child[octant] != NULL && start->child[octant]->isDataNode) {

            // if another object exists at the same point
            if (start->child[octant]->x == x
                && start->child[octant]->y == y
                && start->child[octant]->z == z) {
                Node* node = new Node();
                node->x = x;
                node->y = y;
                node->z = z;
                node->m = m;
                omp_init_lock(&node->lock);
                node->isDataNode = true;
                start->child[octant]->objectsAtSamePoint.push_back(node);

                // [LOCK] : Release lock on start
                releaseLock (start);

                break;
            }

            start->child[octant]->isDataNode = false;

            int oldX = start->child[octant]->x;
            int oldY = start->child[octant]->y;
            int oldZ = start->child[octant]->z;
            int oldM = start->child[octant]->m;

            int octantOfNewPoint = whichOctant(x, y, z, originX, originY, originZ);
            int octantOfOldPoint = whichOctant(oldX, oldY, oldZ,  originX, originY, originZ);

            double centreOfMass = m + oldM;
            double centreOfMassX = ((m * x) + (oldM * oldX)) / centreOfMass;
            double centreOfMassY = ((m * y) + (oldM * oldY)) / centreOfMass;
            double centreOfMassZ = ((m * z) + (oldM * oldZ)) / centreOfMass;

            // Update x y and z of the current node
            // by adding the x,y and z of the new node to the current node
            start->x = ((start->x * start->m) + (x * m)) / (start->m + m);
            start->y = ((start->y * start->m) + (y * m)) / (start->m + m);
            start->z = ((start->z * start->m) + (z * m)) / (start->m + m);

            // Update centre of mass of the current node
            start->m = start->m + m;

            while (octantOfNewPoint == octantOfOldPoint) {

                if (start->child[octantOfNewPoint] == NULL) {
                    start->child[octantOfNewPoint] = new Node();
                    omp_init_lock(&start->child[octantOfNewPoint]->lock);
                }

                start->child[octantOfNewPoint]->isDataNode = false;

                // Update originX, originY and originZ of the new octant that was just created
                setOriginAndRange(originX, originY, originZ, rangeX, rangeY, rangeZ, octantOfNewPoint);

                Node* temp = start->child[octantOfNewPoint];
                // [LOCK] : Release lock on start
                releaseLock (start);
                // [LOCK] : Obtain lock on start->child[octantOfNewPoint]
                getLock (temp);

                start = temp;

                // Update x, y, z and centre of mass of new octant
                start->x = centreOfMassX;
                start->y = centreOfMassY;
                start->z = centreOfMassZ;
                start->m = centreOfMass;

                // Update octant of the new point with respect to the newly created octant
                octantOfNewPoint = whichOctant(x, y, z, originX, originY, originZ);

                // Update octant of the old point with respect to the newly created octant
                octantOfOldPoint = whichOctant(oldX, oldY, oldZ, originX, originY, originZ);
            }

            // Set new point to the correct octant
            start->child[octantOfNewPoint] = new Node();
            start->child[octantOfNewPoint]->x = x;
            start->child[octantOfNewPoint]->y = y;
            start->child[octantOfNewPoint]->z = z;
            start->child[octantOfNewPoint]->m = m;
            omp_init_lock(&start->child[octantOfNewPoint]->lock);
            start->child[octantOfNewPoint]->isDataNode = true;

            // Set old point to the correct octant
            start->child[octantOfOldPoint] = new Node();
            start->child[octantOfOldPoint]->x = oldX;
            start->child[octantOfOldPoint]->y = oldY;
            start->child[octantOfOldPoint]->z = oldZ;
            start->child[octantOfOldPoint]->m = oldM;
            omp_init_lock(&start->child[octantOfOldPoint]->lock);
            start->child[octantOfOldPoint]->isDataNode = true;

            // [LOCK] : Release lock on start
            releaseLock (start);

            break;
        }

            // If the octant is not null and it is not a data node
            // This implies that the octant is not a leaf node
            // So we need to go further down the tree
        else if (start->child[octant] != NULL) {

            // Update x y and z of the current node
            // by adding the x,y and z of the new node to the current node
            start->x = ((start->x * start->m) + (x * m)) / (start->m + m);
            start->y = ((start->y * start->m) + (y * m)) / (start->m + m);
            start->z = ((start->z * start->m) + (z * m)) / (start->m + m);

            start->m = (start->m + m);

            Node* temp = start->child[octant];
            // [LOCK] : Release lock on start
            releaseLock (start);
            // [LOCK] : Obtain lock on start->child[octant]
            getLock (temp);

            start = temp;
        }

            // Else the octant is empty
            // In this case, we can simply create a leaf node in the octant
            // And add out point to that octant
        else {
            start->child[octant] = new Node();
            start->child[octant]->x = x;
            start->child[octant]->y = y;
            start->child[octant]->z = z;
            start->child[octant]->m = m;
            start->child[octant]->isDataNode = true;
            omp_init_lock(&start->child[octant]->lock);

            // Update x y and z of the current node
            // by adding the x,y and z of the new node to the current node
            start->x = ((start->x * start->m) + (x * m)) / (start->m + m);
            start->y = ((start->y * start->m) + (y * m)) / (start->m + m);
            start->z = ((start->z * start->m) + (z * m)) / (start->m + m);
            start->m = start->m + m;

            // [LOCK] : Release lock on start
            releaseLock (start);

            break;
        }
    }
}

void getLockForForce (int i, omp_lock_t* locks) {
    // Get lock
    omp_set_lock(&locks[i]);
}

void releaseForForce (int i, omp_lock_t* locks) {
    // Release lock
    omp_unset_lock(&locks[i]);
}

void calculateNBody (double x, double y, double z, double m,
                     Node* root, double widthOfRegion, int id,
                     double* force_x, double* force_y, double* force_z,
                     omp_lock_t* locks) {
    if (root == NULL)
        return ;

    if (root->isDataNode) {
        // calculate #1

        double d = pow (x - root->x, 2) + pow (y - root->y, 2) + pow (z - root->z, 2);
        d = sqrt (d);
        if (d != 0) {
            getLockForForce (id, locks);
            force_x[id] = force_x[id] + ((g * m * root->m * (root->x - x)) / pow (d, 3));
            force_y[id] = force_y[id] + ((g * m * root->m * (root->y - y)) / pow (d, 3));
            force_z[id] = force_z[id] + ((g * m * root->m * (root->z - z)) / pow (d, 3));
            releaseForForce(id, locks);
        }

        for (int i = 0; i < root->objectsAtSamePoint.size();i++) {
            //calculate #2

            double d = pow (x - root->objectsAtSamePoint[i]->x, 2)
                       + pow (y - root->objectsAtSamePoint[i]->y, 2)
                       + pow (z - root->objectsAtSamePoint[i]->z, 2);

            d = sqrt (d);

            if (d != 0) {
                getLockForForce (id, locks);
                force_x[id] = force_x[id]
                              + ((g * m * root->objectsAtSamePoint[i]->m * (root->objectsAtSamePoint[i]->x - x)) / pow (d, 3));
                force_y[id] = force_y[id]
                              + ((g * m * root->objectsAtSamePoint[i]->m * (root->objectsAtSamePoint[i]->y - y)) / pow (d, 3));
                force_z[id] = force_z[id]
                              + ((g * m * root->objectsAtSamePoint[i]->m * (root->objectsAtSamePoint[i]->z - z)) / pow (d, 3));
                releaseForForce(id, locks);
            }
        }
    } else {
        double d = pow (x - root->x, 2) + pow (y - root->y, 2) + pow (z - root->z, 2);
        d = sqrt (d);
        double s = widthOfRegion;

        if (d == 0)
            return ;
        if (s / d <= theta) {
            // calculate here

            // calculate
            double d = pow (x - root->x, 2) + pow (y - root->y, 2) + pow (z - root->z, 2);
            d = sqrt (d);
            if (d != 0) {
                getLockForForce (id, locks);
                force_x[id] = force_x[id] + ((g * m * root->m * (root->x - x)) / pow (d, 3));
                force_y[id] = force_y[id] + ((g * m * root->m * (root->y - y)) / pow (d, 3));
                force_z[id] = force_z[id] + ((g * m * root->m * (root->z - z)) / pow (d, 3));
                releaseForForce(id, locks);
            }

        } else {
            for (int i = 0; i < 8; i++) {
                if (root->child[i] != NULL) {
                    calculateNBody (x, y, z, m, root->child[i], widthOfRegion/2, id, force_x, force_y, force_z, locks);
                }
            }
        }
    }
}


int main (int argc, char **argv) {

    if (argc < 2) {
        printf ("Insufficient number of arguments. Please try again.\n");
        return 0;
    }

    // OPEN FILE
    FILE* fp = fopen (argv[1], "r");

    if (fp == NULL) {
        printf ("Failed to open file. Please try again.\n");
        return 0;
    } else {
        int n;
        fscanf (fp, "%d\n", &n);

        double* fx = (double *)malloc(n * sizeof(double));
        double* fy = (double *)malloc(n * sizeof(double));
        double* fz = (double *)malloc(n * sizeof(double));
        double* fm = (double *)malloc(n * sizeof(double));

        double* force_x = (double *)malloc(n * sizeof(double));
        double* force_y = (double *)malloc(n * sizeof(double));
        double* force_z = (double *)malloc(n * sizeof(double));

        omp_lock_t* force_locks = (omp_lock_t *)malloc(n * sizeof(omp_lock_t));

        // Read the coordinates and the masses of the n bodies from the input file
        for (int i = 0; i < n; i++) {
            fscanf (fp, "%lf %lf %lf %lf\n", &fx[i], &fy[i], &fz[i], &fm[i]);

            force_x[i] = 0.0;
            force_y[i] = 0.0;
            force_z[i] = 0.0;

            omp_init_lock(&force_locks[i]);

            if (max_x < abs(fx[i]))
                max_x = abs(fx[i]);

            if (max_y < abs(fy[i]))
                max_y = abs(fy[i]);

            if (max_z < abs(fz[i]))
                max_z = abs(fz[i]);
        }

        max_x = 2 * max_x;
        max_y = 2 * max_y;
        max_z = 2 * max_z;

        max_x = max (max_x, max(max_y, max_z));
        max_y = max_x;
        max_z = max_x;

        // Initialize requisite data structures
        // Number of chunks
        int chunks = atoi (argv[2]);
        int t = atoi (argv[3]);
        int m = chunks;
        // Root of the tree
        Node* root = (Node *)malloc(m * sizeof(Node));
        for (int i = 0; i < m; i++) {
            root[i].isDataNode = false;
            root[i].x = 0.0;
            root[i].y = 0.0;
            root[i].z = 0.0;
            root[i].m = 0.0;
            omp_init_lock(&root[i].lock);
            for (int j = 0; j < 8; j++) {
                root[i].child[j] = NULL;
            }
        }

        int chunkSize = n / chunks;

        omp_set_max_active_levels(2);

        double calc_time = 0;
        double insert_time = 0;

        #pragma omp parallel proc_bind(spread) num_threads(chunks)
        {
            int tid = omp_get_thread_num();

            double start = omp_get_wtime();

            // Insert each of the n points into the tree
            #pragma omp parallel proc_bind(master) num_threads(t)
            #pragma omp for
            for (int i = 0; i < chunkSize; i++) {
                int id = (tid * chunkSize) + i;
                if (id < n)
                    insert (&root[tid], fx[id], fy[id], fz[id], fm[id]);
            }

            double end = omp_get_wtime();
            #pragma omp critical
            insert_time = max(insert_time, end - start);

            start = end;

            // Calculate the result for the n bodies
            #pragma omp for
            for (int i = 0; i < n; i++) {
                calculateNBody (fx[i], fy[i], fz[i], fm[i], &root[tid], max_x, i, force_x, force_y, force_z, &force_locks[0]);
            }

            end = omp_get_wtime();
            #pragma omp critical
            calc_time = max(calc_time, end - start);
        }
    }

    // CLOSE FILE
    fclose(fp);

    return 0;
}
