/*
*Student Name: Gökhan UYSAL
*Student Number: 2017400024
*Compile Status: Compiling
*Program Status: Working
*Note: -
*/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "mpi.h"

/*
 * This function is used for calculating manhattan distance between two instance
 * @param A[] this is one of the instances
 * @param B[] this is one of the instances
 * @param size the size of arrays
 * @return double value which is result
 */
double manhattanDistance(double A[], double B[], int size) {
	double result = 0;
	for (int i = 0; i < size; i++) {
		result += fabs(A[i] - B[i]);
	}
	return result;
}

/*
 * This func is uset to calculate relief difference
 * @param target[] this is the target instance in the algorithm
 * @param source[] this is the nearestHit or nearestMiss instance of the algorithm
 * @param max[] this is the max of features named as max in the algorithm
 * @param min[] this is the min of features named as min in the algorithm
 * @param index this is the desired index for difference
 * @return double value which is calculated
 */
double diff(double target[], double source[], double max[], double min[],
		int index) {
	return fabs(target[index] - source[index]) / (max[index] - min[index]);
}

/*
 * this function is used to help qsort() function in the main function.
 * @param *a is one of the compared element in array
 * @param *b is one of the compared element in array
 * @return 1, -1 or 0 according to comparison
 */
int compare(const void *a, const void *b) {
	int f = *((int*) a);
	int s = *((int*) b);
	if (f > s) {
		return 1;
	} else if (f < s) {
		return -1;
	} else {
		return 0;
	}
}

int M; 		// Iteration number of the algorithm
int T;		// Top features count of the algorithm

int main(int argc, char *argv[]) {
	setvbuf(stdout, NULL, _IONBF, 0); // Disables the buffered printing to stdout
	int rank; // rank of the current processor
	int size; // total number of processors

	MPI_Init(&argc, &argv);					// initialization of the MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // gets the rank of the current processor
	MPI_Comm_size(MPI_COMM_WORLD, &size); // gets the total number of processors

	int P;	// Processor count
	int N;	// Instance coun
	int A;	// Feature count

	FILE *cin = fopen(argv[1], "r");	// opening file for input file
	fscanf(cin, "%d", &P); 	// reads number of processors from input file
	fscanf(cin, "%d", &N);	// reads number of instances from input file
	fscanf(cin, "%d", &A);	// reads number of features from input file

	double arr[(N * (A + 1)) + ((A + 1) * (N / (P - 1)))]; // stores preferences of all // root storage on master
	double pref[(A + 1) * (N / (P - 1))]; // stores preferences of each // local disk on processors

	// If it's master processor, reads from input file
	if (rank == 0) {

		fscanf(cin, "%d", &M);// reads number of iteration in the algorithm from input file
		fscanf(cin, "%d", &T);	// reads number of top features from input file
		double num = 0;			// for reading data
		int j = 0, i = ((A + 1) * (N / (P - 1)));	// loop iteration counts
		for (int i = 0; i < ((A + 1) * (N / (P - 1))); i++) {// Loop for filling 0's to arr[] array for scatter
			arr[i] = 0.0;
		}
		while (fscanf(cin, "%lf", &num) == 1) {	// Reading data from input file
			arr[i] = num;
			i++;
		}
		fclose(cin);	// closing file

		for (int i = 1; i < P; i++) {// Sending M and T values to slave processors
			MPI_Send(&M, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&T, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

	}

	// sends data from root array arr to pref array on each processor
	MPI_Scatter(arr, ((N / (P - 1)) * (A + 1)), MPI_DOUBLE, pref,
			((N / (P - 1)) * (A + 1)), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int masterSignal = 1; // Used for termination of the program
	while (masterSignal) {

		if (rank != 0) {		// If it's Slave processor
			int T = -1;		// initializing top features value
			int M = -1;		// initializing number of iterations

			/*
			 * Receiving M and T values from master processor
			 */
			MPI_Recv(&M, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&T, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			double max[A];				// maximum of features
			double min[A];				// minimum of features
			double nearestHit[A];		// nearestHit instance for the algorithm
			double nearestMiss[A];	// nearestMiss instance for the algorithm
			double targetInstance[A];	// Chosen instance for the algorithm
			double W[A];				// Weight array
			int result[T];				// Result indexes
			for (int i = 0; i < A; i++) {	// filling the weight array with 0's
				W[i] = 0;
			}
			int targetID = 0;				// Chosen instance's id
			double targetClass = 0;			// Chosen instance's target class
			for (int i = 0; i < A; i++) {// finding maximum and inimum of features
				double maxTemp = pref[i];// temp variable for storing temp maximums
				double minTemp = pref[i];// temp variable for storing temp minimums
				for (int j = i; j < ((N / (P - 1)) * (A + 1)); j += (A + 1)) {
					if (pref[j] > maxTemp) {
						maxTemp = pref[j];
					}
					if (pref[j] < minTemp) {
						minTemp = pref[j];
					}
				}
				max[i] = maxTemp;
				min[i] = minTemp;
			}

			for (int iteration = 0; iteration < M; iteration++) {// M iterations for the algorithm
				for (int j = 0; j < A; j++) {// Chosing target instance and filling target instance attributes
					if (j == A - 1) {
						targetClass = pref[j + 1 + iteration * (A + 1)];
					}
					targetInstance[j] = pref[j + iteration * (A + 1)];
				}
				targetID = iteration;
				double minNearestHit = INT_MAX;	// temp variable for nearestHit
				double minNearestMiss = INT_MAX;// temp variable for nearestMiss

				for (int i = 0; i < (N / (P - 1)); i++) {// for loop for finding nearestHit and nearestMiss
					if (targetID != i) {// if it is not chosen instance do the algorithm
						if (targetClass == pref[(i + 1) * (A + 1) - 1]) { // if the target class of instance is same wtih chosen instance target class fill nearestHit
							double temp[A];
							for (int j = 0; j < A; j++) {
								temp[j] = pref[j + i * (A + 1)];
							}
							double tempDistance = manhattanDistance(
									targetInstance, temp, A); // find manhattan distance
							if (tempDistance < minNearestHit) {
								minNearestHit = tempDistance;
								for (int k = 0; k < A; k++) {
									nearestHit[k] = temp[k];
								}
							}
						} else {		// if it is not same fill nearestMiss
							double temp[A];
							for (int j = 0; j < A; j++) {
								temp[j] = pref[j + i * (A + 1)];
							}
							double tempDistance = manhattanDistance(
									targetInstance, temp, A);// find manhattan distance
							if (tempDistance < minNearestMiss) {
								minNearestMiss = tempDistance;
								for (int k = 0; k < A; k++) {
									nearestMiss[k] = temp[k];
								}
							}
						}
					}
				}
				for (int i = 0; i < A; i++) {// Find the relief difference and fill the weights
					double leftSide = (diff(targetInstance, nearestHit, max,
							min, i) / M);	// difference for nearestHit
					double rightSide = (diff(targetInstance, nearestMiss, max,
							min, i) / M);	// difference for nearestMiss
					W[i] = W[i] - leftSide + rightSide;
				}
			}
			for (int i = 0; i < T; i++) {// finding the index of top T weights
				double tempMax = INT_MIN;
				int tempIndex = -1;
				for (int j = 0; j < A; j++) {
					if (W[j] > tempMax) {
						tempMax = W[j];
						tempIndex = j;
					}
				}
				result[i] = tempIndex;
				W[tempIndex] = INT_MIN;
			}

			/*
			 * sorting the result[] array for printing with the help of compare function at the above
			 */
			qsort(result, sizeof(result) / sizeof(*result), sizeof(*result),
					compare);

			int sync = 1; // sync message for printing sync
			if (rank != 1) {// if it is not slave 1 wait for message to recieve
				MPI_Recv(&sync, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);
			}
			/*
			 * printing results to stdout
			 */
			printf("Slave P%d :", rank);
			for (int i = 0; i < T; i++) {
				printf(" %d", result[i]);
			}
			printf("\n");
			if (rank != P - 1) {// if it is not last slave send message to next slave for sync
				MPI_Send(&sync, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			}

			/*
			 * sending results to master slave
			 */
			MPI_Send(result, T, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}

		if (rank == 0) {	// if it is master processor
			masterSignal = 0;// changing masterSignal value to 0 for prepearing to terminate

			int size = T * (P - 1);
			int temp[T];	// array for storing the receieved results in loop
			int result[size];// array for storing the receieved result in general.
			for (int i = 1; i < P; i++) {// loop for receiving results messageg from slaves
				MPI_Recv(temp, T, MPI_INT, i, 0, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);
				for (int j = 0; j < T; j++) {
					result[(i - 1) * T + j] = temp[j];
				}
			}

			/*
			 * sorting result[] array for printing the results
			 */
			qsort(result, sizeof(result) / sizeof(*result), sizeof(*result),
					compare);

			for (int i = 0; i < size; i++) {// Eliminating the duplicates in the array
				for (int j = i + 1; j < size; j++) {
					/* If any duplicate found */
					if (result[i] == result[j]) {
						/* Delete the current duplicate element */
						for (int k = j; k < size; k++) {
							result[k] = result[k + 1];
						}

						/* Decrement size after removing duplicate element */
						size--;

						/* If shifting of elements occur then don't increment j */
						j--;
					}
				}
			}

			/*
			 * printing results to stdout
			 */
			printf("Master P0 :");
			for (int i = 0; i < size; i++) {
				printf(" %d", result[i]);
			}
			printf("\n");

		}

		MPI_Bcast(&masterSignal, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast to terminate successfully

	}

	// ****************************************** //

	MPI_Barrier(MPI_COMM_WORLD); // synchronizing processes
	MPI_Finalize();

	return 0;
}
