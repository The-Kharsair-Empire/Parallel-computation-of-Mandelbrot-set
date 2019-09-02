#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


int main(int argc, char* argv[]){

	/* screen ( integer) coordinate */
	int iX,iY;
	const int iXmax = 8000; // default
	const int iYmax = 8000; // default

	/* world ( double) coordinate = parameter plane*/
	double Cx, Cy;
	const double CxMin = -2.5;
	const double CxMax = 1.5;
	const double CyMin = -2.0;
	const double CyMax = 2.0;

	/* */
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 
	FILE * fp;
	// char *filename = "Mandelbrot_Parallel.ppm";
	char filename[256] = {0};
	char *comment = "# ";

	// RGB color array
	static unsigned char color[3];
	unsigned char *colorArray = NULL;
	/* Z = Zx + Zy*i;	Z0 = 0 */
	double Zx, Zy;
	double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	/*  */
	int Iteration;
	const int IterationMax = 2000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;
	
	/* Clock information */
	int rank, size, rows_per_procs, rows_remain, rows_to_divide_evenly;
	const int masterThread = 0;
	const int tag = 0;
	double start, end, time, threadStart, threadEnd, threadTime;
	int responsible_row;
	int rc;

	rc = MPI_Init(&argc, &argv);

	MPI_Status stat;

	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	start = threadStart = MPI_Wtime();
	MPI_Comm_size( MPI_COMM_WORLD, &size);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank);

	if (rank == masterThread) {
		// main master opens the file , and calculate some parameter that will later be broadcast to other threads
		/*create new file,give it a name and open it in binary mode  */
		sprintf(filename, "Mandelbrot_Parallel_with_%d_processes.ppm", size);
		fp = fopen(filename, "wb"); /* b -  binary mode */
		/*write ASCII header to the file (PPM file format)*/
		fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);

		printf("File: %s successfully opened for writing.\n", filename);
		printf("Computing Mandelbrot Set. Please wait...\n");
		rows_per_procs = iYmax / size;
		rows_remain = iYmax % size;
	}

	MPI_Bcast(&rows_per_procs, 1, MPI_INT, masterThread, MPI_COMM_WORLD);
	MPI_Bcast(&rows_remain, 1, MPI_INT, masterThread, MPI_COMM_WORLD);
	// printf("rank %d recieved rows_per_procs: %d, rows_remain: %d\n", rank, rows_per_procs, rows_remain);
	int pixel_pos;
	int cur_row = 0;
	responsible_row = rank;

	if (rank == masterThread) {
		// printf("rank 0 here, i am rank: %d, i will calculate %d rows", rank, rows_per_procs+1);
		colorArray = (unsigned char*)malloc(iXmax*iYmax*3*sizeof(char));
	} else if (rank < rows_remain){
		// printf("rank %d here, i am < %d, i will calculate %d rows\n", rank, rows_remain, rows_per_procs+1);
		colorArray = (unsigned char*)malloc(iXmax*3*(rows_per_procs+1)*sizeof(char));
	} else {
		// printf("rank %d here, i am >= %d, i will calculate %d rows\n", rank, rows_remain, rows_per_procs);
		colorArray = (unsigned char*)malloc(iXmax*3*rows_per_procs*sizeof(char));
	}
	//test
	// int row_cal = 0;
	//test
	while (responsible_row < iYmax) {
		//test
		// printf("rank: %d 's responsible row: %d\n", rank, responsible_row);
		// row_cal += 1;
		//test	
		Cy = CyMin + (responsible_row * PixelHeight);
		if (fabs(Cy) < (PixelHeight / 2)) {
			Cy = 0.0; /* Main antenna */
		}
		pixel_pos = 0;
		for (iX = 0; iX < iXmax; iX++) {
			Cx = CxMin + (iX * PixelWidth);
			/* initial value of orbit = critical point Z= 0 */
			Zx = 0.0;
			Zy = 0.0;
			Zx2 = Zx * Zx;
			Zy2 = Zy * Zy;

			/* do calculation */
			for(Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++) {
				Zy = (2 * Zx * Zy) + Cy;
				Zx = Zx2 - Zy2 + Cx;
				Zx2 = Zx * Zx;
				Zy2 = Zy * Zy;
			};

			if (Iteration == IterationMax) {
				// Point within the set. Mark it as black
				colorArray[cur_row+pixel_pos] = 0;
				colorArray[cur_row+pixel_pos+1] = 0;
				colorArray[cur_row+pixel_pos+2] = 0;
			}
			else {
				// Point outside the set. Mark it as white
				double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0); //how far away it is from inside the set.// aka. escape velocity
				if (c < 1) {
					colorArray[cur_row+pixel_pos] = 0;
					colorArray[cur_row+pixel_pos+1] = 0;
					colorArray[cur_row+pixel_pos+2] = 255*c;
				}
				else if (c < 2) {
					colorArray[cur_row+pixel_pos] = 0;
					colorArray[cur_row+pixel_pos+1] = 255*(c-1);
					colorArray[cur_row+pixel_pos+2] = 255;
				}
				else {
					colorArray[cur_row+pixel_pos] = 255*(c-2);
					colorArray[cur_row+pixel_pos+1] = 255;
					colorArray[cur_row+pixel_pos+2] = 255;
				}
			}
			pixel_pos += 3;
		}
		cur_row += iXmax*3;
		responsible_row += size;
		
	}
	// printf("i am rank: %d, my row_cal is: %d\n", rank, row_cal);
	
	if (rank == masterThread) {
		
		int offset;
		if (masterThread < rows_remain){
			offset = (rows_per_procs+1) * iXmax * 3;
		} else {
			offset = (rows_per_procs) * iXmax * 3;
		}
		 
		for (int i = 1; i < size; i++){
			if (i < rows_remain) {
				// printf("on offset - 1: %c\n", colorArray[offset-1]);
				// printf("on offset: %c\n", colorArray[offset]);
				MPI_Recv(colorArray+offset, (rows_per_procs+1)*3*iXmax, MPI_UNSIGNED_CHAR, i, tag, MPI_COMM_WORLD, &stat);
				offset += (rows_per_procs+1) * iXmax * 3;
			} else {
				// printf("on offset - 1: %c\n", colorArray[offset-1]);
				// printf("on offset: %c\n", colorArray[offset]);
				MPI_Recv(colorArray+offset, (rows_per_procs)*3*iXmax, MPI_UNSIGNED_CHAR, i, tag, MPI_COMM_WORLD, &stat);
				offset += (rows_per_procs) * iXmax * 3;
			}
		}
		end = MPI_Wtime();
		time = end - start;
		printf("Mandelbrot computational time: %lf\n", time);
		int nextRow, jump, nextPixel;
		_Bool notFinished = 1;
		for (int block = 0; block < rows_per_procs+1 && notFinished; block ++) {
			nextRow = block * iXmax * 3;
			for (int each_proc = 0; each_proc < size; each_proc ++){
				nextPixel = 0;
				if (block == rows_per_procs && each_proc >= rows_remain){
					notFinished = 0;
					break;
				}
				if (each_proc < rows_remain) jump = rows_per_procs + 1;
			 	else jump = rows_per_procs;
	
				for (int row = 0; row < iXmax; row ++){
					color[0] = colorArray[nextRow+nextPixel];
					color[1] = colorArray[nextRow+nextPixel+1];
					color[2] = colorArray[nextRow+nextPixel+2];
					fwrite(color, 1, 3, fp);
					nextPixel += 3;
				}
				nextRow += jump * iXmax * 3;
				
			}
		}
		fclose(fp);
	} else {
		if ( rank < rows_remain) {
			MPI_Send(colorArray, (rows_per_procs+1)*3*iXmax, MPI_UNSIGNED_CHAR, masterThread, tag, MPI_COMM_WORLD);
		} else {
			MPI_Send(colorArray, rows_per_procs*3*iXmax, MPI_UNSIGNED_CHAR, masterThread, tag, MPI_COMM_WORLD);
		}
	}
	free(colorArray);
	threadEnd = MPI_Wtime();
	threadTime = threadEnd-threadStart;
	printf("Rank %d took %lf time to finish its task, (if it is rank 0, it will include time of writing files)\n", rank, threadTime);

	MPI_Finalize();
	return 0;

}