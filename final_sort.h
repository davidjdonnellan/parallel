/**
 * MPI Program that performs simple sorting
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

double * merge_array(int n, double * a, int m, double * b);
void     merge_sort(int n, double * a);
void     swap (double * a, double * b);


int MPI_Exchange( int n, double * array, int rank1, int rank2, MPI_Comm comm,double *compTime, double * commTime);
int MPI_Sort_direct(int n, double * array, int root, MPI_Comm comm); 
int MPI_Sort_ranking(int n, double * array, int root, MPI_Comm comm); 
int MPI_Sort_bucket(int n, double * arrary, double max, int root, MPI_Comm comm);
int MPI_Sort_OddEven( int n, double * array, int root, MPI_Comm comm);
int isSorted(int n, double * a,int root, MPI_Comm comm);



int MPI_Exchange( int n, double * array, int rank1, int rank2, MPI_Comm comm,double *compTime, double * commTime ) {
      
      int rank, size, result, i, tag1 = 0, tag2 = 1, tmpTime=0,commTime1=0,compTime1=0;

      double * b = ( double * ) calloc( n, sizeof( double ) );
      double * c;
       
      MPI_Status status;
      MPI_Comm_rank( comm, &rank );
      MPI_Comm_size( comm, &size );
     
      //L8.6
      if( rank == rank1 )
      {
        result = MPI_Send( &array[ 0 ], n, MPI_DOUBLE, rank2, tag1, comm );
        result = MPI_Recv( &b[ 0 ], n, MPI_DOUBLE, rank2, tag2, comm, &status );
        c = merge_array( n, array, n, b );
        for( i = 0; i < n; i++ )
        {
          array[ i ] = c[ i ];
        }
      }
      else if( rank == rank2 )
      {
        result = MPI_Recv( &b[ 0 ], n, MPI_DOUBLE, rank1, tag1, comm, &status );
        result = MPI_Send( &array[ 0 ], n, MPI_DOUBLE, rank1, tag2, comm) ;
        c = merge_array( n, array, n, b );
        for( i =0; i < n; i++ )
        {
          array[ i ] = c[ i + n ];
        }
      }

      *compTime = compTime1;
      *commTime = commTime1;
      //printf("Exchange: Processor %d generates %lf communication time and %lf computation time \n", rank, commTime, compTime);
      return MPI_SUCCESS;
    }


int isSorted(int n, double *a, int root, MPI_Comm comm){
	//find rank and size
	double tmpTime = 0, commTime = 0, compTime = 0;
	int rank, size, i, answer=1;
	MPI_Comm_rank( comm, &rank );
	MPI_Comm_size( comm, &size );

	//allocate space for first and last
	double *first = (double *)calloc(size, sizeof(double));
	double *last = (double *)calloc(size, sizeof(double));
    
    	tmpTime = MPI_Wtime();
	//gather first and last
	MPI_Gather(&a[0], 1, MPI_DOUBLE, first, 1, MPI_DOUBLE, root, comm);
	MPI_Gather(&a[n-1], 1, MPI_DOUBLE, last, 1, MPI_DOUBLE, root, comm);
   	tmpTime = MPI_Wtime() + tmpTime;
   	 commTime += tmpTime;
	
   	 tmpTime = MPI_Wtime();
	//if root, test
	if(rank == root){
		for(i=0; i<size-1;i++){
			if(last[i] > first[i+1])
			{
				answer = 0;			
			}		
		}	
	}
	tmpTime = MPI_Wtime() - tmpTime;
    	compTime += tmpTime;
    
    	tmpTime = MPI_Wtime();
	MPI_Bcast(&answer, 1, MPI_INT, root, comm);
    	tmpTime = MPI_Wtime() + tmpTime;
   	 commTime += tmpTime;
    
	return answer;
}

int MPI_Sort_bucket(int n, double * a, double max, int root, MPI_Comm comm){
   int rank, size, i, j, error; 
	// find rank and size
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    double commTime = 0, compTime = 0, tmpTime;
	
	// allocate the extra memory / arrays needed
    double * bucket = (double *) calloc(n, sizeof(double));
    int bucketCount = 0;
    
    int * bucketCounts = (int *) calloc(size, sizeof(int));
    int * displs = (int *) calloc(size, sizeof(int));
    
    tmpTime = MPI_Wtime();
	// Brodcast the array to the processor
    error = MPI_Bcast(a, n, MPI_DOUBLE, root, comm);
    if(error!=MPI_SUCCESS) return error;
    tmpTime = MPI_Wtime() - tmpTime;
    commTime += tmpTime;
    
    tmpTime = MPI_Wtime();
    //scan a and collect b elements of bucket rank
    for(i=0; i<n;i++){
        if (a[i] >= rank*max/size && a[i]<(rank+1)*max/size){
            bucket[bucketCount++] = a[i];
            
        }
    }
    // merge sort bucket
    merge_sort(bucketCount, bucket);
    tmpTime = MPI_Wtime() - tmpTime;
    compTime += tmpTime;
	

    tmpTime = MPI_Wtime();
    //Gather-v the buckets to a
    //gather the bucketCounta and calculate the displavement
    MPI_Gather(&bucketCount, 1, MPI_INT, bucketCounts, 1, MPI_INT, root, comm);
    tmpTime = MPI_Wtime() - tmpTime;
    commTime += tmpTime;
    

    tmpTime = MPI_Wtime();
    if(rank==root){
        displs[0] = 0;
        for(i=1;i<size;i++){
            displs[i] = displs[i-1]+bucketCounts[i-1];
        }
    }
    tmpTime = MPI_Wtime() - tmpTime;
    compTime += tmpTime;
    

    tmpTime = MPI_Wtime();
    error = MPI_Gatherv(bucket, bucketCount,MPI_DOUBLE, a,
        bucketCounts ,displs,MPI_DOUBLE, root, comm);
    tmpTime = MPI_Wtime() - tmpTime;
    commTime += tmpTime;
    

    if(error != MPI_SUCCESS){
        return error;
    }
   

    printf("Bucket Sort: Processor %d generates %lf communication and %lf computation \n", rank, commTime, compTime);
    
    return MPI_SUCCESS;
	
}
//pasted this in from the SimpleSort file
int MPI_Sort_ranking(int n, double * a, int root, MPI_Comm comm)
{
	int rank,size,i,j;
	
	// find rank and size
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	
	double commTime=0, compTime=0, tmpTime=0;
	
	// allocate the extra memory / arrays needed
	int * ranking = (int *) calloc(n/size, sizeof(int));
	int * finalRanking = (int *) calloc(n, sizeof(int));
	double * b = (double *) calloc(n/size, sizeof(double));
	
	tmpTime = MPI_Wtime();
	// Brodcast the array to the processor
	MPI_Bcast(a, n, MPI_DOUBLE, root, comm);
	tmpTime = MPI_Wtime() - tmpTime;
	commTime += tmpTime;
	
	tmpTime = MPI_Wtime();
	// P rank generates an array ranking with ranking[i] is the rank of a[i+rank*n/size] in the array
	for(i=0;i<(n/size);i++){
		ranking[i] = 0;
		for(j=0;j<n;j++){
			if(a[i+rank*(n/size)]>a[j]){
				ranking[i]++;
			}
		}
	}
	tmpTime = MPI_Wtime() - tmpTime;
	compTime += tmpTime;


	tmpTime = MPI_Wtime();
	// Gather the array ranking to finalRanking
	MPI_Gather(ranking, n/size, MPI_INT, finalRanking, n/size, MPI_INT, root, comm);
	tmpTime = MPI_Wtime() - tmpTime;
	commTime += tmpTime;
	
	tmpTime = MPI_Wtime();
	// if processor 0 then restore the order in the array b and move b back to a
	if(rank==root){
		for(i=0;i<n;i++){
			b[finalRanking[i]] = a[i];
		}
		for(i=0;i<n;i++)
			{a[i] = b[i];}
	
	}

	tmpTime = MPI_Wtime() - tmpTime;
	compTime += tmpTime;
	printf("Rank Process %d generate %lf communication and %lf computation\n", rank, commTime, compTime);
	
	return MPI_SUCCESS;
	
	
}
int MPI_Sort_OddEven( int n, double * array, int root, MPI_Comm comm ) {

      double compTime = 0, commTime = 0, tmpTime1 = 0, tmpTime2 = 0;
      // get rank and size of comm
	 int rank, size, i, j, error;
    
	// find rank and size
    	MPI_Comm_rank(comm, &rank);
   	MPI_Comm_size(comm, &size);
      
      //allocate space for numElements/numProcessors amount of doubles
	double * localA = (double*)calloc(n/size, sizeof(double));
         
      //scatter a to local_a
	tmpTime1 = MPI_Wtime();
	MPI_Scatter(array, n/size, MPI_DOUBLE, localA, n/size, MPI_DOUBLE, root, comm);
        tmpTime1 = MPI_Wtime() - tmpTime1;
	commTime += tmpTime1;


      //sort local_a using mergeSort
	tmpTime2 = MPI_Wtime();
	merge_sort(n/size,localA);
	tmpTime2 = MPI_Wtime() - tmpTime2;
	compTime += tmpTime2;
       
      //odd-even iterations
	for( i=0;i<size;i++){
		// exhcange with rank+1
		if((rank+1)%2==0){
			if(rank<size-1){
				MPI_Exchange(n/size,localA,rank,rank+1,comm,&tmpTime1,&tmpTime2);}
				compTime += tmpTime1; 
				commTime += tmpTime2;
		}
		// exchange with rank-1
		else{
			if(rank>0){
				MPI_Exchange(n/size,localA,rank-1,rank,comm,&tmpTime1,&tmpTime2);}
				compTime += tmpTime1;
				commTime+= tmpTime2;
		
		}
		// you must now call isSorted	
		//isSorted(n, localA,root,comm);
	}
       //gather local_a
    	MPI_Gather(localA,n/size,MPI_DOUBLE,array,n/size,MPI_DOUBLE,root,comm);

	//printf("Odd Even Sort: Processor %d generates %lf communication and %lf computation \n", rank, commTime, compTime);
	return MPI_SUCCESS;
      
}
int MPI_Sort_direct(int n, double * array, int root, MPI_Comm comm){
	
	int rank, size, error;
    double  commTime=0, compTime=0, tmpTime = 0;
	double * localArray;

	//get rank and size
	MPI_Comm_size (comm, &size);
	MPI_Comm_rank (comm, &rank);

	//allocate localArray
	localArray = (double *) calloc(n/size, sizeof(double));

	// scatter array to local array
    tmpTime = MPI_Wtime();
	error = MPI_Scatter(array, n/size, MPI_DOUBLE, localArray, n/size, MPI_DOUBLE, root, comm);
	if(error != MPI_SUCCESS) return error;
    tmpTime = MPI_Wtime() - tmpTime;
	commTime+= tmpTime;

	//proc rank sorts local array
	tmpTime = MPI_Wtime();
	merge_sort(n/size, localArray);
	tmpTime = MPI_Wtime() - tmpTime;
	compTime += tmpTime;

	//gather local arrays
    tmpTime = MPI_Wtime();
	error = MPI_Gather(localArray, n/size, MPI_DOUBLE, array, n/size, MPI_DOUBLE, root, comm);
	if(error != MPI_SUCCESS) return error;
    tmpTime = MPI_Wtime() - tmpTime;
	commTime += tmpTime;
	
	// if proc 0 then merge size - 1 times the chunk in there
    tmpTime = MPI_Wtime();
	if(rank == 0){
		for(int i=1;i<size;i++){
			localArray = merge_array(i*n/size, array, n/size, array+i*n/size);
			double * tmp = merge_array(i*n/size, array, n/size, array+i*n/size);
			for(int j = 0; j < (i+1)*n/size; j++){array[j] = tmp[j];}
		}
	}	
	tmpTime = MPI_Wtime() - tmpTime;
	compTime+= tmpTime;
     printf("Direct Sort: Processor %d generates %lf communication and %lf computation \n", rank, commTime, compTime);

	return MPI_SUCCESS;
}
double * merge_array(int n, double * a, int m, double * b){

   int i,j,k;
   double * c = (double *) calloc(n+m, sizeof(double));

   for(i=j=k=0;(i<n)&&(j<m);)

      if(a[i]<=b[j])c[k++]=a[i++];
      else c[k++]=b[j++];

   if(i==n)for(;j<m;)c[k++]=b[j++];
   else for(;i<n;)c[k++]=a[i++];

return c;
}
// function to merge sort the array a with n elements
void merge_sort(int n, double * a){

   double * c;
   int i;

   if (n<=1) return;

   if(n==2) {

      if(a[0]>a[1])swap(&a[0],&a[1]);
      return;
   }

   merge_sort(n/2,a);merge_sort(n-n/2,a+n/2);

   c=merge_array(n/2,a,n-n/2,a+n/2);

   for(i=0;i<n;i++)a[i]=c[i];

return;
}
// swap two doubles
void swap (double * a, double * b){

   double temp;

   temp=*a;*a=*b;*b=temp;

}

