#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

void First_Bisecting_data(float *process_data_set, int *cluster_assign, int *cluster_size, float *cluster_center, float *local_centers, int cid, float *new_local_centers, int k, int n,int dim,int ktwo, int J, int Rank, int proc);
void Second_Bisecting_data(float *process_data_set, int *cluster_assign, int *cluster_size, float *cluster_center, float *local_centers, int arr[2], float *new_local_centers, int k, int n,int dim,int ktwo, int J, int Rank, int proc, int *global_cluster_size);
void bisectmean(float *process_data_set, int *cluster_assign, float *cluster_center, int *global_cluster_size, int *cluster_size, float *new_local_centers, float *local_centers, int n, int dim, int k , int ktwo, int J, int Rank, int proc);
void generateData(float *data, int dataSize, int dim, int k, int proc, int Rank);
float euclidean_distance(float *p1, float *p2, int dim);
void update_centers(float *data, float *cluster_center, int *cluster_assign, int *cluster_size, int k, int dim, int n);
float sse(float *data, int *cluster_assign, float *cluster_center,  int k, int dim, int n);
int find_max(float *buff, int size);
void global_update(float *global_centroids, int K, int DIM, int proc);
void kmeans(float *data, int *cluster_assign, float *cluster_center, int *cluster_size, int k, int n, int dim);
void kmeans1(float *data, int *cluster_assign, float *cluster_center, int *cluster_size, int k, int n, int dim);
void assign(float *data, float *cluster_center, int *cluster_assign, int *cluster_size, int k, int n, int dim);
int indexing(int *buff, int size , int value);
void max2(int *buff, int size, int arr[2]);



int ROOT = 0;

int main()
{

    int argc;
    char** argv;
	float sse_global;
    float *global_centroids;
    float *prev_global_centroids;
    int *global_cluster_assign;
    int *global_cluster_size;
	int dim = 8;  									 //dimension of the data 
	int n = 120000;									//total number of rows  in the data set
	int k = 32;     									//total number of cluster 
    float ss = 0.0;	
    int Rank;      
    int proc;											// Total Number of proc 
	int J=1;
    int dataSize;
    
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    
    dataSize = n/proc;
    float *data = malloc((dataSize * dim) * sizeof(float));
    int *cluster_assign = malloc(dataSize * sizeof(int));
    float *cluster_center = malloc((k * dim) * sizeof(float));
    int *cluster_size = malloc(k * sizeof(int));

  
    float *local_centers = malloc((2 *dim) * sizeof(float));
    float *new_local_centers = malloc((2 *dim) * sizeof(float));

    
    if(Rank == 0){
        global_centroids = malloc((k * dim) * sizeof(float));
        prev_global_centroids = malloc((k * dim) * sizeof(float));
        global_cluster_assign = malloc(n * sizeof(int));
        global_cluster_size = malloc((k * dim) * sizeof(float));
    }

    
    generateData(data, dataSize, dim, k, proc, Rank);

    
    bisectmean(data, cluster_assign, cluster_center, global_cluster_size, cluster_size, new_local_centers, local_centers, dataSize, dim, k, 2, J, Rank, proc);

    
    kmeans1(data, cluster_assign, cluster_center, cluster_size, k, dataSize, dim);
    ss = sse(data, cluster_assign, cluster_center, k, dim, dataSize);

   
    MPI_Reduce(cluster_center, global_centroids, (k*dim), MPI_FLOAT, MPI_SUM, ROOT, MPI_COMM_WORLD);
    MPI_Reduce(cluster_size, global_cluster_size, k, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&ss, &sse_global, 1, MPI_FLOAT, MPI_SUM, ROOT, MPI_COMM_WORLD);

    if(Rank == 0){
        int i,j;
        global_update(global_centroids, k, dim, proc);
        sse_global = sse_global / proc;

               
            for (i = 0; i < k; ++i) {
				printf("Cluster %d size is %d",i,global_cluster_size[i]);
                printf("\nCluster : %d centers are :", i);
            for(j=0; j<dim; j++){
                printf("%f\n", global_centroids[i*dim +j]);
            if(j != dim-1){
                printf("-> ");
            }
        }
      
        }
        printf("\nSSE: %f\n", sse_global );
    }

    if(Rank == ROOT){
        free(global_centroids);
        free(prev_global_centroids);
        free(global_cluster_assign);
        free(global_cluster_size);
    }
 	free(data);
    free(cluster_assign);
    free(cluster_center);
    free(cluster_size);
    free(local_centers);
    free(new_local_centers);
    MPI_Finalize();
    return 0;
}

void global_update(float *global_centroids, int K, int DIM, int proc) {

    int i, j ,d;
    int flag = 0;

   
    for (i = 0; i < K; ++i){
        float mean = 0.0;
        for (j = 0; j < DIM; ++j){
            mean = global_centroids[i*DIM +j] / (float) proc;
            global_centroids[i*DIM +j] = mean;
        }
    }
}


void generateData(float *data, int dataSize, int dim, int k, int proc, int Rank) {
    int i,p, j,seeuclidean_distance,size, num_of_seeuclidean_distances;
    int offset = 0;

    num_of_seeuclidean_distances = (12/proc);
    for (p = 0; p < num_of_seeuclidean_distances; ++p) {
        seeuclidean_distance = p*proc + (Rank+1); 
        srand(seeuclidean_distance);

        size = ((dataSize*dim)/num_of_seeuclidean_distances) * (p+1);
        for (i = offset; i < size; ++i) {
            data[i] = (float) rand() / (float) RAND_MAX;
        }
        offset = i;
    }
}


int find_max(float *buff, int size) {
    int i;
    int max_index = 0;
    float max = buff[0];
    for (i = 0; i < size; ++i) {
        if (max < buff[i]) {
            max = buff[i];
            max_index = i;
        }
    }
    return max_index;
}

void bisectmean(float *process_data_set, int *cluster_assign, float *cluster_center, int *global_cluster_size, int *cluster_size, float *new_local_centers,
                      float *local_centers, int n, int dim, int k , int ktwo, int J, int Rank, int proc) {
    int cid,i,j;
    cid = -1;
    int arr[2];
    int indices[2];

    while(J < k){

        if(J == 1){
         
            First_Bisecting_data(process_data_set, cluster_assign, cluster_size, cluster_center ,local_centers, cid, new_local_centers, k, n, dim, ktwo, J, Rank, proc);
            J++;
        }else{
            
            MPI_Reduce(cluster_size, global_cluster_size, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            if(Rank == 0){

                max2(global_cluster_size, J, arr);
                indices[0] = indexing(global_cluster_size, J , arr[0]);
                indices[1] = indexing(global_cluster_size, J , arr[1]);
            }
           
            MPI_Bcast(indices, 2, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(arr, 2, MPI_INT, 0, MPI_COMM_WORLD);

            if(arr[0]< 2*arr[1]){

               
                Second_Bisecting_data(process_data_set, cluster_assign, cluster_size, cluster_center ,local_centers, indices, new_local_centers, k, n, dim, ktwo, J, Rank, proc, global_cluster_size);
                J = J + 2;

            }else if(arr[0] >= 2*arr[1]){

               
                First_Bisecting_data(process_data_set, cluster_assign, cluster_size, cluster_center ,local_centers, indices[0], new_local_centers, k, n, dim, ktwo, J, Rank, proc);
                J++;
            }

        }

    }

}




int indexing(int *buff, int size , int value) {
    int i;
    for (i = 0; i < size; ++i) {
        if (buff[i] == value) {
            return i;
        }
    }
}

void Second_Bisecting_data(float *process_data_set, int *cluster_assign, int *cluster_size, float *cluster_center, float *local_centers,
                        int indices[2], float *new_local_centers, int k, int n,int dim,int ktwo, int J, int Rank, int proc,
                         int *global_cluster_size) {

    int i,j, cn1=0, cn2=0, c1, dataSize_size1, dataSize_size2, max_indx;
    int flag = 1;
    float assign1;
    float *data1, *data2;
    float *temp_local_centers;
    int r[2];
    float local_centers1[2*dim];
    float local_centers2[2*dim];

    
    memset(local_centers, 0, (dim*2) * sizeof(float));
    memset(local_centers1, 0, (dim*2) * sizeof(float));
    memset(local_centers2, 0, (dim*2) * sizeof(float));
    temp_local_centers = malloc( ((J+2)*dim) * sizeof(float));

   

    
    dataSize_size1 = cluster_size[indices[0]];
    data1 = malloc( (dataSize_size1* dim) * sizeof(float));
    for (i = 0; i < n; ++i) {
        if (indices[0] == cluster_assign[i]) { 
            memcpy(&data1[cn1*dim], &process_data_set[i*dim], dim*sizeof(float));
            cn1++;
        }
    }

    
    dataSize_size2 = cluster_size[indices[1]];
    data2 = malloc( (dataSize_size2* dim) * sizeof(float));
    for (i = 0; i < n; ++i) {
        if (indices[1] == cluster_assign[i]) { 
            memcpy(&data2[cn2*dim], &process_data_set[i*dim], dim*sizeof(float));
            cn2++;
        }
    }

    if(Rank == 0){
        srand(time(NULL));
        r[0] = rand()% dataSize_size1 + 0;
        memcpy(&local_centers1[0], &data1[r[0]*dim], dim * sizeof(float));
        r[1] = rand()% dataSize_size2 + 0;
        memcpy(&local_centers2[0], &data2[r[1]*dim], dim * sizeof(float));
    }
    
    MPI_Bcast(&local_centers1[0], dim, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_centers2[0], dim, MPI_FLOAT, 0, MPI_COMM_WORLD);

   
    float max = 0.0;
    int indx = 0;
    for (i = 0; i < dataSize_size1; ++i) {
        if((i == r[0]) && (Rank == 0)){
            continue;
        }else{
            assign1 = euclidean_distance(&local_centers1[0], &data1[i*dim], dim);
            if(assign1 > max){
                max = assign1;
                indx = i;
            }
        }
    }
    memcpy(&local_centers1[dim], &data1[indx*dim], dim * sizeof(float));

    max = 0.0;
    for (i = 0; i < dataSize_size2; ++i) {
        if((i == r[1]) && (Rank == 0)){
            continue;
        }else{
            assign1 = euclidean_distance(&local_centers2[0], &data2[i*dim], dim);
            if(assign1 > max){
                max = assign1;
                indx = i;
            }
        }
    }
    memcpy(&local_centers2[dim], &data2[indx*dim], dim * sizeof(float));

    c1 = 0;
    for (i = 0; i < J; ++i) {
        if(indices[0] != i && indices[1] != i){ 
            memcpy(&temp_local_centers[c1*dim], &cluster_center[i*dim], dim*sizeof(float));
            c1++;
        }
    }
  
    kmeans(data1, cluster_assign, local_centers1, cluster_size, 2, dataSize_size1, dim);
    kmeans(data2, cluster_assign, local_centers2, cluster_size, 2, dataSize_size2, dim);

    
    for (i = 0; i < 2; ++i) {
        memcpy(&temp_local_centers[c1*dim], &local_centers1[i*dim], dim * sizeof(float));
        c1++;
    }

    for (i = 0; i < 2; ++i) {
        memcpy(&temp_local_centers[c1*dim], &local_centers2[i*dim], dim * sizeof(float));
        c1++;
    }


    memcpy(&cluster_center[0], &temp_local_centers[0], (c1*dim)*sizeof(float));
    assign(process_data_set, cluster_center, cluster_assign, cluster_size, c1, n, dim);

    free(data1);
    free(data2);
    free(temp_local_centers);
}

void First_Bisecting_data(float *process_data_set, int *cluster_assign, int *cluster_size, float *cluster_center, float *local_centers,
                        int cid, float *new_local_centers, int k, int n,int dim,int ktwo, int J, int Rank, int proc) {

    int i,r,j, c, c1, dataSize_size, max_indx;
    int flag = 1;
    float assign1;
    float *data;
    float *temp_local_centers;
    float *recv_max;
    float maxi[proc];

    dataSize_size = cluster_size[cid];
    data = malloc( (dataSize_size* dim) * sizeof(float));
    if (Rank == 0) recv_max = malloc( (proc* dim) * sizeof(float));
    temp_local_centers = malloc( ((J+1)*dim) * sizeof(float));

   
    memset(local_centers, 0, (dim*2) * sizeof(float));
    memset(new_local_centers, 0, (dim*2) * sizeof(float));
    memset(temp_local_centers, 0, ((J+1)*dim) * sizeof(float));

  
    if(J == 1){
       
        if(Rank == 0){
            srand(time(NULL));
            r = rand()% n + 0;
            memcpy(&local_centers[0], &process_data_set[r*dim], dim * sizeof(float));
        }
        MPI_Bcast(&local_centers[0], dim, MPI_FLOAT, 0, MPI_COMM_WORLD);

      
        float max = 0.0;
        int indx = 0;
        for (i = 0; i < n; ++i) {
            if((i == r) && (Rank == 0)){
                continue;
            }else{
                assign1 = euclidean_distance(&local_centers[0], &process_data_set[i*dim], dim);
                if(assign1 > max){
                    max = assign1;
                    indx = i;
                }
            }
            
        }
        memcpy(&local_centers[dim], &process_data_set[indx*dim], dim * sizeof(float));

   
        if(Rank != 0){
            MPI_Send(&local_centers[dim], dim, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        }else{
            memcpy(&recv_max[0], &process_data_set[indx*dim], dim * sizeof(float));
            for (j = 1; j < proc; ++j) {
                MPI_Recv(&recv_max[j*dim], dim, MPI_FLOAT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            for (j = 0; j < proc; ++j) {
                maxi[j] = euclidean_distance(&local_centers[0], &recv_max[j*dim], dim);
            }

            max_indx = find_max(maxi, proc);
            memcpy(&local_centers[dim], &recv_max[max_indx*dim], dim * sizeof(float));
        }
      
        MPI_Bcast(&local_centers[dim], dim, MPI_FLOAT, 0, MPI_COMM_WORLD);
       
        kmeans(process_data_set, cluster_assign, local_centers, cluster_size, 2, n, dim);
      
        memcpy(&cluster_center[0], &local_centers[0], (2*dim)*sizeof(float));

    } else if(J > 1) {
        c = 0;
        for (i = 0; i < n; ++i) {
            if (cid == cluster_assign[i]) { 
                memcpy(&data[c*dim], &process_data_set[i*dim], dim*sizeof(float));
                c++;
            }
        }

        if(Rank == 0){
            srand(time(NULL));
            r = rand()% dataSize_size + 0;
            memcpy(&local_centers[0], &data[r*dim], dim * sizeof(float));
        }
        MPI_Bcast(&local_centers[0], dim, MPI_FLOAT, 0, MPI_COMM_WORLD);


        float max = 0.0;
        int indx = 0;
        for (i = 0; i < dataSize_size; ++i) {
            if((i == r) && (Rank == 0)){
                continue;
            }else{
                assign1 = euclidean_distance(&local_centers[0], &data[i*dim], dim);
                if(assign1 > max){
                    max = assign1;
                    indx = i;
                }
            }

        }
        memcpy(&local_centers[dim], &data[indx*dim], dim * sizeof(float));

        c1 = 0;
        for (i = 0; i < J; ++i) {
            if(cid != i){ 
                memcpy(&temp_local_centers[c1*dim], &cluster_center[i*dim], dim*sizeof(float));
                c1++;
            }
        }
        
        kmeans(data, cluster_assign, local_centers, cluster_size, 2, dataSize_size, dim);


        for (i = 0; i < 2; ++i) {
            memcpy(&temp_local_centers[c1*dim], &local_centers[i*dim], dim * sizeof(float));
            c1++;
        }
        memcpy(&cluster_center[0], &temp_local_centers[0], (c1*dim)*sizeof(float));
      
        assign(process_data_set, cluster_center, cluster_assign, cluster_size, c1, n, dim);
    }

    free(data);
    free(temp_local_centers);
        if(Rank == 0) free(recv_max);
}


float euclidean_distance(float *p1, float *p2, int dim) {
    float assign1 = 0.0;
    int i;

    for ( i = 0; i < dim; ++i) {
        assign1 += pow(p1[ i ] - p2[ i ], 2);
    }

    return assign1;
}


void assign(float *data, float *cluster_center, int *cluster_assign, int *cluster_size, int k, int n, int dim) {

    int i,j,l,q,p;
    float temp[k];
    int cid;
    float min;

  
    for (p = 0; p < k; ++p){
    	cluster_size[p] = 0;
    }
 
    for (p = 0; p < n; ++p){
        cluster_assign[p] = 0;
    }

    for (i = 0; i < n; i++) {
    
        for (j = 0; j < k; j++) {
            temp[j] = sqrt(euclidean_distance(&data[ i * dim ], &cluster_center[ j * dim], dim));
         
        }
     
        min = temp[0];
        cid= 0;
        for (l = 0; l < k; ++l) {
            if (min > temp[l]) {
                min = temp[l];
                cid = l;
            }
        }
     
        ++(cluster_size[cid]);
        cluster_assign[i] = cid;
    }
}


void update_centers(float *data, float *cluster_center, int *cluster_assign, int *cluster_size, int k, int dim, int n) {

	float *new_centroids = malloc((k * dim) * sizeof(float));  
    float temp[dim];
    float mean;
    int i, p, j, l, m, q;

  
    for (i = 0; i < k; ++i){
		mean = 0.0;
	   
	    for (p = 0; p < dim; ++p) {
	        temp[p] = 0.0;
	    }
	  
	    for (j = 0; j < n; ++j){
	    
	        if (i == cluster_assign[ j ]) {
	            for ( l = 0; l < dim; ++l) {
	                temp[l] += (data[j*dim + l]);
	              
	            }
	        }
	    }
	  
	    for (m = 0; m < dim; ++m){
	        mean = temp[ m ] / (float) cluster_size[ i ];
	        new_centroids[i*dim+m] = mean;
	    }
    }
 
    for (q = 0; q < k; ++q) {
        memcpy(&cluster_center[ q * dim ], &new_centroids[ q * dim ], dim * sizeof(float));
    }
    free(new_centroids);
}

void kmeans(float *data, int *cluster_assign, float *cluster_center, int *cluster_size, int k, int n, int dim){

    int flag = 0;
    int index = 1;
    float threshold = 0.000001;
    float val = 0.0;
    float *temp = malloc((k * dim) * sizeof(float));  

    while (1) {

    
    	++index;
    	int q, p, i, j, x;
    	for (q = 0; q < k; ++q) {
        	memcpy(&temp[q*dim ], &cluster_center[q*dim ], dim * sizeof(float));
    	}
    	
    	
        assign(data, cluster_center, cluster_assign, cluster_size, k, n, dim);

    	update_centers(data, cluster_center, cluster_assign, cluster_size, k, dim, n);
    
	  
    	float val;
        int c = 0;
        for (x = 0; x < k; ++x) {
            val = sqrt(euclidean_distance(&temp[ x * dim ], &cluster_center[ x * dim ], dim));
            if (val < pow(10.0, -6))
                c += 1;
        }

	
		if (c == k){
           
    		break;
		}
    }

    free(temp);
}

void kmeans1(float *data, int *cluster_assign, float *cluster_center, int *cluster_size, int k, int n, int dim){

    int flag = 0;
    int index = 1;
    float threshold = 0.000001;
    float val = 0.0;
    int cnt = 0;
    float *temp = malloc((k * dim) * sizeof(float));  

    while (cnt < 5) {

    
        ++index;
        int q, p, i, j, x;
        for (q = 0; q < k; ++q) {
            memcpy(&temp[q*dim ], &cluster_center[q*dim ], dim * sizeof(float));
        }
        
    
        assign(data, cluster_center, cluster_assign, cluster_size, k, n, dim);
       
        update_centers(data, cluster_center, cluster_assign, cluster_size, k, dim, n);
    
      
        float val;
        int c = 0;
        for (x = 0; x < k; ++x) {
            val = sqrt(euclidean_distance(&temp[ x * dim ], &cluster_center[ x * dim ], dim));
            if (val < pow(10.0, -6))
                c += 1;
        }

 
        cnt++;
    }
    free(temp);
}


float sse(float *data, int *cluster_assign, float *cluster_center,  int k, int dim, int n){

    float sum[k];
    int i1;
    int i;
    int j;
    int l;
    float sse = 0.0;


    for (i1 = 0; i1 < k; ++i1) {
        sum[i1] = 0.0;
    }

    for (i = 0; i < k; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == cluster_assign[j]) {
                sum[i] += euclidean_distance(&data[j*dim], &cluster_center[i*dim], dim);
            }
        }
    }

    for (l = 0; l < k; ++l) {
        sse += sum[l];
    }
    return sse;
}

void max2(int *buff, int size, int arr[2]) {

    int i, j, max, max1;
    int twoMax[2] = {-1,-1};
    int temp[size-1];
    int max_index = 0;
    int flag = -1;
    int c=0;

    max = buff[0];
    for (i = 0; i < size; ++i) {
        if (max < buff[i] ) {
            max = buff[i];
            max_index = i;
        }
    }
    twoMax[0] = max;

    for ( j = 0; j < size; ++j) {
        if(buff[j] != max){
            temp[c] = buff[j];
            c++;
        }
    }

    max = temp[0];
    for (i = 0; i < size-1; ++i) {
        if (max < temp[i] ) {
            max = temp[i];
        }
    }
    twoMax[1] = max;
    memcpy(&arr[0], &twoMax[0], 2 * sizeof(float));
}