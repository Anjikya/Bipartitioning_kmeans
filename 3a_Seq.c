// ConsoleApplication2.cpp : Defines the entry point for the console application.
//


#include<stdio.h>
#include<stdlib.h>
#include<math.h>




void Bi_means(int dim, int n, float *data, int k, int* cluster_assign, float* cluster_center, int* cluster_size, int cluster1, int cluster2);
void cluster_assignment(int dim, int n, int k, float* data, float* cluster_center, int cluster1, int cluster2, int * cluster_assign);
float calc_distance(int dim, float* p1, float* p2, int p, int q);
void recalculate_centroids(int dim, int n, int k, float* data, int* cluster_assign, int* cluster_size, int cluster1, int cluster2,float*);
void calc_centroid(int k, int** points, int* cluster_size, float* data, int dim, int one_cluster, int new_cluster,float*);
float* new_cluster_center(int dim, int n, float* data, float * cluster_center, int cluster_count);

float* calc_farthest_in_cluster(int dim, int cluster, float *data, int * cluster_assign, int n, float* cluster_center);
int calc_ssd(int dim, int n, float* data, int cluster_count, float *cluster_center, int* cluster_assign);



float random()
{

	return (float)rand() / RAND_MAX;
}

float square_error(float* data, int* cluster_assign, float* cluster_center, int n, int k, int dim)
{	int i,j;
	float distance = 0.0;
	for (i = 0; i<k; i++)
	{
		for (j = 0; j<n; j++)
		{
			if (cluster_assign[j] == i)
			{
				distance += sqrt(calc_distance(dim, data, cluster_center, j*dim, i*dim));
			}
		}
	}
	return distance;
}


void main()
{
	int n = 10000;
	int k = 32;
	int dim = 3;
	int i,ii,jj, j;
	float *data = (float*)calloc(n*dim, sizeof(float));
	int *cluster_assign = (int*)calloc(n, sizeof(int));
	float *cluster_center = (float*)calloc(k*dim, sizeof(float));
	int *cluster_size = (int*)calloc(k, sizeof(int));
	float *point = (float*)calloc(dim, sizeof(float));
	float *point1 = (float*)calloc(dim, sizeof(float));
	int cluster1 = 0;
	int cluster2 = 1;
	int cluster_count = 1;
	int flag = 0;
	int val=0;
	srand(1);
	for (i = 0; i<n; i++)
		for ( j = 0; j<dim; j++)
		{
			data[i*dim + j] = random();
			//printf("%f\n", data[i*dim + j]);
		}

	//choosing first center randomly
	memcpy(&cluster_center[0], &data[0], dim * sizeof(float));

	//choosing centers of the cluster that are distant from each other 
	//initial_cluster_centers(data, cluster_center, dim, n, k);

	//kmeans clustering

	for (i = 0; i < n; i++)
		cluster_assign[i] = 0;
	
	cluster_size[0] = n;

	point = calc_farthest_in_cluster(dim,cluster1,data,cluster_assign,n,cluster_center);

	memcpy(&cluster_center[dim], &point[0], dim * sizeof(float));


	while (cluster_count != k)
	{
		Bi_means(dim, n, data, k, cluster_assign, cluster_center, cluster_size, cluster1, cluster2);
		cluster_count++;
		//for (i = 0; i < cluster_count; i++)
			//printf("cluster %d---size %d\n", i, cluster_size[i]);

		cluster1 = calc_ssd(dim, n, data, cluster_count, cluster_center, cluster_assign);

		for ( ii = 0; ii < cluster_count; ii++)
		{
			if (cluster_size[ii] == 1)
			{
				flag = 1;
				cluster_count--;
				point = new_cluster_center(dim, n, data, cluster_center, cluster_count);
				for ( jj = 0; jj < dim; jj++)
				{
					cluster_center[ii*dim + jj] = point[jj];
				}

			}
		}
			cluster2 = cluster_count;
			if (cluster_count == k)
				break;
			if (flag != 1)
			{
				while (1)
				{
					if (cluster_assign[val] == cluster1)
						break;

					val++;
				}

				memcpy(&cluster_center[cluster1*dim], &data[val*dim], dim * sizeof(float));
				point1 = calc_farthest_in_cluster(dim, cluster1, data, cluster_assign, n, cluster_center);
				memcpy(&cluster_center[cluster2*dim], &point1[0], dim * sizeof(float));
			}

	}

	/*for (i = 0; i < k; i++)
		for (j = 0; j < dim; j++)
		{
			printf("Centers of cluster %d------Coordinate %d is %f \n", i, j, cluster_center[i*dim + j]);
		}*/

	for (i = 0; i < k; i++)
		printf("Size of cluster%d is %d\n", i, cluster_size[i]);

	float sse = square_error(data, cluster_assign, cluster_center, n, k, dim);
	printf("SSE = %f\n", sse);

	
	getchar();

}






float* calc_farthest_in_cluster(int dim,int cluster, float *data, int * cluster_assign, int n, float* cluster_center)
{
	float *temp = (float*)calloc(dim, sizeof(float));
	int i, j;
	float dist = 0.0;
	float max = 0.0;
	int point;
	for (i = 0; i < n; i++)
	{
		if (cluster_assign[i] == cluster)
		{
			dist = calc_distance(dim, data, cluster_center, i*dim, cluster*dim);
			if (max < dist)
			{
				max = dist;
				point = i;
			}
		}
	}

	for (i = 0; i < dim; i++)
	{
		temp[i] = data[point*dim+i];
	}

	return temp;

}


int calc_ssd(int dim, int n, float* data, int cluster_count,float *cluster_center,int* cluster_assign)
{

	int i, j;
	int c;
	float *val = (float*)calloc(cluster_count, sizeof(float));
	int max = 0;
	for (i = 0; i < cluster_count; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (cluster_assign[j] == i)
			{
				val[i] += calc_distance(dim, data, cluster_center, j*dim, i*dim);
			}
		}
		if (max < val[i])
		{
			max = val[i];
			c = i;
		}
	}

	return c;
}

float* new_cluster_center(int dim, int n, float* data, float * cluster_center, int cluster_count)
{
	int i, j;
	float dist;
	float max = 0.0;
	int p;
	float* distances = (float*)calloc(n, sizeof(float));
	float* temp = (float*)calloc(dim, sizeof(float));
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < cluster_count; j++)
		{
			dist = sqrt(calc_distance(dim, data, cluster_center, i*dim, j*dim));
			if (max < dist)
			{
				max = dist;
				distances[i] = dist;
			}
		}
		max = 0.0;
	}
	max = 0.0;
	for (i = 0; i < n; i++)
	{
		if (max < distances[i])
		{
			max = distances[i];
			p = i;
		}
	}
	for (i = 0; i < dim; i++)
	{
		temp[i] = data[p*dim + i];
	}

	return temp;
}


void Bi_means(int dim, int n, float *data, int k, int* cluster_assign, float* cluster_center, int* cluster_size,int cluster1,int cluster2)
{
	
	float *center = (float*)calloc(2*dim, sizeof(float));
	int i;
	int ii, jj, kk;
	int rr;

	
	
	
		cluster_assignment(dim, n, k, data, cluster_center, cluster1,cluster2,cluster_assign);
		recalculate_centroids(dim,n,k, data, cluster_assign,cluster_size,cluster1,cluster2,cluster_center);
		
		for ( ii = 0; ii < dim; ii++)
		{
			cluster_center[cluster1*dim + ii] = center[ii];
		}
		for (ii = 0; ii < dim; ii++)
		{
			cluster_center[cluster2*dim + ii] = center[dim+ii];
		}
}


void cluster_assignment(int dim, int n, int k, float* data, float* cluster_center,int cluster1,int cluster2,int * cluster_assign)
{
	int i;
	int j;
	float min;
	int cluster = 0;
	int *cluster_assigned = (int*)calloc(n, sizeof(int));
	float distance_from_center1=0.0;
	float distance_from_center2=0.0;

	for (i = 0; i < n; i++)
	{
		if (cluster_assign[i] == cluster1)
		{
			distance_from_center1 = sqrt(calc_distance(dim, data, cluster_center, i*dim, cluster1*dim));
			distance_from_center2 = sqrt(calc_distance(dim, data, cluster_center, i*dim,cluster2*dim));
			if (distance_from_center1 > distance_from_center2)
				cluster_assign[i] = cluster2;
			else
				cluster_assign[i] = cluster1;
		}
	}



	
}

float calc_distance(int dim, float* p1, float* p2, int p, int q)
{

	int i;
	float value;
	float distance = 0;
	for (i = 0; i < dim; i++)
	{

		distance += pow((p1[p + i] - p2[q + i]), 2);
	}
	return distance;
}

void recalculate_centroids(int dim, int n, int k, float* data, int* cluster_assign, int* cluster_size,int cluster1,int cluster2,float* cluster_center)
{
	int **points = (int**)calloc(2, sizeof(int));
	float *centers = (float*)calloc(2*dim, sizeof(float));
	int i;
	int j;
	int x1 = 0;
	int x2 = 0;

	for (i = 0; i < 2; i++)
	{
		points[i] = (int*)calloc(n, sizeof(int));
	}
	
	
		for (j = 0; j < n; j++)
		{
			if (cluster_assign[j] == cluster1)
			{
				points[0][x1++] = j;
				
			}
			else if (cluster_assign[j] == cluster2)
			{
				points[1][x2++] = j;
				
			}

		}
		

		cluster_size[cluster1] = x1;
		
		cluster_size[cluster2] = x2;
		

		
	

	 calc_centroid(k, points, cluster_size, data, dim,cluster1,cluster2,cluster_center);



}

void calc_centroid(int k, int** points, int* cluster_size, float* data, int dim, int cluster1, int cluster2, float* cluster_center)
{
	//float *centercentroids = (float*)calloc(2*dim, sizeof(float));
	int i;
	int j;
	int p = 0;
	int q = 0;
	float *coordinates = (float*)calloc(2 * dim, sizeof(float));



	for (i = 0; i < dim; i++)
	{
		while (p != cluster_size[cluster1])
		{
			coordinates[i] += data[points[0][p] * dim + i];
			p++;
		}
		cluster_center[cluster1*dim + i] = coordinates[i] / cluster_size[cluster1];
		p = 0;
	}

	for (i = 0; i < dim; i++)
		coordinates[i] = 0.0;

	for (i = 0; i < dim; i++)
	{
		while (q != cluster_size[cluster2])
		{
			coordinates[i] += data[points[1][q] * dim + i];
			q++;
		}
		cluster_center[cluster2*dim + i] = coordinates[i] / cluster_size[cluster2];
		//printf("cluster %d with coordinate %d center is %f ", cluster2, i, cluster_center[cluster2*dim + i]);
		q = 0;
	}


}

