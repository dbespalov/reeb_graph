/******************************************************************************
 *  Topological similarity estimation for 3D models (triangle meshes)         *
 *    using multiresolutional Reeb graphs.                                    *
 *                                                                            *
 *  Copyright (C) 2004  Drexel University                                     *
 *  Implemented by:  Dmitriy Bespalov (bespalov@gmail.com)                    *
 *                                                                            *
 *    This file is part of reeb_graph.                                        *
 *    reeb_graph is free software: you can redistribute it and/or modify      *
 *    it under the terms of the GNU General Public License as published by    *
 *    the Free Software Foundation, either version 2 of the License, or       *
 *    (at your option) any later version.                                     *
 *                                                                            *
 *    reeb_graph is distributed in the hope that it will be useful,           *
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *    GNU General Public License for more details.                            *
 *                                                                            *
 *    You should have received a copy of the GNU General Public License       *
 *    along with reeb_graph.  If not, see <http://www.gnu.org/licenses/>.     *
 *****************************************************************************/

//Class with array manipulation static routines

public class aman{
    
    public static int size_inc = 2000;
    public static int size_inc_sparse = 20;
    public static int size_inc_base_points = 200;
    public static int size_inc_tsets = 50;
    public static int size_inc_queue = 100;
    public static int size_inc_edges = 1000;
    
    public static Point[] cut_ends(Point pts[], int points_size){
	Point temp[], points[];
	points = pts;
	temp = points;
	points = new Point[points_size];
	for(int i = 0; i < points_size; i++)
	    points[i] = temp[i];
	
	return points;
    }
    
    public static Triangle[] cut_ends(Triangle trigs[], int triangles_size){
	Triangle temp2[], triangles[];
	triangles = trigs;
	temp2 = triangles;
	triangles = new Triangle[triangles_size];
	for(int i = 0; i < triangles_size; i++)
	    triangles[i] = temp2[i];
	
	return triangles;
    }
    
    public static Point[] expand(Point pts[], int points_size){
	Point temp[], points[];
	points = pts;
	temp = points;
	points = new Point[temp.length + size_inc];
	for(int i = 0; i < points_size; i++)
	    points[i] = temp[i];
	
	return points;
    }
        
    public static Triangle[] expand(Triangle trigs[], int triangles_size){
	Triangle temp[], triangles[];
	triangles = trigs;
	temp = triangles;
	triangles = new Triangle[temp.length + size_inc];
	for(int i = 0; i < triangles_size; i++)
	    triangles[i] = temp[i];
	
	return triangles;
    }
    
    public static int[] expand(int arr[], int array_size){
	int temp[], array[];
	array = arr;
	temp = array;
	
	array = new int[temp.length + size_inc];
	for(int i = 0; i < array_size; i++)
	    array[i] = temp[i];
		
	return array;
    }
    
    public static double[] expand(double arr[], int array_size){
	double temp[], array[];
	array = arr;
	temp = array;
	
	array = new double[temp.length + size_inc];
	for(int i = 0; i < array_size; i++)
	    array[i] = temp[i];
		
	return array;
    }
    
    public static int[][] expand_sparse(int sprs[][], int sparse_size){
	int temp[][], sparse[][];
	sparse = sprs;
	temp = sparse;
	
	sparse = new int[temp.length + size_inc][];
	for(int i = 0; i < sparse_size; i++)
	    sparse[i] = temp[i];
	
	return sparse;
    }
    
    public static int[] expand_sparse(int sprs[], int sparse_size){
	int temp[], sparse[];
	sparse = sprs;
	temp = sparse;
	sparse = new int[temp.length + size_inc_sparse];
	for(int i = 0; i < temp[0]; i++)
	    sparse[i] = temp[i];
	
	return sparse;
    }
    
    public static Point[] expand_base(Point pts[], int points_size){
	Point temp[], points[];
	points = pts;
	temp = points;
	points = new Point[temp.length + size_inc_base_points];
	for(int i = 0; i < points_size; i++)
	    points[i] = temp[i];
	
	return points;
    }
    
    public static double[] expand_base(double area[], int area_size){
	double temp[], base_areas[];
	base_areas = area;
	temp = base_areas;
	base_areas = new double[temp.length + size_inc_base_points];
	for(int i = 0; i < area_size; i++)
	    base_areas[i] = temp[i];
	
	return base_areas;
    }
    
    public static GValueStorage[][] expand_base(GValueStorage g_vals[][], int base_size){
	
	GValueStorage temp[][], g_values[][];
	g_values = g_vals;
	temp = g_values;
	
	g_values = new GValueStorage[temp.length + size_inc_base_points][];
	for(int i = 0; i < base_size; i++)
	    g_values[i] = temp[i];
	    	
	return g_values;
    }
    
    public static int[][] expand_tsets(int one_range_tsets[][], int length){
	int tsets[][], temp[][];
	
	tsets = one_range_tsets;
	temp = tsets;
	int incr = size_inc_tsets;
	
	if(one_range_tsets.length >= size_inc_tsets*3)
	    incr = size_inc_tsets * 10;
	
	tsets = new int[temp.length+incr][];
	
	for(int i = 0; i < length; i++){
	    tsets[i] = temp[i];
	}
	
	return tsets;	
    }
    
    public static int[] expand_queue(int q[]){
	int queue[], temp[];
	queue = q;
	temp = queue;
	
	int len = (temp[1] - temp[0]) + 3;
	
	queue = new int[len + size_inc_queue];
	queue[0] = 2;
	queue[1] = 2;
	//System.out.println("queue.length = " + queue.length);
	for(int i = temp[0]; i < temp[1]; i++){
	    
	    //System.out.println("queue[1] = " + queue[1]);
	    
	    queue[queue[1]] = temp[i];
	    queue[1]++;	    
	}
	
	return queue;
    }
    
    public static int[][] expand_edges(int eds[][], int length){
	int edges[][], temp[][];
	edges = eds;
	temp = edges;
	
	edges = new int[temp.length + size_inc_edges][2];
	for(int i = 0; i < length; i++){
	    edges[i][0] = temp[i][0];
	    edges[i][1] = temp[i][1];
	}
	    
	return edges;
    }
    
    public static int[][] expand_triangle_edges(int eds[][], int length){
	int edges[][], temp[][];
	edges = eds;
	temp = edges;
	
	edges = new int[temp.length + size_inc_edges][3];
	for(int i = 0; i < length; i++){
	    edges[i][0] = temp[i][0];
	    edges[i][1] = temp[i][1];
	    edges[i][2] = temp[i][2];
	}
	    
	return edges;
    }
    
    public static int[][] trimSparse(int sprs[][], int sparse_size){
	int temp[][], sparse[][];
	sparse = sprs;
	temp = sparse;
	
	for(int i = 0; i < sparse_size; i++){
	    int array[] = sparse[i];
	    int count = 0;
	    for(int j = 0; j < array[0]; j++){
		if(array[j] > -1)
		    count++;	
	    }
	    
	    int new_array[] = new int[count + size_inc_sparse];
	    new_array[0] = 1;
	    for(int j = 1; j < array[0]; j++){
		if(array[j] > -1){
		    new_array[new_array[0]] = array[j];
		    new_array[0]++;
		}	
	    }
	    
	    sparse[i] = new_array;
	}
	
	return sparse;
    }
    
    public static int[][][] trimTsets(int all_tsets[][][], int lens[]){
	int temp[][][], tsets[][][];
	
	tsets = all_tsets;
	temp = tsets;
	
	tsets = new int[temp.length][][];
	
	for(int i = 0; i < temp.length; i++){
	    tsets[i] = new int[lens[i]][];
	    for(int j = 0; j < lens[i]; j++){
		tsets[i][j] = temp[i][j];		
	    }   
	}
	
	return tsets;
    }    
    
    public static int[] trimOneTset(int one_tset[], int length){
	int temp[], tset[];
	
	tset = one_tset;
	temp = tset;
	
	tset = new int[length];
	for(int i = 0; i < length; i++){
	    tset[i] = temp[i];
	}
	
	return tset;
    }  
    
    public static int[][] trimRegions(int rgns[][], int regions_size){
	int temp[][], regions[][];
	regions = rgns;
	temp = regions;
	
	regions = new int[regions_size][];
	for(int i = 0; i < regions_size; i++){
	    regions[i] = temp[i];
	}
	return regions;
    }
}
