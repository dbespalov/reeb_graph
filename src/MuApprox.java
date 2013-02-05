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

//Computes approximate values for mu function

import java.util.*;
import java.io.*;
import java.util.Vector.*;
import java.lang.Math.*;
import java.lang.Object.*;
import java.lang.System.*;

public class MuApprox {

    public static double coefficient;
   
    public double r;
    public double S;
    public static double infinity = 100000000000000.0;

    public int sparse[][];
    public Point points[];
    public int points_length;

    public GValueStorage g_values[][];
    public int base_points_length;
    
    //public Vector base_points;
    public Point base_points[];
    public double base_areas[];
        
    public double mu_values[];

    public Point unvisited_vertices[];

    public double the_area;

    public int found;
    public int not_found;

    public int base_counter;
    public int count0, count1, count2;
    
    //PrintWriter out;
    
    // main method for the procedure
    public double[] doProcess(double coeff, double s, Point pts[], int points_size, int sprs[][]) throws IOException{
	
	/*out = new PrintWriter(new FileWriter("sparse"));
	sparse = sprs;
	points = pts;
	for(int i = 0; i < points_size; i++){
	    int adj[] = sparse[i];
	    
	    out.print(i + " [" + adj[0] + "]: ");
	    for(int j = 1; j < adj[0]; j++){
		out.print(adj[j] + " ");
	    }
	    out.println("");
	}
	
	for(int i = 0; i < points_size; i++){
	    out.print(i + ": ");
	    if(points[i] == null)
		out.println("NULL");
	    else
		out.println("(" + points[i].X + ", " + points[i].Y + ", " + points[i].Z + ")");
	}
	
	out.close();
	*/


	coefficient = coeff;
	S = s;
	the_area = 0;
	
	//System.out.println("coefficient = " + coefficient);
	
	found = 0;
	not_found = 0;
	int counter = 0;
	base_counter = 5;
	count0 = count1 = count2 = 0;
	
	sparse = sprs;
	points = pts;
	points_length = points_size;
	
	//System.out.println("Calculating threshold...");
	calculateThreshold();
	System.out.println("Threshold r = " + r);
	
	//System.out.println("Calculating G approximations...");
	calculatePaths();
	
	System.out.print("Total number of base points: ");
	System.out.println(base_points_length);
		

	//System.out.println("Calculating MU approximations...");
	for(int i = 0; i < points_length; i++){
	    mu_values[i] = mu(i);
	}
	
	//System.out.println("Size of g_values: " + base_points_length);
	
	//System.out.println("Found: " + found);
	//System.out.println("Not Found: " + not_found);
		
	the_area = 0;
	for(int i = 0; i < base_points_length; i++){
	    the_area = the_area + base_areas[i];
	    if(base_areas[i] < 0)
		System.out.println("WARNING in MuApprox: base_areas area is negative for i=" + i);
	}
    
	//System.out.println("Number of areas with 0 points: " + count0);
	//System.out.println("Number of areas with 1 points: " + count1);	
	//System.out.println("Number of areas with 2 points: " + count2);
	//System.out.print("The real area that is taken by base vertices is: ");
	//System.out.println(the_area - count1 - (count2*2));
	
	return mu_values;
    }
    
    public void calculatePaths() throws IOException{
	
	//create VLIST and fill unvisited_vertices
	//MinHeap VLIST = new MinHeap();
	MinHeapElement temp_el;
	Point temp_point;
	unvisited_vertices = new Point[points_length];
	
	base_points_length = 0;
	base_points = new Point[aman.size_inc_base_points];
	base_areas = new double[aman.size_inc_base_points];
	
	//fills up the unvisited_verteces with all vertices
	for(int i = 0; i < points_length; i++){
	    unvisited_vertices[i] = points[i];
	}
	
	g_values = new GValueStorage[aman.size_inc_base_points][points_length];
	mu_values = new double[points_length];
	
	int last_index = 0;
	int counter = 0;
	boolean stopper = false;
	while(stopper == false){
	    	    
	    //takes one vertex from unvisited_vertices and puts it into
	    //base_points. Then calls calculateShortestPath
	    int j = last_index;
	    Point temp = null;
	    while(j < unvisited_vertices.length){
		
		if(unvisited_vertices[j] != null){
		    temp = unvisited_vertices[j];
		    last_index = j;
		    break;
		}
		j++;
	    }
	    
	    if(temp != null){
		int temp_int = j;
		
		if(base_points_length >= base_points.length){
		    base_points = aman.expand_base(base_points, base_points_length);
		    base_areas = aman.expand_base(base_areas, base_points_length);
		    g_values = aman.expand_base(g_values, base_points_length);
		}
	
		base_points[base_points_length] = temp;
		
		MinHeapElement arr[] = new MinHeapElement[points_length+1];
		for(int i = 1; i <= points_length; i++){
		    temp_el = new MinHeapElement();
		    temp_el.key = infinity;
		    temp_el.index = i-1;
		    arr[i] = temp_el;
		}
		MinHeap2 VLIST = new MinHeap2(arr, points_length);
		
		g_values[base_points_length] = calculateShortestPath(temp_int, VLIST);
		base_points_length++;
	    }
	    else if(temp == null)
		stopper = true;
	}
    }
    
    //returns vector that holds the path from the base to every single point
    public GValueStorage[] calculateShortestPath(int base_vertex_index, MinHeap2 VLIST){
	GValueStorage result[] = new GValueStorage[points_length];
			
	//fills up result vector with the infinity values
	for(int z = 0; z < result.length; z++){
	    result[z] = new GValueStorage();
	    result[z].index = z;
	    result[z].value = infinity;
	}
	
	MinHeapElement start_el = new MinHeapElement();
	MinHeapElement smallest;
	GValueStorage temp_storage;
		
	start_el.key = 0.0;
	start_el.index = base_vertex_index;
	
	VLIST.decreaseKey(base_vertex_index, 0.0);
		
	//int counter = VLIST.count - 1;
	double length_VVa;
	double g_V, g_Va;
	while(VLIST.heap_size != 0){
	    smallest = VLIST.extractMin();
	    g_V = smallest.key;
	    
	    //saves the values of g for each vertex removed from VLIST
	    if(result[smallest.index].value > smallest.key){
		result[smallest.index].value = smallest.key;
	    }
	    	    
	    //for each vertex from adj_vers checks g(Va) > g(V) + length(V,Va)
	    int adj[] = sparse[smallest.index];
	    int len = adj[0];
	    for(int j = 1; j < len; j++){
		int index = adj[j];
		
		length_VVa = calculateDistance(points[smallest.index], points[index]);
		
		g_Va = result[index].value;
		
		//Main Check
		if(g_Va > g_V + length_VVa){
		    result[index].value = g_V + length_VVa;
		    
		    //out.println("the point "+ j + " has been updated!");
		    
		    VLIST.decreaseKey(index, g_V + length_VVa);
		}
	    }
	}
		
	Point base_ver = points[base_vertex_index];
	

	//creating the vector of vertices that are in the area	
	int points_in_area_length = 0;
	int points_in_area[] = new int[aman.size_inc];

	for(int i = 0; i < result.length; i++){
	    temp_storage = result[i];
	    //checks if g(V) <= r (i.e., V point is in the area)
	    if(temp_storage.value <= r){
		
		//checks if V was not included in the area
		//then adds it to the area vector
		if(unvisited_vertices[temp_storage.index] != null
		   && temp_storage.index != base_vertex_index){
		    if(points_in_area_length >= points_in_area.length)
			points_in_area = aman.expand(points_in_area, points_in_area_length);
		    
		    points_in_area[points_in_area_length] = temp_storage.index;
		    points_in_area_length++;
		}
		unvisited_vertices[temp_storage.index] = null;		    
	    }
	}
	
	if(points_in_area_length >= points_in_area.length)
	    points_in_area = aman.expand(points_in_area, points_in_area_length);
	
	points_in_area[points_in_area_length] = base_vertex_index;
	points_in_area_length++;
	
	unvisited_vertices[base_vertex_index] = null;
	
	//calculates the area around the base point
	calculateBaseArea(points_in_area, points_in_area_length);
	
    	return result;	
    }
    
    public double mu(int point_index){
	if(base_points == null){
	    System.out.println("ERROR in MuApprox: base_points are not initialized!\n Exiting...");
	    System.exit(1);
	}
	
	double value = 0;
	Point base_point;

	for(int i = 0; i < base_points_length; i++){
	    value = value + (g(point_index, i) * base_areas[i]);
	}
	return value;
    }
    
    public double g(int point_index, int base_index){
	GValueStorage temp = g_values[base_index][point_index];

	if(temp.value != infinity){
	    found++;
	    return temp.value;
	}
	else {
	    not_found++;
	    return 0;
	}
    }
        
    public void calculateBaseArea(int points_in_area[], int points_in_area_length){
	int indexA, indexB, indexC;
	Point A, B, C;
	double area = 0;
	
	for(int i = 0; i < points_in_area_length; i++){
	    indexA = points_in_area[i];
	  
	    int adj1[] = sparse[indexA];
	    int len1 = adj1[0];
	    for(int j = 1; j < len1; j++){
		indexB = adj1[j];
		
		if(indexA > indexB){
		    
		    int adj2[] = sparse[indexB];
		    int len2 = adj2[0];
		    for(int k = 1; k < len2; k++){
			indexC = adj2[k];
			
			if(indexB > indexC && isConnected(indexC, indexA) == true){
			    area = area
				+ calculateTrigArea(points[indexA], points[indexB], points[indexC]);
			}
		    }
		}	    
	    }	    
	}
	
	//	the_area = the_area + area;
	
	if(points_in_area_length == 0){
	    area = 0;
	    count0++;
	}
	
	if(points_in_area_length == 1){
	    area = 1;
	    count1++;
	}
	
	if(points_in_area_length == 2){
	    area = 2;
	    count2++;
	}
	
	the_area = the_area + area;
	
	base_areas[base_points_length] = area;
    }
    
    public double calculateTrigArea(Triangle T){
	double a, b, c, p, area;
	
	a = Math.sqrt(Math.pow(points[T.c].X - points[T.b].X, 2)
		      + Math.pow(points[T.c].Y - points[T.b].Y, 2)
		      + Math.pow(points[T.c].Z - points[T.b].Z, 2));
	
	b = Math.sqrt(Math.pow(points[T.a].X - points[T.c].X, 2)
		      + Math.pow(points[T.a].Y - points[T.c].Y, 2)
		      + Math.pow(points[T.a].Z - points[T.c].Z, 2));
	
	c = Math.sqrt(Math.pow(points[T.a].X - points[T.b].X, 2)
		      + Math.pow(points[T.a].Y - points[T.b].Y, 2)
		      + Math.pow(points[T.a].Z - points[T.b].Z, 2));
	
	p = (a + b + c) / 2;
	
	area = p * (p - a) * (p - b) * (p - c);
	
	if(area < 0)
	    area = 1.0;

	area = Math.sqrt(area);
	
	return area;
    }

    public double calculateTrigArea(Point A, Point B, Point C){
	double a, b, c, p, area;

	a = Math.sqrt(Math.pow(C.X - B.X, 2)
		      + Math.pow(C.Y - B.Y, 2)
		      + Math.pow(C.Z - B.Z, 2));

	b = Math.sqrt(Math.pow(A.X - C.X, 2)
		      + Math.pow(A.Y - C.Y, 2)
		      + Math.pow(A.Z - C.Z, 2));

	c = Math.sqrt(Math.pow(A.X - B.X, 2)
		      + Math.pow(A.Y - B.Y, 2)
		      + Math.pow(A.Z - B.Z, 2));

	p = (a + b + c) / 2;

	area = (p * (p - a) * (p - b) * (p - c));
	
	if(area < 0)
	    area = 1.0;
	
	area = Math.sqrt(area);

	return area;
    }

    public void calculateThreshold(){
	Triangle T;
	
	r = Math.sqrt(coefficient * S);
    }

    public double calculateDistance(Point A, Point B){
	return Math.sqrt(((A.X - B.X)*(A.X - B.X))
			 + ((A.Y - B.Y)*(A.Y - B.Y))
			 + ((A.Z - B.Z)*(A.Z - B.Z)));
    }
    
    public boolean isConnected(int index1, int index2){
	int adj[] = sparse[index1];
	int len = adj[0];
	
	for(int i = 1; i < len; i++){
	    if(adj[i] == index2){
		return true;	
	    }
	}
	
	return false;
    }
    

}


