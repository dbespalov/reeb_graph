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


import java.util.*;
import java.io.*;
import java.util.Vector.*;
import java.lang.Math.*;
import java.lang.Object.*;
import java.lang.System.*;


//Computes attributes for the finest resolution Reeb graph

public class AttributeCalculation {
    

    public double whole_area; //the area that is taken up by the object
    
    public Point points[]; //all points of the object
    public int sparse[][]; //sparse matrix, hold connectivity information of all points
    public double mu_values[]; //vector of values of mu function for each point
    public int points_length;
    
    public int MRG[][][];  //Multiresolutional reeb-graph
    public ReebGraphElement reebs[][];
    public int all_Tsets[][]; //vector of all Tsets of the object. It's a vector of vectors. Each vector in it holds
                             //indecies of the points (in Vector points) as Integer objects
    
    public AttributeElement attributes[];//Vector of attributes. Holds attributes for all nodes (tset) that 
                                         //are in the finest resolution reeb graph
        
    public double the_area; //area that is taken by all nodes in finest resolution
    
    public int count0;
    public int count1;
    public int count2;
    public int count3;
    
    public int INFINITY = 1000000000;
    
    
    // NOT USED!
    // public int pred[][]; //parents in BFS tree
    // public int c_child[][]; //common children in BFS tree
    // public int dist[];


    //main method of the procedure
    public AttributeElement[] doProcess(Point pts[], int pts_size, int sprs[][], double mu[], int Tsets[][],
			  int MR_graph[][][], ReebGraphElement rbs[][], double S){
	count0 = 0;
	count1 = 0;
	count2 = 0;
	count3 = 0;
	the_area = 0;
	
	points = pts;
	all_Tsets = Tsets;
	MRG = MR_graph;
	mu_values = mu;
	whole_area = S;
	sparse = sprs;
	reebs = rbs;
	points_length = pts_size;

	attributes = new AttributeElement[reebs[0].length];
	
	//System.out.println("Constructing BFS tree...");
	//constructBFS();
	
	//System.out.println("Calculating atributes...");
	calculateAttributes();
	
	// System.out.println("Area of the object: " + whole_area);
	// System.out.println("Area that is taken by the nodes in finest resolution: " + the_area);
	
	// System.out.print("Number of Tsets that consist out of 0 point: ");
	// System.out.println(count0);	
	// System.out.print("Number of Tsets that consist out of 1 point: ");
	// System.out.println(count1);	
	// System.out.print("Number of Tsets that consist out of 2 point: ");
	// System.out.println(count2);	
	// System.out.print("Number of times corrections were made for area calculation: ");
	// System.out.println(count3);	
	
	normalizeAttributes();
	
	return attributes;
    }





    
    //multiplies all the values of function a(n) by a coefficient
    //that sets the value SIM(M, M) = 1
    public void normalizeAttributes(){
	double coefficient = whole_area / the_area;
	for(int i = 0; i < attributes.length; i++){
	    attributes[i].a = attributes[i].a * coefficient;
	}	
    }
    
    //Goes through each Tsets and calculates attributes l and a for each one
    public void calculateAttributes(){
	double rnum = (double) MRG.length;
		
	double lens[] = new double[all_Tsets.length];
	double sum_len = 0;
	double len;
	
	for(int i = 0; i < all_Tsets.length; i++){
	    int Tset[] = all_Tsets[i];
	    
	    //computes value of a(m) for each tset
	    AttributeElement temp = new AttributeElement();
	    temp.a = (1 / rnum) * (calculateTsetArea(Tset) / whole_area);
	    
	    attributes[i] = temp;
	    
	    int ptr_ind = Tset[0];
	    double temp_double;
	    double min_m = mu_values[ptr_ind];
	    double max_m = mu_values[ptr_ind];
	 	    
	    //looks for min and maximum values of mu among all points that lie in one tset
	    for(int j = 1; j < Tset.length; j++){
		ptr_ind = Tset[j];
		temp_double = mu_values[Tset[j]];
		
		if(min_m > temp_double)
		    min_m = temp_double;
		
		if(max_m < temp_double)
		    max_m = temp_double;
	    }
	    
	    len = max_m - min_m;
	    sum_len = sum_len + len;
	    
	    //creates a vector of the len values for future calculation of the value of l(m) function
	    lens[i] = len;
	}
	
	//calculates l(m) for each tset
	for(int i = 0; i < all_Tsets.length; i++){
	    attributes[i].l = (1 / rnum) * (lens[i] / sum_len);
	}
    }

        


    //computes the area that is taken by one tset
    public double calculateTsetArea(int points_in_area[]){
	
	int reference[] = new int[points_length];
	for(int i = 0; i < points_length; i++)
	    reference[i] = 0;
	for(int i = 0; i < points_in_area.length; i++)
	    reference[points_in_area[i]] = 1;
	
	int indexA, indexB, indexC;
	Point A, B, C;
	double area = 0;
	
	for(int i = 0; i < points_in_area.length - 2; i++){
	    indexA = points_in_area[i];
	    int adj[] = sparse[indexA];
	    int len = adj[0];
	    
	    for(int j = 1; j < len; j++){
		indexB = adj[j];
		
		if(indexA > indexB && reference[indexB] == 1){
		    
		    int adj2[] = sparse[indexB];
		    int len2 = adj2[0];
		    for(int k = 1; k < len2; k++){
			indexC = adj2[k];
			
			if(indexB > indexC &&
			   reference[indexC] == 1 &&
			   isConnected(indexA, indexC)==true){
			    
			    area = area + calculateTrigArea(points[indexA], points[indexB], points[indexC]);
			}
		    }
		}	    
	    }	    
	}
	
	the_area = the_area + area;
	
	if(points_in_area.length == 0){
	    count0++;
	    return 0;
	}
	
	if(points_in_area.length == 1){
	    count1++;
	    return 0;
	}
	
	if(points_in_area.length == 2){
	    count2++;
	    return 0;
	}
	
	if(area == 0.0)
	    count3++;	
	
	return area;
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
        
    //computes area of a triangle which is defined by 3 points A, B, C
    public double calculateTrigArea(Point A, Point B, Point C){
	double a, b, c, p, area;
	
	a = Math.sqrt((C.X - B.X)*(C.X - B.X)
		      + (C.Y - B.Y)*(C.Y - B.Y)
		      + (C.Z - B.Z)*(C.Z - B.Z));

	b = Math.sqrt((A.X - C.X)*(A.X - C.X)
		      + (A.Y - C.Y)*(A.Y - C.Y)
		      + (A.Z - C.Z)*(A.Z - C.Z));

	c = Math.sqrt((A.X - B.X)*(A.X - B.X)
		      + (A.Y - B.Y)*(A.Y - B.Y)
		      + (A.Z - B.Z)*(A.Z - B.Z));

	p = (a + b + c) / 2;
       
	area = (p * (p - a) * (p - b) * (p - c));
	
	if(area < 0){
	    return 0;
	}
	
	area = Math.sqrt(area);

	return area;
    }

    //calculate distance between two points A and B
    public double calculateDistance(Point A, Point B){
	return Math.sqrt(Math.pow(A.X - B.X, 2)
			 + Math.pow(A.Y - B.Y, 2)
			 + Math.pow(A.Z - B.Z, 2));
    }
    


        //creates lists of cross-edges, parents and children using BFS algorithm
    // public void constructBFS(){
	
    // 	int s = 0; //source node
	
    // 	int color[] = new int[points_length]; //1 - white (not visited), 
    // 	                                      //2 - gray (in the process), 
    // 	                                      //3 - black (finished)
    // 	int queue[] = new int[aman.size_inc_queue];
    // 	queue[0] = 2;
    // 	queue[1] = 2;
	
	
    // 	int d[] = new int[points_length]; //holds distance (level) info
    // 	pred = new int[points_length][aman.size_inc_sparse];
    // 	c_child = new int[points_length][aman.size_inc_sparse];
    // 	//c_edges = new int[aman.size_inc_edges][2];
	
    // 	for(int i = 0; i < points_length; i++){
    // 	    color[i] = 1;
    // 	    d[i] = INFINITY;
    // 	    pred[i][0] = 1;
    // 	    c_child[i][0] = 1;
    // 	}
	
    // 	color[s] = 2;
    // 	d[s] = 0;
	
    // 	queue[queue[1]] = s;
    // 	queue[1]++;
	
    // 	int u, v;
	
    // 	while(queue[0] != queue[1]){
    // 	    u = queue[queue[0]];
    // 	    queue[0]++;
	    
    // 	    int adj[] = sparse[u];
    // 	    int len = adj[0];
	    
    // 	    for(int i = 1; i < len; i++){
    // 		v = adj[i];
    // 		if(color[v] == 1){
    // 		    color[v] = 2;
    // 		    d[v] = d[u] + 1;
		    
    // 		    int temp[] = pred[v];
    // 		    int temp_l = pred[v][0];
    // 		    if(temp_l >= temp.length)
    // 			pred[v] = aman.expand_sparse(temp, temp_l);
		    
    // 		    pred[v][temp_l] = u;
    // 		    pred[v][0]++;
		    
    // 		    if(queue[1] >= queue.length)
    // 			queue = aman.expand_queue(queue);
		    
    // 		    queue[queue[1]] = v;
    // 		    queue[1]++;
    // 		}
    // 	    }
	    
    // 	    color[u] = 3;
    // 	}
	
    // 	for(int i = 0; i < points_length; i++){
	    
    // 	    if(pred[i][0] != 1){
    // 		int pars[] = pred[i];
    // 		int len = pred[i][0];
		
    // 		for(int j = 1; j < len; j++){
    // 		    int index = pars[j];
		    
    // 		    int temp[] = c_child[index];
    // 		    int temp_l = c_child[index][0];
		    
    // 		    if(temp_l >= temp.length)
    // 			c_child[index] = aman.expand_sparse(temp, temp_l);
		    
    // 		    c_child[index][temp_l] = i;
    // 		    c_child[index][0]++;		    
    // 		}		
    // 	    }	    
    // 	}
	
    // 	dist = d;
    // }


}
