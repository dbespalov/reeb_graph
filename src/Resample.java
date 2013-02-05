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

// Resamples the verticies before the construction of MRG until
// the specified number of points (points_number) is reached

import java.util.*;
import java.io.*;
import java.util.Vector.*;
import java.lang.Math.*;


public class Resample {

    public Point points[];
    public int sparse[][];
    
    public double THRESHOLD; //threshold for resampling vertices
    
    public int points_length;

    public Vector doProcess(int points_number, Point pts[], int points_l, int sprs[][]){
	
	points = pts;
	points_length = points_l;
	sparse = sprs;
	
	THRESHOLD = approximateEdgeLength() * 2;
	
	//System.out.println("Number of points before: " + points_length);
	while(points_length < points_number - (points_number * .0333)){
	    resampleOnce(points_number);
	    THRESHOLD = THRESHOLD - (THRESHOLD * 0.05);
	    aman.trimSparse(sparse, points_length);
	}
	
	//System.out.println("Number of points now: " + points_length);
	//System.out.println("Size of matrix now: " + matrix.length + "X" + matrix[0].length);
	
	checkAlgorithm();
	
	int all_points[] = new int[points_length];
	
	for(int i = 0; i < points_length; i++){
	    all_points[i] = i;
	}
	
	Vector result = new Vector();
	result.addElement(points);
	result.addElement(new Integer(points_length));
	result.addElement(sparse);
	return result;
    }

    public void resampleOnce(int points_number){
	boolean stopper = false;
	int i = 0;
	while(stopper == false){
	    int adj[] = sparse[i];
	    for(int j = 1; j < adj[0]; j++){
		int index = adj[j];
		if(index > -1 && distance(i, index) > THRESHOLD){
		    Point A,B,C;
		    		    
		    A = points[i];
		    B = points[index];
		    		    
		    if(points_length >= points.length){
			points = aman.expand(points, points_length);
			sparse = aman.expand_sparse(sparse, points_length);
		    }

		    //create new point
		    Point P = new Point();
		    
		    P.X = (A.X + B.X) / 2.0;
		    P.Y = (A.Y + B.Y) / 2.0;
		    P.Z = (A.Z + B.Z) / 2.0;   
		    
		    points[points_length] = P;
		    points_length++;
		    
		    //removes old connections
		    sparse[i][j] = -1;
		    
		    int k = 1;
		    int len = sparse[index][0];
		    boolean stop = false;
		    while(stop == false){
			if(sparse[index][k] == i){
			    sparse[index][k] = -1;
			    stop = true;
			}
			k++;
		    }
		    
		    //adds new ones
		    len = sparse[i][0];
		    int adj3[] = sparse[i];
		    if(adj3.length <= len)
			sparse[i] = aman.expand_sparse(adj3, len);
		    
		    sparse[i][len] = points_length-1;
		    sparse[i][0]++;
		    
		    len = sparse[index][0];
		    adj3 = sparse[index];
		    if(adj3.length <= len)
			sparse[index] = aman.expand_sparse(adj3, len);
		    
		    sparse[index][len] = points_length-1;
		    sparse[index][0]++;
		    
		    //create new element in sparse matrix
		    sparse[points_length-1] = new int[aman.size_inc_sparse + 1];
		    sparse[points_length-1][0] = 3;
		    sparse[points_length-1][1] = i;
		    sparse[points_length-1][2] = index;
		    
		    //look for extra connections
		    adj = sparse[i];
		    int adj2[] = sparse[index];
		    for(k = 1; k < adj[0]; k++){
			int i1 = adj[k];
			if(i1 > -1){
			    for(int h = 1; h < adj2[0]; h++){
				int i2 = adj2[h];
				if(i2 > -1 && i1 == i2){
				    
				    len = sparse[i1][0];
				    adj3 = sparse[i1];
				    if(adj3.length <= len)
					sparse[i1] = aman.expand_sparse(adj3, len);
				    
				    sparse[i1][len] = points_length-1;
				    sparse[i1][0]++;
				    
				    len = sparse[points_length - 1][0];
				    adj3 = sparse[points_length-1];
				    if(adj3.length <= len)
					sparse[points_length-1] = aman.expand_sparse(adj3, len);
				    
				    sparse[points_length-1][len] = i1;
				    sparse[points_length-1][0]++;
				}
			    }			    
			}			
		    }
		}		
	    }	    
	    
	    i++;
	    
	    if(i >= points_length || points_length >= points_number * 2){
		stopper = true;		
	    }	    
	}	
    }

    public double approximateEdgeLength(){
	
	double AB;
	double max = 0;
	
	for(int i = 0; i < points_length; i++){
	    int adj[] = sparse[i];
	    int len = adj[0];
	    
	    for(int j = 1; j < len; j++){
		if(i > adj[j]){
				    
		    AB = distance(i, adj[j]);
		    
		    if(max < AB)
			max = AB;
		}
	    }
	}
			
	return max;
    }

    public double distance(int a, int b){
	Point A,B;
	A = points[a];
	B = points[b];
	
	double dist = Math.sqrt(Math.pow((A.X - B.X), 2)
				+ Math.pow((A.Y - B.Y), 2)
				+ Math.pow((A.Z - B.Z), 2));
	return dist;
    }
    
    public int placeInPoints(Point P){
	Point P1;
	
	for(int i = 0; i < points_length; i++){
	    P1 = points[i];
	    
	    if(ecomp.eq(P.X,P1.X) && 
	       ecomp.eq(P.Y,P1.Y) &&
	       ecomp.eq(P.Z,P1.Z))
		return i;	    
	}
	
	return points_length;
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

    public void checkAlgorithm(){
	
	for(int i = 0; i < points_length; i++){
	    
	    for(int j = 1; j < sparse[i][0]; j++){
		
		int index = sparse[i][j];
		boolean found = false;
		for(int k = 1; k < sparse[index][0]; k++){
		    
		    if(sparse[index][k] == i){
			found = true;			
		    }	    
		}
		
		if(found == false){
		    System.out.println("WARNING in Resample: can not find point i = " + i);
		}
	    }
	}
    }
        

}

