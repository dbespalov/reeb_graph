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

// We use adjacency list representation for 3D models (triangle meshes)
// (the adjacency list is stored as a sparse matrix)

import java.util.*;
import java.io.*;
import java.util.Vector.*;
import java.lang.Math.*;
import java.lang.Object.*;

public class SparseMatrix {

    public int sparse[][];

    //Main function of the class.
    public Vector createMatrix(Triangle triangles[], int triangles_length,  Point points[], int points_length){
	
	triangles = makeRandom(triangles, triangles_length);
	
	Triangle temp;
	int temp_int = 0;
		
	//System.out.println("The size of triangles: " + triangles_length);
	//System.out.println("The size of points: " + points_length);
	
	sparse = new int[points.length][aman.size_inc_sparse+1];
	for(int i = 0; i < points_length; i++){
	    sparse[i][0] = 1;
	}
	
	//add each triangle with addToMatrix 
	for(int i = 0; i < triangles_length; i++){
	    int a,b,c;
	    a = triangles[i].a;
	    b = triangles[i].b;
	    c = triangles[i].c;
	    
	    if(a != b && isConnected(a, b) == false)
		addConnection(a,b);
	    
	    if(b != c && isConnected(b, c) == false)
		addConnection(b,c);
			    
	    if(a != c && isConnected(a, c) == false)
		addConnection(a,c);
	}
		
	//checkAlgorithm(points_length);
	
	Vector result = new Vector();
	result.addElement(points);
	result.addElement(new Integer(points_length));
	result.addElement(sparse);
	return result;
    }
    
    public void checkAlgorithm(int points_length){
	
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
		    System.out.println("WARNING in SparseMatrix: can not find point i=" + i);	    
		}
	    }
	}
    }
    
    public void addConnection(int index1, int index2){
	int length = sparse[index1][0];
	if(sparse[index1].length <= length)
	    sparse[index1] = aman.expand_sparse(sparse[index1], length);
	
	sparse[index1][length] = index2;
	sparse[index1][0]++;
	
	length = sparse[index2][0];
	if(sparse[index2].length <= length)
	    sparse[index2] = aman.expand_sparse(sparse[index2], length);
	
	sparse[index2][length] = index1;
	sparse[index2][0]++;	
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
    
    public Triangle[] makeRandom(Triangle triangles[], int length) {
	Triangle temp;
	for(int i = 0; i < length; i++){
	    
	    double random1 = Math.random() * (length - 1);
	    double random2 = Math.random() * (length - 1);
	    
	    int n1 = (int) random1;
	    int n2 = (int) random2;
	    
	    temp = triangles[n1];
	    triangles[n1] = triangles[n2];
	    triangles[n2] = temp;	    
	}
			
	return triangles;
    }


}
