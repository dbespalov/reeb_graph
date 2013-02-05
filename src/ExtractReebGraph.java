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
import java.util.Random.*;
import java.io.*;

 /*
 |  The method for estimating similarity of 3D mesh models based on 
 |  multiresolutional Reeb graphs (MRG) was proposed by Hilaga et al. [1].
 |
 |  ExtractReebGraph.java implements construction of MRGs 
 |  for 3D mesh models (see Section 4 in [1]) stored in (pseudo-) VRML format.
 |
 |  CompareReebGraph.java implements matching algorithm 
 |  for a pair of MRGs (see Section 5 in [1])
 |
 |
 |  References:
 |  [1] Topology Matching for Fully Automatic Similarity Estimation of 3D Shapes 
 |      Masaki Hilaga, Yoshihisa Shinagawa, Taku Kohmura and Tosiyasu L. Kunii       
 |      SIGGRAPH, 2001                                                              
 */

 /*
  | WARNING. VRML parser can only handle specially-formatted face-vertex meshes:
  |
  | 1. 3D points are assumed to be enclosed by strings "point [\n" and "]\n".
  |     Coordinates for each point is specified on a separate line in the file
  |         e.g.: string "1.5 3.2 0.2,\n"
  |
  |  2. Face lists are assumed to be enclosed by strings "coordIndex [\n" and "]\n"
  |     Each face must appear on a separate line
  |         e.g.: string "0, 2, 1, -1,\n" encodes a triangle face 
  |               (quad faces are also allowed)
  |
  | 3. Please refer to the sample VRML models for examples of required format
 */


public class ExtractReebGraph {

    public static int points_number;
    public static double mu_coeff;
    public static int MRG_number;
    
    public static Triangle triangles[];
    public static int triangles_length;

    public static int sparse[][];
    public static Point points[];
    public static int points_length;
    public static double mu_values[];
    
    public static int all_Tsets[][];
    public static int MRG[][][];
    public static ReebGraphElement reebs[][];
    public static AttributeElement attributes[];

    public static double whole_area;
    
    /*
    |  usage: java -Xmx512m  ExtractReebGraph  <num_pts>  <mu_coeff>  <mrg_size>  <model_1>.wrl  <model_2>.wrl ... <model_N>.wrl
    |
    |    <num_pts>       -- target number of vertices before MRG is extracted 
    |                       (triangle meshes in each 3D model are resampled to match <num_pts>)
    |
    |    <mu_coeff>      -- coefficient for calculating threshold parameter r=sqrt(mu_coeff * area(S)),
    |                       used to approximate mu values
    |
    |    <mrg_size>      -- number of ranges K used in the finest resolution of MRG
    |	   
    |    <model_i>.wrl   -- i-th VRML model to process, where i=[1,N]. 
    |                       Extracted MRG for this model is stored in a text file <model_i>.mrg 
    */

    public static void main (String arg[]) throws IOException {
		
	points_number = (new Integer(arg[0])).intValue();
	mu_coeff = (new Double(arg[1])).doubleValue();

	MRG_number = (new Integer(arg[2])).intValue();
	
	//System.out.println("points_number = " + points_number);
	//System.out.println("mu_coeff = " + mu_coeff);
	//System.out.println("MRG_number = " + MRG_number);

	//more than one VRML model to process can be specified on the command line
	for(int i = 3; i < arg.length; i++){
	    main_one(arg[i]);	    
	}
    }
      



    // computes MRG for a VRML model (arg is the filename)
    public static void main_one (String arg) throws IOException {
	
	System.out.println("Working on " + arg + " ...");
	String str, filename;
	filename = arg;
	
	

	boolean Stop=false,Stop2 = false;
	Point P,V1,V2,V3,SA,SB,SC,SD;
	int sa,sb,sc,sd;
	int a,b,c,d,m,n,o,p,q, prevsize = 0;
	double Angle1,Angle2,Angle3,length1,length2,length3;
	Triangle T;
	Point Points[] = new Point[aman.size_inc];
	points_length = 0;
	Point Ptemp[];
	Triangle Triangles[] = new Triangle[aman.size_inc];
	triangles_length = 0;
	Triangle Ttemp[];
	
	//get input wrl file
	BufferedReader in = new BufferedReader(new FileReader(filename));
	
	/*
	| WARNING. VRML parser can only handle specially-formatted face-vertex meshes:
	|
	| 1. 3D points are assumed to be enclosed by strings "point [\n" and "]\n".
	|     Coordinates for each point is specified on a separate line in the file
	|         e.g.: string "1.5 3.2 0.2,\n"
	|
	|  2. Face lists are assumed to be enclosed by strings "coordIndex [\n" and "]\n"
	|     Each face must appear on a separate line
	|         e.g.: string "0, 2, 1, -1,\n" encodes a triangle face 
        |               (quad faces are also allowed)
	|
	| 3. Please refer to the sample VRML models for examples of required format
	*/

	//*******************************************************************************
	// VRML parser BEGIN 
	//*******************************************************************************

	//read in point and triangle coordinates from vrml file
	try {
	    while (true) {
	        Stop = false;
	        str = in.readLine();
	        if (str.indexOf("point") != -1 && str.indexOf('[') != -1
		    && str.lastIndexOf(']') < str.lastIndexOf('[')  ) {
		    prevsize = points_length;
		    Stop2 = false;
	            while (!Stop) {
                        str = in.readLine();
			str = str.trim();
		        int temp_int = str.indexOf(']');
			
			if (str.indexOf("#IGNORE") != -1 || 
                            //(str.indexOf(']') != -1 && temp_int == -1)) {
			    temp_int == 0){
		            Stop = true;
			}
	                else if (Stop == false) {
	                    m = str.indexOf(" ");
	                    n = str.indexOf(" ", m+1);
	                    o = str.indexOf(',');
			    			    
			    int temp_int2 = str.indexOf(']');
			    int temp_int3 = str.indexOf(" ", n+1);
			    if(o == -1){
				if(temp_int2 == -1){
				    o = str.length();				    
				}
				else{
				    if(temp_int3 == -1)
					o = temp_int2;
				    else
					o = temp_int3;
				}
			    }
			    
		            if (m != -1) {
				
				if(points_length >= Points.length)
				    Points = aman.expand(Points, points_length);
				
				P = new Point();
                                P.X = Double.valueOf(str.substring(0,m))
                                                     .doubleValue();
                                P.Y = Double.valueOf(str.substring(m,n))
                                                     .doubleValue();
                                P.Z = Double.valueOf(str.substring(n,o))
                                                     .doubleValue();
	                        Points[points_length] = P;
				points_length++;
		            }
			    
			    if(str.indexOf(']') != -1)
				Stop = true;			    
			}
		    }
	        }
		else if (str.indexOf("coordIndex") != -1 && 
                         str.indexOf('[') != -1 && 
			 str.lastIndexOf(']') < str.lastIndexOf('[') &&
			 !Stop2) {

		    //System.out.println("size of Points is " + Points.length);
		    Stop2 = true;
	            while (!Stop) {
                        str = in.readLine();
			str = str.trim();
			
			int temp_int = str.indexOf(']');
						
			if(temp_int == 0){
		            Stop = true;
			}
			else {
		            q = str.indexOf('-');
			    if(q == -1){
				m = -1;				
				n = -1;
				o = -1;
				p = -1;
			    }
			    else{
				String temp = str.substring(0, q);
				
				m = temp.indexOf(',', 0);
				n = temp.indexOf(',', m+1);
				o = temp.indexOf(',', n+1);
				p = temp.indexOf(',', o+1);
			    }
                            			    			    
		            if (str.indexOf(']') != -1){
				Stop = true;
			    }
			    
			    if (m != -1 && p == -1) {
				
				a = Integer.parseInt(str.substring(0,m));
       	                        b = Integer.parseInt(str.substring(m+2,n));
       	                        c = Integer.parseInt(str.substring(n+2,o));
				T = new Triangle();
		                T.a = prevsize+a;
		                T.b = prevsize+b;
		                T.c = prevsize+c;
				
				if(triangles_length >= Triangles.length)
				    Triangles = aman.expand(Triangles, triangles_length);
				
				Triangles[triangles_length] = T;
				triangles_length++;
			    }
			    else if (m != -1 && p != -1) {
       	                        a = Integer.parseInt(str.substring(0,m));
       	                        b = Integer.parseInt(str.substring(m+2,n));
       	                        c = Integer.parseInt(str.substring(n+2,o));
			        d = Integer.parseInt(str.substring(o+2,p));
				
				SA = Points[prevsize+a];
				SB = Points[prevsize+b];
				SC = Points[prevsize+c];
				SD = Points[prevsize+d];

				V1 = new Point();
				V1.X = SB.X-SA.X;
				V1.Y = SB.Y-SA.Y;
				V1.Z = SB.Z-SA.Z;
				V2 = new Point();
				V2.X = SC.X-SA.X;
				V2.Y = SC.Y-SA.Y;
				V2.Z = SC.Z-SA.Z;
				V3 = new Point();
				V3.X = SD.X-SA.X;
				V3.Y = SD.Y-SA.Y;
				V3.Z = SD.Z-SA.Z;

				length1 = Math.sqrt
                                          (V1.X*V1.X+V1.Y*V1.Y+V1.Z*V1.Z);
				length2 = Math.sqrt
				          (V2.X*V2.X+V2.Y*V2.Y+V2.Z*V2.Z);
				length3 = Math.sqrt
				          (V3.X*V3.X+V3.Y*V3.Y+V3.Z*V3.Z);

				Angle1 = Math.acos((V1.X*V2.X+V1.Y*V2.Y+V1.Z*V2.Z)/(length1*length2));
				Angle2 = Math.acos((V1.X*V3.X+V1.Y*V3.Y+V1.Z*V3.Z)/(length1*length3));
				Angle3 = Math.acos((V2.X*V3.X+V2.Y*V3.Y+V2.Z*V3.Z)/(length2*length3));
				if (Angle1 > Angle2) {
				    if (Angle3 > Angle1) {
					
					sa = prevsize+c;
					sc = prevsize+d;
					sb = prevsize+b;
				    }
				    else {
					sa = prevsize+b;
					sc = prevsize+c;
					sb = prevsize+d;
				    }
				}
				else {
				    if (Angle2 > Angle3) {
					sa = prevsize+b;
					sc = prevsize+d;
					sb = prevsize+c;
				    }
				    else {
					sa = prevsize+c;
					sc = prevsize+d;
					sb = prevsize+b;
				    }
				}
				sd = prevsize+a;

				if(triangles_length >= Triangles.length-1)
				    Triangles = aman.expand(Triangles, triangles_length);
								
				T = new Triangle();

		                T.a = sa;
		                T.b = sb;
		                T.c = sc;
				
				Triangles[triangles_length] = T;
				triangles_length++;

				T = new Triangle();
				T.a = sa;
		                T.b = sd;
		                T.c = sc;
				
				Triangles[triangles_length] = T;
				triangles_length++;
			    }
			}
		    }
	        }
	    }
	}
	
	catch(NullPointerException e) {}
	
	//NOTE that the actual sizes of points and triangles are stored in
	//points_length and triangles_length, respectively
	points = Points;      
	triangles = Triangles;  
				
	//*******************************************************************************
	// VRML parser END
	//*******************************************************************************
	
	
	//Calculate epsilon
	double eps = 1000000000;
	//System.out.println("eps: " + eps);
	//System.out.println("points_length = " + points_length);
	//System.out.println("triangles_length = " + triangles_length);
	
	for(int i = 0; i < points_length; i++){
	    Point A;
	    A = points[i];
	    
	    if(eps > Math.abs(A.X) && Math.abs(A.X) != 0.0)
		eps = Math.abs(A.X);
	    if(eps > Math.abs(A.Y) && Math.abs(A.Y) != 0.0)
		eps = Math.abs(A.Y);
	    if(eps > Math.abs(A.Z) && Math.abs(A.Z) != 0.0)
		eps = Math.abs(A.Z);
	}
	
	ecomp.epsilon = Math.abs(eps / 1000.0);
	//System.out.println("epsilon = " + ecomp.epsilon);
	
	//calculates area of the model
	whole_area = calculateWholeArea(triangles, triangles_length);
		
	//create connectivity matrix (i.e., adjacency graph) for mesh model
	SparseMatrix sparsem = new SparseMatrix();
	Vector sparse_res = sparsem.createMatrix(triangles, triangles_length, points, points_length);
	points = (Point[]) sparse_res.elementAt(0);
	Integer points_length_o = (Integer) sparse_res.elementAt(1);
	points_length = points_length_o.intValue();
	sparse = (int[][]) sparse_res.elementAt(2);
	//System.out.println("Area of the object using Triangles: " + whole_area);
	
	//Resample faces until needed number of vertices is reached (specified at command line)
	System.out.println("Resampling faces...");
	Resample resample = new Resample();
	Vector resample_res = resample.doProcess(points_number, points, points_length/*, matrix*/, sparse);
	points = (Point[]) resample_res.elementAt(0);
	points_length_o = (Integer) resample_res.elementAt(1);
	points_length = points_length_o.intValue();
	sparse = (int[][]) resample_res.elementAt(2);
	System.out.println("After resampling, points_length = " + points_length);
	
	//Calculates mu values
	System.out.println("Calculating mu values...");
	MuApprox mu_val = new MuApprox();
	mu_values = mu_val.doProcess(mu_coeff, whole_area, points, points_length, sparse);
	
	//Normalize mu values
	MuNormalization mu_norm = new MuNormalization();
	mu_values = mu_norm.normalize(mu_values);
	
	//Construct MRG
	System.out.println("Constructing MRG...");
	MRGConstrLight mrg = new MRGConstrLight();
	Vector mrg_result = mrg.doProcess(MRG_number, points, points_length, sparse, mu_values, false);
	points = (Point[]) mrg_result.elementAt(0);
	points_length_o = (Integer) mrg_result.elementAt(1);
	points_length = points_length_o.intValue();	
	sparse = (int[][]) mrg_result.elementAt(2);
	mu_values = (double[]) mrg_result.elementAt(3);
	all_Tsets = (int[][]) mrg_result.elementAt(4);
	MRG = (int[][][]) mrg_result.elementAt(5);
	reebs = (ReebGraphElement[][]) mrg_result.elementAt(6);
	
	//Calculate attributes for the nodes in MRG
	System.out.println("Calculating attributes for Tsets...");
	AttributeCalculation a_calc = new AttributeCalculation();
	attributes = a_calc.doProcess(points, points_length, sparse, mu_values, all_Tsets, MRG, reebs, whole_area);
	
	//Write MRG into an '.mrg' text file 
	System.out.println("Saving MRG to a file...");
	SaveInText saver = new SaveInText();
	saver.saveIt(filename, MRG, reebs, attributes, mrg.FINEST_RESOLUTION);
	
	//System.out.println("All done!");
    }
    
    public static double calculateWholeArea(Triangle trigs[], int length){
	double area = 0;
	
	for(int i = 0; i < length; i++){
	    area = area + calculateTrigArea(points[trigs[i].a], points[trigs[i].b], points[trigs[i].c]);
	}
	
	return area;
    }
    
    public static double calculateTrigArea(Point A, Point B, Point C){
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
    
    public static  double calculateSparseArea(){
	
	int points_in_area[] = new int[points_length];
	
	for(int i = 0; i < points_length; i++){
	    points_in_area[i] = i;
	}
		
	int indexA, indexB, indexC;
	Point A, B, C;
	double area = 0;
	
	for(int i = 0; i < points_in_area.length; i++){
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
	
	if(points_in_area.length == 0){
	    area = 0;
	}
	
	if(points_in_area.length == 1){
	    area = 1;
	}
	
	if(points_in_area.length == 2){
	    area = 2;
	}
	
	//the_area = the_area + area;
	
	//base_areas.addElement(new Double(area));
	return area;
    }
    
    public static boolean isConnected(int index1, int index2){
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

