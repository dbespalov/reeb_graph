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

//Normalization of mu values

import java.util.*;
import java.io.*;
import java.util.Vector.*;
import java.lang.Math.*;
import java.lang.Object.*;
import java.lang.System.*;


public class MuNormalization {

    public double[] normalize(double mu_values[]){
	double minimum, maximum, temp;
		
	//looks for min and max
	temp = mu_values[0];
	minimum = mu_values[0];
	maximum = mu_values[0];
	for(int i = 0; i < mu_values.length; i++){
	    temp = mu_values[i];
	    
	    if(temp > maximum)
		maximum = temp;
	    
	    if(temp < minimum)
		minimum = temp;
	}
	
	//System.out.println("Minimum value of mu is: " + minimum);
	//System.out.println("Maximum value of mu is: " + maximum);
	
	//Normalize values using formula
	//MUn(V) = (MU(V) - minimum) / maximum;
	for(int i = 0; i < mu_values.length; i++){
	    temp = mu_values[i];
	    mu_values[i] = (temp - minimum) / maximum;
	}
	
	
	temp = mu_values[0];
	minimum = mu_values[0];
	maximum = mu_values[0];
	for(int i = 0; i < mu_values.length; i++){
	    temp = mu_values[i];
	    
	    if(temp > maximum)
		maximum = temp;
	    
	    if(temp < minimum)
		minimum = temp;
	}
	
	//System.out.println("Minimum value of mu is: " + minimum);
	//System.out.println("Maximum value of mu is: " + maximum);
	
	return mu_values;
    }
}
