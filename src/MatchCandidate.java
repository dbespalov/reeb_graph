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

// MatchCandidate maintains a list of R-Nodes (indicies) from a one MRG 
// that are selected as possible candidates for matching with 
// R-Nodes from the other MRG 
public class MatchCandidate {
    
    public Vector vector;
    
    public void cloneElement(MatchCandidate element){
	vector = (Vector) element.vector.clone();
    }
    
    public void addToElement(MatchCandidate element){
	vector = addIntegerVector(vector, element.vector);
    }
    
    public Vector addIntegerVector(Vector target, Vector source){
	Vector result = null;
	
	if(target == null && source != null)
	    result = new Vector();
	else if(target != null)
	    result = (Vector) target.clone();
	
	if(source != null){
	    for(int i = 0; i < source.size(); i++){
		Integer temp = (Integer) source.elementAt(i);
		
		if(! result.contains(temp))
		    result.addElement(temp);
	    }
	}
	
	return result;	
    }    
}
