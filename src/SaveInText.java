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

// Save MRG into a '.mrg' text file

public class SaveInText {
    
    public void saveIt(String filename, int MRG[][][], ReebGraphElement reebs[][], 
		       AttributeElement attributes[], int res_num)  throws IOException {
	
	String graph_filename = filename.substring(0,filename.indexOf(".wrl")) + ".mrg";
	    
	PrintWriter out = new PrintWriter(new FileWriter(graph_filename));
	
	//writing attributes
	out.println("attributes{");
	out.println(attributes.length);
	
	for(int i = 0; i < attributes.length; i++){
	    AttributeElement element = attributes[i];
	    
	    out.print(element.a);
	    out.print(" ");
	    out.println(element.l);
	}
	    
	out.println("}");
	
	//saving size of MRG
	out.println(MRG.length);
	
	//writing the Reeb Graph Elements
	
	for(int i = 0; i < MRG.length; i++){
	    
	    out.println("elements{");
	    out.println(reebs[i].length);
	    
	    for(int j = 0; j < reebs[i].length; j++){
		ReebGraphElement el = reebs[i][j];
		
		out.println(el.index);
		
		out.print(el.left_bound);
		out.print(" ");
		out.println(el.right_bound);
		
		for(int k = 0; k < el.Tsets.size(); k++){
		    Integer temp_i = (Integer) el.Tsets.elementAt(k);
		    out.print(temp_i.intValue());
		    out.print(" ");		    
		}
		out.println("");
		
		if(el.parents != null){
		    for(int k = 0; k < el.parents.size(); k++){
			Integer temp_i = (Integer) el.parents.elementAt(k);
			out.print(temp_i.intValue());
			out.print(" ");		    
		    }
		}
		else{
		    out.print("NULL");
		}
		
		out.println("");
	    }

	    out.println("}");
	    
	    out.println("connectivity{");
	    out.println(MRG[i].length);
	    
	    for(int j = 0; j < MRG[i].length; j++){
		out.print(j);
		out.print(" ");
		for(int k = 1; k < MRG[i][j][0]; k++){
		    out.print(MRG[i][j][k]);
		    out.print(" ");
		}

		out.println("");
	    }
	    
	    out.println("}");
		
	}
	
	out.close();
    }

}
