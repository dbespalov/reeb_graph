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
 |  Matching a pair of MRGs proceeds as following:
 |  
 |  1. Two MRGs are stored in CompareReebGraph.MRG1 and CompareReebGraph.MRG2 vectors
 |     that contain RNode objects. 
 |
 |  2. Matching is performed starting with coarsest resolution. As such, each RNode
 |     object maintains RNode.children vector to store R-edges between R-nodes
 |     from a coarser to a finear resolution of MRG. 
 |
 |     Note that a list RNode.children is the same as the list stored in
 |     ReebGraphElement.parents, where instances of RNode and ReebGraphElement 
 |     contain data for the same R-node in MRG. 
 |     ReebGraphElement object is used by ExtractReebGraph program 
 |     (see comments in MRGConstrLight.java), while RNode object 
 |     is used by CompareReebGraph program. 
 */

public class CompareReebGraph {

    public static String file_name;
    public static int line_counter;
    
    public static Vector attributes1;
    public static Vector attributes2;
    
    public static Vector MRG1;
    public static Vector MRG2;
        
    public static double w;
    
    public static Vector NLIST;
    public static Vector MPAIR;
    public static double SIM_R_S;
    
    public static PrintWriter out;
    
    
    /*
    |  usage: java -Xmx512m  CompareReebGraph  <num_pts>  <mu_coeff>  <mrg_size>  <sim_weight>  <model_1>.wrl  <model_2>.wrl ... <model_N>.wrl
    |
    |    <num_pts>       -- target number of vertices before MRG is extracted 
    |                       (triangle meshes in each 3D model are resampled to match <num_pts>)
    |
    |    <mu_coeff>      -- coefficient for calculating threshold parameter r=sqrt(mu_coeff * area(S)),
    |                       used to approximate mu values
    |
    |    <mrg_size>      -- number of ranges K used in the finest resolution of MRG
    |   
    |    <sim_weight>    -- weight w used in similarity function (trade-off between attributes a and l)
    |	   
    |    <model_i>.wrl   -- a list of VRML models to compare, where i=[1,N]. 
    |                       It is assumed that each VRML model was processed using ExtractReebGraph program,
    |                       and MRG for model <model_i>.wrl is stored in text file <model_i>.mrg 
    */


    public static void main (String arg[]) throws IOException {
	int points_number = (new Integer(arg[0])).intValue();
	
	double mu_coeff = (new Double(arg[1])).doubleValue();
	
	int MRG_number = (new Integer(arg[2])).intValue();	
	
	w = (new Double(arg[3])).doubleValue();
	
	//the results of comparison will be saved in a file
	out = new PrintWriter(new FileWriter("log_" + points_number 
					     + "_" + mu_coeff + "_" + MRG_number + "_" + w));
		
	for(int i = 4; i < arg.length; i++){
	    for(int j = 4; j < arg.length; j++){
		String argument[] = new String [2];
		argument[0] = arg[i];
		argument[1] = arg[j];
		main_one(argument);
		
		out.println("Similarity between " + arg[i] + " and " + arg[j] + " is " + SIM_R_S);
	    }
	}
	
	out.close();
    }
    
    //function that compares a pair of MRGs
    public static void main_one (String arg[]) throws IOException {
	
	//System.out.println("Comparing " + arg[0] + " and " + arg[1] + " ...");
	
	attributes1 = null;
	attributes2 = null;
	
	MRG1 = null;
	MRG2 = null;
	
	//read both of the objects from mrg files
	readOneFile(arg[0]);
	readOneFile(arg[1]);
	
	if(MRG1.size() != MRG2.size()){
	    System.out.println("ERROR: number of resolutions must match! Exiting...");
	    System.exit(1);	    
	}
	
	//calculate the attributes for the pair of MRGs
	calculateRestAttributes(MRG1, attributes1);
	calculateRestAttributes(MRG2, attributes2);
	
	//computes parents in each MRG
	computeParents();
	
	//compare two MRGs 
	doComparison(arg);
	
	//performs simple tests for both objects
	//tester();
	
	attributes1 = null;
	attributes2 = null;
    }
    
    public static void doComparison(String str[]) throws IOException{
	NLIST = new Vector();
	MPAIR = new Vector();
	Vector list1 = new Vector();
	Vector list2 = new Vector();
	SIM_R_S = 0;
	
	//puts nodes from coarsest resolution graph in to NLIST
	Vector reeb_graph = (Vector) MRG1.elementAt(MRG1.size() - 1);
	for(int i = 0; i < reeb_graph.size(); i++){
	    Vector vec = (Vector) reeb_graph.elementAt(i);
	    RNode node = (RNode) vec.elementAt(0);
	    
	    list1.addElement(node);	    
	}
	
	reeb_graph = (Vector) MRG2.elementAt(MRG2.size() - 1);
	for(int i = 0; i < reeb_graph.size(); i++){
	    Vector vec = (Vector) reeb_graph.elementAt(i);
	    RNode node = (RNode) vec.elementAt(0);
	    
	    list2.addElement(node);	    
	}
	
	NLIST.addElement(list1);
	NLIST.addElement(list2);
	
	
	while(list1.size() != 0 && list2.size() != 0){
	    lookForMatchingPair(str);
	}
	



	
	/*//this code is just for image****************************************
	String image_filename1 = str[0].substring(0,str[0].indexOf(".wrl")) + "_" 
	    + str[1].substring(0,str[1].indexOf(".wrl")) + ".img" ;
	String image_filename2 = str[1].substring(0,str[1].indexOf(".wrl")) + "_" 
	    + str[0].substring(0,str[0].indexOf(".wrl")) + ".img";
	PrintWriter out_image1 = new PrintWriter(new FileWriter(image_filename1));
	PrintWriter out_image2 = new PrintWriter(new FileWriter(image_filename2));
	
	for(int i = 0; i < MPAIR.size(); i++){
	    MPairElement el = (MPairElement) MPAIR.elementAt(i);
	    
	    if(calculateIndexInMRG1(el.node1) == 0
	       && calculateIndexInMRG2(el.node2) == 0){
		
		if(el.vector_nodes1 == null){
		    out_image1.println(el.node1.index);
		}
		else{
		    for(int j = 0; j < el.vector_nodes1.size(); j++){
			RNode temp = (RNode) el.vector_nodes1.elementAt(j);
			out_image1.print(temp.index + " ");			
		    }
		    out_image1.println("");
		}
		
		if(el.vector_nodes2 == null){
		    out_image2.println(el.node2.index);
		}
		else{
		    for(int j = 0; j < el.vector_nodes2.size(); j++){
			RNode temp = (RNode) el.vector_nodes2.elementAt(j);
			out_image2.print(temp.index + " ");
		    }
		    out_image2.println("");
		}
	    }
	}
	
	out_image2.close();
	out_image1.close();
	//this code is just for image****************************************
	*/

	

	//Computes total similarity
	for(int i = 0; i < MPAIR.size(); i++){
	    MPairElement el = (MPairElement) MPAIR.elementAt(i);
	    SIM_R_S = SIM_R_S + sim(el.node1, el.node2);
	}

	System.out.print("Similarity  between "+str[0]+" and "+str[1]+" is: ");
	System.out.println(SIM_R_S);
    }
    
    //takes a node from NLIST that has the largest value of sim(m,m) and find a matching 
    //pair for this node
    public static void lookForMatchingPair(String str[]){
	Vector temp_vector, temp_vector2;
	temp_vector = new Vector();
	temp_vector2 = new Vector();	
	
	Vector list1 = (Vector) NLIST.elementAt(0);
	Vector list2 = (Vector) NLIST.elementAt(1);
	
	if(list1.size() == 0 || list2.size() == 0)
	    return;
	
	boolean is_first_time = true;
	int select_nlist2 = -1, index2 = -1;
	int same_range_node_num;
	
	int select_nlist, index;
	MPairElement pair;
	Vector tmp_v = findMaximumSim(); //returns a node that has the largest value for sim(m,m)
	pair = (MPairElement) tmp_v.elementAt(0);
	Integer select_nlist_o = (Integer) tmp_v.elementAt(1);
	Integer index_o = (Integer) tmp_v.elementAt(2);
	select_nlist = select_nlist_o.intValue();
	index = index_o.intValue();
	
	Vector n_list, vec_m, vec_n;
	RNode node, matched_node, temp_node, m = null, n = null, temp_node2, temp_node3, temp_node4;
	double range_diff, max_mat, temp_double;
	int what_res_m = -1, what_res_n = -1;
	Integer temp_int_m, temp_int_n;
	
	max_mat = -1000000000;
	
	//look for a matching pair

	//if(pair.node1 == null){
	if(select_nlist == 2){
	    n_list = (Vector) NLIST.elementAt(0);
	    node = pair.node2;
	}
	else{
	    n_list = (Vector) NLIST.elementAt(1);
	    node = pair.node1;
	}
	
	Vector vector_nodes;
	boolean first = true;
	int i = 0;
	
	boolean check = false;
	
	for(i = 0; i < n_list.size(); i++){
	    temp_node = (RNode) n_list.elementAt(i);
	    	    
	    //both nodes have to be from the same range
	    if(temp_node.left_bound == node.left_bound
	       && temp_node.right_bound == node.right_bound){
		
		//if(pair.node1 == null){
		if(select_nlist == 2){	
		    m = temp_node;
		    n = node;
		}
		else{
		    m = node;
		    n = temp_node;
		}
		
		range_diff = m.right_bound - m.left_bound;
		what_res_m = calculateIndexInMRG1(m);
		what_res_n = calculateIndexInMRG2(n);
		
		Vector temp_vec = (Vector) MRG1.elementAt(what_res_m);
		vec_m = (Vector) temp_vec.elementAt(m.index);
		
		temp_vec = (Vector) MRG2.elementAt(what_res_n);
		vec_n = (Vector) temp_vec.elementAt(n.index); 
		
		//parents of those nodes have to match
		if(parentsAreMatched(m, n)){
		    
		    //MLISTs for two nodes must be the same
		    if(n.MLIST.equals(m.MLIST)){
			
			int candidates_count = 1;
			Vector match_candidates = new Vector();
			RNode the_node = temp_node;
			
			if(select_nlist == 1)
			    same_range_node_num = howManySameRangeNodesInNList(the_node, 2);
			else
			    same_range_node_num = howManySameRangeNodesInNList(the_node, 1);
			
			// takes into account sequence of R-nodes in order to account for cases when mu values are divided into different ranges

			//select_nlist=1 indicate that CompareReebGraph.MRG2 has the largest sim(m,m) value
			//select_nlist=0 indicate that CompareReebGraph.MRG1 has the largest sim(m,m) value
			while(createMatchCandidates(the_node, n_list, match_candidates, candidates_count, select_nlist)){
			    
			    for(int h = 0; h < match_candidates.size(); h++){
				MatchCandidate candidate_element = (MatchCandidate) match_candidates.elementAt(h);
				
				if(candidates_count != 1){
				    temp_node = new RNode();
				    temp_node.index = the_node.index;
				    temp_node.left_bound = the_node.left_bound;
				    temp_node.right_bound = the_node.right_bound;
				    temp_node.parent = the_node.parent;
				    temp_node.children = copyChildren(null, the_node.children);
				    temp_node.MLIST = new Vector();
				    temp_node.MLIST = (Vector) the_node.MLIST.clone();
				    temp_node.attribute = new AttributeElement();
				    temp_node.attribute.a = the_node.attribute.a;
				    temp_node.attribute.l = the_node.attribute.l;
				    
				    for(int g = 1; g < candidate_element.vector.size(); g++){
					Integer tmp_in = (Integer) candidate_element.vector.elementAt(g);
					temp_node4 = (RNode) n_list.elementAt(tmp_in.intValue());
					temp_node.children = copyChildren(temp_node.children, temp_node4.children);
					
					temp_node.attribute.a = temp_node.attribute.a + temp_node4.attribute.a;
					temp_node.attribute.l = temp_node.attribute.l + temp_node4.attribute.l;
				    }
				    
				    vector_nodes = candidate_element.vector;
				}
				else{
				    vector_nodes = null;
				    temp_node = the_node;
				}
				
				if(select_nlist == 2){	
				    m = temp_node;
				    n = node;
				}
				else{
				    m = node;
				    n = temp_node;
				}
				
				temp_double = mat(m, n);
				
				if(max_mat < temp_double){
				    max_mat = temp_double;
				    
				    if(select_nlist == 2){
					pair.node1 = temp_node;
					
					if(vector_nodes == null)
					    temp_vector = null;
					else
					    temp_vector = (Vector) vector_nodes.clone();
					
					select_nlist2 = 1;
					index2 = i;
				    }
				    else{
					pair.node2 = temp_node;
					
					if(vector_nodes == null)
					    temp_vector = null;
					else
					    temp_vector = (Vector) vector_nodes.clone();
					
					select_nlist2 = 2;
					index2 = i;
				    }
				}
				
				if(select_nlist == select_nlist2)
				    System.out.println("WARNING in CompareReebGraph.lookForMatchingPair(): select_nlist and select_nlist2 point to the same MRG!");
				
				temp_int_m = (Integer) m.MLIST.elementAt(m.MLIST.size() - 1);
				temp_int_n = (Integer) n.MLIST.elementAt(n.MLIST.size() - 1);
				
				//propagate labels
				if(vector_nodes == null){
				    for(int j = 1; j < vec_m.size(); j++){
					temp_node2 = (RNode) vec_m.elementAt(j);
					
					if(temp_node2.left_bound == m.right_bound
					   && temp_node2.right_bound == m.right_bound + range_diff
					   && temp_int_m.intValue() > 0){
					    
					    
					    temp_node2.MLIST.addElement(new Integer(temp_int_m.intValue() + 1));
					}
					
					if(temp_node2.right_bound == m.left_bound
					   && temp_node2.left_bound == m.left_bound - range_diff
					   && temp_int_m.intValue() < 0){
					    
					    temp_node2.MLIST.addElement(new Integer(temp_int_m.intValue() - 1));    
					}
				    }
				    
				    for(int j = 1; j < vec_n.size(); j++){
					temp_node2 = (RNode) vec_n.elementAt(j);
					
					if(temp_node2.left_bound == n.right_bound
					   && temp_node2.right_bound == n.right_bound + range_diff
					   && temp_int_n.intValue() > 0){
					    
					    temp_node2.MLIST.addElement(new Integer(temp_int_n.intValue() + 1));
					}
					
					if(temp_node2.right_bound == n.left_bound
					   && temp_node2.left_bound == n.left_bound - range_diff
					   && temp_int_n.intValue() < 0){
					    
					    temp_node2.MLIST.addElement(new Integer(temp_int_n.intValue() - 1));
					}					
				    }
				}
			    }
			    
			    match_candidates = new Vector();
			    // if(same_range_node_num > 20 && candidates_count > 2 && candidates_count < 6){
			    // 	candidates_count = 6;
			    // }
			    // else
			    candidates_count++;
			}
		    }
		}
	    }
	}
	
    	// if a matching mair can not be found then R-node is removed from NLIST
	// and lookForMatchingPair(str) is called recursively
	if(pair.node1 == null || pair.node2 == null){
	    
	    if(select_nlist == 1){
		RNode nn = (RNode) list1.remove(index);
	    }
	    else{
		RNode nn = (RNode) list2.remove(index);
	    }
	    
	    pair = null;
	    lookForMatchingPair(str);
	}
	
	//removes matching nodes from NLISTs and add this pair to MPAIR list
	else{
	    
	    Vector new_vec = null;
	    if(temp_vector != null){
		new_vec = new Vector();
		for(int r = 0; r < temp_vector.size(); r++){
		    Integer ti = (Integer) temp_vector.elementAt(r);
		    
		    RNode tn = (RNode) n_list.elementAt(ti.intValue());
		    new_vec.addElement(tn);
		}
	    }
	    
	    if(select_nlist == 1){
		pair.vector_nodes1 = null;
		if(temp_vector == null){
		    pair.vector_nodes2 = null;
		}
		else{
		    pair.vector_nodes2 = new_vec;
		}
	    }
	    else{
		pair.vector_nodes2 = null;
		if(temp_vector == null){
		    pair.vector_nodes1 = null;
		}
		else{
		    pair.vector_nodes1 = new_vec;
		}
	    }	    
	    	    
	    MPAIR.addElement(pair);
	    
	    m = pair.node1;
	    n = pair.node2;
	    
	    if(select_nlist == 1)
		list1.removeElementAt(index);
	    else
		list2.removeElementAt(index);
	    
	    if(temp_vector == null){
		if(select_nlist2 == 1)
		    list1.removeElementAt(index2);
		else
		    list2.removeElementAt(index2);
		
	    }
	    else{
		
		for(int h = 0; h < temp_vector.size(); h++){
		    Integer tmp_i2 = (Integer) temp_vector.elementAt(h);
		    
		    RNode ntemp = (RNode) n_list.elementAt(tmp_i2.intValue());
				    
		    n_list.setElementAt(null, tmp_i2.intValue());

		}
				
		int h = 0;
		while(h < n_list.size()){
		    RNode tmp_i2 = (RNode) n_list.elementAt(h);
		    
		    if(tmp_i2 == null)
			n_list.removeElementAt(h);
		    else
			h++;
		}
	    }
	    
	    if(select_nlist == select_nlist2)
		System.out.println("WARNING in CompareReebGraph.lookForMatchingPair(): select_nlist and select_nlist2 point to the same MRG!");
	    	    
	    what_res_m = calculateIndexInMRG1(m);
	    what_res_n = calculateIndexInMRG2(n);
	    
	    if(m.children != null && what_res_m > 0){
		vec_m = (Vector) MRG1.elementAt(what_res_m - 1);
		for(int z = 0; z < m.children.size(); z++){
		    Integer temp_int = (Integer) m.children.elementAt(z);
		    
		    Vector tv = (Vector) vec_m.elementAt(temp_int.intValue());
		    temp_node = (RNode) tv.elementAt(0);
		    
		    if(! list1.contains(temp_node))
			list1.addElement(temp_node);
		}
	    }
	    
	    if(n.children != null && what_res_n > 0){
		vec_n = (Vector) MRG2.elementAt(what_res_n - 1);
		for(int z = 0; z < n.children.size(); z++){
		    Integer temp_int = (Integer) n.children.elementAt(z);
		    Vector tv = (Vector) vec_n.elementAt(temp_int.intValue());
		    temp_node = (RNode) tv.elementAt(0);
		    
		    if(! list2.contains(temp_node))
			list2.addElement(temp_node);
		}
	    }
	}
    }

    //calculates the number of nodes that are connected to R-node node
    public static int howManySameRangeNodesInNList(RNode node, int select_nlist){
	int counter = 1;
	Vector nlist;
	if(select_nlist == 1)
	    nlist = (Vector) NLIST.elementAt(0);
	else
	    nlist = (Vector) NLIST.elementAt(1);
	
	Vector vector_nodes = new Vector();
	for(int i = 0; i < nlist.size(); i++){
	    RNode temp_node = (RNode) nlist.elementAt(i);
	    
	    vector_nodes.addElement(new Integer(i));
	    
	    if(temp_node.left_bound == node.left_bound
	       && temp_node.right_bound == node.right_bound
	       && isConnectedToOneNode(node, vector_nodes, nlist, select_nlist, 5)){
		counter++;
	    }
	    else
		vector_nodes.removeElementAt(vector_nodes.size() - 1);
	}
	
	return counter;
    }
    
    //select_nlist=1 indicate that CompareReebGraph.MRG2 has the largest sim(m,m) value
    //select_nlist=0 indicate that CompareReebGraph.MRG1 has the largest sim(m,m) value
    public static boolean createMatchCandidates(RNode node, Vector n_list, Vector match_candidates, int candidates_count, int select_nlist){
	
	//only takes a look at the sequences that are no longer than 5 nodes
	if(candidates_count > 6)
	    return false;
		
	//this is used so that a sequence with the maximum length can be taken into account
	if(candidates_count == 6){
		    
	    int node_index_in_nlist = -1;
	    MatchCandidate element = new MatchCandidate();
	    RNode temp_node;
	    
	    for(int i = 0; i < n_list.size(); i++){
				
		temp_node = (RNode) n_list.elementAt(i);
		
		if(temp_node.index == node.index 
		   && temp_node.left_bound == node.left_bound
		   && temp_node.right_bound == node.right_bound
		   && temp_node.attribute.a == node.attribute.a
		   && temp_node.attribute.l == node.attribute.l){
		    
		    node_index_in_nlist = i;
		    element.vector = new Vector();
		    element.vector.addElement(new Integer(node_index_in_nlist));
		}
	    }
	    
	    if(node_index_in_nlist == -1){
		System.out.println("ERROR in CompareReebGraph.createMatchCandidates(): can not find R-node in n_list that match ranges of R-node node. Exiting... ");
		System.exit(1);
	    }
		
	    for(int i = 0; i < n_list.size(); i++){
		
		temp_node = (RNode) n_list.elementAt(i);
		
		if(node_index_in_nlist != i 
		   && temp_node.left_bound == node.left_bound
		   && temp_node.right_bound == node.right_bound){
		    
		    int mrg_num;
		    if(select_nlist == 1)
			mrg_num = 2;
		    else
			mrg_num = 1;
		    
		    Integer node_index = new Integer(i);
		    if(! element.vector.contains(node_index)
		       && isConnectedToOneNode(temp_node, element.vector, n_list, mrg_num, candidates_count)){
			element.vector.addElement(node_index);
		    }
		}
	    }
	    
	    match_candidates.addElement(element);
	    	    
	    if(element.vector.size() > 5)
		return true;
	    else
		return false;
	}
	
	//computes sequences that are 1 to 5 nodes long
	int counter = 0;
	boolean stopper = false;
	int node_index_in_nlist = -1;
	
	MatchCandidate element = new MatchCandidate();
	RNode temp_node;
	
	for(int i = 0; i < n_list.size(); i++){
	    temp_node = (RNode) n_list.elementAt(i);
	    
	    if(temp_node.index == node.index 
	       && temp_node.left_bound == node.left_bound
	       && temp_node.right_bound == node.right_bound
	       && temp_node.attribute.a == node.attribute.a
	       && temp_node.attribute.l == node.attribute.l){
		
		node_index_in_nlist = i;
		element.vector = new Vector();
		element.vector.addElement(new Integer(node_index_in_nlist));
	    }
	}
	
	if(node_index_in_nlist == -1){
	    System.out.println("ERROR in CompareReebGraph.createMatchCandidates(): can not find R-node in n_list that match ranges of R-node node. Exiting... ");
	    System.exit(1);
	}
	
	counter++;
	match_candidates.addElement(element);
	
	if(counter == candidates_count){
	    return true;
	}
		
	boolean found_one = false;
	
	while(stopper != true){
	    match_candidates.addElement(null);
	    element = (MatchCandidate) match_candidates.remove(0);
	    
	    found_one = false;
	    while(element != null){
		for(int i = 0; i < n_list.size(); i++){
		    temp_node = (RNode) n_list.elementAt(i);
		    
		    Integer tmp_int;
		    if(element.vector.size() > 1)
			tmp_int = (Integer) element.vector.elementAt(element.vector.size() - 1);
		    else
			tmp_int = new Integer(-1);
		    
		    if(tmp_int.intValue() < i && node_index_in_nlist != i 
		       && temp_node.left_bound == node.left_bound
		       && temp_node.right_bound == node.right_bound){
			
			int mrg_num;
			if(select_nlist == 1)
			    mrg_num = 2;
			else
			    mrg_num = 1;
			
			Integer node_index = new Integer(i);
			if(! element.vector.contains(node_index)
			   && isConnectedToOneNode(temp_node, element.vector, n_list, mrg_num, candidates_count)){
			    
			    
			    MatchCandidate new_el = new MatchCandidate();
			    new_el.cloneElement(element);
			    new_el.vector.addElement(node_index);
						    
			    found_one = true;
			    
			    match_candidates.addElement(new_el);
			}
		    }
		}		    
		element = (MatchCandidate) match_candidates.remove(0);
	    }
	    
	    if(found_one == true){
		counter++;
	    }
	    else if(found_one == false){
		stopper = true;
	    }
	    
	    if(counter == candidates_count)
		return true;
	}
	
	return false;
    }
    
    // returns union of two sets of integers (target and new_children) 
    public static Vector copyChildren(Vector target, Vector new_children){
	Vector result = null;
	if(target == null && new_children != null)
	    result = new Vector();
	else if(target != null)
	    result = (Vector) target.clone();
	
	if(new_children != null){
	    for(int i = 0; i < new_children.size(); i++){
		Integer temp = (Integer) new_children.elementAt(i);
		
		if(! result.contains(temp))
		    result.addElement(temp);
	    }
	}
	
	return result;	
    }
    
    //checks that node and all R-nodes whos indices are stored in vector_nodes are connected to one R-node
    public static boolean isConnectedToOneNode(RNode node, Vector vector_nodes, Vector n_list, int what_MRG, int candidates_count){
	if(vector_nodes == null)
	    return false;
	
	Integer temp_int = (Integer) vector_nodes.elementAt(0);
	Vector reeb_graph, vec, temp_vec;
	RNode main_node = (RNode) n_list.elementAt(temp_int.intValue());
	RNode temp_node, temp_node2;
	boolean is_connected = true;
						   
	if(what_MRG == 1)
	    reeb_graph = (Vector) MRG1.elementAt(calculateIndexInMRG1(main_node));
	else
	    reeb_graph = (Vector) MRG2.elementAt(calculateIndexInMRG2(main_node));
	
	vec = (Vector) reeb_graph.elementAt(main_node.index);
	
	for(int i = 1; i < vec.size(); i++){
	    temp_node = (RNode) vec.elementAt(i);
	    
	    for(int j = 1; j < vector_nodes.size(); j++){
		temp_int = (Integer) vector_nodes.elementAt(j);
		temp_node2 = (RNode) n_list.elementAt(temp_int.intValue());
		
		temp_vec = (Vector) reeb_graph.elementAt(temp_node2.index);
		
		if(! temp_vec.contains(temp_node)){
		    is_connected = false;
		    break;
		}	
	    }
	    
	    if(is_connected == true){
		Vector tmp = (Vector) reeb_graph.elementAt(node.index);
		if(tmp.contains(temp_node) && ! node.equals(temp_node))
		    return true;
	    }
	    
	    is_connected = true;
	}
	
	return false;	
    }

    public static boolean allAlone(RNode node, Vector vector_nodes, Vector n_list, int what_MRG){
	for(int i = 0; i < vector_nodes.size(); i++){
	    Integer temp_int = (Integer) vector_nodes.elementAt(i);
	    
	    RNode temp = (RNode) n_list.elementAt(temp_int.intValue());
	    
	    if(! isAlone(temp, what_MRG))
		return false;
	}
	
	return isAlone(node, what_MRG);
    }
    
    public static boolean isAlone(RNode node, int what_MRG){
	Vector mrg;
	int index_in_mrg;
	if(what_MRG == 1){
	    mrg = MRG1;
	    index_in_mrg = calculateIndexInMRG1(node);
	}
	else{
	    mrg = MRG2;
	    index_in_mrg = calculateIndexInMRG2(node);
	}
	
	Vector reeb_graph = (Vector) mrg.elementAt(index_in_mrg);
	Vector adj_vec = (Vector) reeb_graph.elementAt(node.index);
	
	RNode temp = (RNode) adj_vec.elementAt(0);
	
	if(adj_vec.size() == 1 && temp.index == node.index
	   && (node.children == null || node.children.size() == 1))
	    return true;
	else
	    return false;
    }
        
    public static Vector findMaximumSim(){
	Vector result = new Vector();
	int select_nlist, index;
	
	Vector NLIST1 = (Vector) NLIST.elementAt(0);
	Vector NLIST2 = (Vector) NLIST.elementAt(1);
	RNode temp_node, maximum_node;
	double maximum_sim, temp;
	int MRG_no = 1;
	
	temp_node = (RNode) NLIST1.elementAt(0);
	maximum_sim = sim(temp_node, temp_node);
	maximum_node = temp_node;
	index = 0;
	
	for(int i = 1; i < NLIST1.size(); i++){
	    temp_node = (RNode) NLIST1.elementAt(i);
	    temp = sim(temp_node, temp_node);
	    
	    if(temp > maximum_sim){
		maximum_sim = temp;
		maximum_node = temp_node;
		MRG_no = 1;
		index = i;
	    }	    
	}
	
	for(int i = 0; i < NLIST2.size(); i++){
	    temp_node = (RNode) NLIST2.elementAt(i);
	    temp = sim(temp_node, temp_node);
	    
	    if(temp > maximum_sim){
		maximum_sim = temp;
		maximum_node = temp_node;
		MRG_no = 2;
		index = i;
	    }
	}
	
	MPairElement element = new MPairElement();
	
	if(MRG_no == 1){
	    element.node1 = maximum_node;	    
	    element.node2 = null;
	}
	else{
	    element.node1 = null;	    
	    element.node2 = maximum_node;
	}
	
	select_nlist = MRG_no;
	
	result.addElement(element);
	result.addElement(new Integer(select_nlist));
	result.addElement(new Integer(index));
	
	return result;
    }

    //R-node m must be from MRG1, and R-node n is from MRG2
    public static boolean parentsAreMatched(RNode m, RNode n){
	if(m.parent == null || n.parent == null)
	    return true;
	
	for(int i = 0; i < MPAIR.size(); i++){
	    MPairElement el = (MPairElement) MPAIR.elementAt(i);
	    
	    if(el.vector_nodes1 == null && el.vector_nodes2 == null){
		if(m.parent.equals(el.node1) && n.parent.equals(el.node2))
		    return true;
	    }
	    else{
		if(el.vector_nodes1 != null){
		    if(el.vector_nodes1.contains(m.parent) && n.parent.equals(el.node2))
			return true;
		}
		else{
		    if(el.vector_nodes2.contains(n.parent) && m.parent.equals(el.node1))
			return true;
		}
	    }
	}
	
	return false;
    }
    
    public static void computeParents(){
	RNode element, element2;
	Integer temp_int;
	Vector reeb_graph, next_reeb_graph;
	
	//calculating parent for each RNode in MRG1
	for(int i = 1; i < MRG1.size(); i++){
	    reeb_graph = (Vector) MRG1.elementAt(i);
	    
	    for(int j = 0; j < reeb_graph.size(); j++){
		Vector vec = (Vector) reeb_graph.elementAt(j);
		element = (RNode) vec.elementAt(0);
		
		if(i == MRG1.size() - 1)
		    element.parent = null;
		
		if(element.children != null){
		    
		    for(int k = 0; k < element.children.size(); k++){
			temp_int = (Integer) element.children.elementAt(k);
			
			next_reeb_graph = (Vector) MRG1.elementAt(i - 1);
			Vector vec2 = (Vector) next_reeb_graph.elementAt(temp_int.intValue());
			element2 = (RNode) vec2.elementAt(0);
			
			element2.parent = element;
		    }		    
		}		
	    }	    
	}

	//calculating parent for each RNode in MRG2
	for(int i = 1; i < MRG2.size(); i++){
	    reeb_graph = (Vector) MRG2.elementAt(i);
	    
	    for(int j = 0; j < reeb_graph.size(); j++){
		Vector vec = (Vector) reeb_graph.elementAt(j);
		element = (RNode) vec.elementAt(0);
		
		if(element.children != null){
		    
		    for(int k = 0; k < element.children.size(); k++){
			temp_int = (Integer) element.children.elementAt(k);
			
			next_reeb_graph = (Vector) MRG2.elementAt(i - 1);
			Vector vec2 = (Vector) next_reeb_graph.elementAt(temp_int.intValue());
			element2 = (RNode) vec2.elementAt(0);
			
			element2.parent = element;
		    }		    
		}		
	    }	    
	}
    }

    public static double mat(RNode m, RNode n){
	double loss_m_n = - loss(m, n);
	
	double range_length_m = m.right_bound - m.left_bound;
	double range_length_n = n.right_bound - n.left_bound;
	
	RNode adj_m1 = adj(m, 1, m.right_bound, m.right_bound + range_length_m);
	RNode adj_m2 = adj(m, 1, m.left_bound - range_length_m, m.left_bound);
	
	RNode adj_n1 = adj(n, 2, n.right_bound, n.right_bound + range_length_n);
	RNode adj_n2 = adj(n, 2, n.left_bound - range_length_n, n.left_bound);
	
	double loss_adj1 = loss(adj_m1, adj_n1);
	double loss_adj2 = loss(adj_m2, adj_n2);
	
	return (loss_m_n - (loss_adj1 + loss_adj2));	
    }
    
    public static double loss(RNode m, RNode n){

	return ((0.5 * (sim(m, m) + sim(n, n))) - sim(m, n));
    }
    
    public static RNode adj(RNode m, int what_MRG, double left_b, double right_b){
	Vector reeb_graph = null;
	
	if(what_MRG == 1){
	    reeb_graph = (Vector) MRG1.elementAt(calculateIndexInMRG1(m));
	}
	
	else if(what_MRG == 2){
	    reeb_graph = (Vector) MRG2.elementAt(calculateIndexInMRG2(m));
	}
	
	else{
	    System.out.println("ERROR in CompareReebGraph.adj(): wrong MRG index! Exiting...");
	    System.exit(1);
	}
	
	Vector adj_vers = (Vector) reeb_graph.elementAt(m.index);
	RNode result = new RNode();
	result.attribute = new AttributeElement();
	result.attribute.a = 0.0;
	result.attribute.a = 0.0;

	for(int i = 1; i < adj_vers.size(); i++){
	    RNode node = (RNode) adj_vers.elementAt(i);
	    
	    if(node.left_bound == left_b || node.right_bound == right_b){
		result.attribute.a = result.attribute.a + node.attribute.a;
		result.attribute.l = result.attribute.l + node.attribute.l;
	    }	    
	}
	
	return result;	
    }
    
    public static int calculateIndexInMRG1(RNode m){
	double range_length = m.right_bound - m.left_bound;
	double n = - (Math.log(range_length)) / (Math.log(2));
		
	int result = (int) n;
	return (MRG1.size() - 1 - result);
    }
    
    public static int calculateIndexInMRG2(RNode m){
	double range_length = m.right_bound - m.left_bound;	
	double n = - (Math.log(range_length)) / (Math.log(2));
	
	int result = (int) n;
	return (MRG2.size() - 1 - result);
    }
            
    public static double sim(RNode m, RNode n){
	double min_a, min_l;
	if(m.attribute.a > n.attribute.a)
	    min_a = n.attribute.a;
	else
	    min_a = m.attribute.a;
	
	if(m.attribute.l > n.attribute.l)
	    min_l = n.attribute.l;
	else
	    min_l = m.attribute.l;
	
	return (w * min_a + (1 - w) * min_l);
    }
        
    public static void calculateRestAttributes(Vector MRG, Vector attributes){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
	
	for(int i = 0; i < reeb_graph.size(); i++){
	    Vector elements = (Vector) reeb_graph.elementAt(i);
	    RNode node = (RNode) elements.elementAt(0);
	    
	    node.attribute = new AttributeElement();
	    AttributeElement temp = (AttributeElement) attributes.elementAt(i);
	    node.attribute.a = temp.a;
	    node.attribute.l = temp.l;	    
	}
		
	for(int i = 1; i < MRG.size(); i++){
	   reeb_graph = (Vector) MRG.elementAt(i);
	    
	    for(int j = 0; j < reeb_graph.size(); j++){
		Vector elements = (Vector) reeb_graph.elementAt(j);
		RNode node = (RNode) elements.elementAt(0);
		
		node.attribute = new AttributeElement();
		node.attribute.a = 0.0;
		node.attribute.l = 0.0;
		Vector prev_graph = (Vector) MRG.elementAt(i - 1);
		
		for(int k = 0; k < node.children.size(); k++){
		    Integer temp_int = (Integer) node.children.elementAt(k);
		    Vector adj_vec = (Vector) prev_graph.elementAt(temp_int.intValue());
		    RNode temp_node = (RNode) adj_vec.elementAt(0);		   
		    
		    node.attribute.a = node.attribute.a + temp_node.attribute.a;
		    node.attribute.l = node.attribute.l + temp_node.attribute.l;		    
		}
	    }
	}
    }
       
    public static void readOneFile(String filename) throws IOException {
	Error error = new Error();

	try{
	    filename = filename.substring(0,filename.indexOf(".wrl")) + ".mrg";
	    
	    BufferedReader in = new BufferedReader(new FileReader(filename));
	    String str;
	    StringTokenizer st;
	    
	    line_counter = 1;
	    file_name = filename;
	    
	    str = in.readLine();

	    if(! str.equalsIgnoreCase("attributes{")){
		throw(error);
	    }
	    
	    Vector attributes = new Vector();
	    line_counter++;
	    str = in.readLine();
	    
	    int att_num = Integer.valueOf(str).intValue();
	    for(int i = 0; i < att_num; i++){
		line_counter++;
		str = in.readLine();
		
		st = new StringTokenizer(str, " ");
		double d1, d2;
		String st1, st2;
		
		st1 = st.nextToken();
		st2 = st.nextToken();
		
		if(st1.equalsIgnoreCase("NaN"))
		    st1 = "0.0";
		if(st2.equalsIgnoreCase("NaN"))
		    st2 = "0.0";
		
		d1 = Double.valueOf(st1).doubleValue();
		d2 = Double.valueOf(st2).doubleValue();
		
		AttributeElement a_element = new AttributeElement();
		
		a_element.a = d1;
		a_element.l = d2;
		
		attributes.addElement(a_element);

		st = null;
	    }
	    
	    line_counter++;
	    str = in.readLine();
	    
	    if(! str.equalsIgnoreCase("}")){
		throw(error);
	    }
	    
	    line_counter++;
	    str = in.readLine();
	    int MRG_num = Integer.valueOf(str).intValue();
	    Vector MRG = new Vector();
	    
	    for(int i = 0; i < MRG_num; i++){
		Vector reeb_graph = new Vector();
		
		line_counter++;
		str = in.readLine();
		if(! str.equalsIgnoreCase("elements{")){
		    throw(error);
		}
		
		line_counter++;
		str = in.readLine();
		int graph_size = Integer.valueOf(str).intValue();
		
		for(int j = 0; j < graph_size; j++){
		    //reads index
		    line_counter++;
		    str = in.readLine();
		    int index = Integer.valueOf(str).intValue();
		    
		    //reads range
		    line_counter++;
		    str = in.readLine();
		    st = new StringTokenizer(str, " ");
		    double left_b, right_b;
		    
		    left_b = Double.valueOf(st.nextToken()).doubleValue();
		    right_b = Double.valueOf(st.nextToken()).doubleValue();
		    
		    //reads Tsets
		    line_counter++;
		    str = in.readLine();
		    st = new StringTokenizer(str, " ");
		    
		    Vector Tsets = new Vector();
		    while(st.hasMoreTokens()){
			Integer temp_integer = new Integer(Integer.valueOf(st.nextToken()).intValue());
			Tsets.addElement(temp_integer);			
		    }

		    //reads children (in SaveGraph.java it's called parents,
		    //here it has a better common sense to call children)
		    line_counter++;
		    str = in.readLine();
		    Vector children;
		    
		    if(str.equalsIgnoreCase("NULL")){
			children = null;
		    }
		    else{
			st = new StringTokenizer(str, " ");
			children = new Vector();
			while(st.hasMoreTokens()){
			    Integer temp_integer2 = new Integer(Integer.valueOf(st.nextToken()).intValue());
			    children.addElement(temp_integer2);			
			}			
		    }
		    
		    RNode el = new RNode();
		    el.index = index;
		    el.left_bound = left_b;
		    el.right_bound = right_b;
		    el.Tsets = Tsets;
		    el.children = children;
		    
		    el.MLIST = new Vector();
		    el.MLIST.addElement(new Integer(0));
		    
		    Vector vec = new Vector();
		    
		    vec.addElement(el);
		    reeb_graph.addElement(vec);
		}
		
		line_counter++;
		str = in.readLine();
		
		if(! str.equalsIgnoreCase("}")){
		    throw(error);
		}
		
		line_counter++;
		str = in.readLine();
		
		if(! str.equalsIgnoreCase("connectivity{")){
		    throw(error);
		}
		
		line_counter++;
		str = in.readLine();
		graph_size = Integer.valueOf(str).intValue();
		
		for(int j = 0; j < graph_size; j++){    
		    line_counter++;
		    str = in.readLine();
		    st = new StringTokenizer(str, " ");
		    
		    int main_ver_ind = Integer.valueOf(st.nextToken()).intValue();
		    
		    Vector adj_ver = (Vector) reeb_graph.elementAt(main_ver_ind);
		    
		    while(st.hasMoreTokens()){
			int ver_to_add = Integer.valueOf(st.nextToken()).intValue();
			
			Vector vec_el_to_add = (Vector) reeb_graph.elementAt(ver_to_add);
			RNode el_to_add = (RNode) vec_el_to_add.elementAt(0);
			adj_ver.addElement(el_to_add);
		    }
		    
		}
		
		line_counter++;
		str = in.readLine();
		
		if(! str.equalsIgnoreCase("}")){
		    throw(error);
		}
		
		MRG.addElement(reeb_graph);
	    }
	    
	    if(attributes1 == null && MRG1 == null){
		attributes1 = attributes;
		MRG1 = MRG;
	    }

	    else if(attributes2 == null && MRG2 == null){
		attributes2 = attributes;
		MRG2 = MRG;
	    }
	    else{
		System.out.println("WARNING: can only compare two MRGs at a time!");
	    }
	    	    
	    in.close();
	}
	
	catch(Error e){
	    System.out.print("The input file ");
	    System.out.print(file_name);
	    System.out.print(" contains wrong input at line ");
	    System.out.print(line_counter);
	    System.out.println("! Program terminating.....");
	    System.exit(1);
	}	
    }




    public static String getNodeName(RNode element){
	String node_name;
	Integer temp_integer = new Integer(element.index);
	Double temp_double = new Double(element.left_bound * Math.pow(2, MRG1.size() - 1));
	
	Integer temp = new Integer(temp_double.intValue() + 1);
		
	node_name = "v" + temp_integer.toString() + "_" + temp.toString() + "_" + calculateIndexInMRG1(element);

	return node_name;
    }

    
    // public static void tester(){
    // 	Vector nlist1 = (Vector) NLIST.elementAt(0);
    // 	Vector nlist2 = (Vector) NLIST.elementAt(1);
	
    // 	System.out.print("The size of NLIST1 is:");
    // 	System.out.println(nlist1.size());
	
    // 	System.out.print("The size of NLIST2 is:");
    // 	System.out.println(nlist2.size());
		
    // 	double total_sim1 = 0;
    // 	double total_sim2 = 0;

    // 	for(int i = 0; i < MRG1.size(); i++){
    // 	    Vector reeb_graph = (Vector) MRG1.elementAt(i);
	    
    // 	    for(int j = 0; j < reeb_graph.size(); j++){
    // 		Vector vec = (Vector) reeb_graph.elementAt(j);
		
    // 		RNode node = (RNode) vec.elementAt(0);
		
    // 		total_sim1 = total_sim1 + sim(node, node);		
    // 	    }
    // 	}
	
    // 	for(int i = 0; i < MRG2.size(); i++){
    // 	    Vector reeb_graph = (Vector) MRG2.elementAt(i);
	    
    // 	    for(int j = 0; j < reeb_graph.size(); j++){
    // 		Vector vec = (Vector) reeb_graph.elementAt(j);
		
    // 		RNode node = (RNode) vec.elementAt(0);
		
    // 		total_sim2 = total_sim2 + sim(node, node);		
    // 	    }
    // 	}
	
    // 	System.out.print("Total sim for MRG1: ");
    // 	System.out.println(total_sim1);
    // 	System.out.print("Total sim for MRG2: ");
    // 	System.out.println(total_sim2);
    // }
        

}

