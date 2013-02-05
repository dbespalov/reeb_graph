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

    /*
    |  Construction of multiresolutional Reeb graph (MRG) proceeds as following:
    |
    |  1. Creation of a reeb graph at the finest resolution 
    |     with MRGConstrLight.createFinestResolutionReebGraph() method.
    |
    |  2. Reeb graphs at coarser levels are created by MRGConstrLight.createMRG() method
    |     that recursively merges neighboring ranges of the normalized mu values 
    |
    |  3. Data for each R-node in MRG is maintained by a ReebGraphElement object, stored
    |     in MRGConstrLight.reebs[res_idx][rnode_index] array, where res_idx is 
    |     resolution index in MRG (res_idx=0 is the finest resolution), and rnode_index
    |     is the R-node index. 
    |
    |  4. R-edges between R-nodes at the same resolution of MRG are stored via 
    |     adjacency lists in MRGConstrLight.MRG[res_idx][rnode_index][] array, where
    |     res_idx is resolution index in MRG (res_idx=0 is the finest resolution), and
    |     rnode_index is the R-node index. 
    |     Integer rnode2_index=MRGConstrLight.MRG[res_idx][rnode_index][adj_idx] 
    |     for adj_idx > 0 indicates that R-node 
    |     MRGConstrLight.reebs[res_idx][rnode_index] is connected to R-node
    |     MRGConstrLight.reebs[res_idx][rnode2_index]. 
    |     MRGConstrLight.MRG[res_idx][rnode_index][0] stores the length of the 
    |     adjacency list for R-node at rnode_index at resolution res_idx.
    |
    |  5. R-edges that connect R-nodes from different resolutions of MRG are
    |     stored in ReebGraphElement.parents. For each R-node at resolution
    |     res_idx, ReebGraphElement.parents vector contains indices to R-nodes
    |     at a previous (finer) resolution (indexed by res_idx-1) -- i.e., 
    |     ReebGraphElement objects from MRGConstrLight.reebs[res_idx-1][] array.
    */


//Constructs MRG for a given mesh model and its mu values
public class MRGConstrLight {

    public double percentage_threshold = 0.01;
    
    public Point points[];
    public int points_length;
    public int sparse[][];
    public double mu_values[];
    
    public int all_Tsets[][][];  //used to hold all Tsets that are created
    public int tsets[][];

    public double ranges[];    
    //public Vector other_resolutions; //holds info how Tsets of diff res can be constructed
    //from all_Tsets vector
    
    public static int FINEST_RESOLUTION;
    public int MRG[][][];      //vector of Reeb Graph of different resolutions
    public ReebGraphElement reebs[][];
    public int Rnumber;
    
    public int info[];
    public Vector stacker;
    
    public int cons[][];

    public Vector doProcess(int mrg_num, Point pts[], int points_size, 
			    int sprs[][], double mu[], boolean chooser){
	
	//System.out.println("mrg_num = " + mrg_num);
	
	double mrg_num_d = (double) mrg_num;
	ecomp.epsilon2 = 1 / (mrg_num_d * mrg_num_d);

	//System.out.println("epsilon2 = " + ecomp.epsilon2);
	
	FINEST_RESOLUTION = mrg_num;
	points = pts;
	mu_values = mu;
	points_length = points_size;
	sparse = sprs;
	
	divideRanges();
	
	//System.out.println("Resampling edges that lie in different ranges...");
	doResampling();
	sparse = aman.trimSparse(sparse, points_length);
	if(checker() == false)
	    System.out.println("WARNING in MRGConstrLight: Resampling wasn't successful!");
	
	//System.out.println("Creating T-sets...");
	createTSets();
	
	//System.out.println("Creating finest resolution reeb graph...");
	createFinestResolutionReebGraph();
	
	// if chooser is TRUE, subroutine reduceNumberNodes() will be called -- not implemented, see commented-out code below
	if(chooser == true){
	    System.out.println("ERROR in MRGConstrLight: this function is not implemented (see commented-out code at the end of MRGConstrLight.java)\n Exiting...");
	    System.exit(1);
	    
	    /*System.out.println("The number of nodes before reduction: " + all_Tsets.size());
	    System.out.println("Starting reducing the number of nodes...");
	    reduceNumberNodes();
	    System.out.println("The number of nodes after reduction: " + all_Tsets.size());*/
	}
	    
	//System.out.println("Creating reeb graphs with lower resolutions...");
	createMRG();
	
	System.out.println("Done creating MRG of size: " + MRG.length);
	
	//special subroutine for testing
	//testing();
	
	Vector result = new Vector();
	result.addElement(points);
	result.addElement(new Integer(points_length));
	result.addElement(sparse);
	result.addElement(mu_values);
	result.addElement(tsets);
	result.addElement(MRG);
	result.addElement(reebs);
	
	return result;
    }
    
    public void divideRanges(){
	ranges = new double[FINEST_RESOLUTION + 1];
	
	ranges[0] = 0.0;
	
	for(int i = 1; i < FINEST_RESOLUTION; i++){
	    double i_double = (double) i;
	    double finest_double = (double) FINEST_RESOLUTION;
	    ranges[i] = i_double / finest_double;
	}
	
	ranges[FINEST_RESOLUTION] = 1.0;
	
	double temp[];
	temp = mu_values;
	
	mu_values = new double[points.length];
	for(int i = 0; i < points_length; i++){
	    mu_values[i] = temp[i];	    
	}	
    }

    public void doResampling(){
	int range1, range2;
	int index;
	for(int i = 0; i < points_length; i++){
	    
	    int adj[] = sparse[i];
	    int len = sparse[i][0];
	    range1 = whatRange(i);
	    
	    for(int j = 1; j < len; j++){
		index = adj[j];		    
		if(index > -1){
		    range2 = whatRange(index);
		    
		    if(! isInOneRange(i, index)){
			//create new vertex
			double mu1, mu2;
			double new_mu;
			int new_index;
			
			mu1 = mu_values[i];
			mu2 = mu_values[index];
			
			//calculate new_mu
			if(range1 > range2)
			    new_mu = ranges[range2 + 1];
			else
			    new_mu = ranges[range1 + 1];
			
			createNewVertex(i, index, mu1, mu2, new_mu);
		    }
		}
	    }
	}
    }
    
    public boolean checker(){
	for(int i = 0; i < points_length; i++){
	    int range1 = whatRange(i);
	    int adj[] = sparse[i];
	    for(int j = 1; j < adj[0]; j++){
		int range2 = whatRange(adj[j]);
		
		if(! isInOneRange(i, adj[j])){
		    return false;
		}
	    }
	}
	return true;
    }
    
    public boolean isInOneRange(int index1, int index2){
	double muV1 = mu_values[index1];
	double muV2 = mu_values[index2];
	
	for(int i = 0; i < ranges.length - 1; i++){
	    double left_range = ranges[i];
	    double right_range = ranges[i+1];
	    
	    if(left_range <= muV1 && muV1 <= right_range){
		if(left_range <= muV2 && muV2 <= right_range){
		    return true;
		}
	    }	    
	}
	return false;
    }
    

    
    public void createTSets(){
	int range;
	
	all_Tsets = new int[FINEST_RESOLUTION][][];
	
	info = new int[points_length];
	
	int counter = 0;
	for(int i = 0; i < info.length; i++){
	
	    range = whatRange(i);	    
	    double rng = ranges[range];
	    double muV = mu_values[i];
	    
	    if(rng == muV && range != 0){
		info[i] = 2;
		counter++;
	    }
	    else
		info[i] = 1;
	}
	
	//System.out.print("Number of vertices that lie exactly on the border: ");
	//System.out.println(counter);
	
	//holds lengths of one_range_tsets
	int lens[] = new int[FINEST_RESOLUTION];
	for(int i = 0; i < FINEST_RESOLUTION; i++){
	    lens[i] = 0;
	    all_Tsets[i] = new int[aman.size_inc_tsets][];
	}
	
	int i = 0;
	while(i < info.length){
	    if(info[i] == 0)
		i++;
	    else if(info[i] == 1){
		int rang = whatRange(i);
		int one_set[] = createOneTset(i, rang);
		
		int one_range_tsets[][] = all_Tsets[rang];
		
		if(lens[rang] >= one_range_tsets.length)
		    one_range_tsets = aman.expand_tsets(one_range_tsets, lens[rang]);
		
		one_range_tsets[lens[rang]] = one_set;
		lens[rang]++;
		
		all_Tsets[rang] = one_range_tsets;
	    }
	    else if(info[i] == 2 || info[i] == 3){
		int rang = whatRange(i) - 1;
		int one_set[] = createOneTset(i, rang);
		
		int one_range_tsets[][] = all_Tsets[rang];
		
		if(lens[rang] >= one_range_tsets.length){
		    one_range_tsets = aman.expand_tsets(one_range_tsets, lens[rang]);
		}
		
		one_range_tsets[lens[rang]] = one_set;
		lens[rang]++;
		
		all_Tsets[rang] = one_range_tsets;
	    }
	}
	
	all_Tsets = aman.trimTsets(all_Tsets, lens);
	
	//for(i = 0; i < all_Tsets.length; i++){
	    //int tmp[][] = all_Tsets[i];
	    //System.out.println("The number of nodes in the range " + i + " : " + tmp.length);
	    //System.out.println("*******************************************");
	    
	    // for(int j = 0; j < tmp.length; j++){
	    // 	int tmp2[] = tmp[j];
	    // 	System.out.print(tmp2.length + " ");
	    // }
	    //System.out.println(" ");
	    //System.out.println("*******************************************");
	//}	
    }
        
    public int[] createOneTset(int point_index, int range){
	
	int tset[] = new int[aman.size_inc];
	int tset_length = 0;
	
	stacker = new Vector();
	
	int real_range;
	Integer temp;
	int index = point_index;
	
	if(info[point_index] > 0){
	    stacker.addElement(new Integer(index));
	}
	    	
	while(stacker.size() != 0){
	    
	    temp = (Integer) stacker.remove(0);
	    index = temp.intValue();
	    
	    if(createOnePoint(index,  range) == true){
	    
		if(tset_length >= tset.length)
		    tset = aman.expand(tset, tset_length);
		
		tset[tset_length] = index;
		tset_length++;
	    }
	}
	
	tset = aman.trimOneTset(tset, tset_length);
	return tset;
    }
        
    public boolean createOnePoint(int point_index, int range){
	
	boolean result = false;
	int real_range = whatRange(point_index);
	
	if(info[point_index] > 0){
	    result = true;
	    
	    if(info[point_index] == 1){
		info[point_index] = 0;
		
		if(real_range != range)
		    System.out.println("WARNING in MRGConstrLight: error code 0");	
	    }
	    
	    else if(info[point_index] == 2){
		
		if(real_range == range)
		    info[point_index] = 3;
		else if(real_range == range + 1)
		    info[point_index] = 1;
		else
		    System.out.println("WARNING in MRGConstrLight: error code 1");	
	    }
	    
	    else if(info[point_index] == 3){
		info[point_index] = 0;
		
		if(real_range != range + 1)
		    System.out.println("WARNING in MRGConstrLight: error code 2");	
	    }
	    
	    int adj[] = sparse[point_index];
	    int len = adj[0];
	    for(int i = 1; i < len; i++){
		int index = adj[i];
		if(info[index] > 0){
		    real_range = whatRange(index);
		    
		    if(info[index] == 1 && real_range == range){
			if(! stacker.contains(new Integer(index)))
			    stacker.addElement(new Integer(index));
		    }
		    
		    else if(info[index] == 3 && real_range == range + 1){
			if(! stacker.contains(new Integer(index)))
			    stacker.addElement(new Integer(index));
		    }
		    
		    else if(info[index] == 2){
			
			if(real_range == range + 1 || real_range == range){
			    if(! stacker.contains(new Integer(index)))
				stacker.addElement(new Integer(index));
			}
		    }
		}
	    }
	}

	return result;	
    }
        
    public void createFinestResolutionReebGraph(){
	int MRGsize = 0;
	int tm = ranges.length - 1;
	while(tm != 1){
	    MRGsize++;
	    tm = tm / 2;
	}
	MRGsize++;
		
	int R_number = 0;
	//calculate how many R-nodes (T-sets) we have
	for(int i = 0; i < all_Tsets.length; i++){
	    int temp[][] = all_Tsets[i];
	    R_number = R_number + temp.length;
	}
	
	Rnumber = R_number;
	
	//System.out.println("number of nodes: " + Rnumber);
	
	MRG = new int[MRGsize][][];
	MRG[0] = new int[Rnumber][aman.size_inc_sparse];
	
	for(int i = 0; i < Rnumber; i++){
	    MRG[0][i][0] = 1;
	}
	
	reebs = new ReebGraphElement[MRGsize][];
	reebs[0] = new ReebGraphElement[Rnumber];
	
	for(int i = 0; i < R_number; i++){
	    ReebGraphElement rel = new ReebGraphElement();
	    rel.index = i;
	    
	    int temp_int = calculateRange(i);
	    rel.left_bound = ranges[temp_int];
	    rel.right_bound = ranges[temp_int + 1];
	    
	    rel.parents = null;

	    rel.Tsets = new Vector();
	    rel.Tsets.addElement(new Integer(i));
	    
	    reebs[0][i] = rel;
	}
	
	//System.out.println("Number of intervals: " + all_Tsets.length);
	
	int index1 = 0;
	int index2 = 0;
	int main_index = 0;
	int count = 0;
	
	for(int i = 0; i < all_Tsets.length - 1; i++){
	    int current[][] = all_Tsets[i];
	    int following[][] = all_Tsets[i + 1];
	    
	    //System.out.println("Working on interval: " + i);
	    
	    for(int j = 0; j < current.length; j++){
		int Tset1[] = current[j];
		
		for(int k = 0; k < following.length; k++){
		    int Tset2[] = following[k];
		    
		    if(isConnectedToTset(Tset1, Tset2)){
			index2 = main_index + current.length + k;
			index1 = main_index + j;
			
			int len1 = MRG[0][index1][0];
			int len2 = MRG[0][index2][0];
			
			int adj1[] = MRG[0][index1];
			int adj2[] = MRG[0][index2];
			
			if(len1 >= MRG[0][index1].length)
			    MRG[0][index1] = aman.expand_sparse(adj1, len1);
			    
			MRG[0][index1][len1] = index2;
			MRG[0][index1][0]++;
			
			if(len2 >= MRG[0][index2].length)
			    MRG[0][index2] = aman.expand_sparse(adj2, len2);
			
			MRG[0][index2][len2] = index1;
			MRG[0][index2][0]++;
		    }
		}
	    }
	    
	    main_index = main_index + current.length;
	}
		
	//puts all Tsets into a single vector
	putTsetsInVector();
    }
    
    public int calculateRange(int index){
	int temp = 0;
	
	for(int i = 0; i < all_Tsets.length; i++){
	    int vec[][] = all_Tsets[i];
	    temp = temp + vec.length;
	    
	    if(temp > index){
		return i;		    
	    }
	}
	
	System.out.println("WARNING in MRGConstrLight: problem with calculating range. Returning -1");
	return -1;
    }
        
    public void putTsetsInVector(){
	int tcount = 0;
	for(int i = 0; i < all_Tsets.length; i++){
	    int temp[][] = all_Tsets[i];
	    tcount = tcount + temp.length;
	}
	
	tsets = new int[tcount][];
	int k = 0;
	for(int i = 0; i < all_Tsets.length; i++){
	    int temp[][] = all_Tsets[i];
	    
	    for(int j = 0; j < temp.length; j++){
		tsets[k] = temp[j];
		k++;
	    }
	}
    }
    
    
    
    public void createMRG(){
	
	//Main loop that creates all other resolutions
	int prev_MRG_count = 0;
	int curr_MRG_count;
	ReebGraphElement element, copy_element;
	while(ranges.length > 2){    
	    
	    curr_MRG_count = prev_MRG_count+1;
	    //create a clone of the RG
	    reebs[curr_MRG_count] = new ReebGraphElement[reebs[prev_MRG_count].length];
	    
	    for(int i = 0; i < reebs[curr_MRG_count].length; i++){
		element = reebs[prev_MRG_count][i];
		copy_element = new ReebGraphElement();
		
		copy_element.index = element.index;
		copy_element.left_bound = element.left_bound;
		copy_element.right_bound = element.right_bound;
		
		copy_element.Tsets = (Vector) element.Tsets.clone();
		
		copy_element.parents = new Vector();
		copy_element.parents.addElement(new Integer(element.index));
		
		reebs[curr_MRG_count][i] = copy_element;
	    }
	    
	    //holds info on which nodes to be unified
	    cons = new int[reebs[prev_MRG_count].length][];
	    for(int i = 0; i < cons.length; i++){
		int vec[] = MRG[prev_MRG_count][i];
		cons[i] = new int[vec.length];
		
		for(int j = 0; j < vec[0]; j++){
		    cons[i][j] = vec[j];		    
		}	
	    }
	    	    
	    for(int i = 0; i < ranges.length - 2; i = i + 2){
		unifyTwoRanges(i, i + 1, i + 2, prev_MRG_count);
	    }
	    
	    //removes every other element from ranges
	    double temp_ranges[];
	    temp_ranges = ranges;
	    ranges = new double[(temp_ranges.length/2) + 1];
	    
	    int j = 0;
	    for(int i = 0; i < temp_ranges.length; i = i + 2){
		ranges[j] = temp_ranges[i];
		j++;
	    }
	    ranges[0] = temp_ranges[0];
	    ranges[ranges.length-1] = temp_ranges[temp_ranges.length-1];
	    
	    updateReebs(prev_MRG_count);
	    updateMRG(prev_MRG_count);
	    
	    int counter = 0;
	    for(int i = 0; i < MRG[curr_MRG_count].length; i++){
		ReebGraphElement el1 = reebs[curr_MRG_count][i];
		int adj[] = MRG[curr_MRG_count][i];
		
		for(j = 1; j < adj[0]; j++){
		    int index = adj[j];
		    ReebGraphElement el2 = reebs[curr_MRG_count][index];
		    
		    if(el1.left_bound == el2.left_bound && 
		       el1.right_bound == el2.right_bound){
			counter++;
		    }	    
		}
	    }
	    
	    //System.out.println("we have " + counter + " pairs of nodes that are connected and lie in one interval!");
	    
	    prev_MRG_count++;
	    
	    //System.out.println("In resolution " + (prev_MRG_count+1) + " there are " + MRG[prev_MRG_count].length + " nodes!");
	}
    }
    
    public void updateReebs(int prev_MRG_index){
	int curr_MRG_index = prev_MRG_index+1;
	ReebGraphElement temp_reeb[];
	temp_reeb = reebs[curr_MRG_index];
	int null_count = 0;
	for(int i = 0; i < temp_reeb.length; i++){
	    if(temp_reeb[i] == null){
		null_count++;
	    }
	}
	
	reebs[curr_MRG_index] = new ReebGraphElement[temp_reeb.length - null_count];
	null_count = 0;
	for(int i = 0; i < temp_reeb.length; i++){
	    if(temp_reeb[i] != null){
		temp_reeb[i].index = null_count;
		reebs[curr_MRG_index][null_count] = temp_reeb[i];
		null_count++;
	    }
	}
    }
    
    public void updateMRG(int prev_MRG_index){
	int curr_MRG_index = prev_MRG_index+1;
	int size_reeb = reebs[curr_MRG_index].length;
	int parents[][];
	parents = new int[size_reeb][];
	
	for(int i = 0; i < size_reeb; i++){
	    Vector pars = reebs[curr_MRG_index][i].parents;
	    int temp = pars.size();
	    parents[i] = new int[temp];
	    
	    for(int j = 0; j < temp; j++){
		Integer temp_int = (Integer) pars.elementAt(j);
		parents[i][j] = temp_int.intValue();
	    }
	}
	
	MRG[curr_MRG_index] = new int[size_reeb][aman.size_inc_sparse];
	
	for(int i = 0; i < size_reeb; i++)
	    MRG[curr_MRG_index][i][0] = 1;
	
	for(int i = 0; i < size_reeb; i++){
	    
	    int j = i+1;
	    while(j < size_reeb){
		
		for(int k1 = 0; k1 < parents[i].length; k1++){
		    for(int k2 = 0; k2 < parents[j].length; k2++){
			if(areNodesConnected(prev_MRG_index, parents[i][k1], parents[j][k2]) == true &&
			   areNodesConnected(curr_MRG_index, i, j) == false){
			    
			    int index1 = i;
			    int index2 = j;
			    
			    int len1 = MRG[curr_MRG_index][index1][0];
			    int len2 = MRG[curr_MRG_index][index2][0];
			    
			    int adj1[] = MRG[curr_MRG_index][index1];
			    int adj2[] = MRG[curr_MRG_index][index2];
			    
			    if(len1 >= MRG[curr_MRG_index][index1].length)
				MRG[curr_MRG_index][index1] = aman.expand_sparse(adj1, len1);
			    
			    MRG[curr_MRG_index][index1][len1] = index2;
			    MRG[curr_MRG_index][index1][0]++;
			    
			    if(len2 >= MRG[curr_MRG_index][index2].length)
				MRG[curr_MRG_index][index2] = aman.expand_sparse(adj2, len2);
			    
			    MRG[curr_MRG_index][index2][len2] = index1;
			    MRG[curr_MRG_index][index2][0]++;
			}
		    }
		}
				
		j++;
	    }
	}
    }
    
    public boolean areNodesConnected(int mrg_index, int index1, int index2){
	int adj[] = MRG[mrg_index][index1];
	int len = adj[0];
	
	for(int i = 1; i < len; i++){
	    
	    if(adj[i] == index2)
		return true;	    
	}
	
	return false;	
    }

    public void unifyTwoRanges(int left_range_index, int mid_range_index, int right_range_index, int prev_MRG_index){
	double left_range, mid_range, right_range;
	ReebGraphElement element = null;
	ReebGraphElement element2 = null;
	int curr_MRG_index = prev_MRG_index+1;
	
	left_range = ranges[left_range_index];
	mid_range = ranges[mid_range_index];
	right_range = ranges[right_range_index];
	
	int i = 0;
	while(i < reebs[curr_MRG_index].length){
	    if(reebs[curr_MRG_index][i] != null){
		element = reebs[curr_MRG_index][i];
		if((element.left_bound == left_range && element.right_bound == mid_range) 
		   || (element.left_bound == mid_range && element.right_bound == right_range)
		   || (element.left_bound == left_range && element.right_bound == right_range)){
		    
		    int j = 1;
		    int temp_i = i;
		    while(j < cons[element.index][0]){
			int index = cons[element.index][j];
			if(reebs[curr_MRG_index][index] != null){
			    element2 = reebs[curr_MRG_index][index];
			    
			    if((element2.left_bound == left_range && element2.right_bound == mid_range) 
			       || (element2.left_bound == mid_range && element2.right_bound == right_range)
			       || (element2.left_bound == left_range && element2.right_bound == right_range)){
				
				unifyTwoNodes(element.index, element2.index, left_range, right_range, prev_MRG_index);
				temp_i = -1;
			    }
			}
			j++;
		    }
		    i = temp_i;
		}
	    }
	    i++;
	}
	
	for(i = 0; i < reebs[curr_MRG_index].length; i++){
	    if(reebs[curr_MRG_index][i] != null){
		ReebGraphElement elt = reebs[curr_MRG_index][i];
		
		if((elt.left_bound == left_range && elt.right_bound == mid_range)
		   || (elt.left_bound == mid_range && elt.right_bound == right_range)){
		    elt.left_bound = left_range;
		    elt.right_bound = right_range;
		    
		    reebs[curr_MRG_index][i] = elt;
		}
	    }
	}
	
    }
    
    public void unifyTwoNodes(int node_index1, int node_index2, double left_b, double right_b, int prev_MRG_index){
	int curr_MRG_index = prev_MRG_index + 1;
	int bigger, smaller;
	
	if(node_index1 == node_index2){
	    System.out.println("WARNING in MRGConstrLight: node_index1 and node_index2 are the same!");
	    return;
	}
	
	ReebGraphElement el1 = reebs[curr_MRG_index][node_index1];
	ReebGraphElement el2 = reebs[curr_MRG_index][node_index2];
	
	for(int i = 0; i < el2.Tsets.size(); i++){
	    Integer temp_integer = (Integer) el2.Tsets.elementAt(i);
	    if(! el1.Tsets.contains(temp_integer))
		el1.Tsets.addElement(temp_integer);
	}
	
	for(int i = 0; i < el2.parents.size(); i++){
	    Integer int_t = (Integer) el2.parents.elementAt(i);
	    
	    if(! el1.parents.contains(int_t)){
		el1.parents.addElement(int_t);
	    }
	}
	
	el1.left_bound = left_b;
	el1.right_bound = right_b;
		
	el2.left_bound = left_b;
        el2.right_bound = right_b;
	
	el2.index = -1;
	
	reebs[curr_MRG_index][node_index2] = null;
	
	int adj1[];
	int len1;
	int adj2[] = cons[node_index2];
	int len2 = adj2[0];
	int index1, index2;
	boolean found;
	
	for(int i = 1; i < len2; i++){
	    index2 = adj2[i];
	    
	    adj1 = cons[node_index1];
	    len1 = adj1[0];
	    
	    found = false;
	    for(int j = 1; j < len1; j++){
		index1 = adj1[j];
		if(index2 == index1)
		    found = true;
	    }
	    
	    if(found == false && index2 != node_index1){
		
		if(len1 >= adj1.length)
		    cons[node_index1] = aman.expand_sparse(adj1, len1);
		
		cons[node_index1][len1] = index2;
		cons[node_index1][0]++;		 
	    }	    
	}
    }

    public boolean isConnectedToTset(int TSet1[], int TSet2[]){
	for(int i = 0; i < TSet1.length; i++){
	    int temp = TSet1[i];
	    if(isConnectedToTSet(temp, TSet2))
	       return true;
	}
		    
	return false;
    }
    
    public boolean isConnectedToTSet(int ptr_index, int TSet[]){
	if(TSet == null || TSet.length == 0)
	    return true;
	
	int adj[] = sparse[ptr_index];
	int len = adj[0];
	
	for(int i = 0; i < TSet.length; i++){
	    int index = TSet[i];
	    
	    for(int j = 1; j < len; j++){
		if(adj[j] == index)
		    return true;		
	    }
	}
	
	return false;
    }
    
    public int whatRange(int ptr_index){
	double muV = mu_values[ptr_index];
	
	if(Math.abs(muV - 1.0) < ecomp.epsilon2){
	    return (ranges.length - 2);
	}
	
	for(int i = 0; i < ranges.length - 1; i++){
	    double left_range = ranges[i];
	    double right_range = ranges[i+1];
	    
	    if(left_range <= muV && muV < right_range){
		return i;		
	    }	    
	}
	
	System.out.println("WARNING in MRGConstrLight: algorithm in whatRange is wrong! Returning -1");

	//System.out.println("number of points: " + points_length);
	//System.out.println("point: " + ptr_index);
	//System.out.println("mu value: " + muV);
	//System.out.println("epsilon: " + ecomp.epsilon2);
	
	return -1;
    }
    
    //creates new vertex with the value new_mu, puts it in the points vector and mu_values
    public void createNewVertex(int point1, int point2, double mu1, double mu2, double new_mu){
	int i = point1;
	int index = point2;
	int adj[] = sparse[i];
	int adj2[] = sparse[index];
	
	//System.out.println("p1 = " + point1 + ", p2 = " + point2 + ", mu1 = " + mu1 + ", mu2 = " + mu2 + ", new_mu = " + new_mu);
	
	Point P1 = points[point1];
	Point P2 = points[point2];

	Point P = new Point();
	
	if(mu1 > mu2){
	    P.X = (P1.X*(new_mu - mu2) + P2.X*(mu1 - new_mu)) / (mu1 - mu2);
	    P.Y = (P1.Y*(new_mu - mu2) + P2.Y*(mu1 - new_mu)) / (mu1 - mu2);
	    P.Z = (P1.Z*(new_mu - mu2) + P2.Z*(mu1 - new_mu)) / (mu1 - mu2);
	}
	else{
	    P.X = (P1.X*(mu2 - new_mu) + P2.X*(new_mu - mu1)) / (mu2 - mu1);
	    P.Y = (P1.Y*(mu2 - new_mu) + P2.Y*(new_mu - mu1)) / (mu2 - mu1);
	    P.Z = (P1.Z*(mu2 - new_mu) + P2.Z*(new_mu - mu1)) / (mu2 - mu1);
	}
	
	if(points_length >= points.length){
	    points = aman.expand(points, points_length);
	    mu_values = aman.expand(mu_values, points_length);
	    sparse = aman.expand_sparse(sparse, points_length);
	}
	
	points[points_length] = P;
	mu_values[points_length] = new_mu;
	
	points_length++;
		    
	//removes old connections
	int k;
	int len;
	boolean stop;
		
	k = 1;
	len = sparse[i][0];
	stop = false;
	while(stop == false){
	    if(sparse[i][k] == index){
		sparse[i][k] = -1;
		stop = true;
	    }
	    k++;
	}
	
	k = 1;
	len = sparse[index][0];
	stop = false;
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
	
	//looks for extra connections
	adj = sparse[i];
	adj2 = sparse[index];
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

/*public void testing(){
	for(int i = 0; i < MRG.size(); i++){
	    Vector reeb_graph = (Vector) MRG.elementAt(i);
	    
	    System.out.print("reeb_graph, resolution: ");
	    System.out.println(i + 1);
	    
	    for(int j = 0; j < reeb_graph.size(); j++){
		Vector reeb_elements = (Vector) reeb_graph.elementAt(j);
		
		ReebGraphElement temp_el = (ReebGraphElement) reeb_elements.elementAt(0);
				
		System.out.print(temp_el.index);
		
		System.out.print(":[");
		System.out.print(temp_el.left_bound);
		System.out.print(", ");
		System.out.print(temp_el.right_bound);
		System.out.print("]");

		System.out.print(" (");
		for(int k = 0; k < temp_el.Tsets.size(); k++){
		    Integer temp_int = (Integer) temp_el.Tsets.elementAt(k);

		    System.out.print(temp_int.intValue());
		    System.out.print(", ");
		}
		System.out.print(")");
		System.out.print("  ");
		
		for(int l = 1; l < reeb_elements.size(); l++){
		    temp_el = (ReebGraphElement) reeb_elements.elementAt(l);
		    System.out.print(temp_el.index);
		    System.out.print("  ");		    
		}
		
		System.out.println("");
	    }
	    
	    System.out.println("");
	}
	
	System.out.println("The Tsets are: ");
	
	for(int i = 0; i < all_Tsets.size(); i++){
	    Vector temp_vec = (Vector) all_Tsets.elementAt(i);
	    
	    for(int j = 0; j < temp_vec.size(); j++){
		Integer temp_i = (Integer) temp_vec.elementAt(j);
		
		System.out.print(temp_i.intValue());
		System.out.print("  ");
	    }
	    
	    System.out.println("");
	}
	}*/


    
/*public void reduceNumberNodes(){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
	int el_degree;
	int i = 0;
	boolean checker = false;
	
	//for(int i = 0; i < reeb_graph.size(); i++){
	while(i < reeb_graph.size()){
	    if(isNodeNotImportant(i)){
		el_degree = connectionDegree(i);
		
		Vector adj_vec = (Vector) reeb_graph.elementAt(i);
		ReebGraphElement element = (ReebGraphElement) adj_vec.elementAt(0);	    
		
		checker = false;
		
		if(el_degree == 1){
		    checker = lookForNodesToReduce(element, 1);
		}
		else if(el_degree == 2){
		    checker = lookForNodesToReduce(element, 2);
		}
				
		if(checker == false)
		    i++;
	    }
	    else
		i++;
	}
    }   
    
    public boolean lookForNodesToReduce(ReebGraphElement element, int number_connection){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
		
	for(int i = 0; i < reeb_graph.size(); i++){
	    if(isNodeNotImportant(i)){
		Vector adj_vec = (Vector) reeb_graph.elementAt(i);
		ReebGraphElement element2 = (ReebGraphElement) adj_vec.elementAt(0);
		if(ifQualified(element, element2, number_connection)){
		    reduceConnection(element, element2);
		    return true;
		}		
	    }
	    }
	
	return false;
    }
    
    public void reduceConnection(ReebGraphElement element1, ReebGraphElement element2){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
	Vector adj_vec1 = (Vector) reeb_graph.elementAt(element1.index);
	//ReebGraphElement element1 = (ReebGraphElement) adj_vec1.elementAt(0);
	Vector adj_vec2 = (Vector) reeb_graph.elementAt(element2.index);
	//ReebGraphElement element2 = (ReebGraphElement) adj_vec2.elementAt(0);
	
	if(adj_vec2.size() == 1){
	    copyTsets(element1.index, element2.index);
	    reeb_graph.removeElementAt(element2.index);
	    all_Tsets.removeElementAt(element2.index);
	    
	    reduceIndecies(element2.index);
	    all_Tsets.setElementAt(element1.Tsets, element1.index);
	}
	else if(adj_vec2.size() == 2){
	    ReebGraphElement element3 = (ReebGraphElement) adj_vec2.elementAt(1);
	    Vector adj = (Vector) reeb_graph.elementAt(element3.index);
	    for(int i = 1; i < adj.size(); i++){
		ReebGraphElement element4 = (ReebGraphElement) adj.elementAt(i);
		
		if(element4.index == element2.index){
		    adj.removeElementAt(i);
		    break;
		}
	    }
	    
	    copyTsets(element1.index, element2.index);
	    reeb_graph.removeElementAt(element2.index);
	    all_Tsets.removeElementAt(element2.index);
	    
	    reduceIndecies(element2.index);
	    all_Tsets.setElementAt(element1.Tsets, element1.index);
	}
	else
	    System.out.println("size of adj_vec2 has to be 1 or 2!");
    }
    
    public void reduceIndecies(int index){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
	for(int i = 0; i < reeb_graph.size(); i++){
	    Vector adj = (Vector) reeb_graph.elementAt(i);
	    ReebGraphElement el = (ReebGraphElement) adj.elementAt(0);
	    
	    if(el.index > index){
		el.index--;		
	    }	    
	}	
    }
    
    public void copyTsets(int index1, int index2){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
	Vector adj_vec1 = (Vector) reeb_graph.elementAt(index1);
	ReebGraphElement element1 = (ReebGraphElement) adj_vec1.elementAt(0);
	Vector adj_vec2 = (Vector) reeb_graph.elementAt(index2);
	ReebGraphElement element2 = (ReebGraphElement) adj_vec2.elementAt(0);
	
	for(int i = 0; i < element2.Tsets.size(); i++){
	    Integer temp = (Integer) element2.Tsets.elementAt(i);
	    
	    if(! element1.Tsets.contains(temp)){
		element1.Tsets.addElement(temp);		
	    }
	}	
    }
    
    public boolean ifQualified(ReebGraphElement element1, ReebGraphElement element2, int number_connection){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
	Vector adj_vec1 = (Vector) reeb_graph.elementAt(element1.index);
	//ReebGraphElement element1 = (ReebGraphElement) adj_vec1.elementAt(0);
	Vector adj_vec2 = (Vector) reeb_graph.elementAt(element2.index);
	//ReebGraphElement element2 = (ReebGraphElement) adj_vec2.elementAt(0);
	
	if(number_connection == 1){
	    if(element1.left_bound == element2.left_bound
	       && element1.right_bound == element2.right_bound
	       && element1.index != element2.index
	       && adj_vec1.size() == 1 && adj_vec2.size() == 1){
		return true;		
	    }
	    
	    return false;
	}
	else if(number_connection == 2){
	    double lb1, lb2, rb1, rb2;
	    int i1, i2;
	    	    
	    if(element1.left_bound == element2.left_bound
	       && element1.right_bound == element2.right_bound
	       && element1.index != element2.index
	       && adj_vec1.size() == 2 && adj_vec2.size() == 2){
		
		ReebGraphElement element3 = (ReebGraphElement) adj_vec1.elementAt(1);
		ReebGraphElement element4 = (ReebGraphElement) adj_vec2.elementAt(1);
		
		if(element3.left_bound == element4.left_bound
		   && element3.right_bound == element4.right_bound
		   && element3.index == element4.index){
		    
		    return true;
		}		
	    }
	    return false;
	}
	
	else
	    System.out.println("connection_number has to be 1 or 2!");
	
	return false;
    }
    
    public boolean isNodeNotImportant(int index){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
	Vector adj_vec = (Vector) reeb_graph.elementAt(index);
	ReebGraphElement element = (ReebGraphElement) adj_vec.elementAt(0);
	
	double size_tsets = (double) element.Tsets.size();
	double all_points = (double) points.size();
	
	if((size_tsets / all_points) <= percentage_threshold)
	    return true;
	else
	    return false;
    }
    
    public int connectionDegree(int index){
	Vector reeb_graph = (Vector) MRG.elementAt(0);
	Vector adj_vec = (Vector) reeb_graph.elementAt(index);
	ReebGraphElement element = (ReebGraphElement) adj_vec.elementAt(0);
	if(adj_vec.size() == 1)
	    return 1;
	else if(adj_vec.size() == 2)
	    return 2;
	else
	    return 0;
	    }    
*/    

}
