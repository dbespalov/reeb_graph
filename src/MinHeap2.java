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


//Min-Heap implementation

import java.util.*;
import java.io.*;
import java.util.Vector.*;
import java.lang.Math.*;
import java.lang.Object.*;

public class MinHeap2 {
    
    public int iconverter[];
    public MinHeapElement heap[];
    public int heap_size;

    
    public MinHeap2(MinHeapElement arr[], int length){
	iconverter = new int[length];
	
	heap = arr;
	heap[0] = null;
	heap_size = length;
	
	for(int i = 1; i <= heap_size; i++){
	    iconverter[heap[i].index] = i;
	}
	
	buildMinHeap();
    }
    
    public int parent(int index){
	return (index/2);
    }
    
    public int left(int index){
	return (2*index);
    }
    
    public int right(int index){
	return (2*index + 1);
    }
    
    public void minHeapify(int index){
	int temp = min_heapify(index);

	while(temp != -1)
	    temp = min_heapify(temp);
    }
    
    public int min_heapify(int index){
	int l, r, smallest;
	MinHeapElement temp;
	
	l = 2*index;
	r = 2*index + 1;
	
	if(l <= heap_size && heap[l].key < heap[index].key)
	    smallest = l;
	else
	    smallest = index;
	
	if(r <= heap_size && heap[r].key < heap[smallest].key)
	    smallest = r;
	
	if(smallest != index){
	    MinHeapElement el1 = heap[index];
	    MinHeapElement el2 = heap[smallest];
	    
	    heap[index] = el2;
	    heap[smallest] = el1;
	    
	    iconverter[el2.index] = index;
	    iconverter[el1.index] = smallest;
	    	    
	    return smallest;
	}
	
	return -1;
    }
    
    public void buildMinHeap(){
	
	for(int i = (heap_size / 2); i >= 1; i--){
	    minHeapify(i);	    
	}
    }
    
    public MinHeapElement extractMin(){
	if(heap_size < 1){
	    System.out.println("WARNING in MinHeap2: min-heap is empty! Returning null...");
	    return null;
	}
	
	MinHeapElement min = heap[1];
	
	iconverter[heap[1].index] = -1;
	heap[1] = heap[heap_size];
	heap_size--;
	iconverter[heap[1].index] = 1;
	
	minHeapify(1);
	
	return min;
    }
    
    public void decreaseKey(int index_in_points, double key){
	int index = iconverter[index_in_points];
	
	if(key >= heap[index].key){
	    System.out.println("WARNING in MinHeap2: new key is bigger than the current key! Ignoring...");
	    return;
	}
	
	heap[index].key = key;
	
	int i = index;
	int parent_i = index/2;
	while(i > 1 && heap[parent_i].key > heap[i].key){
	    MinHeapElement el1 = heap[i];
	    MinHeapElement el2 = heap[parent_i];
	    heap[i] = el2;
	    heap[parent_i] = el1;
	    
	    iconverter[el2.index] = i;
	    iconverter[el1.index] = parent_i;

	    i = parent_i;
	    parent_i = i/2;
	}	
    }
    
}
