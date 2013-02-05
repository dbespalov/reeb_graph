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

//class for epsilon-based floating num. comparisons

public class ecomp{
    
    public static boolean eq(double r, double q){
	if(Math.abs(r - q) < epsilon)
	    return true;
	else
	    return false;
    }
    
    public static boolean gr(double r, double q){
	if((r - epsilon) > q)
	    return true;
	else
	    return false;
    }
    
    public static boolean le(double r, double q){
	if((r + epsilon) < q)
	    return true;
	else
	    return false;
    }
    
    public static boolean ge(double r, double q){
	if((r + epsilon) > q)
	    return true;
	else
	    return false;
    }
    
    public static boolean lq(double r, double q){
	if((r - epsilon) < q)
	    return true;
	else
	    return false;
    }
    
    public static double epsilon;
    public static double epsilon2;
}
