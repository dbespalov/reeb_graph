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

//contains properties of a 3D point

public class Point {

    public double X,Y,Z;
    public int TriangleLocation;       //number triangle the point is in
    public static double epsilon2 = .000001;

    public boolean IsEqualTo(Point P) {
	if (Math.abs(X-P.X)<epsilon2) {
            if (Math.abs(Y- P.Y)<epsilon2) {
    	        if (Math.abs(Z-P.Z)<epsilon2) {
		    return(true);
		}
	    }
	}
	return(false);
    }

    //returns the distance between 2 points in 3D space
    public double GetDistance(Point P) {
	
        double Distance = Math.sqrt(Math.pow(X-P.X, 2) +
       				    Math.pow(Y-P.Y, 2) +
       				    Math.pow(Z-P.Z, 2));
        return(Distance);
    }

	public boolean IsBetween(Point A, Point B) {
	return(((X <= A.X && X >= B.X)||(X >= A.X && X <= B.X)) &&
               ((Y <= A.Y && Y >= B.Y)||(Y >= A.Y && Y <= B.Y)) &&
               ((Z <= A.Z && Z >= B.Z)||(Z >= A.Z && Z <= B.Z)) &&
                !IsEqualTo(A) && !IsEqualTo(B));
    }

    //returns true if a 3D point P (this.X, this.Y, this.Z) 
    // is in between the points A and B with a tolerance of epsilon2
    public boolean IsBetween2(Point A, Point B) {
	double GX,LX,GY,LY,GZ,LZ;
        if (A.X > B.X) {
            GX = A.X;
            LX = B.X;
        }
	else {
            GX = B.X;
            LX = A.X;
        }
        if (A.Y > B.Y) {
            GY = A.Y;
            LY = B.Y;
        }
        else {
            GY = B.Y;
            LY = A.Y;
        }
        if (A.Z > B.Z) {
            GZ = A.Z;
            LZ = B.Z;
        }
        else {
            GZ = B.Z;
            LZ = A.Z;
        }
        return(((X > LX || Math.abs(X-LX)<epsilon2) &&
                (X < GX || Math.abs(X-GX)<epsilon2)) &&
               ((Y > LY || Math.abs(Y-LY)<epsilon2) &&
                (Y < GY || Math.abs(Y-GY)<epsilon2)) &&
               ((Z > LZ || Math.abs(Z-LZ)<epsilon2) &&
                (Z < GZ || Math.abs(Z-GZ)<epsilon2)));
    } 


    // Computes the dot product between 
    //        the coplanar normal to a side of a triangle (AB) 
    //   and 
    //        vector AP in the plane of the triangle

    // Assumes that A and B are vertices of some triangle, 
    // and point P (this.X, this.Y, this.Z) is in the plane of that triangle.
    public double dotProduct(Point A, Point B, int which) {
	if (which == 1) {
            return((A.Z-B.Z)*(Y-A.Y) + (B.Y-A.Y)*(Z-A.Z));
        }
        else if (which == 2) {
            return((A.Z-B.Z)*(X-A.X) + (B.X-A.X)*(Z-A.Z));
        }
        else {
            return((A.Y-B.Y)*(X-A.X) + (B.X-A.X)*(Y-A.Y));
        }
    }

}
