package com.graduate.model;

import com.graduate.util.MatrixOperator;

public class Camera {

    private double[][] intrinsicParameterMatrix;
    private double[][] extrinsicParameterMatrix;
    private double[][] rotateMatrix;
    private double[][] position;

    //    private double focalLength;         //   Meters.

//    private double fieldOfView;         //   Degrees.
//
//    private int resX;                   //   Horizontal resolution in pixels.
//    private int resY;                   //   Vertical resolution in pixels.
//
//    private double dX;                  //   Center-to-center horizontal pixel spacing.
//    private double dY;                  //   Center-to-center vertical pixel spacing.
//
//    private int cX;                     //   X coordinate of center pixel.
//    private int cY;                     //   Y coordinate of center pixel.

    public Camera(double[][]intrinsicParameterMatrix, double[][]extrinsicParameterMatrix,
                  double[][]rotateMatrix,double[][]position){
        this.intrinsicParameterMatrix=intrinsicParameterMatrix;
        this.extrinsicParameterMatrix=extrinsicParameterMatrix;
        this.rotateMatrix=rotateMatrix;
        this.position=position;
    }

    public Camera(){}

    public int[][]convertWorldToPixel(double[][]worldCoordinate){
        worldCoordinate=MatrixOperator.jointByRow(worldCoordinate,new double[][]{{1}});
        double[][]temp= MatrixOperator.multiply(
                MatrixOperator.multiply(intrinsicParameterMatrix,extrinsicParameterMatrix),worldCoordinate);
        return new int[][]{{new Double(temp[0][0]/temp[2][0]).intValue()},{new Double(temp[1][0]/temp[2][0]).intValue()}};
    }

    public double[][] getDirectionInWorld(int[][]pixelCoordinate){
        double[][]direction={{(double)pixelCoordinate[0][0]},{(double)pixelCoordinate[1][0]},{1}};
        direction=MatrixOperator.multiply(
                MatrixOperator.reverse(intrinsicParameterMatrix),direction);
        direction=MatrixOperator.toUnitVector(direction);
        direction=MatrixOperator.multiply(
                MatrixOperator.transposition(rotateMatrix),direction);
        return new double[][]{{position[0][0]},{position[1][0]},{position[2][0]},
                {direction[0][0]},{direction[1][0]},{direction[2][0]}};
    }
}
