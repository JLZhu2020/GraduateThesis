package com.graduate.model;

import com.graduate.util.MatrixOperator;

public class Camera {

    private double[][] intrinsicParameterMatrix;
    private double[][] extrinsicParameterMatrix;
    private double[][] rotateMatrix;
    private double[][] position;
    private double K1;      //radial distortion coefficient
    private double K2;      //radial distortion coefficient
    private double focalLength;

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

    //By default, the constructor will return a pinhole camera without lens distortion
    public Camera(double[][]intrinsicParameterMatrix, double[][]extrinsicParameterMatrix,
                  double[][]rotateMatrix, double[][]position, double focalLength, double K1, double K2){
        this.intrinsicParameterMatrix=intrinsicParameterMatrix;
        this.extrinsicParameterMatrix=extrinsicParameterMatrix;
        this.rotateMatrix=rotateMatrix;
        this.position=position;
        this.focalLength=focalLength;
        this.K1=K1;
        this.K2=K2;
    }

    public Camera(){}

    //only pinhole camera model without lens distortion could do this process
    public int[][]convertWorldToPixel(double[][]worldCoordinate){
        worldCoordinate=MatrixOperator.jointByRow(worldCoordinate,new double[][]{{1}});
        double[][]temp= MatrixOperator.multiply(
                MatrixOperator.multiply(intrinsicParameterMatrix,extrinsicParameterMatrix),worldCoordinate);
        return new int[][]{{new Double(temp[0][0]/temp[2][0]).intValue()},{new Double(temp[1][0]/temp[2][0]).intValue()}};
    }

    //this function will use given K1 and K2 to imitate lens distortion
    public int[][]convertWorldToPixelWithDistortion(double[][]worldCoordinate){
        worldCoordinate=MatrixOperator.jointByRow(worldCoordinate,new double[][]{{1}});
        double[][]temp= MatrixOperator.multiply(
                MatrixOperator.multiply(intrinsicParameterMatrix,extrinsicParameterMatrix),worldCoordinate);
        return new int[][]{{new Double(temp[0][0]/temp[2][0]).intValue()},{new Double(temp[1][0]/temp[2][0]).intValue()}};
    }

    public double[][] getDirectionInWorld(int[][]pixelCoordinate){
        double[][]direction={{(double)pixelCoordinate[0][0]},{(double)pixelCoordinate[1][0]},{1}};
        double[][]ImageToCCS={{1/focalLength,0,0},{0,1/focalLength,0},{0,0,1}};
        double[][]DistToPixel={
                intrinsicParameterMatrix[0].clone(),
                intrinsicParameterMatrix[1].clone(),
                intrinsicParameterMatrix[2].clone()};
        DistToPixel[0][0]/=focalLength;
        DistToPixel[1][1]/=focalLength;
        direction=MatrixOperator.multiply(
                MatrixOperator.reverse(DistToPixel),direction);
        direction=eliminateLensDistortion(direction);
        direction=MatrixOperator.multiply(ImageToCCS,direction);
        direction=MatrixOperator.toUnitVector(direction);
        direction=MatrixOperator.multiply(
                MatrixOperator.transposition(rotateMatrix),direction);
        return new double[][]{{position[0][0]},{position[1][0]},{position[2][0]},
                {direction[0][0]},{direction[1][0]},{direction[2][0]}};
    }

    //if camera's parameter is all known, could apply linear regression to determine K1 and K2
    //input WCS coordinate and pixel coordinate are augment by last line all 1
    public void calibrateK1K2(double[][]WCS, double[][]pixel){
        double[][]CCSToImage={{focalLength,0,0},{0,focalLength,0},{0,0,1}};
        double[][]DistToPixel={
                intrinsicParameterMatrix[0].clone(),
                intrinsicParameterMatrix[1].clone(),
                intrinsicParameterMatrix[2].clone()};
        DistToPixel[0][0]/=focalLength;
        DistToPixel[1][1]/=focalLength;
        double[][]PixelToDist=MatrixOperator.reverse(DistToPixel);
        double[][]Dist=MatrixOperator.multiply(PixelToDist,pixel);
        double[][]Image=MatrixOperator.multiply(CCSToImage,
                MatrixOperator.multiply(extrinsicParameterMatrix,WCS));
        int length=Dist[0].length;
        double[][]X=new double[2*length][2];
        double[][]Y=new double[2*length][1];
        for(int i=0;i<length;i++){
            double rSquare=Dist[0][i]*Dist[0][i]+Dist[1][i]*Dist[1][i];
            X[2*i][0]=1;
            X[2*i+1][0]=1;
            X[2*i][1]=rSquare;
            X[2*i+1][1]=rSquare;
            Y[2*i][0]=(Image[0][i]-Dist[0][i])/rSquare;
            Y[2*i+1][0]=(Image[1][i]-Dist[1][i])/rSquare;
        }
        double[][]K=MatrixOperator.linearRegression(X,Y);
        this.K1=K[0][0];
        this.K2=K[1][0];
    }

    /*The radial lens distortion equation is: Ximage=Xdist(1+K1*r^2+K2*r^4)
    *                                         Yimage=Ydist(1+K1*r^2+K2*r^4)
    * */
    public double[][] eliminateLensDistortion(double[][]distortedCoordinate){
        double xdist=distortedCoordinate[0][0];
        double ydist=distortedCoordinate[1][0];
        double rSquare=xdist*xdist+ydist*ydist;
        return new double[][]{
                {xdist*(1+K1*rSquare+K2*rSquare*rSquare)},
                {ydist*(1+K1*rSquare+K2*rSquare*rSquare)},
                {1}};
    }
}
