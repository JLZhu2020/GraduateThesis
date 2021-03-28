package com.graduate.model;

import com.graduate.util.MatrixOperator;

public class Camera {

    private double[][] intrinsicParameterMatrix;
    private double[][] extrinsicParameterMatrix;
    private double[][] rotateMatrix;
    private double[][] position;
    private double K1;      //radial distortion coefficient
    private double K2;      //radial distortion coefficient
    private double K1P;
    private double K2P;
    private double K1D;
    private double K2D;
    private double focalLength;
    private double[][] imageToPixel;
    private double[][] cameraToImage;

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
        this.K1P=0.0;
        this.K2P=0.0;
        this.K1D=0.0;
        this.K2D=0.0;
        this.cameraToImage=new double[][]{{focalLength,0,0},{0,focalLength,0},{0,0,1}};
        double[][]imageToPixel={intrinsicParameterMatrix[0].clone(),
                        intrinsicParameterMatrix[1].clone(),intrinsicParameterMatrix[2].clone()};
        imageToPixel[0][0]/=focalLength;
        imageToPixel[1][1]/=focalLength;
        this.imageToPixel=imageToPixel;
    }

    public void setK1PandK2P(double K1P, double K2P){
        this.K1P=K1P;
        this.K2P=K2P;
    }

    public void setK1DandK2D(double K1D, double K2D){
        this.K1D=K1D;
        this.K2D=K2D;
    }

    public Camera(){}

    //only pinhole camera model without lens distortion could do this process
    public int[][]convertWorldToPixel(double[][]worldCoordinate){
        int n=worldCoordinate[0].length;
        double[][]imageCoordinate= worldToImage(worldCoordinate);
        double[][]pixelCoordinate=MatrixOperator.multiply(imageToPixel,imageCoordinate);
        int[][]res=new int[2][n];
        for(int i=0;i<n;i++){
            res[0][i]=new Double(pixelCoordinate[0][i]).intValue();
            res[1][i]=new Double(pixelCoordinate[1][i]).intValue();
        }
        return res;
    }

    //this function will use given K1 and K2 to imitate lens distortion
    public int[][]convertWorldToPixelWithDistortion(double[][]worldCoordinate){
        int n=worldCoordinate[0].length;
        double[][]imageCoordinate= worldToImage(worldCoordinate);
        for(int i=0;i<n;i++){
            double r2=imageCoordinate[0][i]*imageCoordinate[0][i]+imageCoordinate[1][i]*imageCoordinate[1][i];
            double rate=1+K1*r2+K2*r2*r2;
            imageCoordinate[0][i]*=rate;
            imageCoordinate[1][i]*=rate;
        }
        double[][]pixelCoordinate=MatrixOperator.multiply(imageToPixel,imageCoordinate);
        int[][]res=new int[2][n];
        for(int i=0;i<n;i++){
            res[0][i]=(int)Math.round(pixelCoordinate[0][i]);//new Double(pixelCoordinate[0][i]).intValue();
            res[1][i]=(int)Math.round(pixelCoordinate[1][i]);//new Double(pixelCoordinate[1][i]).intValue();
        }
        return res;
    }

    public double[][] worldToImage(double[][]worldCoordinate){
        int n=worldCoordinate[0].length;
        double[][]lastRow=new double[1][n];
        for(int i=0;i<n;i++)lastRow[0][i]=1;
        worldCoordinate=MatrixOperator.jointByRow(worldCoordinate,lastRow);
        double[][]cameraCoordinate= MatrixOperator.multiply(extrinsicParameterMatrix,worldCoordinate);
        double[][]imageCoordinate= MatrixOperator.multiply(cameraToImage,cameraCoordinate);
        for(int i=0;i<n;i++){
            imageCoordinate[0][i]/=imageCoordinate[2][i];
            imageCoordinate[1][i]/=imageCoordinate[2][i];
            imageCoordinate[2][i]=1.0;
        }
        return imageCoordinate;
    }

    public double[][] pixelToDist(int[][]pixelCoordinate){
        int n=pixelCoordinate[0].length;
        double[][]augPixel=new double[3][n];
        for(int i=0;i<n;i++){
            augPixel[0][i]=(double)pixelCoordinate[0][i];
            augPixel[1][i]=(double)pixelCoordinate[1][i];
            augPixel[2][i]=1.0;
        }
        double[][] pixelToDist=MatrixOperator.reverse(imageToPixel);
        return MatrixOperator.multiply(pixelToDist,augPixel);
    }

    public double[][] getDirectionInWorldP(int[][]pixelCoordinate){
        double[][]augPixel={{(double)pixelCoordinate[0][0]},{(double)pixelCoordinate[1][0]},{1}};
        double[][]pixelToImage=MatrixOperator.reverse(imageToPixel);
        double[][]imageCoordinate=MatrixOperator.multiply(pixelToImage,augPixel);
        double r2=imageCoordinate[0][0]*imageCoordinate[0][0]+imageCoordinate[1][0]*imageCoordinate[1][0];
        double rate=1+K1P*r2+K2P*r2*r2;
        imageCoordinate[0][0]*=rate;
        imageCoordinate[1][0]*=rate;
        double[][]imageToCamera=MatrixOperator.reverse(cameraToImage);
        double[][]cameraCoordinate=MatrixOperator.multiply(imageToCamera,imageCoordinate);
        double[][]rotateT=MatrixOperator.transposition(rotateMatrix);
        double[][]directInWorld=MatrixOperator.multiply(rotateT,cameraCoordinate);
        directInWorld=MatrixOperator.toUnitVector(directInWorld);
        return new double[][]{{position[0][0]},{position[1][0]},{position[2][0]},
                {directInWorld[0][0]},{directInWorld[1][0]},{directInWorld[2][0]}};
    }

    public double[][] getDirectionInWorldD(int[][]pixelCoordinate){
        double[][]augPixel={{(double)pixelCoordinate[0][0]},{(double)pixelCoordinate[1][0]},{1}};
        double[][]pixelToImage=MatrixOperator.reverse(imageToPixel);
        double[][]imageCoordinate=MatrixOperator.multiply(pixelToImage,augPixel);
        double r2=imageCoordinate[0][0]*imageCoordinate[0][0]+imageCoordinate[1][0]*imageCoordinate[1][0];
        double rate=1+K1D*r2+K2D*r2*r2;
        imageCoordinate[0][0]/=rate;
        imageCoordinate[1][0]/=rate;
        double[][]imageToCamera=MatrixOperator.reverse(cameraToImage);
        double[][]cameraCoordinate=MatrixOperator.multiply(imageToCamera,imageCoordinate);
        double[][]rotateT=MatrixOperator.transposition(rotateMatrix);
        double[][]directInWorld=MatrixOperator.multiply(rotateT,cameraCoordinate);
        directInWorld=MatrixOperator.toUnitVector(directInWorld);
        return new double[][]{{position[0][0]},{position[1][0]},{position[2][0]},
                {directInWorld[0][0]},{directInWorld[1][0]},{directInWorld[2][0]}};
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
