package com.graduate.model;

import com.graduate.util.MatrixOperator;

public class Camera {

    private double[][] intrinsicParameterMatrix;
    private double[][] extrinsicParameterMatrix;
    private double[][] rotateMatrix;
    private double[][] position;
    private double K1;      //radial distortion coefficient
    private double K2;      //radial distortion coefficient
    private double[][] Kinverse;
    private double focalLength;
    private double[][] imageToPixel;
    private double[][] cameraToImage;

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
        this.cameraToImage=new double[][]{{focalLength,0,0},{0,focalLength,0},{0,0,1}};
        double[][]imageToPixel={intrinsicParameterMatrix[0].clone(),
                        intrinsicParameterMatrix[1].clone(),intrinsicParameterMatrix[2].clone()};
        imageToPixel[0][0]/=focalLength;
        imageToPixel[1][1]/=focalLength;
        this.imageToPixel=imageToPixel;
    }

    public void setKinverse(double[][]Kinverse){
        this.Kinverse=Kinverse;
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
        int m2Tomm2=1000000;
        for(int i=0;i<n;i++){
            double r2=imageCoordinate[0][i]*imageCoordinate[0][i]*m2Tomm2+imageCoordinate[1][i]*imageCoordinate[1][i]*m2Tomm2;
            double rate=1+K1*r2+K2*r2*r2;
            imageCoordinate[0][i]*=rate;
            imageCoordinate[1][i]*=rate;
        }
        double[][]pixelCoordinate=MatrixOperator.multiply(imageToPixel,imageCoordinate);
        int[][]res=new int[2][n];
        for(int i=0;i<n;i++){
            res[0][i]=new Double(Math.floor(pixelCoordinate[0][i])).intValue();
            res[1][i]=new Double(Math.floor(pixelCoordinate[1][i])).intValue();
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
            augPixel[0][i]=(double)pixelCoordinate[0][i]+0.5;
            augPixel[1][i]=(double)pixelCoordinate[1][i]+0.5;
            augPixel[2][i]=1.0;
        }
        double[][] pixelToDist=MatrixOperator.reverse(imageToPixel);
        return MatrixOperator.multiply(pixelToDist,augPixel);
    }

    public double[][] getDirectionInWorld(int[][]pixelCoordinate){
        double[][]augPixel={{(double)pixelCoordinate[0][0]+0.5},{(double)pixelCoordinate[1][0]+0.5},{1}};
        double[][]pixelToImage=MatrixOperator.reverse(imageToPixel);
        double[][]imageCoordinate=MatrixOperator.multiply(pixelToImage,augPixel);
        if(outOfScreen(imageCoordinate))return null;
        double rate=1;
        if(this.Kinverse!=null){
            int m2Tomm2=1000000;
            double r2=imageCoordinate[0][0]*imageCoordinate[0][0]*m2Tomm2+imageCoordinate[1][0]*imageCoordinate[1][0]*m2Tomm2;
            double r2n=r2;
            for(int i=0;i<Kinverse.length;i++){
                rate+=Kinverse[i][0]*r2n;
                r2n*=r2;
            }
        }
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

    public boolean outOfScreen(double[][]imageCoordinate){
        if(Math.abs(imageCoordinate[0][0])<=0.00768&&Math.abs(imageCoordinate[1][0])<=0.006144){
            return false;
        }else return true;
    }
}
