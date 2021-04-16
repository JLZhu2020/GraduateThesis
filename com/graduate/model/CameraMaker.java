package com.graduate.model;

import com.graduate.util.MatrixOperator;

import java.util.Arrays;

public class CameraMaker {

    public static Camera makeOneCamera(double[]data){
        if(data==null||data.length!=10){
            return null;
        }
        double[][]position=MatrixOperator.makeVector(new double[]{data[0],data[1],data[2]});
        double[][]rotateMatrix=calculateRotateMatrix(data[3],data[4]);
        double[][]transferVector=MatrixOperator.multiply(rotateMatrix,position);
        MatrixOperator.negative(transferVector);
        double[][]extrinsicParameterMatrix=MatrixOperator.jointByCol(rotateMatrix,transferVector);
//        double[][]lastLine={{0},{0},{0},{1}};
//        extrinsicParameterMatrix=MatrixOperator.jointByRow(extrinsicParameterMatrix,lastLine);
        double fieldOfView=data[5];
        int resX= (int) data[6];
        int resY= (int) data[7];
        double dX=data[8];
        double dY=data[9];
        double iMin=dX*(double)Math.min(resX,resY);
        double focalLength = iMin / (2.0 * Math.tan(0.5 * fieldOfView * Math.PI/180.0));
        double[][]intrinsicParameterMatrix={
                {focalLength/dX,0,resX/2},
                {0,focalLength/dY,resY/2},
                {0,0,1}};
        return new Camera(intrinsicParameterMatrix,extrinsicParameterMatrix,rotateMatrix,
                position,focalLength,-1.28E-4, 1.61E-6);
    }

    public static Camera[] makeCameras(double[][]data){
        if(data==null||data.length==0||data[0].length!=10){
            return null;
        }
        Camera[] cameras=new Camera[data.length];
        for(int i=0;i<data.length;i++){
            cameras[i]=makeOneCamera(data[i]);
        }
        return cameras;
    }

    public static double[][] calculateRotateMatrix(double heading, double elevation, double...another){
        if(another.length>1){
            System.out.println("to much input parameter");
            return null;
        }
        heading=heading*Math.PI/180;
        elevation=elevation*Math.PI/180;
        double cosHeading   = Math.cos(heading);
        double sinHeading   = Math.sin(heading);
        double cosElevation = Math.cos(elevation);
        double sinElevation = Math.sin(elevation);
        double[][]rotateMatrix={
                {-cosHeading,sinHeading,0.0},
                {-sinHeading * sinElevation, -cosHeading * sinElevation, cosElevation},
                {sinHeading * cosElevation, cosHeading * cosElevation, sinElevation}};
        return rotateMatrix;
    }


    public static void main(String[] args) {
        calculateRotateMatrix(1.0,1.0);

    }
}
