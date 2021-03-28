package com.graduate.model;

import com.graduate.util.FileOperator;
import com.graduate.util.MatrixOperator;

public class TheWorld {
    Camera[]cameras;
    BallTrack ballTrack;

    //do the first method(linear regression way) and polynomial calibration equation
    public double[][] calibration1P(double[][]distCoordinate, double[][]imageCoordinate){
        int n=distCoordinate[0].length;
        double[][]X=new double[2*n][2];
        double[][]Y=new double[2*n][1];
        for(int i=0;i<n;i++){
            X[2*i][0]=1;
            X[2*i+1][0]=1;
            double r2=distCoordinate[0][i]*distCoordinate[0][i]+distCoordinate[1][i]*distCoordinate[1][i];
            X[2*i][1]=r2;
            X[2*i+1][1]=r2;
            Y[2*i][0]=(imageCoordinate[0][i]-distCoordinate[0][i])/(distCoordinate[0][i]*r2);
            Y[2*i+1][0]=(imageCoordinate[1][i]-distCoordinate[1][i])/(distCoordinate[1][i]*r2);
        }
        return MatrixOperator.linearRegression(X,Y);
    }

    public double[][] rayIntersection(double[][]rayFunctions){
        int length=rayFunctions.length;
        double[][]augmentationMatrix=new double[6][length];
        for(int i=0;i<length;i++){
            augmentationMatrix[0][i]=1-rayFunctions[3][i]*rayFunctions[3][i];   //1-ai^2
            augmentationMatrix[1][i]=1-rayFunctions[4][i]*rayFunctions[4][i];   //1-bi^2
            augmentationMatrix[2][i]=1-rayFunctions[5][i]*rayFunctions[5][i];   //1-ci^2
            augmentationMatrix[3][i]=rayFunctions[3][i]*rayFunctions[4][i];     //ai*bi
            augmentationMatrix[4][i]=rayFunctions[3][i]*rayFunctions[5][i];     //ai*ci
            augmentationMatrix[5][i]=rayFunctions[4][i]*rayFunctions[5][i];     //bi*ci
        }
        double[][]A=new double[3][3];
        double[][]b=new double[3][1];
        A[0][0]=MatrixOperator.innerProduct(augmentationMatrix[0]);
        A[0][1]=-MatrixOperator.innerProduct(augmentationMatrix[3]);
        A[0][2]=-MatrixOperator.innerProduct(augmentationMatrix[4]);
        A[1][0]=-MatrixOperator.innerProduct(augmentationMatrix[3]);
        A[1][1]=MatrixOperator.innerProduct(augmentationMatrix[1]);
        A[1][2]=-MatrixOperator.innerProduct(augmentationMatrix[5]);
        A[2][0]=-MatrixOperator.innerProduct(augmentationMatrix[4]);
        A[2][1]=-MatrixOperator.innerProduct(augmentationMatrix[5]);
        A[2][2]=MatrixOperator.innerProduct(augmentationMatrix[2]);
        b[0][0]=MatrixOperator.innerProduct(augmentationMatrix[0],rayFunctions[0])
                -MatrixOperator.innerProduct(augmentationMatrix[3],rayFunctions[1])
                -MatrixOperator.innerProduct(augmentationMatrix[4],rayFunctions[2]);
        b[1][0]=-MatrixOperator.innerProduct(augmentationMatrix[3],rayFunctions[0])
                +MatrixOperator.innerProduct(augmentationMatrix[1],rayFunctions[1])
                -MatrixOperator.innerProduct(augmentationMatrix[5],rayFunctions[2]);
        b[2][0]=-MatrixOperator.innerProduct(augmentationMatrix[4],rayFunctions[0])
                -MatrixOperator.innerProduct(augmentationMatrix[5],rayFunctions[1])
                +MatrixOperator.innerProduct(augmentationMatrix[2],rayFunctions[2]);
        return MatrixOperator.solveLinearEquation(A,b);
    }

    public static void main(String[] args) {
        String[]cameraFile  = FileOperator.readFromDatFile("cameraData.dat");
        double[][]cameraData= FileOperator.dealTheData(cameraFile);
        TheWorld theWorld=new TheWorld();
        theWorld.cameras=CameraMaker.makeCameras(cameraData);
        String[] ballFile = FileOperator.readFromDatFile("ball.dat");
        double[][]ballData=FileOperator.dealTheData(ballFile);
        theWorld.ballTrack=new BallTrack(ballData);
        double[][]court={
                {0.0,0.0,-4.115,4.115,-4.115,4.115,-4.115,-4.115,4.115,4.115},
                {6.4,-6.4,6.4,6.4,-6.4,-6.4,11.885,-11.885,11.885,-11.885},
                {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01}};
        double K1P=0;
        double K2P=0;
        for(int i=0;i<theWorld.cameras.length;i++){
            Camera camerai=theWorld.cameras[i];
            double[][]distCoordinates=camerai.pixelToDist(camerai.convertWorldToPixelWithDistortion(court));
            double[][]imageCoordinates=camerai.worldToImage(court);
            double[][]KP=theWorld.calibration1P(distCoordinates,imageCoordinates);
            MatrixOperator.nicePrint(KP);
            MatrixOperator.nicePrint(imageCoordinates);
            MatrixOperator.nicePrint(distCoordinates);
            System.out.println("==============================================================");
            K1P+=KP[0][0];
            K2P+=KP[1][0];
        }
        K1P/=theWorld.cameras.length;
        K2P/=theWorld.cameras.length;
        System.out.println(K1P+" "+K2P);
        for(Camera camera:theWorld.cameras)camera.setK1PandK2P(K1P,K2P);
//        double[][]position=theWorld.ballTrack.getPosition(0.677);
//        double[][]rayFunctions=null;
//        for(Camera camera:theWorld.cameras){
//            int[][]pixelCoordinate=camera.convertWorldToPixel(position);
//            double[][]rayFunction=camera.getDirectionInWorldP(pixelCoordinate);
//            rayFunctions= MatrixOperator.jointByCol(rayFunctions,rayFunction);
//        }
//        MatrixOperator.nicePrint(position);
//        MatrixOperator.nicePrint(theWorld.rayIntersection(rayFunctions));
//        MatrixOperator.nicePrint(rayFunctions);
    }



}
