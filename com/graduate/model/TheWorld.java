package com.graduate.model;

import com.graduate.util.FileOperator;
import com.graduate.util.MatrixOperator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class TheWorld {
    Camera[]cameras;
    BallTrack ballTrack;

    //do the first method(linear regression way) and polynomial calibration equation
    public double[][] calibration1P(double[][]distCoordinate, double[][]imageCoordinate){
        int n=distCoordinate[0].length;
        double[][]X=new double[2*n][2];
        double[][]Y=new double[2*n][1];
        int m2Tomm2=1000000;
        for(int i=0;i<n;i++){
            double r2=distCoordinate[0][i]*distCoordinate[0][i]*m2Tomm2+distCoordinate[1][i]*distCoordinate[1][i]*m2Tomm2;
            X[2*i][0]=r2*distCoordinate[0][i];
            X[2*i+1][0]=r2*distCoordinate[1][i];
            X[2*i][1]=r2*r2*distCoordinate[0][i];
            X[2*i+1][1]=r2*r2*distCoordinate[1][i];
            Y[2*i][0]=imageCoordinate[0][i]-distCoordinate[0][i];
            Y[2*i+1][0]=imageCoordinate[1][i]-distCoordinate[1][i];
        }
        return MatrixOperator.linearRegression(X,Y);
    }

    public double[][] getTrueK(double[][]distCoordinate, double[][]imageCoordinate, int order){
        int n=distCoordinate[0].length;
        double[][]X=new double[2*n][order];
        double[][]Y=new double[2*n][1];
        int m2Tomm2=1000000;
        for(int i=0;i<n;i++){
            double r2=imageCoordinate[0][i]*imageCoordinate[0][i]*m2Tomm2+imageCoordinate[1][i]*imageCoordinate[1][i]*m2Tomm2;
            double r2n=r2;
            for(int j=0;j<order;j++){
                X[2*i][j]=imageCoordinate[0][i]*r2n;
                X[2*i+1][j]=imageCoordinate[1][i]*r2n;
                r2n*=r2;
            }
            Y[2*i][0]=distCoordinate[0][i]-imageCoordinate[0][i];
            Y[2*i+1][0]=distCoordinate[1][i]-imageCoordinate[1][i];
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

    public void trackError(double timeMax, String fileName, boolean dist) throws IOException{
        File file=new File(fileName);
        FileWriter fileWriter=new FileWriter(file.getName());
        BufferedWriter bufferedWriter=new BufferedWriter(fileWriter);
        double frameRateHz=340.0;
        double h=1/frameRateHz;
        double time=0.0;

        double maxXError=0.0, maxYError=0.0, maxZError=0.0, maxTotalError=0.0;
        double minXError=10000.0, minYError=10000.0, minZError=10000.0, minTotalError=10000.0;
        double meanXError=0.0, meanYError=0.0, meanZError=0.0, meanTotalError=0.0;
        double devXError=0.0, devYError=0.0, devZError=0.0, devTotalError=0.0;
        double[][]errors=null;

        while (time<=timeMax){
            double[][]worldCoordinate=this.ballTrack.getPosition(time);
            double[][]rayFunctions=null;
            for(Camera camera : this.cameras){
                int[][]pixelCoordinate=null;
                if(dist) {
                    pixelCoordinate=camera.convertWorldToPixelWithDistortion(worldCoordinate);
                }else pixelCoordinate=camera.convertWorldToPixel(worldCoordinate);
                double[][]rayFunction=camera.getDirectionInWorld(pixelCoordinate);
                rayFunctions=MatrixOperator.jointByCol(rayFunctions,rayFunction);
            }
            double[][]WCS=this.rayIntersection(rayFunctions);
            double r2=Math.pow(worldCoordinate[0][0]-WCS[0][0],2)+
                    Math.pow(worldCoordinate[1][0]-WCS[1][0],2)+
                    Math.pow(worldCoordinate[2][0]-WCS[2][0],2);
            double r=Math.sqrt(r2);
            double[][]error={{WCS[0][0]-worldCoordinate[0][0]},
                    {WCS[1][0]-worldCoordinate[1][0]},
                    {WCS[2][0]-worldCoordinate[2][0]},{r}};
            maxXError=Math.max(Math.abs(error[0][0]),maxXError);
            maxYError=Math.max(Math.abs(error[1][0]),maxYError);
            maxZError=Math.max(Math.abs(error[2][0]),maxZError);
            maxTotalError=Math.max(maxTotalError,r);
            minXError=Math.min(Math.abs(error[0][0]),minXError);
            minYError=Math.min(Math.abs(error[1][0]),minYError);
            minZError=Math.min(Math.abs(error[2][0]),minZError);
            minTotalError=Math.min(minTotalError,r);
            errors=MatrixOperator.jointByCol(errors,error);
            String data=""+time+" "+worldCoordinate[0][0]+" "+worldCoordinate[1][0]+" "+worldCoordinate[2][0]+" "+
                    WCS[0][0]+" "+WCS[1][0]+" "+WCS[2][0]+" "+
                    1000*error[0][0]+" "+1000*error[1][0]+" "+1000*error[2][0]+" "+1000*r+"\n";
            bufferedWriter.write(data);
            time+=h;
        }

        int length=errors[0].length;
        for(int i=0;i<length;i++){
            meanXError+=errors[0][i];
            meanYError+=errors[1][i];
            meanZError+=errors[2][i];
            meanTotalError+=errors[3][i];
            devXError+=errors[0][i]*errors[0][i];
            devYError+=errors[1][i]*errors[1][i];
            devZError+=errors[2][i]*errors[2][i];
            devTotalError+=errors[3][i]*errors[3][i];
        }

        meanXError/=length;
        meanYError/=length;
        meanZError/=length;
        meanTotalError/=length;
        devXError=Math.sqrt(devXError/length-meanXError*meanXError);
        devYError=Math.sqrt(devYError/length-meanYError*meanYError);
        devZError=Math.sqrt(devZError/length-meanZError*meanZError);
        devTotalError=Math.sqrt(devTotalError/length-meanTotalError*meanTotalError);
        System.out.println("");
        System.out.println("maxXError: "+1000*maxXError+" maxYError: "+1000*maxYError+" maxZError: "+1000*maxZError+" maxTotalError :"+1000*maxTotalError);
        System.out.println("minXError: "+1000*minXError+" minYError: "+1000*minYError+" minZError: "+1000*minZError+" minTotalError :"+1000*minTotalError);
        System.out.println("meanXError: "+1000*meanXError+" meanYError: "+1000*meanYError+" meanZError: "+1000*meanZError+" meanTotalError :"+1000*meanTotalError);
        System.out.println("devXError: "+1000*devXError+" devYError: "+1000*devYError+" devZError: "+1000*devZError+" devTotalError :"+1000*devTotalError);
        System.out.println("");
        bufferedWriter.close();
    }

    public static void main(String[] args) throws IOException {
        String[]cameraFile  = FileOperator.readFromDatFile("cameraData.dat");
        double[][]cameraData= FileOperator.dealTheData(cameraFile);
        TheWorld theWorld=new TheWorld();
        theWorld.cameras=CameraMaker.makeCameras(cameraData);
        String[] ballFile = FileOperator.readFromDatFile("ball.dat");
        double[][]ballData=FileOperator.dealTheData(ballFile);
        theWorld.ballTrack=new BallTrack(ballData);
        String[] courtFile = FileOperator.readFromDatFile("court.dat");
        double[][]court=FileOperator.dealTheData(courtFile);

        double frameRateHz=340.0;
        double h=1/frameRateHz;
        double time=0.0;
        double timeMax=ballData[ballData.length-1][0];

        System.out.println("=====Camera without lens distortion=========");
        theWorld.trackError(timeMax,"track.dat", false);
        System.out.println("=====Camera with lens distortion=========");
        theWorld.trackError(timeMax,"trackD.dat", true);

        double[][]Ks=null;
        double[][]Kinverses=null;
        for(int i=0;i<theWorld.cameras.length;i++){
            Camera camera=theWorld.cameras[i];
            double[][]distCoordinates=camera.pixelToDist(camera.convertWorldToPixelWithDistortion(court));
            double[][]imageCoordinates=camera.worldToImage(court);
            double[][]K=theWorld.getTrueK(distCoordinates,imageCoordinates,1);
            double K1=K[0][0];
            double[][]Kinverse=new double[1][1];
            Kinverse[0][0]=-K1;
            camera.setKinverse(Kinverse);
            Ks=MatrixOperator.jointByCol(Ks,K);
            Kinverses=MatrixOperator.jointByCol(Kinverses,Kinverse);
        }

        System.out.println("=====Camera calibrated by fourth order inverse lens distortion=========");
        System.out.println("estimated K:");
        MatrixOperator.nicePrint(Ks);
        System.out.println("inverse series K:");
        MatrixOperator.nicePrint(Kinverses);
        theWorld.trackError(timeMax, "trackC1.dat", true);

        Ks=null;
        Kinverses=null;
        for(int i=0;i<theWorld.cameras.length;i++){
            Camera camera=theWorld.cameras[i];
            double[][]distCoordinates=camera.pixelToDist(camera.convertWorldToPixelWithDistortion(court));
            double[][]imageCoordinates=camera.worldToImage(court);
            double[][]K=theWorld.getTrueK(distCoordinates,imageCoordinates,2);
            double K1=K[0][0];
            double K2=K[1][0];
            double[][]Kinverse=new double[2][1];
            Kinverse[0][0]=-K1;
            Kinverse[1][0]=3*K1*K1-K2;
            camera.setKinverse(Kinverse);
            Ks=MatrixOperator.jointByCol(Ks,K);
            Kinverses=MatrixOperator.jointByCol(Kinverses,Kinverse);
        }
        System.out.println("=====Camera calibrated by second order inverse lens distortion=========");
        System.out.println("estimated K:");
        MatrixOperator.nicePrint(Ks);
        System.out.println("inverse series K:");
        MatrixOperator.nicePrint(Kinverses);
        theWorld.trackError(timeMax, "trackC2.dat", true);

        Ks=null;
        Kinverses=null;
        for(int i=0;i<theWorld.cameras.length;i++){
            Camera camera=theWorld.cameras[i];
            double[][]distCoordinates=camera.pixelToDist(camera.convertWorldToPixelWithDistortion(court));
            double[][]imageCoordinates=camera.worldToImage(court);
            double[][]K=theWorld.getTrueK(distCoordinates,imageCoordinates,3);
            double K1=K[0][0];
            double K2=K[1][0];
            double K3=K[2][0];
            double[][]Kinverse=new double[3][1];
            Kinverse[0][0]=-K1;
            Kinverse[1][0]=3*K1*K1-K2;
            Kinverse[2][0]=-12*K1*K1*K1+8*K1*K2-K3;
            camera.setKinverse(Kinverse);
            Ks=MatrixOperator.jointByCol(Ks,K);
            Kinverses=MatrixOperator.jointByCol(Kinverses,Kinverse);
        }
        System.out.println("=====Camera calibrated by third order inverse lens distortion=========");
        System.out.println("estimated K:");
        MatrixOperator.nicePrint(Ks);
        System.out.println("inverse series K:");
        MatrixOperator.nicePrint(Kinverses);
        theWorld.trackError(timeMax, "trackC3.dat", true);
    }
}
