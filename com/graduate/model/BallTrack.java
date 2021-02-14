package com.graduate.model;

public class BallTrack {

    private double[][]sourceData;

    private double[]tBall;
    private double[]xBall;
    private double[]yBall;
    private double[]zBall;
    private double[]vxBall;
    private double[]vyBall;
    private double[]vzBall;

    public BallTrack(double[][]sourceData){
        this.setSourceData(sourceData);
    }

    public BallTrack(){}

    public double[][] getSourceData() {
        return sourceData;
    }

    public void setSourceData(double[][] sourceData) {
        this.sourceData = sourceData;
        if(sourceData==null||sourceData.length==0||sourceData[0].length!=7)return;
        int len=sourceData.length;
        tBall=new double[len];
        xBall=new double[len];
        yBall=new double[len];
        zBall=new double[len];
        vxBall=new double[len];
        vyBall=new double[len];
        vzBall=new double[len];
        for(int i=0;i<len;i++){
           tBall[i]=sourceData[i][0];
           xBall[i]=sourceData[i][1];
           yBall[i]=sourceData[i][2];
           zBall[i]=sourceData[i][3];
           vxBall[i]=sourceData[i][4];
           vyBall[i]=sourceData[i][5];
           vzBall[i]=sourceData[i][6];
        }
    }

    //This function uses cubic Hermite interpolation to generate in between positions.
    public double[][] getPosition(double t){
        int len=tBall.length;
        if(t<=tBall[0]) {
            return new double[][]{{xBall[0]}, {yBall[0]}, {zBall[0]}};
        }else if(t>=tBall[len-1]){
            return new double[][]{{xBall[len-1]}, {yBall[len-1]}, {zBall[len-1]}};
        }else{
            int i=findT(t);
            if(tBall[i]==t){
                return new double[][]{{xBall[len-1]}, {yBall[len-1]}, {zBall[len-1]}};
            }
            double[][]res=new double[3][1];
            double tA = tBall[i];
            double tB = tBall[i+1];

            double deltaT        = tB - tA;
            double deltaTSquared = deltaT * deltaT;

            double fA    = (1.0 + 2.0 * (t - tA) / deltaT) * (tB - t) * (tB - t) / deltaTSquared;
            double fB    = (1.0 + 2.0 * (tB - t) / deltaT) * (t - tA) * (t - tA) / deltaTSquared;
            double fADot =  (t - tA) * (tB - t) * (tB - t) / deltaTSquared;
            double fBDot = -(t - tA) * (t - tA) * (tB - t) / deltaTSquared;

            res[0][0] = xBall[i] * fA + xBall[i+1] * fB + vxBall[i] * fADot + vxBall[i+1] * fBDot;
            res[1][0] = yBall[i] * fA + yBall[i+1] * fB + vyBall[i] * fADot + vyBall[i+1] * fBDot;
            res[2][0] = zBall[i] * fA + zBall[i+1] * fB + vzBall[i] * fADot + vzBall[i+1] * fBDot;
            return res;
        }


    }

    public int findT(double t){
        int i=0;
        int j=tBall.length-1;
        while(i<j-1){
            int mid=(i+j)/2;
            if(tBall[mid]>t){
                j=mid;
            }else i=mid;
        }
        return i;
    }

}
