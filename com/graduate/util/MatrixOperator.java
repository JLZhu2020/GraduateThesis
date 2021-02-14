package com.graduate.util;

import java.util.Arrays;

public class MatrixOperator {

    //return matrix is L, input will be the U
    public static double[][] luDecomposition(double[][]matrix){
        if(matrix==null||matrix.length!=matrix[0].length||det(matrix)==0){
            System.out.println("Could not reverse");
            return null;
        }
        int len=matrix.length;
        double[][] matrixI=matrixI(len);
        for(int i=0;i<len;i++){
//            if(matrix[i][i]==0.0){
//                int j=i+1;
//                while(matrix[j][i]==0.0)j++;
//                exchangeRow(matrix,i,j);
//                exchangeRow(matrixI,i,j);
//            }
            for(int j=i+1;j<len;j++){
                rowJminusrowI(matrixI,i,j,matrix[j][i]/matrix[i][i]);
                rowJminusrowI(matrix,i,j,matrix[j][i]/matrix[i][i]);
            }
        }
        return reverse(matrixI);
    }

    public static void choleskyFactoriation(double[][]matrix){
        if(matrix==null||matrix.length!=matrix[0].length||det(matrix)==0||!isSymmetric(matrix)){
            System.out.println("Could not factor");
            return;
        }
        int len=matrix.length;
        for(int i=0;i<len;i++){
            if(matrix[i][i]==0.0){
                int j=i+1;
                while(matrix[j][i]==0.0)j++;
                exchangeRow(matrix,i,j);
            }
            for(int j=i+1;j<len;j++){
                double t=matrix[j][i]/matrix[i][i];
                rowJminusrowI(matrix,i,j,t);
            }
        }
        for(int i=0;i<len;i++){
            double k=Math.sqrt(matrix[i][i]);
            rowMultiply(matrix,i,1/k);
        }
    }

    public static double[] getEigenValue(double[][]matrix){

        return null;
    }

    public static double[][] solveLinearEquation(double[][]matrixA,double[][]vectorb){
        if(matrixA==null||matrixA.length!=matrixA[0].length||vectorb==null||matrixA.length!=vectorb.length){
            System.out.println("Unsolvable function!");
            return null;
        }else if(det(matrixA)==0.0){
            return new double[matrixA.length][1];
        }
        int len=matrixA.length;
        double[][] x=new double[len][1];
        for(int i=0;i<len;i++){
            if(matrixA[i][i]==0.0){
                int j=i+1;
                while(matrixA[j][i]==0.0)j++;
                exchangeRow(matrixA,i,j);
                exchangeRow(vectorb,i,j);
            }
            rowMultiply(vectorb,i,1/matrixA[i][i]);
            rowMultiply(matrixA,i,1/matrixA[i][i]);
            for(int j=i+1;j<matrixA.length;j++){
                rowJminusrowI(vectorb,i,j,matrixA[j][i]);
                rowJminusrowI(matrixA,i,j,matrixA[j][i]);
            }
        }
        for(int i=len-1;i>=0;i--){
            for(int j=len-1;j>i;j--){
                vectorb[i][0]-=matrixA[i][j]*x[j][0];
            }
            x[i][0]=vectorb[i][0]/matrixA[i][i];
        }
        return x;
    }

    //the vector b or x should be a column vector
    public static double[][] makeVector(double[]array){
        if(array==null)return null;
        double[][]vector=new double[array.length][1];
        for(int i=0;i<array.length;i++){
            vector[i][0]=array[i];
        }
        return vector;
    }

    public static double[][] multiply(double[][] matrixA,double [][] matrixB){
        if(matrixA==null||matrixB==null||matrixA[0].length!=matrixB.length){
            System.out.println("Could not multiply");
            return null;
        }
        double[][]matrixC=new double[matrixA.length][matrixB[0].length];
        for(int i=0;i<matrixC.length;i++){
            for(int j=0;j<matrixC[0].length;j++){
                double val=0.0;
                for(int k=0;k<matrixA[0].length;k++){
                    val+=matrixA[i][k]*matrixB[k][j];
                }
                matrixC[i][j]=val;
            }
        }
        return matrixC;
    }

    public static void exchangeRow(double[][]matrix,int i,int j){
        if(i<0||j<0||matrix==null||i>=matrix.length||j>=matrix.length){
            System.out.println("Could not exchangeRow");
            return;
        }
        for(int k=0;k<matrix[0].length;k++){
            double temp=matrix[i][k];
            matrix[i][k]=matrix[j][k];
            matrix[j][k]=temp;
        }
    }

    //its for calculating det; remove ith row and jth col from matrix and return.
    public static double[][] subMatrix(double[][]matrix,int i,int j){
        if(i<0||j<0||matrix==null||i>=matrix.length||j>=matrix[0].length){
            System.out.println("Do not have such subMatrix");
            return matrix;
        }
        double[][]subMatrix=new double[matrix.length-1][matrix[0].length-1];
        int rowPlus=0;
        int colPlus=0;
        for(int l=0;l<subMatrix.length;l++){
            if(l==i){
                rowPlus=1;
            }
            colPlus=0;
            for(int m=0;m<subMatrix[0].length;m++){
                if(m==j){
                    colPlus=1;
                }
                subMatrix[l][m]=matrix[l+rowPlus][m+colPlus];
            }
        }
        return subMatrix;
    }

    public static double det(double[][]matrix){
        if(matrix==null||matrix.length!=matrix[0].length){
            System.out.println("Do not have Det");
            return 0.0;
        }
        if(matrix.length==1)return matrix[0][0];
        double det=0.0;
        for(int i=0;i<matrix.length;i++){
            if(i%2==0){
                det+=matrix[0][i]*det(subMatrix(matrix,0,i));
            }else det-=matrix[0][i]*det(subMatrix(matrix,0,i));
        }
        return det;
    }

    public static double[][] matrixI(int i){
        if(i<=0){
            return null;
        }
        double[][] matrixI=new double[i][i];
        for(int j=0;j<i;j++){
            matrixI[j][j]=1;
        }
        return matrixI;
    }

    //the ith row of matrix times num
    public static void rowMultiply(double[][]matrix,int i,double num){
        if(matrix==null||i>=matrix.length||i<0){
            System.out.println("could not do row multiply");
            return;
        }
        for(int j=0;j<matrix[0].length;j++){
            matrix[i][j]*=num;
        }
    }

    public static void rowJminusrowI(double [][]matrix,int i,int j,double num){
        if(i<0||j<0||matrix==null||i>=matrix.length||j>=matrix.length){
            System.out.println("Could not exchangeRow");
            return;
        }
        for(int k=0;k<matrix[0].length;k++){
            matrix[j][k]-=num*matrix[i][k];
        }
    }

    public static double[][] reverse(double[][]matrix){
        if(matrix==null||matrix.length!=matrix[0].length||det(matrix)==0){
            System.out.println("not reversible");
            return null;
        }
        int len=matrix.length;
        double[][] matrixI=matrixI(len);
        for(int i=0;i<len;i++){
            if(matrix[i][i]==0.0){
                int j=i+1;
                while(matrix[j][i]==0.0)j++;
                exchangeRow(matrix,i,j);
                exchangeRow(matrixI,i,j);
            }
            rowMultiply(matrixI,i,1/matrix[i][i]);
            rowMultiply(matrix,i,1/matrix[i][i]);
            for(int j=i+1;j<len;j++){
                rowJminusrowI(matrixI,i,j,matrix[j][i]);
                rowJminusrowI(matrix,i,j,matrix[j][i]);
            }
        }
        for(int i=len-1;i>=0;i--){
            for(int j=i-1;j>=0;j--){
                rowJminusrowI(matrixI,i,j,matrix[j][i]);
                rowJminusrowI(matrix,i,j,matrix[j][i]);
            }
        }
        return matrixI;
    }

    public static double[][] transposition(double[][]matrix){
        if(matrix==null)return null;
        double[][]transposition=new double[matrix[0].length][matrix.length];
        for(int i=0;i<transposition.length;i++){
            for(int j=0;j<transposition[0].length;j++){
                transposition[i][j]=matrix[j][i];
            }
        }
        return transposition;
    }

    public static void negative(double[][]matrix){
        for(int i=0;i<matrix.length;i++){
            for(int j=0;j<matrix[0].length;j++){
                matrix[i][j]*=-1;
            }
        }
    }

    //int put m1Xn and m2Xn matrix, return (m1+m2)Xn matrix, matrixB is below matrixA
    public static double[][] jointByRow(double[][]matrixA,double[][]matrixB){
        if(matrixA==null||matrixB==null||matrixA[0].length!=matrixB[0].length){
            System.out.println("could not joint");
            return matrixA==null?(matrixB==null?null:matrixB):matrixA;
        }
        double[][]matrixC=new double[matrixA.length+matrixB.length][matrixA[0].length];
        for(int i=0;i<matrixA.length;i++){
            for(int j=0;j<matrixA[0].length;j++){
                matrixC[i][j]=matrixA[i][j];
            }
        }
        for(int i=0;i<matrixB.length;i++){
            for(int j=0;j<matrixB[0].length;j++){
                matrixC[i+matrixA.length][j]=matrixB[i][j];
            }
        }
        return matrixC;
    }

    public static double[] sumVector(int len){
        if(len<1){
            return null;
        }
        double[]vector1=new double[len];
        for(int i=0;i<len;i++){
            vector1[i]=1;
        }
        return vector1;
    }

    public static double innerProduct(double[]vector,double[]...optional){
        double[]vector2;
        if(optional.length==0){
            vector2=sumVector(vector.length);
        }else vector2=optional[0];
        double sum=0;
        for(int i=0;i<vector.length;i++){
            sum+=vector[i]*vector2[i];
        }
        return sum;
    }

    //int put mXn1 and mXn2 matrix, return mX(n1+n2) matrixC, matrixB is on the right of matrixA
    public static double[][] jointByCol(double[][]matrixA,double[][]matrixB){
        if(matrixA==null||matrixB==null||matrixA.length!=matrixB.length){
            System.out.println("could not joint");
            return matrixA==null?(matrixB==null?null:matrixB):matrixA;
        }
        double[][]matrixC=new double[matrixA.length][matrixA[0].length+matrixB[0].length];
        for(int i=0;i<matrixA.length;i++){
            for(int j=0;j<matrixA[0].length;j++){
                matrixC[i][j]=matrixA[i][j];
            }
        }
        for(int i=0;i<matrixB.length;i++){
            for(int j=0;j<matrixB[0].length;j++){
                matrixC[i][j+matrixA[0].length]=matrixB[i][j];
            }
        }
        return matrixC;
    }

    public static double vectorSquareSum(double[][]vector){
        if(vector==null||vector.length<1||vector[0].length!=1){
            System.out.println("input not a vector");
            return 0;
        }
        double norm=0.0;
        for(int i=0;i<vector.length;i++){
            norm+=vector[i][0]*vector[i][0];
        }
        return norm;
    }

    public static double[][] toUnitVector(double[][]vector){
        if(vector==null||vector.length<1||vector[0].length!=1){
            System.out.println("input not a vector");
            return null;
        }
        double norm=vectorSquareSum(vector);
        norm=Math.sqrt(norm);
        double[][]res=new double[vector.length][1];
        for(int i=0;i<vector.length;i++){
            res[i][0]=vector[i][0]/norm;
        }
        return res;
    }

    //if the input is a symmetric matrix it will return true.
    public static boolean isSymmetric(double[][]matrix){
        if(matrix==null||matrix.length!=matrix[0].length){
            return false;
        }
        for(int i=0;i<matrix.length-1;i++){
            for(int j=i;j<matrix.length;j++){
                if(matrix[i][j]!=matrix[j][i])return false;
            }
        }
        return true;
    }

    public static void niceToString(double[][]matrix){
        if(matrix==null){
            System.out.println(matrix);
            return;
        }else{
            for(int i=0;i<matrix.length;i++){
                System.out.println(Arrays.toString(matrix[i]));
            }
        }
    }

}
