package com.graduate.util;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FileOperator {

    //read the .dat file by line
    public static String[] readFromDatFile(String fileName){
        String root = System.getProperty("user.dir");
        fileName=root+ File.separator+fileName;
        List<String> res=new ArrayList<>();
        try {
            FileReader fileReader=new FileReader(fileName);
            BufferedReader bufferedReader=new BufferedReader(fileReader);
            String line;
            while((line=bufferedReader.readLine())!=null){
                res.add(line);
            }
            fileReader.close();
        }catch (Exception e){
            e.printStackTrace();
        }
        return res.toArray(new String[res.size()]);
    }

    //convert the string[]file to double[][]matrix, each line of the file should have the same amount of number
    public static double[][] dealTheData(String[]file){
        if(file==null||file.length==0){
            return null;
        }
        List<Double>helper=new ArrayList<>();
        String s=file[0];
        while(s.length()>0){
            s=s.trim();
            int i=s.indexOf(' ');
            if(i==-1){
                helper.add(new Double(s));
                break;
            }
            helper.add(new Double(s.substring(0,i)));
            s=s.substring(i);
        }
        double[][]res=new double[file.length][helper.size()];
        for(int i=0;i<file.length;i++){
            s=file[i];
            for(int j=0;j<helper.size();j++){
                s=s.trim();
                int k=s.indexOf(' ');
                if(k==-1){
                    res[i][j]=new Double(s);
                }else{
                    res[i][j]=new Double(s.substring(0,k));
                    s=s.substring(k);
                }
            }
        }
        return res;
    }

    public static File makeFile(String fileName, String[] data){
        try{
            File file=new File(fileName);
            FileWriter fileWriter=new FileWriter(fileName);
            BufferedWriter bufferedWriter=new BufferedWriter(fileWriter);
            for(int i=0;i<data.length;i++){
                bufferedWriter.write(data[i]);
            }
            bufferedWriter.close();
            return file;
        }catch (Exception e){
            e.printStackTrace();
        }
        return null;
    }

}
