package gene_project5;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class aligment {
  
  ArrayList<String> reads = new ArrayList<String>();
  ArrayList<String> chromRef = new ArrayList<String>();
  HashMap<String , value> hash = new HashMap<String , value>();
  
  public static void main(String[] args) {
    aligment al = new aligment();
    al.readFromFile();
    System.out.println(al.reads.size());
    //al.printReads();
    al.getHash();
    //al.randomAL();
    al.hashAL();
  }
  
  public void getHash() {
    int length;
    String s;
    String ref;
    int chromNum = chromRef.size();
    for (int j = 0;j < chromNum;j++) {
      ref = chromRef.get(j);
      length = ref.length();
      for (int i = 0;i < length - 11;i++) {
        s = ref.substring(i,i + 12);
        value va = new value(j, i);
        hash.put(s, va);
      }
    }
  }
  public void hashAL() {
    int length;
    String s;
    String read;
    String ref;
    int readLength = reads.get(0).length();
    ArrayList<alResult> results = new ArrayList<alResult>();
    int maxScore;
    int chromNum = chromRef.size();
    
    length = reads.size();
    for (int i = 0;i < length;i++) {
      read = reads.get(i);
      readLength = reads.get(i).length();
      alResult re = new alResult(read);    
      for (int j = 0;j < readLength - 11;j++) {
        String seed;
        seed = read.substring(j, j + 12);
        if (hash.get(seed) != null) { 
          value va = hash.get(seed);
          int k = va.chromNum;
          ref = chromRef.get(k);
          int position = va.position;
          if (position - j + readLength > ref.length() || position - j < 0) {
            continue;
          }
          s = ref.substring(position - j, position - j + readLength);
          smith sm = new smith(read, s);
          sm.compute();
          String cigar = sm.getCigar();
          maxScore = sm.getMaxScore();
          va = new value(k, position - j);
          samData da = new samData(va, maxScore, cigar);
          re.resultsForEachRead.add(da);
        } else {
          System.out.println(i + " " + j + seed);
        }
        
        
      }
      results.add(re);
    }
    
    try { File writeName = new File("lib/HashALoutput.txt"); // 相对路径，如果没有则要建立一个新的output.txt文件 
      writeName.createNewFile(); // 创建新文件,有同名的文件的话直接覆盖
      try (FileWriter writer = new FileWriter(writeName);
          BufferedWriter out = new BufferedWriter(writer)
      ) { 
        int size = results.size();
        for (int i = 0;i < size;i++) {
          alResult re = results.get(i);
          int sizeRe = re.resultsForEachRead.size();
          int maxIndex = 0;
          int maxscore = 0;
          for (int j = 0;j < sizeRe;j++) {
            if (re.resultsForEachRead.get(j).maxScore > maxscore) {
              maxscore = re.resultsForEachRead.get(j).maxScore;
              maxIndex = j;
            }
          }
          out.write("readID:" + i + "\t");
          out.write("chromID:" + re.resultsForEachRead.get(maxIndex).va.chromNum + "\t");
          out.write("position:" + re.resultsForEachRead.get(maxIndex).va.position + "\t");
          out.write("socre:" + re.resultsForEachRead.get(maxIndex).maxScore);
          out.write("cigar:" + re.resultsForEachRead.get(maxIndex).cigar + "\r\n");
        }
      
        out.flush(); // 把缓存区内容压入文件 
        } 
      } catch (IOException e) { 
        e.printStackTrace(); 
      }

  }
  
  public void randomAL() {
    int length;
    int randomPosition;
    String s;
    String read;
    String ref;
    int readLength = reads.get(0).length();
    ArrayList<alResult> results = new ArrayList<alResult>();
    int maxScore;
    int chromNum = chromRef.size();
    
    length = reads.size();
    for (int i = 0;i < length;i++) {
      read = reads.get(i);
      readLength = reads.get(i).length();
      read = read.substring(13, readLength);
      alResult re = new alResult(read);    
      for (int j = 0;j < 10;j++) {
        int k = (int)(0+Math.random()*(chromNum - 1 - 0 + 1));
        ref = chromRef.get(k);
        randomPosition = (int)(0+Math.random()*(length - 1 - 0 + 1));
        if (randomPosition + readLength > ref.length()) {
          randomPosition = ref.length() - readLength;
        }
        s = ref.substring(randomPosition, randomPosition + readLength);
        smith sm = new smith(read, s);
        sm.compute();
        String cigar = sm.getCigar();
        maxScore = sm.getMaxScore();
        value va = new value(k, randomPosition);
        samData da = new samData(va, maxScore, cigar);
        re.resultsForEachRead.add(da);
      }
      results.add(re);
    }
    
    try { File writeName = new File("lib/output.txt"); // 相对路径，如果没有则要建立一个新的output.txt文件 
      writeName.createNewFile(); // 创建新文件,有同名的文件的话直接覆盖
      try (FileWriter writer = new FileWriter(writeName);
          BufferedWriter out = new BufferedWriter(writer)
      ) { 
        int size = results.size();
        for (int i = 0;i < size;i++) {
          alResult re = results.get(i);
          int sizeRe = re.resultsForEachRead.size();
          int maxIndex = 0;
          int maxscore = 0;
          for (int j = 0;j < sizeRe;j++) {
            if (re.resultsForEachRead.get(j).maxScore > maxscore) {
              maxscore = re.resultsForEachRead.get(j).maxScore;
              maxIndex = j;
            }
          }
          out.write("readID:" + i + "\t");
          out.write("chromID:" + re.resultsForEachRead.get(maxIndex).va.chromNum + "\t");
          out.write("position:" + re.resultsForEachRead.get(maxIndex).va.position + "\t");
          out.write("socre:" + re.resultsForEachRead.get(maxIndex).maxScore);
          out.write("cigar:" + re.resultsForEachRead.get(maxIndex).cigar + "\r\n");
        }
      
        out.flush(); // 把缓存区内容压入文件 
        } 
      } catch (IOException e) { 
        e.printStackTrace(); 
      }

  }
  
  public void readFromFile() {
    File f=new File("lib/DRR076693_1.1k.fastq");
    File f2=new File("lib/sacCer3.fa");
        FileReader fre;
        FileReader fre2;
        try {
          fre = new FileReader(f);
          BufferedReader bre=new BufferedReader(fre);
          fre2 = new FileReader(f2);
          BufferedReader bre2=new BufferedReader(fre2);
          String str="";
          int i = 0;
          while((str=bre.readLine())!=null) //●判断最后一行不存在，为空
          {
            if(i % 4 == 1) {
              reads.add(str);
            }
            i++;
          }
          i = 0;
          StringBuilder ref = new StringBuilder();
          while((str=bre2.readLine())!=null) //●判断最后一行不存在，为空
          {
            if (str.contains(">chr") && i == 0) {
              ref = new StringBuilder();
            } else if (str.contains(">chr") && i != 0) {
              chromRef.add(ref.toString());
              System.out.println(ref.length());
              ref = new StringBuilder();        
            } else {
              ref.append(str);
            }
            i++;
          }
          
          bre.close();
          fre.close();
          bre2.close();
          fre2.close();
        } catch (FileNotFoundException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }      
  }
  
  public void printReads() {
    int length = reads.size();
    for (int i = 0;i < length;i++) {
      System.out.println(reads.get(i));
      if (i == 20) {
        break;
      }
    }
  }
}

class alResult {
  String read;
  ArrayList<samData> resultsForEachRead = new ArrayList<samData>();
  
  public alResult(String read) {
    this.read = read;
  }
}

class samData {
  value va;
  int maxScore;
  String cigar;
  
  public samData(value va, int max, String ci) {
    this.va = va;
    this.maxScore = max;
    this.cigar = ci;
  }
}

class value {
  int chromNum;
  int position;
  
  public value(int chromNum, int position) {
    this.chromNum = chromNum;
    this.position = position;
  }
}