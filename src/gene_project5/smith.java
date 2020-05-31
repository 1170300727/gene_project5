package gene_project5;

import java.util.Arrays;
import java.util.HashMap;

public class smith {
  private String s1; // the reference string
  private String s2; // the read
  private int penaltyInsertion = 7;
  private int penaltyDeletion = 7;
  private int score[][];
  private int maxScore;
  // 0--end(not include),1--up,2--left,3--leftup,4--end(include)
  private int back[][];
  private HashMap<String , Integer> penaltyMatrix = new HashMap<String , Integer>();
  
  public static void main(String[] args) {
    String s1 = "ACGAACGT";
    String s2 = "ACGACGT";
    String cigar;
    smith sm = new smith(s1, s2);
    sm.compute();
    sm.printScoreMatirx();
    System.out.println();
    sm.printBack();
    cigar = sm.getCigar();
    System.out.println(cigar);
  }
  
  public smith(String s1, String s2) {
    this.s1 = s1;
    this.s2 = s2;
    penaltyMatrix.put("AA", -5);
    penaltyMatrix.put("CC", -5);
    penaltyMatrix.put("GG", -5);
    penaltyMatrix.put("TT", -5);
    penaltyMatrix.put("AC", 4);
    penaltyMatrix.put("AG", 4);
    penaltyMatrix.put("AT", 4);
    penaltyMatrix.put("CA", 4);
    penaltyMatrix.put("CG", 4);
    penaltyMatrix.put("CT", 4);
    penaltyMatrix.put("TA", 4);
    penaltyMatrix.put("TC", 4);
    penaltyMatrix.put("TG", 4);
    penaltyMatrix.put("GA", 4);
    penaltyMatrix.put("GC", 4);
    penaltyMatrix.put("GT", 4);
  }
  
  public String getCigar() {
    int max = 0;
    int maxI = 0;
    int maxJ = 0;
    int m = s1.length();
    int n = s2.length();
    for (int i = 0;i < m;i++) { 
      for (int j = 0;j < n;j++) {
        if(score[i][j] > max) {
          max = score[i][j];
          maxI = i;
          maxJ = j;
        }
      }
    }
    String result = "";
    this.maxScore = max;
    while (back[maxI][maxJ] != 0 && back[maxI][maxJ] != 4) {
      if (back[maxI][maxJ] == 1) {
        result = result + "I";
        maxI--;
      } else if (back[maxI][maxJ] == 2) {
        result = result + "D";
        maxJ--;
      } else if (back[maxI][maxJ] == 3) {
        result = result + "M";
        maxI--;
        maxJ--;
      }
    }
    result = new StringBuffer(result).reverse().toString();
    return result;
  }
  
  public void compute() {
    int m = s1.length();
    int n = s2.length();
    score = new int[m][n];
    back = new int[m][n];
    int[] scores = new int[4];
    String connect;
    int tempScore;
    int max;
    for (int i = 0;i < n;i++) { //column
      for (int j = 0;j < m;j++) { //crow
        connect = String.valueOf(s1.charAt(j)) + String.valueOf(s2.charAt(i));
        tempScore = penaltyMatrix.get(connect);
        if (i == 0 && j == 0) { 
          if (tempScore > 0) {
            score[j][i] = tempScore;
            back[j][i] = 4;
          } else {
            score[j][i] = 0;
            back[j][i] = 0;
          }
        } else if (i == 0 && j != 0) { //the first column
          if (score[j - 1][i] - penaltyInsertion > 0) {
            score[j][i] = score[j - 1][i] - penaltyInsertion;
            back[j][i] = 4;
          } else {
            score[j][i] = 0;
            back[j][i] = 0;
          }
        } else if (i != 0 && j == 0) { //the first row
          if (score[j][i - 1] - penaltyDeletion > 0) {
            score[j][i] = score[j][i - 1] - penaltyDeletion;
            back[j][i] = 4;
          } else {
            score[j][i] = 0;
            back[j][i] = 0;
          }
        } else { //normal situation
          scores[0] = 0;
          scores[1] = score[j - 1][i] - penaltyInsertion;
          scores[2] = score[j][i - 1] - penaltyDeletion;
          scores[3] = score[j - 1][i - 1] - tempScore;
          Arrays.sort(scores); // ascending sort
          max = scores[3];
          score[j][i] = max;
          if (max == scores[0]) {
            back[j][i] = 0;
          } else if (max == scores[1]) {
            back[j][i] = 1;
          } else if (max == scores[2]) {
            back[j][i] = 2;
          } else if (max == scores[3]) {
            back[j][i] = 3;
          }
        }
      }
    }
  }
  
  public void printScoreMatirx() {
    int m = s1.length();
    int n = s2.length();
    for (int i = 0;i < m;i++) { 
      for (int j = 0;j < n;j++) {
        System.out.print(score[i][j] + "\t");
      }
      System.out.println();
    }
  }
  
  public void printBack() {
    int m = s1.length();
    int n = s2.length();
    for (int i = 0;i < m;i++) { 
      for (int j = 0;j < n;j++) {
        System.out.print(back[i][j] + "\t");
      }
      System.out.println();
    }
  }
  
  public int getMaxScore() {
    return this.maxScore;
  }
}
