package ttp;

import newrep.Item;
import utils.Deb;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class TTPSolution {
    private int[] tour;
    private int[] pickingPlan;

    public boolean isTspImprove() {
        return tspImprove;
    }

    public void setTspImprove(boolean tspImprove) {
        this.tspImprove = tspImprove;
    }

    public boolean tspImprove;

    public long[] getL() {
        return L;
    }

    public void setL(long[] l) {
        L = l;
    }

    public long[] L;
    public Double[] getScoresArr() {
        return scoresArr;
    }

    public void setScoresArr(Double[] scoresArr) {
        this.scoresArr = scoresArr;
    }

    public Double[] scoresArr;
    public ArrayList<Item> getBoundaryItemsArr() {
        return boundaryItemsArr;
    }

    public void setBoundaryItemsArr(ArrayList<Item> boundaryItemsArr) {
        this.boundaryItemsArr = boundaryItemsArr;
    }

    public ArrayList<Item> boundaryItemsArr;

    public float[] getLeastPickedPR() {
        return leastPickedPR;
    }

    public void setLeastPickedPR(float[] leastPickedPR) {
        this.leastPickedPR = leastPickedPR;
    }

    public float[] getHighestUnpickedPR() {
        return highestUnpickedPR;
    }

    public void setHighestUnpickedPR(float[] highestUnpickedPR) {
        this.highestUnpickedPR = highestUnpickedPR;
    }

    public float[] getPrefixMinimumValue() {
        return prefixMinimumValue;
    }

    public void setPrefixMinimumValue(float[] prefixMinimumValue) {
        this.prefixMinimumValue = prefixMinimumValue;
    }

    public float[] getPostfixMaximumValue() {
        return postfixMaximumValue;
    }

    public void setPostfixMaximumValue(float[] postfixMaximumValue) {
        this.postfixMaximumValue = postfixMaximumValue;
    }

    private float[] leastPickedPR;
    private float[] highestUnpickedPR;
    private float[] prefixMinimumValue;
    private float[] postfixMaximumValue;

    public long fp;
    public double ft;
    public double ob;
    public long wend;

    public double standScore;
    public double averageWeight;

    // time accumulator
    public double[] timeAcc;
    // time record
    public double[] timeRec;
    // weight accumulator
    public long[] weightAcc;
    // weight record at each iteration
    public long[] weightRec;
    // tour mapper
    public int[] mapCI;


    //S5部分
    public int[] tspTour;
    public int[] packingPlan;
    public double finalProfit   = Double.NEGATIVE_INFINITY;
    public double finalTravelTime   = Double.POSITIVE_INFINITY;
    public long finalDistance   = Long.MAX_VALUE;
    public double objectiveScore   = Double.NEGATIVE_INFINITY;
    public double finalCapacityFree = Double.POSITIVE_INFINITY;
    public double finalWeight = Double.POSITIVE_INFINITY;
    public long computationTime = Long.MAX_VALUE;

    private void initSolution(int[] tour, int[] pickingPlan) {
        this.tour = tour;
        this.pickingPlan = pickingPlan;

        this.leastPickedPR = new float[tour.length];
        this.highestUnpickedPR = new float[tour.length];
        this.prefixMinimumValue = new float[tour.length];
        this.postfixMaximumValue = new float[tour.length];
        // records
        this.timeAcc = new double[this.tour.length];
        this.timeRec = new double[this.tour.length];
        this.weightAcc = new long[this.tour.length];
        this.weightRec = new long[this.tour.length];
        this.mapCI = new int[this.tour.length];
        this.boundaryItemsArr = new ArrayList<Item>();

        //S5部分
        this.tspTour = this.tour;
        this.packingPlan = this.pickingPlan;
    }

    public TTPSolution() {

    }

    public TTPSolution(int[] tour, int[] pickingPlan) {
        initSolution(tour, pickingPlan);
    }

    public TTPSolution(int m, int n) {

        this.tour = new int[m];
        this.pickingPlan = new int[n];

        this.leastPickedPR = new float[m];
        this.highestUnpickedPR = new float[m];
        this.prefixMinimumValue = new float[m];
        this.postfixMaximumValue = new float[m];
        this.boundaryItemsArr = new ArrayList<Item>();
        // records
        timeAcc = new double[tour.length];
        timeRec = new double[tour.length];
        weightAcc = new long[tour.length];
        weightRec = new long[tour.length];
        mapCI = new int[tour.length];
    }

    public TTPSolution(TTPSolution s2) {

        this.tour = Arrays.copyOf(s2.tour, s2.tour.length);
        this.pickingPlan = Arrays.copyOf(s2.pickingPlan, s2.pickingPlan.length);
        this.scoresArr = Arrays.copyOf(s2.scoresArr, s2.scoresArr.length);
        this.L = Arrays.copyOf(s2.L, s2.L.length);
        this.leastPickedPR = Arrays.copyOf(s2.leastPickedPR, s2.leastPickedPR.length);
        this.highestUnpickedPR = Arrays.copyOf(s2.highestUnpickedPR, s2.highestUnpickedPR.length);
        this.prefixMinimumValue = Arrays.copyOf(s2.prefixMinimumValue, s2.prefixMinimumValue.length);
        this.postfixMaximumValue = Arrays.copyOf(s2.postfixMaximumValue, s2.postfixMaximumValue.length);
        ArrayList<Item> list1 = new ArrayList<Item>(s2.boundaryItemsArr);
        this.setBoundaryItemsArr(list1);



        this.fp = s2.fp;
        this.ft = s2.ft;
        this.ob = s2.ob;
        this.objectiveScore = s2.objectiveScore;
        this.wend = s2.wend;
        this.standScore = s2.standScore;
        this.averageWeight = s2.averageWeight;

        this.timeAcc = Arrays.copyOf(s2.timeAcc, s2.timeAcc.length);
        this.timeRec = Arrays.copyOf(s2.timeRec, s2.timeRec.length);
        this.weightAcc = Arrays.copyOf(s2.weightAcc, s2.weightAcc.length);
        this.weightRec = Arrays.copyOf(s2.weightRec, s2.weightRec.length);
        this.mapCI = Arrays.copyOf(s2.mapCI, s2.mapCI.length);
    }

    public TTPSolution(String filePath) {
        File solFile = new File(filePath);
        BufferedReader br = null;

        int nbCities = 0, nbItems = 0;

        try {
            br = new BufferedReader(new FileReader(solFile));
            String line;

            // scan tour
            while ((line = br.readLine()) != null) {

                // number of cities
                if (line.startsWith("DIMENSION")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    nbCities = Integer.parseInt(line);
                }

                // number of items
                if (line.startsWith("NUMBER OF ITEMS")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    nbItems = Integer.parseInt(line);
                }

                if (line.startsWith("TOUR_SECTION")) {
                    this.tour = new int[nbCities];
                    for (int j = 0; j < nbCities; j++) {
                        line = br.readLine();
                        tour[j] = Integer.parseInt(line);
                    }
                }

                if (line.startsWith("PP_SECTION")) {
                    this.pickingPlan = new int[nbItems];
                    for (int j = 0; j < nbItems; j++) {
                        line = br.readLine();
                        pickingPlan[j] = Integer.parseInt(line);
                    }
                }
            } // end while

            br.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        this.initSolution(tour, pickingPlan);
    }

    @Override
    public String toString() {
        // the tour
        String s = "tsp tour    : (";
        for (int i = 0; i < tour.length; i++) {
            s += tour[i] + " ";
        }
        s = s.substring(0, s.length() - 1) + ")\n";

        // the picking plan
        s += "picking plan: (";
        for (int i = 0; i < pickingPlan.length; i++) {
            int pp = pickingPlan[i];
            s += pp + " ";
        }
        s = s.substring(0, s.length() - 1) + ")";

        return s;
    }

    @Override
    public TTPSolution clone() {
        return new TTPSolution(this);
    }

    @Override
    public boolean equals(Object o2) {

        TTPSolution s2 = (TTPSolution) o2;

        for (int i = 0; i < this.tour.length; i++) {
            if (this.tour[i] != s2.tour[i]) return false;
        }
        for (int i = 0; i < this.pickingPlan.length; i++) {
            if (this.pickingPlan[i] != s2.pickingPlan[i]) return false;
        }

        return true;
    }

    // getters
    public int[] getTour() {
        return tour;
    }

    public int[] getPickingPlan() {
        return pickingPlan;
    }

    // setters
    public void setTour(int[] tour) {
        this.tour = tour;
    }

    public void setPickingPlan(int[] pickingPlan) {
        this.pickingPlan = pickingPlan;
    }


    public String output() {
        String s =
                "DIMENSION : " + tour.length + "\n" +
                        "NUMBER OF ITEMS : " + pickingPlan.length + "\n" +
                        "\n";

        s += "TOUR_SECTION\n";
        for (int x : tour) {
            s += x + "\n";
        }
        s += "\n";

        s += "PP_SECTION\n";
        for (int x : pickingPlan) {
            s += x + "\n";
        }

        s += "EOF";

        return s;
    }

    public void printStats() {

        Deb.echo("============");
        Deb.echo(" STATISTICS ");
        Deb.echo("============");
        Deb.echo("objective   : " + this.ob);
        Deb.echo("final time  : " + this.fp);
        Deb.echo("final weight: " + this.wend);
        Deb.echo("final profit: " + this.fp);

        int cmpItems = 0;
        int nbItems = this.pickingPlan.length;
        for (int x : this.pickingPlan) if (x != 0) cmpItems++;
        Deb.echo("percent inserted: " + cmpItems + "/" + nbItems + "(" +
                String.format("%.2f", (cmpItems * 100.0) / nbItems) + "%)");
        Deb.echo("============");
    }




    //S5部分
    public double getObjective() {
        return this.objectiveScore;
    }

    public void reset() {
        finalProfit = Double.NEGATIVE_INFINITY;
        finalTravelTime = Double.POSITIVE_INFINITY;
        finalDistance = Long.MAX_VALUE;
        objectiveScore = Double.NEGATIVE_INFINITY;
        finalCapacityFree = Double.POSITIVE_INFINITY;
        finalWeight = Double.POSITIVE_INFINITY;
        computationTime = Long.MAX_VALUE;
    }
    public String answer() {
        int[] tourOut = new int[tspTour.length - 1];
        for (int i = 0; i < tspTour.length - 1; i++){
            tourOut[i] = tspTour[i] + 1;
        }

        int itemsPerCity = packingPlan.length / (tspTour.length  -2);

        List<Integer> packingPlanList = new ArrayList<Integer>();
        int packingPlanIndex = 0;
        for (int i = 1; i < tspTour.length - 1; i++){
            int city = tspTour[i];
            for (int j = 0; j < itemsPerCity; j++){
                if (packingPlan[packingPlanIndex] == 1){
                    int itemIndex = j * (tspTour.length  -2) + city - 1;
                    int itemIndexFrom1 = itemIndex + 1;
                    packingPlanList.add(itemIndexFrom1 );
                }
                packingPlanIndex++;
            }
        }
        Collections.sort(packingPlanList);

        int[] packingOut = new int[packingPlanList.size()];
        for (int i = 0; i < packingPlanList.size(); i++){
            packingOut[i] = packingPlanList.get(i);
        }
        return Arrays.toString(tourOut) + "\r\n" + Arrays.toString(packingOut) +  "\r\n";
    }

    public void writeResult(String title){
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(title, false));
            writer.write(this.answer());
            writer.flush();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();

        }
    }


    // the following print functions print varying levels of detail
    public void print() {
        System.out.print(finalCapacityFree+" "+finalWeight+" "+finalProfit+" "+finalDistance+" "+finalTravelTime+" "+objectiveScore +" "+computationTime);
    }
    public void println() {
        this.print();
        System.out.println();
    }
    public void altPrint(){
        System.out.println();
        System.out.println(String.format("%-30s= %.2f" , "Leftover Knapsack Capacity", finalCapacityFree));
        System.out.println(String.format("%-30s= %.2f" , "Total weight of items", finalWeight));
        System.out.println(String.format("%-30s= %.2f" , "Total profit", finalProfit));
        System.out.println(String.format("%-30s= %d" , "Total time (raw)", finalDistance));
        System.out.println(String.format("%-30s= %.2f" , "Total time", finalTravelTime));
        System.out.println(String.format("%-30s= %.2f" , "ob", objectiveScore));
        System.out.println(String.format("%-30s= %d" , "Computation time (ms)", computationTime));
    }
    public void printFull() {
        this.println();
        System.out.println("tspTour "+Arrays.toString(tspTour));
        System.out.println("packingPlan "+Arrays.toString(packingPlan));
    }
    public void printFullForCCode() {
        System.out.print("[");
        for (int i=0; i<tspTour.length-2; i++) {
            System.out.print((tspTour[i]+1)+",");
        }
        System.out.println((tspTour[tspTour.length-2]+1)+"]");

        System.out.print("[");
        for (int i=0; i<packingPlan.length-1; i++) {
            if (packingPlan[i]==1) {
                System.out.print((i+1)+",");
            }
        }
        if (packingPlan[packingPlan.length-1]==1) {
            System.out.println(packingPlan.length+"]");
        } else {
            System.out.println("]");
        }
    }
}
