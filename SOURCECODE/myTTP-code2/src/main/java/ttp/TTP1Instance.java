package ttp;

import newrep.Item;
import utils.CityCoordinates;
import utils.ConfigHelper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

public class TTP1Instance extends TTPInstance {
    /**
     * knapsack renting ratio per time unit
     */
    protected double rent;

    /**
     * instance to string
     */
    @Override
    public String toString() {
        return super.toString() + "rent     : $" + this.rent;
    }

    /**
     * reads the instance from a .ttp file using file's name
     *
     * @param fileName
     */
    public TTP1Instance(String fileName) {
        //this.name = fileName;
        String[] sp = fileName.split("/", 2);
        this.directory = sp[0];
        this.name = sp[1];
        String[] sp2 = directory.split("-", 2);
        this.tspName = sp2[0];

        String ttpData = ConfigHelper.getProperty("ttpdata");

        this.ttpFile = new File(ttpData + fileName);
        this.file = ttpFile;
        BufferedReader br = null;                       //BufferedReader用于加快读取字符的速度

        try {
            br = new BufferedReader(new FileReader(this.ttpFile));
            String line;

            while ((line = br.readLine()) != null) {

                // instance name
                if (line.startsWith("PROBLEM NAME")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    //this.name = line;
                    this.problemName = line;
                }

                // KP data type
                if (line.startsWith("KNAPSACK DATA TYPE")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    this.knapsackDataType = line;
                }

                // number of cities
                if (line.startsWith("DIMENSION")) {
                    // if (line.startsWith("NUMBER OF NODES")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    this.nbCities = Integer.parseInt(line);
                    this.numberOfNodes = this.nbCities;
                }

                // number of items
                if (line.startsWith("NUMBER OF ITEMS")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    this.nbItems = Integer.parseInt(line);
                    this.numberOfItems = this.nbItems;
                }

                // knapsack capacity
                if (line.startsWith("CAPACITY OF KNAPSACK")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    this.capacity = Long.parseLong(line);
                    this.capacityOfKnapsack = this.capacity;
                }

                // minimum velocity
                if (line.startsWith("MIN SPEED")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    this.minSpeed = Double.parseDouble(line);
                }

                // maximum velocity
                if (line.startsWith("MAX SPEED")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    this.maxSpeed = Double.parseDouble(line);
                }

                // rent
                if (line.startsWith("RENTING RATIO")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    this.rent = Double.parseDouble(line);
                    this.rentingRatio = this.rent;
                }

                // edge weight
                if (line.startsWith("EDGE_WEIGHT_TYPE")) {
                    line = line.substring(line.indexOf(":") + 1);
                    line = line.replaceAll("\\s+", "");
                    this.edgeWeightType = line;
                }

                // nodes
                if (line.startsWith("NODE_COORD_SECTION")) {
                    this.items = new int[this.numberOfItems][4];
                    // coordinates
                    this.coordinates = new CityCoordinates[this.nbCities];

                    this.nodes = new double[this.numberOfNodes][3];//S5部分
                    for (int i = 0; i < this.nbCities; i++) {
                        line = br.readLine();
                        String[] parts = line.split("\\s+");
                        this.coordinates[i] = new CityCoordinates(Double.parseDouble(parts[1]), Double.parseDouble(parts[2]));
                        //S5部分

                        for (int j = 0; j < parts.length; j++) {
                            double temp = Double.parseDouble(parts[j]);
                            // adjust city number by 1
                            if (j == 0) temp = temp - 1;
                            this.nodes[i][j] = temp;
                        }
                    }

                    // distance matrix
                    if (nbCities < 10000) {
                        this.setDist(new long[this.nbCities][this.nbCities]);
                        for (int i = 0; i < nbCities; i++) {
                            for (int j = 0; j < nbCities; j++) {
//                getDist()[i][j] = Math.round(this.coordinates[i].distanceEuclid(this.coordinates[j]));
                                getDist()[i][j] = (long) Math.ceil(this.coordinates[i].distanceEuclid(this.coordinates[j]));
                                //System.out.println(this.coord[i] + "&" + this.coord[j] + "->" + dist[i][j]);
                            }
                        }
                    }
                }

                // items
                if (line.startsWith("ITEMS SECTION")) {
                    this.profits = new int[this.nbItems];
                    this.weights = new int[this.nbItems];
                    this.availability = new int[this.nbItems];

                    this.items = new int[this.numberOfItems][4];//S5
                    int itemsPerCity = this.numberOfItems / (this.numberOfNodes-1);//S5
                    this.profitabilityRatio = new double[this.numberOfNodes][itemsPerCity];//S5

                    for (int i = 0; i < this.nbItems; i++) {
                        line = br.readLine();
                        String[] splittedLine = line.split("\\s+");

                        this.profits[i] = Integer.parseInt(splittedLine[1]);
                        this.weights[i] = Integer.parseInt(splittedLine[2]);
                        this.availability[i] = Integer.parseInt(splittedLine[3]);



                        //S5部分

                        for (int j = 0; j < splittedLine.length; j++) {
                            int temp = Integer.parseInt(splittedLine[j]);
                            // adjust city number by 1
                            if (j==0) temp =  temp - 1;  // item numbers start here with 0 --> in TTP files with 1
                            if (j==3) {
                                int k=0;
                                while (this.profitabilityRatio[temp-1][k] > 0.000001){
                                    k++;
//                                    if (k==7){
//                                        System.out.println("222");
//                                    }
                                }
                                this.profitabilityRatio[temp-1][k] = (float)this.items[i][1] / (float)this.items[i][2];

                                temp =  temp - 1;  // city numbers start here with 0 --> in TTP files with 1
                            }
                            this.items[i][j] = temp;
                        }
                    }

                }


            } // end while

            br.close();

        } catch (IOException ex) {
            ex.printStackTrace();
        }
        clusterItems();               //将物品聚类到相应的城市中
    }


    public double getRent() {
        return rent;
    }


    /**
     * objective function
     *
     * @param s the TTP solution
     */
    public void objective(TTPSolution s) {

        int[] x = s.getTour(); //city Array
        int[] z = s.getPickingPlan();

        long[][] D = getDist();
        double C = (maxSpeed - minSpeed) / capacity; // velocity const
        double velocity;

        long acc;       // iteration weight accumulator 迭代权累加器
        long wc = 0;    // current weight
        long fp = 0;    // final profit
        double ft = 0;  // tour time
        double ob;      // objective value


        // visit all cities
        for (int i = 0; i < this.nbCities; i++) {
            acc = 0;
            // check only items contained in current city
            for (int j : clusters[x[i] - 1]) { //根据tour处在当前位置城市中的物品下标:j
                if (z[j] != 0) {
                    fp += profits[j];
                    acc += weights[j];
                }
            }

            wc += acc;
            velocity = maxSpeed - wc * C;

            int h = (i + 1) % nbCities; //next城市
            ft += distFor(x[i] - 1, x[h] - 1) / velocity; //去往下个城市需要花费的时间

            // record important data for future use
            s.timeAcc[i] = ft;//until now time
            s.timeRec[i] = distFor(x[i] - 1, x[h] - 1) / velocity;//当前城市到下一个城市花费的时间
            s.weightAcc[i] = wc;//到当前城市时累积的重量
            s.weightRec[i] = acc;//在当前城市拾取累积的重量

            // map indices to their associated cities 将城市和索引对调
            s.mapCI[x[i] - 1] = i;//当前城市x[i]在tour中的下标为i
        }

        ob = fp - ft * rent;

        // solution properties
        s.fp = fp;
        s.ft = ft;
        s.wend = capacity - wc;
        s.ob = ob;
        s.objectiveScore=ob;
    }

    public void calculateScoreObjective(TTPSolution s,TTP1Instance ttp) {

        int[] x = s.getTour(); //city Array
        int[] z = s.getPickingPlan();

        long[][] D = getDist();
        double C = (maxSpeed - minSpeed) / capacity; // velocity const
        double velocity;

        long acc;       // iteration weight accumulator 迭代权累加器
        long wc = 0;    // current weight
        long fp = 0;    // final profit
        double ft = 0;  // tour time
        double ob;      // objective value


        float prefixMinimumValue_inCity = (float) 0.0;
        float postfixMaximumValue_inCity = (float) 10000.0;
        float[] LPPR = s.getLeastPickedPR();;
        float[] HUPR = s.getHighestUnpickedPR();


        // visit all cities
        for (int i = 0; i < this.nbCities; i++) {
            acc = 0;
            // check only items contained in current city
            for (int j : clusters[x[i] - 1]) { //根据tour处在当前位置城市中的物品下标:j
                int itemIndex = j;
                if (z[j] != 0) {
                    fp += profits[j];
                    acc += weights[j];
                }
            }


            //原本代码
            wc += acc;
            velocity = maxSpeed - wc * C;

            int h = (i + 1) % nbCities; //next城市
            ft += distFor(x[i] - 1, x[h] - 1) / velocity; //去往下个城市需要花费的时间

            // record important data for future use
            s.timeAcc[i] = ft;//until now time
            s.timeRec[i] = distFor(x[i] - 1, x[h] - 1) / velocity;//当前城市到下一个城市花费的时间
            s.weightAcc[i] = wc;//到当前城市时累积的重量
            s.weightRec[i] = acc;//在当前城市拾取累积的重量

            // map indices to their associated cities 将城市和索引对调
            s.mapCI[x[i] - 1] = i;//当前城市x[i]在tour中的下标为i

        }
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        int[] A = ttp.getAvailability();
        int[] tour = s.getTour();
        double maxSpeed = ttp.getMaxSpeed();
        // neighbor solution
        int origBF,k,i;

        // distances of all tour cities (city -> end)
        long[] L = new long[nbCities];

        // time approximations
        double t1;
        L[nbCities-1] = ttp.distFor(tour[nbCities-1] - 1,0);
        for (i=nbCities-2; i >= 0; i--) {
            L[i] = L[i+1] + ttp.distFor(tour[i+1]-1, tour[i]-1);
        }
        for (k = 0; k < nbItems; k++) {
            // index where Bit-Flip happened
            origBF = s.mapCI[A[k] - 1];  //计算当前item的城市index
            // calculate time approximations
            t1 = L[origBF]*(1/(maxSpeed-C*ttp.weightOf(k)) - 1/maxSpeed);//计算当前物品选取与否的时间差
            // affect score to item
            s.scoresArr[k] = (Math.pow(ttp.profitOf(k),7.4)*(L[0]-L[origBF])/maxSpeed)/(Math.pow(ttp.weightOf(k),7.4)*t1);
        }
        double leastPickPR_inCity = Float.POSITIVE_INFINITY;
        double highestUnPickPR_inCity = (float) 0.0;
        boolean noItemsSelectedInCity = true;
        for ( i = 0; i < this.nbCities; i++) {
            for (int j : clusters[x[i] - 1]) { //根据tour处在当前位置城市中的物品下标:j
                int itemIndex = j;
                double PR =  s.scoresArr[itemIndex];
                if (z[j] != 0) {
                    noItemsSelectedInCity = false;
                    if (PR < leastPickPR_inCity) {
                        leastPickPR_inCity = PR;
                    }
                }else {
                    if (PR > highestUnPickPR_inCity) {
                        highestUnPickPR_inCity = PR;
                    }
                }
            }
            if (noItemsSelectedInCity) {
                leastPickPR_inCity = highestUnPickPR_inCity + 1;
            }
            //在当前城市中对应的最少拾取的价值和最大未拾取的价值
            LPPR[x[i] - 1] = (float) leastPickPR_inCity;
            HUPR[x[i] - 1] = (float) highestUnPickPR_inCity;
        }
        s.setLeastPickedPR(LPPR);
        s.setHighestUnpickedPR(HUPR);
//找到最小拾取和最大未拾取的preMin和postMax
            float[] preMinV = s.getPrefixMinimumValue();
            float[] posMaxV = s.getPostfixMaximumValue();
            float startPreMin = LPPR[0];
            for ( k = 0; k < preMinV.length; k++) {
                if (LPPR[k] < startPreMin) {
                    startPreMin = LPPR[k];
                }
                preMinV[k] = startPreMin;
            }
            float startPostMax = HUPR[posMaxV.length - 1];
            for (int l = posMaxV.length - 1; l > -1; l--) {
                if (HUPR[l] > startPostMax) {
                    startPostMax = HUPR[l];
                }
                posMaxV[l] = startPostMax;
            }


            s.setPrefixMinimumValue(preMinV);
            s.setPostfixMaximumValue(posMaxV);


        //原来代码
        ob = fp - ft * rent;

        // solution properties
        s.fp = fp;
        s.ft = ft;
        s.wend = capacity - wc;
        s.ob = ob;
        s.objectiveScore=ob;
    }

    public void CocoObjective(TTPSolution s) {

        int[] x = s.getTour(); //city Array
        int[] z = s.getPickingPlan();

        long[][] D = getDist();
        double C = (maxSpeed - minSpeed) / capacity; // velocity const
        double velocity;

        long acc;       // iteration weight accumulator 迭代权累加器
        long wc = 0;    // current weight
        long fp = 0;    // final profit
        double ft = 0;  // tour time
        double ob;      // objective value


        float prefixMinimumValue_inCity = (float) 0.0;
        float postfixMaximumValue_inCity = (float) 10000.0;
        float[] LPPR = s.getLeastPickedPR();;
        float[] HUPR = s.getHighestUnpickedPR();

        // visit all cities
        for (int i = 0; i < this.nbCities; i++) {
            acc = 0;
            float leastPickPR_inCity = (float) 100000.0;
            float highestUnPickPR_inCity = (float) 0.0;

//            boolean noItemsSelectedInCity = true;
//            for (int k = 0; k < (PR_clusters[x[i] - 1]).size(); k++) {
//                float PR = PR_clusters[x[i] - 1].get(k);
//                int itemIndex = clusters[x[i] - 1].get(k);
//                if (z[itemIndex] != 0) {//当前物品被拾取
//                    noItemsSelectedInCity = false;
//                    if (PR < leastPickPR_inCity) {
//                        leastPickPR_inCity = PR;
//                    }
//                } else {
//                    if (PR > highestUnPickPR_inCity) {
//                        highestUnPickPR_inCity = PR;
//                    }
//                }
//            }
//            if (noItemsSelectedInCity) {
//                leastPickPR_inCity = highestUnPickPR_inCity + 1;
//            }



            // check only items contained in current city
            for (int j : clusters[x[i] - 1]) { //根据tour处在当前位置城市中的物品下标:j
                if (z[j] != 0) {//当前物品被拾取
                    fp += profits[j];
                    acc += weights[j];
                }
            }

                wc += acc;
                velocity = maxSpeed - wc * C;

                int h = (i + 1) % nbCities; //next城市
                ft += distFor(x[i] - 1, x[h] - 1) / velocity; //去往下个城市需要花费的时间

                // record important data for future use
                s.timeAcc[i] = ft;//until now time
                s.timeRec[i] = distFor(x[i] - 1, x[h] - 1) / velocity;//当前城市到下一个城市花费的时间
                s.weightAcc[i] = wc;//到当前城市时累积的重量
                s.weightRec[i] = acc;//在当前城市拾取累积的重量

                // map indices to their associated cities 将城市和索引对调
                s.mapCI[x[i] - 1] = i;//当前城市x[i]在tour中的下标为i


            //在当前城市中对应的最少拾取的价值和最大未拾取的价值
//            LPPR[x[i] - 1] = leastPickPR_inCity;
//            HUPR[x[i] - 1] = highestUnPickPR_inCity;
            }
//            s.setLeastPickedPR(LPPR);
//            s.setHighestUnpickedPR(HUPR);
//
//            float[] preMinV = s.getPrefixMinimumValue();
//            float[] posMaxV = s.getPostfixMaximumValue();
//            float startPreMin = LPPR[0];
//            for (int k = 0; k < preMinV.length; k++) {
//                if (LPPR[k] < startPreMin) {
//                    startPreMin = LPPR[k];
//                }
//                preMinV[k] = startPreMin;
//            }
//            float startPostMax = HUPR[posMaxV.length - 1];
//            for (int l = posMaxV.length - 1; l > -1; l--) {
//                if (HUPR[l] > startPostMax) {
//                    startPostMax = HUPR[l];
//                }
//                posMaxV[l] = startPostMax;
//            }
//
//
//            s.setPrefixMinimumValue(preMinV);
//            s.setPostfixMaximumValue(posMaxV);

            ob = fp - ft * rent;

            // solution properties
            s.fp = fp;
            s.ft = ft;
            s.wend = capacity - wc;
            s.ob = ob;
            s.objectiveScore = ob;
    }

    public void MarkUncheckedAll(TTPSolution s) {
        ArrayList<Item> boundaryItemsArr = s.boundaryItemsArr;
        Iterator<Item> iter = boundaryItemsArr.iterator();
        while (iter.hasNext()){
            Item item = iter.next();
            item.isChecked = false;
        }
    }

    public void MarkChecked(TTPSolution s,int itemIndex) {
        ArrayList<Item> boundaryItemsArr = s.boundaryItemsArr;
        Iterator<Item> iter = boundaryItemsArr.iterator();
        while (iter.hasNext()){
            Item item = iter.next();
            if (item.itemId == itemIndex){
                item.isChecked = true;
            }
        }
    }

    public boolean AllChecked(TTPSolution s){
        ArrayList<Item> boundaryItemsArr = s.boundaryItemsArr;
        Iterator<Item> iter = boundaryItemsArr.iterator();
        boolean isAllChedked = true;
        while (iter.hasNext()){
            Item item = iter.next();
            if (item.isChecked == false){
                isAllChedked = false;
                break;
            }
        }
        return isAllChedked;
    }

    public void UpdateItemsInBag(TTPSolution s) {

    }

    public void CocoUpdateIndexFunction(TTPSolution s) {
        int[] x = s.getTour(); //city Array
        int[] z = s.getPickingPlan();
        float[] LPPR = s.getLeastPickedPR();;
        float[] HUPR = s.getHighestUnpickedPR();

        float[] preMinV = s.getPrefixMinimumValue();
        float[] posMaxV = s.getPostfixMaximumValue();

        ArrayList<Item> boundaryItemsArr = s.boundaryItemsArr;
        boundaryItemsArr.clear();
        int[] leastPackedIndex=new int[this.nbCities]; int[] highestUnpackedIndex=new int[this.nbCities];
        for (int i = 0; i < this.nbCities; i++) {
            float leastPickPR_inCity = (float) 100000.0;
            float highestUnPickPR_inCity = (float) 0.0;

            boolean noItemsSelectedInCity = true;
            for (int k = 0; k < (PR_clusters[x[i] - 1]).size(); k++) {
                float PR = PR_clusters[x[i] - 1].get(k);
                int itemIndex = clusters[x[i] - 1].get(k);
                if (z[itemIndex] != 0) {//当前物品被拾取
                    noItemsSelectedInCity = false;
                    if (PR < leastPickPR_inCity) {
                        leastPickPR_inCity = PR;

                        leastPackedIndex[i] = itemIndex;
                    }
                } else {
                    if (PR > highestUnPickPR_inCity) {
                        highestUnPickPR_inCity = PR;
                        highestUnpackedIndex[i]=itemIndex;
                    }
                }
            }
            if (noItemsSelectedInCity) {
                leastPickPR_inCity = highestUnPickPR_inCity + 1;
            }

            LPPR[x[i] - 1] = leastPickPR_inCity;
            HUPR[x[i] - 1] = highestUnPickPR_inCity;
        }

        float startPreMin = LPPR[0];
        for (int k = 0; k < preMinV.length; k++) {
            if (LPPR[k] < startPreMin) {
                startPreMin = LPPR[k];
                Item currentItem = new Item(leastPackedIndex[k],false);
                boundaryItemsArr.add(currentItem);
            }
            preMinV[k] = startPreMin;
        }
        float startPostMax = HUPR[posMaxV.length - 1];
        for (int l = posMaxV.length - 1; l > -1; l--) {
            if (HUPR[l] > startPostMax) {
                startPostMax = HUPR[l];
                Item currentItem = new Item(highestUnpackedIndex[l],false);
                boundaryItemsArr.add(currentItem);
            }
            posMaxV[l] = startPostMax;
        }

        s.setBoundaryItemsArr(boundaryItemsArr);
    }



    //S5部分
    public int getItemIndex(int pPIndex, int[] tour){
        int itemIndex=0;
        int itemsPerCity = this.numberOfItems / (tour.length - 2);

        int itemNum=(pPIndex%itemsPerCity);
        int tourIndex=(int)Math.ceil(pPIndex/itemsPerCity)+1;
        int cityID=tour[tourIndex];
        itemIndex=(itemNum*(tour.length - 2))+(cityID-1);

        //System.out.println("IN: "+itemNum+" TI:"+tourIndex+" CID: "+cityID+" II:"+itemIndex);
        return itemIndex;
    }

    /**
     * Translated code of the original "TTP1Objective.m".
     *
     * Important note: in contrast to the MATLAB code, city numbers start from 0
     * and item numbers start from 0.
     *
//     * @param distances         a n by n matrix that shows the distances between the cities (there are n cities)
//     * @param weights           the weight of each item (1 by m)
//     * @param values            the profit of each item (1 by m)
//     * @param av                a m by n matrix showing if the ith item is available in the jth city.
//     * @param tour              a 1 by n+1 array showing the tour (a complete tour)
//     * @param z                 a 1 by m array, showing which item from which city (if z(i)==j, it means item i from city j)  -->>>> the other way around:
//     * @param weightofKnapsack  maximum weight of the knapsack
//     * @param vmax              maximum velocity
//     * @param vmin              minimum velocity
//     * @param rentRate          the rent rate of the knapsack
     * @return TTP object:
     *          "fp" final profit gained form the picked items,
     *          "ft" the time takes to finish the tour (including changes of the speed),
     *          "ob" objective value,
     *          "wend" weight of the knapsack at the end of the tour
     */
    public void evaluate(TTPSolution solution, boolean printDetails) {
        boolean debugPrint = !false;

        int[] tour = solution.tspTour;
        int[] z = solution.packingPlan;
        long weightofKnapsack = this.capacityOfKnapsack;
        double rentRate = this.rentingRatio;
        double vmin = this.minSpeed;
        double vmax = this.maxSpeed;
        solution.finalDistance = 0;

        // correctness check: does the tour start and end in the same city
        //S5
//        if(tour[0]!=tour[tour.length-1]) {
//            System.out.println("ERROR: The last city must be the same as the first city");
//            solution.reset();
//            return;
//        }

        double weightUntilNow = 0;
        solution.finalTravelTime = 0;
        solution.finalProfit = 0;

        // the following is used for a different interpretation of "packingPlan"
        int itemsPerCity = solution.packingPlan.length / (solution.tspTour.length-2);
        if (debugPrint)
            System.out.println("itemsPerCity="+itemsPerCity+" solution.tspTour.length="+solution.tspTour.length);

        if (printDetails) System.out.println("i,tour[i],tour[h],distance[ih],ftraw,wItemsPickedHere,pItemsPickedHere,wEndFree,wEndUsed,fp,ft,ob");

        for (int i=0; i<tour.length-1; i++) {

            double wItemsPickedHere = 0;
            double pItemsPickedHere = 0;

            // important: nothing to be picked at the first city!
            if (debugPrint)
                System.out.print("\ni="+i+" checking packing: ");

            int currentCityTEMP = tour[i]; // what's the current city? --> but the items start at city 2 in the TTP file, so I have to take another 1 off!

            int currentCity = currentCityTEMP - 1;

            if (debugPrint) {
                System.out.println();
                System.out.print("city "+currentCityTEMP+" cityIndexForItem[][] "+currentCity+" (this.numberOfNodes="+this.numberOfNodes+"): ");
            }
            //查找当前城市中有的物品下标
            if (i>0){
                for (int itemNumber = 0; itemNumber < itemsPerCity; itemNumber++) {
                    int indexOfPackingPlan = (i-1)*itemsPerCity + itemNumber;
                    if (debugPrint) {
                        System.out.print("indexOfPackingPlan="+indexOfPackingPlan+" ");
                    }

                    // what is the next item's index in items-array?
                    int itemIndex = currentCity + itemNumber*(this.numberOfNodes-1);
//                    if(i<3)
                    // System.out.println(indexOfPackingPlan + " CI " + (i) + " id: " + itemIndex);
                    if (debugPrint)
                        System.out.print("itemIndex="+itemIndex+" ");

                    if (z[indexOfPackingPlan]>0) {

                        int currentWC = this.items[itemIndex][2];
                        weightUntilNow = weightUntilNow + currentWC;            // for debugging purposes only
                        wItemsPickedHere = wItemsPickedHere + currentWC;        // for debugging purposes only

                        //System.out.println("i: "+i+" ppI: "+indexOfPackingPlan + " CI " + (currentCity+1) + " id: " + itemIndex+ " IN: "+itemNumber+" WTC: "+currentWC+" WT: "+wc);
                        int currentFP = this.items[itemIndex][1];
                        solution.finalProfit = solution.finalProfit+currentFP;
                        pItemsPickedHere += currentFP;                          // for debugging purposes only

                        if (debugPrint)
                            System.out.print("[fp="+currentFP+",wc="+currentWC+"] ");
                    }
                }
            }
            if (debugPrint)
                System.out.println();

            int h = (i+1)%(tour.length-1); //h: next tour city index
            if (debugPrint)
                System.out.println("  i="+i+" h="+h + " tour[i]="+tour[i]+" tour[h]="+tour[h]);

            long distance = (long)Math.ceil(distances(tour[i],tour[h]));

            // compute the raw distance
            solution.finalDistance += distance;

            // compute the adjusted (effective) distance
            solution.finalTravelTime = solution.finalTravelTime + (distance / (vmax-weightUntilNow*(vmax-vmin)/weightofKnapsack)); // "1-" replaced by "vmax-" on 11/09/2016 based on comment by Christopher Crouch

            if (debugPrint)
                System.out.println("i="+i+" tour[i]="+tour[i]+" tour[h]="+tour[h]+" distance="+distance+" fp="+solution.finalProfit + " ft=" + solution.finalTravelTime);

            if (printDetails) {
                System.out.println(i+","+tour[i]+","+tour[h]+","+distance+","+
                        solution.finalDistance+","+
                        wItemsPickedHere+","+pItemsPickedHere+","+
                        (weightofKnapsack - weightUntilNow)+","+weightUntilNow+","+
                        solution.finalProfit+","+
                        solution.finalTravelTime+","+(solution.finalProfit - solution.finalTravelTime * rentRate) );
            }
        }
        //System.out.println("wc:" + wc);
        solution.finalWeight = weightUntilNow;
        solution.finalCapacityFree = weightofKnapsack - weightUntilNow;
        solution.objectiveScore = solution.finalProfit - solution.finalTravelTime * rentRate;
    }





    // used to simulate the distance matrix
    public double distances(int i, int j) {
        boolean debugPrint = false;
        double result = 0;
        result = Math.sqrt(
                (this.nodes[i][1]-this.nodes[j][1]) * (this.nodes[i][1]-this.nodes[j][1])
                        + (this.nodes[i][2]-this.nodes[j][2]) * (this.nodes[i][2]-this.nodes[j][2])
        );

        if (debugPrint)
            System.out.println(" distance="+this.nodes[i][1]+ " " +this.nodes[j][1]+" "+this.nodes[i][2]+" "+this.nodes[j][2]+"->"+result);

        return result;
    }

    public void printInstance() {
        printInstance(true);
    }

    /*
     * prints the details
     * shortSummary:
     *   true   prints a short version in one line
     *   false  prints the node and item data as well
     */
    public void printInstance(boolean shortSummary) {
        if (shortSummary) {
            System.out.print("TTP Instance: ");
        } else {
            System.out.println("---- TTP Instance START ----");
        }

        System.out.println(
                this.problemName +
                        " " + this.knapsackDataType +
                        " " + this.numberOfNodes +
                        " " + this.numberOfItems +
                        " " + this.capacityOfKnapsack +
                        " " + this.minSpeed +
                        " " + this.maxSpeed +
                        " " + this.rentingRatio
        );

        if (shortSummary) {
        } else {
            for (double[] i:this.nodes) {
                System.out.println(Arrays.toString(i));
            }
            for (int[] i:this.items) {
                System.out.println(Arrays.toString(i));
            }
            System.out.println("---- TTP Instance END ----");
        }
    }

    public void adjustS5Evaluate(TTPSolution s, boolean printDetails) {

        int[] x = s.getTour(); //city Array
        int[] z = s.getPickingPlan();

        int[] CLKresult = new int[s.tspTour.length-1];
        for (int i = 0; i < x.length-1; i++) {
            CLKresult[i] = x[i]+1;
        }

        long[][] D = getDist();
        double C = (maxSpeed - minSpeed) / capacity; // velocity const
        double velocity;

        long acc;       // iteration weight accumulator 迭代权累加器
        long wc = 0;    // current weight
        long fp = 0;    // final profit
        double ft = 0;  // tour time
        double ob;      // objective value


        // visit all cities
        for (int i = 0; i < this.nbCities; i++) {
            acc = 0;
            // check only items contained in current city
            for (int j : clusters[CLKresult[i] - 1]) { //根据tour处在当前位置城市中的物品下标:j
                if (z[j] != 0) {
                    fp += profits[j];
                    acc += weights[j];
                }
            }

            wc += acc;
            velocity = maxSpeed - wc * C;

            int h = (i + 1) % nbCities; //next城市
            ft += distFor(CLKresult[i] - 1, CLKresult[h] - 1) / velocity; //去往下个城市需要花费的时间

            // record important data for future use
            s.timeAcc[i] = ft;//until now time
            s.timeRec[i] = distFor(CLKresult[i] - 1, CLKresult[h] - 1) / velocity;//当前城市到下一个城市花费的时间
            s.weightAcc[i] = wc;//到当前城市时累积的重量
            s.weightRec[i] = acc;//在当前城市拾取累积的重量

            // map indices to their associated cities 将城市和索引对调
            s.mapCI[CLKresult[i] - 1] = i;//当前城市xCLKresult[i]在tour中的下标为i
        }

        ob = fp - ft * rent;

        // solution properties
        s.fp = fp;
        s.ft = ft;
        s.wend = capacity - wc;
        s.ob = ob;

        s.objectiveScore = ob;
    }
}
