package solver;

import ttp.TTP1Instance;
import ttp.TTPSolution;
import utils.Deb;
import utils.GraphHelper;
import utils.Quicksort;
import utils.TwoOptHelper;

import java.text.DecimalFormat;
import java.util.HashSet;

public abstract class LocalSearch extends SearchHeuristic  {

    protected TTPSolution s0;
    protected HashSet<Integer>[] candidates;

    // initial solution
    public void setS0(TTPSolution s0) {
        this.s0 = s0;
    }
    public TTPSolution getS0() {
        return s0;
    }

    /**
     * first fit or best fit
     */
    public boolean firstfit;


    public LocalSearch() {
        super();
    }

    public LocalSearch(TTP1Instance ttp) {
        super(ttp);
        // generate Delaunay triangulation
        if (! System.getProperty("os.name").contains("Windows")) {    //返回了"Linux",so result is false
//      candidates = GraphHelper.delaunayKNN(ttp,10);
            candidates = GraphHelper.delaunay(ttp);

            if (debug) Deb.echo("Delaunay triangulation: OK");
        }
    }



    @Override
    public void setTTP(TTP1Instance ttp) {
        super.setTTP(ttp);
        // generate Delaunay triangulation
        if (! System.getProperty("os.name").contains("Windows")) {
//    candidates = GraphHelper.delaunayKNN(ttp,5);
            candidates = GraphHelper.delaunay(ttp);
        }
    }

    public void setDelaunayLevel(int L) {
        // generate Delaunay triangulation
        if (! System.getProperty("os.name").contains("Windows")) {
            candidates = GraphHelper.delaunayKNN(ttp, L);
        }
    }


    /**
     * use first fit strategy
     */
    public void firstfit() {
        firstfit = true;
    }


    /**
     * use best fit strategy
     */
    public void bestfit() {
        firstfit = false;
    }

    @Override
    public String getName() {
        String suf = this.firstfit ? "-FF" : "-BF";
        return this.name + suf;
    }




    /**
     * KP pre-processing
     *
     * base on the item insertion heuristic in MATLS
     * described in "improving efficiency of heuristics..."
     *
     * uses only delta_1 and delta_2 approximations
     * in order to limit the number of inserted items
     */
    public TTPSolution insertT2(TTPSolution sol) {

        // TTP data
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        long[][] D = ttp.getDist();
        int[] A = ttp.getAvailability();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;
        double R = ttp.getRent();

        // initial solution data
        int[] tour = sol.getTour();
        int[] pickingPlan = sol.getPickingPlan();

        // neighbor solution
        int origBF;
        int i, k, itr;

        // distances of all tour cities (city -> end)
        long[] L = new long[nbCities];
        // current weight
        long wCurr;
        // time approximations
        double t1, t2, t3, a, b1, b2;

        // store `distance to end` of each tour city
        L[nbCities - 1] = ttp.distFor(tour[nbCities - 1] - 1, 0);
        for (i = nbCities - 2; i >= 0; i--) {
            L[i] = L[i + 1] + ttp.distFor(tour[i + 1] - 1, tour[i] - 1);
        }

        // sort item according to score
        Double[] scores = new Double[nbItems];
        int[] insertedItems = new int[nbItems];

        for (k = 0; k < nbItems; k++) {
            // index where Bit-Flip happened
            origBF = sol.mapCI[A[k] - 1];
            // calculate time approximations
            t1 = L[origBF] * (1 / (maxSpeed - C * ttp.weightOf(k)) - 1 / maxSpeed);
            // affect score to item
            scores[k] = (ttp.profitOf(k) - R * t1) / ttp.weightOf(k);
            // empty the knapsack
            pickingPlan[k] = 0;
        }

        // evaluate solution after emptying knapsack
        ttp.objective(sol);

        // sort items according to score
        Quicksort<Double> qs = new Quicksort<>(scores);
        qs.sort();
        int[] sortedItems = qs.getIndices();

        // loop & insert items
        int nbInserts = 0;
        wCurr = 0;
        int v2 = 0, v3 = 0;
        for (itr = 0; itr < nbItems; itr++) {

            k = sortedItems[itr];

            // check if new weight doesn't exceed knapsack capacity
            if (wCurr + ttp.weightOf(k) > capacity) {
                continue;
            }

            // index where Bit-Flip happened
            origBF = sol.mapCI[A[k] - 1];

            /* insert item if it has a potential gain */
            // time approximations t2 (worst-case time)
            t2 = L[origBF] * (1 / (maxSpeed - C * (wCurr + ttp.weightOf(k))) - 1 / (maxSpeed - C * wCurr));
            if (ttp.profitOf(k) > R * t2) {
                v2++;
                pickingPlan[k] = A[k];
                wCurr += ttp.weightOf(k);
                insertedItems[nbInserts++] = k;
            }
        } // END FOR k

        // evaluate solution & update vectors
        ttp.objective(sol);

        // debug msg
        if (this.debug) {
            Deb.echo(">> item insertion: best=" + sol.ob);
            Deb.echo("   wend: " + sol.wend);

            Deb.echo("==> nb t2: " + v2 + " | nb t3: " + v3);
            Deb.echo("==> nb inserted: " + nbInserts + "/" + nbItems + "(" +
                    String.format("%.2f", (nbInserts * 100.0) / nbItems) + "%)");
            Deb.echo("==> w_curr: " + wCurr);
        }

        return sol;
    }


    /**
     * KP pre-processing
     *
     * base on the item insertion heuristic in MATLS
     * described in "improving efficiency of heuristics..."
     */
    public TTPSolution insertAndEliminate(TTPSolution sol) {

        // TTP data
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        long[][] D = ttp.getDist();
        int[] A = ttp.getAvailability();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;
        double R = ttp.getRent();

        // initial solution data
        int[] tour = sol.getTour();
        int[] pickingPlan = sol.getPickingPlan();

        // neighbor solution
        int origBF;
        int i, k, itr;

        // distances of all tour cities (city -> end)
        long[] L = new long[nbCities];
        // current weight
        long wCurr;
        // time approximations
        double t1, t2, t3, a, b1, b2;

        // store `distance to end` of each tour city
        L[nbCities-1] = ttp.distFor(tour[nbCities-1] - 1,0);
        for (i=nbCities-2; i >= 0; i--) {
            L[i] = L[i+1] + ttp.distFor(tour[i+1]-1, tour[i]-1);
        }
        sol.L = L;
        // sort item according to score
        Double[] scores = new Double[nbItems];
        int[] insertedItems = new int[nbItems];

        for (k = 0; k < nbItems; k++) {
            // index where Bit-Flip happened
            origBF = sol.mapCI[A[k] - 1];  //计算当前item的城市index
            // calculate time approximations
            t1 = L[origBF]*(1/(maxSpeed-C*ttp.weightOf(k)) - 1/maxSpeed);//计算当前物品选取与否的时间差
            // affect score to item
//            scores[k] = (ttp.profitOf(k)-R*t1) / ttp.weightOf(k);//计算性价比
            scores[k] = (ttp.profitOf(k)-R*t1) / (Math.pow(ttp.weightOf(k),1.5));//计算性价比
            // empty the knapsack
            pickingPlan[k] = 0;
        }

        // evaluate solution after emptying knapsack
        ttp.objective(sol);

        // sort items according to score
        Quicksort<Double> qs = new Quicksort<>(scores);
        qs.sort();
        int[] sortedItems = qs.getIndices();
        Double[] sortedData = qs.getData();
        Double standData = 0.0;
        Double sumData=0.0;int count=0;
        for (int g = 0; g < sortedData.length; g++) {
            if (sortedData[g]>0){
                sumData+=sortedData[g];
                count++;
            }
        }
        standData = sumData/count;
//        standData = sortedData[0]/10;
//        standData = (standData+sortedData[0])/2;

        count=0;
        for (int g = 0; g < sortedData.length; g++) {
            if (sortedData[g]>standData){
                sumData+=sortedData[g];
                count++;
            }
        }

        DecimalFormat df = new DecimalFormat("0.0000") ;
        Double rato = Double.parseDouble(df.format((float)count/(1*nbItems))) ;
        // loop & insert items
        int nbInserts = 0;
        wCurr = 0;
        int v2 = 0, v3 = 0;


        //----------------------我的代码-----
        standData = standData+(sortedData[0]-standData)*0.65;
//        standData = standData*0.5;
        for (int s = nbCities - 1; s > 0; s--) {
            int cityNode = tour[s] - 1;
            for (int j : ttp.clusters[tour[s] - 1]) {
                if (wCurr + ttp.weightOf(j) > capacity*rato) {
                    continue;
                }
                if (scores[j] > standData) {
                    pickingPlan[j] = A[j];
                    wCurr += ttp.weightOf(j);
                    insertedItems[nbInserts++] = j;
                }else break;
            }
        }


        //原始代码

        for (itr = 0; itr < nbItems; itr++) {

            k = sortedItems[itr];

            // check if new weight doesn't exceed knapsack capacity
            if (wCurr + ttp.weightOf(k) > capacity) {
                continue;
            }

            // index where Bit-Flip happened
            origBF = sol.mapCI[A[k] - 1];//计算当前item的城市index
            if (pickingPlan[k] != 0) continue;
            /* insert item if it has a potential gain */
            // time approximations t2 (worst-case time)
            t2 = L[origBF] * (1/(maxSpeed-C*(wCurr+ttp.weightOf(k))) - 1/(maxSpeed-C*wCurr));
            if (ttp.profitOf(k) > R*t2) {
                v2++;
                pickingPlan[k] = A[k];//第k个物品所在的城市
                wCurr += ttp.weightOf(k);
                insertedItems[nbInserts++] = k;//nbInserts代表插入物品的数量，insertedItems代表选取物品数组
            }
            else {
                // time approximations t3 (expected time)
                a = wCurr / L[0];
                b1 = maxSpeed - C * (wCurr + ttp.weightOf(k));//拾取k物品后的速度
                b2 = maxSpeed - C * wCurr;//不拾取k物品的速度
                t3 = (1 / a) * Math.log(
                        ( (a * L[0] + b1) * (a * (L[0] - L[origBF]) + b2) ) /
                                ( (a * (L[0] - L[origBF]) + b1) * (a * L[0] + b2) )
                );
                if (ttp.profitOf(k) > R*t3) {v3++;
                    pickingPlan[k] = A[k];
                    wCurr += ttp.weightOf(k);
                    insertedItems[nbInserts++] = k;
                }
                else continue;
            }
        } // END FOR k

        // evaluate solution & update vectors
        ttp.objective(sol);

        // debug msg
        if (this.debug) {
            Deb.echo(">> item insertion: best=" + sol.ob);
            Deb.echo("   wend: "+sol.wend);

            Deb.echo("==> nb t2: "+v2+" | nb t3: "+v3);
            Deb.echo("==> nb inserted: "+nbInserts+"/"+nbItems+"("+
                    String.format("%.2f", (nbInserts * 100.0) / nbItems)+"%)");
            Deb.echo("==> w_curr: "+wCurr);
        }

        if (nbItems > 100000 || nbCities > 50000) return sol;



        //===========================
        // elimination
        //===========================
        // best solution
        int kBest = 0;
        double GBest = sol.ob;

        // neighbor solution
        long fp;
        double ft, G;
        long wc;
        int r;



        // improvement indicator
        boolean improved = false;

        // browse items in the new order...
        for (itr = 0; itr < nbInserts; itr++) {
            k = insertedItems[itr];//insertedItems代表选取物品数组

            fp = sol.fp - ttp.profitOf(k);//去掉该物品后的价值

            // index where Bit-Flip happened
            origBF = sol.mapCI[A[k] - 1];//表示已经走过的城市个数，nbCities-OrigBF表示剩下没走的城市数mCities-OrigBF表示剩下没走的城市数

            // starting time
            ft = origBF == 0 ? 0 : sol.timeAcc[origBF - 1];

            // recalculate velocities from bit-flip city
            for (r = origBF; r < nbCities; r++) {
                wc = sol.weightAcc[r] - ttp.weightOf(k);;
                ft += ttp.distFor(tour[r]-1, tour[(r + 1) % nbCities]-1) / (maxSpeed - wc * C);
            }//去掉该物品后所花费的时间

            G = Math.round(fp - ft * R);

            // update best
            if (G > GBest) {

                kBest = k;
                GBest = G;
                improved = true;
                if (firstfit) break;
            }

        } // END FOR k

        /* update if improvement */
        if (improved) {

            // bit-flip
            pickingPlan[kBest] = 0;

            // evaluate & update vectors
            ttp.objective(sol);

            // debug msg
            if (this.debug) {
                Deb.echo(">> item elimination: best=" + sol.ob);
            }
        }

        return sol;
    }

    public TTPSolution DeleteAndInsert(TTPSolution sol){
        TTPSolution sBest = sol.clone();
        int[] tour = sol.getTour();
        int[] pickingPlan = sol.getPickingPlan();
        int nbItems = ttp.getNbItems();
        int[] A = ttp.getAvailability();

        int k,l;
//        CalculateScoreArr(sol);
        Double [] scores = sol.scoresArr;
        Quicksort<Double> qs = new Quicksort<>(scores);
        qs.sort();
        int[] sortedItems = qs.getIndices();
        int DeleteNum = (int) Math.floor(Math.sqrt(nbItems));
        int InsertNum = (int) Math.floor(Math.sqrt(nbItems));
        while (DeleteNum>0){
            for (int i = scores.length-1; i > 0; i--) {
                if (DeleteNum<=0) break;
                k = sortedItems[i];
                if (pickingPlan[k]!=0){
                    pickingPlan[k]=0;
                    sol.wend+=ttp.weightOf(k);

                }
                DeleteNum--;
            }
        }
//        CalculateScoreArr(sol);
        while (InsertNum>0){
            for (int i = scores.length-1; i > 0; i--) {
                if (InsertNum<=0 ) break;
                k = sortedItems[i];
                if (pickingPlan[k] == 0 && ttp.weightOf(k) > sol.wend) continue;
                if (pickingPlan[k]==0){
                    pickingPlan[k] =  A[k];
                    sol.wend-=ttp.weightOf(k);

                }
                InsertNum--;
            }
        }
        ttp.objective(sol);
        if (sol.ob > sBest.ob) {
            sBest = sol.clone();
        }
        return sBest;
    }
    public  void CalculateScoreArr(TTPSolution s){
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        int[] A = ttp.getAvailability();
        int[] tour = s.getTour();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity; // velocity const
        double R = ttp.getRent();
        // neighbor solution
        int origBF,k,i;

        // distances of all tour cities (city -> end)
        long[] L = s.L;

        // time approximations
        double t2;

        for (k = 0; k < nbItems; k++) {
            // index where Bit-Flip happened
            origBF = s.mapCI[A[k] - 1];  //计算当前item的城市index
            // calculate time approximations
            t2 = L[origBF]*(1/(maxSpeed-C*(ttp.weightOf(k)+s.weightAcc[origBF])) - 1/maxSpeed);//计算当前物品选取与否的时间差
            // affect score to item
            s.scoresArr[k] = ttp.profitOf(k)-R*t2;
        }

        //计算均值
        Quicksort<Double> qs = new Quicksort<>(s.scoresArr);
        qs.sort();
        int[] sortedItems = qs.getIndices();
        Double[] sortedData = qs.getData();
        Double standData = 0.0,sumWeight = 0.0,averageWeight=0.0;
        Double sumData=0.0;int count=0;
        for (int g = 0; g < sortedData.length; g++) {
            if (sortedData[g]>0){
                sumData+=sortedData[g];
                sumWeight+=ttp.weightOf(sortedItems[g]);
                count++;
            }
        }
        standData = sumData/count;
        averageWeight = sumWeight/count;
        s.averageWeight = averageWeight;
        s.standScore = standData;
    }
    public TTPSolution HybridInsert(TTPSolution sol) {
        // TTP data
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        long[][] D = ttp.getDist();
        int[] A = ttp.getAvailability();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;
        double R = ttp.getRent();

        // initial solution data
        int[] tour = sol.getTour();
        int[] pickingPlan = sol.getPickingPlan();

        // neighbor solution
        int origBF;
        int i, k, itr;

        // distances of all tour cities (city -> end)
        long[] L = new long[nbCities];
        // current weight
        long wCurr;
        // time approximations
        double t1, t2, t3, a, b1, b2;

        // store `distance to end` of each tour city
        L[nbCities-1] = ttp.distFor(tour[nbCities-1] - 1,0);
        for (i=nbCities-2; i >= 0; i--) {
            L[i] = L[i+1] + ttp.distFor(tour[i+1]-1, tour[i]-1);
        }

        // sort item according to score
        Double[] scores = new Double[nbItems];
        int[] insertedItems = new int[nbItems];

        for (k = 0; k < nbItems; k++) {
            // index where Bit-Flip happened
            origBF = sol.mapCI[A[k] - 1];  //计算当前item的城市index
            // calculate time approximations
            t1 = L[origBF]*(1/(maxSpeed-C*ttp.weightOf(k)) - 1/maxSpeed);//计算当前物品选取与否的时间差
            // affect score to item
//            scores[k] = (ttp.profitOf(k)-R*t1) / ttp.weightOf(k);//计算性价比
            scores[k] = (Math.pow(ttp.profitOf(k),7.4)*(L[0]-L[origBF])/maxSpeed)/(Math.pow(ttp.weightOf(k),7.4)*t1);
            // empty the knapsack
            pickingPlan[k] = 0;

//            ttp.PR_clusters[ttp.availability[k]-1].add(scores[k].floatValue());
        }

        // evaluate solution after emptying knapsack
        ttp.objective(sol);


        // sort items according to score
        Quicksort<Double> qs = new Quicksort<>(scores);
        qs.sort();
        int[] sortedItems = qs.getIndices();

        // loop & insert items
        int nbInserts = 0;
        wCurr = 0;
        int v2 = 0, v3 = 0;

        double GStar = Double.NEGATIVE_INFINITY;
        for (itr = 0; itr < nbItems; itr++) {

            k = sortedItems[itr];

            // check if new weight doesn't exceed knapsack capacity
            if (wCurr + ttp.weightOf(k) > capacity) {
                continue;
            }
            // index where Bit-Flip happened
            origBF = sol.mapCI[A[k] - 1];//计算当前item的城市index

            /* insert item if it has a potential gain */
            // time approximations t2 (worst-case time)
            t2 = L[origBF] * (1/(maxSpeed-C*(wCurr+ttp.weightOf(k))) - 1/(maxSpeed-C*wCurr));
            if (ttp.profitOf(k) > R*t2) {
                v2++;
                pickingPlan[k] = A[k];//第k个物品所在的城市
                wCurr += ttp.weightOf(k);
                insertedItems[nbInserts++] = k;//nbInserts代表插入物品的数量，insertedItems代表选取物品数组

                ttp.objective(sol);
                GStar = sol.ob;
            }else {
                wCurr += ttp.weightOf(k);
                pickingPlan[k] = A[k];
                ttp.objective(sol);
                if (sol.ob<GStar){
                    pickingPlan[k] = 0;
                    wCurr -= ttp.weightOf(k);
                }else {
                    GStar = sol.ob;
                }
            }

        }

        return sol;
    }



    public TTPSolution Eliminate(TTPSolution sol) {

        // TTP data
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        long[][] D = ttp.getDist();
        int[] A = ttp.getAvailability();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;
        double R = ttp.getRent();

        // initial solution data
        int[] tour = sol.getTour();
        int[] pickingPlan = sol.getPickingPlan();

        // neighbor solution
        int origBF;
        int i, k, itr;






        //===========================
        // elimination
        //===========================
        // best solution
        int kBest = 0;
        double GBest = sol.ob;

        // neighbor solution
        long fp;
        double ft, G;
        long wc;
        int r;
        // delta parameters
        int deltaP, deltaW;

        // best solution


        // improvement indicator
        boolean improved = false;

        // browse items in the new order...
        for (k = 0; k < nbItems; k++) {

            if (pickingPlan[k] == 0 && ttp.weightOf(k) > sol.wend) continue;

            // calculate deltaP and deltaW
            if (pickingPlan[k] != 0) {
                deltaP = -ttp.profitOf(k);
                deltaW = -ttp.weightOf(k);
                fp = sol.fp + deltaP;

                // index where Bit-Flip happened
                origBF = sol.mapCI[A[k] - 1];//表示已经走过的城市个数，nbCities-OrigBF表示剩下没走的城市数mCities-OrigBF表示剩下没走的城市数

                // starting time
                ft = origBF == 0 ? 0 : sol.timeAcc[origBF - 1];

                // recalculate velocities from bit-flip city
                for (r = origBF; r < nbCities; r++) {
                    wc = sol.weightAcc[r] - ttp.weightOf(k);;
                    ft += ttp.distFor(tour[r]-1, tour[(r + 1) % nbCities]-1) / (maxSpeed - wc * C);
                }//去掉该物品后所花费的时间

                G = (fp - ft * R);
                // update best
                if (G > GBest) {

                    kBest = k;
                    GBest = G;
                    improved = true;
                }
            }

        } // END FOR k

        /* update if improvement */
        if (improved) {

            // bit-flip
            pickingPlan[kBest] = 0;

            // evaluate & update vectors
            ttp.objective(sol);

            // debug msg
            if (this.debug) {
                Deb.echo(">> item elimination: best=" + sol.ob);
            }
        }

        return sol;
    }


    /**
     * 2-opt search
     *
     * deal with the TSKP sub-problem
     * 2-opt heuristic with Delaunay candidate generator
     */
    public TTPSolution fast2opt(TTPSolution sol) {

        // TTP data
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;

        // initial solution data
        int[] tour;

        // delta parameters
        double deltaT;

        // improvement indicator
        boolean improved;

        // best solution
        ttp.objective(sol);
        int iBest=0, jBest=0;
        double ftBest = sol.ft;

        // neighbor solution
        double ft;
        long wc;
        int i, j, c1, c2, q;
        int nbIter = 0;

        // current tour
        tour = sol.getTour();

        // search params
        double threshold = -0.1;
        if (nbItems >= 100000) {
            threshold = -10;
        }
        if (nbCities >= 50000) { // ex. pla85000 based instances
            threshold = -1000;
        }

        // search
        do {
            improved = false;
            nbIter++;

            // cleanup and stop execution if interrupted
            if (Thread.currentThread().isInterrupted()) break;

            // fast 2-opt
            for (i = 1; i < nbCities - 1; i++) {
                int node1 = tour[i] - 1;
                for (int node2 : candidates[node1]) {
                    j = sol.mapCI[node2];
                    //if (j<=i) continue;

                    // calculate final time with partial delta
                    ft = sol.ft;
                    wc = i - 2 < 0 ? 0 : sol.weightAcc[i - 2]; // fix index...
                    deltaT = 0;
                    for (q = i - 1; q <= j; q++) {//两点之间有几个城市

                        wc += TwoOptHelper.get2optValue(q, sol.weightRec, i, j);
                        c1 = TwoOptHelper.get2optValue(q, tour, i, j) - 1;
                        c2 = TwoOptHelper.get2optValue((q + 1) % nbCities, tour, i, j) - 1;

                        deltaT += -sol.timeRec[q] + ttp.distFor(c1,c2) / (maxSpeed - wc * C);
                    }

                    // retrieve neighbor's final time
                    ft = ft + deltaT;

                    // update best
                    if (ft - ftBest < threshold) { // soft condition
                        iBest = i;
                        jBest = j;
                        ftBest = ft;
                        improved = true;

                        if (firstfit) break;
                    }

                    //if (firstfit && improved) break;
                } // END FOR j
                if (firstfit && improved) break;
            } // END FOR i


            //===================================
            // update if improvement
            //===================================
            if (improved) {

                // apply 2-opt move
                TwoOptHelper.do2opt(tour, iBest, jBest);

                // evaluate & update vectors
                ttp.objective(sol);
            }

            // debug msg
            if (this.debug) {
                Deb.echo(">> TSKP " + nbIter +
                        ": ob=" + String.format("%.0f", sol.ob) +
                        " | ft=" + String.format("%.0f", sol.ft));
            }

        } while (improved && nbIter<maxIterTSKP);

        if (debug) Deb.echo("==> 2-opt :" + nbIter + " iterations");

        // in order to compute sol.timeAcc
        // we need to use objective function
        ttp.objective(sol);

        return sol;
    }


    /**
     * 2-opt search
     *
     * deal with the TSKP sub-problem
     * 2-opt heuristic without neighborhood reduction techniques
     */
    public TTPSolution slow2opt(TTPSolution sol) {

        // TTP data
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;

        // initial solution data
        int[] tour;

        // delta parameters
        double deltaT;

        // improvement indicator
        boolean improved;

        // best solution
        ttp.objective(sol);
        int iBest=0, jBest=0;
        double ftBest = sol.ft;

        // neighbor solution
        double ft;
        long wc;
        int i, j, c1, c2, q;
        int nbIter = 0;

        // Delaunay triangulation
        //ArrayList<Integer>[] candidates = ttp.delaunay();

        // current tour
        tour = sol.getTour();

        // search params
        double threshold = -0.1;
        if (nbItems >= 100000) {
            threshold = -10;
        }
        if (nbCities >= 50000) { // ex. pla85000 based instances
            threshold = -1000;
        }

        // search
        do {
            improved = false;
            nbIter++;

            // cleanup and stop execution if interrupted
            if (Thread.currentThread().isInterrupted()) break;

            // fast 2-opt
            for (i = 1; i < nbCities - 1; i++) {
                for (j=i+1; j < nbCities; j++) {

                    // calculate final time with partial delta
                    ft = sol.ft;
                    wc = i - 2 < 0 ? 0 : sol.weightAcc[i - 2]; // fix index...
                    deltaT = 0;
                    for (q = i - 1; q <= j; q++) {

                        wc += TwoOptHelper.get2optValue(q, sol.weightRec, i, j);
                        c1 = TwoOptHelper.get2optValue(q, tour, i, j) - 1;
                        c2 = TwoOptHelper.get2optValue((q + 1) % nbCities, tour, i, j) - 1;

                        deltaT += -sol.timeRec[q] + ttp.distFor(c1,c2) / (maxSpeed - wc * C);
                    }

                    // retrieve neighbor's final time
                    ft = ft + deltaT;

                    // update best
                    if (ft - ftBest < threshold) { // soft condition
                        iBest = i;
                        jBest = j;
                        ftBest = ft;
                        improved = true;

                        if (firstfit) break;
                    }

                    //if (firstfit && improved) break;
                } // END FOR j
                if (firstfit && improved) break;
            } // END FOR i


            //===================================
            // update if improvement
            //===================================
            if (improved) {

                // apply 2-opt move
                TwoOptHelper.do2opt(tour, iBest, jBest);

                // evaluate & update vectors
                ttp.objective(sol);
            }

            // debug msg
            if (this.debug) {
                Deb.echo(">> TSKP " + nbIter +
                        ": ob=" + String.format("%.0f", sol.ob) +
                        " | ft=" + String.format("%.0f", sol.ft));
            }

        } while (improved && nbIter<maxIterTSKP);


        // in order to compute sol.timeAcc
        // we need to use objective function
        ttp.objective(sol);
        return sol;
    }


    /**
     * bit-flip search
     *
     * deal with the KRP sub-problem
     * this function applies a simple bit-flip
     */
    public TTPSolution lsBitFlip(TTPSolution sol) {

        // TTP data
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        int[] A = ttp.getAvailability();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;
        double R = ttp.getRent();

        // initial solution data
        int[] tour = sol.getTour();
        int[] pickingPlan = sol.getPickingPlan();

        // delta parameters
        int deltaP, deltaW;

        // best solution
        double GBest = sol.ob;

        // neighbor solution
        long fp;
        double ft, G;
        long wc;
        int origBF;
        int k, r, kBest=0;
        int nbIter = 0;

        boolean improved;

        // start search
        do {
            improved = false;
            nbIter++;

            // browse items in the new order...
            for (k = 0; k < nbItems; k++) {

                // cleanup and stop execution if interrupted
                if (Thread.currentThread().isInterrupted()) break;

                // check if new weight doesn't exceed knapsack capacity
                if (pickingPlan[k] == 0 && ttp.weightOf(k) > sol.wend) continue;

                // calculate deltaP and deltaW
                if (pickingPlan[k] == 0) {
                    deltaP = ttp.profitOf(k);
                    deltaW = ttp.weightOf(k);
                } else {
                    deltaP = -ttp.profitOf(k);
                    deltaW = -ttp.weightOf(k);
                }
                fp = sol.fp + deltaP;


                // index where Bit-Flip happened
                origBF = sol.mapCI[A[k] - 1];

                // starting time
                ft = origBF == 0 ? 0 : sol.timeAcc[origBF - 1];

                // recalculate velocities from bit-flip city
                for (r = origBF; r < nbCities; r++) {
                    wc = sol.weightAcc[r] + deltaW;
                    ft += ttp.distFor(tour[r]-1,tour[(r + 1) % nbCities] - 1) / (maxSpeed - wc * C);
                }

                G = fp - ft * R;

                // update best
                if (G > GBest) {
                    kBest = k;
                    GBest = G;

                    improved = true;
                    if (firstfit) break;
                }

            } // END FOR k


            //=====================================
            // update if improvement
            //=====================================
            if (improved) {

                // bit-flip
                pickingPlan[kBest] = pickingPlan[kBest] != 0 ? 0 : A[kBest];


                //===========================================================
                // recover accumulation vectors
                //===========================================================
                if (pickingPlan[kBest] != 0) {
                    deltaP = ttp.profitOf(kBest);
                    deltaW = ttp.weightOf(kBest);
                } else {
                    deltaP = -ttp.profitOf(kBest);
                    deltaW = -ttp.weightOf(kBest);
                }
                fp = sol.fp + deltaP;
                origBF = sol.mapCI[A[kBest] - 1];
                ft = origBF == 0 ? 0 : sol.timeAcc[origBF - 1];
                for (r = origBF; r < nbCities; r++) {
                    // recalculate velocities from bit-flip city
                    wc = sol.weightAcc[r] + deltaW;
                    ft += ttp.distFor(tour[r] - 1, tour[(r + 1) % nbCities] - 1) / (maxSpeed - wc * C);
                    // recover wacc and tacc
                    sol.weightAcc[r] = wc;
                    sol.timeAcc[r] = ft;
                }
                G = fp - ft * R;
                sol.ob = G;
                sol.fp = fp;
                sol.ft = ft;
                sol.wend = capacity - sol.weightAcc[nbCities - 1];
                //===========================================================

                // debug msg
                if (!this.debug) {
                    Deb.echo(">> KRP: " + nbIter +
                            " | ob=" + String.format("%.2f",sol.ob) +
                            " | ft=" + String.format("%.2f",sol.ft)
                    );
                }
            }

        } while (improved && nbIter<maxIterKRP);

        if (debug) Deb.echo("==> bitflip :" + nbIter + " iterations");
        if (this.debug) {
            Deb.echo(">> KRP: " + nbIter +
                    " | ob=" + String.format("%.2f",sol.ob) +
                    " | ft=" + String.format("%.2f",sol.ft)
            );
        }
        // in order to recover all history vectors
        ttp.objective(sol);

        return sol;
    }
}
