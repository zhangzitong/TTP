package solver;

import Optimisation.Optimisation;
import ea.Initialization;
import newrep.Item;
import ttp.TTP1Instance;
import ttp.TTPSolution;
import utils.Deb;
import utils.RandGen;
import utils.TwoOptHelper;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class CocoSolver extends LocalSearch {

    public double T_abs;       // absolute temperature
    public double T0;           // initial temperature
    public double alpha;       // cooling rate
    public double trialFactor; // number of trials (per temperature)

    public int[] tourS5;
    public int[] CLKresult;

    public Constructive construct;
    public Initialization init;


    public double GTSPBest;
    public boolean tspImprove;

    public CocoSolver() {
        super();
        // use default config
        SAConfig();
    }

    public CocoSolver(TTP1Instance ttp) {
        super(ttp);
        // use default config
        SAConfig();
    }


    // SA params config
    // default config
    void SAConfig() {

        int nbItems = ttp.getNbItems();

        T_abs = 1;
//    T0 = 100.0;
//    alpha = 0.95;
        alpha = 0.9578;
        T0 = 98;

        trialFactor = generateTFLinFit(nbItems);
    }



    /**
     * simulated annealing
     *
     * deal with the KRP sub-problem
     * this function applies a simple bit-flip
     */
    public TTPSolution simulatedAnnealing(TTPSolution sol) {

        // copy initial solution into improved solution
        TTPSolution sBest = sol.clone();

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
        int k, r;
        int nbIter = 0;

        double T = T0;

        long trials = Math.round(nbItems*trialFactor);

        if (debug) Deb.echo(">>>> TRIAL FACTOR: "+trialFactor);

        //===============================================
        // start simulated annealing process
        //===============================================
        long[] L = sol.L;
        do {
            nbIter++;

            // cleanup and stop execution if interrupted
            if (Thread.currentThread().isInterrupted()) break;

            for (int u=0; u<trials; u++) {

                // browse items randomly
                k = RandGen.randInt(0, nbItems - 1);

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

                // handle velocity constraint
                // index where Bit-Flip happened
                origBF = sol.mapCI[A[k] - 1];
                // starting time
                ft = origBF == 0 ? .0 : sol.timeAcc[origBF - 1];
                // recalculate velocities from bit-flip city
                // to recover objective value
                for (r = origBF; r < nbCities; r++) {
                    wc = sol.weightAcc[r] + deltaW;
                    ft += ttp.distFor(tour[r] - 1, tour[(r + 1) % nbCities] - 1) / (maxSpeed - wc * C);
                }
                // compute recovered objective value
                G = fp - ft * R;

                //=====================================
                // update if improvement or
                // Boltzmann condition satisfied
                //=====================================
                double mu = Math.random();
                double energy_gap = G - GBest;


                Double t2 = L[origBF]*(1/(maxSpeed-C*(ttp.weightOf(k)+sol.weightAcc[origBF])) - 1/maxSpeed);//计算当前物品选取与否的时间差
                Double tempP= ttp.profitOf(k)-R*t2;//当前物品的增值

                boolean acceptance = energy_gap > 0 || Math.exp(energy_gap / T) > mu;

                if (acceptance == true && energy_gap<0 ){
                    if (origBF<nbCities/2 && pickingPlan[k] == 0
                            && ttp.weightOf(k)>sol.averageWeight
//                            && tempP > sol.standScore
//                            && Math.random()>0.1
                    ){
                        acceptance = false;
                        continue;
                    }
                }

//                if (energy_gap>0){
//                    //如果在前面的较重价高物品被选取,则倾向于把该物品扔掉
//                    if (origBF<nbCities/2 && pickingPlan[k] == 0
//                            && ttp.weightOf(k)>sol.averageWeight
//                            && tempP > sol.standScore
//                            && Math.random()>0.1){
//                        acceptance = false;
//                    }
//                }
//                else {
//                    //如果在后面的较轻价高物品未被选取,则倾向于抛弃该物品
//                    if(origBF>nbCities/2 && pickingPlan[k] != 0
//                            && ttp.weightOf(k)<sol.averageWeight
//                            && tempP > sol.standScore
//                            && Math.random()>0.1){
//                        acceptance = true;
//
//                    }
//                }


                if (acceptance) {

                    GBest = G;

                    // bit-flip
                    pickingPlan[k] = pickingPlan[k] != 0 ? 0 : A[k];

                    //===========================================================
                    // recover accumulation vectors
                    //===========================================================
                    if (pickingPlan[k] != 0) {
                        deltaP = ttp.profitOf(k);
                        deltaW = ttp.weightOf(k);
                    } else {
                        deltaP = -ttp.profitOf(k);
                        deltaW = -ttp.weightOf(k);
                    }
                    fp = sol.fp + deltaP;
                    origBF = sol.mapCI[A[k] - 1];
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

                }

            }

            // update best if improvement
            if (sol.ob > sBest.ob) {
                sBest = sol.clone();
            }

            if (this.debug) {
                Deb.echo(">> KRP " + nbIter + ": ob=" +
                        String.format("%.0f",sol.ob));
            }

            // cool down temperature
            T = T * alpha;

            // stop when temperature reach absolute value
        } while (T > T_abs);


        //我的代码
       sBest = DeleteAndInsert(sBest);

        if (this.debug) {
            Deb.echo(">> DeleteAndInsertKRP " + nbIter + ": ob=" +
                    String.format("%.0f",sBest.ob));
        }



        //原始代码
        // in order to recover all history vector
        ttp.objective(sBest);

        return sBest;
    }

    public void initTourAndPackPlan(int iterN){
        int nbCities = ttp.getNbCities();
        if (tourS5 == null){
            tourS5 = new int[nbCities];
        }

        tourS5= Optimisation.linkernTour(ttp.file.getPath(), ttp.numberOfNodes+1,ttp.getTspName());
        if (construct==null){
            construct = new Constructive(ttp);
        }
        if (init==null){
            init = new Initialization(ttp);
        }

        // use Lin-Kernighan to initialize the tour

        // pre-process the knapsack
        // insert and eliminate items
//            if (ttp.getNbCities() < 30000) s0 = insertAndEliminate(s0);
//            else s0 = insertT2(s0);

        if (CLKresult == null){
            CLKresult = new int[ttp.numberOfNodes];
        }

        for (int i = 0; i < tourS5.length-1; i++) {
            CLKresult[i] = tourS5[i]+1;
        }

        if (iterN == 0){
            s0 = new TTPSolution(
                    CLKresult,
                    construct.zerosPickingPlan()
            );
            //Insert算法
//            s0 = insertT2(s0);
        }else {

            s0 = new TTPSolution(
                init.rlinkern(),
//                    CLKresult,
//                construct.linkernTour(),
//                s0.getPickingPlan()
                    construct.zerosPickingPlan()
            );
//            s0 = insertT2(s0);
        }
//            if (ttp.getNbCities() < 30000)
                s0 = insertAndEliminate(s0);
        if (s0.ob<0){
//            s0 = insertT2(s0);
        }
//            else s0 = insertT2(s0);

//        int[] Insertpp = new int[ttp.getNbItems()];
//        Insertpp = s0.getPickingPlan();

        //PackIterative算法
//        s0 = new TTPSolution(
//                tourS5,
//                s0.getPickingPlan()
//        );
//        ttp.adjustS5Evaluate(s0, !false);
//        s0 = Optimisation.HT(ttp, 600*1000, true, s0);


//        int[] pickingPlan = s0.getPickingPlan().clone();


//        Initialization init = new Initialization(ttp);
//        s0 = new TTPSolution(
//                init.rlinkern(),
////                construct.linkernTour(),
//                s0.getPickingPlan()
//        );
        ttp.CocoObjective(s0);
    }



    @Override
    public TTPSolution search() {

        //===============================================
        // generate initial solution
        //===============================================
//        Constructive construct = new Constructive(ttp);
//        Initialization init = new Initialization(ttp);

        tspImprove = false;

        int nbCities = ttp.getNbCities();
        if (tourS5 == null) {
            tourS5 = new int[nbCities];
        }

        tourS5 = Optimisation.linkernTour(ttp.file.getPath(), ttp.numberOfNodes + 1, ttp.getTspName());
        if (construct == null) {
            construct = new Constructive(ttp);
        }
        if (init == null) {
            init = new Initialization(ttp);
        }

        // use Lin-Kernighan to initialize the tour

        // pre-process the knapsack
        // insert and eliminate items
//            if (ttp.getNbCities() < 30000) s0 = insertAndEliminate(s0);
//            else s0 = insertT2(s0);

        if (CLKresult == null) {
            CLKresult = new int[ttp.numberOfNodes];
        }

        for (int i = 0; i < tourS5.length - 1; i++) {
            CLKresult[i] = tourS5[i] + 1;
        }

        List<String> list = new ArrayList<String>();
        Boolean isGoThroughAll=false,isRepeat=false;
        s0 = new TTPSolution(
//                    init.rlinkern(),
//                    CLKresult,
                construct.linkernTour(),
//                s0.getPickingPlan()
                construct.zerosPickingPlan()
        );
        StringBuilder builder = new StringBuilder();
        for (int value : s0.getTour()){
            builder.append(value);
        }
        list.add(builder.toString());

        s0.scoresArr = construct.zerosScores();
        s0.L = construct.zerosDistance();
        s0 = insertAndEliminate(s0);

        ttp.objective(s0);
        if (this.debug) {
            Deb.echo("STARTING SOL >> " + s0.ob);
        }
        // copy initial solution into improved solution
        TTPSolution sol = s0.clone();
        // best found solution
        TTPSolution sBest = sol.clone();
        // best found
        double GBest = sol.ob;
        // number of iterations
        int nbIter = 0;
        // improvement tag
        boolean improved;

        //===============================================
        // start cosolver search
        //===============================================
        do {


            nbIter++;
            do {
                nbIter++;
                improved = false;

                // 2-opt heuristic on TSKP
                sol = fast2optPGCH(sol);
                tspImprove = sol.tspImprove;
                if (!tspImprove) {
//                    sol = fast2optPGCHNoImproveTSP(sol);
                }
                // simple bit-flip on KRP
//                sol = lsBitFlip(sol);
                sol = simulatedAnnealing(sol);
                // update best if improvement
                if (sol.ob > GBest) {
                    GBest = sol.ob;
                    improved = true;
                    Deb.echo("innerOb-best: " + sol.ob);
                }

                if (Thread.currentThread().isInterrupted()) {
                    if (sol.ob > sBest.ob) sBest = sol.clone();
                    return sBest;
                }
            } while (improved);
            //===============================================

            // tag idle step
            if (sol.ob > sBest.ob) {
                sBest = sol.clone();


                // debug msg
                if (this.debug) {
                    Deb.echo("Best " + nbIter + ":");
                    Deb.echo("outob-best: " + sol.ob);
                    Deb.echo("wend   : " + sol.wend);
                    Deb.echo("---");
                }
            }

            // stop execution if interrupted
            if (Thread.currentThread().isInterrupted()) return sBest;

            // restart if no improvement
            if (!Thread.currentThread().isInterrupted()) {
                if (debug) {
                    Deb.echo("===> RESTART");
                }

                // restart
                int[] newTour = init.rlinkern();
                StringBuilder builder2 = new StringBuilder();
                for (int value : newTour){
                    builder2.append(value);
                }
                String strTour = builder2.toString();
                while (list.contains(strTour)){
                     newTour = init.rlinkern();
                    StringBuilder builder3 = new StringBuilder();
                    for (int value : newTour){
                        builder3.append(value);
                    }
                    strTour = builder3.toString();
                }
                list.add(strTour);

                sol = new TTPSolution(
                        newTour,
//                        construct.linkernTour(),
//                        sol.getTour(),

                        construct.zerosPickingPlan()
//                        sol.getPickingPlan()
                );
                sol = insertAndEliminate(sol);
            }
            sol.scoresArr = construct.zerosScores();

            sol = insertAndEliminate(sol);

            ttp.objective(sol);



            // stop when time expires
        } while (!Thread.currentThread().isInterrupted());
        //===============================================

        return sBest;
    }

    public TTPSolution fast2optPGCHMyImproved(TTPSolution sol) {

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
        boolean improved,TotalImproved;

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
            TTPSolution sInnerBest;
            sInnerBest = sol.clone();
            improved = false;
            nbIter++;
            TotalImproved=false;
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

                        deltaT += -sol.timeRec[q] + ttp.distFor(c1, c2) / (maxSpeed - wc * C);
                    }
                    // retrieve neighbor's final time
                    ft = ft + deltaT;

                    // update best
                    if (ft - ftBest < threshold) { // soft condition
                        iBest = i;
                        jBest = j;
                        ftBest = ft;
                        improved = true;
//                        ifImproved = true;

                    }else {
//                        sol = PGCH(sol, i, j);
//                        tour = sol.getTour();
//                        // update best
//                        if (sInnerBest.ob<sol.ob){
//                            TotalImproved = true;
//                        }
                    }



                    //if (firstfit && improved) break;
                } // END FOR j
            } // END FOR i

            if (improved) {

                // apply 2-opt move
                TwoOptHelper.do2opt(tour, iBest, jBest);

                // evaluate & update vectors
                ttp.objective(sol);
            }
            if (TotalImproved && !improved){
//                sol = PGCH(sol,iTemp,jTemp);
                ttp.objective(sol);
            }


            // debug msg
            if (this.debug) {
                Deb.echo(">> TSKP " + nbIter +
                        ": ob=" + String.format("%.0f", sol.ob) +
                        " | ft=" + String.format("%.0f", sol.ft));
            }

        } while (improved || TotalImproved);

        if (debug) Deb.echo("==> 2-opt :" + nbIter + " iterations");

        // in order to compute sol.timeAcc
        // we need to use objective function
        ttp.objective(sol);

        return sol;
    }

    public TTPSolution fast2optPGCH(TTPSolution sol) {

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
        boolean improved,TotalImproved,ifImproved;

        // best solution
        ttp.objective(sol);

        int iBest=0, jBest=0;
        int iTemp=0, jTemp=0;
        double ftBest = sol.ft;

        // neighbor solution
        double ft;
        long wc;
        int i, j, c1, c2, q;
        int nbIter = 0;

        // current tour

        int[] pickingPlan = sol.getPickingPlan();


        TTPSolution sBest = sol.clone();
        tour = sol.getTour();
        // best solution
        double GBest = sol.ob;
        GTSPBest = sol.ob;

        double threshold = -0.1;
        if (nbItems >= 100000) {
            threshold = -10;
        }
        if (nbCities >= 50000) { // ex. pla85000 based instances
            threshold = -1000;
        }
        TTPSolution sInnerBest;
        ifImproved = false;
        // search
        do {
             sInnerBest = sol.clone();
            improved = false;
            nbIter++;

            // cleanup and stop execution if interrupted
            if (Thread.currentThread().isInterrupted()) {
                break;
            }

            // fast 2-opt+PGCH
            for (i = 1; i < nbCities - 1; i++) {
                int node1 = tour[i] - 1;
                for (int node2 : candidates[node1]) {
                    j = sol.mapCI[node2];

//                    for (int z : ttp.clusters[node1]) { //根据tour处在当前位置城市中的物品下标
//                        int itemIndex = j;
//                        float PR = (float) sol.scoresArr[itemIndex];
//                    }


                    // calculate final time with partial delta
                    ft = sol.ft;
                    wc = i - 2 < 0 ? 0 : sol.weightAcc[i - 2]; // fix index...
                    deltaT = 0;
//                    if (i>=j) continue;
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
                        ifImproved = true;

                    }

                } // END FOR j

            } // END FOR i

//            sol.setTour(sBest.getTour());
//            sol.setPickingPlan(sBest.getPickingPlan());
//            sol.ob = sBest.ob;
//            ttp.CocoObjective(sol);
            //===================================
            // update if improvement
            //===================================
            if (improved) {

                // apply 2-opt move
                TwoOptHelper.do2opt(tour, iBest, jBest);

                // evaluate & update vectors
                ttp.objective(sol);
                GBest = sol.ob;
            }
//            if (TotalImproved && !improved){
//                sol = PGCH(sol,iTemp,jTemp);
//                ttp.CocoObjective(sol);
//            }

            // debug msg
            if (this.debug) {
                Deb.echo(">> TSKP " + nbIter +
                        ": ob=" + String.format("%.0f", sol.ob) +
                        " | ft=" + String.format("%.0f", sol.ft));
            }



        } while (improved) ;

//        if (debug) Deb.echo("==> 2-opt :" + nbIter + " iterations");

        // in order to compute sol.timeAcc
        // we need to use objective function
        sol.tspImprove = ifImproved;
        ttp.objective(sol);
        CalculateScoreArr(sol);
        return sol;
    }
    public void CalculateDistance(TTPSolution s){
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        int[] A = ttp.getAvailability();
        int[] tour = s.getTour();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity; // velocity const
        // neighbor solution
        int origBF,k,i;

        // distances of all tour cities (city -> end)
        long[] L = new long[nbCities];

        // time approximations
        double t2;
        L[nbCities-1] = ttp.distFor(tour[nbCities-1] - 1,0);
        for (i=nbCities-2; i >= 0; i--) {
            L[i] = L[i+1] + ttp.distFor(tour[i+1]-1, tour[i]-1);
        }
        s.L = L;
    }



    public TTPSolution fast2optPGCHNoImproveTSP(TTPSolution sol) {

        // TTP data
        int nbCities = ttp.getNbCities();
        int nbItems = ttp.getNbItems();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;

        // initial solution data
        int[] tour;


        CalculateDistance(sol);
        CalculateScoreArr(sol);

        // delta parameters
        double deltaT;

        // improvement indicator
        boolean improved,TotalImproved,ifImproved;

        // best solution
        ttp.objective(sol);

        int iBest=0, jBest=0;
        int iTemp=0, jTemp=0;
        double ftBest = sol.ft;

        // neighbor solution
        double ft;
        long wc;
        int i, j, c1, c2, q;
        int nbIter = 0;

        // current tour

        int[] pickingPlan = sol.getPickingPlan();


        TTPSolution sBest = sol.clone();
        tour = sol.getTour();
        // best solution
        double GBest = sol.ob;
        GTSPBest = sol.ob;

        double threshold = -0.1;
        if (nbItems >= 100000) {
            threshold = -10;
        }
        if (nbCities >= 50000) { // ex. pla85000 based instances
            threshold = -1000;
        }
        TTPSolution sInnerBest;
        ifImproved = false;
        float[] LPPR = sol.getLeastPickedPR();;
        float[] HUPR = sol.getHighestUnpickedPR();
        float[] preMinV = sol.getPrefixMinimumValue();
        float[] posMaxV = sol.getPostfixMaximumValue();

        // search
        do {
            sInnerBest = sol.clone();
            improved = false;
            TotalImproved=false;
            nbIter++;

            // cleanup and stop execution if interrupted
            if (Thread.currentThread().isInterrupted()) {
                break;
            }

            // fast 2-opt+PGCH
            for (i = 1; i < nbCities - 1; i++) {
                int node1 = tour[i] - 1;

                int restoreNode1 = 0,restoreNode2=0;
                float absValue1= (float) 0.0;float absValue2= (float) 0.0;
                List<Integer> list = new ArrayList<Integer>();

                for (int node22 : candidates[node1]) { //根据tour处在当前位置城市中的物品下标
                    j = sol.mapCI[node22];
//                    if (absValue1 < (HUPR[mu] - LPPR[i]) ){
//                        absValue1 = (HUPR[mu] - LPPR[i]);
//                        restoreNode1 = mu;
//                    }
//                }
//
////                for (int node2 : candidates[node1]) {//----------------------
//                j = sol.mapCI[restoreNode1];
                    // calculate final time with partial delta
                    ft = sol.ft;
                    wc = i - 2 < 0 ? 0 : sol.weightAcc[i - 2]; // fix index...
                    deltaT = 0;
                    if (i >= j) continue;


                    for (q = i - 1; q <= j; q++) {//两点之间有几个城市

                        wc += TwoOptHelper.get2optValue(q, sol.weightRec, i, j);
                        c1 = TwoOptHelper.get2optValue(q, tour, i, j) - 1;
                        c2 = TwoOptHelper.get2optValue((q + 1) % nbCities, tour, i, j) - 1;

                        deltaT += -sol.timeRec[q] + ttp.distFor(c1, c2) / (maxSpeed - wc * C);
                    }
                    // retrieve neighbor's final time
                    ft = ft + deltaT;

                    // update best
                    if (ft - ftBest < threshold) { // soft condition
                        iBest = i;
                        jBest = j;
                        ftBest = ft;
                        improved = true;
                        ifImproved = true;

                    } else {
                        if ((ft - ftBest) < (int) Math.floor(Math.sqrt(ftBest))) {
                            sol = PGCH(sol, i, j);
                        }
                    }
                }
                 // END FOR j

                if (sInnerBest.ob<sol.ob){
                    TotalImproved = true;
                }

//                }//------------
            } // END FOR i

            //===================================
            // update if improvement
            //===================================
            if (improved) {

                // apply 2-opt move
                TwoOptHelper.do2opt(tour, iBest, jBest);

                // evaluate & update vectors
                ttp.objective(sol);
                GBest = sol.ob;
            }
            if (TotalImproved && !improved){
//                sol = PGCH(sol,iTemp,jTemp);
                ttp.objective(sol);
            }

            // debug msg
            if (this.debug) {
                Deb.echo(">> TSKPP " + nbIter +
                        ": ob=" + String.format("%.0f", sol.ob) +
                        " | ft=" + String.format("%.0f", sol.ft));
            }



        } while (improved || TotalImproved) ;

//        if (debug) Deb.echo("==> 2-opt :" + nbIter + " iterations");

        // in order to compute sol.timeAcc
        // we need to use objective function
        sol.tspImprove = ifImproved;
        ttp.objective(sol);
        CalculateDistance(sol);
        CalculateScoreArr(sol);
        return sol;
    }

    public Boolean IFPGCH(TTPSolution s, int m, int n) {
//        ttp.CocoObjective(sol);

        int[] A = ttp.getAvailability();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;
        double R = ttp.getRent();

        // best solution
//        double GBest = sol.ob;
        TTPSolution sol = s.clone();
//        ttp.objective(sBest);

        int nbCities = ttp.getNbCities();
        int[] tour = sol.getTour();
        int [] myTour = TwoOptHelper.mydo2opt(tour, m, n);
        // evaluate & update vectors
        sol.setTour(myTour);

        // initial solution data
        int[] x = sol.getTour();
        int[] z = sol.getPickingPlan();
        int[] pickingPlan = sol.getPickingPlan();
        float[] LPPR = sol.getLeastPickedPR();;
        float[] HUPR = sol.getHighestUnpickedPR();
        float[] preMinV = sol.getPrefixMinimumValue();
        float[] posMaxV = sol.getPostfixMaximumValue();


        // best solution
        int kBest = 0;
        double GBest = sol.ob;

        // neighbor solution
        int origBF;
        int  k, itr;

        // neighbor solution
        long fp;
        double ft, G;
        long wc;
        int r;
        int deltaP, deltaW;
        boolean improved = false;

        // visit part cities
        for (int i = m; i < n+1; i++){
            for (int y = 0; y < (ttp.PR_clusters[x[i] - 1]).size(); y++) {    //根据tour处在当前位置城市中的物品下标:j
//                float PR = ttp.PR_clusters[x[i] - 1].get(y);
                 k = ttp.clusters[x[i] - 1].get(y);
//                if (z[k] != 0 && (PR < preMinV[i])){    //必须丢弃
//                    z[k] = 0;
//                    sol.wend += ttp.weightOf(k);
//                }


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
                origBF = sol.mapCI[A[k] - 1];//表示已经走过的城市个数，nbCities-OrigBF表示剩下没走的城市数mCities-OrigBF表示剩下没走的城市数

                // starting time
                ft = origBF == 0 ? 0 : sol.timeAcc[origBF - 1];

                // recalculate velocities from bit-flip city
                for (r = origBF; r < nbCities; r++) {
                    wc = sol.weightAcc[r] - ttp.weightOf(k);;
                    ft += ttp.distFor(tour[r]-1, tour[(r + 1) % nbCities]-1) / (maxSpeed - wc * C);
                }//去掉该物品后所花费的时间

                G = fp - ft * R;

                // update best
                if (G > GBest) {

                    kBest = k;
                    GBest = G;
                    // bit-flip
                    pickingPlan[k] = pickingPlan[k] != 0 ? 0 : A[k];

                    //===========================================================
                    // recover accumulation vectors
                    //===========================================================
                    if (pickingPlan[k] != 0) {
                        deltaP = ttp.profitOf(k);
                        deltaW = ttp.weightOf(k);
                    } else {
                        deltaP = -ttp.profitOf(k);
                        deltaW = -ttp.weightOf(k);
                    }
                    fp = sol.fp + deltaP;
                    origBF = sol.mapCI[A[k] - 1];
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
                }



            }
        }


        sol.setPickingPlan(z);
        ttp.CocoObjective(sol);

        if (sol.ob>s.ob){
            GTSPBest = sol.ob;
        }
        return sol.ob>s.ob;
    }


    public TTPSolution PGCH(TTPSolution s, int m, int n) {
//        ttp.CocoObjective(sol);

        int[] A = ttp.getAvailability();
        double maxSpeed = ttp.getMaxSpeed();
        double minSpeed = ttp.getMinSpeed();
        long capacity = ttp.getCapacity();
        double C = (maxSpeed - minSpeed) / capacity;
        double R = ttp.getRent();

        // best solution
//        double GBest = sol.ob;
        TTPSolution sol = s.clone();
//        ttp.objective(sBest);

        int nbCities = ttp.getNbCities();
        int[] tour = sol.getTour();
        int [] myTour = TwoOptHelper.mydo2opt(tour, m, n);
        // evaluate & update vectors
        sol.setTour(myTour);

        // initial solution data
        int[] x = sol.getTour();
        int[] z = sol.getPickingPlan();
        int[] pickingPlan = sol.getPickingPlan();
        float[] LPPR = sol.getLeastPickedPR();;
        float[] HUPR = sol.getHighestUnpickedPR();
        float[] preMinV = sol.getPrefixMinimumValue();
        float[] posMaxV = sol.getPostfixMaximumValue();


        // visit part cities
        for (int i = m; i < n+1; i++){
            for (int k = 0; k < (ttp.clusters[x[i] - 1]).size(); k++) {    //根据tour处在当前位置城市中的物品下标:j

                int itemIndex = ttp.clusters[x[i] - 1].get(k);
                if (z[itemIndex] != 0 ){    //必须丢弃
                    z[itemIndex] = 0;
                    sol.wend += ttp.weightOf(itemIndex);
                }
            }
        }
        ttp.objective(sol);
        CalculateScoreArr(sol);
//        sol.setPickingPlan(z);

        for (int i = n; i >m-1; i--){
            for (int k = (ttp.clusters[x[i] - 1]).size()-1; k >-1 ; k--) {    //根据tour处在当前位置城市中的物品下标:j
                int itemIndex = ttp.clusters[x[i] - 1].get(k);
                double PR =  sol.scoresArr[itemIndex];

                if ((ttp.weightOf(itemIndex) <= sol.wend) && z[itemIndex] == 0 && (PR > sol.standScore)){ //可以拾取
                    z[itemIndex] = A[itemIndex];
                    sol.wend -= ttp.weightOf(itemIndex);
                }
            }
        }

        sol.setPickingPlan(z);
        ttp.objective(sol);
        if (sol.ob>s.ob){
            s = sol.clone();
//            s.setTour(x);
//            s.setPickingPlan(z);
//            s.ob=sol.ob;

        }
        return s;
    }
    public TTPSolution BBFS(TTPSolution sol) {

        // copy initial solution into improved solution
        TTPSolution sBest = sol.clone();

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
        int k, r;
//        int nbIter = 0;

//        double T = T0;

        long trials = Math.round(nbItems*trialFactor);

//        if (debug) Deb.echo(">>>> TRIAL FACTOR: "+trialFactor);

        //===============================================
        // start simulated annealing process
        //===============================================
        ttp.CocoUpdateIndexFunction(sol);
        ttp.MarkUncheckedAll(sol);
        Random rm = new Random();
        List<Integer> list = new ArrayList<Integer>();
        while (!ttp.AllChecked(sol)){
            int tempRandomIndex = 0;

            pickingPlan = sol.getPickingPlan();

            int boundItemCount = sol.boundaryItemsArr.size();

            tempRandomIndex = rm.nextInt(boundItemCount);
            Item tempItem = sol.boundaryItemsArr.get(tempRandomIndex);

            //不重复选取物品
            if (!list.contains(tempItem.itemId) && (!tempItem.isChecked)){
                list.add(tempItem.itemId);
                ttp.MarkChecked(sol,tempItem.itemId);


                // check if new weight doesn't exceed knapsack capacity
                if (pickingPlan[tempItem.itemId] == 0 && ttp.weightOf(tempItem.itemId) > sol.wend) continue;

                // calculate deltaP and deltaW
                if (pickingPlan[tempItem.itemId] == 0) {
                    deltaP = ttp.profitOf(tempItem.itemId);
                    deltaW = ttp.weightOf(tempItem.itemId);
                } else {
                    deltaP = -ttp.profitOf(tempItem.itemId);
                    deltaW = -ttp.weightOf(tempItem.itemId);
                }
                fp = sol.fp + deltaP;

                // handle velocity constraint
                // index where Bit-Flip happened
                origBF = sol.mapCI[A[tempItem.itemId] - 1];
                // starting time
                ft = origBF == 0 ? .0 : sol.timeAcc[origBF - 1];
                // recalculate velocities from bit-flip city
                // to recover objective value
                for (r = origBF; r < nbCities; r++) {
                    wc = sol.weightAcc[r] + deltaW;
                    ft += ttp.distFor(tour[r] - 1, tour[(r + 1) % nbCities] - 1) / (maxSpeed - wc * C);
                }
                // compute recovered objective value
                G = fp - ft * R;

                //=====================================
                // update if improvement or
                //=====================================
                double energy_gap = G - GBest;
                boolean acceptance = energy_gap > 0 ;

                if (acceptance) {

                    GBest = G;
                    // bit-flip
                    pickingPlan[tempItem.itemId] = pickingPlan[tempItem.itemId] != 0 ? 0 : A[tempItem.itemId];

                    //===========================================================
                    // recover accumulation vectors
                    //===========================================================
                    if (pickingPlan[tempItem.itemId] != 0) {
                        deltaP = ttp.profitOf(tempItem.itemId);
                        deltaW = ttp.weightOf(tempItem.itemId);
                    } else {
                        deltaP = -ttp.profitOf(tempItem.itemId);
                        deltaW = -ttp.weightOf(tempItem.itemId);
                    }
                    fp = sol.fp + deltaP;
                    origBF = sol.mapCI[A[tempItem.itemId] - 1];
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
                    sol.objectiveScore = G;
                    sol.fp = fp;
                    sol.ft = ft;
                    sol.wend = capacity - sol.weightAcc[nbCities - 1];
                    //===========================================================

                    sol.setPickingPlan(pickingPlan);
//                    ttp.objective(sol);
                    ttp.CocoUpdateIndexFunction(sol);
                    list.clear();
                }
            }
        }
        if (this.debug) {
            Deb.echo(">> KRP "  + ": ob=" +
                    String.format("%.0f",sol.ob));
        }

        // in order to recover all history vector
        ttp.CocoObjective(sol);

        return sol;
    }






    public static double generateTFLinFitx(int xi) {
        int i=-1;
//    int[] x = new int[]{         50,    204,  609,  1147,  8034, 38875, 105318,  253568, 338090};
        int[] x = new int[]{         1,   75,    375,  790,  2102, 15111, 70250, 140500,  338090};
        double[] y = new double[]{57872,  13896,  700,   350,    16,     1,   0.16,  0.0493,   0.03};

//    int[] x = new int[]{         1,   75,    375,  790,  2102, 15111, 70250, 140500,  338090};
//    double[] y = new double[]{
//      57872, 13896, 350, 16, 1, .16, .0493, .03
//    };

        int n=y.length;
        for (int k=0; k<n-1; k++) {
            if (x[k] <= xi && xi < x[k + 1]) {
                i = k;
                break;
            }
        }
        if (xi <= x[0]) {
            return 57872.0;
        }
        if (xi >= x[n-1]) {
            return 0.03;
        }

        double m = ( y[i]-y[i+1] ) / ( x[i]-x[i+1] );
        double b = y[i]-m*x[i];

        double yi = m*xi + b;

        return yi;
    }

    public static double generateTFLinFit(int xi) {
        int i=-1;
//    int[] z = new int[]{ 1,  75, 375, 790, 2102, 15111, 70250, 140500, 338090};

//    int[] x = new int[]{         1,   75,    375,  790,  2102, 15111, 70250, 140500,  338090};
//    double[] y = new double[]{57872,  13896,  700,   350,    16,     1,   0.16,  0.0493,   0.03};

        int[] x = new int[]{ 1, 130, 496, 991, 3038, 18512, 75556, 169046, 338090};
        double[] y = new double[]{57872,  13896,  700,   350,    16,     1,   0.16,  0.0493,   0.03};


        int n=y.length;
        for (int k=0; k<n-1; k++) {
            if (x[k] <= xi && xi < x[k + 1]) {
                i = k;
                break;
            }
        }
        if (xi <= x[0]) {
            return 57872.0;
        }
        if (xi >= x[n-1]) {
            return 0.03;
        }

        double m = ( y[i]-y[i+1] ) / ( x[i]-x[i+1] );
        double b = y[i]-m*x[i];

        double yi = m*xi + b;

        return yi;
    }

    public static double generateTFExpFit(int X) {
        double a =   1.614e+05;
        double b =   -0.006915;
        double c =       145.4;
        double d =  -5.191e-05;

        double A = a * Math.exp(b * X) + c * Math.exp(d * X);

        return A;
    }

    public static double generateTFManualFit(int X) {
        int nbItems = X;
        double trialFactor;
        if (nbItems < 500)
            trialFactor = 1000;// reduce... (time)

        else if (nbItems < 1000)
            trialFactor = 100;
        else if (nbItems < 5000)
            trialFactor = 50;

        else if (nbItems < 20000)
            trialFactor = 10; //was 5... retest others
        else if (nbItems < 100000)
            trialFactor = 1;
        else if (nbItems < 200000)
            trialFactor = .04;
        else
            trialFactor = .03;
        return trialFactor;
    }

}
