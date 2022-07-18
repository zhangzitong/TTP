package ttp;

import utils.CityCoordinates;

import java.io.File;
import java.util.ArrayList;


public abstract class TTPInstance {

    protected String name;
    protected String directory;
    protected String tspName;

    protected String knapsackDataType;
    protected int nbCities;
    protected int nbItems;
    protected long capacity;
    public double minSpeed;
    public double maxSpeed;
    protected String edgeWeightType;
    protected CityCoordinates[] coordinates;
    protected long[][] dist = null;
    public int[] availability;
    protected int[] profits;
    protected int[] weights;

    protected File ttpFile;

    // item clusters per city
    public ArrayList<Integer>[] clusters;


    public ArrayList<Float>[] getPR_clusters() {
        return PR_clusters;
    }

    public void setPR_clusters(ArrayList<Float>[] PR_clusters) {
        this.PR_clusters = PR_clusters;
    }

    public ArrayList<Float>[] PR_clusters;
/*
S5部分
 */
    public String problemName;
    public int numberOfNodes;
    public int numberOfItems;
    public long capacityOfKnapsack;
    public double rentingRatio;
    public double[][] nodes;
    public int[][] items;

    ArrayList<Float>pR;
    public double[][] profitabilityRatio;
    public File file;

    public void setName(String name) {
        this.name = name;
    }

    public void setDirectory(String directory) {
        this.directory = directory;
    }

    public void setTspName(String tspName) {
        this.tspName = tspName;
    }

    public String getKnapsackDataType() {
        return knapsackDataType;
    }

    public void setKnapsackDataType(String knapsackDataType) {
        this.knapsackDataType = knapsackDataType;
    }

    public void setNbCities(int nbCities) {
        this.nbCities = nbCities;
    }

    public void setNbItems(int nbItems) {
        this.nbItems = nbItems;
    }

    public void setCapacity(long capacity) {
        this.capacity = capacity;
    }

    public void setMinSpeed(double minSpeed) {
        this.minSpeed = minSpeed;
    }

    public void setMaxSpeed(double maxSpeed) {
        this.maxSpeed = maxSpeed;
    }

    public String getEdgeWeightType() {
        return edgeWeightType;
    }

    public void setEdgeWeightType(String edgeWeightType) {
        this.edgeWeightType = edgeWeightType;
    }

    public void setCoordinates(CityCoordinates[] coordinates) {
        this.coordinates = coordinates;
    }

    public void setAvailability(int[] availability) {
        this.availability = availability;
    }

    public void setProfits(int[] profits) {
        this.profits = profits;
    }

    public void setWeights(int[] weights) {
        this.weights = weights;
    }

    public File getTtpFile() {
        return ttpFile;
    }

    public void setTtpFile(File ttpFile) {
        this.ttpFile = ttpFile;
    }

    public void setClusters(ArrayList<Integer>[] clusters) {
        this.clusters = clusters;
    }

    public String getProblemName() {
        return problemName;
    }

    public void setProblemName(String problemName) {
        this.problemName = problemName;
    }

    public int getNumberOfNodes() {
        return numberOfNodes;
    }

    public void setNumberOfNodes(int numberOfNodes) {
        this.numberOfNodes = numberOfNodes;
    }

    public int getNumberOfItems() {
        return numberOfItems;
    }

    public void setNumberOfItems(int numberOfItems) {
        this.numberOfItems = numberOfItems;
    }

    public long getCapacityOfKnapsack() {
        return capacityOfKnapsack;
    }

    public void setCapacityOfKnapsack(long capacityOfKnapsack) {
        this.capacityOfKnapsack = capacityOfKnapsack;
    }

    public double getRentingRatio() {
        return rentingRatio;
    }

    public void setRentingRatio(double rentingRatio) {
        this.rentingRatio = rentingRatio;
    }

    public double[][] getNodes() {
        return nodes;
    }

    public void setNodes(double[][] nodes) {
        this.nodes = nodes;
    }

    public int[][] getItems() {
        return items;
    }

    public void setItems(int[][] items) {
        this.items = items;
    }

    public ArrayList<Float> getpR() {
        return pR;
    }

    public void setpR(ArrayList<Float> pR) {
        this.pR = pR;
    }

    public double[][] getProfitabilityRatio() {
        return profitabilityRatio;
    }

    public void setProfitabilityRatio(double[][] profitabilityRatio) {
        this.profitabilityRatio = profitabilityRatio;
    }

    public File getFile() {
        return file;
    }

    public void setFile(File file) {
        this.file = file;
    }



    @Override
    public String toString() {
        String s = "";
        if (nbCities < 20 && nbItems < 20) {
            // coordinates
            s += "cities | coordinates:\n";
            for (int i=0; i<this.nbCities; i++) {
                s += String.format("%6d | ", i+1) + this.coordinates[i] + "\n";
            }
            s+="\n";
            // distance matrix
            s += "distance matrix:\n";
            for (int i=0; i<this.nbCities; i++) {
                for (int j=0; j<this.nbCities; j++) {
                    s += String.format("%5d", this.getDist()[i][j]);
                }
                s += "\n";
            }
            s+="\n";

            // items
            s += "items   : ";
            for (int i=0; i<this.nbItems; i++) {
                s +=  String.format("%5d", i+1);
            }
            s += "\n";
            s += "values  : ";
            for (int i=0; i<this.nbItems; i++) {
                s +=  String.format("%5d", this.profits[i]);
            }
            s += "\n";
            s += "weights : ";
            for (int i=0; i<this.nbItems; i++) {
                s +=  String.format("%5d", this.weights[i]);
            }
            s += "\n";
            s += "city ref: ";
            for (int i=0; i<this.nbItems; i++) {
                s +=  String.format("%5d", this.availability[i]);
            }
            s+="\n\n";

        }

        s += "name     : " + this.name + "\n";
        s += "#cities  : " + this.nbCities + "\n";
        s += "#items   : " + this.nbItems + "\n";
        s += "capacity : " + this.capacity + "\n";
        s += "min speed: " + this.getMinSpeed() + "\n";
        s += "max speed: " + this.getMaxSpeed() + "\n";

        return s;
    }

    public void setDist(long[][] dist) {
        this.dist = dist;
    }

    public String getTspName() {
        return tspName;
    }
    public String getDirectory() {
        return directory;
    }
    public String getName() {
        return name;
    }
    public long[][] getDist() {
        return dist;
    }
    public int[] getAvailability() {
        return availability;
    }
    public int[] getWeights() {
        return weights;
    }
    public int[] getProfits() {
        return profits;
    }
    public int getNbCities() {
        return nbCities;
    }
    public int getNbItems() {
        return nbItems;
    }
    public double getMaxSpeed() {
        return maxSpeed;
    }
    public double getMinSpeed() {
        return minSpeed;
    }
    public long getCapacity() {
        return capacity;
    }
    public ArrayList<Integer>[] getClusters() {
        return clusters;
    }
    public CityCoordinates[] getCoordinates() {
        return coordinates;
    }

    public int profitOf(int i) {
        return this.profits[i];
    }
    public int weightOf(int i) {
        return this.weights[i];
    }

    // used to simulate the distance matrix
    public long distFor(int i, int j) {
        if (dist==null) {
//      return Math.round(this.coordinates[i].distanceEuclid(this.coordinates[j]));
            return (long)Math.ceil(this.coordinates[i].distanceEuclid(this.coordinates[j]));
        }
        return dist[i][j];
    }


    /**
     * organize items per city
     */
    public void clusterItems() {

        clusters = new ArrayList[nbCities];
        PR_clusters = new ArrayList[nbCities];
        int i;

        // init cluster arrays
        for (i=0; i<nbCities; i++) {
            clusters[i] = new ArrayList<>();
            PR_clusters[i] = new ArrayList<>();
        }
        // fill clusters
        for (i=0; i<nbItems; i++) {
            clusters[ availability[i]-1 ].add(i);
//            float currentItemPR = (float) this.profits[i] / this.weights[i] ;
//            PR_clusters[availability[i]-1].add(currentItemPR);
        }
    }
}
