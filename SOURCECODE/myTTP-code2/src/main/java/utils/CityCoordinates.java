package utils;

/*
 * city coordinates
 * */
public class CityCoordinates {
    private double x;
    private double y;

    public CityCoordinates(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public CityCoordinates() {
    }

    @Override
    public String toString() {
        return "(" + String.format("%.2f", this.x) + "," +
                     String.format(".2f", this.y) + ")";
    }
// getters
    public double getX() {
        return x;
    }
    public double getY() {
        return y;
    }
//  setters
    public void setX(double x) {
        this.x = x;
    }
    public void setY(double y) {
        this.y = y;
    }

    public double distanceEuclid(CityCoordinates c2){
        double P = this.x - c2.x;
        double Q = this.y - c2.y;
        return Math.sqrt(P*P + Q*Q);
    }
}
