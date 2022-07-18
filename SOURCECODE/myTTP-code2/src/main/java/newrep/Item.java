package newrep;

public class Item {
	public int weight;
	public int profit;
	public boolean isSelected;
	public int itemId;
	public double[] profitPerWc;
	public boolean isChecked;

	public Item(int profit, int weight, int itemId){
		this.profit = profit;
		this.weight = weight;
		this.itemId = itemId;
		this.isSelected = false;
		this.isChecked = false;
	}
	public Item(int itemId,boolean isChecked){

		this.itemId = itemId;
//		this.isSelected = false;
		this.isChecked = isChecked;
	}
}
