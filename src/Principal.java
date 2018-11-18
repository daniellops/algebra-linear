
public class Principal {
	public static void main(String[] args) {
		Complex [][] C = {{new Complex(1, 0), new Complex(2, 0), new Complex(3, 0)}, {new Complex(2, 0), new Complex(2, 0), new Complex(0, 0)}, {new Complex(3, 0), new Complex(0, 0), new Complex(3, 0)}};
		//xxdouble [][] C = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
//		HouseHolder hh = AL.householder(C);
//		AL.fatoracaoQR(hh.getMatrizTridiagonalizada());
		AL.autoVetorAutoValorHouseHolderQR(C);
	}
}
