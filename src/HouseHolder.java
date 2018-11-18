
public class HouseHolder {
	
	private Complex[][] h;
	private Complex[][] matrizTridiagonalizada;
	
	public HouseHolder(Complex[][] matrizTridiagonalizada) {
		super();
		this.matrizTridiagonalizada = matrizTridiagonalizada;
		Complex[][] h = new Complex[matrizTridiagonalizada.length][matrizTridiagonalizada[0].length];
		
		this.h = Matrix.identidade(h.length);
	}

	public Complex[][] getH() {
		return h;
	}

	public void setH(Complex[][] h) {
		this.h = h;
	}

	public Complex[][] getMatrizTridiagonalizada() {
		return matrizTridiagonalizada;
	}

	public void setMatrizTridiagonalizada(Complex[][] matrizTridiagonalizada) {
		this.matrizTridiagonalizada = matrizTridiagonalizada;
	}
	
}
