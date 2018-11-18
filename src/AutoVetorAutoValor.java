
public class AutoVetorAutoValor {

	private Complex[][] autoVetor;
	private Complex autoValor;

	public AutoVetorAutoValor() {
		super();
	}

	public AutoVetorAutoValor(Complex[][] autoVetor, Complex autoValor) {
		this.autoValor = autoValor;
		this.autoVetor = autoVetor;
	}

	public Complex[][] getAutoVetor() {
		return autoVetor;
	}

	public void setAutoVetor(Complex[][] autoVetor) {
		this.autoVetor = autoVetor;
	}

	public Complex getAutoValor() {
		return autoValor;
	}

	public void setAutoValor(Complex autoValor) {
		this.autoValor = autoValor;
	}

}
