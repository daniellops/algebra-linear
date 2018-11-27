
public class Matrix {

	public static Complex[][] identidade(int ordem) {
		Complex[][] c = new Complex[ordem][ordem];
		for (int i = 0; i < ordem; i++) {
			for (int j = 0; j < ordem; j++) {
				double re = 0;
				if (i == j)
					re = 1;
				c[i][j] = new Complex(re, 0);
			}
		}
		return c;
	}

	public static Complex[][] clone(Complex[][] a) {
		Complex[][] c = new Complex[a.length][a[0].length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {
				c[i][j] = new Complex(a[i][j].re(), a[i][j].im());
			}
		}
		return c;
	}

	public static Complex[][] criarMatrizZerada(int linhas, int colunas) {
		Complex[][] c = new Complex[linhas][colunas];
		for (int i = 0; i < c.length; i++) {
			for (int j = 0; j < c[0].length; j++) {
				c[i][j] = new Complex(0, 0);
			}
		}
		return c;
	}

	public static Complex[][] soma(Complex[][] a, Complex[][] b) {
		Complex[][] c = new Complex[a.length][a[0].length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {
				c[i][j] = a[i][j].plus(b[i][j]);
			}
		}
		return c;
	}

	public static Complex[][] subtracao(Complex[][] a, Complex[][] b) {
		Complex[][] c = new Complex[a.length][a[0].length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {
				c[i][j] = a[i][j].minus(b[i][j]);
			}
		}
		return c;
	}

	public static Complex[][] multiplicacao(Complex[][] a, Complex[][] b) {
		Complex[][] c = new Complex[a.length][b[0].length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b[0].length; j++) {
				for (int k = 0; k < a[0].length; k++) {
//					Cij = Cij + Aik × Bkj
					if (c[i][j] == null)
						c[i][j] = new Complex(0, 0);
					c[i][j] = c[i][j].plus(a[i][k].times(b[k][j]));
				}
			}
		}
		return c;
	}

	public static Complex[][] multiplicacao(Complex[][] a, double b) {
		Complex[][] c = new Complex[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {
				c[i][j] = a[i][j].times(new Complex(b, 0));
			}
		}
		return c;
	}

	public static Complex[][] transposta(Complex[][] a) {
		Complex[][] c = new Complex[a[0].length][a.length];
		for (int i = 0; i < c.length; i++) {
			for (int j = 0; j < c[0].length; j++) {
				c[i][j] = new Complex(a[j][i].re(), a[j][i].im());
			}
		}
		return c;
	}

	public static Complex[][] getColuna(Complex[][] a, int j) {
		Complex[][] c = new Complex[a.length][1];
		for (int i = 0; i < c.length; i++) {
			c[i][0] = new Complex(a[i][j].re(), a[i][j].im());
		}
		return c;
	}

	public static Complex[][] getLinha(Complex[][] a, int i) {
		Complex[][] c = new Complex[1][a[0].length];
		for (int j = 0; j < c.length; j++) {
			c[0][j] = new Complex(a[i][j].re(), a[i][j].im());
		}
		return c;
	}
	
	public static Complex[][] getVetorDiagonal(Complex[][] a){
		Complex[][] c = new Complex[a.length][1];
		for (int i = 0; i < a.length; i++) {
			c[i][0] = new Complex(a[i][i].re(), a[i][i].im());
		}
		return c;
	}

	public static Complex[][] getSubMatriz(Complex[][] a, int i0, int i1, int j0, int j1) {
		Complex[][] c = new Complex[i1 - i0 + 1][j1 - j0 + 1];
		for (int i = 0; i < c.length; i++) {
			for (int j = 0; j < c.length; j++) {
				c[i][j] = new Complex(a[i0 + i][j0 + j].re(), a[i0 + i][j0 + j].im());
			}
		}
		return c;
	}
	
	public static boolean ehTransposta(Complex[][] a) {
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < i; j++) {
				if (a[i][j].re() != a[i][j].re() || a[i][j].im() != a[i][j].im())
					return false;
			}
		}
		return true;
	}

	public static double calcularModuloVetor(double[] vetor) {
		double modulo = 0;
		for (int i = 0; i < vetor.length; i++) {
			modulo += vetor[i] * vetor[i];
		}
		return Math.sqrt(modulo);
	}

	public static double calcularModuloVetor(Complex[][] vetor) {
		double modulo = 0;
		for (int i = 0; i < vetor.length; i++) {
			for (int j = 0; j < vetor[0].length; j++) {
				modulo += vetor[i][j].re() * vetor[i][j].re();
			}
		}
		return Math.sqrt(modulo);
	}

	public static Complex[][] normalizar(Complex[][] vetor) {
		return multiplicacao(vetor, 1.0 / calcularModuloVetor(vetor));
	}

}
