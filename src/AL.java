import java.util.List;

public class AL {

	public static int getNumeroDeLinhas(Complex[][] a) {
		return a.length;
	}

	public static int getNumeroDeColunas(Complex[][] a) {
		return a[0].length;
	}

	public static HouseHolder householder(Complex[][] matriz) {
		Complex[][] matrizTridiagonalizada = Matrix.clone(matriz.clone());
		HouseHolder hh = new HouseHolder(matrizTridiagonalizada);

		int ordem = getNumeroDeLinhas(matrizTridiagonalizada);

		for (int k = 0; k < ordem - 2; k++) {
			Complex[][] p = Matrix.getColuna(matrizTridiagonalizada, k);
			for (int i = 0; i <= k; i++) {
				p[i][0] = new Complex(0, 0);
			}

			double np = Matrix.calcularModuloVetor(p);

			Complex[][] plinha = Matrix.criarMatrizZerada(ordem, 1);
			int alfa = (int) Math.signum(matrizTridiagonalizada[k][k].re());
			plinha[k + 1][0] = new Complex(np * alfa, 0);

			Complex[][] n = Matrix.subtracao(p, plinha);

			Complex[][] nzinho = Matrix.normalizar(n);

			Complex[][] matrizIdentidade = Matrix.identidade(ordem);
			Complex[][] x = Matrix.multiplicacao(nzinho, Matrix.multiplicacao(Matrix.transposta(nzinho), 2.0));
			Complex[][] matrizH = Matrix.subtracao(matrizIdentidade, x);

			hh.setH(Matrix.multiplicacao(matrizH, hh.getH()));
			matrizTridiagonalizada = Matrix.multiplicacao(matrizH,
					Matrix.multiplicacao(matrizTridiagonalizada, Matrix.transposta(matrizH)));
		}
		hh.setMatrizTridiagonalizada(matrizTridiagonalizada);
		return hh;
	}

	public static QR fatoracaoQR(Complex[][] matriz) {
		Complex[][] r = Matrix.clone(matriz);
		int ordem = r.length;
		Complex[][] q = Matrix.identidade(ordem);

		for (int c = 0; c < ordem - 1; c++) {
			Complex[][] v = Matrix.getColuna(r, c);
			for (int i = 0; i < c; i++) {
				v[i][0] = new Complex(0, 0);
			}

			double nv = Matrix.calcularModuloVetor(v);

			v[c][0] = new Complex(v[c][0].re() - nv, 0);

			double nnzinho = Matrix.calcularModuloVetor(v);
			Complex[][] nzinho = Matrix.multiplicacao(v, 1.0 / nnzinho);

			Complex[][] matrizIdentidade = Matrix.identidade(ordem);
			Complex[][] x = Matrix.multiplicacao(nzinho, Matrix.multiplicacao(Matrix.transposta(nzinho), 2.0));
			Complex[][] Hc = Matrix.subtracao(matrizIdentidade, x);

			r = Matrix.multiplicacao(Hc, r);
			q = Matrix.multiplicacao(Hc, q);
		}

		QR qr = new QR();
		qr.setQ(q);
		qr.setR(r);
		return qr;
	}

	public static List<AutoVetorAutoValor> autoVetorAutoValorHouseHolderQR(Complex[][] matriz) {
		boolean ehTransposta = Matrix.ehTransposta(matriz);
		
		HouseHolder h = householder(matriz);
		QR qr = new QR();
		qr.setQ(Matrix.identidade(h.getMatrizTridiagonalizada().length));
		Complex[][] matrizTransformada = h.getMatrizTridiagonalizada();
		
		double modulo = 0;
		
		do {
			QR aux = fatoracaoQR(matrizTransformada);
			matrizTransformada = Matrix.multiplicacao(Matrix.transposta(aux.getQ()), Matrix.multiplicacao(matrizTransformada, aux.getQ()));
			qr.setQ(Matrix.multiplicacao(aux.getQ(), qr.getQ()));
			
			Complex[][] vetorDiagonal = Matrix.getVetorDiagonal(matrizTransformada);
			modulo = Matrix.calcularModuloVetor(vetorDiagonal);
		} while (modulo > 0.0000000001);
		for (int i = 0; i < matrizTransformada.length; i++) {
			for (int j = 0; j < matrizTransformada.length; j++) {
				System.out.println(matrizTransformada[i][j]);
			}
		}
		
		if(ehTransposta) {
			
		} else {
			
		}
		return null;
	}
	
	public static List<AutoVetorAutoValor> autoVetorAutoValorMatrizDiagonal(Complex[][] a){
		for (int i = 0; i < a.length; i++) {
			AutoVetorAutoValor auto = new AutoVetorAutoValor();
			Complex autoValor = new Complex(a[i][i].re(), a[i][i].im());
			Complex[][] autoVetor = Matrix.criarMatrizZerada(a.length, 1);
			autoVetor[i][0] = new Complex(1, 0);
			auto.setAutoValor(autoValor);
			auto.setAutoValor(autoValor);
		}
		return null;
	}
	
	public static List<AutoVetorAutoValor> autoVetorAutoValorMatrizComDentes(Complex[][] a){
		return null;
	}

}
