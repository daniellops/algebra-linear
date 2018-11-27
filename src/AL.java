import java.util.ArrayList;
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
		
		double moduloInicial = 0;
		double moduloFinal = 0;
		
		do {
			moduloInicial = moduloFinal;
			QR aux = fatoracaoQR(matrizTransformada);
			matrizTransformada = Matrix.multiplicacao(aux.getR(), aux.getQ());
			qr.setQ(Matrix.multiplicacao(aux.getQ(), qr.getQ()));
			
			Complex[][] vetorDiagonal = Matrix.getVetorDiagonal(matrizTransformada);
			moduloFinal = Matrix.calcularModuloVetor(vetorDiagonal);
		} while (Math.abs(moduloFinal - moduloInicial) > 0.0000000001);
		
		List<AutoVetorAutoValor> autos;
		if(ehTransposta) {
			autos = autoVetorAutoValorMatrizDiagonal(matrizTransformada);
		} else {
			autos = autoVetorAutoValorMatrizComDentes(matrizTransformada);
		}
		autos = autoVetorAutoValorMatrizComDentes(matrizTransformada);
		for (AutoVetorAutoValor auto : autos) {
			Complex[][] autoVetor = auto.getAutoVetor();
			autoVetor = Matrix.multiplicacao(qr.getQ(), autoVetor);
			autoVetor = Matrix.multiplicacao(h.getH(), autoVetor);
			auto.setAutoVetor(autoVetor);
		}
		return autos;
	}
	
	public static List<AutoVetorAutoValor> autoVetorAutoValorMatrizDiagonal(Complex[][] a){
		List<AutoVetorAutoValor> autos = new ArrayList<>();
		for (int i = 0; i < a.length; i++) {
			AutoVetorAutoValor auto = new AutoVetorAutoValor();
			Complex autoValor = new Complex(a[i][i].re(), a[i][i].im());
			Complex[][] autoVetor = Matrix.criarMatrizZerada(a.length, 1);
			autoVetor[i][0] = new Complex(1, 0);
			auto.setAutoValor(autoValor);
			auto.setAutoValor(autoValor);
			autos.add(auto);
		}
		return autos;
	}
	
	public static List<AutoVetorAutoValor> autoVetorAutoValorMatrizComDentes(Complex[][] a){
		int i = 0;
		List<AutoVetorAutoValor> autos = new ArrayList<>();
		while (i < a.length) {
			List<AutoVetorAutoValor> avav = new ArrayList<>();
			
			if (i < a.length - 1 && a[i + 1][i].re() > 0.001) {
				Complex[][] sub = Matrix.getSubMatriz(a, i, i + 1, i, i + 1);
				List<AutoVetorAutoValor> avavDois = autoVetorAutoValorMatrizDoisPorDois(sub);
				
				AutoVetorAutoValor auto1 = new AutoVetorAutoValor();
				AutoVetorAutoValor auto2 = new AutoVetorAutoValor();
				
				auto1.setAutoValor(avavDois.get(0).getAutoValor());
				auto2.setAutoValor(avavDois.get(1).getAutoValor());
				
				Complex[][] autoVetor1 = Matrix.criarMatrizZerada(a.length, 1);
				Complex[][] autoVetor2 = Matrix.criarMatrizZerada(a.length, 1);

				autoVetor1[i][0] = avavDois.get(0).getAutoVetor()[0][0];
				autoVetor1[i + 1][0] = avavDois.get(0).getAutoVetor()[1][0];
				autoVetor2[i][0] = avavDois.get(1).getAutoVetor()[0][0];
				autoVetor2[i + 1][0] = avavDois.get(1).getAutoVetor()[1][0];
				
				auto1.setAutoVetor(autoVetor1);
				auto2.setAutoVetor(autoVetor2);
				
				avav.add(auto1);
				avav.add(auto2);
				i++;
			} else {
				AutoVetorAutoValor auto = new AutoVetorAutoValor();
				auto.setAutoValor(new Complex(a[i][i].re(), a[i][i].im()));
				
				Complex[][] autoVetor = Matrix.criarMatrizZerada(a.length, 1);
				autoVetor[i][0] = new Complex(1, 0);

				auto.setAutoVetor(autoVetor);
				avav.add(auto);
			}
			
			for (AutoVetorAutoValor av : avav) {
				int j = i - 1;
				while(j >= 0) {
					if (j != 0 && a[j][j - 1].re() > 0.001) {
						int k = j + 1;
						Complex[][] valor = Matrix.criarMatrizZerada(2, 1);
						while(k <= i) {
							valor[0][0] = valor[0][0].plus(a[j - 1][k].times(av.getAutoVetor()[k][0]));
							valor[1][0] = valor[1][0].plus(a[j][k].times(av.getAutoVetor()[k][0]));
							k++;
						}
						j--;
					} else {
						Complex valor = new Complex(1, 0);
						int k = j + 1;
						while(k <= i) {
							valor = valor.plus(a[j][k].times(av.getAutoVetor()[k][0]));
							k++;
						}
						valor = valor.divides(a[j][j].minus(a[i][i]));
						valor = valor.times(new Complex(-1, 0));
						av.getAutoVetor()[j][0] = valor;
					}
					j--;
				}
				autos.add(av);
			}
			
			i++;
		}
		return autos;
	}
	
	public static List<AutoVetorAutoValor> autoVetorAutoValorMatrizDoisPorDois(Complex[][] m){
		double a = 1.0;
		double b = (m[0][0].plus(m[1][1]).re()) * -1;
		double c = (m[0][0].times(m[1][1]).re()) - (m[0][1].times(m[1][0]).re());
		List<AutoVetorAutoValor> autos = new ArrayList<AutoVetorAutoValor>();
		double delta = Math.pow(b, 2) - (4 * a * c);
		
		if(delta < 0) {
			double raiz = Math.sqrt(Math.abs(delta));
			double real = Math.pow((m[0][0].plus(m[1][1]).re()), 2);
			Complex autoValor1 = new Complex(real, raiz);
			Complex autoValor2 = new Complex(real, -raiz);
			autoValor1 = autoValor1.divides(new Complex(2, 0));
			autoValor2 = autoValor1.divides(new Complex(2, 0));
			AutoVetorAutoValor auto1 = new AutoVetorAutoValor();
			AutoVetorAutoValor auto2 = new AutoVetorAutoValor();
			auto1.setAutoValor(autoValor1);
			auto2.setAutoValor(autoValor2);
			autos.add(auto1);
			autos.add(auto2);
		} else {
			double raiz = Math.sqrt(delta);
			double real = Math.pow((m[0][0].plus(m[1][1]).re()), 2);
			double autoValor1 = (real + raiz) / 2;
			double autoValor2 = (real - raiz) / 2;
			AutoVetorAutoValor auto1 = new AutoVetorAutoValor();
			AutoVetorAutoValor auto2 = new AutoVetorAutoValor();
			auto1.setAutoValor(new Complex(autoValor1, 0));
			auto2.setAutoValor(new Complex(autoValor2, 0));
			autos.add(auto1);
			autos.add(auto2);
		}
		
		for (AutoVetorAutoValor auto : autos) {
			Complex[][] autoVetor = new Complex[2][1];
			
			autoVetor[0][1] = new Complex(1, 0);
			autoVetor[0][0] = new Complex(-1, 0);
			auto.setAutoVetor(autoVetor);
		}
		
		return autos;
	}


	public static SVD svd(Complex[][] a) {
		Complex[][] aux1 = Matrix.multiplicacao(a, Matrix.transposta(a));
		Complex[][] aux2 = Matrix.multiplicacao(Matrix.transposta(a), a);
		Complex[][] abarra;
		if(aux1.length * aux1[0].length < aux2.length * aux2[0].length) {
			abarra = aux1;
		} else {
			abarra = aux2;
		}
		Complex[][] u = Matrix.criarMatrizZerada(a.length, a.length);
		Complex[][] s = Matrix.criarMatrizZerada(a.length, a[0].length);
		Complex[][] v = Matrix.criarMatrizZerada(a[0].length, a[0].length);
		
		List<AutoVetorAutoValor> autos = autoVetorAutoValorHouseHolderQR(abarra);
		int indice = 0;
		for (AutoVetorAutoValor auto : autos) {
			Complex[][] autoVetor = auto.getAutoVetor();
			for (int i = 0; i < autoVetor.length; i++) {
				v[i][indice] = new Complex(autoVetor[i][0].re(), autoVetor[i][0].im());
			}
			s[indice][indice] = new Complex(auto.getAutoValor().re(), auto.getAutoValor().im());
		}
		
		for (int i = 0; i < a[0].length; i++) {
			double rho = Math.sqrt(s[i][i].re());
			rho = 1 / rho;
			
			Complex[][] av = Matrix.multiplicacao(a, Matrix.getColuna(v, i));
			Complex[][] coluna = Matrix.multiplicacao(av, rho);
			
			for (int j = 0; j < coluna.length; j++) {
				u[j][i] = coluna[j][0];
			}
		}
		
		SVD svd = new SVD();
		svd.setU(u);
		svd.setS(s);
		svd.setV(v);
		return svd;
	}

}
