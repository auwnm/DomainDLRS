package dlModel;

import java.util.ArrayList;
import phylogeny.Node;
import phylogeny.Tree;

// Birth-death Model
// Code adapted from Rasmussen and Kellis et. al. (2011) 

public class BirthDeathProbs {

	private static final double INFINITY = Double.POSITIVE_INFINITY;

	// Rasmussen and Kellis et. al. (2011)
	// P(0,t) Probability that there will be no survivors after time t
	public static double P0(double birthRate, double deathRate, double time) {
		double l = birthRate;
		double u = deathRate;
		double t = time;
		double r = l - u;
		double p0;

		if (time == 0.0)
			return 0.0;

		if (birthRate == deathRate) {
			double lt = l * t;
			double lt1 = 1.0 + lt;
			p0 = lt / lt1;
		} else {
			double ert = Math.exp(-r * t);
			double luert = l - u * ert;
			p0 = (u - u * ert) / luert;
		}
		return p0;
	}

	// Rasmussen and Kellis et. al. (2011)
	// P(1,t) Probability that there will be only 1 survivor after time t
	public static double P1(double birthRate, double deathRate, double time) {
		double l = birthRate;
		double u = deathRate;
		double t = time;
		double r = l - u;
		double p1;

		if (time == 0.0)
			return 1.0;

		if (birthRate == deathRate) {
			double lt = l * t;
			double lt1 = 1.0 + lt;
			p1 = 1.0 / lt1 / lt1;
		} else {
			double ert = Math.exp(-r * t);
			double luert = l - u * ert;
			p1 = r * r * ert / luert / luert;
		}
		return p1;
	}

	public static double[] calcExtinctionProb(Tree tree, double[] BDRates) {

		double p0, p1;
		double birthRate = BDRates[0];
		double deathRate = BDRates[1];
		double[] extinctionProb = new double[tree.nnodes];

		// get nodes in post order
		ArrayList<Node> nodes = tree.getTreePostList();

		for (int i = 0; i < tree.nnodes; i++) {
			Node node = nodes.get(i);

			if (node.isLeaf()) {
				extinctionProb[node.id] = -INFINITY;
			} else {
				double prod = 0.0;
				for (int j = 0; j < node.nchildren; j++) {
					Node child = node.children.get(j);

					// Compute Po and P1
					// double t = tree.at[child.id];
					double t = tree.vt[child.parent.id] - tree.vt[child.id];

					double dc = Math.exp(extinctionProb[child.id]);

					p0 = BirthDeathProbs.P0(birthRate, deathRate, t);
					p1 = BirthDeathProbs.P1(birthRate, deathRate, t);

					// Probability of Extinction
					prod += Math.log(p0 + dc * p1 / (1.0 - (birthRate / deathRate) * p0 * dc));
				}

				extinctionProb[node.id] = prod;
			}
		}
		return extinctionProb;
	}

	public static void showExtinctionProb(Tree htree, double[] extinctionPrb) {
		String str;
		ArrayList<Node> orderedNodes = new ArrayList<Node>();
		htree.getTreePostOrder(orderedNodes, htree.root);

		System.out.println("Death Probabilities:");
		for (int i = 0; i < htree.nnodes; i++) {
			if (orderedNodes.get(i).isLeaf())
				str = String.format("@ Leaf %-3s %-25s", orderedNodes.get(i).id, orderedNodes.get(i).name);
			else
				str = String.format("@ Node %-3s parent(%-3s,%-3s)\t", orderedNodes.get(i).id,
						orderedNodes.get(i).children.get(0).id, orderedNodes.get(i).children.get(1).id);
			System.out.println(str + "\tP(Extinction) = " + Math.exp(extinctionPrb[orderedNodes.get(i).id]));
		}

	}

	/*
	 * // Probability function for computing probability of death of a lineage // at
	 * a certain point of time in an sliced part of the host tree. public static
	 * double deathAttime(double p0, double p1, double lu, double dy, double to) {
	 * return p0 + ( (dy * p1) / (1.0 - lu * p0 * dy) ); }
	 * 
	 * // Returns the probability of 1 gene giving rise to ngenes after time 'time'
	 * double birthDeathCount(int ngenes, float time, double birthRate, float
	 * deathRate) { final double l = birthRate; final double u = deathRate;
	 * 
	 * if (birthRate == deathRate) { final double ut = time / (1.0 / birthRate +
	 * time); if (ngenes == 0) return ut; return (1.0 - ut)*(1.0 - ut) * ipow(ut,
	 * ngenes-1); }
	 * 
	 * final double r = l - u; final double a = u / l;
	 * 
	 * final double ut = (1.0 - exp(-r*time)) / (1.0 - a * exp(-r*time)); final
	 * double p0 = a*ut;
	 * 
	 * if (ngenes == 0) return p0;
	 * 
	 * // p(0, t) = ut // p(1, t) = ... // (1.0 - p0)*(1.0 - ut) * ut^{ngenes-1}
	 * return (1.0 - p0)*(1.0 - ut) * ipow(ut, ngenes-1); }
	 * 
	 * // returns the probability of 'start' genes giving rise to 'end' genes after
	 * // time 'time' // slower more stable computation double
	 * birthDeathCountsSlow(int start, int end, float time, float birth, float
	 * death) { final double l = birth; final double u = death; final double r = l -
	 * u; final double a = u / l;
	 * 
	 * final double ertime = Math.exp(-r*time); final double ut = (1.0 - ertime) /
	 * (1.0 - a * ertime); final double p0 = a*ut;
	 * 
	 * // all 'start' genes die out if (end == 0) { return ipow(p0, start); }
	 * 
	 * final int iter = (start < end) ? start : end; double p = 0.0;
	 * 
	 * for (int j=0; j<=iter; j++) { p += fchoose(start, j) * fchoose(start + end -
	 * j - 1, start - 1) * ipow(p0, start-j) * ipow(ut, end-j) * ipow(1 - p0 - ut,
	 * j); }
	 * 
	 * // do not allow invalid values to propogate if (isnan(p) || isinf(p) || p >
	 * 1.0) { printf("p=%e genes=(%d, %d) b=%f d=%f t=%f\n", p, start, end, birth,
	 * death, time); fflush(stdout); assert(0); } return p; }
	 * 
	 * 
	 * // returns the probability of 'start' genes giving rise to 'end' genes after
	 * // time 'time' // much faster computation than birthDeathCounts2 double
	 * birthDeathCounts(int start, int end, float time, float birth, float death) {
	 * if (start == 0) { if (end == 0) return 1.0; else return 0.0; }
	 * 
	 * final double ertime = exp((birth-death)*time); final double tmp =
	 * (ertime-1.0) / (birth*ertime - death); final double a = death * tmp; final
	 * double b = birth * tmp;
	 * 
	 * // all 'start' genes die out if (end == 0) { return ipow(a, start); }
	 * 
	 * // compute base case double f = ipow(a, start) * ipow(b, end); if (start > 1)
	 * f *= (start + end - 1); for (int k=2; k<start; k++) f *= (start + end - k) /
	 * double(k);
	 * 
	 * 
	 * double p = f; double x = start; double y = end; double z = start + end - 1;
	 * final double oneab = 1.0 - a - b; final int iter = (start < end) ? start :
	 * end; for (int j=1; j<=iter; j++) { f *= (oneab * x * y / (j * a * b * z)); p
	 * += f; x--; y--; z--; }
	 * 
	 * if (p < 0.0) p = 0.0;
	 * 
	 * 
	 * if (p > 1.0) { // resort to a slower more stable function return
	 * birthDeathCountsSlow(start, end, time, birth, death); }
	 * 
	 * return p; }
	 * 
	 * 
	 * 
	 * // Probability of no birth from 'n' lineages starting at time 0, // evolving
	 * until time 'T' with 'birth' and 'death' rates // for a reconstructed process.
	 * double probNoBirth(int n, float T, float birth, float death) { if (birth ==
	 * 0.0) return 1.0; else if (birth == death) { return 1.0 / ipow(1.0 + birth *
	 * T, n); }
	 * 
	 * final double r = birth - death; final double exprt = exp(-r * T);
	 * 
	 * return pow(1.0 - (birth*(1.0 - exprt)) / (birth - death * exprt), n); }
	 * 
	 * 
	 * 
	 * 
	 * 
	 * //===========================================================================
	 * == // sampling
	 * 
	 * 
	 * // Probability density for for next birth at time 't' given // 'n' lineages
	 * starting at time 0, evolving until time 'T' with a // 'birth' and 'death'
	 * rates for a reconstructed process. double birthWaitTime(float t, int n, float
	 * T, float birth, float death) { if (birth == death) { final double t2 = t - T;
	 * final double nl = 1.0 / birth; return birth * n * ipow(-nl+t2, n) /
	 * ipow(-nl-T, n) / (1.0-birth*t2); }
	 * 
	 * final double r = birth - death; final double a = death / birth;
	 * 
	 * return n * r * exp(-n*r*t) * \ pow(1.0 - a * exp(-r * (T - t)), n-1) / \
	 * pow(1.0 - a * exp(-r * T), n); }
	 * 
	 * // Probability density for for next birth at time 't' given // 'n'=1 lineages
	 * starting at time 0, evolving until time 'T' with a // 'birth' and 'death'
	 * rates for a reconstructed process. double birthWaitTime1(float t, float T,
	 * float birth, float death) { // special case if (birth == death) { final
	 * double t2 = t - T; final double nl = 1.0 / birth; return birth * (-nl+t2) /
	 * (-nl-T) / (1.0-birth*t2); }
	 * 
	 * final double r = birth - death; final double a = death / birth;
	 * 
	 * return r * exp(-r*t) / (1.0 - a * exp(-r * T)); }
	 * 
	 * 
	 * // numerator for birthWaitTime double birthWaitTimeNumer(float t, int n,
	 * float T, float birth, float death, double denom) { const double r = birth -
	 * death; const double a = death / birth;
	 * 
	 * return n * r * exp(-n*r*t) * pow(1.0 - a * exp(-r * (T - t)), n-1) / denom; }
	 * 
	 * // denominator for birthWaitTime double birthWaitTimeDenom(int n, float T,
	 * float birth, float death) { const double r = birth - death; const double a =
	 * death / birth;
	 * 
	 * return pow(1.0 - a * exp(-r * T), n); }
	 * 
	 * 
	 * 
	 * // Probability density for for next birth at time 't' given // 'n'=1 lineages
	 * starting at time 0, evolving until time 'T' with a // 'birth' and 'death'
	 * rates for a reconstructed process. // The denominator 'denom' must be
	 * precomputed double birthWaitTimeNumer1(float t, float T, float birth, float
	 * death, float denom) { const double r = birth - death; return r * exp(-r*t) /
	 * denom; }
	 * 
	 * // Computes the denominator for birthWaitTime1 double
	 * birthWaitTimeDenom1(float T, float birth, float death) { const double r =
	 * birth - death; const double a = death / birth; return 1.0 - a * exp(-r * T);
	 * }
	 * 
	 * 
	 * 
	 * 
	 * // Sample the next birth event from a reconstructed birthdeath process. //
	 * Let there be 'n' lineages at time 0 that evolve until time 'T' with //
	 * 'birth' and 'death' rates. // Conditioned that a birth will occur double
	 * sampleBirthWaitTime(int n, float T, float birth, float death) {
	 * 
	 * // TODO: could make this more efficient
	 * 
	 * if (birth == death) { double start_y = birthWaitTime(0, n, T, birth, death);
	 * double end_y = birthWaitTime(T, n, T, birth, death); double M = max(start_y,
	 * end_y);
	 * 
	 * while (true) { double t = frand(T); double f = birthWaitTime(t, n, T, birth,
	 * death);
	 * 
	 * if (frand() <= f / M) return t; } } else { // uses rejection sampling double
	 * denom = birthWaitTimeDenom(n, T, birth, death); double start_y =
	 * birthWaitTimeNumer(0, n, T, birth, death, denom); double end_y =
	 * birthWaitTimeNumer(T, n, T, birth, death, denom); double M = max(start_y,
	 * end_y);
	 * 
	 * while (true) { double t = frand(T); double f = birthWaitTimeNumer(t, n, T,
	 * birth, death, denom);
	 * 
	 * if (frand() <= f / M) return t; } } }
	 * 
	 * 
	 * 
	 * // Sample the next birth event from a reconstructed birthdeath process. //
	 * Let there be 'n'=1 lineages at time 0 that evolve until time 'T' with //
	 * 'birth' and 'death' rates. // Conditioned that a birth will occur double
	 * sampleBirthWaitTime1(float T, float birth, float death) { // TODO: could make
	 * this much more efficient
	 * 
	 * if (birth == death) { // uses rejection sampling double start_y =
	 * birthWaitTime1(0, T, birth, death); double end_y = birthWaitTime1(T, T,
	 * birth, death); double M = max(start_y, end_y);
	 * 
	 * while (true) { double t = frand(T); double f = birthWaitTime1(t, T, birth,
	 * death);
	 * 
	 * if (frand() <= f / M) return t; }
	 * 
	 * } else { // uses rejection sampling double denom = birthWaitTimeDenom1(T,
	 * birth, death); double start_y = birthWaitTimeNumer1(0, T, birth, death,
	 * denom); double end_y = birthWaitTimeNumer1(T, T, birth, death, denom); double
	 * M = max(start_y, end_y);
	 * 
	 * while (true) { double t = frand(T); double f = birthWaitTimeNumer1(t, T,
	 * birth, death, denom);
	 * 
	 * if (frand() <= f / M) return t; } } }
	 */

}
