package dlModel;

import java.util.Arrays;

import phylogeny.Reconciliation;

public class Realization {

	public double[] vertexTimesCache;
	public double[] impVertexTimesCache;

	public void Cache(Reconciliation R) {

		if (R.guest.vt != null)
			this.vertexTimesCache = Arrays.copyOf(R.guest.vt, R.guest.vt.length);

		if (R.guestPrime.vt != null)
			this.impVertexTimesCache = Arrays.copyOf(R.guestPrime.vt, R.guestPrime.vt.length);
	}

	public void RestoreCache(Reconciliation R) {

		if (vertexTimesCache != null)
			R.guest.vt = this.vertexTimesCache;
		if (impVertexTimesCache != null)
			R.guestPrime.vt = this.impVertexTimesCache;
	}

	public void clearCache() {
		this.vertexTimesCache = null;
		this.impVertexTimesCache = null;
	}

}
