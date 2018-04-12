package proposers;

import common.LogDouble;

public interface Proposer {

	// Executes an actual perturbation of the state parameter(s) and cache it for
	// restoring purpose
	public boolean cacheAndPerturb();

	// Clears the cached (previous) state of perturbed parameters when a proposed
	// state has been accepted.
	public void clearCache();

	// Restores the cached (previous) state of perturbed parameters when a proposed
	// state has been rejected.
	public void restoreCache();

	// Returns the "forward" probability density Q(x';x) for obtaining the new value
	// x' given the old value x.
	public LogDouble getForwardDensity();

	// Returns the "backward" probability density Q(x;x') for obtaining the old
	// value x given the new value x'.
	public LogDouble getBackwardDensity();

	// The ratio between the "backward" and "forward" proposal densities
	// respectively.
	// It will returns the ratio Q(x;x')/Q(x';x) for the old state x and the new
	// state x',
	public LogDouble getDensityRatio();

	// This method will return the proposer name
	public String getProposerName();

	// This method will return the object of statistics associated with this
	// proposer
	public ProposerStatistics getProposerStatistics();

	// Returns whether this proposer is active or not; return true if enabled; false
	// if disabled.
	public boolean isEnabled();

	// Controls whether this proposer is active or not. isActive true to enable;
	// false to disable.
	public void setEnabled(boolean isActive);

	// Returns true if the proposal is indeed valid.
	public boolean hasValidProposal();

}
